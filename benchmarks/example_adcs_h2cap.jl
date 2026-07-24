# example_adcs_h2cap.jl — EXPLORATORY variant of the ADCS SDP co-design that swaps
# the H∞ control side for an H₂ control SDP WITH A HARD CONTROL-EFFORT CAP.
#
# THIS FILE CHANGES NOTHING ELSE (example_adcs.jl, example_adcs_hinf.jl, the paper,
# and src/ are all untouched).  It is a standalone experiment answering: "on the SDP
# path, can we bound the actuator effort with a hard LMI cap that a Riccati cannot
# express, and lower the effort weight R accordingly?"
#
# Control side:  H₂ state-feedback SDP in (Σ, Z₀, L), K = L Σ⁻¹:
#     min tr(QΣ) + tr(R Z₀)
#     s.t. [Z₀ L; Lᵀ Σ] ⪰ 0                         (Z₀ ⪰ K Σ Kᵀ, input second moment)
#          [Σ−W  AΣ+BL; (AΣ+BL)ᵀ  Σ] ⪰ 0            (Lyapunov LMI, closed-loop cov)
#          tr(Z₀) ≤ u_cap                            (HARD control-effort cap — SDP only)
#   with R = 0.1·I (softened, because the hard cap now regulates the effort).  The
#   Riccati/DARE can only *weight* effort through R; it cannot *cap* E[uᵀu] — the cap
#   is a constraint on the feasible set, so no choice of R reproduces it.  The
#   envelope gradients are the standard covariance-SDP duals ∂J/∂A = −2 S₁₂ Σᵀ,
#   ∂J/∂B = −2 S₁₂ Lᵀ (the cap carries no θ, so it adds no gradient term).
#
# Estimation side, budget, box, weights, single-origin linearization: identical to
# example_adcs_hinf.jl (Q=diag(4,4,4,1,1,1), c_b=5, dense H₂ estimation SDP).
#
# Run:  julia --project benchmarks/example_adcs_h2cap.jl

include("../src/LQR.jl")     # dare, dlyap (Riccati cross-checks)
include("../src/BFGS.jl")    # cached_objective, verify_gradient, multistart_design
using JuMP, Clarabel, LinearAlgebra, ForwardDiff, Printf

# ── Spacecraft parameters (identical to example_adcs_hinf.jl) ──────────────────
const Jx = 4.0; const Jy = 6.0; const Jz = 5.0
const Jmat = Matrix(Diagonal([Jx, Jy, Jz])); const Jinv = inv(Jmat)
const dt_s = 0.1
const n = 6; const m = 3; const ry = 6

const v_st0 = 0.20; const v_gyro0 = 0.20
const Σ_w  = Matrix(Diagonal([1e-6, 1e-6, 1e-6, 1e-4, 1e-4, 1e-4]))
const Qmat = Matrix(Diagonal([4.0, 4.0, 4.0, 1.0, 1.0, 1.0]))
const Rw   = 0.1                              # softened effort weight (cap does the bounding)
const Rmat = Rw * Matrix(I, m, m)
const Cmat = Matrix(1.0I, ry, n)
Vmat(θ) = Matrix(Diagonal([v_st0/θ[2]^2, v_st0/θ[2]^2, v_st0/θ[2]^2,
                           v_gyro0/θ[3]^2, v_gyro0/θ[3]^2, v_gyro0/θ[3]^2]))

function gyro_jac(ω)
    ωx, ωy, ωz = ω
    [0.0          (Jz-Jy)*ωz   (Jz-Jy)*ωy;
     (Jx-Jz)*ωz   0.0          (Jx-Jz)*ωx;
     (Jy-Jx)*ωy   (Jy-Jx)*ωx   0.0]
end

function model_lin(θ, ω_op)
    T = eltype(θ); e_rw = θ[1]
    A = Matrix{T}(I, n, n)
    A[1,4] = dt_s; A[2,5] = dt_s; A[3,6] = dt_s
    A[4:6, 4:6] .= I(3) .- dt_s .* (Jinv * gyro_jac(ω_op))
    B = zeros(T, n, m); B[4:6, 1:3] .= dt_s .* e_rw .* Jinv
    return A, B
end

const ω_ens = [[0.0, 0.0, 0.0]]              # single linearization about the origin
const M_ens = length(ω_ens)

# ── CONTROL inner value: H₂ SDP with hard effort cap tr(Z₀) ≤ u_cap ────────────
# Returns (J, ∇θJ, J_state, J_effort, effort=tr(Z₀), status). u_cap=Inf ⇒ plain H₂.
function ctrl_h2(θ, ω_op; u_cap = Inf)
    A, B = model_lin(θ, ω_op)
    mdl = Model(Clarabel.Optimizer); set_silent(mdl)
    @variable(mdl, Σ[1:n, 1:n], Symmetric)
    @variable(mdl, Z0[1:m, 1:m], Symmetric)
    @variable(mdl, L[1:m, 1:n])
    @constraint(mdl, Symmetric(Σ .- 1e-7 .* Matrix(I, n, n)) in PSDCone())
    @constraint(mdl, Symmetric([Z0 L; permutedims(L) Σ]) in PSDCone())
    clyap = @constraint(mdl, Symmetric([Σ .- Σ_w  (A*Σ .+ B*L);
                                        permutedims(A*Σ .+ B*L)  Σ]) in PSDCone())
    isfinite(u_cap) && @constraint(mdl, tr(Z0) <= u_cap)
    @objective(mdl, Min, tr(Qmat * Σ) + tr(Rmat * Z0))
    optimize!(mdl)
    st = termination_status(mdl)
    (st in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.NUMERICAL_ERROR)) ||
        return (1e8, fill(0.0, length(θ)), 1e8, 1e8, NaN, st)
    Σv = value.(Σ); Lv = value.(L); Z0v = value.(Z0)
    J = objective_value(mdl)
    Jstate = tr(Qmat * Σv); Jeff = tr(Rmat * Z0v); eff = tr(Z0v)
    S = dual(clyap); S12 = S[1:n, n+1:2n]
    gA = -2 .* (S12 * permutedims(Σv))
    gB = -2 .* (S12 * permutedims(Lv))
    Aj = ForwardDiff.jacobian(t -> vec(model_lin(t, ω_op)[1]), θ)
    Bj = ForwardDiff.jacobian(t -> vec(model_lin(t, ω_op)[2]), θ)
    ∇ = [tr(gA' * reshape(Aj[:, k], n, n)) + tr(gB' * reshape(Bj[:, k], n, m))
         for k in eachindex(θ)]
    return J, ∇, Jstate, Jeff, eff, st
end

# Unconstrained H₂ control cost via the DARE (cross-check for the SDP).
function ctrl_dare(θ, ω_op)
    A, B = model_lin(θ, ω_op)
    P = dare(A, B, Qmat, Rmat); tr(P * Σ_w)
end

# ── ESTIMATION inner value: dense H₂ covariance SDP (P^ε form) ─────────────────
function est_sdp_dense(A, V; withgrad = true)
    mdl = Model(Clarabel.Optimizer); set_silent(mdl)
    @variable(mdl, P[1:n, 1:n], Symmetric)
    @variable(mdl, F[1:n, 1:ry])
    @variable(mdl, Z[1:ry, 1:ry], Symmetric)
    @constraint(mdl, Symmetric([Z permutedims(F); F P]) in PSDCone())
    lmi2 = @constraint(mdl, Symmetric([P .- Qmat  (permutedims(A)*P .- permutedims(Cmat)*permutedims(F));
                                       (P*A .- F*Cmat)  P]) in PSDCone())
    @objective(mdl, Min, tr(Σ_w * P) + tr(V * Z))
    optimize!(mdl)
    st = termination_status(mdl)
    (st in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.NUMERICAL_ERROR)) ||
        return (1e8, zeros(n, n), zeros(ry, ry), st)
    J = objective_value(mdl); Zv = value.(Z)
    S = dual(lmi2); S12 = S[1:n, n+1:2n]
    (J, -2 .* (value.(P) * permutedims(S12)), Matrix(Zv), st)
end

function est_vg(θ, ω_op)
    A, _ = model_lin(θ, ω_op); V = Vmat(θ)
    J, gA, gV, st = est_sdp_dense(A, V)
    J >= 1e8 && return (1e8, fill(0.0, length(θ)))
    Aj = ForwardDiff.jacobian(t -> vec(model_lin(t, ω_op)[1]), θ)
    Vj = ForwardDiff.jacobian(t -> vec(Matrix(Vmat(t))), θ)
    ∇ = [tr(gA' * reshape(Aj[:, k], n, n)) + tr(gV' * reshape(Vj[:, k], ry, ry))
         for k in eachindex(θ)]
    (J, ∇)
end
kalman_cost(A, V) = tr(Qmat * dare(Matrix(A'), Matrix(Cmat'), Σ_w, Matrix(V)))

# ── Averaged inner value (control H₂+cap  +  estimation H₂) ────────────────────
function inner_vg(θ; u_cap = Inf)
    Jd = 0.0; Je = 0.0; g = zeros(length(θ))
    for ω in ω_ens
        jd, gd, _, _, _, _ = ctrl_h2(θ, ω; u_cap = u_cap)
        je, ge = est_vg(θ, ω)
        Jd += jd; Je += je; g .+= gd .+ ge
    end
    (Jd/M_ens, Je/M_ens, g ./ M_ens)
end

# ── Shared power/mass budget (identical to example_adcs_hinf.jl) ───────────────
const θ_nom = [1.0, 1.0, 1.0]; const B_res = sum(θ_nom); const c_b = 5.0
J_des(θ)  = c_b * (sum(θ) - B_res)^2
∇J_des(θ) = fill(2c_b * (sum(θ) - B_res), 3)
const θ_lb = [0.3, 0.3, 0.3]; const θ_ub = [3.0, 5.0, 5.0]
const θ_names = ["e_rw", "α_st", "α_gyro"]
const θ_roles = ["f/B (wheel authority)", "V (star tracker)", "V (rate gyro)"]

function adcs_h2cap_report(; seed = 20240624, n_starts = 5, iterations = 100, cap_frac = 0.6)
    bar = "="^78
    println(bar)
    println("  ContEst (H₂ control SDP + effort cap  +  H₂ estimation SDP)  —  ADCS")
    println("  single linearization about origin, R = $Rw·I, $n states")
    println(bar)

    # Set the effort cap from the UNCONSTRAINED baseline effort so it binds.
    _, _, _, _, eff0, _ = ctrl_h2(θ_nom, ω_ens[1]; u_cap = Inf)
    u_cap = cap_frac * eff0
    @printf("\n  unconstrained baseline effort tr(Z₀) = %.4g  →  hard cap u_cap = %.4g (%.0f%%)\n",
            eff0, u_cap, 100*cap_frac)

    # ── Validity: unconstrained H₂ control SDP vs DARE; estimation SDP vs Kalman ──
    Jsdp0, _, _, _, _, _ = ctrl_h2(θ_nom, ω_ens[1]; u_cap = Inf)
    Jdare0 = ctrl_dare(θ_nom, ω_ens[1])
    A0, _ = model_lin(θ_nom, ω_ens[1]); V0 = Vmat(θ_nom)
    Jest_sdp, _, _, _ = est_sdp_dense(A0, V0); Jk = kalman_cost(A0, V0)
    println("\n── Inner-solver validity (baseline, unconstrained) ─────────")
    @printf("  H₂ control : SDP tr(QΣ)+tr(RZ₀) = %.5f   DARE tr(PW) = %.5f   gap %+.2f%%\n",
            Jsdp0, Jdare0, 100*(Jsdp0-Jdare0)/Jdare0)
    @printf("  H₂ estim.  : SDP tr(QΣ) = %.5f   Kalman = %.5f   gap %+.2f%%\n",
            Jest_sdp, Jk, 100*(Jest_sdp-Jk)/Jk)

    # ── Objective (cap active) + gradient check ──────────────────────────────
    f, g! = cached_objective(θ -> begin
        Jd, Je, gg = inner_vg(θ; u_cap = u_cap); (Jd + Je + J_des(θ), gg .+ ∇J_des(θ))
    end)
    ∇a, ∇fd, rel = verify_gradient(f, g!, θ_nom)
    println("\n── Gradient verification at θ_nom (cap active) ─────────────")
    @printf("  analytic : [%s]\n", join((@sprintf("%9.4f", v) for v in ∇a), " "))
    @printf("  fin.diff : [%s]\n", join((@sprintf("%9.4f", v) for v in ∇fd), " "))
    @printf("  relative error: %.2e  %s\n", rel, rel < 5e-2 ? "(OK)" : "(WARNING)")

    # ── Multi-start co-design under the cap ──────────────────────────────────
    println("\n── Multi-start BFGS ($n_starts starts) ──")
    best_J, best_θ, results = multistart_design(f, g!, θ_lb, θ_ub;
                                                n_starts = n_starts, θ_nom = θ_nom,
                                                seed = seed, g_tol = 1e-5, iterations = iterations)
    Js = sort([r.J for r in results]); nd = 1
    for k in 2:length(Js); Js[k]-Js[k-1] > 1e-4 && (nd += 1); end
    @printf("  starts run: %d   distinct minima: %d   best J_tot = %.5f\n", n_starts, nd, best_J)

    # ── Cost decomposition (baseline vs optimum, cap active) ─────────────────
    function bd(θ)
        jd, gd, Js, Jef, eff, _ = ctrl_h2(θ, ω_ens[1]; u_cap = u_cap)
        je, _ = est_vg(θ, ω_ens[1])
        (; J_est = je, J_state = Js, J_eff = Jef, effort = eff,
           J_cont = jd, J_c = jd + je, J_des = J_des(θ), J_tot = jd + je + J_des(θ))
    end
    base = bd(θ_nom); opt = bd(best_θ)
    pct(b, o) = abs(b) < 1e-12 ? "    — " : @sprintf("%+6.2f%%", 100*(o-b)/b)
    println("\n── Cost components: baseline θ_nom → optimal θ* ────────────")
    @printf("  %-30s %12s %12s %10s\n", "component", "baseline", "optimal", "Δ%")
    for (nm, b, o) in (("Estimation (H₂)   J_est", base.J_est,  opt.J_est),
                       ("Control state     tr(QΣ)", base.J_state, opt.J_state),
                       ("Control effort    tr(RZ₀)", base.J_eff, opt.J_eff),
                       ("Control total     J_cont", base.J_cont, opt.J_cont),
                       ("Inner  J_c=J_cont+J_est",  base.J_c,    opt.J_c),
                       ("Design cost       J_des",  base.J_des,  opt.J_des),
                       ("TOTAL   J_tot",            base.J_tot,  opt.J_tot))
        @printf("  %-30s %12.5f %12.5f  %s\n", nm, b, o, pct(b, o))
    end
    @printf("\n  effort tr(Z₀):  baseline %.4g  optimal %.4g   (cap = %.4g%s)\n",
            base.effort, opt.effort, u_cap,
            abs(opt.effort - u_cap) < 1e-3*u_cap ? ", ACTIVE" : "")
    @printf("  Net total reduction: %.5f  (%.2f%%)\n",
            base.J_tot-opt.J_tot, 100*(base.J_tot-opt.J_tot)/base.J_tot)

    println("\n── Design parameters: baseline → optimal ───────────────────")
    @printf("  %-10s %10s %10s %11s   %s\n", "param", "baseline", "optimal", "Δ", "enters/role")
    for i in 1:3
        @printf("  %-10s %10.3f %10.3f %+11.3f   %s\n",
                θ_names[i], θ_nom[i], best_θ[i], best_θ[i]-θ_nom[i], θ_roles[i])
    end
    @printf("  budget:  Σθ_nom = %.3f → Σθ* = %.3f  (B_res = %.1f)\n", sum(θ_nom), sum(best_θ), B_res)
    println(bar)
    return best_θ, best_J, base, opt, u_cap
end

if abspath(PROGRAM_FILE) == @__FILE__
    adcs_h2cap_report()
end
