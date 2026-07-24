# example_adcs_hinf.jl — EXPLORATORY reformulation of the spacecraft ADCS co-design
# as a LINEAR robust/stochastic SDP problem, instead of the nonlinear constrained
# eKF–MPC study in example_adcs.jl.  (Companion to example_pll_hinf.jl.)
#
# THIS FILE CHANGES NOTHING ELSE. It answers: "can the ADCS shared-budget co-design
# be posed on the ContEst SDP path (H∞ for control, H₂ for estimation) rather than
# the eKF–MPC path?"
#
# What is done differently from example_adcs.jl (deliberately, per the request):
#   • SKIP the input constraints (|τ_i| ≤ 0.8 dropped) — required to leave the MPC
#     path; a steady-state state-feedback design has no horizon to enforce them.
#   • LINEARISE the attitude dynamics about the pointing equilibrium ω=0.  The only
#     nonlinearity is the gyroscopic coupling ω×(Jω); its Jacobian VANISHES at ω=0,
#     so linearising there gives three decoupled double integrators.  (An ensemble of
#     nonzero body-rate points was tested and gave a byte-identical optimum: the
#     design θ enters B and V only, never the operating-point-dependent block of A,
#     so a single linearization about the origin suffices — ω_ens can hold more.)
#
# Design axes, SAME meaning / SAME shared budget as example_adcs.jl:
#   • θ = [e_rw (wheel authority, enters B) ; α_st, α_gyro (sensor precisions, in V)]
#   • all three draw on ONE budget  e_rw + α_st + α_gyro ≈ B_res  (soft equality) —
#     the competition is at the DESIGN-cost level, preserved exactly here.
#
# Inner solvers (both ContEst SDP-path values, gradients from the LMI duals):
#   • CONTROL    H∞ bounded-real-lemma (BRL) SDP in (Y=P⁻¹, L=KY), fixed γ²:
#                  J_det = min tr(Y⁻¹W) = tr(X W);  e_rw enters B;  the BRL dual gives
#                  the exact envelope gradients ∂J/∂A = 2 S₃₁ Y, ∂J/∂B = 2 S₃₁ Lᵀ
#                  (dense, r_x=6; equals the game-Riccati value tr(XW) at a full cert).
#   • ESTIMATION H₂ covariance SDP in the observability-gramian (P^ε) form (dense):
#                  J_est = tr(Q Σ);  α_st,α_gyro enter V(θ);  ∂J/∂V = Z*.
#
# Run:  julia --project benchmarks/example_adcs_hinf.jl

include("../src/LQR.jl")     # dare (kalman_cost reference)
include("../src/BFGS.jl")    # cached_objective, verify_gradient, multistart_design
using JuMP, Clarabel, LinearAlgebra, ForwardDiff, Printf, Random

# ── Spacecraft parameters (identical to example_adcs.jl) ──────────────────────
const Jx = 4.0; const Jy = 6.0; const Jz = 5.0
const Jmat = Matrix(Diagonal([Jx, Jy, Jz])); const Jinv = inv(Jmat)
const dt_s = 0.1

const n = 6      # states [ϕ,θ_a,ψ, ωx,ωy,ωz]
const m = 3      # inputs [τx,τy,τz]
const ry = 6     # outputs: full state (star tracker attitude + rate gyro)

# ── Noise / cost weights (identical to example_adcs.jl) ───────────────────────
const v_st0 = 0.20; const v_gyro0 = 0.20
const Σ_w  = Matrix(Diagonal([1e-6, 1e-6, 1e-6, 1e-4, 1e-4, 1e-4]))
const Qmat = Matrix(Diagonal([4.0, 4.0, 4.0, 1.0, 1.0, 1.0]))   # attitude 4× rate
# Control-effort weight. With the input constraints (|τ|≤0.8) dropped on the SDP
# path, R is the only thing bounding actuation; raise it to keep the H∞ controller
# from spending unlimited torque.  Non-const so a driver can sweep it.
Rmat = Matrix(Diagonal([1.0, 1.0, 1.0]))
const Cmat = Matrix(1.0I, ry, n)                 # full-state readout (C = I)
Vmat(θ) = Matrix(Diagonal([v_st0/θ[2]^2, v_st0/θ[2]^2, v_st0/θ[2]^2,
                           v_gyro0/θ[3]^2, v_gyro0/θ[3]^2, v_gyro0/θ[3]^2]))

const γ² = 16.0         # H∞ attenuation level (admissible at baseline; matches the
                        # paper's convention).  The ADCS plant is benign (no sharp
                        # resonant mode), so the worst-case premium over H₂ is small;
                        # as γ→∞ the control cost collapses to the H₂/LQG value.

# Jacobian of the gyroscopic term  gyro(ω)=ω×(Jω)  at ω_op (couples the rate axes).
#   gyro = [(Jz−Jy)ωyωz, (Jx−Jz)ωzωx, (Jy−Jx)ωxωy]
function gyro_jac(ω)
    ωx, ωy, ωz = ω
    [0.0            (Jz-Jy)*ωz     (Jz-Jy)*ωy;
     (Jx-Jz)*ωz     0.0            (Jx-Jz)*ωx;
     (Jy-Jx)*ωy     (Jy-Jx)*ωx     0.0]
end

# ── Linearised discrete model about a body-rate operating point ω_op ──────────
# A depends on ω_op (gyroscopic coupling), NOT on θ; e_rw enters B only.
function model_lin(θ, ω_op)
    T = eltype(θ)
    e_rw = θ[1]
    A = Matrix{T}(I, n, n)
    A[1,4] = dt_s; A[2,5] = dt_s; A[3,6] = dt_s            # attitude ← rate
    A[4:6, 4:6] .= I(3) .- dt_s .* (Jinv * gyro_jac(ω_op)) # rate coupling at ω_op
    B = zeros(T, n, m)
    B[4:6, 1:3] .= dt_s .* e_rw .* Jinv                    # wheel authority → rates
    return A, B
end

# ── Linearization point: the pointing equilibrium (origin) ────────────────────
# A single linearization about ω=0.  At the origin the gyroscopic Jacobian
# ∂(ω×Jω)/∂ω vanishes, so A is three decoupled double integrators.  (An ensemble
# of nonzero body-rate points was tested and gave a byte-identical optimum: the
# design θ enters B and V only, never the operating-point-dependent block of A, so
# averaging over operating points does not move the optimum — one point suffices.)
const ω_ens = [[0.0, 0.0, 0.0]]
const M_ens = length(ω_ens)

# ── CONTROL inner value (H∞ bounded-real-lemma SDP) at one operating point ────
#   J = min tr(T)  s.t.  [T  W^{1/2}; W^{1/2}  Y] ⪰ 0  and  the BRL LMI
#       (in Y=P⁻¹, L=KY at fixed γ²); at optimum J = tr(Y⁻¹W) = tr(XW) (game value).
# Envelope gradient from the BRL dual S: ∂J/∂A = 2 S₃₁ Y, ∂J/∂B = 2 S₃₁ Lᵀ, chained
# through A(θ),B(θ).  Dense (r_x=6) — no ring/block structure needed at this scale.
function ctrl_vg(θ, ω_op)
    A, B = model_lin(θ, ω_op)
    Wh = Matrix(sqrt(Symmetric(Σ_w)))
    Qh = Matrix(sqrt(Symmetric(Qmat))); Rh = Matrix(sqrt(Symmetric(Rmat)))
    nw = n; nz = n + m; μ = γ²
    mdl = Model(Clarabel.Optimizer); set_silent(mdl)
    @variable(mdl, T[1:n, 1:n], Symmetric)
    @variable(mdl, Y[1:n, 1:n], Symmetric)
    @variable(mdl, L[1:m, 1:n])
    @constraint(mdl, Symmetric([T Wh; Wh Y]) in PSDCone())
    AYBL = A * Y .+ B * L
    CzYDzL = [Qh * Y; Rh * L]
    Mblk = [ -Y            zeros(n, nw)            permutedims(AYBL)   permutedims(CzYDzL);
             zeros(nw, n)  -μ * Matrix(I, nw, nw)  permutedims(Wh)     zeros(nw, nz);
             AYBL          Wh                      -Y                  zeros(n, nz);
             CzYDzL        zeros(nz, nw)           zeros(nz, n)        -Matrix(I, nz, nz) ]
    @constraint(mdl, brl, Symmetric(-Mblk) in PSDCone())
    @objective(mdl, Min, tr(T))
    optimize!(mdl)
    st = termination_status(mdl)
    (st in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.NUMERICAL_ERROR)) ||
        return (1e8, fill(0.0, length(θ)))
    J = objective_value(mdl)
    S = dual(brl); b3 = (n + nw + 1):(2n + nw)
    S31 = Matrix(S[b3, 1:n])
    Yv = value.(Y); Lv = value.(L)
    gA = 2 .* (S31 * Yv)                    # ∂J/∂A   (BRL dual, plant block)
    gB = 2 .* (S31 * permutedims(Lv))       # ∂J/∂B
    Aj = ForwardDiff.jacobian(t -> vec(model_lin(t, ω_op)[1]), θ)
    Bj = ForwardDiff.jacobian(t -> vec(model_lin(t, ω_op)[2]), θ)
    ∇ = [tr(gA' * reshape(Aj[:, k], n, n)) + tr(gB' * reshape(Bj[:, k], n, m))
         for k in eachindex(θ)]
    return J, ∇
end

# ── ESTIMATION inner value (H₂) at one operating point: dense P^ε covariance SDP ─
function est_sdp_dense(A, V; withgrad = true)
    mdl = Model(Clarabel.Optimizer); set_silent(mdl)
    @variable(mdl, P[1:n, 1:n], Symmetric)
    @variable(mdl, F[1:n, 1:ry])
    @variable(mdl, Z[1:ry, 1:ry], Symmetric)
    @constraint(mdl, Symmetric([Z F'; F P]) in PSDCone())
    lmi2 = @constraint(mdl, Symmetric([P .- Qmat  (A'*P .- Cmat'*F');
                                       (P*A .- F*Cmat)  P]) in PSDCone())
    @objective(mdl, Min, tr(Σ_w * P) + tr(V * Z))
    optimize!(mdl)
    st = termination_status(mdl)
    (st in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.NUMERICAL_ERROR)) ||
        return (1e8, zeros(n, n), zeros(ry, ry), st)
    J = objective_value(mdl)
    withgrad || return (J, zeros(n, n), zeros(ry, ry), st)
    Zv = value.(Z)
    S = dual(lmi2); S12 = S[1:n, n+1:2n]
    gA = -2 .* (value.(P) * permutedims(S12))
    gV = Matrix(Zv)
    (J, gA, gV, st)
end

function est_vg(θ, ω_op)
    A, _ = model_lin(θ, ω_op); V = Vmat(θ)
    J, gA, gV, st = est_sdp_dense(A, V)
    J >= 1e8 && return (1e8, fill(0.0, length(θ)))
    # A is θ-independent ⇒ only V(θ) contributes; keep the A-term for generality.
    Aj = ForwardDiff.jacobian(t -> vec(model_lin(t, ω_op)[1]), θ)
    Vj = ForwardDiff.jacobian(t -> vec(Matrix(Vmat(t))), θ)
    ∇ = [tr(gA' * reshape(Aj[:, k], n, n)) + tr(gV' * reshape(Vj[:, k], ry, ry))
         for k in eachindex(θ)]
    return J, ∇
end

kalman_cost(A, V) = tr(Qmat * dare(Matrix(A'), Matrix(Cmat'), Σ_w, Matrix(V)))

function inner_vg(θ)
    Jd = 0.0; Je = 0.0; g = zeros(length(θ))
    for ω_op in ω_ens
        jd, gd = ctrl_vg(θ, ω_op)
        je, ge = est_vg(θ, ω_op)
        Jd += jd; Je += je; g .+= gd .+ ge
    end
    return (Jd/M_ens, Je/M_ens, g ./ M_ens)
end

# ── Design cost: shared power/mass budget (identical to example_adcs.jl) ───────
const θ_nom = [1.0, 1.0, 1.0]
const B_res = sum(θ_nom); const c_b = 5.0
J_des(θ)  = c_b * (sum(θ) - B_res)^2
∇J_des(θ) = fill(2c_b * (sum(θ) - B_res), 3)

const θ_lb = [0.3, 0.3, 0.3]; const θ_ub = [3.0, 5.0, 5.0]
const θ_names = ["e_rw", "α_st", "α_gyro"]
const θ_roles = ["f/B (wheel authority)", "V (star tracker)", "V (rate gyro)"]

function eval_J(θ)
    Jd, Je, g = inner_vg(θ)
    (Jd + Je + J_des(θ), g .+ ∇J_des(θ))
end
const contest_f, contest_g! = cached_objective(eval_J)

function adcs_hinf_report(; seed = 20240624, n_starts = 5, iterations = 100, R_w = 1.0)
    global Rmat = Matrix(Diagonal(fill(float(R_w), m)))   # control-effort weight
    bar = "="^76
    println(bar)
    println("  ContEst (H∞ control + H₂ estimation SDP)  —  spacecraft ADCS co-design")
    println("  linearised about $M_ens operating point(s), γ²=$γ², $n states")
    println("  control-effort weight  R = $(R_w)·I")
    println(bar)

    ∇a, ∇fd, rel = verify_gradient(contest_f, contest_g!, θ_nom)
    println("\n── Gradient verification at θ_nom ──────────────────────────")
    @printf("  analytic : [%s]\n", join((@sprintf("%9.4f", v) for v in ∇a), " "))
    @printf("  fin.diff : [%s]\n", join((@sprintf("%9.4f", v) for v in ∇fd), " "))
    @printf("  relative error: %.2e  %s\n", rel, rel < 5e-2 ? "(OK)" : "(WARNING)")

    A0, _ = model_lin(θ_nom, ω_ens[1]); V0 = Vmat(θ_nom)
    Jsdp, _, _, _ = est_sdp_dense(A0, V0); Jk = kalman_cost(A0, V0)
    println("\n── H₂ estimation SDP validity (baseline, equilibrium op point) ──")
    @printf("  dense P^ε SDP tr(QΣ) = %.4f   centralised Kalman = %.4f   gap %+.2f%%\n",
            Jsdp, Jk, 100*(Jsdp-Jk)/Jk)

    println("\n── Multi-start BFGS ($n_starts starts: θ_nom + $(n_starts-1) random) ──")
    best_J, best_θ, results = multistart_design(contest_f, contest_g!, θ_lb, θ_ub;
                                                n_starts = n_starts, θ_nom = θ_nom,
                                                seed = seed, g_tol = 1e-5, iterations = iterations)
    Js = sort([r.J for r in results]); nd = 1
    for k in 2:length(Js); Js[k]-Js[k-1] > 1e-4 && (nd += 1); end
    @printf("  starts run: %d   distinct minima: %d   best J_tot = %.4f\n", n_starts, nd, best_J)
    println(nd == 1 ? "  → all starts agree: optimum is (numerically) global." :
                      "  → multimodal: reporting the best minimum only.")

    function bd(θ)
        Jd, Je, _ = inner_vg(θ)
        (; J_est = Je, J_det = Jd, J_c = Jd+Je, J_des = J_des(θ), J_tot = Jd+Je+J_des(θ))
    end
    base = bd(θ_nom); opt = bd(best_θ)
    pct(b, o) = abs(b) < 1e-9 ? "    — " : @sprintf("%+6.2f%%", 100*(o-b)/b)
    println("\n── Cost components: baseline θ_nom → optimal θ* (ensemble mean) ──")
    @printf("  %-28s %12s %12s %12s %8s\n", "component", "baseline", "optimal", "Δ (incr +)", "%")
    for (nm, b, o) in (("Estimation (H₂)  J_est", base.J_est, opt.J_est),
                       ("Control    (H∞)  J_det", base.J_det, opt.J_det),
                       ("Inner  J_c=J_det+J_est", base.J_c,   opt.J_c),
                       ("Design cost      J_des", base.J_des, opt.J_des),
                       ("TOTAL   J_tot=J_c+J_des", base.J_tot, opt.J_tot))
        @printf("  %-28s %12.4f %12.4f %12.4f  %s\n", nm, b, o, o-b, pct(b, o))
    end
    @printf("\n  Net total reduction: %.4f  (%.2f%%)\n",
            base.J_tot-opt.J_tot, 100*(base.J_tot-opt.J_tot)/base.J_tot)

    println("\n── Design parameters: baseline → optimal (best minimum) ────")
    @printf("  %-10s %10s %10s %11s   %s\n", "param", "baseline", "optimal", "Δ", "enters/role")
    for i in 1:length(θ_nom)
        @printf("  %-10s %10.3f %10.3f %+11.3f   %s\n",
                θ_names[i], θ_nom[i], best_θ[i], best_θ[i]-θ_nom[i], θ_roles[i])
    end
    @printf("  budget:  Σθ_nom = %.3f → Σθ* = %.3f  (B_res = %.1f)\n",
            sum(θ_nom), sum(best_θ), B_res)
    println(bar)
    return best_θ, best_J, base, opt
end

if abspath(PROGRAM_FILE) == @__FILE__
    adcs_hinf_report()
end
