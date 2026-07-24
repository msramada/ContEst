# example_adcs_hinfest.jl — EXPLORATORY ADCS SDP co-design that combines
#   • CONTROL    : H₂ state-feedback SDP WITH A HARD CONTROL-EFFORT CAP (from
#                  example_adcs_h2cap.jl), R = 0.1·I, cap tr(Z₀) ≤ u_cap on the
#                  steady-state input second-moment (control-effort covariance);
#   • ESTIMATION : H∞ (robust / worst-case) filter SDP — the transpose-dual of the
#                  H∞ bounded-real-lemma control SDP — replacing the H₂ estimation
#                  SDP of example_adcs_h2cap.jl / example_adcs_hinf.jl.
#
# THIS FILE CHANGES NOTHING ELSE (example_adcs.jl, example_adcs_hinf.jl,
# example_adcs_h2cap.jl, the paper, and src/ are all untouched).
#
# H∞ ESTIMATION recipe (as requested):
#   1. find the MINIMUM feasible γ² of the robust-filter game (bisection on the
#      filter GARE admissibility flag, `hinf_filter_gare`);
#   2. FIX γ² = 5 · γ²_min (a comfortable robustness margin above the boundary);
#   3. read the design gradient from the DUALS of the fixed-γ filter BRL SDP
#      (envelope contraction ⟨S, ∂M/∂θ⟩ with the LMI variables held at optimum).
#
# The filter BRL SDP is the exact transpose-dual of the control BRL SDP in
# example_adcs_hinf.jl under the swap  (A,B,√W,√Q,√R) → (Aᵀ,Cᵀ,√Q,√W,√V):
#   value  tr(T) = tr(Y⁻¹Q) = tr(Q Σe)   (worst-case guaranteed estimation cost),
# validated against tr(Q Σe) from the filter GARE.  The estimation-error weight is
# M = Q (Le = √Q).  On this study θ_h = (α_st,α_gyro) enters V only, so the whole
# estimation gradient is the ∂/∂V block of the BRL dual, obtained generically by
# ForwardDiff of the constraint matrix contracted with the dual.
#
# Run:  julia --project benchmarks/example_adcs_hinfest.jl

include("../src/LQR.jl")     # dare, dlyap  (Riccati cross-checks)
include("../src/Hinf.jl")    # hinf_filter_gare  (γ²_min bisection + GARE cross-check)
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
# (identical to example_adcs_h2cap.jl).  Returns (J, ∇θJ, J_state, J_effort,
#  effort=tr(Z₀), status).  u_cap=Inf ⇒ plain H₂.  Z₀ ⪰ K Σ Kᵀ is the stationary
#  input second-moment (control-effort covariance); tr(Z₀) ≤ u_cap caps it.
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

function ctrl_dare(θ, ω_op)
    A, B = model_lin(θ, ω_op)
    tr(dare(A, B, Qmat, Rmat) * Σ_w)
end

# ── ESTIMATION inner value: H∞ robust-filter BRL SDP (transpose-dual) ──────────
# Builds the (3n+ry)×(3n+ry) BRL LMI in (T, Y=Σe⁻¹, L=filter-gain·Y) for a FIXED
# attenuation γ².  Value tr(T) = tr(Q Σe) (worst-case guaranteed estimation cost).
# The block layout is example_adcs_hinf.jl's control BRL under the transpose-dual
# swap: A→Aᵀ, control-input B→measurement-injection Cᵀ, disturbance/trace channel
# √W→√Q, state-performance √Q→√W, input-performance √R→√V.
function est_hinf_brl_matrices(A, C, V, γ²; Yv = nothing, Lv = nothing)
    Ad = permutedims(A); Bd = permutedims(C)          # dual system (Aᵀ, Cᵀ)
    Qh = Matrix(sqrt(Symmetric(Qmat)))                # √Q : trace-obj + disturbance channel
    Wsh = Matrix(sqrt(Symmetric(Σ_w)))                # √W : state-performance weight
    # √V by elementwise diagonal sqrt (V is diagonal here): ForwardDiff-safe, whereas
    # sqrt(Symmetric(V)) is an eigendecomposition that is non-differentiable when V has
    # repeated eigenvalues (e.g. V(θ_nom) ∝ I) — that would make ∂J/∂V return NaN.
    Vh  = Matrix(Diagonal(sqrt.(diag(V))))            # √V : input-(injection-)performance weight
    nw = n; nz = n + ry
    # When Yv,Lv are supplied (envelope gradient), treat them as fixed data so the
    # matrix depends on θ only through Ad(θ) and Vh(θ).
    Y = Yv === nothing ? nothing : Yv
    L = Lv === nothing ? nothing : Lv
    return (Ad, Bd, Qh, Wsh, Vh, nw, nz, Y, L)
end

function est_hinf_brl_M(A, C, V, γ², Yv, Lv)
    Ad, Bd, Qh, Wsh, Vh, nw, nz, _, _ = est_hinf_brl_matrices(A, C, V, γ²)
    AYBL   = Ad * Yv .+ Bd * Lv
    CzYDzL = [Wsh * Yv; Vh * Lv]
    Mblk = [ -Yv                 zeros(eltype(Yv), n, nw)  permutedims(AYBL)   permutedims(CzYDzL);
             zeros(eltype(Yv), nw, n)  -γ² * Matrix(I, nw, nw)  permutedims(Qh)  zeros(eltype(Yv), nw, nz);
             AYBL                Qh                        -Yv                 zeros(eltype(Yv), n, nz);
             CzYDzL              zeros(eltype(Yv), nz, nw)  zeros(eltype(Yv), nz, n)  -Matrix(I, nz, nz) ]
    return Mblk
end

function est_hinf_sdp(A, C, V, γ²)
    Qh = Matrix(sqrt(Symmetric(Qmat)))
    Wsh = Matrix(sqrt(Symmetric(Σ_w))); Vh = Matrix(sqrt(Symmetric(Matrix(V))))
    nw = n; nz = n + ry
    Ad = permutedims(A); Bd = permutedims(C)
    mdl = Model(Clarabel.Optimizer); set_silent(mdl)
    @variable(mdl, T[1:n, 1:n], Symmetric)
    @variable(mdl, Y[1:n, 1:n], Symmetric)
    @variable(mdl, L[1:ry, 1:n])
    @constraint(mdl, Symmetric([T Qh; Qh Y]) in PSDCone())
    AYBL = Ad * Y .+ Bd * L
    CzYDzL = [Wsh * Y; Vh * L]
    Mblk = [ -Y            zeros(n, nw)            permutedims(AYBL)   permutedims(CzYDzL);
             zeros(nw, n)  -γ² * Matrix(I, nw, nw) permutedims(Qh)     zeros(nw, nz);
             AYBL          Qh                      -Y                  zeros(n, nz);
             CzYDzL        zeros(nz, nw)           zeros(nz, n)        -Matrix(I, nz, nz) ]
    @constraint(mdl, brl, Symmetric(-Mblk) in PSDCone())
    @objective(mdl, Min, tr(T))
    optimize!(mdl)
    st = termination_status(mdl)
    (st in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.NUMERICAL_ERROR)) ||
        return (1e8, nothing, nothing, nothing, st)
    return (objective_value(mdl), Matrix(dual(brl)), value.(Y), value.(L), st)
end

# Worst-case estimation cost from the filter GARE (cross-check for the SDP value).
function est_hinf_gare_cost(A, C, W, V, γ²)
    Le = Matrix(sqrt(Symmetric(Qmat)))
    Σe, _, _, ok = hinf_filter_gare(A, C, W, V, Le, γ²)
    ok ? (tr(Qmat * Σe), true) : (1e8, false)
end

# Minimum feasible γ² of the robust-filter game (bisection on GARE admissibility).
function est_hinf_gamma_min(A, C, W, V; lo = 1e-2, hi = 1e6, iters = 40)
    Le = Matrix(sqrt(Symmetric(Qmat)))
    feas(g²) = hinf_filter_gare(A, C, W, V, Le, g²)[4]
    feas(hi) || return hi                      # even huge γ² infeasible → give up high
    while !feas(lo); lo *= 2; lo > hi && return hi; end
    for _ in 1:iters
        mid = sqrt(lo * hi)                    # geometric bisection (γ² spans decades)
        feas(mid) ? (hi = mid) : (lo = mid)
    end
    return hi
end

# H∞ estimation value + design gradient at fixed γ², from the BRL dual (envelope).
function est_hinf_vg(θ, ω_op, γ²)
    A, _ = model_lin(θ, ω_op); V = Vmat(θ)
    J, S, Yv, Lv, st = est_hinf_sdp(A, Cmat, V, γ²)
    (J >= 1e8 || S === nothing) && return (1e8, fill(0.0, length(θ)))
    # Envelope: dJ/dθ_k = ⟨S, ∂(-Mblk)/∂θ_k⟩ with (T,Y,L) held at optimum.  The
    # constraint is -Mblk⪰0 with dual S; θ enters only through Aᵀ(θ) and √V(θ).
    # Match the validated control-BRL sign convention (∂J/∂A = 2 S₃₁Y): there
    # dJ/dθ = ⟨S, ∂Mblk/∂θ⟩, so we contract S with ∂Mblk/∂θ directly.
    g = ForwardDiff.gradient(
        t -> begin
            A_t, _ = model_lin(t, ω_op); V_t = Vmat(t)
            M = est_hinf_brl_M(A_t, Cmat, V_t, γ², Yv, Lv)
            sum(S .* M)
        end, θ)
    return J, g
end

kalman_cost(A, V) = tr(Qmat * dare(Matrix(A'), Matrix(Cmat'), Σ_w, Matrix(V)))

# ── Averaged inner value (control H₂+cap  +  estimation H∞) ────────────────────
function inner_vg(θ; u_cap = Inf, γ²_est = 1.0)
    Jd = 0.0; Je = 0.0; g = zeros(length(θ))
    for ω in ω_ens
        jd, gd, _, _, _, _ = ctrl_h2(θ, ω; u_cap = u_cap)
        je, ge = est_hinf_vg(θ, ω, γ²_est)
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

function adcs_hinfest_report(; seed = 20240624, n_starts = 5, iterations = 100,
                              cap_frac = 0.6, γ_margin = 5.0)
    bar = "="^80
    println(bar)
    println("  ContEst (H₂ control SDP + effort cap  +  H∞ robust-filter SDP)  —  ADCS")
    println("  single linearization about origin, R = $Rw·I, $n states")
    println(bar)

    A0, _ = model_lin(θ_nom, ω_ens[1]); V0 = Vmat(θ_nom)

    # ── (1) minimum γ²; (2) fix γ² = margin · γ²_min ──────────────────────────
    γ²_min = est_hinf_gamma_min(A0, Cmat, Σ_w, V0)
    γ²_est = γ_margin * γ²_min
    @printf("\n  minimum feasible γ² (robust filter) = %.4g   →   fixed γ² = %.1f·γ²_min = %.4g\n",
            γ²_min, γ_margin, γ²_est)

    # ── Effort cap from the UNCONSTRAINED baseline effort so it binds ──────────
    _, _, _, _, eff0, _ = ctrl_h2(θ_nom, ω_ens[1]; u_cap = Inf)
    u_cap = cap_frac * eff0
    @printf("  unconstrained baseline effort tr(Z₀) = %.4g  →  hard cap u_cap = %.4g (%.0f%%)\n",
            eff0, u_cap, 100*cap_frac)

    # ── Inner-solver validity ─────────────────────────────────────────────────
    Jsdp0, _, _, _, _, _ = ctrl_h2(θ_nom, ω_ens[1]; u_cap = Inf)
    Jdare0 = ctrl_dare(θ_nom, ω_ens[1])
    Jest_sdp, _, _, _, _ = est_hinf_sdp(A0, Cmat, V0, γ²_est)
    Jest_gare, okg = est_hinf_gare_cost(A0, Cmat, Σ_w, V0, γ²_est)
    Jk = kalman_cost(A0, V0)
    println("\n── Inner-solver validity (baseline) ────────────────────────")
    @printf("  H₂ control  : SDP tr(QΣ)+tr(RZ₀) = %.5f   DARE tr(PW) = %.5f   gap %+.2f%%\n",
            Jsdp0, Jdare0, 100*(Jsdp0-Jdare0)/Jdare0)
    @printf("  H∞ estim.   : BRL SDP tr(QΣe) = %.5f   filter GARE = %.5f   gap %+.2f%%\n",
            Jest_sdp, Jest_gare, okg ? 100*(Jest_sdp-Jest_gare)/Jest_gare : NaN)
    @printf("  H∞ vs H₂ estim.: H∞ tr(QΣe) = %.5f   Kalman (H₂) = %.5f   robust premium %+.2f%%\n",
            Jest_sdp, Jk, 100*(Jest_sdp-Jk)/Jk)

    # ── Objective (cap active, fixed γ²) + gradient check ─────────────────────
    f, g! = cached_objective(θ -> begin
        Jd, Je, gg = inner_vg(θ; u_cap = u_cap, γ²_est = γ²_est)
        (Jd + Je + J_des(θ), gg .+ ∇J_des(θ))
    end)
    ∇a, ∇fd, rel = verify_gradient(f, g!, θ_nom)
    println("\n── Gradient verification at θ_nom (cap active, fixed γ²) ────")
    @printf("  analytic : [%s]\n", join((@sprintf("%9.4f", v) for v in ∇a), " "))
    @printf("  fin.diff : [%s]\n", join((@sprintf("%9.4f", v) for v in ∇fd), " "))
    @printf("  relative error: %.2e  %s\n", rel, rel < 5e-2 ? "(OK)" : "(WARNING)")

    # ── Multi-start co-design ─────────────────────────────────────────────────
    println("\n── Multi-start BFGS ($n_starts starts) ──")
    best_J, best_θ, results = multistart_design(f, g!, θ_lb, θ_ub;
                                                n_starts = n_starts, θ_nom = θ_nom,
                                                seed = seed, g_tol = 1e-5, iterations = iterations)
    Js = sort([r.J for r in results]); nd = 1
    for k in 2:length(Js); Js[k]-Js[k-1] > 1e-4 && (nd += 1); end
    @printf("  starts run: %d   distinct minima: %d   best J_tot = %.5f\n", n_starts, nd, best_J)

    # ── Cost decomposition (baseline vs optimum, cap active) ──────────────────
    function bd(θ)
        jd, _, Js, Jef, eff, _ = ctrl_h2(θ, ω_ens[1]; u_cap = u_cap)
        je, _ = est_hinf_vg(θ, ω_ens[1], γ²_est)
        (; J_est = je, J_state = Js, J_eff = Jef, effort = eff,
           J_cont = jd, J_c = jd + je, J_des = J_des(θ), J_tot = jd + je + J_des(θ))
    end
    base = bd(θ_nom); opt = bd(best_θ)
    pct(b, o) = abs(b) < 1e-12 ? "    — " : @sprintf("%+6.2f%%", 100*(o-b)/b)
    println("\n── Cost components: baseline θ_nom → optimal θ* ────────────")
    @printf("  %-30s %12s %12s %10s\n", "component", "baseline", "optimal", "Δ%")
    for (nm, b, o) in (("Estimation (H∞)   J_est", base.J_est,  opt.J_est),
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
    return best_θ, best_J, base, opt, γ²_min, γ²_est, u_cap
end

if abspath(PROGRAM_FILE) == @__FILE__
    adcs_hinfest_report()
end
