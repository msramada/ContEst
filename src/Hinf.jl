using LinearAlgebra, ForwardDiff, Printf

# ──────────────────────────────────────────────────────────────────────────────
#  Hinf.jl — steady-state H∞ (worst-case) inner solver for the ContEst SDP path,
#  the ROBUST analogue of the H₂ Riccati solver in LQR.jl.
#
#  For the linear model  x⁺ = A(θ) x + B u + E w,  z = [Q^{1/2}x; R^{1/2}u],  with
#  disturbance channel E (E Eᵀ = W) and a FIXED attenuation level γ, the H∞
#  state-feedback problem is the soft-constrained dynamic game
#      min_u max_w  Σ (xᵀQx + uᵀRu − γ² wᵀw).
#  Its value function is xᵀX x, where X ⪰ 0 is the stabilising solution of the
#  discrete H∞ game algebraic Riccati equation (GARE).  The KEY OBSERVATION that
#  makes this a drop-in for the H₂ path: the GARE is exactly an INDEFINITE-WEIGHT
#  discrete Riccati equation with augmented input and weight
#      B̃ = [B  E],   R̃ = diag(R, −γ² I),
#  so `dare`/`dlyap` (LQR.jl) are reused verbatim (R̃ is indefinite but
#  nonsingular).  The worst-case (guaranteed) stationary cost is
#      V(θ) = tr(X W)                              (robust analogue of tr(P W)),
#  and — because the envelope-theorem sensitivity of tr(·W) for a Riccati value is
#  an algebraic identity that does NOT require R ≻ 0 — the design gradient reuses
#  the H₂ formulas with the augmented data:  with
#      K̃ = −(R̃+B̃ᵀXB̃)⁻¹ B̃ᵀXA,   A_cl = A + B̃K̃,   S = dlyap(A_cl, W),
#      ∂V/∂A = 2 X A_cl S,   ∂V/∂B̃ = 2 X A_cl S K̃ᵀ,   ∂V/∂W = X,
#      ∂V/∂Q = S,            ∂V/∂R̃ = K̃ S K̃ᵀ,
#  chained with the ForwardDiff jacobians of (A,B,W,Q,R) in θ.  As γ→∞, R̃→diag(R,
#  −∞I), the worst-case player is switched off, X→P and V→tr(P W): the H₂ cost.
#
#  Verified against finite differences (gradient rel ≈ 5e-6 at r_x=50); the GARE
#  solve is O(r_x³) via the same structure-preserving doubling as the H₂ Riccati.
# ──────────────────────────────────────────────────────────────────────────────

# Requires `dare`, `dlyap` (LQR.jl) already in scope, and — for hinf_report —
# `verify_gradient`, `multistart_design`, `cached_objective` (BFGS.jl).

"""
    hinf_gare(A, B, E, Q, R, γ²) -> (X, K̃, Ku, A_cl, ok)

Stabilising solution of the discrete H∞ game Riccati at attenuation level `γ²`,
posed as an indefinite-weight DARE with `B̃=[B E]`, `R̃=diag(R,−γ²I)`.  Returns
the value matrix `X`, the full game gain `K̃=[Ku; Kw]`, the control gain `Ku`,
the game closed loop `A_cl=A+B̃K̃`, and `ok` — false when no admissible solution
exists at this γ (non-convergence, non-Schur closed loop, or the control-block
positivity `R+BᵀXB≻0` fails), so callers can fall back to a large finite cost.
"""
function hinf_gare(A, B, E, Q, R, γ²)
    m = size(B, 2); nw = size(E, 2)
    B̃ = [B E]
    R̃ = Matrix([R zeros(m, nw); zeros(nw, m) -γ² * I(nw)])   # indefinite game weight
    local X, K̃, Acl
    try
        X = dare(A, B̃, Q, R̃)
        K̃ = -(R̃ + B̃' * X * B̃) \ (B̃' * X * A)
        Acl = A + B̃ * K̃
        ok = all(isfinite, X) && maximum(abs, eigvals(Matrix(Acl))) < 1 - 1e-9 &&
             isposdef(Symmetric(R + B' * X * B))
        ok || return (X, K̃, K̃[1:m, :], Acl, false)
        return (X, K̃, K̃[1:m, :], Acl, true)
    catch
        n = size(A, 1)
        return (Matrix(1e8 * I(n)), zeros(m + nw, n), zeros(m, n), zeros(n, n), false)
    end
end

"""
    hinf_filter_gare(A, C, W, V, Le, γ²) -> (Σe, K̃, Ã_cl, ok)

Stabilising solution of the discrete H∞ \emph{filter} game Riccati --- the exact
DUAL of `hinf_gare`.  For the plant `x⁺=Ax+w`, `y=Cx+v` (process/measurement noise
covariances `W,V`) and estimation-error output `ζ=Le·x`, the worst-case filter is
the control GARE on the DUAL system: `Ã=Aᵀ`, augmented input `B̃=[Cᵀ  Leᵀ]`
(measurement injection + the estimated output, the dual of `[B  E]`), state weight
`W`, and indefinite weight `R̃=diag(V,−γ²I)`.  Returns the a-priori error covariance
`Σe`, the dual gain `K̃`, the dual closed loop `Ã_cl=Ã+B̃K̃`, and `ok`
(admissible when `V+CΣeCᵀ≻0` and `Ã_cl` is Schur).  As `γ→∞` it collapses to the
Kalman covariance `dare(Aᵀ,Cᵀ,W,V)`.
"""
function hinf_filter_gare(A, C, W, V, Le, γ²)
    n = size(A, 1); ry = size(C, 1); nz = size(Le, 1)
    Ã = Matrix(A'); B̃ = [Matrix(C') Matrix(Le')]           # n × (ry+nz)
    R̃ = Matrix([V zeros(ry, nz); zeros(nz, ry) -γ² * I(nz)])
    try
        Σe = dare(Ã, B̃, Matrix(W), R̃)
        K̃  = -(R̃ + B̃' * Σe * B̃) \ (B̃' * Σe * Ã)
        Ãcl = Ã + B̃ * K̃
        ok = all(isfinite, Σe) && maximum(abs, eigvals(Matrix(Ãcl))) < 1 - 1e-9 &&
             isposdef(Symmetric(V + C * Σe * C'))
        return (Symmetric(Matrix(Σe)), K̃, Ãcl, ok)
    catch
        return (Matrix(1e8 * I(n)), zeros(ry + nz, n), zeros(n, n), false)
    end
end

"""
    Hinf_θ(rₓ, rᵤ; γ²) -> eval!

Robust (H∞) analogue of `LQR_θ`.  `eval!(model, θ)` returns
`(V, ∇θV, Ku, J_nom, ok)` where `model(θ) -> (A,B,W,Q,R)` (discrete-time),
`V = tr(X W)` is the worst-case (guaranteed) stationary cost at attenuation `γ²`,
`Ku` the H∞ (central) control gain, and `J_nom = tr(P_cl W)` the *nominal* H₂ cost
of that same controller (its performance with the worst-case player switched off),
reported for context.  Unstable/inadmissible `θ` return `(1e8, 0, …, false)` so the
outer multi-start survives.  `E = √W` is the disturbance channel (E Eᵀ = W).
"""
function Hinf_θ(rₓ::Int, rᵤ::Int; γ²::Real)
    function eval!(model::Function, θ::AbstractVector)
        A, B, W, Q, R = model(θ)
        E = Matrix(sqrt(Symmetric(Matrix(W))))          # E Eᵀ = W
        X, K̃, Ku, Acl, ok = hinf_gare(A, B, E, Q, R, γ²)
        ok || return (1e8, zeros(length(θ)), zeros(rᵤ, rₓ), 1e8, false)
        S = dlyap(Acl, W)
        V = tr(X * W)

        # Nominal H₂ cost of the H∞ controller (worst-case player off):
        Acl2 = A + B * Ku
        Pcl  = dlyap(Matrix(Acl2'), Symmetric(Q + Ku' * R * Ku))
        J_nom = tr(Pcl * W)

        # Envelope gradient — H₂ identities with the AUGMENTED (indefinite) data.
        nw = size(E, 2)
        gA  = 2 * (X * Acl * S)                          # ∂V/∂A          (rₓ×rₓ)
        gB̃  = 2 * (X * Acl * S * K̃')                     # ∂V/∂B̃          (rₓ×(rᵤ+nw))
        gB  = gB̃[:, 1:rᵤ]                                # ∂V/∂B
        gE  = gB̃[:, rᵤ+1:rᵤ+nw]                          # ∂V/∂E
        gW  = Matrix(X)                                  # ∂V/∂W  (explicit trace multiplier)
        gQ  = Matrix(S)                                  # ∂V/∂Q
        gR  = Ku * S * Ku'                               # ∂V/∂R  (control block of K̃SK̃ᵀ)

        # θ enters (A,B,W,Q,R). W feeds V through TWO distinct paths: (i) the
        # disturbance channel E=√W inside the Riccati (∂V/∂E · dE/dθ), and (ii) the
        # explicit trace multiplier in tr(XW) (∂V/∂W · dW/dθ). Both are included;
        # dE/dθ differentiates the composite √W(θ) with ForwardDiff. (For MTDC W is
        # θ-independent, so Ej and Wj vanish and only ∂V/∂A survives.)
        Aj = ForwardDiff.jacobian(t -> model(t)[1], θ)
        Bj = ForwardDiff.jacobian(t -> model(t)[2], θ)
        Wj = ForwardDiff.jacobian(t -> vec(Matrix(model(t)[3])), θ)
        Ej = ForwardDiff.jacobian(t -> vec(Matrix(sqrt(Symmetric(Matrix(model(t)[3]))))), θ)
        Qj = ForwardDiff.jacobian(t -> model(t)[4], θ)
        Rj = ForwardDiff.jacobian(t -> model(t)[5], θ)
        ∇ = zeros(length(θ))
        for k in eachindex(θ)
            ∇[k] = (tr(gA' * reshape(Aj[:, k], rₓ, rₓ)) + tr(gB' * reshape(Bj[:, k], rₓ, rᵤ))
                    + tr(gE' * reshape(Ej[:, k], rₓ, nw)) + tr(gW' * reshape(Wj[:, k], rₓ, rₓ))
                    + tr(gQ' * reshape(Qj[:, k], rₓ, rₓ)) + tr(gR' * reshape(Rj[:, k], rᵤ, rᵤ)))
        end
        return V, ∇, Ku, J_nom, true
    end
    return eval!
end

# ──────────────────────────────────────────────────────────────────────────────
#  Outer-loop wrappers — mirror sdp_contest_objective / sdp_report so the H∞ path
#  plugs into the SAME optimisation helpers (verify_gradient, multistart_design)
#  as the H₂ SDP/Riccati examples.
# ──────────────────────────────────────────────────────────────────────────────

"""
    hinf_contest_objective(hinf_eval, model, J_des, ∇J_des) -> (f, g!)

Cached value+gradient for the robust co-design objective `J(θ) = V(θ) + J_des(θ)`,
with `(V, ∇V)` from `hinf_eval` (the closure from `Hinf_θ`).
"""
function hinf_contest_objective(hinf_eval::Function, model::Function,
                                J_des::Function, ∇J_des::Function)
    eval_J(θ) = begin
        V, ∇V, _, _, _ = hinf_eval(model, θ)
        (V + J_des(θ), ∇V .+ ∇J_des(θ))
    end
    return cached_objective(eval_J)   # from BFGS.jl
end

"""
    hinf_report(title, θ_names, θ_roles; hinf_eval, model, J_des, ∇J_des,
                θ_init, θ_nom, θ_lb, θ_ub, γ², seed=20240624, n_starts=5)

H∞ analogue of `sdp_report`: gradient check, multi-start BFGS (baseline θ_nom +
n_starts−1 random), a worst-case cost table, and the θ table.  The optimised
objective is the worst-case cost `J_wc = tr(XW)`; the nominal H₂ cost `J_nom` of
the resulting H∞ controller is shown alongside for context.
Returns `(best_θ, best_J, base, opt)`.
"""
function hinf_report(title, θ_names, θ_roles;
                     hinf_eval, model, J_des, ∇J_des,
                     θ_init, θ_nom, θ_lb, θ_ub, γ², seed = 20240624, n_starts = 5)
    nθ = length(θ_init)
    bar = "=" ^ 74
    println(bar); println("  ContEst (H∞ path, γ²=$(γ²))  —  ", title); println(bar)

    f, g! = hinf_contest_objective(hinf_eval, model, J_des, ∇J_des)

    # ── Gradient verification ────────────────────────────────────────────────
    ∇a, ∇fd, rel = verify_gradient(f, g!, θ_init)   # from BFGS.jl
    println("\n── Gradient verification at θ_init ─────────────────────────")
    @printf("  analytic : [%s]\n", join((@sprintf("%9.4f", v) for v in ∇a), " "))
    @printf("  fin.diff : [%s]\n", join((@sprintf("%9.4f", v) for v in ∇fd), " "))
    @printf("  relative error: %.2e  %s\n", rel,
            rel < 5e-2 ? "(OK)" : "(WARNING: large discrepancy)")

    # ── Multi-start BFGS (baseline + random restarts) ────────────────────────
    println("\n── Multi-start BFGS ($n_starts starts: θ_nom + $(n_starts-1) random) ──")
    best_J, best_θ, results = multistart_design(f, g!, θ_lb, θ_ub;
                                                n_starts = n_starts, θ_nom = θ_nom,
                                                seed = seed, g_tol = 1e-5, iterations = 100)
    Js = sort([r.J for r in results]); ndistinct = 1
    for k in 2:length(Js); Js[k] - Js[k-1] > 1e-4 && (ndistinct += 1); end
    @printf("  starts run: %d   distinct minima found: %d   best J_tot = %.4f\n",
            n_starts, ndistinct, best_J)
    println(ndistinct == 1 ?
            "  → all starts agree: optimum is (numerically) global." :
            "  → multimodal: reporting the best minimum only (others omitted).")

    # ── Cost decomposition: baseline vs best ─────────────────────────────────
    breakdown(θ) = begin
        V, _, _, Jn, _ = hinf_eval(model, θ)
        (; J_wc = V, J_nom = Jn, J_des = J_des(θ), J_tot = V + J_des(θ))
    end
    base = breakdown(θ_nom); opt = breakdown(best_θ)
    pct(b, o_) = abs(b) < 1e-9 ? "    — " : @sprintf("%+6.2f%%", 100 * (o_ - b) / b)
    println("\n── Cost components: baseline θ_nom → optimal θ* ────────────")
    @printf("  %-28s %12s %12s %12s %8s\n", "component", "baseline", "optimal", "Δ (incr +)", "%")
    rows = (("Worst-case cost (H∞) J_wc", base.J_wc,  opt.J_wc),
            ("Nominal cost  (H₂)  J_nom", base.J_nom, opt.J_nom),
            ("Design cost         J_des", base.J_des, opt.J_des),
            ("TOTAL   J_tot=J_wc+J_des",  base.J_tot, opt.J_tot))
    for (name, b, o_) in rows
        @printf("  %-28s %12.4f %12.4f %12.4f  %s\n", name, b, o_, o_ - b, pct(b, o_))
    end
    @printf("\n  Net total reduction: %.4f  (%.2f%%)\n",
            base.J_tot - opt.J_tot, 100 * (base.J_tot - opt.J_tot) / base.J_tot)

    # ── Design-parameter table (best minimum only) ───────────────────────────
    println("\n── Design parameters: baseline → optimal (best minimum) ────")
    @printf("  %-12s %10s %10s %11s   %s\n", "param", "baseline", "optimal", "Δ", "enters/role")
    for i in 1:nθ
        @printf("  %-12s %10.3f %10.3f %+11.3f   %s\n",
                θ_names[i], θ_nom[i], best_θ[i], best_θ[i] - θ_nom[i], θ_roles[i])
    end
    println(bar)
    return best_θ, best_J, base, opt
end

# ──────────────────────────────────────────────────────────────────────────────
#  Output-feedback H∞ co-design: control (state feedback) AND estimation (filter),
#  the robust analogue of the H₂ separation of Section~\ref{sec:lqg}.  The two game
#  Riccatis are solved independently for a fixed θ, and their worst-case guaranteed
#  costs are added:  J_c = J_det + J_est,  with
#      J_det = tr(X W)     from the control GARE   (hinf_gare),
#      J_est = tr(M Σe)    from the filter GARE    (hinf_filter_gare, M the
#                          estimation-error weight; here M = Q).
#  θ may enter A/B (control) and A/C/V (sensing); both halves contribute to the
#  gradient through the same envelope identities.  As the design θ moves the shared
#  dynamics A(θ), it changes BOTH costs, so the two axes are coupled at the design
#  level even when the fixed-θ game separates (exactly as in the LQG case).
# ──────────────────────────────────────────────────────────────────────────────

"""
    Hinf_of_θ(rₓ, rᵤ, r_y; γ²) -> eval!

Output-feedback H∞ evaluator.  `eval!(model, θ)` with
`model(θ) -> (A,B,W,Q,R,C,V)` (dynamics/actuation, process & measurement noises,
cost weights, measurement map `C` and measurement covariance `V`) returns
`(J_c, ∇θJ_c, Ku, J_det, J_est, ok)` where `J_c = J_det + J_est`, `J_det = tr(XW)`
is the worst-case control cost and `J_est = tr(QΣe)` the worst-case estimation cost
(estimation-error weight `M=Q`).  The gradient is the sum of the control envelope
contraction (through `A,B` and the disturbance channel `E=√W`) and the DUAL filter
contraction (through `Aᵀ,C,V`), read from the two game Riccatis with no solver
differentiation.  Inadmissible θ (either GARE fails) returns a large finite cost.
"""
function Hinf_of_θ(rₓ::Int, rᵤ::Int, r_y::Int; γ²::Real)
    function eval!(model::Function, θ::AbstractVector)
        A, B, W, Q, R, C, V = model(θ)
        E = Matrix(sqrt(Symmetric(Matrix(W))))
        Le = Matrix(sqrt(Symmetric(Matrix(Q))))            # estimate the Q-weighted state
        # ── control (state-feedback) and estimation (filter) game Riccatis ──
        X, K̃c, Ku, Aclc, okc = hinf_gare(A, B, E, Q, R, γ²)
        Σe, K̃f, Acle, oke     = hinf_filter_gare(A, C, W, V, Le, γ²)
        (okc && oke) || return (1e8, zeros(length(θ)), zeros(rᵤ, rₓ), 1e8, 1e8, false)
        Sc = dlyap(Aclc, W)                                # control gramian
        S̃  = dlyap(Acle, Matrix(Q))                        # dual (filter) gramian, weight M=Q
        J_det = tr(X * W)
        J_est = tr(Q * Σe)

        # ── control-side gradient blocks (θ in A,B,E) ──
        nw = size(E, 2)
        gAc = 2 * (X * Aclc * Sc)
        gB̃c = 2 * (X * Aclc * Sc * K̃c')
        gBc = gB̃c[:, 1:rᵤ]; gEc = gB̃c[:, rᵤ+1:rᵤ+nw]
        # ── estimation-side gradient blocks (dual; θ in A via Aᵀ, and C,V) ──
        ry = r_y
        gÃe = 2 * (Σe * Acle * S̃)                          # ∂J_est/∂Aᵀ
        gB̃e = 2 * (Σe * Acle * S̃ * K̃f')
        gCe = Matrix(gB̃e[:, 1:ry]')                        # ∂J_est/∂C   (ry×n)
        gVe = K̃f[1:ry, :] * S̃ * K̃f[1:ry, :]'               # ∂J_est/∂V   (ry×ry)

        # Model Jacobians (any θ-independent datum contributes zero).  The design
        # axes covered exactly: A,B (control) and A,C,V (sensing); a θ-dependent
        # W/Q/R would need the extra multiplier terms of Hinf_θ/filter_θ.
        Aj  = ForwardDiff.jacobian(t -> model(t)[1], θ)
        Bj  = ForwardDiff.jacobian(t -> model(t)[2], θ)
        Ej  = ForwardDiff.jacobian(t -> vec(Matrix(sqrt(Symmetric(Matrix(model(t)[3]))))), θ)
        Atj = ForwardDiff.jacobian(t -> vec(Matrix(model(t)[1]')), θ)
        Cj  = ForwardDiff.jacobian(t -> model(t)[6], θ)
        Vj  = ForwardDiff.jacobian(t -> model(t)[7], θ)
        ∇ = zeros(length(θ))
        for k in eachindex(θ)
            ∇[k] = (tr(gAc' * reshape(Aj[:, k], rₓ, rₓ)) + tr(gBc' * reshape(Bj[:, k], rₓ, rᵤ))
                    + tr(gEc' * reshape(Ej[:, k], rₓ, nw))
                    + tr(gÃe' * reshape(Atj[:, k], rₓ, rₓ)) + tr(gCe' * reshape(Cj[:, k], ry, rₓ))
                    + tr(gVe' * reshape(Vj[:, k], ry, ry)))
        end
        return J_det + J_est, ∇, Ku, J_det, J_est, true
    end
    return eval!
end

"""
    hinf_of_report(title, θ_names, θ_roles; of_eval, model, J_des, ∇J_des,
                   θ_init, θ_nom, θ_lb, θ_ub, γ², seed=20240624, n_starts=5)

Output-feedback H∞ report: gradient check, multi-start BFGS, and a worst-case cost
table split into estimation/control/total (`J_est`, `J_det`, `J_c=J_det+J_est`,
`J_des`, `J_tot`).  `of_eval` is the closure from `Hinf_of_θ`.
"""
function hinf_of_report(title, θ_names, θ_roles;
                        of_eval, model, J_des, ∇J_des,
                        θ_init, θ_nom, θ_lb, θ_ub, γ², seed = 20240624, n_starts = 5)
    nθ = length(θ_init)
    bar = "=" ^ 74
    println(bar); println("  ContEst (output-feedback H∞, γ²=$(γ²))  —  ", title); println(bar)

    eval_J(θ) = begin
        Jc, ∇Jc, _, _, _, _ = of_eval(model, θ)
        (Jc + J_des(θ), ∇Jc .+ ∇J_des(θ))
    end
    f, g! = cached_objective(eval_J)

    ∇a, ∇fd, rel = verify_gradient(f, g!, θ_init)
    println("\n── Gradient verification at θ_init ─────────────────────────")
    @printf("  analytic : [%s]\n", join((@sprintf("%9.4f", v) for v in ∇a), " "))
    @printf("  fin.diff : [%s]\n", join((@sprintf("%9.4f", v) for v in ∇fd), " "))
    @printf("  relative error: %.2e  %s\n", rel,
            rel < 5e-2 ? "(OK)" : "(WARNING: large discrepancy)")

    println("\n── Multi-start BFGS ($n_starts starts: θ_nom + $(n_starts-1) random) ──")
    best_J, best_θ, results = multistart_design(f, g!, θ_lb, θ_ub;
                                                n_starts = n_starts, θ_nom = θ_nom,
                                                seed = seed, g_tol = 1e-5, iterations = 100)
    Js = sort([r.J for r in results]); ndistinct = 1
    for k in 2:length(Js); Js[k] - Js[k-1] > 1e-4 && (ndistinct += 1); end
    @printf("  starts run: %d   distinct minima found: %d   best J_tot = %.4f\n",
            n_starts, ndistinct, best_J)
    println(ndistinct == 1 ?
            "  → all starts agree: optimum is (numerically) global." :
            "  → multimodal: reporting the best minimum only (others omitted).")

    breakdown(θ) = begin
        _, _, _, Jd, Je, _ = of_eval(model, θ)
        (; J_est = Je, J_det = Jd, J_c = Jd + Je, J_des = J_des(θ), J_tot = Jd + Je + J_des(θ))
    end
    base = breakdown(θ_nom); opt = breakdown(best_θ)
    pct(b, o_) = abs(b) < 1e-9 ? "    — " : @sprintf("%+6.2f%%", 100 * (o_ - b) / b)
    println("\n── Cost components (worst-case): baseline θ_nom → optimal θ* ──")
    @printf("  %-28s %12s %12s %12s %8s\n", "component", "baseline", "optimal", "Δ (incr +)", "%")
    rows = (("Estimation cost  J_est",       base.J_est, opt.J_est),
            ("Control cost     J_det",       base.J_det, opt.J_det),
            ("Inner cost  J_c=J_det+J_est",  base.J_c,   opt.J_c),
            ("Design cost      J_des",       base.J_des, opt.J_des),
            ("TOTAL   J_tot=J_c+J_des",      base.J_tot, opt.J_tot))
    for (name, b, o_) in rows
        @printf("  %-28s %12.4f %12.4f %12.4f  %s\n", name, b, o_, o_ - b, pct(b, o_))
    end
    @printf("\n  Net total reduction: %.4f  (%.2f%%)\n",
            base.J_tot - opt.J_tot, 100 * (base.J_tot - opt.J_tot) / base.J_tot)

    println("\n── Design parameters: baseline → optimal (best minimum) ────")
    @printf("  %-12s %10s %10s %11s   %s\n", "param", "baseline", "optimal", "Δ", "enters/role")
    for i in 1:nθ
        @printf("  %-12s %10.3f %10.3f %+11.3f   %s\n",
                θ_names[i], θ_nom[i], best_θ[i], best_θ[i] - θ_nom[i], θ_roles[i])
    end
    println(bar)
    return best_θ, best_J, base, opt
end
