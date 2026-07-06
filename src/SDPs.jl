using JuMP, Clarabel, LinearAlgebra, ForwardDiff, Printf

# ──────────────────────────────────────────────────────────────────────────────
#  SDPs.jl — steady-state covariance/H₂ LMI inner solver for the ContEst SDP path
#
#  Infinite-horizon, control-only (state-feedback) co-design: for a linear model
#      x⁺ = A(θ) x + B(θ) u + w,   w ~ N(0,W),
#  the stationary LQG/H₂ cost is the value of the semidefinite program
#      min_{Σ⪰0, Z₀, L}  tr(Q Σ) + tr(R Z₀)
#      s.t.  [Z₀ L; Lᵀ Σ] ⪰ 0                          (input-cost Schur: Z₀ ⪰ LΣ⁻¹Lᵀ)
#            [Σ−W  AΣ+BL; (AΣ+BL)ᵀ  Σ] ⪰ 0             (discrete Lyapunov LMI)
#  with stabilising gain K = L Σ⁻¹.  θ enters only the dynamics (A,B) — this is
#  the control side of ContEst (CCD); the estimation side (θ in the measurement
#  map h) requires the MPC-eKF path instead.
#
#  Gradient by the envelope theorem from the Lyapunov-LMI dual S (blocks S11,S12):
#      ∂J/∂A = −2 S12 Σᵀ,  ∂J/∂B = −2 S12 Lᵀ,  ∂J/∂W = S11,  ∂J/∂Q = Σ,  ∂J/∂R = Z₀
#  chained with the ForwardDiff jacobians of (A,B,W,Q,R) w.r.t. θ.
# ──────────────────────────────────────────────────────────────────────────────

function SDP_θ(rₓ::Int, rᵤ::Int)
    model = Model(Clarabel.Optimizer)
    set_silent(model)
    # Decision variables.  Σ and Z₀ are declared Symmetric so *all* their entries
    # (not just the upper triangle) are constrained — required for a non-diagonal
    # Q/R, otherwise the lower triangle is free and tr(QΣ) can be unbounded.
    @variable(model, Σ[1:rₓ, 1:rₓ], Symmetric)
    @variable(model, Z₀[1:rᵤ, 1:rᵤ], Symmetric)
    @variable(model, L[1:rᵤ, 1:rₓ])
    @constraint(model, Symmetric(Σ - 1e-5 * I(rₓ)) in PSDCone())
    @constraint(model, Symmetric([Z₀ L; L' Σ]) in PSDCone())
    c_lyapunov = @constraint(model,
        Symmetric([Σ zeros(rₓ, rₓ); zeros(rₓ, rₓ) Σ]) in PSDCone())

    function eval!(myLinearModel::Function, θval::AbstractVector)
        A0, B0, W0, Q0, R0 = myLinearModel(θval)
        @objective(model, Min, tr(Q0 * Σ) + tr(R0 * Z₀))
        delete(model, c_lyapunov)
        c_lyapunov = @constraint(model,
            Symmetric([(Σ - W0)        (A0 * Σ + B0 * L);
                       (A0 * Σ + B0 * L)'       Σ]) in PSDCone())
        optimize!(model)

        # Guard infeasibility: a random BFGS start can produce an unstable A(θ)
        # for which the Lyapunov LMI is infeasible. Return a large finite cost
        # and zero gradient instead of asserting, so the multi-start survives.
        if primal_status(model) != FEASIBLE_POINT
            return (1e8, zeros(length(θval)), zeros(rᵤ, rₓ), 1e8, 0.0)
        end

        J = objective_value(model)
        Σstar = value.(Σ)
        Lstar = value.(L)
        Z0star = value.(Z₀)
        Jstate = tr(Q0 * Σstar)          # state-regulation share  tr(QΣ)
        Jeffort = tr(R0 * Z0star)        # control-effort share    tr(RZ₀)

        S = JuMP.dual(c_lyapunov)        # size (2rₓ)×(2rₓ)
        S11 = S[1:rₓ, 1:rₓ]
        S12 = S[1:rₓ, rₓ+1:2rₓ]
        gradA = -2 .* (S12 * Σstar')     # ∂J/∂A   (rₓ×rₓ)
        gradB = -2 .* (S12 * Lstar')     # ∂J/∂B   (rₓ×rᵤ)
        gradW = S11                      # ∂J/∂W
        gradQ = Σstar                    # ∂J/∂Q
        gradR = Z0star                   # ∂J/∂R

        # Chain with the sensitivities of (A,B,W,Q,R) to θ (differentiate the
        # PASSED model closure, not any global).
        A_jac = ForwardDiff.jacobian(θ -> myLinearModel(θ)[1], θval)
        B_jac = ForwardDiff.jacobian(θ -> myLinearModel(θ)[2], θval)
        W_jac = ForwardDiff.jacobian(θ -> myLinearModel(θ)[3], θval)
        Q_jac = ForwardDiff.jacobian(θ -> myLinearModel(θ)[4], θval)
        R_jac = ForwardDiff.jacobian(θ -> myLinearModel(θ)[5], θval)
        ∇θJ = zeros(length(θval))
        for k in 1:length(θval)
            A_θk = reshape(A_jac[:, k], rₓ, rₓ)
            B_θk = reshape(B_jac[:, k], rₓ, rᵤ)
            W_θk = reshape(W_jac[:, k], rₓ, rₓ)
            Q_θk = reshape(Q_jac[:, k], rₓ, rₓ)
            R_θk = reshape(R_jac[:, k], rᵤ, rᵤ)
            ∇θJ[k] = (tr(gradA' * A_θk) + tr(gradB' * B_θk)
                      + tr(gradW' * W_θk) + tr(gradQ' * Q_θk) + tr(gradR' * R_θk))
        end
        K_LMI = Lstar / Σstar
        return J, ∇θJ, K_LMI, Jstate, Jeffort
    end
    return eval!
end

# ──────────────────────────────────────────────────────────────────────────────
#  Outer-loop wrappers — mirror BFGS.jl's contest_objective / GUIDE §6 so the SDP
#  path plugs into the SAME optimisation helpers (verify_gradient, bfgs_design,
#  multistart_design) as the MPC-eKF examples.  `cached_objective`,
#  `verify_gradient`, `multistart_design` come from BFGS.jl (included after this).
# ──────────────────────────────────────────────────────────────────────────────

"""
    sdp_contest_objective(sdp_eval, model, J_des, ∇J_des) -> (f, g!)

Cached value+gradient pair for the SDP co-design objective
`J(θ) = J_c(θ) + J_des(θ)`, where `(J_c, ∇J_c)` come from the compiled SDP
evaluator `sdp_eval` (the closure returned by `SDP_θ`) applied to the linear
model `model(θ) -> (A,B,W,Q,R)`.
"""
function sdp_contest_objective(sdp_eval::Function, model::Function,
                               J_des::Function, ∇J_des::Function)
    eval_J(θ) = begin
        J, ∇J, _, _, _ = sdp_eval(model, θ)
        (J + J_des(θ), ∇J .+ ∇J_des(θ))
    end
    return cached_objective(eval_J)   # from BFGS.jl
end

"""
    sdp_cost_breakdown(sdp_eval, model, J_des, θ)
        -> (; J_state, J_effort, J_c, J_des, J_tot)

Control-only cost decomposition (analogue of GUIDE §6). `J_c = tr(QΣ)+tr(RZ₀)`
splits into the state-regulation share `J_state` and the control-effort share
`J_effort`; `J_des` is the hardware cost; `J_tot = J_c + J_des`.
"""
function sdp_cost_breakdown(sdp_eval::Function, model::Function, J_des::Function, θ)
    _, _, _, Jstate, Jeffort = sdp_eval(model, θ)
    Jc = Jstate + Jeffort
    return (; J_state = Jstate, J_effort = Jeffort, J_c = Jc,
              J_des = J_des(θ), J_tot = Jc + J_des(θ))
end

"""
    sdp_report(title, θ_names, θ_roles; sdp_eval, model, J_des, ∇J_des,
               θ_init, θ_nom, θ_lb, θ_ub, seed=20240624, n_starts=5)

SDP analogue of `report_contest.jl`: gradient check, multi-start BFGS
(baseline θ_nom + n_starts−1 random), control-only cost table and θ table.
Returns `(best_θ, best_J, base, opt)`.
"""
function sdp_report(title, θ_names, θ_roles;
                    sdp_eval, model, J_des, ∇J_des,
                    θ_init, θ_nom, θ_lb, θ_ub, seed = 20240624, n_starts = 5)
    nθ = length(θ_init)
    bar = "=" ^ 74
    println(bar); println("  ContEst (SDP path)  —  ", title); println(bar)

    sdp_f, sdp_g! = sdp_contest_objective(sdp_eval, model, J_des, ∇J_des)

    # ── Gradient verification ────────────────────────────────────────────────
    ∇a, ∇fd, rel = verify_gradient(sdp_f, sdp_g!, θ_init)   # from BFGS.jl
    println("\n── Gradient verification at θ_init ─────────────────────────")
    @printf("  analytic : [%s]\n", join((@sprintf("%9.4f", v) for v in ∇a), " "))
    @printf("  fin.diff : [%s]\n", join((@sprintf("%9.4f", v) for v in ∇fd), " "))
    @printf("  relative error: %.2e  %s\n", rel,
            rel < 5e-2 ? "(OK)" : "(WARNING: large discrepancy)")

    # ── Multi-start BFGS (baseline + random restarts) ────────────────────────
    println("\n── Multi-start BFGS ($n_starts starts: θ_nom + $(n_starts-1) random) ──")
    best_J, best_θ, results = multistart_design(sdp_f, sdp_g!, θ_lb, θ_ub;
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
    base = sdp_cost_breakdown(sdp_eval, model, J_des, θ_nom)
    opt  = sdp_cost_breakdown(sdp_eval, model, J_des, best_θ)
    pct(b, o_) = abs(b) < 1e-9 ? "    — " : @sprintf("%+6.2f%%", 100 * (o_ - b) / b)
    println("\n── Cost components: baseline θ_nom → optimal θ* ────────────")
    @printf("  %-28s %12s %12s %12s %8s\n", "component", "baseline", "optimal", "Δ (incr +)", "%")
    rows = (("State regulation  J_state", base.J_state, opt.J_state),
            ("Control effort    J_effort", base.J_effort, opt.J_effort),
            ("Control cost (H₂) J_c",      base.J_c,     opt.J_c),
            ("Design cost       J_des",    base.J_des,   opt.J_des),
            ("TOTAL             J_tot",    base.J_tot,   opt.J_tot))
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
