# BFGS.jl — reusable outer-loop optimisation for ContEst co-design
#
# The examples define a plant (`ControlModel` → `mpc_eval`) and a design cost
# (`J_des`, `∇J_des`); everything needed to turn those into an optimal design is
# here, so an example only has to *call* these helpers rather than re-implement
# the cached objective, the gradient check, and the BFGS / multi-start loops.
#
# Pipeline:
#   contest_objective(...)  →  (f, g!)     cached value+gradient for Optim
#   verify_gradient(f, g!)  →  rel. error  (analytic vs finite differences)
#   bfgs_design(...)        →  one box-constrained BFGS solve
#   multistart_design(...)  →  many starts, keep the best (escape local minima)

using Optim, LinearAlgebra, Random, Printf

# ── Cached objective ──────────────────────────────────────────────────────────

"""
    cached_objective(eval_J) -> (f, g!)

Wrap a user cost `eval_J(θ) -> (J, ∇J)` (value and gradient computed together)
into the `(f, g!)` pair `Optim` expects, caching the last `θ` so the underlying
solve runs **once** per point even though `Optim` queries value and gradient
separately.
"""
function cached_objective(eval_J::Function)
    _θ = Ref(Float64[])
    _J = Ref(0.0)
    _g = Ref(Float64[])
    function refresh!(θ)
        if _θ[] != θ
            J, g = eval_J(θ)
            _θ[] = copy(θ); _J[] = J; _g[] = copy(g)
        end
    end
    f(θ)     = (refresh!(θ); _J[])
    g!(G, θ) = (refresh!(θ); G .= _g[])
    return f, g!
end

"""
    contest_objective(mpc_eval, x0, Σ0, u_lin, J_des, ∇J_des) -> (f, g!)

Convenience wrapper for the standard ContEst upper-level objective

    J(θ)  = J_c(θ) + J_des(θ),   ∇J(θ) = ∇J_c(θ) + ∇J_des(θ)

where `(J_c, ∇J_c)` come from the compiled MPC evaluator `mpc_eval` (the closure
returned by `nonlinear_mpc_θ`). Returns a cached `(f, g!)` ready for
[`bfgs_design`](@ref) / [`multistart_design`](@ref).
"""
function contest_objective(mpc_eval::Function, x0, Σ0, u_lin, J_des::Function, ∇J_des::Function)
    eval_J(θ) = begin
        _, Jc, ∇Jc, _ = mpc_eval(x0, Σ0, u_lin, θ)
        (Jc + J_des(θ), ∇Jc .+ ∇J_des(θ))
    end
    return cached_objective(eval_J)
end

# ── Gradient verification ─────────────────────────────────────────────────────

"""
    verify_gradient(f, g!, θ; h = 1e-5) -> (∇_analytic, ∇_fd, rel_err)

Compare the analytic gradient (`g!`) against forward finite differences of `f`
at `θ`. `rel_err = ‖∇_analytic − ∇_fd‖ / ‖∇_fd‖`; it should be `< 5e-2`
(typically `~1e-5`) before the optimiser output can be trusted.
"""
function verify_gradient(f::Function, g!::Function, θ; h = 1e-5)
    n  = length(θ)
    ∇a = zeros(n); g!(∇a, copy(θ))
    f0 = f(copy(θ))
    ∇fd = zeros(n)
    for i in 1:n
        θp = copy(θ); θp[i] += h
        ∇fd[i] = (f(θp) - f0) / h
    end
    return ∇a, ∇fd, norm(∇a - ∇fd) / (norm(∇fd) + 1e-10)
end

# ── BFGS (single start) ───────────────────────────────────────────────────────

"""
    bfgs_design(f, g!, θ_lb, θ_ub, θ_init; iterations=80, g_tol=1e-4,
                show_trace=false, kwargs...) -> Optim result

Box-constrained BFGS (`Fminbox(BFGS())`) from a single initial point, using the
analytic gradient `g!`. Extra keyword arguments are forwarded to
`Optim.Options`. The returned object is a standard `Optim` result, so
`result.minimizer`, `Optim.minimum(result)`, `Optim.converged(result)`,
`Optim.iterations(result)` all apply.
"""
function bfgs_design(f::Function, g!::Function, θ_lb, θ_ub, θ_init;
                     iterations = 80, g_tol = 1e-4, show_trace = false, kwargs...)
    return optimize(f, g!, θ_lb, θ_ub, θ_init,
                    Fminbox(BFGS()),
                    Optim.Options(; iterations = iterations, g_tol = g_tol,
                                  show_trace = show_trace, kwargs...))
end

# ── Multi-start BFGS ──────────────────────────────────────────────────────────

"""
    multistart_design(f, g!, θ_lb, θ_ub; n_starts=8, θ_nom=(θ_lb.+θ_ub)./2,
                      seed=nothing, verbose=false, kwargs...)
        -> (best_J, best_θ, results)

Run [`bfgs_design`](@ref) from `n_starts` initial points — `θ_nom` plus
`n_starts-1` samples drawn uniformly inside the box `[θ_lb, θ_ub]` — and keep the
best. `J(θ)` is generally non-convex in `θ`, so this both escapes local minima
and, when every start returns the same point, certifies the optimum as global.

Pass `seed` for reproducible random starts. `results` is a vector of
`(; x0, θ, J, iters)` named tuples, one per start. Remaining keyword arguments
are forwarded to `bfgs_design`.
"""
function multistart_design(f::Function, g!::Function, θ_lb, θ_ub;
                           n_starts = 8, θ_nom = (θ_lb .+ θ_ub) ./ 2,
                           seed = nothing, verbose = false, kwargs...)
    seed === nothing || Random.seed!(seed)
    d = length(θ_nom)
    starts = [copy(θ_nom)]
    for _ in 1:(n_starts - 1)
        push!(starts, θ_lb .+ (θ_ub .- θ_lb) .* rand(d))
    end

    best_J, best_θ = Inf, copy(θ_nom)
    results = NamedTuple[]
    for (i, x0) in enumerate(starts)
        r  = bfgs_design(f, g!, θ_lb, θ_ub, copy(x0); kwargs...)
        Js = Optim.minimum(r); θs = r.minimizer; it = Optim.iterations(r)
        push!(results, (x0 = x0, θ = θs, J = Js, iters = it))
        if Js < best_J
            best_J, best_θ = Js, θs
        end
        verbose && @printf("  start %2d:  J = %.6g   (%d iters)\n", i, Js, it)
    end
    return best_J, best_θ, results
end
