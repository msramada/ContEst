# Simulate.jl — Monte-Carlo closed-loop evaluation of the realised ContEst cost.
#
# The reported estimation/control comparison comes from SAMPLED TRAJECTORIES of the
# ORIGINAL cost, not the open-loop MPC prediction. For each sample we draw a true
# initial state, run the true nonlinear plant forward under a receding-horizon
# controller, run the eKF on the noisy measurements, and accumulate the realised
# cost. Both components are evaluated on the TRUE state:
#
#   J_det = Σ_t ‖x_true_t − s_x‖²_Q + ‖u_t − s_u‖²_R      (prob.runningcost, true state)
#   J_est = Σ_t (x_true_t − x_est_t)' Q (x_true_t − x_est_t)   (realised estimation error)
#   J_c   = J_det + J_est
#
# averaged over `n_samples` independent noise/IC realisations. The `controller` is
# a unified closure (x_est, Σ_est, θ) -> u, so the caller selects the info-state
# MPC (ContEst design) or the certainty-equivalence mean-only MPC (baseline).
#
# Requires src/eKF.jl and src/MPCs.jl (update, resolve_covars, ControlModel).

using LinearAlgebra, Random

# Cholesky factor for sampling N(0, S); small jitter keeps borderline-PSD Σ safe.
_chol_L(S) = cholesky(Symmetric(Matrix(S)) + 1e-12 * I).L

function simulate_mc(prob::ControlModel, controller::Function, θ::Vector,
                     Q::AbstractMatrix, x_ic::Vector, Σ_ic::AbstractMatrix;
                     N_sim::Int = prob.N, n_samples::Int = 200, seed::Int = 20240624)
    Random.seed!(seed)
    W, V = resolve_covars(prob.Covars, θ)
    W = Matrix(W); V = Matrix(V)
    Lw  = _chol_L(W)
    Lv  = _chol_L(V)
    Lic = _chol_L(Σ_ic)

    Jdet = 0.0
    Jest = 0.0
    for _ in 1:n_samples
        x_true = x_ic .+ Lic * randn(prob.n)
        x_est  = copy(x_ic)
        Σ_est  = Matrix(Σ_ic)
        for _t in 1:N_sim
            u = controller(x_est, Σ_est, θ)
            e = x_true .- x_est
            Jdet += prob.runningcost(x_true, u)
            Jest += dot(e, Q * e)
            # advance the TRUE plant and take a noisy measurement
            x_true = prob.f(x_true, u, θ) .+ Lw * randn(prob.n)
            y = prob.h(x_true, u, θ) .+ Lv * randn(prob.o)
            # eKF measurement-mode update on the realised measurement
            (x_est, Σ_est) = update((x_est, Σ_est), u, y, θ, prob.f, prob.h, (W, V); mode = "update")
            Σ_est = (Σ_est + Σ_est') / 2      # numerical symmetrisation
        end
    end
    return (; J_det = Jdet / n_samples, J_est = Jest / n_samples)
end
