using JuMP, Clarabel, LinearAlgebra, ForwardDiff


mutable struct ControlModel
    f::Function # state dynamics
    h::Function # output dynamics
    # Noise covariances (W_process, V_measurement). Either a fixed tuple of
    # matrices, OR a Function θ -> (W, V) when the noise levels are themselves
    # design variables (e.g. a PLL bandwidth that sets the frequency-estimate
    # variance). In the latter case the θ-dependence flows through the eKF
    # covariance recursion and the envelope-theorem gradient picks it up
    # automatically via ForwardDiff — no change to the gradient machinery.
    Covars
    n::Int # original state dimension
    m::Int # control dimension
    o::Int # output dimension
    N::Int # prediction horizon
    runningcost::Function # original deterministic running cost
    stoch_runningcost::Function # stochastic running cost
    constraint_function::Function
end

# Resolve the covariance spec (fixed tuple or θ-dependent function) into a pair
# of dense matrices, keeping the eltype (Float64 or ForwardDiff.Dual) intact.
resolve_covars(Covars, θ₀) = let C = (Covars isa Function ? Covars(θ₀) : Covars)
    (Matrix(C[1]), Matrix(C[2]))
end

function info_state_dynamics(ζ::Vector, u₀::Vector, θ₀::Vector, f::Function, h::Function, Covars)
    C = resolve_covars(Covars, θ₀)
    n = size(C[1], 1)
	x = ζ[1:n]
	Σ_up_flat = ζ[n+1:end]
	Σ = vec_to_mat(Σ_up_flat)
	infoState = (x, Σ)
	(x⁺, Σ⁺) = update(infoState, u₀, [0.0], θ₀, f, h, C)
	return [x⁺; mat_uptriang_to_vec(Σ⁺)]
end

function mat_uptriang_to_vec(Σ_up::Matrix)
	n = size(Σ_up, 1)
	Σ_up_flat = Vector{eltype(Σ_up)}(undef, n * (n + 1) ÷ 2)
	ind = 0
	for j = 1:n
		for i = 1:j
			ind += 1
			Σ_up_flat[ind] = Σ_up[i, j]
		end
	end
	return Σ_up_flat
end

function vec_to_mat(Σ_up_flat::AbstractVector)
    L = length(Σ_up_flat)
    n = Int((sqrt(8L + 1) - 1) / 2)

    Σ = Matrix{eltype(Σ_up_flat)}(undef, n, n)
    ind = 1
    for j in 1:n
        for i in 1:j
            v = Σ_up_flat[ind]
            Σ[i, j] = v
            Σ[j, i] = v
            ind += 1
        end
    end
    return Σ
end


function nonlinear_mpc_θ(prob::ControlModel)
    n = prob.n
    n_info = n + n * (n + 1) ÷ 2
    # Define the optimization problem
    model = Model(Clarabel.Optimizer)
    set_silent(model)
    @variable(model, ζ[i = 1:n_info, t = 1:prob.N])
    @variable(model, u[i = 1:prob.m, t = 1:prob.N-1])
    objective = sum(prob.stoch_runningcost(ζ[:,t], u[:,t]) for t in 1:prob.N-1)
                    + prob.stoch_runningcost(ζ[:,prob.N], zeros(prob.m))
    @objective(model, Min, objective)

    
    @constraint(model, Symmetric(vec_to_mat(ζ[prob.n+1:end, 1])) -
        1e-6 * Matrix{Float64}(I, prob.n, prob.n) in PSDCone())
    @constraint(model,
        [t in 1:prob.N-1],
        prob.constraint_function(ζ[1:prob.n ,t+1], u[:,t]) .<= 0.0)

    # place holders
    c_dynamics = @constraint(model,
        [t in 1:prob.N-1],
        ζ[:,t+1] == ζ[:,t])
    c_initial_state = @constraint(model, ζ[:,1] == zeros(n_info))
    

    # Find optimal solution
    function eval(x₀₀::Vector{Float64}, Σ₀₀::Matrix{Float64}, u₀::Vector{Float64}, θ::Vector; grad::Bool = true)
        for t = 1:prob.N-1
            delete(model, c_dynamics[t])
        end
        delete(model, c_initial_state)

        ζ₀ = zeros(n_info)
        ζ₀[1:n] = x₀₀
        ζ₀[n+1:end] = mat_uptriang_to_vec(Σ₀₀)
        function myModel(θval::Vector)
            A_info = ForwardDiff.jacobian(ζvar -> info_state_dynamics(ζvar, u₀, θval, prob.f, prob.h, prob.Covars), ζ₀)
            B_info = ForwardDiff.jacobian(uvar -> info_state_dynamics(ζ₀, uvar, θval, prob.f, prob.h, prob.Covars), u₀)
            return A_info, B_info
        end
        A_info, B_info = myModel(θ)
        c_dynamics = @constraint(model,
        [t in 1:prob.N-1],
        ζ[:,t+1] == A_info * ζ[:,t] + B_info * u[:,t])
        c_initial_state = @constraint(model, ζ[:,1] == ζ₀)
        optimize!(model)
        #@assert termination_status(model) == MOI.OPTIMAL "Solver did not reach OPTIMAL."
        J = objective_value(model)
        # First receding-horizon input. Keep this separate from u₀: the parameter
        # gradient below must re-linearise about the SAME operating point (u₀) the
        # QP was built with, not the optimised input.
        u_first = value.(u)[:,1]

        feasibility = (primal_status(model) == FEASIBLE_POINT)
        ∇θJ = zeros(length(θ))
        # The dual/Jacobian contraction is only needed for the design gradient.
        # Closed-loop simulation only needs u_first, so callers pass grad=false to
        # skip this (the expensive part) entirely.
        if grad
            jacob_wrt_A = zeros(size(A_info))
            jacob_wrt_B = zeros(size(B_info))
            for k = 1:prob.N-1
                λ = JuMP.dual(c_dynamics[k])
                jacob_wrt_A .+= λ * value.(ζ[:,k])'
                jacob_wrt_B .+= λ * value.(u[:,k])'
            end
            A_jacobians_wrt_θ = ForwardDiff.jacobian(θvar -> myModel(θvar)[1], θ)
            B_jacobians_wrt_θ = ForwardDiff.jacobian(θvar -> myModel(θvar)[2], θ)
            for k in 1:length(θ)
                A_θk = reshape(A_jacobians_wrt_θ[:,k], n_info, n_info)
                B_θk = reshape(B_jacobians_wrt_θ[:,k], n_info, prob.m)
                ∇θJ[k] = tr(jacob_wrt_A' * A_θk) + tr(jacob_wrt_B' * B_θk)
            end
        end
        return u_first, J, ∇θJ, feasibility
    end
    return eval
end

# ── Certainty-equivalence (mean-only) MPC ─────────────────────────────────────
# A baseline controller that does NOT use the information state: it plans on the
# linearised MEAN dynamics x_{t+1} = A x_t + B u_t (A = ∂f/∂x, B = ∂f/∂u at the
# current eKF mean, no offset — matching the origin-through convention of the
# info-state MPC) and minimises the DETERMINISTIC running cost `prob.runningcost`
# only (no tr(QΣ) covariance term). It still consumes the eKF mean estimate as its
# initial state, but is blind to the covariance. Returns the first receding-horizon
# input; no gradient is produced (reporting/simulation use only).
function certainty_equivalence_mpc(prob::ControlModel)
    n = prob.n
    model = Model(Clarabel.Optimizer)
    set_silent(model)
    @variable(model, x[i = 1:n, t = 1:prob.N])
    @variable(model, u[i = 1:prob.m, t = 1:prob.N-1])
    @objective(model, Min,
        sum(prob.runningcost(x[:,t], u[:,t]) for t in 1:prob.N-1)
            + prob.runningcost(x[:,prob.N], zeros(prob.m)))
    @constraint(model,
        [t in 1:prob.N-1],
        prob.constraint_function(x[:,t+1], u[:,t]) .<= 0.0)
    c_dynamics = @constraint(model, [t in 1:prob.N-1], x[:,t+1] == x[:,t])
    c_initial_state = @constraint(model, x[:,1] == zeros(n))

    function eval(x₀::Vector{Float64}, u₀::Vector{Float64}, θ::Vector)
        for t = 1:prob.N-1
            delete(model, c_dynamics[t])
        end
        delete(model, c_initial_state)
        A = ∇ₓf(x₀, u₀, θ, prob.f)
        B = ∇ᵤf(x₀, u₀, θ, prob.f)
        c_dynamics = @constraint(model,
            [t in 1:prob.N-1],
            x[:,t+1] == A * x[:,t] + B * u[:,t])
        c_initial_state = @constraint(model, x[:,1] == x₀)
        optimize!(model)
        return value.(u)[:,1]
    end
    return eval
end

