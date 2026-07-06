using ForwardDiff

function ∇ₓf(x₀::Vector, u₀::Vector, θ₀::Vector, f::Function)
	ForwardDiff.jacobian(x -> f(x, u₀, θ₀), x₀)
end

function ∇ᵤf(x₀::Vector, u₀::Vector, θ₀::Vector, f::Function)
	ForwardDiff.jacobian(u -> f(x₀, u, θ₀), u₀)
end

function ∇ₓh(x₀::Vector, u₀::Vector, θ₀::Vector, h::Function)
	ForwardDiff.jacobian(x -> h(x, u₀, θ₀), x₀)
end

function ∇ᵤh(x₀::Vector, u₀::Vector, θ₀::Vector, h::Function)
	ForwardDiff.jacobian(u -> h(x₀, u, θ₀), u₀)
end


function time_update(infoState::Tuple{Vector, Matrix}, u₀::Vector, θ₀::Vector, f::Function, Covars::Tuple{Matrix, Matrix})
	F = ∇ₓf(infoState[1], u₀, θ₀, f) # jacobian of state dynamics
	x = f(infoState[1], u₀, θ₀)
	Σ = F * infoState[2] * F' + Covars[1]
	return (x, Σ)
end

function measurement_update(infoState::Tuple{Vector, Matrix}, u₀::Vector, y₁::Vector, θ₀::Vector, h::Function, Covars::Tuple{Matrix, Matrix}; mode = "predict")
	H = ∇ₓh(infoState[1], u₀, θ₀, h) # jacobian of measurement dynamics
	L = infoState[2] * H' / (H * infoState[2] * H' + Covars[2] + 0.0001*I) #Kalman Gain
	if mode == "predict"
	x = infoState[1]
	else
	x = infoState[1] + L * (y₁ - h(infoState[1], u₀, θ₀))
	end
	Σ = (I - L * H) * infoState[2]
	return (x, Σ)
end

function update(infoState::Tuple{Vector, Matrix}, u₀::Vector, y₁::Vector, θ₀::Vector, f, h, Covars::Tuple{Matrix, Matrix}; mode = "predict")
	infoState = time_update(infoState, u₀, θ₀, f, Covars)
	infoState = measurement_update(infoState, u₀, y₁, θ₀, h, Covars; mode=mode)
	return infoState
end

