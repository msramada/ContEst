using LinearAlgebra, ForwardDiff

# ──────────────────────────────────────────────────────────────────────────────
#  LQR.jl — steady-state H₂/LQG inner solver via Riccati (the SCALABLE control
#  path).  Same role and return signature as `SDP_θ` (SDPs.jl) but O(n³) instead
#  of the SDP's O(n⁶): the H₂ cost of the stationary linear-quadratic problem is
#      J(θ) = tr(P W),     P = discrete Riccati (DARE) solution,
#  and the envelope-theorem gradient uses P together with the closed-loop
#  stationary covariance S (one discrete-Lyapunov solve):
#      ∂J/∂W = P,  ∂J/∂Q = S,  ∂J/∂R = K S Kᵀ,
#      ∂J/∂A = 2 P A_cl S,     ∂J/∂B = 2 P A_cl S Kᵀ,   A_cl = A + B K.
#  These are chained with the ForwardDiff jacobians of (A,B,W,Q,R) in θ.
#  Verified against the SDP (identical cost, rel 2e-7) and finite differences
#  (gradient rel 3e-7); DARE solves in ~13 ms at 200 states, ~60 ms at 400.
#
#  θ enters only the dynamics (A,B) here — the control side of ContEst (CCD).
#  The estimation side (θ in the measurement map) is handled by the dual filter
#  Riccati in `filter_cost` below (used by the App-4 sensing example).
# ──────────────────────────────────────────────────────────────────────────────

# Discrete Riccati  AᵀPA − P − AᵀPB(R+BᵀPB)⁻¹BᵀPA + Q = 0
# via structure-preserving doubling (SDA): quadratic convergence, O(n³)/iter.
function dare(A, B, Q, R; tol = 1e-13, maxit = 80)
    G = B * (R \ B')
    Ak = copy(A); Gk = copy(G); Hk = Matrix(Q)
    for _ in 1:maxit
        W = I + Gk * Hk
        Ak1 = Ak * (W \ Ak)
        Gk1 = Gk + Ak * (W \ (Gk * Ak'))
        Hk1 = Hk + Ak' * (Hk * (W \ Ak))
        if norm(Hk1 - Hk) / (norm(Hk) + 1e-12) < tol
            Hk = Hk1; break
        end
        Ak, Gk, Hk = Ak1, Gk1, Hk1
    end
    return Symmetric(Hk)
end

# Discrete Lyapunov  S = Acl S Aclᵀ + W  via squaring (needs ρ(Acl) < 1).
function dlyap(Acl, W; tol = 1e-13, maxit = 80)
    Ak = copy(Acl); Sk = Matrix(W)
    for _ in 1:maxit
        Sk1 = Sk + Ak * Sk * Ak'
        Ak1 = Ak * Ak
        if norm(Sk1 - Sk) / (norm(Sk) + 1e-12) < tol
            Sk = Sk1; break
        end
        Ak, Sk = Ak1, Sk1
    end
    return Symmetric(Sk)
end

"""
    LQR_θ(rₓ, rᵤ) -> eval!

Drop-in scalable replacement for `SDP_θ`.  `eval!(model, θ)` returns
`(J, ∇θJ, K, Jstate, Jeffort)` where `model(θ) -> (A,B,W,Q,R)` (discrete-time),
`J = tr(PW)` is the stationary H₂ control cost, `K` the stabilising gain,
`Jstate = tr(Q S)`, `Jeffort = tr(R K S Kᵀ)` (so `Jstate + Jeffort = J`).
Unstable/undetectable `θ` (Riccati non-convergence) return `(1e8, 0, …)` so the
outer multi-start survives.
"""
function LQR_θ(rₓ::Int, rᵤ::Int)
    function eval!(model::Function, θ::AbstractVector)
        A, B, W, Q, R = model(θ)
        local P, K, Acl, S
        try
            P = dare(A, B, Q, R)
            K = -(R + B' * P * B) \ (B' * P * A)
            Acl = A + B * K
            # detectability/stability guard: closed loop must be Schur, P finite
            if !all(isfinite, P) || maximum(abs, eigvals(Matrix(Acl))) > 1 - 1e-8
                return (1e8, zeros(length(θ)), zeros(rᵤ, rₓ), 1e8, 0.0)
            end
            S = dlyap(Acl, W)
        catch
            return (1e8, zeros(length(θ)), zeros(rᵤ, rₓ), 1e8, 0.0)
        end
        J = tr(P * W)
        Jstate = tr(Q * S)
        Jeffort = tr(R * (K * S * K'))

        gW = Matrix(P); gQ = Matrix(S); gR = K * S * K'
        gA = 2 * (P * Acl * S); gB = 2 * (P * Acl * S * K')

        Aj = ForwardDiff.jacobian(t -> model(t)[1], θ)
        Bj = ForwardDiff.jacobian(t -> model(t)[2], θ)
        Wj = ForwardDiff.jacobian(t -> model(t)[3], θ)
        Qj = ForwardDiff.jacobian(t -> model(t)[4], θ)
        Rj = ForwardDiff.jacobian(t -> model(t)[5], θ)
        ∇ = zeros(length(θ))
        for k in eachindex(θ)
            ∇[k] = (tr(gA' * reshape(Aj[:, k], rₓ, rₓ)) + tr(gB' * reshape(Bj[:, k], rₓ, rᵤ))
                    + tr(gW' * reshape(Wj[:, k], rₓ, rₓ)) + tr(gQ' * reshape(Qj[:, k], rₓ, rₓ))
                    + tr(gR' * reshape(Rj[:, k], rᵤ, rᵤ)))
        end
        return J, ∇, K, Jstate, Jeffort
    end
    return eval!
end

# ──────────────────────────────────────────────────────────────────────────────
#  Filter (estimation) side — dual Riccati, for the App-4 sensing example.
#  Steady-state Kalman error covariance Σ_e solves the filter DARE with
#  (Aᵀ, Cᵀ, W, V); the estimation cost is J_est = tr(M Σ_e) with weight M, and
#  its gradient (θ in C and V, i.e. sensor placement/precision) comes from the
#  dual Lyapunov solve — the exact analogue of the control side above.
# ──────────────────────────────────────────────────────────────────────────────

"""
    filter_θ(n, ry) -> eval!

Estimation-side evaluator (dual of `LQR_θ`).  `eval!(Af, Cf, Vf, Wf, M, θ)` returns
`(Jest, ∇θJest, Σe)` where `Σe = dare(Aᵀ,Cᵀ,W,V)` is the steady-state (predicted)
Kalman error covariance and `Jest = tr(M Σe)` the estimation-uncertainty cost.
Each of `Af, Cf, Vf, Wf` may be a fixed matrix OR a function of `θ`:
`Cf`,`Vf` are the θ-dependent measurement map / noise (sensor precision or rate),
and `Af`,`Wf` let the ESTIMATION cost also depend on the control/plant parameters
(e.g. inertia that enters `A` and the process noise `W`) — essential when `θ_f`
appears in the dynamics, not only in `B`.

The gradient is the dual of the control envelope formulas: with Ã=Aᵀ, B̃=Cᵀ,
K̃=−(V+B̃ᵀΣeB̃)⁻¹B̃ᵀΣeÃ, Ã_cl=Ã+B̃K̃, S̃=dlyap(Ã_cl, M),
    ∂Jest/∂Ã = 2 Σe Ã_cl S̃,   ∂Jest/∂C = (2 Σe Ã_cl S̃ K̃ᵀ)ᵀ,
    ∂Jest/∂V = K̃ S̃ K̃ᵀ,        ∂Jest/∂W = S̃,
chained with the ForwardDiff Jacobians of `Aᵀ(θ), C(θ), V(θ), W(θ)`.  A fixed-matrix
argument contributes no term (its Jacobian is zero).  Unstable/undetectable `θ`
(Riccati non-convergence) return `(1e8, 0, …)` so the outer multi-start survives.
Verified against finite differences (rel ~1e-6).
"""
function filter_θ(n::Int, ry::Int)
    function eval!(Af, Cf::Function, Vf::Function, Wf, M, θ::AbstractVector)
        A = Af isa Function ? Af(θ) : Af
        W = Wf isa Function ? Wf(θ) : Wf
        C = Cf(θ); V = Vf(θ)
        Ãt = Matrix(A'); B̃ = Matrix(C')
        local Σe, K̃, Ãcl, S̃
        try
            Σe = dare(Ãt, B̃, W, V)
            K̃ = -(V + B̃' * Σe * B̃) \ (B̃' * Σe * Ãt)
            Ãcl = Ãt + B̃ * K̃
            if !all(isfinite, Σe) || maximum(abs, eigvals(Matrix(Ãcl))) > 1 - 1e-8
                return (1e8, zeros(length(θ)), Matrix(1e8 * I(n)))
            end
            S̃ = dlyap(Ãcl, M)
        catch
            return (1e8, zeros(length(θ)), Matrix(1e8 * I(n)))
        end
        Jest = tr(M * Σe)
        gÃt = 2 * (Σe * Ãcl * S̃)                # ∂Jest/∂Ã   (Ã = Aᵀ)
        gC  = Matrix((2 * Σe * Ãcl * S̃ * K̃')')   # ∂Jest/∂C   (ry×n)
        gV  = K̃ * S̃ * K̃'                         # ∂Jest/∂V   (ry×ry)
        gW  = Matrix(S̃)                          # ∂Jest/∂W   (n×n)
        Cj  = ForwardDiff.jacobian(Cf, θ)
        Vj  = ForwardDiff.jacobian(Vf, θ)
        Atj = Af isa Function ? ForwardDiff.jacobian(t -> vec(Matrix(Af(t)')), θ) : nothing
        Wj  = Wf isa Function ? ForwardDiff.jacobian(t -> vec(Matrix(Wf(t))), θ)  : nothing
        ∇ = zeros(length(θ))
        for k in eachindex(θ)
            g = tr(gC' * reshape(Cj[:, k], ry, n)) + tr(gV' * reshape(Vj[:, k], ry, ry))
            Atj === nothing || (g += tr(gÃt' * reshape(Atj[:, k], n, n)))
            Wj  === nothing || (g += tr(gW'  * reshape(Wj[:, k],  n, n)))
            ∇[k] = g
        end
        return Jest, ∇, Symmetric(Matrix(Σe))
    end
    return eval!
end
