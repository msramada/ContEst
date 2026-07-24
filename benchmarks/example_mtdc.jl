# example_mtdc.jl — ContEst App 5: multi-terminal HVDC (MTDC) DC-voltage droop
# coordination + sparse voltage-sensor allocation.  BASE LIBRARY (model + inner
# solvers); the full co-design study + figure live in example_mtdc_sparse.jl.
#
# MERGED two-axis co-design on ONE robust inner value (30 states, NT=15 terminals):
#
#   control axis  θ_f = k   droop gains (enter A)   — H∞ game (tr(XW)) / block-diag BRL SDP
#   sensing axis  θ_h = α   sensor gains (enter C)   — H₂ estimation SDP, J_est = tr(Q Σ^ε)
#
# The droop enters A(θ), so it reshapes BOTH the H∞ worst-case control value AND the
# estimation error dynamics; the two axes therefore couple and the co-design does not
# separate.  The sensing side is priced by the H₂ estimation SDP in the observability-
# gramian "P^ε" form (paper eq. Jest):
#     min tr(W P^ε) + tr(V Z)  s.t. [Z Fᵀ; F P^ε]⪰0, [P^ε−Q  AᵀP^ε−CᵀFᵀ; P^εA−FC  P^ε]⪰0,
# with F = P^ε G (G the filter gain) and continuous sensor gains α_i∈[0,1] scaling the
# rows of C(α).  Exactly as on the control side we restrict P^ε and F to be
# BLOCK-DIAGONAL (one 2×2 block per terminal — a decentralised filter): its true cost is
# within ~2% of the centralised Kalman value (validated below), it solves in O(n) instead
# of the dense SDP's O(n⁶), and its LMI dual gives the EXACT envelope gradients
#     ∂J_est/∂A = −2 P^ε S₁₂ᵀ,   ∂J_est/∂C = 2 Fᵀ S₁₂ᵀ.
# α is an OUTER design variable (like the droop): in the P^ε form C enters bilinearly
# with the gain F, so α cannot be an SDP variable; it is optimised by BFGS through the
# envelope gradient, with an ℓ1 penalty selecting the sensor set.  There are NO
# covariance caps: every reported estimation number is the TRUE cost tr(Q Σ^ε_true),
# Σ^ε_true = dlyap(A−GC, W+GVGᵀ) — the steady-state error covariance of the SDP's gain.
#
# Why H∞ (control): the failure mode is a RESONANT PEAK — converter lag makes excessive
# droop provoke a lightly-damped inter-terminal DC-voltage oscillation whose worst-case
# (not RMS) amplification threatens the grid.  H∞ prices that.
#
# States  x = [ΔV_i, ΔI_i] per terminal (DC-bus voltage, converter current)
# Inputs  u = supplementary power-reference modulation per terminal
# Outputs y = ΔV_i (DC-bus voltages telemetered; converter currents estimated)
# Design  θ = (droop k_i ; sensor gains α_i)
#
# Run:  julia --project benchmarks/example_mtdc.jl        (quick validity + droop sanity)

include("../src/LQR.jl")     # dare, dlyap
include("../src/BFGS.jl")    # cached_objective, multistart_design, verify_gradient
include("../src/Hinf.jl")    # hinf_gare (H∞ control game Riccati, dense reference)
using JuMP, Clarabel, LinearAlgebra, Printf, ForwardDiff

# ── Reduced DC grid: NT=15 terminals ⇒ 2·NT = 30 states ───────────────────────
const NT  = 15
const dt  = 0.005                 # fast DC dynamics ⇒ small step [s]
const γ²  = 16.0                  # H∞ disturbance-attenuation level (exact envelope gradient)
nb(i) = (i == 1 ? NT : i - 1, i == NT ? 1 : i + 1)      # ring neighbours
idx(i) = (2i - 1, 2i)                                    # per-terminal 2×2 block indices

const Cdc  = [0.8 + 0.4sin(0.9i) for i in 1:NT]         # DC-node capacitances
const τcnv = [0.02 + 0.01cos(1.1i) for i in 1:NT]       # converter current lag [s]
const Gdc  = 5.0                                         # DC line conductance (ring)
const k0   = 12.0                                        # droop scale (leverage)
const σinf = [1.0 + 1.5 * (0.5 + 0.5sin(1.7i + 2)) for i in 1:NT]   # infeed intensity

const Qmat = Matrix(Diagonal([iseven(s) ? 0.2 : 60.0 for s in 1:2NT]))  # penalise ΔV ≫ ΔI
const Rmat = Matrix(2.0 * I(NT))
const v_meas = 1.0e-3             # DC-voltage measurement variance (per sensor)

const n = 2NT; const m = NT
const Iidx = [2i for i in 1:NT]              # converter-current states (estimated)

# measurement map C(α): sensor i reads ΔV_i = state 2i-1, scaled by gain α_i
Cmat(α) = (C = zeros(eltype(α), m, n); for i in 1:m; C[i, 2i-1] = α[i]; end; C)
const Vmat = Matrix(v_meas * I, m, m)        # per-sensor measurement covariance

# ── Linear model  θ ↦ (A,B,W)  (discrete-time); droop k0·θ_i enters A ─────────
function model(θ)
    T = eltype(θ)
    Ac = zeros(T, n, n); Bc = zeros(T, n, m)
    for i in 1:NT
        v = 2i - 1; c = 2i; (l, r) = nb(i)
        Ac[v, v] = -2Gdc / Cdc[i]
        Ac[v, 2l-1] += Gdc / Cdc[i]
        Ac[v, 2r-1] += Gdc / Cdc[i]
        Ac[v, c] = 1.0 / Cdc[i]
        Ac[c, c] = -1.0 / τcnv[i]
        Ac[c, v] = -(k0 * θ[i]) / τcnv[i]          # droop feedback (θ enters A)
        Bc[c, i] = 1.0 / τcnv[i]
    end
    Ad = Matrix(1.0I, n, n) + dt * Ac
    Bd = dt * Bc
    W = zeros(n, n)
    for i in 1:NT
        W[2i-1, 2i-1] = (dt * σinf[i] / Cdc[i])^2      # infeed disturbance on DC node
        W[2i, 2i] = 1e-8
    end
    return Ad, Bd, W
end

# ── Design cost weights and box ──────────────────────────────────────────────
const c_k = 5e-3                 # converter-stress price per unit droop
const λs_default = 0.04          # ℓ1 sensor-sparsity weight
const θ_lb = fill(0.5, NT); const θ_ub = fill(4.0, NT); const θ_nom = fill(3.0, NT)

# ── Control inner value: H∞ worst-case cost via the game Riccati (dense reference) ──
function Jdet_ric(θ; μ = γ²)
    A, B, W = model(θ); E = Matrix(sqrt(Symmetric(Matrix(W))))
    X, _, _, _, ok = hinf_gare(A, B, E, Qmat, Rmat, μ); ok ? tr(X * Matrix(W)) : 1e6
end

# ── Estimation inner value: H₂ SDP in the observability-gramian (P^ε) form ─────
#   min tr(W P) + tr(V Z)  s.t.  [Z Fᵀ; F P] ⪰ 0,  [P−Q  AᵀP−CᵀFᵀ; PA−FC  P] ⪰ 0
# pattern ∈ (:full, :blockdiag, :banded); P,F restricted to the ring block structure.
# Returns (J, G, gA, gC, status): SDP value, filter gain G=P⁻¹F, and the exact envelope
# gradients gA=∂J/∂A (n×n), gC=∂J/∂C (m×n) from the dual of the second LMI.
function est_sdp(A, C, W, V, Q; pattern = :blockdiag, withgrad = true)
    nn = size(A, 1); ry = size(C, 1)
    mdl = Model(Clarabel.Optimizer); set_silent(mdl)
    if pattern == :full
        Pv = @variable(mdl, [1:nn, 1:nn], Symmetric); P = Matrix(Pv)
        Fv = @variable(mdl, [1:nn, 1:ry]); F = Matrix(Fv)
    else
        Pe = zeros(AffExpr, nn, nn); Fe = zeros(AffExpr, nn, ry)
        for i in 1:NT
            Pi = @variable(mdl, [1:2, 1:2], Symmetric); r = idx(i)
            for a in 1:2, b in 1:2; Pe[r[a], r[b]] = Pi[a, b]; end
            nbrs = pattern == :blockdiag ? (i,) : (i, i == 1 ? NT : i - 1, i == NT ? 1 : i + 1)
            for j in nbrs
                Fij = @variable(mdl, [1:2, 1:1]); r2 = idx(i)
                Fe[r2[1], j] = Fij[1, 1]; Fe[r2[2], j] = Fij[2, 1]
            end
        end
        if pattern == :banded
            for i in 1:NT
                j = i == NT ? 1 : i + 1; Pij = @variable(mdl, [1:2, 1:2]); ri = idx(i); rj = idx(j)
                for a in 1:2, b in 1:2; Pe[ri[a], rj[b]] = Pij[a, b]; Pe[rj[b], ri[a]] = Pij[a, b]; end
            end
        end
        P = Matrix(Pe); F = Matrix(Fe)
    end
    @variable(mdl, Z[1:ry, 1:ry], Symmetric)
    @constraint(mdl, Symmetric([Z permutedims(F); F P]) in PSDCone())
    lmi2 = @constraint(mdl, Symmetric([P .- Q  (permutedims(A) * P .- permutedims(C) * permutedims(F));
                                       (P * A .- F * C)  P]) in PSDCone())
    @objective(mdl, Min, tr(W * P) + tr(V * Z))
    optimize!(mdl)
    st = termination_status(mdl)
    (st in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.NUMERICAL_ERROR)) ||
        return (1e6, zeros(nn, ry), zeros(nn, nn), zeros(ry, nn), st)
    J = objective_value(mdl)
    Pval = value.(P); Fval = value.(F)
    G = Symmetric(Matrix(Pval)) \ Matrix(Fval)          # filter gain G = (P^ε)⁻¹ F
    withgrad || return (J, G, zeros(nn, nn), zeros(ry, nn), st)
    S = dual(lmi2); S12 = S[1:nn, nn+1:2nn]
    gA = -2 .* (Matrix(Pval) * permutedims(S12))                 # ∂J/∂A = −2 P S₁₂ᵀ
    gC = 2 .* (permutedims(Matrix(Fval)) * permutedims(S12))     # ∂J/∂C =  2 Fᵀ S₁₂ᵀ
    (J, G, gA, gC, st)
end

# True steady-state estimation cost of the SDP filter gain on the FULL plant:
#   ε⁺ = (A−GC)ε + w − Gv  ⇒  Σ^ε_true = dlyap(A−GC, W + G V Gᵀ),  J_est = tr(Q Σ^ε_true).
function est_true_cost(A, C, W, V, Q, G)
    Acl = A - G * C
    Σ = dlyap(Acl, Matrix(W) + G * V * permutedims(G))
    tr(Q * Σ)
end

# Centralised (dense) Kalman predictor cost — reference for the block-diagonal validity.
function kalman_cost(A, C, W, V, Q)
    Σ = dare(Matrix(A'), Matrix(C'), Matrix(W), Matrix(V)); tr(Q * Σ)
end

# Estimation value + gradients as a function of (θ, α): builds A(θ), C(α), solves est_sdp.
function est_val_grad(θ, α; pattern = :blockdiag)
    A, _, W = model(θ); C = Cmat(α)
    J, G, gA, gC, st = est_sdp(A, C, W, Vmat, Qmat; pattern = pattern)
    (J, G, gA, gC, st)
end

# ── Standalone sanity: block-diagonal estimation validity + game-Riccati droop ──
function mtdc_report()
    bar = "="^74
    println(bar); println("  ContEst App 5 — MTDC base library sanity ($n states, NT=$NT)"); println(bar)
    A, B, W = model(θ_nom); C = Cmat(ones(m))
    Jk = kalman_cost(A, C, W, Vmat, Qmat)
    println("\n── H₂ estimation SDP validity at baseline θ=3 (all $NT sensors) ──")
    @printf("  %-12s %12s %14s %12s\n", "pattern", "SDP tr(WP)", "true tr(QΣ)", "vs Kalman")
    for pat in (:full, :blockdiag, :banded)
        J, G, _, _, st = est_sdp(A, C, W, Vmat, Qmat; pattern = pat)
        Jt = est_true_cost(A, C, W, Vmat, Qmat, G)
        @printf("  %-12s %12.4f %14.4f %+11.1f%%\n", pat, J, Jt, 100*(Jt-Jk)/Jk)
    end
    @printf("  centralised Kalman reference tr(QΣ) = %.4f\n", Jk)
    println(bar)
end

if abspath(PROGRAM_FILE) == @__FILE__
    mtdc_report()
end
