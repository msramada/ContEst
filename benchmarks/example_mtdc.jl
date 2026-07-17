# example_mtdc.jl — ContEst App 5: multi-terminal HVDC (MTDC) DC-voltage droop
# coordination + sparse voltage-sensor allocation with estimation-quality caps.
#
# MERGED two-axis co-design on ONE robust inner value (30 states, NT=15 terminals):
#
#   control axis  θ_f = k   droop gains (enter A)   — H∞ game Riccati, J_det = tr(XW)
#   sensing axis  θ_h = α   sensor gains (enter C)   — ℓ1 covariance-SDP, J_est = tr(QΣ)
#
# The droop enters A(θ), so it reshapes BOTH the H∞ worst-case control value AND the
# estimation prior Σ_pred(θ) = dlyap(A(θ),W); the two axes therefore couple and the
# co-design does not separate. Estimation/allocation is posed as the convex
# information-form estimation SDP (paper eq. Jest): the posterior error-covariance
# bound X ⪰ (Σ_pred(θ)^{-1} + Σ_j α_j C_jᵀV_j^{-1}C_j)^{-1}, with sparse ℓ1 sensor
# selection and PER-CONVERTER-CURRENT estimation-error caps  [X]_{2i,2i} ≤ tol.
# Those caps + discrete selection make the estimation half genuinely SDP-only
# (no algebraic Riccati expresses per-state covariance caps with an active set);
# the cap duals λ_i are the shadow prices that localize the essential sensors.
#
# Why H∞ (control): the failure mode is a RESONANT PEAK — converter lag makes
# excessive droop provoke a lightly-damped inter-terminal DC-voltage oscillation
# whose worst-case (not RMS) amplification threatens the grid. H∞ prices that.
#
# States  x = [ΔV_i, ΔI_i] per terminal (DC-bus voltage, converter current)
# Inputs  u = supplementary power-reference modulation per terminal
# Outputs y = ΔV_i (DC-bus voltages telemetered; converter currents estimated)
# Design  θ = (droop k_i ; sensor gains α_i)
#
# Run:  julia --project benchmarks/example_mtdc.jl

include("../src/LQR.jl")     # dare, dlyap
include("../src/BFGS.jl")    # cached_objective, multistart_design, verify_gradient
include("../src/Hinf.jl")    # hinf_gare (H∞ control game Riccati)
using JuMP, Clarabel, LinearAlgebra, Printf, ForwardDiff

# ── Reduced DC grid: NT=15 terminals ⇒ 2·NT = 30 states ───────────────────────
const NT  = 15
const dt  = 0.005                 # fast DC dynamics ⇒ small step [s]
const γ²  = 8.0                   # H∞ disturbance-attenuation level
nb(i) = (i == 1 ? NT : i - 1, i == NT ? 1 : i + 1)      # ring neighbours

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
# per-sensor information (sensor j measures ΔV_j = state 2j-1); θ-independent
const Gsens = [ (e = zeros(n); e[2j-1] = 1.0; (e * e') / v_meas) for j in 1:NT ]

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

# ── Inner values ─────────────────────────────────────────────────────────────
# H∞ worst-case control cost via the game Riccati (indefinite-weight DARE).
function Jdet(θ)
    A, B, W = model(θ); E = Matrix(sqrt(Symmetric(Matrix(W))))
    X, _, _, _, ok = hinf_gare(A, B, E, Qmat, Rmat, γ²)
    ok ? tr(X * W) : 1e6
end
# Estimation prior information = inverse stationary predicted covariance Σ_pred(θ)^{-1}.
Sinfo(α) = sum(α[j] * Gsens[j] for j in 1:m)          # measurement information (θ-indep)
Yprior(θ) = begin A, B, W = model(θ); inv(Symmetric(Matrix(dlyap(A, Matrix(W))))) end
# Posterior error-covariance (information form) for sensor gains α.
estX(θ, α) = inv(Symmetric(Matrix(Yprior(θ) + Sinfo(α))))
Jest(θ, α) = tr(Qmat * estX(θ, α))
curv(θ, α) = (X = estX(θ, α); [X[j, j] for j in Iidx])

# ── Droop objective value + EXACT envelope gradient ──────────────────────────
# J(θ) = J_det(θ) + J_est(θ,α) + c_k Σθ, with θ entering ONLY through A(θ):
#   ∂J_det/∂A = 2 X A_cl S           (H∞ game-Riccati envelope, no solver diff.)
#   ∂J_est/∂A = 2 Λ A Σ_pred,  Λ = dlyap(Aᵀ, Σ_pred⁻¹ X Q X Σ_pred⁻¹)
# both contracted against the model Jacobian ∂A/∂θ (ForwardDiff of A only).
function droop_val_grad(α, θ)
    A, B, W = model(θ); Wm = Matrix(W)
    E = Matrix(sqrt(Symmetric(Wm)))
    X, K̃, Ku, Acl, ok = hinf_gare(A, B, E, Qmat, Rmat, γ²)
    ok || return (1e6, zeros(NT))
    Sc = dlyap(Acl, Wm); Jdet_ = tr(X * Wm)
    gA_det = 2 * (X * Acl * Sc)                         # ∂J_det/∂A
    P  = dlyap(A, Wm)                                   # Σ_pred(θ)
    Pi = inv(Symmetric(Matrix(P)))
    Xe = inv(Symmetric(Matrix(Pi + Sinfo(α))))
    Jest_ = tr(Qmat * Xe)
    GP = Symmetric(Pi * Xe * Qmat * Xe * Pi)           # ∂J_est/∂Σ_pred
    Λ  = dlyap(Matrix(A'), Matrix(GP))                 # Aᵀ Λ A - Λ + GP = 0
    gA_est = 2 * (Λ * A * P)                            # ∂J_est/∂A
    gA = gA_det .+ gA_est
    Aj = ForwardDiff.jacobian(t -> vec(model(t)[1]), θ)
    g  = [ tr(gA' * reshape(Aj[:, k], n, n)) + c_k for k in 1:NT ]
    return (Jdet_ + Jest_ + c_k * sum(θ), g)
end

# ── Droop step: multi-start projected BFGS on the EXACT envelope gradient ─────
function droop_opt(α; n_starts = 5)
    f, g! = cached_objective(θ -> droop_val_grad(α, θ))
    _, θbest, _ = multistart_design(f, g!, θ_lb, θ_ub; n_starts = n_starts, θ_nom = θ_nom, seed = 1)
    return θbest
end

# ── Sensor step: convex ℓ1 covariance-SDP with per-current caps ───────────────
#   min tr(QX) + λs Σα   s.t.  [Yprior+ΣαⱼGⱼ  I; I  X] ⪰ 0,  [X]_{2i,2i} ≤ tol,  0≤α≤1
function sensor_sdp(θ; λs = λs_default, tol = Inf, capped = false)
    Y = Matrix(Yprior(θ))
    mdl = Model(Clarabel.Optimizer); set_silent(mdl)
    @variable(mdl, 0 <= α[1:m] <= 1)
    @variable(mdl, X[1:n, 1:n], Symmetric)
    Yinfo = Y + sum(α[j] * Gsens[j] for j in 1:m)
    @constraint(mdl, Symmetric([Yinfo Matrix(I, n, n); Matrix(I, n, n) X]) in PSDCone())
    caprefs = capped ? [@constraint(mdl, X[j, j] <= tol) for j in Iidx] : ConstraintRef[]
    @objective(mdl, Min, tr(Qmat * X) + λs * sum(α))
    optimize!(mdl)
    st = termination_status(mdl); ok = st in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL)
    (ok ? value.(α) : fill(NaN, m)),
    (ok && capped ? [-dual(c) for c in caprefs] : Float64[]), st
end

# ── Report ───────────────────────────────────────────────────────────────────
function mtdc_report()
    bar = "="^74
    println(bar); println("  ContEst App 5 — MTDC droop + sparse sensing co-design (H∞, γ²=$γ², $n states)"); println(bar)
    allon = ones(NT)
    # exact envelope gradient vs finite differences (droop outer loop)
    fchk, gchk! = cached_objective(θ -> droop_val_grad(allon, θ))
    _, _, rel = verify_gradient(fchk, gchk!, θ_nom)
    @printf("\n  droop envelope-gradient check (analytic vs finite diff): rel err = %.2e\n", rel)
    Jd0, Je0 = Jdet(θ_nom), Jest(θ_nom, allon)
    θ1 = droop_opt(allon); Jd1, Je1 = Jdet(θ1), Jest(θ1, allon)
    println("\n── A. Droop co-design (all $NT sensors) ──")
    @printf("  baseline θ=3 : J_det=%.4f  J_est=%.4f  J_c=%.4f\n", Jd0, Je0, Jd0 + Je0)
    @printf("  ContEst θ*∈[%.2f,%.2f]: J_det=%.4f (-%.1f%%)  J_est=%.4f (-%.1f%%)  J_c=%.4f (-%.1f%%)\n",
            minimum(θ1), maximum(θ1), Jd1, 100*(Jd0-Jd1)/Jd0, Je1, 100*(Je0-Je1)/Je0,
            Jd1+Je1, 100*((Jd0+Je0)-(Jd1+Je1))/(Jd0+Je0))
    Jdes0 = c_k*sum(θ_nom); Jdes1 = c_k*sum(θ1)
    @printf("  J_tot (incl. converter stress J_des): %.4f → %.4f  (-%.1f%%)\n",
            Jd0+Je0+Jdes0, Jd1+Je1+Jdes1, 100*((Jd0+Je0+Jdes0)-(Jd1+Je1+Jdes1))/(Jd0+Je0+Jdes0))

    base  = maximum(curv(θ1, allon))
    prior = maximum(diag(inv(Symmetric(Matrix(Yprior(θ1)))))[Iidx])
    @printf("  achievable current est-var band at θ*: floor=%.3e … no-sensor prior=%.3e\n", base, prior)

    println("\n── B. ℓ1 sensor frontier at θ* (covariance-SDP) ──")
    @printf("  %-8s %-10s %-9s %s\n", "λs", "#sensors", "J_est", "max curr var")
    for λ in (0.0, 0.01, 0.02, 0.04)
        α, _, _ = sensor_sdp(θ1; λs = λ); on = Float64.(α .> 0.5)
        @printf("  %-8.2f %2d/%-7d %-9.4f %.3e\n", λ, count(>(0.5), α), NT, Jest(θ1, on), maximum(curv(θ1, on)))
    end

    println("\n── C. Capped ℓ1-SDP at θ* (λs=0.04): per-current caps [Σ]_{2i,2i} ≤ tol ──")
    @printf("  %-11s %-9s %-9s %-13s %s\n", "tol", "#sensors", "J_est", "max curr var", "binding (λ_i)")
    for frac in (0.75, 0.70, 0.65, 0.62, 0.58)
        tol = frac * prior
        α, λd, st = sensor_sdp(θ1; λs = 0.04, tol = tol, capped = true)
        if st in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL)
            on = Float64.(α .> 0.5)
            act = [(i, round(λd[i], digits = 1)) for i in 1:NT if λd[i] > 1e-2]
            @printf("  %-11.3e %2d/%-7d %-9.4f %.3e     %s\n", tol, count(>(0.5), α), NT,
                    Jest(θ1, on), maximum(curv(θ1, on)), isempty(act) ? "none" : string(act))
        else
            @printf("  %-11.3e INFEASIBLE (%s)\n", tol, st)
        end
    end

    println("\n── D. Merged joint optimum (capped menu, then re-optimize droop) ──")
    αc, _, _ = sensor_sdp(θ1; λs = 0.04, tol = 0.65 * prior, capped = true); menu = Float64.(αc .> 0.5)
    θ2 = droop_opt(menu); Jd2, Je2 = Jdet(θ2), Jest(θ2, menu)
    @printf("  kept %d/%d sensors: %s\n", count(>(0), menu), NT, string(findall(>(0), menu)))
    @printf("  droop θ**∈[%.2f,%.2f]:  J_det=%.4f  J_est=%.4f  J_c=%.4f\n",
            minimum(θ2), maximum(θ2), Jd2, Je2, Jd2 + Je2)
    @printf("  vs baseline (θ=3, all %d sensors): J_c %.4f→%.4f  (-%.1f%%),  sensors %d→%d\n",
            NT, Jd0 + Je0, Jd2 + Je2, 100*((Jd0+Je0)-(Jd2+Je2))/(Jd0+Je0), NT, count(>(0), menu))
    Jdes0 = c_k*sum(θ_nom); Jdes2 = c_k*sum(θ2)
    @printf("  J_tot (incl. J_des): %.4f → %.4f  (-%.1f%%)\n",
            Jd0+Je0+Jdes0, Jd2+Je2+Jdes2, 100*((Jd0+Je0+Jdes0)-(Jd2+Je2+Jdes2))/(Jd0+Je0+Jdes0))
    println(bar)
    return (; θ1, Jd0, Je0, Jd1, Je1, θ2, Jd2, Je2, menu)
end

if abspath(PROGRAM_FILE) == @__FILE__
    mtdc_report()
end
