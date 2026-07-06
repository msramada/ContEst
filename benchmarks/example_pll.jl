# example_pll.jl – ContEst: co-designing the ESTIMATOR itself (PLL / dynamic
# state estimation) in a low-inertia, inverter-based power grid.
#
# Why this example
# ────────────────
# The other studies move the SENSING axis θ_h through what/where one senses: PMU
# reporting rate (WAMS), a shared sensor budget (ADCS), sensor location
# (distillation tray), or discrete placement (sensing). Here θ_h instead tunes a
# *dynamic estimator* — the phase-locked loop (PLL) that
# every grid-following (GFL) inverter uses to estimate grid frequency/phase for
# its fast-frequency / damping response. That is an estimator design variable
# classical control co-design cannot represent and sensor placement cannot
# capture: it is not WHERE we sense, but HOW the on-board estimator is tuned.
#
# The PLL bandwidth b_i is a genuine bias–variance knob:
#   • low  b_i : smooth, low-noise frequency estimate, but SLUGGISH — it lags the
#                true frequency, so the estimate the eKF/controller receives is
#                delayed → worse estimation and damping;
#   • high b_i : fast tracking, but the loop passes more measurement noise into
#                the frequency estimate → the reported frequency is noisier
#                (larger measurement covariance) → worse estimation.
# So b_i enters BOTH the dynamics f (the PLL tracking state) AND the measurement
# covariance V(θ) (the frequency-readout noise) — the first study to exercise the
# design-of-noise-levels extension (paper §"Co-designing the weights and noise
# levels"): V is here a function of θ, and the envelope-theorem gradient flows
# through it unchanged.
#
# Two design axes (ContEst):
#   • θ_f : per-bus inverter damping-actuator effectiveness  e_i ∈ [0.3,3]  (in f)
#   • θ_h : per-bus PLL bandwidth  b_i ∈ [2,20] rad/s                       (in f AND V)
#
# The disturbance is concentrated at bus 1 (large post-fault kick + prior), while
# buses 2,3 are quiet — so the informative-vs-quiet trade-off should push the PLL
# at the disturbed bus toward fast tracking and the quiet buses toward low noise:
# a spatially heterogeneous ESTIMATOR co-design.
#
# States  x = [Δδ_1..Δδ_B, Δω_1..Δω_B, Δωᵖ_1..Δωᵖ_B]  (angle, freq, PLL freq est.)
# Inputs  u = [ΔP_1..ΔP_B]                              (inverter damping power)
# Outputs y = [Δδ_1..Δδ_B (PMU angle), Δωᵖ_1..Δωᵖ_B (PLL frequency)]
#
# Run:  julia --project benchmarks/example_pll.jl

include("../src/eKF.jl")
include("../src/MPCs.jl")
include("../src/Simulate.jl")
include("../src/BFGS.jl")
include("../report_contest.jl")

using LinearAlgebra, ForwardDiff, Printf

# ── Network size ──────────────────────────────────────────────────────────────
const NB   = 3           # inverter buses (each with its own PLL)
const dt_s = 0.05        # sampling / control period [s]

# ── Heterogeneous, LOW-inertia machine/inverter parameters ────────────────────
const Hgen = [1.6, 2.2, 2.0]       # low inertia constants [s] (inverter-dominated)
const Dgen = [0.4, 0.6, 0.5]       # natural damping [pu]
const δeq  = [0.15, -0.05, 0.10]   # equilibrium angles [rad]

# Line synchronising coefficients (open chain 1–2–3, weak 1–3 tie)
const Kc = let K = zeros(NB, NB)
    K[1,2] = K[2,1] = 1.2
    K[2,3] = K[3,2] = 1.1
    K[1,3] = K[3,1] = 0.4
    K
end

# index helpers into the stacked state
δrange() = 1:NB
ωrange() = NB+1:2NB
prange() = 2NB+1:3NB

# ── Nonlinear dynamics f(x,u,θ) (deviation form ⇒ f(0,0,θ)=0 ∀θ) ──────────────
# PLL state Δωᵖ_i first-order-tracks the true frequency Δω_i with bandwidth b_i.
function f_pll(x, u, θ)
    Δδ = x[δrange()]; Δω = x[ωrange()]; Δωp = x[prange()]
    e  = θ[1:NB];     b  = θ[NB+1:2NB]
    Pe(i) = sum(Kc[i,j] * (sin(δeq[i]-δeq[j]+Δδ[i]-Δδ[j]) - sin(δeq[i]-δeq[j]))
                for j in 1:NB if j != i)
    ωdot(i)  = (-Pe(i) - Dgen[i]*Δω[i] + e[i]*u[i]) / (2*Hgen[i])
    return vcat([Δδ[i]  + dt_s*Δω[i]              for i in 1:NB],
                [Δω[i]  + dt_s*ωdot(i)            for i in 1:NB],
                [Δωp[i] + dt_s*b[i]*(Δω[i]-Δωp[i]) for i in 1:NB])   # PLL tracking
end

# ── Output map h(x,u,θ): PMU angle (direct) + PLL frequency (via the loop) ─────
function h_pll(x, ::AbstractVector, θ)
    Δδ = x[δrange()]; Δωp = x[prange()]
    return vcat([Δδ[i]  for i in 1:NB],      # bus-angle PMU channels
                [Δωp[i] for i in 1:NB])      # PLL frequency-estimate channels
end

# ── θ-DEPENDENT noise: PLL frequency-readout variance grows with bandwidth ─────
# A wider-bandwidth PLL passes more input noise to its output, so the effective
# measurement variance of the frequency channel is σ_ω²·(1 + κ·b_i). The angle
# PMU channels have fixed variance. This is the design-of-V extension.
const σ_δ2 = 4.0e-3      # angle-PMU measurement variance
const σ_ω2 = 2.0e-3      # base PLL-frequency measurement variance
const κ_pll = 0.30       # noise–bandwidth coefficient (per rad/s)

const Σ_w = Diagonal([fill(2.0e-5, NB);   # process noise: angle
                      fill(5.0e-4, NB);   #                frequency
                      fill(5.0e-4, NB)])|> Matrix   #        PLL state
Covars_θ(θ) = (Σ_w,
               Matrix(Diagonal(vcat(fill(σ_δ2, NB),
                                    [σ_ω2*(1 + κ_pll*θ[NB+i]) for i in 1:NB]))))

# ── Cost: regulate true angle & (heavily) true frequency; PLL states unweighted ─
const Q  = Diagonal([fill(15.0, NB); fill(120.0, NB); fill(0.0, NB)])  # perf. on TRUE x
const R  = Diagonal(fill(1.0, NB))
const s_x = zeros(3NB); const s_u = zeros(NB)

function stoch_cost(ζ, u)
    n = 3NB; xk = ζ[1:n]; Σk = vec_to_mat(ζ[n+1:end])
    return (xk - s_x)'*Q*(xk - s_x) + (u - s_u)'*R*(u - s_u) + tr(Q*Σk)
end
det_cost(x, u) = (x - s_x)'*Q*(x - s_x) + (u - s_u)'*R*(u - s_u)
constraint_fn(::AbstractVector, u) = vcat([u[i]-1.5 for i in 1:NB], [-u[i]-1.5 for i in 1:NB])

const N_hor = 12
prob = ControlModel(f_pll, h_pll, Covars_θ, 3NB, NB, 2NB, N_hor, det_cost, stoch_cost, constraint_fn)
const mpc_eval = nonlinear_mpc_θ(prob)
const ce_mpc   = certainty_equivalence_mpc(prob)   # mean-only baseline controller

# ── Post-disturbance initial condition: fault at BUS 1 ────────────────────────
# Bus 1 takes a large angle/frequency kick and a large prior; buses 2,3 quiet.
const x_ic = vcat([0.18, 0.03, 0.05],       # Δδ
                  [0.30, 0.05, 0.08],       # Δω
                  [0.0,  0.0,  0.0])        # Δωᵖ (PLL not yet locked to transient)
const σ_ic = [i == 1 ? 3.0e-2 : 6.0e-3 for i in 1:NB]
const Σ_ic = Diagonal(vcat(σ_ic, σ_ic, σ_ic)) |> Matrix
const u_lin = zeros(NB)

# ── Design cost: actuator upgrade + mild PLL-bandwidth regularisation ──────────
const λ_e = 0.5          # inverter damping-authority upgrade weight
const b_nom = 6.0        # naive default PLL bandwidth
const λ_b = 0.02         # mild regularisation of the PLL bandwidth about nominal

J_des(θ)  = λ_e*sum((θ[1:NB] .- 1).^2) + λ_b*sum((θ[NB+1:2NB] .- b_nom).^2)
∇J_des(θ) = vcat(2λ_e .* (θ[1:NB] .- 1), 2λ_b .* (θ[NB+1:2NB] .- b_nom))

const θ_lb  = [fill(0.3, NB); fill(2.0,  NB)]
const θ_ub  = [fill(3.0, NB); fill(20.0, NB)]
const θ_nom = [ones(NB); fill(b_nom, NB)]     # baseline: nominal actuators, naive PLL

# ── Sampled closed-loop cost breakdown (realised cost on the true trajectory) ──
# Baseline runs the certainty-equivalence controller (eKF mean only); the ContEst
# optimum runs the information-state MPC. Both J_det and J_est use the true state.
ctrl_info(xe, Σe, θ) = mpc_eval(xe, Σe, u_lin, θ; grad = false)[1]
ctrl_mean(xe, Σe, θ) = ce_mpc(xe, u_lin, θ)
function cost_breakdown(θ, role::Symbol)
    ctrl = role == :baseline ? ctrl_mean : ctrl_info
    r  = simulate_mc(prob, ctrl, θ, Matrix(Q), x_ic, Σ_ic)
    Jc = r.J_det + r.J_est
    return (; J_est = r.J_est, J_det = r.J_det, J_c = Jc, J_des = J_des(θ), J_tot = Jc + J_des(θ))
end

# ── ContEst objective (cached value+grad) and report ──────────────────────────
const contest_f, contest_g! = contest_objective(mpc_eval, x_ic, Σ_ic, u_lin, J_des, ∇J_des)

const θ_names = vcat(["e$i"    for i in 1:NB], ["b_pll$i" for i in 1:NB])
const θ_roles = vcat(["f (inverter damping)"    for _ in 1:NB],
                     ["h,V (PLL bandwidth)"     for _ in 1:NB])

if abspath(PROGRAM_FILE) == @__FILE__
    run_contest_report("Low-inertia grid: PLL / DSE estimator co-design",
                       θ_names, θ_roles;
                       contest_f = contest_f, contest_g! = contest_g!,
                       θ_init = θ_nom, θ_nom = θ_nom, θ_lb = θ_lb, θ_ub = θ_ub,
                       cost_breakdown = cost_breakdown, seed = 20240624, n_starts = 5)
end
