# example_wams.jl – ContEst: wide-area monitoring & damping of inter-area
# oscillations in a two-area power network.
#
# Non-trivial reformulation
# ─────────────────────────
# Both design axes are SCARCE-RESOURCE ALLOCATIONS with saturating benefit, so the
# optimum is a genuine interior, spatially heterogeneous trade — not "buy the most
# of everything" (every parameter pinned at its bound), the trivial outcome the
# original precision/effectiveness formulation produced.
#
#   • θ_h : per-bus PMU REPORTING RATE  r_i  under a shared COMMUNICATION BUDGET.
#       A faster PMU stream averages down the phasor readout variance,
#       V_i = v0 / r_i (more samples ⇒ less noise), but the wide-area link carries
#       a fixed total rate Σ_i r_i ≈ B (a soft budget in J_des). Rate is therefore
#       SCARCE and must be ALLOCATED across buses: spending it on a quiet,
#       well-known bus starves an informative one. The benefit saturates (∝1/r), so
#       the optimum pours the budget onto the disturbed area — a spatial sensing
#       allocation, the estimation-forward ContEst axis.
#
#   • θ_f : per-machine PRIMARY DAMPING GAIN  k_i  under a PROVISIONING COST.
#       Each machine augments its natural damping by a designed k_i (extra
#       converter/PSS damping). The closed-loop benefit SATURATES (diminishing
#       returns once the local mode is well damped) while provisioning damping
#       headroom costs c_k·Σ k_i, so the optimum is an interior profile that buys
#       damping where the disturbance actually lives (area 2) and little where it
#       does not — the classic non-monotone spatial trade of the kept
#       inertia/BESS studies, here on the nonlinear eKF–MPC path.
#
# Nonlinear (synchronising-power sin(δ) terms) ⇒ eKF–MPC inner solver.
#
# States  x = [Δδ_1..Δδ_N, Δω_1..Δω_N]   (rotor angle / speed deviations)
# Inputs  u = [ΔP_1..ΔP_N]               (supplementary damping power)
# Outputs y = [Δδ_1..Δδ_N, Δω_1..Δω_N]   (PMU angle/speed, all buses)
# Design  θ = [k_1..k_N (primary damping, in f) ; r_1..r_N (PMU rate, in V(θ))]
#
# Run:  julia --project benchmarks/example_wams.jl

include("../src/eKF.jl")
include("../src/MPCs.jl")
include("../src/Simulate.jl")
include("../src/BFGS.jl")
include("../report_contest.jl")

using LinearAlgebra, ForwardDiff, Printf

# ── Network: 4-machine, two-area (Kundur inter-area benchmark) ────────────────
const NM = 4            # generators (each a PMU bus and a damping actuator)
const NA = NM ÷ 2       # generators per area (two areas)
const dt_s = 0.05       # PMU / control sampling period [s]

area(i) = i ≤ NA ? 1 : 2

# ── Heterogeneous machine parameters ──────────────────────────────────────────
const Hgen = [4.0, 4.6, 3.7, 4.2][1:NM]      # inertia constants [s]
const Dgen = [1.0, 1.3, 0.8, 1.1][1:NM]      # natural damping  [pu]
const δeq  = [0.22, 0.08, -0.05, 0.30][1:NM] # equilibrium rotor angles [rad]

# Synchronising-coefficient matrix (strong intra-area, weak inter-area)
const Kc = let K = zeros(NM, NM), Kin = 2.6, Kout = 0.55
    for i in 1:NM, j in 1:NM
        i == j && continue
        K[i, j] = area(i) == area(j) ? Kin : Kout
    end
    K
end

# index helpers into the stacked state
δrange() = 1:NM
ωrange() = NM+1:2NM

# ── Nonlinear swing dynamics  f(x,u,θ)  (deviation form ⇒ f(0,0,θ)=0 ∀θ) ───────
function f_wams(x, u, θ)
    Δδ = x[δrange()]; Δω = x[ωrange()]; k = θ[1:NM]
    Pe(i) = sum(Kc[i, j] * (sin(δeq[i] - δeq[j] + Δδ[i] - Δδ[j]) - sin(δeq[i] - δeq[j]))
                for j in 1:NM if j != i)
    # speed dynamics: synchronising power, natural + designed PRIMARY DAMPING
    # (Dgen + k_i)·Δω, and the SUPPLEMENTARY wide-area MPC input u.
    ωdot(i) = (-Pe(i) - (Dgen[i] + k[i]) * Δω[i] + u[i]) / (2 * Hgen[i])
    return vcat([Δδ[i] + dt_s * Δω[i]   for i in 1:NM],
                [Δω[i] + dt_s * ωdot(i) for i in 1:NM])
end

# ── PMU output  h(x,u,θ): angle & speed of every bus (rate enters V, not h) ────
function h_wams(x, ::AbstractVector, θ)
    Δδ = x[δrange()]; Δω = x[ωrange()]
    return vcat([Δδ[i] for i in 1:NM], [Δω[i] for i in 1:NM])
end

# ── θ-DEPENDENT measurement noise: rate r_i sets the PMU readout variance ──────
# A faster PMU stream averages down the readout noise, V_i = v0 / r_i, on both the
# angle and the speed channel of bus i.  This is the design-of-V extension: r_i
# (θ_h) enters the measurement covariance rather than the output map h.
const v_δ0 = 0.15        # base angle-PMU variance (at unit rate) — noisy, so rate matters
const v_ω0 = 0.15        # base speed-PMU variance (at unit rate)
const Σ_w  = Diagonal([fill(1e-6, NM); fill(2e-3, NM)]) |> Matrix
Covars_θ(θ) = (Σ_w,
               Matrix(Diagonal(vcat([v_δ0 / θ[NM+i] for i in 1:NM],
                                    [v_ω0 / θ[NM+i] for i in 1:NM]))))

# ── Cost: penalise speeds (oscillation) ≫ angles ──────────────────────────────
const Q  = Diagonal([fill(20.0, NM); fill(100.0, NM)])
const R  = Diagonal(fill(1.0, NM))
const s_x = zeros(2NM); const s_u = zeros(NM)

function stoch_cost(ζ, u)
    n = 2NM; xk = ζ[1:n]; Σk = vec_to_mat(ζ[n+1:end])
    return (xk - s_x)' * Q * (xk - s_x) + (u - s_u)' * R * (u - s_u) + tr(Q * Σk)
end
det_cost(x, u) = (x - s_x)' * Q * (x - s_x) + (u - s_u)' * R * (u - s_u)

# Supplementary MPC authority is BOUNDED, so it cannot fully cancel a mis-tuned
# primary droop — the droop gain k_i must itself be designed.
const u_max = 0.5
constraint_fn(::AbstractVector, u) = vcat([u[i] - u_max for i in 1:NM], [-u[i] - u_max for i in 1:NM])

const N_hor = 12
prob = ControlModel(f_wams, h_wams, Covars_θ, 2NM, NM, 2NM, N_hor, det_cost, stoch_cost, constraint_fn)
const mpc_eval = nonlinear_mpc_θ(prob)
const ce_mpc   = certainty_equivalence_mpc(prob)   # mean-only baseline controller

# ── Post-disturbance initial condition: fault strikes AREA 2 ──────────────────
# Machines NA+1..NM take a large angle/speed kick and a large prior; area 1 quiet.
const x_ic = vcat([area(i) == 2 ? 0.10 : 0.03 for i in 1:NM] .* [(-1.0)^i for i in 1:NM],
                  [area(i) == 2 ? 0.15 : 0.04 for i in 1:NM] .* [(-1.0)^i for i in 1:NM])
const σ_ic = [area(i) == 2 ? 4e-2 : 3e-3 for i in 1:NM]   # prior uncertainty concentrated in area 2
const Σ_ic = Diagonal(vcat(σ_ic, σ_ic)) |> Matrix
const u_lin = zeros(NM)

# ── Design cost: damping provisioning (convex) + shared PMU COMMUNICATION BUDGET ─
# The damping-headroom cost is CONVEX (converter stress/losses grow faster than
# linearly with the damping gain), c_k·Σ k_i². Balanced against the mildly
# saturating damping benefit this yields a robust INTERIOR optimum k*_i ∝ (local
# marginal benefit), so damping concentrates where the disturbance is (area 2).
const k_nom  = 0.5                 # naive uniform primary-damping gain (baseline)
const c_k    = 0.5                 # convex damping-provisioning price
const r_nom  = 0.7                 # naive uniform PMU rate (scarce link)
const B_comm = NM * r_nom          # shared communication budget: total PMU rate
const c_b    = 12.0                # budget-tightness price (soft equality Σr ≈ B)

J_des(θ)  = c_k * sum(θ[1:NM] .^ 2) + c_b * (sum(θ[NM+1:2NM]) - B_comm)^2
∇J_des(θ) = vcat(2c_k .* θ[1:NM],
                 fill(2c_b * (sum(θ[NM+1:2NM]) - B_comm), NM))

const θ_lb  = [fill(0.0, NM); fill(0.3, NM)]
const θ_ub  = [fill(6.0, NM); fill(4.0, NM)]
const θ_nom = [fill(k_nom, NM); fill(r_nom, NM)]   # uniform damping + uniform PMU rate

# ── ContEst objective (cached value+grad) ─────────────────────────────────────
const contest_f, contest_g! = contest_objective(mpc_eval, x_ic, Σ_ic, u_lin, J_des, ∇J_des)

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

const θ_names = vcat(["k$i" for i in 1:NM], ["r$i" for i in 1:NM])
const θ_roles = vcat(["f: primary damping @ area $(area(i))" for i in 1:NM],
                     ["V: PMU rate @ area $(area(i))"        for i in 1:NM])

if abspath(PROGRAM_FILE) == @__FILE__
    run_contest_report("Two-area WAMS: comms-budgeted PMU rate + damping provisioning",
                       θ_names, θ_roles;
                       contest_f = contest_f, contest_g! = contest_g!,
                       θ_init = θ_nom, θ_nom = θ_nom, θ_lb = θ_lb, θ_ub = θ_ub,
                       cost_breakdown = cost_breakdown, seed = 20240624, n_starts = 5)
end
