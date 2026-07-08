# example_adcs.jl – ContEst: spacecraft 3-axis attitude determination & control
# (ADCS) co-design under a shared spacecraft resource budget.
#
# Non-trivial reformulation
# ─────────────────────────
# The original "raise every precision / authority" formulation was trivial: with an
# independent quadratic price on each axis the optimizer just walks each parameter
# partway up its own bound. Here the sensing and actuation hardware instead draw on
# a SINGLE shared spacecraft POWER/MASS BUDGET, so they COMPETE: a better star
# tracker, a better gyro and a stronger reaction wheel each consume the same scarce
# resource, and buying one means starving another. The benefit of each saturates
# (diminishing returns), so the optimum is an interior ALLOCATION of the fixed
# budget, not a corner — and because attitude error is weighted far above rate, the
# budget flows to the star tracker (attitude sensing) and wheel authority while the
# rate gyro is deliberately starved. Which resource wins is not obvious a priori;
# ContEst computes the split.
#
# The dynamics are nonlinear through the gyroscopic coupling ω×(Jω) (zero at ω=0,
# so x=0,u=0 is a fixed point ∀θ) ⇒ the eKF–MPC inner solver.
#
# States   x = [ϕ, θ_a, ψ, ωx, ωy, ωz]   (small-angle attitude rad; body rates rad/s)
# Inputs   u = [τx, τy, τz]               (reaction-wheel torques, N·m)
# Outputs  y = [ϕ,θ_a,ψ (star tracker), ωx,ωy,ωz (rate gyro)]   (noise set by precision)
# Design   θ = [e_rw (wheel authority, in f) ; α_st, α_gyro (sensor precisions, in V)]
#            all three draw on one shared budget e_rw + α_st + α_gyro ≈ B.
#
# Run:  julia --project benchmarks/example_adcs.jl

include("../src/eKF.jl")
include("../src/MPCs.jl")
include("../src/Simulate.jl")
include("../src/BFGS.jl")
include("../report_contest.jl")

using LinearAlgebra, ForwardDiff, Printf

# ── Spacecraft parameters ─────────────────────────────────────────────────────
const Jx = 4.0      # principal moment of inertia, x-axis  [kg·m²]
const Jy = 6.0      # principal moment of inertia, y-axis  [kg·m²]
const Jz = 5.0      # principal moment of inertia, z-axis  [kg·m²]
const Jmat  = Diagonal([Jx, Jy, Jz]) |> Matrix
const Jinv  = inv(Jmat)
const dt_s  = 0.1   # star-tracker/gyro/control sampling period  [s]

# ── Nonlinear attitude dynamics  f(x, u, θ) ──────────────────────────────────
# Reaction-wheel torque enters scaled by the wheel-authority design parameter
# θ[1]=e_rw (same scalar on all three axes). Gyroscopic term ω×(Jω) is the
# nonlinearity (zero at ω=0), so x=0,u=0 is a fixed point for every θ.
function f_adcs(x, u, θ)
    ϕ, θa, ψ   = x[1], x[2], x[3]
    ωx, ωy, ωz = x[4], x[5], x[6]
    e_rw       = θ[1]
    ω   = [ωx, ωy, ωz]
    gyro = cross(ω, Jmat * ω)                  # genuine nonlinearity, = 0 at ω=0
    ω̇   = Jinv * (e_rw .* u .- gyro)
    return [ϕ  + dt_s * ωx,
            θa + dt_s * ωy,
            ψ  + dt_s * ωz,
            ωx + dt_s * ω̇[1],
            ωy + dt_s * ω̇[2],
            ωz + dt_s * ω̇[3]]
end

# ── Measurement  h(x,u,θ): star tracker (attitude) + rate gyro (rates) ────────
# The precision now enters the measurement NOISE V(θ), not the output gain, so h is
# the plain sensor readout.
h_adcs(x, ::AbstractVector, ::AbstractVector) = [x[1], x[2], x[3], x[4], x[5], x[6]]

# ── θ-DEPENDENT sensor noise: precision α sets the readout variance V = v0/α² ──
# A more precise (heavier/more-power) star tracker or gyro reports lower variance,
# with diminishing returns (∝ 1/α²).  α_st (θ[2]) covers the 3 attitude channels,
# α_gyro (θ[3]) the 3 rate channels.
const v_st0   = 0.20      # base star-tracker variance (at unit precision) — sensing matters
const v_gyro0 = 0.20      # base rate-gyro variance
const Σ_w = Diagonal([1e-6, 1e-6, 1e-6, 1e-4, 1e-4, 1e-4]) |> Matrix
Covars_θ(θ) = (Σ_w,
               Matrix(Diagonal([v_st0 / θ[2]^2, v_st0 / θ[2]^2, v_st0 / θ[2]^2,
                                v_gyro0 / θ[3]^2, v_gyro0 / θ[3]^2, v_gyro0 / θ[3]^2])))

# ── Cost matrices: penalise attitude error 4× the body rates ──────────────────
const Q  = Diagonal([200.0, 200.0, 200.0, 50.0, 50.0, 50.0])   # [ϕ,θ_a,ψ, ωx,ωy,ωz]
const R  = Diagonal([1.0, 1.0, 1.0])                           # [τx,τy,τz]
const s_x = zeros(6)
const s_u = zeros(3)

function stoch_cost(ζ, u)
    n = 6; xk = ζ[1:n]; Σk = vec_to_mat(ζ[n+1:end])
    return (xk - s_x)' * Q * (xk - s_x) + (u - s_u)' * R * (u - s_u) + tr(Q * Σk)
end
det_cost(x, u) = (x - s_x)' * Q * (x - s_x) + (u - s_u)' * R * (u - s_u)

# ── Input limits: |τ_i| ≤ 0.8 N·m (loose enough that wheel-authority benefit
# saturates before its bound — so e_rw is an interior trade, not pinned high) ──
constraint_fn(::AbstractVector, u) =
    [u[1] - 0.8; -u[1] - 0.8; u[2] - 0.8; -u[2] - 0.8; u[3] - 0.8; -u[3] - 0.8]

const N_hor = 15   # horizon (15 × 0.1 s = 1.5 s)
prob = ControlModel(f_adcs, h_adcs, Covars_θ, 6, 3, 6, N_hor, det_cost, stoch_cost, constraint_fn)
const mpc_eval = nonlinear_mpc_θ(prob)
const ce_mpc   = certainty_equivalence_mpc(prob)   # mean-only baseline controller

# ── Post-disturbance initial condition ────────────────────────────────────────
const x_ic = [0.15, -0.10, 0.12, 0.05, -0.04, 0.03]    # attitude + rate upset
const Σ_ic = Diagonal([2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2]) |> Matrix
const u_lin = zeros(3)

# ── Design cost: shared spacecraft POWER/MASS BUDGET ──────────────────────────
# Wheel authority and the two sensor precisions all draw on ONE budget; the soft
# equality (Σ resource ≈ B) makes them compete for it. Baseline spends the budget
# uniformly (unit on each); ContEst re-allocates it.
const θ_nom  = [1.0, 1.0, 1.0]          # e_rw, α_st, α_gyro (uniform baseline)
const B_res  = sum(θ_nom)               # shared resource budget
const c_b    = 30.0                     # budget-tightness price (soft equality Σθ ≈ B)
J_des(θ)  = c_b * (sum(θ) - B_res)^2
∇J_des(θ) = fill(2c_b * (sum(θ) - B_res), 3)

const θ_lb   = [0.3, 0.3, 0.3]
const θ_ub   = [3.0, 5.0, 5.0]
const θ_init = copy(θ_nom)

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

# ── Run ───────────────────────────────────────────────────────────────────────
if abspath(PROGRAM_FILE) == @__FILE__
    run_contest_report("Spacecraft ADCS — shared power/mass budget allocation",
                       ["e_rw", "α_st", "α_gyro"],
                       ["f (wheel authority)", "V (star tracker)", "V (rate gyro)"];
                       contest_f = contest_f, contest_g! = contest_g!,
                       θ_init = θ_init, θ_nom = θ_nom, θ_lb = θ_lb, θ_ub = θ_ub,
                       cost_breakdown = cost_breakdown, seed = 20240624, n_starts = 5,
                       surrogate = θ -> mpc_eval(x_ic, Σ_ic, u_lin, θ; grad = false)[2])
end
