# example_distillation.jl – ContEst: binary distillation column  (candidate #4)
#
# Physical setup
# ─────────────
# An Nst-stage binary column (stage 1 = reboiler … stage Nst = condenser)
# separates a light/heavy mixture. The plant is nonlinear through the vapour–
# liquid equilibrium  y(x) = α x / (1 + (α−1) x)  and the bubble-point temperature
# T(x), both of which are nonlinear in the stage composition.
#
# Medium-scale ContEst reformulation (feed and sensor placement)
# ───────────────────────────────────────────────────────────────
#   • θ_f : Nf FEED-STAGE LOCATIONS  f_1…f_Nf  (enter f). The total feed is
#     distributed over stages by a SUM of normalised Gaussians centred at the
#     feed locations; the cumulative feed above stage i sets that stage's
#     internal liquid traffic  L_i(θ_f) = L + F·Σ_{k>i} w_k(θ_f).  Feed
#     placement therefore reshapes the whole convective coupling of the column
#     nonlinearly — a genuine interior allocation of the feed streams, not a
#     "more is better" knob.
#   • θ_h : Ns TEMPERATURE-SENSOR LOCATIONS  s_1…s_Ns  (enter h). Each sensor is a
#     normalised-Gaussian weighted stage temperature; because T(x) is nonlinear,
#     the most informative trays are the steep mid-column stages — an interior
#     optimum, not the end stages.  Sensors compete to cover the profile.
#
# States   x = [Δx1 … ΔxNst]     (stage light-component mole-fraction deviations)
# Inputs   u = [ΔV]              (reboiler boil-up deviation)
# Outputs  y = [ T(s_1) … T(s_Ns) ]   (Ns movable temperature sensors)
#
# Design parameters  θ ∈ ℝ^{Nf + Ns}:
#   θ[1:Nf]              feed-stage locations   (enter f — placement)
#   θ[Nf+1 : Nf+Ns]      sensor-stage locations (enter h — placement)
#
# Run:  julia --project benchmarks/example_distillation.jl

include("../src/eKF.jl")
include("../src/MPCs.jl")
include("../src/Simulate.jl")
include("../src/BFGS.jl")

using LinearAlgebra, Optim, Printf

# ── Column parameters ─────────────────────────────────────────────────────────
const Nst   = 12                     # stages (1 = reboiler, Nst = condenser)
const Nf    = 3                      # feed streams (feed-location design axis)
const Ns    = 3                      # temperature sensors (sensing design axis)
const αrel  = 2.5                    # relative volatility
const Lref  = 2.0                    # external reflux (liquid) flow
const Vboil = 2.5                    # nominal boil-up (vapour) flow
const Ffeed = 2.5                    # total feed flow (large enough that feed
                                     # location materially reshapes liquid traffic)
const Mhold = 1.0                    # stage molar holdup
const dt_d  = 0.05                   # discretisation step
const ctr   = (Nst + 1) / 2          # column mid-point (where the upset lives)

# Nominal steady composition profile (light enriches up the column) and the
# corresponding VLE / temperature steady values.  Pure deviation form: every
# term below carries a Δ or u factor, so the origin is an equilibrium for all θ.
const x_ss = collect(range(0.15, 0.85, length = Nst))
vle(x)  = αrel * x / (1 + (αrel - 1) * x)                 # equilibrium curve
const y_ss = vle.(x_ss)
const TB0 = 380.0; const ΔTB = 40.0                      # bubble-point scale [K]
temp(x) = TB0 - ΔTB * vle(x)                             # nonlinear T(x)
const T_ss = temp.(x_ss)

# Feed distribution over stages (SUM of normalised Gaussians, one per feed
# location) and the resulting internal liquid-traffic profile L_i(θ_f).
const σ_f = 1.2
function feed_weights(flocs)
    g = zeros(eltype(flocs), Nst)
    for f in flocs
        g .+= [exp(-((i - f)^2) / (2 * σ_f^2)) for i in 1:Nst]
    end
    return g ./ sum(g)
end
function liquid_flow(flocs)                # L_i = L + F·(feed entering above stage i)
    w = feed_weights(flocs)
    return [Lref + Ffeed * sum(w[i+1:Nst]) for i in 1:Nst]
end

# ── Nonlinear column dynamics  f(x, u, θ)  (forward Euler) ─────────────────────
function f_dist(x, u, θ)
    Δx = x
    flocs = θ[1:Nf]
    Li = liquid_flow(flocs)
    Δy = [vle(x_ss[i] + Δx[i]) - y_ss[i] for i in 1:Nst]     # nonlinear VLE deviation
    out = map(1:Nst) do i
        xin  = i == Nst ? zero(eltype(Δx)) : Δx[i+1]         # liquid from stage above
        yin  = i == 1   ? zero(eltype(Δx)) : Δy[i-1]         # vapour from stage below
        ydsg = i == 1   ? zero(eltype(y_ss)) : y_ss[i-1]
        liq  = Li[i] * (xin - Δx[i])
        vap  = (Vboil) * (yin - Δy[i]) + u[1] * (ydsg - y_ss[i]) + u[1] * (yin - Δy[i])
        Δx[i] + dt_d / Mhold * (liq + vap)
    end
    return collect(out)
end

# ── Parametric output  h(x, u, θ) ─────────────────────────────────────────────
# Ns movable temperature sensors (inferential control) — each a normalised-
# Gaussian weighted stage temperature at its location s_j.  Placement sets how
# well the composition profile is estimated; because T(x) is nonlinear and the
# upset is mid-column, the best trays are interior, and the Ns sensors must
# spread to cover the profile rather than pile onto one tray.
const σ_sens = 1.2
function h_dist(x, ::AbstractVector, θ)
    Δx = x
    slocs = θ[Nf+1 : Nf+Ns]
    ΔTstage = [temp(x_ss[i] + Δx[i]) - T_ss[i] for i in 1:Nst]
    out = map(1:Ns) do j
        w = [exp(-((i - slocs[j])^2) / (2 * σ_sens^2)) for i in 1:Nst]
        ŵ = w ./ sum(w)
        sum(ŵ .* ΔTstage)
    end
    return collect(out)
end

# ── Fixed process and measurement noise ────────────────────────────────────────
const Σ_w    = Diagonal(fill(2e-5, Nst)) |> Matrix
const v0_sen = 0.6                                # temperature-sensor variance [K²]
const Σ_v    = Matrix(Diagonal(fill(v0_sen, Ns)))
const Covars = (Σ_w, Σ_v)

# ── Cost matrices ─────────────────────────────────────────────────────────────
# Penalise composition deviation across the column, weighted toward mid-column
# where the upset concentrates.
const Q  = Diagonal([20.0 + 30.0 * exp(-((i - ctr)^2) / (2 * (Nst / 5)^2)) for i in 1:Nst])
const R  = Diagonal([2.0])
const s_x = zeros(Nst)
const s_u = zeros(1)

det_cost(x, u) = (x - s_x)' * Q * (x - s_x) + (u - s_u)' * R * (u - s_u)
function stoch_cost(ζ, u)
    n  = Nst
    xk = ζ[1:n]; Σk = vec_to_mat(ζ[n+1:end])
    return (xk - s_x)' * Q * (xk - s_x) + (u - s_u)' * R * (u - s_u) + tr(Q * Σk)
end

# ── Input limit: |ΔV| ≤ 1.5 ───────────────────────────────────────────────────
constraint_fn(::AbstractVector, u) = [u[1] - 1.5; -u[1] - 1.5]

# ── Assemble the MPC evaluator ────────────────────────────────────────────────
const N_hor = 12
prob = ControlModel(f_dist, h_dist, Covars, Nst, 1, Ns, N_hor, det_cost, stoch_cost, constraint_fn)
const mpc_eval = nonlinear_mpc_θ(prob)

# ── Feed-composition upset initial condition (peaks mid-column) ───────────────
const x_ic  = [0.02 + 0.08 * exp(-((i - ctr)^2) / (2 * (Nst / 6)^2)) for i in 1:Nst]
const Σ_ic  = Diagonal(fill(1e-2, Nst)) |> Matrix
const u_lin = zeros(1)

# ── Design cost  J_des(θ) ─────────────────────────────────────────────────────
# Naive baseline: feeds clustered low (a common default), sensors placed low near
# the bottoms product. Location carries only mild regularisation (relocating a
# feed or sensor tray is nearly free), so the placement win is driven by column
# physics rather than by a hardware-allocation budget.
const f_nom = collect(range(3.0, 6.0, length = Nf))     # feeds clustered below middle
const s_nom = collect(range(2.0, 5.0, length = Ns))     # sensors low near bottoms
const θ_nom = vcat(f_nom, s_nom)

const λ_f  = 0.05                                       # feed-placement regularisation
const λ_s  = 0.05                                       # sensor-placement regularisation
J_des(θ) = λ_f * sum((θ[1:Nf] .- f_nom) .^ 2) +
           λ_s * sum((θ[Nf+1:Nf+Ns] .- s_nom) .^ 2)
∇J_des(θ) = vcat(2λ_f .* (θ[1:Nf] .- f_nom),
                 2λ_s .* (θ[Nf+1:Nf+Ns] .- s_nom))

const θ_lb  = vcat(fill(2.0, Nf), fill(2.0, Ns))
const θ_ub  = vcat(fill(Nst - 1.0, Nf), fill(Nst - 1.0, Ns))
const θ_init = copy(θ_nom)

# ── Objective + gradient (cached) ─────────────────────────────────────────────
const contest_f, contest_g! = contest_objective(mpc_eval, x_ic, Σ_ic, u_lin, J_des, ∇J_des)

# ── Sampled closed-loop cost breakdown (realised cost on the true trajectory) ──
# Both the baseline (θ_nom) and the ContEst optimum (θ*) run the SAME
# information-state MPC; the ONLY difference is the design θ, so the reported
# reduction isolates the pure co-design gain (no controller upgrade). Costs are
# Monte-Carlo averages of the realised cost over M = 100 sample trajectories, each
# a T = 30-step closed-loop rollout of the true nonlinear plant (the MPC still
# plans over its N = 12 horizon). Both J_det and J_est use the true state.
const M_mc  = 100        # Monte-Carlo sample trajectories
const T_sim = 30         # closed-loop simulation horizon (steps)
ctrl_info(xe, Σe, θ) = mpc_eval(xe, Σe, u_lin, θ; grad = false)[1]
function cost_breakdown(θ, ::Symbol)
    r  = simulate_mc(prob, ctrl_info, θ, Matrix(Q), x_ic, Σ_ic;
                     n_samples = M_mc, N_sim = T_sim)
    Jc = r.J_det + r.J_est
    return (; J_est = r.J_est, J_det = r.J_det, J_c = Jc, J_des = J_des(θ), J_tot = Jc + J_des(θ))
end

# ── Run ───────────────────────────────────────────────────────────────────────
if abspath(PROGRAM_FILE) == @__FILE__
    include("../report_contest.jl")
    θ_names = vcat(["f$k" for k in 1:Nf], ["s$j" for j in 1:Ns])
    θ_roles = vcat(["f (feed stage)"      for _ in 1:Nf],
                   ["h (temp sensor loc)" for _ in 1:Ns])
    run_contest_report("Binary distillation — multi-feed and sensor placement",
                       θ_names, θ_roles;
                       contest_f = contest_f, contest_g! = contest_g!,
                       θ_init = θ_init, θ_nom = θ_nom, θ_lb = θ_lb, θ_ub = θ_ub,
                       cost_breakdown = cost_breakdown, seed = 20240624, n_starts = 5)
end
