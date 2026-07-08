# make_paper_figs.jl — closed-loop simulation + figures for the ContEst paper.
#
#   julia --project make_paper_figs.jl adcs   # (or distillation / pll)
#
# Includes one example as a library (its report is guarded), recomputes the
# optimal design θ* with BFGS, simulates the closed loop (true nonlinear plant +
# eKF + receding-horizon MPC) under the BASELINE and the OPTIMAL design with the
# SAME noise realization, and writes a 2-panel figure to Paper_LaTeX/figs/.
# It also prints the cost breakdown and θ* used to fill the LaTeX tables.

ENV["GKSwstype"] = "100"   # headless GR

const CASE = isempty(ARGS) ? "adcs" : ARGS[1]
const SPEC = Dict(
    "adcs" => (file = "example_adcs.jl",
               title = "Spacecraft ADCS",
               sig   = "pointing error norm  [rad]",
               sidx  = nothing,
               tunit = "t  [s]", dt = 0.10),
    "distillation" => (file = "example_distillation.jl",
               title = "Binary distillation column",
               sig   = "mid-column composition dev. dx10",
               sidx  = 10,
               tunit = "t  [s]", dt = 0.05),
    "pll" => (file = "example_pll.jl",
               title = "Low-inertia grid: PLL/DSE co-design",
               sig   = "bus-1 frequency dev. dw1",
               sidx  = 4,                       # Δω_1 (disturbed bus), state layout [Δδ;Δω;Δωᵖ]
               tunit = "t  [s]", dt = 0.05),
)[CASE]

include(joinpath("benchmarks", SPEC.file))   # defines prob, mpc_eval, x_ic, Σ_ic, u_lin, Q, θ_nom, θ_lb, θ_ub, contest_f/g!, cost_breakdown
using Plots, Random, LinearAlgebra, Printf

# ── designs to compare ────────────────────────────────────────────────────────
# Multi-start (baseline θ_nom + 4 random restarts), keep the best — matches the
# reported table (report_contest.jl) and escapes local minima.
_, θstar, _ = multistart_design(contest_f, contest_g!, θ_lb, θ_ub;
                                n_starts = 5, θ_nom = θ_nom, seed = 20240624,
                                g_tol = 1e-5, iterations = 100)
θ_base = θ_nom
lab_base, lab_opt = "baseline", "ContEst (optimal)"

# ── closed-loop rollout: true plant + eKF + receding-horizon MPC ──────────────
function closed_loop(θ; T = 80, seed = 7)
    Random.seed!(seed)
    n, m, o = prob.n, prob.m, prob.o
    # Covariances may be a fixed tuple or a θ-dependent function (e.g. PLL noise).
    covars = prob.Covars isa Function ? prob.Covars(θ) : prob.Covars
    W, V = covars
    Lw = cholesky(Symmetric(Matrix(W)) + 1e-12I).L
    Lv = cholesky(Symmetric(Matrix(V)) + 1e-12I).L
    Lic = cholesky(Symmetric(Matrix(Σ_ic)) + 1e-12I).L

    xh = copy(x_ic)                 # estimate mean (prior)
    xt = x_ic .+ Lic * randn(n)     # true state drawn from the prior
    Σ  = copy(Matrix(Σ_ic))
    up = copy(u_lin)
    X  = zeros(n, T + 1); X[:, 1] = xt
    for k in 1:T
        u, _, _, _ = mpc_eval(xh, Σ, up, θ)                 # MPC on current estimate
        xt = prob.f(xt, u, θ) .+ Lw * randn(n)              # true nonlinear plant + process noise
        y  = prob.h(xt, u, θ) .+ Lv * randn(o)              # measurement
        (xh, Σ) = update((xh, Σ), u, y, θ, prob.f, prob.h, covars; mode = "correct")
        X[:, k + 1] = xt; up = u
    end
    return X
end

# representative closed-loop signal (ADCS: pointing-error norm; others: one state)
sigtrace(X) = CASE == "adcs" ? vec(sqrt.(sum(abs2, X[1:3, :]; dims = 1))) : X[SPEC.sidx, :]

# ── cost breakdowns ───────────────────────────────────────────────────────────
# Role-aware: the eKF–MPC examples define cost_breakdown(θ, role) (baseline =
# certainty-equivalence, optimum = info-state MPC); distillation has the 1-arg form.
bd(θ, role) = applicable(cost_breakdown, θ, role) ? cost_breakdown(θ, role) : cost_breakdown(θ)
b = bd(θ_base, :baseline); o_ = bd(θstar, :optimum)

default(fontfamily = "Computer Modern", legendfontsize = 8, guidefontsize = 9,
        tickfontsize = 8, titlefontsize = 10, framestyle = :box, grid = true)

# ── panel (a): closed-loop response, baseline vs ContEst (same noise draw) ─────
Xb = closed_loop(θ_base); Xo = closed_loop(θstar)
tb = (0:size(Xb, 2)-1) .* SPEC.dt
p1 = plot(tb, sigtrace(Xb); label = lab_base, lw = 2, ls = :dash, color = :gray35,
          xlabel = SPEC.tunit, ylabel = SPEC.sig, title = "(a) closed-loop response")
plot!(p1, tb, sigtrace(Xo); label = lab_opt, lw = 2, color = :firebrick)
hline!(p1, [0.0]; color = :black, lw = 0.6, label = "")

labels = ["J_est", "J_det", "J_c", "J_tot"]
bw = [b.J_est, b.J_det, b.J_c, b.J_tot]
ow = [o_.J_est, o_.J_det, o_.J_c, o_.J_tot]
# manual dodged bars (no StatsPlots dependency)
xb = (1:4) .- 0.18; xo = (1:4) .+ 0.18
p2 = bar(xb, bw; bar_width = 0.34, label = lab_base, color = :gray60,
         xticks = (1:4, labels), ylabel = "cost", title = "(b) cost breakdown")
bar!(p2, xo, ow; bar_width = 0.34, label = lab_opt, color = :firebrick)

plt = plot(p1, p2; layout = (1, 2), size = (760, 300), left_margin = 4Plots.mm,
           bottom_margin = 4Plots.mm, dpi = 200)
out = joinpath("Paper_LaTeX", "figs", "contest_$(CASE).pdf")
savefig(plt, out)
println("saved figure: ", out)

# ── print numbers for the LaTeX tables ────────────────────────────────────────
@printf("\n[%s] θ* = [%s]\n", CASE, join((@sprintf("%.3f", v) for v in θstar), ", "))
red(a, c) = 100 * (a - c) / a
@printf("  J_est : %10.4f -> %10.4f  (%.2f%%)\n", b.J_est, o_.J_est, red(b.J_est, o_.J_est))
@printf("  J_det : %10.4f -> %10.4f  (%.2f%%)\n", b.J_det, o_.J_det, red(b.J_det, o_.J_det))
@printf("  J_c   : %10.4f -> %10.4f  (%.2f%%)\n", b.J_c,   o_.J_c,   red(b.J_c,   o_.J_c))
@printf("  J_des : %10.4f -> %10.4f\n", b.J_des, o_.J_des)
@printf("  J_tot : %10.4f -> %10.4f  (%.2f%%)\n", b.J_tot, o_.J_tot, red(b.J_tot, o_.J_tot))
