# experiment_pll_fault_amplified.jl — EXPLORATORY, NOT part of the paper.
#
# Investigates a MORE SEVERE bus-1 fault for the PLL/DSE example: instead of
# the paper's fault (Σ_0(1,1)=3e-2, W(1,1)=2e-5), we scale ONLY the (1,1)
# entries by 100x:
#     Σ_0(1,1) : 3.0e-2 -> 3.0        (bus-1 angle prior variance)
#     W(1,1)   : 2.0e-5 -> 2.0e-3     (bus-1 angle process noise)
# (all other entries of Σ_0 and W, and the post-fault kick x_ic, are exactly
# as in example_pll.jl — unchanged).
#
# Does NOT modify example_pll.jl, experiment_pll_before_after_fault.jl, or any
# src/*.jl file. Everything new lives in this script: a second ControlModel
# `prob_fault` that reuses f_pll/h_pll/det_cost/stoch_cost/constraint_fn from
# example_pll.jl unchanged, with only Covars replaced.
#
# What this produces:
#   1. A Table-4-style parameter table: θ_nom, θ_peace (peacetime, unaffected
#      by the fault severity), θ_star_amp (post-fault design re-solved under
#      the amplified-noise fault, started from θ_peace) -- plus the associated
#      J_est/J_det/J_c/J_tot cost breakdown (M=100 samples, T=30 steps, same
#      convention as example_pll.jl's cost_breakdown).
#   2. A 5-second (T=100 steps @ dt=0.05s) time-domain figure: 100 Monte-Carlo
#      trajectories of the TRUE plant under the amplified post-fault scenario
#      (ICs drawn from the new Σ_0, process noise from the new W), for the
#      PEACE design (θ_peace) vs. the NEW post-fault design (θ_star_amp),
#      showing sample means and 2σ bands for Δδ_1, Δω_1, Δωᵖ_1 (bus 1, the
#      faulted bus).
#
# Run:  julia --project benchmarks/experiment_pll_fault_amplified.jl

ENV["GKSwstype"] = "100"   # headless GR

include("example_pll.jl")   # prob, mpc_eval, f_pll, h_pll, Q, R, det_cost,
                            # stoch_cost, constraint_fn, N_hor, x_ic, Σ_ic,
                            # Σ_w, σ_δ2, σ_ω2, κ_pll, θ_nom, θ_lb, θ_ub,
                            # J_des, ∇J_des, u_lin, M_mc, T_sim, θ_names
using Plots, Random, Printf, LaTeXStrings

# ── amplified fault: scale ONLY the (1,1) entries by 100x ─────────────────────
const Σ_ic_amp = copy(Σ_ic); Σ_ic_amp[1, 1] = 100 * Σ_ic[1, 1]      # 3e-2 -> 3.0
const Σ_w_amp  = copy(Σ_w);  Σ_w_amp[1, 1]  = 100 * Σ_w[1, 1]       # 2e-5 -> 2e-3

@printf("\nΣ_0(1,1): %.4g -> %.4g\n", Σ_ic[1,1], Σ_ic_amp[1,1])
@printf("W(1,1)  : %.4g -> %.4g\n\n", Σ_w[1,1], Σ_w_amp[1,1])

# same θ-dependent V(θ) as example_pll.jl, but process noise now uses Σ_w_amp
Covars_θ_amp(θ) = (Σ_w_amp,
                   Matrix(Diagonal(vcat(fill(σ_δ2, NB),
                                        [σ_ω2 * (1 + κ_pll * θ[NB+i]) for i in 1:NB]))))

const prob_fault = ControlModel(f_pll, h_pll, Covars_θ_amp, 3NB, NB, 2NB, N_hor,
                                det_cost, stoch_cost, constraint_fn)
const mpc_eval_fault = nonlinear_mpc_θ(prob_fault)

const contest_f_fault, contest_g_fault! =
    contest_objective(mpc_eval_fault, x_ic, Σ_ic_amp, u_lin, J_des, ∇J_des)

ctrl_info_fault(xe, Σe, θ) = mpc_eval_fault(xe, Σe, u_lin, θ; grad = false)[1]
function cost_breakdown_fault(θ)
    r = simulate_mc(prob_fault, ctrl_info_fault, θ, Matrix(Q), x_ic, Σ_ic_amp;
                    n_samples = M_mc, N_sim = T_sim)
    Jc = r.J_det + r.J_est
    return (; J_est = r.J_est, J_det = r.J_det, J_c = Jc, J_des = J_des(θ), J_tot = Jc + J_des(θ))
end

# ── TRIGGER 1 (peacetime, unaffected by the fault) — recompute θ_peace ────────
# Same scenario as experiment_pll_before_after_fault.jl's trigger 1 (system at
# rest, uniform quiet-level prior); the fault only concerns the post-fault
# scenario, so peacetime design is unchanged from the paper.
const x_ic_before = zeros(3NB)
const Σ_ic_before = Diagonal(fill(6.0e-3, 3NB)) |> Matrix
const contest_f_before, contest_g_before! =
    contest_objective(mpc_eval, x_ic_before, Σ_ic_before, u_lin, J_des, ∇J_des)

println("="^84)
println("  TRIGGER 1: peacetime design θ_peace (unaffected by fault severity)")
println("="^84)
_, θ_peace, _ = multistart_design(contest_f_before, contest_g_before!, θ_lb, θ_ub;
                                  n_starts = 5, θ_nom = θ_nom, seed = 20240624,
                                  g_tol = 1e-5, iterations = 100)
@printf("  θ_peace = [%s]\n", join((@sprintf("%.3f", v) for v in θ_peace), ", "))

# ── TRIGGER 2 (post-fault, AMPLIFIED noise): re-solve started from θ_peace ────
println("\n" * "="^84)
println("  TRIGGER 2: post-fault design θ_star_amp under 100x-amplified fault")
println("  (Σ_0(1,1) = $(Σ_ic_amp[1,1]), W(1,1) = $(Σ_w_amp[1,1]))")
println("="^84)
∇a, ∇fd, rel = verify_gradient(contest_f_fault, contest_g_fault!, θ_peace)
@printf("  gradient check at θ_peace: rel.err = %.2e  %s\n", rel, rel < 5e-2 ? "(OK)" : "(WARNING)")

_, θ_star_amp, results = multistart_design(contest_f_fault, contest_g_fault!, θ_lb, θ_ub;
                                           n_starts = 5, θ_nom = θ_peace, seed = 20240624,
                                           g_tol = 1e-5, iterations = 100)
let Js = sort([r.J for r in results]), nd = 1
    for k in 2:length(Js); Js[k] - Js[k-1] > 1e-4 && (nd += 1); end
    @printf("  starts: 5   distinct minima: %d\n", nd)
end
@printf("  θ_star_amp = [%s]\n", join((@sprintf("%.3f", v) for v in θ_star_amp), ", "))

# ── Table-4-style parameter table ──────────────────────────────────────────────
println("\n" * "─"^74)
println("  Table (like Table 4): θ_nom, θ_peace, θ_star_amp  [100x-amplified fault]")
println("─"^74)
@printf("  %-8s %-24s %10s %14s %14s\n", "param", "axis", "θ_nom", "θ_peace", "θ_star_amp")
for i in eachindex(θ_nom)
    @printf("  %-8s %-24s %10.3f %14.3f %14.3f\n",
            θ_names[i], θ_roles[i], θ_nom[i], θ_peace[i], θ_star_amp[i])
end

# ── cost breakdown: θ_peace vs θ_star_amp, both evaluated under the amplified fault
println("\n" * "─"^74)
println("  Cost breakdown under the amplified fault (M=$(M_mc), T=$(T_sim) steps)")
println("─"^74)
b_peace = cost_breakdown_fault(θ_peace)
b_star  = cost_breakdown_fault(θ_star_amp)
pct(b, o) = @sprintf("%+6.2f%%", 100*(o-b)/b)
@printf("  %-10s %14s %14s %10s\n", "component", "θ_peace", "θ_star_amp", "Δ%")
for (nm, b, o) in (("J_est", b_peace.J_est, b_star.J_est), ("J_det", b_peace.J_det, b_star.J_det),
                   ("J_c",   b_peace.J_c,   b_star.J_c),
                   ("J_tot", b_peace.J_tot, b_star.J_tot))
    @printf("  %-10s %14.4f %14.4f  %s\n", nm, b, o, pct(b,o))
end

# also, for reference: naive θ_nom under the same amplified fault
b_naive = cost_breakdown_fault(θ_nom)
println("\n  (for reference) θ_nom under the same amplified fault:")
@printf("  %-10s %14.4f\n  %-10s %14.4f\n  %-10s %14.4f\n  %-10s %14.4f\n",
        "J_est", b_naive.J_est, "J_det", b_naive.J_det, "J_c", b_naive.J_c, "J_tot", b_naive.J_tot)

# ── 5-second time-domain Monte Carlo rollout, mean +/- 2σ ─────────────────────
# Custom rollout (simulate_mc only returns aggregate costs): draws n_samples
# true ICs from Σ_ic_amp, simulates the true nonlinear plant with process
# noise from Σ_w_amp and the θ-dependent V(θ), runs the eKF + info-state MPC
# for N_sim steps, and records the full TRUE-state trajectory.
function rollout_trajectories(prob::ControlModel, controller::Function, θ, x_ic, Σ_ic;
                              n_samples = 100, N_sim = 100, seed = 20240624)
    Random.seed!(seed)
    n = prob.n
    W, V = resolve_covars(prob.Covars, θ)
    Lw  = cholesky(Symmetric(Matrix(W)) + 1e-12I).L
    Lv  = cholesky(Symmetric(Matrix(V)) + 1e-12I).L
    Lic = cholesky(Symmetric(Matrix(Σ_ic)) + 1e-12I).L
    X = zeros(n_samples, N_sim + 1, n)
    for s in 1:n_samples
        x_true = x_ic .+ Lic * randn(n)
        x_est  = copy(x_ic)
        Σ_est  = Matrix(Σ_ic)
        X[s, 1, :] = x_true
        for t in 1:N_sim
            u = controller(x_est, Σ_est, θ)
            x_true = prob.f(x_true, u, θ) .+ Lw * randn(n)
            y = prob.h(x_true, u, θ) .+ Lv * randn(prob.o)
            (x_est, Σ_est) = update((x_est, Σ_est), u, y, θ, prob.f, prob.h, (W, V); mode = "update")
            Σ_est = (Σ_est + Σ_est') / 2
            X[s, t + 1, :] = x_true
        end
    end
    return X
end

const N_5s = round(Int, 5.0 / dt_s)   # 100 steps @ dt_s = 0.05s
const n_traj = 100
println("\nRunning $(n_traj) x 5s ($(N_5s)-step) Monte-Carlo rollouts under the amplified " *
        "post-fault scenario for the peace design and the new post-fault design ...")

X_peace = rollout_trajectories(prob_fault, ctrl_info_fault, θ_peace,    x_ic, Σ_ic_amp; n_samples = n_traj, N_sim = N_5s, seed = 20240624)
X_star  = rollout_trajectories(prob_fault, ctrl_info_fault, θ_star_amp, x_ic, Σ_ic_amp; n_samples = n_traj, N_sim = N_5s, seed = 20240624)

tgrid = (0:N_5s) .* dt_s
mean_series(S) = vec(sum(S; dims = 1)) ./ size(S, 1)          # S: n_samples x (N_sim+1)
std_series(S)  = vec(sqrt.(sum((S .- reshape(mean_series(S), 1, :)).^2; dims = 1) ./ (size(S,1) - 1)))

default(fontfamily = "Computer Modern", guidefontsize = 10,
        tickfontsize = 8, framestyle = :box, grid = true, legend = false)

ω1_idx, ωp1_idx = NB + 1, 2NB + 1

# panel data: Δω_1 (true bus-1 frequency) and the PLL's estimation error Δω_1 - Δω̂_1^p
series_peace = (X_peace[:, :, ω1_idx], X_peace[:, :, ω1_idx] .- X_peace[:, :, ωp1_idx])
series_star  = (X_star[:, :, ω1_idx],  X_star[:, :, ω1_idx]  .- X_star[:, :, ωp1_idx])
labels = (L"$\Delta\omega_1$  [pu]", L"$\Delta\omega_1-\Delta\hat\omega_1^{\mathrm{p}}$  [pu]")

panels = []
for i in 1:2
    mP, sP = mean_series(series_peace[i]), std_series(series_peace[i])
    mS, sS = mean_series(series_star[i]),  std_series(series_star[i])
    is_bottom = i == 2
    p = plot(tgrid, mP; ribbon = 2 .* sP, ls = :solid,
             color = :forestgreen, lw = 2, fillalpha = 0.2,
             xlabel = is_bottom ? L"$t$  [s]" : "", ylabel = labels[i],
             xticks = is_bottom ? :auto : (0:5, fill("", 6)))
    plot!(p, tgrid, mS; ribbon = 2 .* sS, ls = :dash,
          color = :red, lw = 2, fillalpha = 0.2)
    hline!(p, [0.0]; color = :black, lw = 0.6)
    push!(panels, p)
end

plt = plot(panels...; layout = (2, 1), size = (620, 560), link = :x,
          left_margin = 4Plots.mm, bottom_margin = 0Plots.mm, top_margin = 0Plots.mm,
          dpi = 200)
outdir = joinpath("ContEst_TeX", "figs")
isdir(outdir) || mkpath(outdir)
out = joinpath(outdir, "pll_fault_amplified_timedomain.pdf")
savefig(plt, out)
println("saved figure: ", out)

out_png = joinpath(outdir, "pll_fault_amplified_timedomain.png")
savefig(plt, out_png)
println("saved figure: ", out_png)
