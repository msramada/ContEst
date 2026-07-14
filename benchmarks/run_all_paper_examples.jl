# run_all_paper_examples.jl — run every example used in the ContEst paper and
# print a consolidated, paper-style results table.
#
# The four studies in the manuscript are:
#   • example_adcs.jl         — spacecraft ADCS, shared power/mass budget   (eKF–MPC, H₂)
#   • example_distillation.jl — feed-stage / sensor-tray placement          (eKF–MPC, H₂)
#   • example_pll.jl          — low-inertia grid PLL/DSE co-design          (eKF–MPC, H₂)
#   • example_mtdc.jl         — multi-terminal HVDC droop coordination      (output-fb H∞)
#
# Each example file is `include`d into its OWN module so their many top-level
# `const`s (θ_nom, Q, θ_lb, …) do not clash. Their `if abspath(PROGRAM_FILE)==…`
# guards keep the per-example reports from auto-running, so here we drive the
# multi-start BFGS ourselves and collect the baseline→optimum cost breakdown.
#
# The final two tables mirror the paper's Table (consolidated formulations) and
# Table (summary of cost reductions).
#
# Run:  julia --project benchmarks/run_all_paper_examples.jl

using Printf

const HERE = @__DIR__

# ── Load each example into an isolated module ─────────────────────────────────
# `include` inside the module resolves the example's own relative `include`s
# (../src/…) against the example file, so each module gets its own machinery.
function load_example(modname::Symbol, file::String)
    @eval Main module $modname
        include(joinpath($HERE, $file))
    end
    # `invokelatest` avoids the Julia ≥1.12 world-age warning for a binding
    # accessed in the same top-level call that defined it.
    return Base.invokelatest(getfield, Main, modname)
end

println("Loading examples (this compiles the eKF–MPC and H∞ machinery)…")
ADCS   = load_example(:Ex_adcs,   "example_adcs.jl")
DISTIL = load_example(:Ex_distil, "example_distillation.jl")
PLL    = load_example(:Ex_pll,    "example_pll.jl")
MTDC   = load_example(:Ex_mtdc,   "example_mtdc.jl")

# ── Common optimiser settings (identical to the per-example reports) ──────────
const SEED     = 20240624
const N_STARTS = 5

# A uniform result record for the summary tables.
struct StudyResult
    name   :: String
    dims   :: NTuple{3,Int}   # (r_x, r_u, r_y)
    inner  :: String          # inner-value description
    theta  :: String          # design-vector description (n_f + n_h, where it enters)
    base   :: NamedTuple      # baseline cost breakdown
    opt    :: NamedTuple      # optimum  cost breakdown
end

# ── eKF–MPC studies (adcs / distillation / pll) ───────────────────────────────
# Reuse the module's own contest objective + multi-start BFGS, then read the
# sampled closed-loop cost breakdown (role-aware where the example supports it).
function run_ekf_study(mod, name, inner, theta)
    @printf("  ▸ %-28s optimising (multi-start BFGS, %d starts)…\n", name, N_STARTS)
    _, best_θ, _ = mod.multistart_design(mod.contest_f, mod.contest_g!,
                                         mod.θ_lb, mod.θ_ub;
                                         n_starts = N_STARTS, θ_nom = mod.θ_nom,
                                         seed = SEED, g_tol = 1e-5, iterations = 100)
    cb = mod.cost_breakdown
    base = applicable(cb, mod.θ_nom, :baseline) ? cb(mod.θ_nom, :baseline) : cb(mod.θ_nom)
    opt  = applicable(cb, best_θ,   :optimum)   ? cb(best_θ,   :optimum)   : cb(best_θ)
    p = mod.prob
    return StudyResult(name, (p.n, p.m, p.o), inner, theta, base, opt)
end

# ── MTDC study (output-feedback H∞) ───────────────────────────────────────────
# The H∞ path has no closed-loop sampler; the worst-case guaranteed costs come
# straight from the two game Riccatis via `of_eval`. Rebuild the same objective
# `hinf_of_report` optimises and read the baseline/optimum breakdown.
function run_mtdc_study(mod, name, inner, theta)
    @printf("  ▸ %-28s optimising (multi-start BFGS, %d starts)…\n", name, N_STARTS)
    eval_J(θ) = begin
        Jc, ∇Jc, _, _, _, _ = mod.of_eval(mod.model, θ)
        (Jc + mod.J_des(θ), ∇Jc .+ mod.∇J_des(θ))
    end
    f, g! = mod.cached_objective(eval_J)
    _, best_θ, _ = mod.multistart_design(f, g!, mod.θ_lb, mod.θ_ub;
                                         n_starts = N_STARTS, θ_nom = mod.θ_nom,
                                         seed = SEED, g_tol = 1e-5, iterations = 100)
    breakdown(θ) = begin
        _, _, _, Jd, Je, _ = mod.of_eval(mod.model, θ)
        (; J_est = Je, J_det = Jd, J_c = Jd + Je,
           J_des = mod.J_des(θ), J_tot = Jd + Je + mod.J_des(θ))
    end
    base = breakdown(mod.θ_nom); opt = breakdown(best_θ)
    r_x = 2 * mod.NT
    return StudyResult(name, (r_x, mod.NT, mod.NT), inner, theta, base, opt)
end

println("\nRunning all four studies…")
results = StudyResult[
    run_ekf_study(ADCS,   "ADCS (spacecraft)", "eKF–MPC (H₂)",       "1+2: e_rw∈f, α∈V"),
    run_ekf_study(DISTIL, "Distillation",      "eKF–MPC (H₂)",       "3+3: feeds∈f, sensors∈h"),
    run_ekf_study(PLL,    "PLL/DSE (grid)",    "eKF–MPC (H₂)",       "3+3: e∈f, b∈f,V"),
    run_mtdc_study(MTDC,  "MTDC (HVDC)",       "H∞ output-fb.",      "25+0: k∈f (droop)"),
]

# ── Reporting helpers ─────────────────────────────────────────────────────────
red(b, o) = b == 0 ? 0.0 : 100 * (b - o) / b               # % reduction (larger is better)
bar = "=" ^ 92

# ── Table 1: consolidated view of the studies (mirrors tab:formulations) ──────
println("\n\n", bar)
println("  TABLE 1 — Consolidated view of the studies (dimensions, design vector, reductions)")
println(bar)
@printf("  %-20s %-11s %-24s %-16s %8s %8s %8s\n",
        "study", "(rx,ru,ry)", "θ: nf+nh (enters)", "inner value",
        "Jest↓", "Jdet↓", "Jtot↓")
println("  " * "-"^88)
for r in results
    @printf("  %-20s %-11s %-24s %-16s %7.1f%% %7.1f%% %7.1f%%\n",
            r.name, "($(r.dims[1]),$(r.dims[2]),$(r.dims[3]))", r.theta, r.inner,
            red(r.base.J_est, r.opt.J_est),
            red(r.base.J_det, r.opt.J_det),
            red(r.base.J_tot, r.opt.J_tot))
end
println(bar)

# ── Table 2: summary of cost reductions, baseline → optimum (mirrors tab:summary)
println("\n", bar)
println("  TABLE 2 — Summary of cost reductions (baseline θ_nom → ContEst optimum θ*)")
println("  J_est: estimation · J_det: control · J_c=J_det+J_est · J_tot=J_c+J_des")
println(bar)
@printf("  %-20s | %-18s %6s | %-18s %6s | %-18s %6s | %7s | %-20s %6s\n",
        "study", "J_est base→opt", "red.", "J_det base→opt", "red.",
        "J_c base→opt", "red.", "J_des", "J_tot base→opt", "red.")
println("  " * "-"^140)
for r in results
    b, o = r.base, r.opt
    @printf("  %-20s | %8.1f→%-8.1f %5.1f%% | %8.1f→%-8.1f %5.1f%% | %8.1f→%-8.1f %5.1f%% | %7.1f | %8.1f→%-8.1f %5.1f%%\n",
            r.name,
            b.J_est, o.J_est, red(b.J_est, o.J_est),
            b.J_det, o.J_det, red(b.J_det, o.J_det),
            b.J_c,   o.J_c,   red(b.J_c,   o.J_c),
            o.J_des,
            b.J_tot, o.J_tot, red(b.J_tot, o.J_tot))
end
println(bar)
println("\nNote: eKF–MPC costs are sampled closed-loop (simulate_mc, seed=$SEED);")
println("MTDC costs are the exact H∞ worst-case guaranteed values. Numbers match the")
println("per-example reports and the paper's summary tables.")
