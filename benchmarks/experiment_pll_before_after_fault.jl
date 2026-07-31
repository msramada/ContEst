# experiment_pll_before_after_fault.jl — reproduces the two-stage
# (peacetime / post-fault) ContEst results now cited in
# ContEst_TeX/main.tex, Section~\ref{sec:sim:pll} and Table~\ref{tab:pll}
# (theta_peace, theta_star, and the trigger-1/trigger-2 percentages). Kept
# in the repo (not deleted) so those paper numbers stay reproducible, per the
# project convention of backing every reported number with runnable code.
# Does not modify example_pll.jl / src/*.jl in any way (read-only include).
#
# Two ContEst applications, same dynamics/cost/box/reporting convention as
# example_pll.jl throughout:
#   TRIGGER 1 (peacetime) — system at rest, no disturbance: x_ic=0
#     (equilibrium), Σ_ic UNIFORM at the paper's quiet-bus level (no bus-1
#     uncertainty boost). Produces θ_peace.
#   TRIGGER 2 (post-fault) — exactly the paper's existing fault scenario
#     (x_ic, Σ_ic as coded in example_pll.jl), but re-solved STARTING FROM
#     θ_peace instead of θ_nom. Produces θ*.
#
# Run:  julia --project benchmarks/experiment_pll_before_after_fault.jl

include("example_pll.jl")   # prob, mpc_eval, θ_nom, θ_lb, θ_ub, J_des, ∇J_des,
                            # Q, M_mc, T_sim, θ_names, u_lin, x_ic, Σ_ic,
                            # contest_f, contest_g!, cost_breakdown (all the
                            # paper's existing AFTER-fault/trigger-2 objects,
                            # unmodified; its own report is guarded, so
                            # including it here runs nothing beyond
                            # definitions)
using Printf

# ── TRIGGER 1 scenario: system at rest, no disturbance ────────────────────
const x_ic_before = zeros(3NB)
const Σ_ic_before = Diagonal(fill(6.0e-3, 3NB)) |> Matrix   # uniform, quiet-bus level

const contest_f_before, contest_g_before! =
    contest_objective(mpc_eval, x_ic_before, Σ_ic_before, u_lin, J_des, ∇J_des)

ctrl_info_before(xe, Σe, θ) = mpc_eval(xe, Σe, u_lin, θ; grad = false)[1]
function cost_breakdown_before(θ)
    r = simulate_mc(prob, ctrl_info_before, θ, Matrix(Q), x_ic_before, Σ_ic_before;
                    n_samples = M_mc, N_sim = T_sim)
    Jc = r.J_det + r.J_est
    return (; J_est = r.J_est, J_det = r.J_det, J_c = Jc, J_des = J_des(θ), J_tot = Jc + J_des(θ))
end

function run_trigger(name, contest_f_, contest_g_!, cost_bd, θ_start; seed = 20240624, n_starts = 5)
    println("\n" * "─"^74)
    println("  $name")
    println("─"^74)
    ∇a, ∇fd, rel = verify_gradient(contest_f_, contest_g_!, θ_start)
    @printf("  gradient check at θ_start: rel.err = %.2e  %s\n", rel, rel < 5e-2 ? "(OK)" : "(WARNING)")

    best_J, best_θ, results = multistart_design(contest_f_, contest_g_!, θ_lb, θ_ub;
                                                n_starts = n_starts, θ_nom = θ_start,
                                                seed = seed, g_tol = 1e-5, iterations = 100)
    Js = sort([r.J for r in results]); nd = 1
    for k in 2:length(Js); Js[k]-Js[k-1] > 1e-4 && (nd += 1); end
    @printf("  starts: %d   distinct minima: %d   best surrogate J = %.5f\n", n_starts, nd, best_J)

    base = cost_bd(θ_start); opt = cost_bd(best_θ)
    pct(b, o) = @sprintf("%+6.2f%%", 100*(o-b)/b)
    @printf("  %-10s %12s %12s %10s\n", "component", "start θ", "optimal", "Δ%")
    for (nm, b, o) in (("J_est", base.J_est, opt.J_est), ("J_det", base.J_det, opt.J_det),
                       ("J_c",   base.J_c,   opt.J_c),
                       ("J_tot", base.J_tot, opt.J_tot))
        @printf("  %-10s %12.4f %12.4f  %s\n", nm, b, o, pct(b,o))
    end
    return (; best_θ, base, opt)
end

bar = "="^84
println(bar)
println("  ContEst TRIGGER 1 (peacetime, started from θ_nom) vs.")
println("  TRIGGER 2 (post-fault, started from θ_peace) — PLL/DSE co-design")
println(bar)

trig1 = run_trigger("TRIGGER 1: peacetime (x_ic=0, uniform quiet-level Σ_0), start=θ_nom",
                    contest_f_before, contest_g_before!, cost_breakdown_before, θ_nom)
θ_peace = trig1.best_θ

trig2 = run_trigger("TRIGGER 2: post-fault (paper's scenario, unmodified), start=θ_peace",
                    contest_f, contest_g!, θ -> cost_breakdown(θ, :x), θ_peace)
θ_star = trig2.best_θ

println("\n" * bar)
println("  Design comparison: θ_nom -> θ_peace -> θ_star")
println(bar)
@printf("  %-10s %10s %14s %14s\n", "param", "θ_nom", "θ_peace", "θ_star")
for i in eachindex(θ_nom)
    @printf("  %-10s %10.3f %14.3f %14.3f\n", θ_names[i], θ_nom[i], θ_peace[i], θ_star[i])
end

println("\n── Combined effect vs. the original naive baseline θ_nom (both triggers) ──")
base_after = cost_breakdown(θ_nom, :x)
pct(b, o) = @sprintf("%+6.2f%%", 100*(o-b)/b)
@printf("  %-10s %12s %12s %10s\n", "component", "θ_nom", "θ_star", "Δ%")
for (nm, b, o) in (("J_est", base_after.J_est, trig2.opt.J_est), ("J_det", base_after.J_det, trig2.opt.J_det),
                   ("J_c",   base_after.J_c,   trig2.opt.J_c),
                   ("J_tot", base_after.J_tot, trig2.opt.J_tot))
    @printf("  %-10s %12.4f %12.4f  %s\n", nm, b, o, pct(b,o))
end
println(bar)
