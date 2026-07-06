# report_contest.jl – shared reporting for the ContEst candidate examples.
# Produces the GUIDE.md §5–§6 output: gradient check, multi-start BFGS summary
# (best minimum only), cost-component table, and design-parameter table.
#
# Requires src/BFGS.jl to be already included (verify_gradient, multistart_design).

using Printf, LinearAlgebra

function run_contest_report(title, θ_names, θ_roles;
                            contest_f, contest_g!, θ_init, θ_nom, θ_lb, θ_ub,
                            cost_breakdown, seed = 20240624, n_starts = 5)
    nθ = length(θ_init)
    bar = "=" ^ 74
    println(bar); println("  ContEst  —  ", title); println(bar)

    # ── Gradient verification ────────────────────────────────────────────────
    ∇a, ∇fd, rel = verify_gradient(contest_f, contest_g!, θ_init)
    println("\n── Gradient verification at θ_init ─────────────────────────")
    @printf("  analytic : [%s]\n", join((@sprintf("%8.4f", v) for v in ∇a), " "))
    @printf("  fin.diff : [%s]\n", join((@sprintf("%8.4f", v) for v in ∇fd), " "))
    @printf("  relative error: %.2e  %s\n", rel,
            rel < 5e-2 ? "(OK)" : "(WARNING: large discrepancy)")

    # ── Multi-start BFGS (baseline + random restarts) ────────────────────────
    println("\n── Multi-start BFGS ($n_starts starts: θ_nom + $(n_starts-1) random) ──")
    best_J, best_θ, results = multistart_design(contest_f, contest_g!, θ_lb, θ_ub;
                                                n_starts = n_starts, θ_nom = θ_nom,
                                                seed = seed, g_tol = 1e-5, iterations = 100)
    # count distinct minima (by objective value, 1e-4 tol)
    Js = sort([r.J for r in results])
    ndistinct = 1
    for k in 2:length(Js)
        Js[k] - Js[k-1] > 1e-4 && (ndistinct += 1)
    end
    @printf("  starts run: %d   distinct minima found: %d   best J_tot = %.4f\n",
            n_starts, ndistinct, best_J)
    println(ndistinct == 1 ?
            "  → all starts agree: optimum is (numerically) global." :
            "  → multimodal: reporting the best minimum only (others omitted).")

    # ── Cost decomposition: baseline vs best ─────────────────────────────────
    # eKF–MPC examples supply a role-aware `cost_breakdown(θ, role)` that reports
    # SAMPLED closed-loop costs — the baseline (:baseline) is driven by a
    # certainty-equivalence controller (eKF mean only), the optimum (:optimum) by
    # the information-state MPC. Legacy 1-arg breakdowns still work.
    base = applicable(cost_breakdown, θ_nom, :baseline) ? cost_breakdown(θ_nom, :baseline) : cost_breakdown(θ_nom)
    opt  = applicable(cost_breakdown, best_θ, :optimum) ? cost_breakdown(best_θ, :optimum) : cost_breakdown(best_θ)
    pct(b, o_) = abs(b) < 1e-9 ? "    — " : @sprintf("%+6.2f%%", 100 * (o_ - b) / b)

    println("\n── Cost components: baseline θ_nom → optimal θ* ────────────")
    @printf("  %-26s %12s %12s %12s %8s\n", "component", "baseline", "optimal", "Δ (incr +)", "%")
    rows = (("Estimation cost  J_est",          base.J_est, opt.J_est),
            ("Deterministic ctrl J_det",        base.J_det, opt.J_det),
            ("Stochastic ctrl (both) J_c",      base.J_c,   opt.J_c),
            ("Design cost  J_des",              base.J_des, opt.J_des),
            ("TOTAL  J_tot",                    base.J_tot, opt.J_tot))
    for (name, b, o_) in rows
        @printf("  %-26s %12.4f %12.4f %12.4f  %s\n", name, b, o_, o_ - b, pct(b, o_))
    end
    @printf("\n  Net total reduction: %.4f  (%.2f%%)\n",
            base.J_tot - opt.J_tot, 100 * (base.J_tot - opt.J_tot) / base.J_tot)

    # ── Design-parameter table (best minimum only) ───────────────────────────
    println("\n── Design parameters: baseline → optimal (best minimum) ────")
    @printf("  %-10s %10s %10s %10s   %s\n", "param", "baseline", "optimal", "Δ", "enters")
    for i in 1:nθ
        @printf("  %-10s %10.3f %10.3f %+10.3f   %s\n",
                θ_names[i], θ_nom[i], best_θ[i], best_θ[i] - θ_nom[i], θ_roles[i])
    end
    println(bar)
    return best_θ, best_J, base, opt
end
