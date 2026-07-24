# example_mtdc_sparse.jl — OFFICIAL structured (block-diagonal) MTDC co-design study.
#
# The device network is a DC RING: terminal i couples only to its two neighbours, so
# A(θ) is block-tridiagonal-plus-corners.  BOTH inner values are posed as structured
# semidefinite programs whose per-terminal 2×2 block-diagonal (decentralised) restriction
# cuts the O(n²) variables to O(n) and keeps the dual strictly complementary, so the
# LMI dual returns the EXACT envelope gradient with no solver differentiation:
#
#   CONTROL  H∞ bounded-real-lemma SDP in Y=P⁻¹, L=KY, at fixed μ=γ²=16:
#            J_det = min tr(Y⁻¹W) s.t. BRL LMI;  ∂J_det/∂A = 2 S₃₁ Y  (dense = game Riccati).
#   SENSING  H₂ estimation SDP in the P^ε form (example_mtdc.jl:est_sdp):
#            J_est = min tr(WP)+tr(VZ);  ∂J_est/∂A = −2 P S₁₂ᵀ,  ∂J_est/∂C = 2 Fᵀ S₁₂ᵀ.
#
# The droop θ enters A(θ) and so reshapes BOTH values → the two axes couple and the
# co-design does not separate.  The sensor gains α∈[0,1] scale the rows of C(α) and are
# optimised jointly with θ (outer BFGS on the envelope gradient) under an ℓ1 penalty that
# selects the sensor set.  Every reported cost is the TRUE value at the design:
# J_det = tr(XW) (game Riccati) and J_est = tr(Q Σ^ε_true) (Lyapunov of the SDP gain).
#
# The figure (contest_mtdc.pdf) contrasts, over the number of retained sensors:
#   (A) sensor pruning at the fixed baseline droop θ_nom  vs
#   (B) joint (α,θ) co-design — showing the estimation cost the co-design buys.
#
# Run:  julia --project benchmarks/example_mtdc_sparse.jl
#
include("example_mtdc.jl")     # model, consts, idx, Cmat/Vmat, est_sdp, est_true_cost,
                               # kalman_cost, est_val_grad, Jdet_ric, hinf_gare, dlyap …
ENV["GKSwstype"] = "100"       # headless GR
using JuMP, Clarabel, LinearAlgebra, ForwardDiff, Printf, Plots, LaTeXStrings

const GAMMA2 = γ²              # 16.0 — exact block-diagonal envelope gradient

# Noise-tolerant box optimizer: the block-diagonal control SDP value carries ~1e-5
# numerical noise, which makes the default HagerZhang line search thrash.  A
# BackTracking line search with a capped iteration budget is robust and cheap; the
# frontier only needs good (not razor-tight) designs.  Reused by droop/joint opts.
const LS = Optim.LineSearches.BackTracking(order = 2)
function opt_box(f, g!, lb, ub, z0; iters = 15, outer = 3)
    z0i = clamp.(z0, lb .+ 1e-6, ub .- 1e-6)      # Fminbox needs a strictly interior start
    r = optimize(f, g!, lb, ub, z0i, Fminbox(BFGS(linesearch = LS)),
                 Optim.Options(iterations = iters, outer_iterations = outer,
                               g_tol = 1e-3, f_reltol = 1e-4))
    r.minimizer, Optim.minimum(r)
end

# ── H∞ control SDP (bounded-real lemma), pattern ∈ (:full, :blockdiag, :banded) ──
function hinf_sdp(θ; pattern = :blockdiag, μ = GAMMA2, withgrad = true)
    A, B, W = model(θ); nn = size(A, 1); mm = size(B, 2)
    Wm = Matrix(W); Wh = Matrix(sqrt(Symmetric(Wm))); E = Wh
    Qh = Matrix(sqrt(Symmetric(Matrix(Qmat)))); Rh = Matrix(sqrt(Symmetric(Matrix(Rmat))))
    nw = nn; nz = nn + mm
    mdl = Model(Clarabel.Optimizer); set_silent(mdl)
    @variable(mdl, T[1:nn, 1:nn], Symmetric)
    if pattern == :full
        Yv = @variable(mdl, [1:nn, 1:nn], Symmetric); Ym = Matrix(Yv)
        Lv = @variable(mdl, [1:mm, 1:nn]); Lm = Matrix(Lv)
    else
        Y = zeros(AffExpr, nn, nn); L = zeros(AffExpr, mm, nn)
        for i in 1:NT
            Yi = @variable(mdl, [1:2, 1:2], Symmetric); r = idx(i)
            for a in 1:2, b in 1:2; Y[r[a], r[b]] = Yi[a, b]; end
            nbrs = pattern == :blockdiag ? (i,) : (i, i == 1 ? NT : i - 1, i == NT ? 1 : i + 1)
            for j in nbrs
                Lij = @variable(mdl, [1:1, 1:2]); c = idx(j); L[i, c[1]] = Lij[1, 1]; L[i, c[2]] = Lij[1, 2]
            end
        end
        if pattern == :banded
            for i in 1:NT
                j = i == NT ? 1 : i + 1; Yij = @variable(mdl, [1:2, 1:2]); ri = idx(i); rj = idx(j)
                for a in 1:2, b in 1:2; Y[ri[a], rj[b]] = Yij[a, b]; Y[rj[b], ri[a]] = Yij[a, b]; end
            end
        end
        Ym = Matrix(Y); Lm = Matrix(L)
    end
    @constraint(mdl, Symmetric([Matrix(T) Wh; Wh Ym]) in PSDCone())
    AYBL = A * Ym .+ B * Lm; CzYDzL = [Qh * Ym; Rh * Lm]
    M = [ -Ym            zeros(nn, nw)          permutedims(AYBL)  permutedims(CzYDzL);
          zeros(nw, nn)  -μ * Matrix(I, nw, nw) permutedims(E)     zeros(nw, nz);
          AYBL           E                      -Ym                zeros(nn, nz);
          CzYDzL         zeros(nz, nw)          zeros(nz, nn)      -Matrix(I, nz, nz) ]
    @constraint(mdl, brl, Symmetric(-M) in PSDCone())
    @objective(mdl, Min, tr(Matrix(T)))
    optimize!(mdl)
    st = termination_status(mdl)
    (st in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.NUMERICAL_ERROR)) || return (1e6, zeros(nn, nn), st, solve_time(mdl), num_variables(mdl))
    J = objective_value(mdl)
    withgrad || return (J, zeros(nn, nn), st, solve_time(mdl), num_variables(mdl))
    S = dual(brl); b3 = (nn + nw + 1):(2nn + nw)
    (J, 2 .* (Matrix(S[b3, 1:nn]) * Matrix(value.(Ym))), st, solve_time(mdl), num_variables(mdl))
end

# ── True (reported) inner values at a design ─────────────────────────────────
# control: game-Riccati worst-case tr(XW) (Jdet_ric, example_mtdc.jl)
# estimation: tr(Q Σ^ε_true) of the block-diagonal SDP filter on the selected sensor set
function est_true_at(θ, sel; pattern = :blockdiag)
    A, _, W = model(θ); C = Cmat(Float64.(sel))
    J, G, _, _, st = est_sdp(A, C, W, Vmat, Qmat; pattern = pattern)
    est_true_cost(A, C, W, Vmat, Qmat, G)
end

# ── Value+gradient builders (block-diagonal SDPs, exact envelope gradients) ────
# droop-only (α fixed): J = J_cont(θ) + J_est(θ,α) + c_k Σθ,  grad over θ (through A).
function droop_vg(θ, α; pattern = :blockdiag)
    Jc, gA_c, _, _, _ = hinf_sdp(θ; pattern = pattern); Jc >= 1e6 && return (1e6, zeros(m))
    Je, _, gA_e, _, _ = est_val_grad(θ, α; pattern = pattern); Je >= 1e6 && return (1e6, zeros(m))
    gA = gA_c .+ gA_e
    Aj = ForwardDiff.jacobian(t -> vec(model(t)[1]), θ)
    (Jc + Je + c_k * sum(θ), [tr(gA' * reshape(Aj[:, k], n, n)) + c_k for k in 1:m])
end
droop_opt(α; pattern = :blockdiag, iters = 25) = begin
    f, g! = cached_objective(θ -> droop_vg(θ, α; pattern = pattern))
    θ, _ = opt_box(f, g!, θ_lb, θ_ub, copy(θ_nom); iters = iters)
    θ
end

# α-step (θ fixed): minimize J_est(θ,α) + λs Σα over α∈[0,1] — clean exact gradient
# ∂/∂α_i = [∂J_est/∂C]_{i,2i-1} + λs from the estimation SDP dual (no control-SDP noise).
function alpha_opt(θ; λs = λs_default, pattern = :blockdiag, α0 = ones(m))
    f, g! = cached_objective(function (α)
        Je, _, _, gC, _ = est_val_grad(θ, α; pattern = pattern)
        (Je + λs * sum(α), [gC[i, 2i-1] + λs for i in 1:m])
    end)
    α, _ = opt_box(f, g!, zeros(m), ones(m), copy(α0); iters = 25)
    α
end

# Joint (α,θ) co-design by ALTERNATING minimization of J_cont(θ)+J_est(θ,α)+λsΣα+c_kΣθ:
# the θ-step (droop_opt, block-diagonal H∞+H₂) and the α-step (alpha_opt) are each
# well-conditioned, so this is far more robust than one 30-var solve on the noisy joint.
function joint_opt(; λs = λs_default, pattern = :blockdiag, z0 = vcat(θ_nom, ones(m)),
                   rounds = 3, iters = 15)
    θ = z0[1:m]; α = z0[m+1:2m]
    for _ in 1:rounds
        α = alpha_opt(θ; λs = λs, pattern = pattern, α0 = α)
        θ = droop_opt(α; pattern = pattern, iters = iters)
    end
    vcat(θ, α)
end

# ── Curve A: greedy sensor pruning at the fixed baseline droop θ_nom ──────────
# Returns costs[k] = est cost with k sensors, and subset[k] = the retained k-set.
function prune_baseline(; pattern = :blockdiag)
    sel = trues(m); costs = fill(NaN, m); subset = Vector{Vector{Bool}}(undef, m)
    costs[m] = est_true_at(θ_nom, sel; pattern = pattern); subset[m] = collect(sel)
    for k in (m-1):-1:1
        best = Inf; bestj = 0
        for j in findall(sel)
            trial = copy(sel); trial[j] = false
            c = est_true_at(θ_nom, trial; pattern = pattern)
            (c < best) && (best = c; bestj = j)
        end
        sel[bestj] = false; costs[k] = best; subset[k] = collect(sel)
    end
    (costs, subset)                          # costs[k], subset[k] (k = #sensors)
end

# ── Curve B: joint co-design on the SAME per-count sensor sets — co-design the droop
# for each retained-sensor budget (identical sensors as curve A; baseline vs co-designed
# droop isolates the co-design gain). Records estimation and control cost per count.
function joint_frontier(subset; pattern = :blockdiag)
    estB = fill(NaN, m); contB = fill(NaN, m)
    for k in m:-1:1
        sel = subset[k]
        θk = droop_opt(Float64.(sel); pattern = pattern, iters = 20)
        estB[k] = est_true_at(θk, sel; pattern = pattern); contB[k] = Jdet_ric(θk)
        @printf("    k=%2d → J_est(co-design)=%.4f  J_cont=%.4f\n", k, estB[k], contB[k])
    end
    (estB, contB)
end

# ── Figure: two panels vs number of sensors — (a) estimation cost, (b) control
# cost — each contrasting baseline-droop pruning with joint α+θ co-design. Panel
# (b) shows the joint control cost stays low (near/below the fixed-droop baseline)
# while estimation improves: co-design buys sensing without sacrificing control.
function mtdc_figure(estA, estB, JcontA_const, contB;
                     path = joinpath("ContEst_TeX", "figs", "contest_mtdc.pdf"))
    default(fontfamily = "Computer Modern", legendfontsize = 12, guidefontsize = 15,
            tickfontsize = 13, titlefontsize = 16, framestyle = :box, grid = true)
    ks = 1:m
    fin(v) = [k for k in ks if isfinite(v[k])]
    kA = fin(estA); kB = fin(estB)
    prune_col = RGB(0.13, 0.55, 0.13)   # forest green
    joint_col = RGB(0.80, 0.0, 0.0)     # red
    # identical styling in both panels: green dashed + circles (baseline/pruning),
    # red solid + diamonds (joint). Legends removed; series described in the caption.
    ms_pts = 6; lw_ln = 3
    p1 = plot(kA, [estA[k] for k in kA]; color = prune_col, ls = :dash, lw = lw_ln, marker = :circle, ms = ms_pts,
              xlabel = "number of voltage sensors",
              ylabel = L"estimation cost  $\mathrm{tr}(Q\Sigma)$", title = "(a) sensing", legend = false)
    plot!(p1, kB, [estB[k] for k in kB]; color = joint_col, ls = :solid, lw = lw_ln, marker = :diamond, ms = ms_pts)
    kBc = fin(contB)
    p2 = plot(ks, fill(JcontA_const, m); color = prune_col, ls = :dash, lw = lw_ln, marker = :circle, ms = ms_pts,
              xlabel = "number of voltage sensors",
              ylabel = L"control cost  $\mathrm{tr}(XW)$", title = "(b) control", legend = false)
    plot!(p2, kBc, [contB[k] for k in kBc]; color = joint_col, ls = :solid, lw = lw_ln, marker = :diamond, ms = ms_pts)
    plt = plot(p1, p2; layout = (1, 2), size = (900, 380), dpi = 200,
               left_margin = 7Plots.mm, bottom_margin = 7Plots.mm)
    mkpath(dirname(path)); savefig(plt, path)
    @printf("  figure written: %s\n", path)
    plt
end

# ── Report ───────────────────────────────────────────────────────────────────
function mtdc_sparse_report(; λs_grid = 0.0:0.01:0.10)
    bar = "="^84
    println(bar)
    println("  MTDC — block-diagonal H∞ control + H₂ sensing co-design (γ²=$(GAMMA2), $(2NT) states)")
    println(bar)
    allon = ones(m)

    # ---- control: full-dense vs block-diagonal vs game-Riccati (at baseline) ----
    hinf_sdp(θ_nom; pattern = :blockdiag, withgrad = false)                      # warmup
    Jf, _, _, tf, nvf = hinf_sdp(θ_nom; pattern = :full)
    Jd, _, _, td, nvd = hinf_sdp(θ_nom; pattern = :blockdiag)
    Jric = Jdet_ric(θ_nom)
    println("\n── Control SDP: full dense vs block-diagonal (baseline θ=3, γ²=$(GAMMA2)) ──")
    @printf("  %-26s %10s %12s %12s\n", "formulation", "#variables", "solve time", "J_det=tr(XW)")
    @printf("  %-26s %10d %10.1f s %12.4f\n", "full dense (Y,L,T)", nvf, tf, Jf)
    @printf("  %-26s %10d %10.2f s %12.4f\n", "block-diagonal (official)", nvd, td, Jd)
    @printf("  game-Riccati check tr(XW)=%.4f (== full dense);  block-diag +%.1f%%, %.0f× fewer vars\n",
            Jric, 100*(Jd-Jf)/Jf, nvf/nvd)

    # ---- estimation validity: block-diagonal vs centralised Kalman (at baseline) ----
    A0, _, W0 = model(θ_nom); C0 = Cmat(allon); Jk = kalman_cost(A0, C0, W0, Vmat, Qmat)
    println("\n── Estimation SDP validity (baseline θ=3, all $NT sensors) ──")
    @printf("  %-14s %12s %14s %12s\n", "pattern", "SDP tr(WP)", "true tr(QΣ)", "vs Kalman")
    for pat in (:full, :blockdiag, :banded)
        Je, G, _, _, _ = est_sdp(A0, C0, W0, Vmat, Qmat; pattern = pat)
        @printf("  %-14s %12.4f %14.4f %+11.1f%%\n", pat, Je,
                est_true_cost(A0, C0, W0, Vmat, Qmat, G), 100*(est_true_cost(A0,C0,W0,Vmat,Qmat,G)-Jk)/Jk)
    end
    @printf("  centralised Kalman reference tr(QΣ)=%.4f\n", Jk)

    # ---- (A) droop co-design (all sensors) ----
    println("\n── (A) Droop co-design, block-diagonal driven (all $NT sensors) ──")
    θ1 = droop_opt(allon)
    Jd0, Je0 = Jdet_ric(θ_nom), est_true_at(θ_nom, trues(m))
    Jd1, Je1 = Jdet_ric(θ1),   est_true_at(θ1, trues(m))
    @printf("  baseline θ=3 : J_det=%.4f  J_est=%.4f  J_c=%.4f\n", Jd0, Je0, Jd0+Je0)
    @printf("  ContEst θ*∈[%.2f,%.2f]: J_det=%.4f (-%.1f%%)  J_est=%.4f (-%.1f%%)  J_c=%.4f (-%.1f%%)\n",
            minimum(θ1),maximum(θ1), Jd1,100*(Jd0-Jd1)/Jd0, Je1,100*(Je0-Je1)/Je0, Jd1+Je1,100*((Jd0+Je0)-(Jd1+Je1))/(Jd0+Je0))
    Jdes0, Jdes1 = c_k*sum(θ_nom), c_k*sum(θ1)
    @printf("  J_tot (incl. J_des): %.4f → %.4f  (-%.1f%%)\n",
            Jd0+Je0+Jdes0, Jd1+Je1+Jdes1, 100*((Jd0+Je0+Jdes0)-(Jd1+Je1+Jdes1))/(Jd0+Je0+Jdes0))

    # ---- (B) sensor frontier: pruning (curve A) vs joint co-design (curve B) ----
    println("\n── (B) Sensor frontier: baseline pruning vs joint α+θ co-design ──")
    estA, subsetA = prune_baseline()
    JcontA = Jdet_ric(θ_nom)                    # control cost fixed under baseline-droop pruning
    println("    [curve A] baseline-droop pruning done")
    estB, contB = joint_frontier(subsetA)
    @printf("  %-9s %-14s %-14s %-14s %-14s\n", "#sensors", "prune J_est", "joint J_est", "prune J_cont", "joint J_cont")
    for k in m:-1:1
        a  = isfinite(estA[k])  ? @sprintf("%.4f", estA[k])  : "—"
        b  = isfinite(estB[k])  ? @sprintf("%.4f", estB[k])  : "—"
        bc = isfinite(contB[k]) ? @sprintf("%.4f", contB[k]) : "—"
        @printf("  %-9d %-14s %-14s %-14.4f %-14s\n", k, a, b, JcontA, bc)
    end
    mtdc_figure(estA, estB, JcontA, contB)

    # ---- (C) merged optimum: pure joint α+θ selection at a mid-sparsity ℓ1 weight ----
    println("\n── (C) Merged optimum (pure joint α+θ selection) ──")
    zc = joint_opt(; λs = 0.05); θ2 = zc[1:m]; sel2 = zc[m+1:2m] .> 0.5
    Jd2, Je2 = Jdet_ric(θ2), est_true_at(θ2, sel2)
    @printf("  kept %d/%d sensors: %s\n", count(sel2), NT, string(findall(sel2)))
    @printf("  droop θ**∈[%.2f,%.2f]: J_det=%.4f  J_est=%.4f  J_c=%.4f\n",
            minimum(θ2),maximum(θ2), Jd2, Je2, Jd2+Je2)
    @printf("  vs baseline (θ=3, all %d sensors): J_c %.4f→%.4f (-%.1f%%),  sensors %d→%d\n",
            NT, Jd0+Je0, Jd2+Je2, 100*((Jd0+Je0)-(Jd2+Je2))/(Jd0+Je0), NT, count(sel2))
    Jdes2 = c_k*sum(θ2)
    @printf("  J_tot (incl. J_des): %.4f → %.4f  (-%.1f%%)\n",
            Jd0+Je0+Jdes0, Jd2+Je2+Jdes2, 100*((Jd0+Je0+Jdes0)-(Jd2+Je2+Jdes2))/(Jd0+Je0+Jdes0))
    println(bar)
    return (; θ1, Jd0, Je0, Jd1, Je1, θ2, sel2, Jd2, Je2, estA, estB, contB, JcontA)
end

if abspath(PROGRAM_FILE) == @__FILE__
    mtdc_sparse_report()
end
