# example_mtdc_alloc.jl — MTDC droop coordination WITH sparse voltage-sensor
# allocation via the ℓ1 epigraph transformation.  (Variant of example_mtdc.jl;
# neither the paper nor the existing sources are modified — this file only
# *includes* them read-only.)
#
# ─────────────────────────────────────────────────────────────────────────────
# What is new relative to example_mtdc.jl
# ─────────────────────────────────────────────────────────────────────────────
# The paper's MTDC study designs only the control-side droop θ_f=k (25+0: no
# sensing axis), and the ℓ1 sensor-allocation mechanism of the paper's Section 7
# is described but never exercised numerically.  Here we ADD a genuine sensing
# axis and demonstrate that mechanism on the robust H∞ output-feedback inner
# value:
#
#   • Not every DC terminal need carry a telemetered voltage sensor: the meshed
#     ring couples neighbouring bus voltages, so an un-sensed terminal can be
#     reconstructed by the H∞ filter from its neighbours.  WHICH voltage sensors
#     to install is a discrete allocation, traded against estimation quality and
#     a per-sensor telemetry/installation cost.
#
#   • Continuous relaxation (Section 7): give each measurement row a gain
#     αᵢ ≥ 0, so C(α) = diag(α)·C_base (αᵢ=0 removes sensor i; larger αᵢ ⇒
#     higher SNR, since the per-channel noise V is fixed).  α joins the design θ
#     and enters the estimation half through C only — the exact envelope gradient
#     ∂J_est/∂C is already supplied by Hinf_of_θ, so NO new gradient machinery.
#
#   • ℓ1 epigraph transformation.  We penalise the number of active sensors by
#     λ_s‖α‖₁ inside J_des.  The epigraph lift introduces slacks t with
#     λ_s‖α‖₁ ↦ λ_s 1ᵀt and −t ⪯ α ⪯ t; here α ≥ 0, so the optimum is t = α and
#     the lifted objective is simply the SMOOTH linear term λ_s 1ᵀα on the
#     nonnegative box (exactly the reduction stated in Section 7).  Smoothness is
#     what lets the bound-constrained (L-)BFGS of BFGS.jl keep its convergence
#     guarantees instead of degrading to a subgradient method.  Its design
#     gradient is the constant λ_s on the α-block — added to the envelope
#     gradient with nothing else changed.
#
# We JOINTLY co-design the droop k (control, enters A) and the sparse sensing α
# (estimation, enters C), turning the paper's "25+0" study into a genuine
# two-axis ContEst problem (25+25) with a discrete sensor menu.
#
# Run:  julia --project benchmarks/example_mtdc_alloc.jl
# ─────────────────────────────────────────────────────────────────────────────

include("example_mtdc.jl")   # brings in NT, dt, Cdc, τcnv, Gdc, k0, γ², σinf,
                             # Qmat, Rmat, Cmat, Vmeas, nb, of_eval, c_droop,
                             # and the whole H∞/BFGS machinery.  Its bottom-of-
                             # file report is guarded by `PROGRAM_FILE==@__FILE__`,
                             # so including it runs NOTHING.
using LinearAlgebra, Printf

# ── Design layout: θ = [ k(1:NT) ; α(NT+1:2NT) ] ─────────────────────────────
const Cbase = copy(Cmat)                       # NT×2NT voltage-selection base map
const τ_on  = 0.05                             # threshold: αᵢ ≤ τ_on ⇒ sensor OFF

# Model with per-terminal droop kᵢ (in A) AND per-terminal sensor gain αᵢ (in C).
# Identical physics to example_mtdc.jl; the ONLY change is C ↦ diag(α)·C_base.
function model_alloc(θ)
    T = eltype(θ)
    Ac = zeros(T, 2NT, 2NT); Bc = zeros(T, 2NT, NT)
    for i in 1:NT
        v = 2i - 1; c = 2i
        (l, r) = nb(i)
        Ac[v, v] = -2Gdc / Cdc[i]
        Ac[v, 2l-1] += Gdc / Cdc[i]
        Ac[v, 2r-1] += Gdc / Cdc[i]
        Ac[v, c]  = 1.0 / Cdc[i]
        Ac[c, c]  = -1.0 / τcnv[i]
        Ac[c, v]  = -(k0 * θ[i]) / τcnv[i]     # droop kᵢ = θ[i] enters A
        Bc[c, i]  = 1.0 / τcnv[i]
    end
    Ad = Matrix(1.0I, 2NT, 2NT) + dt * Ac
    Bd = dt * Bc
    W = zeros(2NT, 2NT)
    for i in 1:NT
        W[2i-1, 2i-1] = (dt * σinf[i] / Cdc[i])^2
        W[2i, 2i] = 1e-8
    end
    α = @view θ[NT+1:2NT]
    C = Diagonal(α) * Cbase                     # sensing gains enter C (Section 7)
    return Ad, Bd, W, Qmat, Rmat, Matrix(C), Vmeas
end

# ── Design cost:  droop stress (existing) + ℓ1 sensor-allocation penalty ──────
# J_des(θ) = c_droop·Σkᵢ + λ_s·Σαᵢ .  The second term is the ℓ1 epigraph penalty
# specialised to α≥0 (t=α): smooth and linear, gradient constant λ_s on α.
make_costs(λs) = (
    θ -> c_droop * sum(@view θ[1:NT]) + λs * sum(@view θ[NT+1:2NT]),
    θ -> vcat(fill(c_droop, NT), fill(λs, NT)),
)

# ── Boxes and baseline (naive: over-provisioned droop, ALL sensors deployed) ──
const θnom_a = vcat(fill(3.0, NT), fill(1.0, NT))       # k=3 (over-provisioned), α=1 (all on)
const lb_a   = vcat(fill(0.5, NT), fill(0.0, NT))       # αᵢ≥0 (0 ⇒ sensor removed)
const ub_a   = vcat(fill(4.0, NT), fill(3.0, NT))

n_active(α) = count(>(τ_on), α)

# ── Objective builder reusing the output-feedback H∞ evaluator (of_eval) ──────
function build_obj(λs)
    Jdes, ∇Jdes = make_costs(λs)
    eval_J(θ) = begin
        Jc, ∇Jc, _, _, _, _ = of_eval(model_alloc, θ)
        (Jc + Jdes(θ), ∇Jc .+ ∇Jdes(θ))
    end
    f, g! = cached_objective(eval_J)
    return f, g!, Jdes
end

# ── Cost breakdown at a design ───────────────────────────────────────────────
function costs_at(θ, Jdes)
    _, _, _, Jd, Je, ok = of_eval(model_alloc, θ)
    (; ok, J_est = Je, J_det = Jd, J_c = Jd + Je, J_des = Jdes(θ), J_tot = Jd + Je + Jdes(θ),
       nsens = n_active(@view θ[NT+1:2NT]))
end

if abspath(PROGRAM_FILE) == @__FILE__
    bar = "=" ^ 78
    println(bar)
    println("  MTDC + sparse voltage-sensor allocation (output-feedback H∞, γ²=$(γ²))")
    println("  design = [droop k(1:$NT) ∈ A] + [sensor gains α(1:$NT) ∈ C], ℓ1 penalty on α")
    println(bar)

    λs_headline = 0.06         # tuned below via the sweep; clear partial sparsity, filter stays admissible

    # ── 1. Gradient verification (confirms exactness through the NEW C(α) axis) ──
    f0, g0!, Jdes0 = build_obj(λs_headline)
    ∇a, ∇fd, rel = verify_gradient(f0, g0!, θnom_a)
    @printf("\n[1] Gradient check at baseline (50 params): rel.err = %.2e  %s\n",
            rel, rel < 5e-2 ? "(OK — envelope gradient exact through C(α))" : "(WARNING)")

    # ── 2. λ_s sparsity/performance frontier (single BFGS from the naive baseline) ──
    println("\n[2] ℓ1 sparsity/performance frontier (single-start BFGS from baseline)")
    @printf("    %-8s %8s %10s %10s %10s %10s %10s\n",
            "λ_s", "#sens", "J_est", "J_det", "J_c", "J_des", "J_tot")
    base0 = costs_at(θnom_a, Jdes0)
    @printf("    %-8s %8d %10.4f %10.4f %10.4f %10.4f %10.4f   (naive baseline)\n",
            "—", base0.nsens, base0.J_est, base0.J_det, base0.J_c, base0.J_des, base0.J_tot)
    for λs in (0.0, 0.01, 0.03, 0.06, 0.10, 0.20, 0.40)
        f, g!, Jdes = build_obj(λs)
        r = bfgs_design(f, g!, lb_a, ub_a, copy(θnom_a); iterations = 200, g_tol = 1e-5)
        c = costs_at(r.minimizer, Jdes)
        @printf("    %-8.2f %8d %10.4f %10.4f %10.4f %10.4f %10.4f\n",
                λs, c.nsens, c.J_est, c.J_det, c.J_c, c.J_des, c.J_tot)
    end

    # ── 3. Headline joint co-design at λ_s = λs_headline (multi-start) ──────────
    println("\n[3] Headline joint co-design (droop + sparse sensing), multi-start")
    Jdes, ∇Jdes = make_costs(λs_headline)
    best_θ, best_J, base, opt = hinf_of_report(
        "MTDC droop + sparse voltage-sensor allocation (λ_s=$(λs_headline))",
        vcat(["k[$i]" for i in 1:NT], ["α[$i]" for i in 1:NT]),
        vcat(["A: droop @ term $i (σ_inf=$(round(σinf[i],digits=1)))" for i in 1:NT],
             ["C: V-sensor @ term $i" for i in 1:NT]);
        of_eval = of_eval, model = model_alloc, J_des = Jdes, ∇J_des = ∇Jdes,
        θ_init = θnom_a, θ_nom = θnom_a, θ_lb = lb_a, θ_ub = ub_a,
        γ² = γ², seed = 20240624, n_starts = 5)

    # ── 4. Sensor-menu summary + threshold/re-solve ────────────────────────────
    αopt = best_θ[NT+1:2NT]
    on   = findall(>(τ_on), αopt)
    off  = findall(<=(τ_on), αopt)
    println("\n[4] Deployed sensor menu at λ_s=$(λs_headline)  (threshold τ=$(τ_on))")
    @printf("    sensors kept: %d / %d      pruned: %s\n",
            length(on), NT, isempty(off) ? "(none)" : join(off, ", "))
    @printf("    active gains α: [%s]\n",
            join((@sprintf("%.2f", αopt[i]) for i in on), ", "))

    # Threshold to a binary deployment, then RE-SOLVE the droop with the sensor
    # set fixed (αᵢ∈{0, re-optimised}); report the realised worst-case cost of the
    # deployed (integer-sensor) design.
    keep_mask = αopt .> τ_on
    Jdes_fixed = θ -> c_droop * sum(@view θ[1:NT])          # sensor cost now sunk
    ∇Jdes_fixed = θ -> vcat(fill(c_droop, NT), zeros(NT))
    function eval_fixed(θ)
        θf = copy(θ); θf[NT+1:2NT] .= ifelse.(keep_mask, θ[NT+1:2NT], 0.0)
        Jc, ∇Jc, _, _, _, _ = of_eval(model_alloc, θf)
        (Jc + Jdes_fixed(θf), ∇Jc .+ ∇Jdes_fixed(θf))
    end
    ff, gf! = cached_objective(eval_fixed)
    θ0f = copy(best_θ); θ0f[NT+1:2NT] .= ifelse.(keep_mask, best_θ[NT+1:2NT], 0.0)
    rf = bfgs_design(ff, gf!, lb_a, ub_a, θ0f; iterations = 200, g_tol = 1e-5)
    θdep = copy(rf.minimizer); θdep[NT+1:2NT] .= ifelse.(keep_mask, rf.minimizer[NT+1:2NT], 0.0)
    cdep = costs_at(θdep, Jdes)
    # Report the PERFORMANCE quantities (J_est, J_det, J_c) of the deployed
    # integer-sensor design; these are what the sparse menu actually buys.  We do
    # NOT roll J_des into a total here, because after committing to a menu the
    # per-sensor cost is a fixed installation count, not the λ_s‖α‖₁ surrogate
    # (whose value depends on the re-grown gains of the kept sensors).
    @printf("\n[5] Thresholded + re-solved deployment (%d of %d sensors, integer menu):\n",
            length(on), NT)
    @printf("    J_est=%.4f  J_det=%.4f  J_c=J_det+J_est=%.4f\n",
            cdep.J_est, cdep.J_det, cdep.J_c)
    @printf("    vs naive 25-sensor baseline:  J_est %.4f→%.4f (%+.1f%%),  J_c %.4f→%.4f (%.1f%% lower)\n",
            base0.J_est, cdep.J_est, 100*(cdep.J_est-base0.J_est)/base0.J_est,
            base0.J_c, cdep.J_c, 100*(base0.J_c-cdep.J_c)/base0.J_c)
    @printf("    ⇒ %d fewer voltage sensors (%.0f%% of the telemetry removed) with estimation cost no worse.\n",
            NT - length(on), 100*(NT-length(on))/NT)
    println(bar)
end
