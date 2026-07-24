# example_pll_hinf.jl — EXPLORATORY reformulation of the PLL / DSE co-design as a
# LINEAR robust/stochastic SDP problem, instead of the nonlinear constrained
# eKF–MPC study in example_pll.jl.
#
# THIS FILE CHANGES NOTHING ELSE. It is a standalone experiment answering the
# question: "can the PLL estimator co-design be posed on the ContEst SDP path
# (H∞ for control, H₂ for estimation) rather than the eKF–MPC path?"
#
# What is done differently from example_pll.jl (deliberately, per the request):
#   • SKIP the input/state constraints (|u_i| ≤ 1.5 dropped) — required to leave
#     the MPC path; a pure state-feedback / steady-state design has no horizon
#     over which to enforce inequalities.
#   • LINEARISE the nonlinear swing dynamics about an OPERATING POINT.  The only
#     nonlinearity is the sinusoidal electrical power Pe(i); linearised about a
#     loading angle set δ_op it gives constant synchronising coefficients
#     K̃_ij = Kc_ij·cos(δ_op_i − δ_op_j).  To avoid over-fitting one operating
#     point we build an ENSEMBLE of operating points {δ_op^(1..M)} and average the
#     inner value and its gradient (a scenario/averaged-gradient robustification):
#         J_c(θ) = (1/M) Σ_s [ J_det^(s)(θ) + J_est^(s)(θ) ],   ∇ averaged likewise.
#
# Two design axes, SAME physical meaning and SAME box as example_pll.jl:
#   • θ_f = e_i ∈ [0.3,3]     inverter damping-actuator effectiveness  (enters B)
#   • θ_h = b_i ∈ [2,20]      per-bus PLL bandwidth                    (enters A AND V)
#
# Inner solvers (the two ContEst SDP-path values):
#   • CONTROL   H∞ state-feedback game Riccati (src/Hinf.jl `hinf_gare`, fixed γ²):
#                 J_det = tr(X W),  worst-case guaranteed cost, exact envelope grad.
#   • ESTIMATION H₂ covariance SDP in the observability-gramian (P^ε) form
#                 (the same LMI used by example_mtdc.jl, here DENSE / centralised):
#                 J_est = tr(Q Σ) = min tr(W P)+tr(V Z) s.t. two LMIs,
#                 exact envelope gradients ∂J/∂A (through the PLL rows of A(b)) and
#                 ∂J/∂V = Z* (through the θ-dependent readout noise V(b)).
#
# Run:  julia --project benchmarks/example_pll_hinf.jl

include("../src/LQR.jl")     # dare, dlyap
include("../src/Hinf.jl")    # hinf_gare  (needs dare/dlyap in scope)
include("../src/BFGS.jl")    # cached_objective, verify_gradient, multistart_design
using JuMP, Clarabel, LinearAlgebra, ForwardDiff, Printf, Random

# ── Network size & physics (identical to example_pll.jl) ──────────────────────
const NB   = 3
const dt_s = 0.05
const Hgen = [1.6, 2.2, 2.0]
const Dgen = [0.4, 0.6, 0.5]
const δeq  = [0.15, -0.05, 0.10]
const Kc = let K = zeros(NB, NB)
    K[1,2] = K[2,1] = 1.2; K[2,3] = K[3,2] = 1.1; K[1,3] = K[3,1] = 0.4; K
end

const n = 3NB          # states  [Δδ; Δω; Δωᵖ]
const m = NB           # inputs   ΔP
const ry = 2NB         # outputs  [Δδ (PMU angle); Δωᵖ (PLL freq est)]
δr() = 1:NB; ωr() = NB+1:2NB; pr() = 2NB+1:3NB

# ── Noise / cost weights (identical to example_pll.jl) ────────────────────────
const σ_δ2 = 4.0e-3; const σ_ω2 = 2.0e-3; const κ_pll = 0.30
const Σ_w  = Matrix(Diagonal([fill(2.0e-5, NB); fill(5.0e-4, NB); fill(5.0e-4, NB)]))
const Qmat = Matrix(Diagonal([fill(15.0, NB); fill(120.0, NB); fill(0.0, NB)]))
const Rmat = Matrix(Diagonal(fill(1.0, NB)))

# Fixed measurement map C: PMU reads Δδ, PLL channel reads Δωᵖ (no θ in h).
const Cmat = let C = zeros(ry, n)
    for i in 1:NB; C[i, i] = 1.0; C[NB+i, 2NB+i] = 1.0; end; C
end
# θ-dependent measurement covariance: PLL freq-readout variance grows with bandwidth.
Vmat(θ) = Matrix(Diagonal(vcat(fill(σ_δ2, NB),
                               [σ_ω2*(1 + κ_pll*θ[NB+i]) for i in 1:NB])))

const γ²   = 60.0        # H∞ attenuation level (admissible at baseline; see report)

# ── Linearised discrete model about an operating point δ_op ───────────────────
# θ = [e (NB); b (NB)].  e enters B (actuator), b enters the PLL rows of A.
# K̃_ij = Kc_ij·cos(δ_op_i − δ_op_j) are the synchronising coefficients at δ_op.
function model_lin(θ, δ_op)
    T = eltype(θ)
    e = θ[1:NB]; b = θ[NB+1:2NB]
    A = Matrix{T}(I, n, n); B = zeros(T, n, m)
    for i in 1:NB
        # Δδ_i⁺ = Δδ_i + dt·Δω_i
        A[i, NB+i] += dt_s
        # Δω_i⁺ = Δω_i + dt·(−Pe_lin − D·Δω + e·u)/(2H)
        c = dt_s / (2Hgen[i])
        A[NB+i, NB+i] += -c * Dgen[i]
        for j in 1:NB
            j == i && continue
            K̃ = Kc[i,j] * cos(δ_op[i] - δ_op[j])
            A[NB+i, i] += -c * K̃            # ∂Pe_i/∂Δδ_i
            A[NB+i, j] +=  c * K̃            # ∂Pe_i/∂Δδ_j
        end
        B[NB+i, i] = c * e[i]
        # Δωᵖ_i⁺ = Δωᵖ_i + dt·b·(Δω_i − Δωᵖ_i)
        A[2NB+i, NB+i]  += dt_s * b[i]
        A[2NB+i, 2NB+i] += -dt_s * b[i]
    end
    return A, B
end

# ── Ensemble of operating points (loading-angle scenarios) ────────────────────
# Nominal equilibrium plus signed loading perturbations → different K̃ → different A.
const δ_ens = [δeq,
               δeq .+ [0.25, -0.15,  0.20],
               δeq .+ [-0.20, 0.25, -0.15],
               δeq .+ [0.30,  0.10, -0.25],
               δeq .+ [-0.10,-0.20,  0.30]]
const M_ens = length(δ_ens)

# ── CONTROL inner value (H∞) at one operating point: tr(XW) + envelope grad ───
# θ enters A (b, via PLL rows) and B (e). W,Q,R are θ-independent here.
function ctrl_vg(θ, δ_op)
    A, B = model_lin(θ, δ_op)
    E = Matrix(sqrt(Symmetric(Σ_w)))
    X, K̃, Ku, Acl, ok = hinf_gare(A, B, E, Qmat, Rmat, γ²)
    ok || return (1e8, fill(0.0, length(θ)))
    S = dlyap(Acl, Σ_w)
    J = tr(X * Σ_w)
    nw = size(E, 2)
    gA  = 2 * (X * Acl * S)
    gB̃  = 2 * (X * Acl * S * K̃')
    gB  = gB̃[:, 1:m]
    Aj = ForwardDiff.jacobian(t -> vec(model_lin(t, δ_op)[1]), θ)
    Bj = ForwardDiff.jacobian(t -> vec(model_lin(t, δ_op)[2]), θ)
    ∇ = [tr(gA' * reshape(Aj[:, k], n, n)) + tr(gB' * reshape(Bj[:, k], n, m))
         for k in eachindex(θ)]
    return J, ∇
end

# ── ESTIMATION inner value (H₂) at one operating point: dense P^ε covariance SDP ─
#   min tr(W P) + tr(V Z)  s.t.  [Z Fᵀ; F P] ⪰ 0,
#                                [P−Q  AᵀP−CᵀFᵀ; PA−FC  P] ⪰ 0
# value = tr(Q Σ) (centralised); envelope grads ∂J/∂A = −2 P S₁₂ᵀ, ∂J/∂V = Z*.
function est_sdp_dense(A, V; withgrad = true)
    mdl = Model(Clarabel.Optimizer); set_silent(mdl)
    @variable(mdl, P[1:n, 1:n], Symmetric)
    @variable(mdl, F[1:n, 1:ry])
    @variable(mdl, Z[1:ry, 1:ry], Symmetric)
    @constraint(mdl, Symmetric([Z F'; F P]) in PSDCone())
    lmi2 = @constraint(mdl, Symmetric([P .- Qmat  (A'*P .- Cmat'*F');
                                       (P*A .- F*Cmat)  P]) in PSDCone())
    @objective(mdl, Min, tr(Σ_w * P) + tr(V * Z))
    optimize!(mdl)
    st = termination_status(mdl)
    (st in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.NUMERICAL_ERROR)) ||
        return (1e8, zeros(n, n), zeros(ry, ry), st)
    J = objective_value(mdl)
    withgrad || return (J, zeros(n, n), zeros(ry, ry), st)
    Zv = value.(Z)
    S = dual(lmi2); S12 = S[1:n, n+1:2n]
    gA = -2 .* (value.(P) * permutedims(S12))     # ∂J/∂A
    gV = Matrix(Zv)                               # ∂J/∂V = Z*
    (J, gA, gV, st)
end

function est_vg(θ, δ_op)
    A, _ = model_lin(θ, δ_op); V = Vmat(θ)
    J, gA, gV, st = est_sdp_dense(A, V)
    J >= 1e8 && return (1e8, fill(0.0, length(θ)))
    Aj = ForwardDiff.jacobian(t -> vec(model_lin(t, δ_op)[1]), θ)
    Vj = ForwardDiff.jacobian(t -> vec(Matrix(Vmat(t))), θ)
    ∇ = [tr(gA' * reshape(Aj[:, k], n, n)) + tr(gV' * reshape(Vj[:, k], ry, ry))
         for k in eachindex(θ)]
    return J, ∇
end

# Reference: centralised steady-state Kalman estimation cost (validity check).
kalman_cost(A, V) = tr(Qmat * dare(Matrix(A'), Matrix(Cmat'), Σ_w, Matrix(V)))

# ── Averaged inner value over the operating-point ensemble ────────────────────
function inner_vg(θ)
    Jd = 0.0; Je = 0.0; g = zeros(length(θ))
    for δ_op in δ_ens
        jd, gd = ctrl_vg(θ, δ_op)
        je, ge = est_vg(θ, δ_op)
        Jd += jd; Je += je; g .+= gd .+ ge
    end
    return (Jd/M_ens, Je/M_ens, g ./ M_ens)
end

# ── Design cost (identical to example_pll.jl) ─────────────────────────────────
const λ_e = 0.5; const b_nom = 6.0; const λ_b = 0.02
J_des(θ)  = λ_e*sum((θ[1:NB] .- 1).^2) + λ_b*sum((θ[NB+1:2NB] .- b_nom).^2)
∇J_des(θ) = vcat(2λ_e .* (θ[1:NB] .- 1), 2λ_b .* (θ[NB+1:2NB] .- b_nom))

const θ_lb  = [fill(0.3, NB); fill(2.0,  NB)]
const θ_ub  = [fill(3.0, NB); fill(20.0, NB)]
const θ_nom = [ones(NB); fill(b_nom, NB)]
const θ_names = vcat(["e$i" for i in 1:NB], ["b_pll$i" for i in 1:NB])
const θ_roles = vcat(["f/B (inverter damping)" for _ in 1:NB],
                     ["h,A,V (PLL bandwidth)"  for _ in 1:NB])

# ── Total co-design objective (cached value+grad) ─────────────────────────────
function eval_J(θ)
    Jd, Je, g = inner_vg(θ)
    (Jd + Je + J_des(θ), g .+ ∇J_des(θ))
end
const contest_f, contest_g! = cached_objective(eval_J)

function pll_hinf_report(; seed = 20240624, n_starts = 5)
    bar = "="^76
    println(bar)
    println("  ContEst (H∞ control + H₂ estimation SDP)  —  PLL / DSE co-design")
    println("  linearised over $M_ens operating points, averaged gradient, γ²=$γ², $n states")
    println(bar)

    # ── Gradient verification ────────────────────────────────────────────────
    ∇a, ∇fd, rel = verify_gradient(contest_f, contest_g!, θ_nom)
    println("\n── Gradient verification at θ_nom ──────────────────────────")
    @printf("  analytic : [%s]\n", join((@sprintf("%8.3f", v) for v in ∇a), " "))
    @printf("  fin.diff : [%s]\n", join((@sprintf("%8.3f", v) for v in ∇fd), " "))
    @printf("  relative error: %.2e  %s\n", rel, rel < 5e-2 ? "(OK)" : "(WARNING)")

    # ── Estimation-SDP validity vs centralised Kalman (baseline, op 1) ────────
    A0, _ = model_lin(θ_nom, δeq); V0 = Vmat(θ_nom)
    Jsdp, _, _, _ = est_sdp_dense(A0, V0); Jk = kalman_cost(A0, V0)
    println("\n── H₂ estimation SDP validity (baseline, nominal operating point) ──")
    @printf("  dense P^ε SDP tr(QΣ) = %.4f   centralised Kalman = %.4f   gap %+.2f%%\n",
            Jsdp, Jk, 100*(Jsdp-Jk)/Jk)

    # ── Multi-start BFGS ─────────────────────────────────────────────────────
    println("\n── Multi-start BFGS ($n_starts starts: θ_nom + $(n_starts-1) random) ──")
    best_J, best_θ, results = multistart_design(contest_f, contest_g!, θ_lb, θ_ub;
                                                n_starts = n_starts, θ_nom = θ_nom,
                                                seed = seed, g_tol = 1e-5, iterations = 100)
    Js = sort([r.J for r in results]); nd = 1
    for k in 2:length(Js); Js[k]-Js[k-1] > 1e-4 && (nd += 1); end
    @printf("  starts run: %d   distinct minima: %d   best J_tot = %.4f\n", n_starts, nd, best_J)
    println(nd == 1 ? "  → all starts agree: optimum is (numerically) global." :
                      "  → multimodal: reporting the best minimum only.")

    # ── Cost decomposition ───────────────────────────────────────────────────
    function bd(θ)
        Jd, Je, _ = inner_vg(θ)
        (; J_est = Je, J_det = Jd, J_c = Jd+Je, J_des = J_des(θ), J_tot = Jd+Je+J_des(θ))
    end
    base = bd(θ_nom); opt = bd(best_θ)
    pct(b, o) = abs(b) < 1e-9 ? "    — " : @sprintf("%+6.2f%%", 100*(o-b)/b)
    println("\n── Cost components: baseline θ_nom → optimal θ* (ensemble mean) ──")
    @printf("  %-28s %12s %12s %12s %8s\n", "component", "baseline", "optimal", "Δ (incr +)", "%")
    for (nm, b, o) in (("Estimation (H₂)  J_est", base.J_est, opt.J_est),
                       ("Control    (H∞)  J_det", base.J_det, opt.J_det),
                       ("Inner  J_c=J_det+J_est", base.J_c,   opt.J_c),
                       ("Design cost      J_des", base.J_des, opt.J_des),
                       ("TOTAL   J_tot=J_c+J_des", base.J_tot, opt.J_tot))
        @printf("  %-28s %12.4f %12.4f %12.4f  %s\n", nm, b, o, o-b, pct(b, o))
    end
    @printf("\n  Net total reduction: %.4f  (%.2f%%)\n",
            base.J_tot-opt.J_tot, 100*(base.J_tot-opt.J_tot)/base.J_tot)

    # ── Design-parameter table ───────────────────────────────────────────────
    println("\n── Design parameters: baseline → optimal (best minimum) ────")
    @printf("  %-10s %10s %10s %11s   %s\n", "param", "baseline", "optimal", "Δ", "enters/role")
    for i in 1:length(θ_nom)
        @printf("  %-10s %10.3f %10.3f %+11.3f   %s\n",
                θ_names[i], θ_nom[i], best_θ[i], best_θ[i]-θ_nom[i], θ_roles[i])
    end
    println(bar)
    return best_θ, best_J, base, opt
end

if abspath(PROGRAM_FILE) == @__FILE__
    pll_hinf_report()
end
