# example_mtdc.jl — ContEst App 5: multi-terminal HVDC (MTDC) DC-voltage droop
# coordination  (X_applications.md §5).  ROBUST (H∞) inner value.
#
# Physical setup
# ──────────────
# NT converter terminals share a meshed DC grid. Each terminal regulates its DC
# voltage with a droop: the converter current responds to the local voltage
# deviation with gain k_i (through a first-order converter lag). Droop gains are
# a coupled, non-monotone knob — too little droop gives poor voltage regulation,
# but the converter lag means excessive droop provokes inter-terminal DC-voltage
# oscillation, and raising one terminal's droop shifts burden onto its
# neighbours. The optimal droop profile is an interior, network-coupled trade.
#
# Why H∞ (not H₂) here.  The failure mode is a RESONANT PEAK: renewable-infeed
# steps excite a lightly-damped inter-terminal oscillation whose worst-case
# amplification, not its noise-averaged variance, is what threatens the DC grid.
# H∞ prices exactly that worst case, so we replace the H₂/LQG inner value with the
# steady-state H∞ (worst-case) guaranteed cost V(θ)=tr(XW), where X solves the
# discrete H∞ game Riccati at a fixed disturbance-attenuation level γ (Hinf.jl).
# This is an indefinite-weight DARE, so it stays on the SAME scalable O(n³)
# Riccati path as the H₂ solver, with the same envelope gradient.
#
# States  x = [ΔV_i, ΔI_i]  per terminal   (DC-bus voltage, converter current)
# Inputs  u = supplementary power-reference modulation per terminal
# Outputs y = ΔV_i (DC-bus voltages telemetered; converter currents estimated)
# Perf.   z = [Q^{1/2}x; R^{1/2}u]; disturbance channel E=√W (renewable infeed)
# Design  θ = per-terminal droop coefficient k_i   (enters A: droop feedback)
#
# Output feedback ⇒ TWO worst-case game Riccatis: a control GARE (J_det=tr(XW))
# and a dual filter GARE (J_est=tr(QΣe)). The droop θ enters A, so it moves BOTH,
# coupling control and estimation at the design level (as in the LQG separation).
#
# Run:  julia --project benchmarks/example_mtdc.jl

include("../src/LQR.jl")     # dare, dlyap (shared Riccati/Lyapunov machinery)
include("../src/BFGS.jl")    # verify_gradient, multistart_design, cached_objective
include("../src/Hinf.jl")    # Hinf_θ, hinf_report (H∞ game-Riccati inner solver)
using LinearAlgebra, Printf

# ── DC grid ──────────────────────────────────────────────────────────────────
const NT = 25                     # terminals ⇒ 2·NT = 50 states, NT design params
const dt = 0.005                  # fast DC dynamics ⇒ small step [s]
nb(i) = (i == 1 ? NT : i - 1, i == NT ? 1 : i + 1)      # ring neighbours

const Cdc  = [0.8 + 0.4sin(0.9i) for i in 1:NT]         # DC-node capacitances
const τcnv = [0.02 + 0.01cos(1.1i) for i in 1:NT]       # converter current lag [s]
const Gdc  = 5.0                                         # DC line conductance (ring)
const k0   = 12.0                                        # droop scale (leverage): the
# worst-case cost V(θ)=tr(XW) is a pronounced U-shape in each θ_i — the interior
# robust optimum is near θ≈1, and pushing droop higher climbs a steep converter-lag
# resonance wall that H∞ penalises directly (the worst-case peak, not its average).

const γ²   = 8.0   # H∞ disturbance-attenuation level (smaller ⇒ more robust/conservative).
# Chosen so the worst-case game Riccati is admissible over the WHOLE box [0.5,4]^NT
# (incl. the θ_nom=3 baseline) with margin; at this γ the resonance wall is steep
# but every design remains feasible. As γ→∞ the value collapses to the H₂ tr(PW).

# Renewable-infeed disturbance intensity on each DC node, heterogeneous.
const σinf = [1.0 + 1.5 * (0.5 + 0.5sin(1.7i + 2)) for i in 1:NT]

const Qmat = Matrix(Diagonal([iseven(s) ? 0.2 : 60.0 for s in 1:2NT]))  # penalise ΔV ≫ ΔI
const Rmat = Matrix(2.0 * I(NT))            # moderately expensive supplementary control

# ── Measurement model (output feedback): DC-bus VOLTAGES are telemetered, the
# converter CURRENTS are not — so the estimator must reconstruct the current states
# from the voltage measurements. This makes MTDC a genuine output-feedback H∞
# co-design with an estimation half (worst-case H∞ filter), not full-state control.
const meas_rows = [2i - 1 for i in 1:NT]                 # measured states = ΔV_i
const Cmat  = Matrix(1.0I, 2NT, 2NT)[meas_rows, :]       # NT × 2NT voltage-selection map
const v_meas = 1.0e-3                                    # DC-voltage measurement variance
const Vmeas  = Matrix(v_meas * I(NT))

# ── Linear model  θ ↦ (A,B,W,Q,R,C,V)  (discrete-time) ───────────────────────
# θ_i = droop-gain scaling: k_i = k0·θ_i enters the converter-current equation
#   τ İ_i = −I_i − k_i ΔV_i + u_i
# and couples to the DC-node capacitor balance
#   C_i V̇_i = −Σ_j G_ij(V_i−V_j) + I_i + w_i .
function model(θ)
    T = eltype(θ)
    Ac = zeros(T, 2NT, 2NT); Bc = zeros(T, 2NT, NT)
    for i in 1:NT
        v = 2i - 1; c = 2i
        (l, r) = nb(i)
        Ac[v, v] = -2Gdc / Cdc[i]                 # C V̇ = -ΣG(V-Vnb) + I
        Ac[v, 2l-1] += Gdc / Cdc[i]
        Ac[v, 2r-1] += Gdc / Cdc[i]
        Ac[v, c]  = 1.0 / Cdc[i]
        Ac[c, c]  = -1.0 / τcnv[i]                # τ İ = -I - k ΔV + u
        Ac[c, v]  = -(k0 * θ[i]) / τcnv[i]        # droop feedback (θ enters A)
        Bc[c, i]  = 1.0 / τcnv[i]
    end
    Ad = Matrix(1.0I, 2NT, 2NT) + dt * Ac
    Bd = dt * Bc
    W = zeros(2NT, 2NT)
    for i in 1:NT
        W[2i-1, 2i-1] = (dt * σinf[i] / Cdc[i])^2   # infeed disturbance on the DC node
        W[2i, 2i] = 1e-8
    end
    return Ad, Bd, W, Qmat, Rmat, Cmat, Vmeas
end

# ── Design cost: droop provisioning / converter stress ∝ Σ k_i ───────────────
# Baseline = the NAIVE over-provisioned droop: an engineer chasing tight DC-voltage
# regulation sets a uniformly HIGH droop, not accounting for the converter lag that
# makes excessive droop provoke inter-terminal resonance. Because the worst-case cost
# V(θ) is U-shaped in each k_i, this baseline sits on the far (resonant) wall of the
# bowl, so ContEst must PULL most terminals' droop DOWN toward the interior optimum
# (and push a few under-provisioned ones up) — a genuinely non-monotone re-allocation,
# and one that also lowers converter stress (design cost), for a material reduction.
const k_nom   = 3.0                        # naive uniform (over-)provisioned droop scaling
const θ_nom   = fill(k_nom, NT)
const c_droop = 0.005                      # small converter-stress price (U-shape in J_c dominates)
J_des(θ)  = c_droop * sum(θ)
∇J_des(θ) = fill(c_droop, NT)

const θ_lb = fill(0.5, NT)
const θ_ub = fill(4.0, NT)

const of_eval = Hinf_of_θ(2NT, NT, NT; γ² = γ²)

if abspath(PROGRAM_FILE) == @__FILE__
    @printf("App 5 — MTDC droop coordination (output-feedback H∞, γ²=%.1f):  %d terminals, %d states, %d meas, %d design params\n",
            γ², NT, 2NT, NT, NT)
    hinf_of_report("App 5: multi-terminal HVDC DC-voltage droop coordination",
                   ["k[$i]" for i in 1:NT],
                   ["A: droop @ terminal (σ_inf=$(round(σinf[i],digits=1)))" for i in 1:NT];
                   of_eval = of_eval, model = model, J_des = J_des, ∇J_des = ∇J_des,
                   θ_init = θ_nom, θ_nom = θ_nom, θ_lb = θ_lb, θ_ub = θ_ub,
                   γ² = γ², seed = 20240624, n_starts = 5)
end
