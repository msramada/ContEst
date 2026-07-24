# MTDC example — archived from the paper (formulation + results)

Snapshot of the **multi-terminal HVDC (MTDC)** study **as it appeared in
`ContEst_TeX/main.tex`** (subsection `\label{sec:sim:mtdc}`, "Multi-terminal
HVDC: robust droop co-design with $H_2$ sparse sensing") **before it was removed
from the manuscript** and the linear-robust example role was taken over by the
reformulated ADCS study (`benchmarks/example_adcs_hinf.jl`, H∞ bounded-real-lemma
SDP control + H₂ estimation SDP, both gradients from the LMI duals).

Captured here so the formulation, numbers, figure, and tables remain documented
and recoverable. The code is untouched and still runs:

- `benchmarks/example_mtdc.jl` — base library (model, `est_sdp`, `est_true_cost`,
  `kalman_cost`, `Jdet_ric`, `hinf_gare`).
- `benchmarks/example_mtdc_sparse.jl` — the official study + figure
  (`ContEst_TeX/figs/contest_mtdc.pdf`).
- Run: `julia --project benchmarks/example_mtdc_sparse.jl`.

(See also `MTDC_old_example_summary.md` for the even-earlier 50-state droop-only
version that this study itself replaced.)

---

## 1. Role in the paper

MTDC was the paper's **linear, robust, large-scale** study — the single case that
exercised the steady-state SDP path on a robust inner value, in contrast to the
three nonlinear eKF–MPC studies (ADCS, distillation, PLL). It carried three
demonstrations that no other study did:

1. **H∞ control (worst-case) + H₂ estimation on one shared design.**
2. **Structured (block-diagonal / decentralized) SDP scalability** — the ring
   sparsity cutting certificate/gain DOF from `O(r_x²)` to `O(r_x)`.
3. **ℓ₁ sparse sensor selection** co-optimized with the plant parameter.
4. **Shared-dynamics coupling**: the droop `θ` enters `A(θ)`, so it reshapes
   *both* the H∞ control value and the estimation error dynamics — a design-level
   coupling a place-then-control pipeline cannot exploit.

> Note on what its removal costs the paper: demonstrations (2)–(4) above are no
> longer exercised by any example. The reformulated ADCS study reproduces the
> H∞-control + H₂-estimation SDP pairing (demonstration 1), but couples the two
> axes through a **shared power/mass budget** (a design-cost coupling) rather than
> through `A(θ)`, is **dense** (6 states, no structured/scalable SDP), and does
> **no sensor selection**. The block-diagonal-SDP and ℓ₁-selection theory
> (Remarks in `sec:hinf`/`sec:envelope`) therefore remains as theory without a
> worked example.

---

## 2. System and formulation

Reduced DC grid of `n_t = 15` terminals on a meshed DC ring, each with a DC-bus
voltage and converter-current state `[ΔV_i, ΔI_i]`, for `r_x = 2 n_t = 30`
states; `r_u = 15` inputs; `r_y = 15` candidate voltage sensors. `Δt = 5 ms`.

Linearized, Euler-discretized state space:

```
x_{k+1} = A(θ) x_k + B u_k + E w_k
y_k     = C(α) x_k + D v_k,     z_k = [Q^{1/2} x_k ; R^{1/2} u_k]
```
with `A(θ) = I + Δt·A_c(θ)`, `B = Δt·B_c`, disturbance channel `E = W^{1/2}`
(renewable infeed), performance output `z_k`.

Continuous-time device model (per terminal `i`, neighbours `j ∼ i` on the ring):
```
C_i^dc ΔV̇_i = -Σ_{j∼i} G_ij (ΔV_i - ΔV_j) + ΔI_i + infeed
τ_i    Δİ_i = -ΔI_i - (k_0 θ_i) ΔV_i + u_i
y_{k,i}     = α_i ΔV_i + v_{k,i},     v_k ~ N(0, v I_15)
```
- `C_i^dc` DC-bus capacitance, `G_ij` ring line conductances, `τ_i` converter lag.
- The **droop** law feeds `-(k_0 θ_i) ΔV_i` into the current with scale `k_0 = 12`;
  this is the **only** place `θ` enters `A`, through converter-current feedback.
- The **sensor gains** `α_i` enter `C`. Converter currents `ΔI_i` are unmeasured.

Weights / noise / box:
- `Q` penalizes `ΔV` (and weights the estimation error); `R = 2 I`.
- Heterogeneous renewable-infeed disturbance `σ_inf,i` per node drives the DC nodes.
- Measurement covariance `V = 10^{-3} I_15`, `D = V^{1/2}`.
- Design cost (converter-stress): `J_des = c_k Σ_i θ_i`, `c_k = 5×10^{-3}`.
- Design box `θ ∈ [0.5, 4]^15`, naive over-provisioned baseline `θ_nom = 3`
  (tight regulation ignoring the resonance).

**Design axes:** `θ_f = k` = droop gains (enter `A`) — control; `θ_h = α` = sensor
gains (enter `C`) — sensing.

### Control side — H∞ game at fixed attenuation `γ² = 16`

Worst-case value `J_cont★(θ) = tr(X W)`, `X` the stabilizing solution of the
discrete H∞ game Riccati (GARE, paper eq. `\eqref{eq:gare}`).

**Official (structured) form — block-diagonal bounded-real-lemma SDP.** In
variables `Y = P⁻¹`, `L = K Y`, `J_cont★ = min tr(Y⁻¹ W)` s.t. the BRL LMI at
fixed γ. A full `(Y,L)` recovers the game/central solution (dense optimum equals
`tr(XW)`). Restricting `Y, L` to **block-diagonal** (one 2×2 block per terminal —
a decentralized certificate/controller) cuts DOF `O(r_x²) → O(r_x)`, is strictly
complementary once `γ² ≥ 16`, and its BRL dual returns the exact envelope
gradient `∂J_cont★/∂A_θ` from the plant block (no solver differentiation).

### Sensing side — H₂ estimation SDP in the observability-gramian (P^ε) form

```
J_est★(θ,α) = min_{P^ε, F, Z}  tr(W P^ε) + tr(V Z)
   s.t.  [Z  Fᵀ ; F  P^ε] ⪰ 0
         [P^ε − Q   AᵀP^ε − CᵀFᵀ ; P^ε A − F C   P^ε] ⪰ 0
```
filter gain `G = (P^ε)⁻¹ F`. `P^ε, F` restricted **block-diagonal** (one 2×2
storage block + one 2×1 gain block per terminal — a decentralized filter). Exact
envelope gradients from the Lyapunov-LMI dual `S'`:
```
∂J_est★/∂A_θ = -2 P^ε S'_12ᵀ ,   ∂J_est★/∂C_θ = 2 Fᵀ S'_12ᵀ
```
Sensor gains `α_i ∈ [0,1]` scale rows `C_i(α) = α_i e_{2i-1}ᵀ`. In the P^ε form
`C_θ` enters bilinearly with `F`, so `α` is an **outer** variable optimized by the
projected method through `∂J_est★/∂α_i = [∂J_est★/∂C_θ]_{i,2i-1}`, with an ℓ₁
penalty `λ_s 1ᵀα` selecting the sensor set.

No covariance caps: every reported estimation number is the **true** steady-state
cost `J_est = tr(Q Σ^ε_true)`, `Σ^ε_true = dlyap(A_θ − G C_θ, W + G V Gᵀ)` (error
covariance of the SDP's decentralized gain on the full plant).

**Coupling.** Because the droop enters `A(θ)` it reshapes both `X` (game Riccati)
and the estimation error dynamics `A_θ − G C_θ`, so `J_cont★` and `J_est★` vary
over a common `θ` and the co-design does not separate.

**Why H∞ for the droop.** The failure mode of a droop grid is a *resonant peak*:
converter lag makes excessive droop provoke a lightly damped inter-terminal
DC-voltage oscillation; an infeed step near resonance is amplified far beyond its
RMS level. H∞ prices that worst-case amplification (H₂ averages it away), so
`J_cont★` is sharply U-shaped in the droop and the tight-regulation baseline
over-provisions onto the resonance wall.

---

## 3. Results

### Structured vs. dense (Table `tab:mtdc-sparse`, at `γ² = 16`, baseline droop)

| formulation | # variables (Y,L) | solve time |
|---|---|---|
| full dense `(Y,L)` | 915 | 118 s |
| block-diagonal (official) | 75 | 0.37 s |

- Block-diagonal uses ≈12× fewer variables and solves ≈320× faster.
- Its price is a conservative worst-case value: the decentralized restriction
  over-bounds `J_cont★` by ≈56% (`1.61 → 2.50` at baseline).
- The **design** it returns is nearly lossless: co-optimizing the droop on the
  block-diagonal surrogate and scoring on the true model gives
  `J_c` within **+1.1%** of the full-dense (game-Riccati) co-design
  (full-dense `1.362` vs. block-diagonal `1.377`), both ≈40% below baseline.
- Estimation validity: at baseline the block-diagonal filter's true estimation
  cost is only **+2.1%** above the centralized Kalman value (+0.1% for the banded
  certificate). The banded / nearest-neighbour certificate matches the dense value
  to ~1% but stays dual-degenerate at every γ² (fast evaluator, not a gradient
  source).

### (A) Droop co-design (all 15 sensors)

From baseline (`J_cont★ = 1.61`, `J_est★ = 0.66`), ContEst lowers the droop to an
interior, disturbance-dependent profile `θ★ ∈ [0.93, 1.44]` (none at a bound):

| quantity | baseline | co-design | reduction |
|---|---|---|---|
| `J_cont★` | 1.61 | 0.98 | −39.1% |
| `J_est★`  | 0.66 | 0.40 | −39.9% |
| `J_c`     | 2.27 | 1.38 | −39.3% |

The single control-side axis improves the estimation half purely by reshaping the
shared dynamics `A(θ)` — the design-level coupling a place-then-control pipeline
cannot exploit.

### (B) Sensor frontier — pruning vs. co-design (Figure `fig:mtdc`)

Sweeping the retained-sensor count: greedy pruning at the fixed baseline droop
degrades gracefully (ring coupling lets the filter reconstruct un-sensed terminals
from neighbours), but the **joint (α, θ) co-design lies far below it at every
count** — co-designing the droop for each sensor budget cuts estimation cost by
35–40% vs. pruning alone:

- at 7 of 15 sensors: `0.53` (joint) vs. `0.88` (prune);
- even at a single sensor: `0.86` (joint) vs. `1.43` (prune).

This estimation gain is bought at **no cost to control**: joint control cost holds
at its robust optimum `tr(XW) ≈ 0.98` across the whole sweep (far below the
over-provisioned baseline `1.61`). Reshaping `A(θ)` improves sensing while
preserving robustness.

Figure `figs/contest_mtdc.pdf`: two panels vs. number of retained voltage sensors
— (a) estimation cost `tr(QΣ^ε_true)`, (b) worst-case control cost `tr(XW)`.
Forest-green dashed + circles = baseline-droop pruning; red solid + diamonds =
joint (α,θ) co-design.

### (C) Merged optimum

Letting the ℓ₁ selection and droop co-optimize freely, ContEst keeps **7 of 15
sensors, {3,4,6,7,11,12,14}**, with droop `θ★★ ∈ [0.92, 1.67]`:

| quantity | value |
|---|---|
| `J_cont★` | 0.98 |
| `J_est★`  | 0.53 |
| `J_c`     | 1.51 (−33.6% vs. full-telemetry baseline) |

— removing over half the voltage sensors, in one envelope loop that tunes the
droop and allocates sensors together.

### Consolidated (as in `tab:formulations`)

| study | (r_x, r_u, r_y) | θ: n_f+n_h (enters) | inner value | J_est↓ | J_cont↓ | J_tot↓ |
|---|---|---|---|---|---|---|
| MTDC | (30, 15, 15) | 15+15: f_θ, C(α) | H∞ + H₂ | 39.9% (15→7 sensors) | 39.1% | **41.4%** |

---

## 4. Novelty framing (from the paper)

The novelty was **not** the sensor selection itself (fixed-plant sensor selection
is mature: greedy/submodular, convex relaxation, global MISDP). What none of those
represents is that here the control-side droop `θ_f` **reshapes** the estimation
problem — it moves the error dynamics `A(θ) − GC` on which the filter runs — so
ContEst co-designs the plant parameter and the sensor set together from a single
envelope gradient, whereas every fixed-plant selector optimizes sensors over a
*frozen* plant. It is the shared-dynamics coupling, not the allocation, that is
beyond their reach.
