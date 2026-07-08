# Example A — AGC + wide-area ACE metering: enlarging the control/estimation coupling gap

**Code:** `benchmarks/codesign_grid_examples.jl` (`build_AGC`, baseline A) and
`benchmarks/example_A_variants.jl` (`build_AGC_v`, the one-factor-at-a-time changes).
Solved on the **LQG-SDP path** (`src/LQR.jl`: `LQR_θ` control Riccati + `filter_θ`
dual filter Riccati), the envelope-theorem design gradient read from the LMI duals.

This note records example A, the realistic changes that were tried to widen the
**joint-vs-sequential coupling gap**, and their measured effect. *(Results table
filled from the runs; see the last section for which changes actually helped.)*

---

## 1. Formulation

A multi-area **automatic generation control (AGC) / load-frequency control**
problem co-designed with the wide-area **Area Control Error (ACE) metering**
architecture. The reason this is a promising ContEst target: AGC acts on the
**estimated** ACE, `ACE_i = ΔP_tie,i + B_i Δf_i`, so estimation error feeds directly
into the control loop — control and estimation are coupled through the plant, not
merely through a shared parameter.

### Plant (per area `i`, `N = 12` areas ⇒ `r_x = 36`)

State `x_i = [Δf_i, ΔP_tie,i, ΔP_m,i]` (frequency, net tie-line flow, governor/turbine
mechanical power); input `u_i` (AGC command); `Δt = 0.1 s`. Discrete-time swing +
tie-line + governor dynamics, with a ring tie topology and a small tie-line loss
term (`−Δt·0.5·ΔP_tie`) that breaks the otherwise-conserved `Σ ΔP_tie` mode (which
would be uncontrollable):

```
Δf_i⁺     = Δf_i + Δt/(2H_i)·( ΔP_m,i − ΔP_tie,i − (D_i + d_i)·Δf_i )
ΔP_tie,i⁺ = ΔP_tie,i + Δt·2π·Σ_j T_ij (Δf_i − Δf_j) − Δt·0.5·ΔP_tie,i
ΔP_m,i⁺   = ΔP_m,i − Δt/T_g·( ΔP_m,i + droop·Δf_i − u_i )
```

`H_i ∈ [3,4.5]`, `D_i = 1`, `T_g = 0.4`, `droop = 0.8`, `T_ij = 0.5` (ring). Load
fluctuations (process noise `W`) are concentrated in areas 1–5 (`3e-3` vs `3e-4`).

### Measurements

Each area reports frequency `Δf_i` (cheap, fixed noise `v_f = 2e-4`) and tie-flow
`ΔP_tie,i` (the scarce, designed channel, noise `V_i = v0_tie / p_i`,
`v0_tie = 4e-3`). `r_y = 2N = 24`.

### Design axes and inner value

- `θ_f = (d_1,…,d_N)` — per-area **fast-reserve damping** (enters `A`).
- `θ_h = (p_1,…,p_N)` — per-area **tie-flow metering precision** (enters `V`).

Inner value = stationary LQG cost `J_c = tr(PW) + tr(MΣᵉ)`, with **(#1)** the
control-induced estimation weight `M(θ) = Kᵀ(R+BᵀPB)K` (the "cost of not knowing the
state", so sensing value is set by the controller), and **(#2)** a single **shared
budget** `Σ d_i + η Σ p_i ≈ B` (`η = 0.4`, `B = 8`) so actuation and sensing compete.
`Q` weights `Δf` (100) and the ACE component `ΔP_tie` (20). Baselines are kept
**on-budget** (uniform `d_i = 0.3`, `p_i ≈ 0.92`), so both the reductions vs baseline
and the joint-vs-sequential gap are clean.

### Metrics

- `J_est ↓ / J_det ↓ / J_c ↓ / J_tot ↓` — reductions vs the naive on-budget baseline.
- **`joint − seq`** — the joint optimum's advantage over a smart, estimation-aware
  **sequential** pipeline (fast reserve for control first, then metering for
  estimation). This is the coupling gap we want to enlarge.

Baseline A: `r_x = 36`, gradient check `1.6e-5`, `J_est −11.8%`, `J_det −7.5%`,
`J_c −8.1%`, `J_tot −7.8%`, **`joint − seq ≈ 4.66%`** (4.68% at 3 starts).

---

## 2. Realistic changes tried (one-factor-at-a-time)

Each change is applied alone to baseline A (`build_AGC_v(...)`), keeping everything
else fixed and the baseline on-budget.

1. **Noisy tie telemetry** (`v0_tie: 4e-3 → 3e-2`). Real tie-line MW is SCADA-grade
   (noisier/slower than PMU frequency). Inflates the control-critical estimation cost
   `tr(MΣᵉ)`, so sensing allocation matters more.
2. **Aggressive AGC** (`R → 0.1·R`). Cheap/high-gain control ⇒ larger `K` ⇒ larger
   `M`, so estimation error is very costly and the controller must be de-tuned where
   sensing is poor — a control-first pipeline over-drives.
3. **Tight shared budget** (`B: 8 → 5`, `η: 0.4 → 0.6`). Fiercer competition between
   fast reserve and metering, so the order of spending matters more.
4. **Disturbance vs poor observability** (disturbed areas 1–5 get 8× worse tie
   metering). Where the control puts reserve is exactly where sensing is hard, so the
   control choice strongly changes where sensing is needed.
5. **Telemetry latency** (first-order lag on the tie channel, augmented state,
   `r_x = 48`). AGC acts on a delayed ACE; the controller must be designed around the
   sensing latency.
6. **Tie measurement bias** (random-walk bias per tie channel, augmented state,
   `r_x = 48`). Biased ACE causes integral-loop error; good metering removes a
   control-critical error the control-first design is blind to.

---

## 3. Results

Solved with `LQR_θ` + `filter_θ` (envelope gradient), multi-start (2 starts),
BLAS threads pinned to 1. `↓` = reduction vs the naive on-budget baseline;
**`joint − seq`** is the coupling gap (higher = more value from joint co-design).

| # | Change | `r_x` | grad err | `J_est ↓` | `J_det ↓` | `J_c ↓` | `J_tot ↓` | **`joint − seq`** |
|---|---|---|---|---|---|---|---|---|
| 0 | baseline A | 36 | 1.6e-5 | 11.8% | 7.5% | 8.1% | 7.8% | **4.66%** |
| 1 | noisy tie telemetry (`v0_tie ×7.5`) | 36 | 1.6e-5 | 11.8% | 7.5% | 8.1% | 7.8% | 4.67% |
| 2 | aggressive AGC (`R ×0.1`) | 36 | 4.6e-5 | 7.1% | 4.3% | 5.0% | 4.9% | 3.15% |
| 3 | tight shared budget (`B 8→5`, `η 0.4→0.6`) | 36 | 1.8e-5 | 4.8% | 3.1% | 3.4% | 3.0% | 0.03% |
| 4 | disturbance ⟂ poor observability | 36 | 1.6e-5 | 11.8% | 7.5% | 8.1% | 7.8% | 4.66% |
| 5 | telemetry latency (lag state) | 48 | 1.6e-5 | 11.8% | 7.5% | 8.1% | 7.8% | 4.67% |
| 6 | tie measurement bias (random walk) | 48 | 4.3e-3 | — | — | — | — | infeasible* |

\*Change 6 hit the filter-Riccati detectability guard (`J ≈ 1e8`) and the gradient
degraded (`4.3e-3`): an additive random-walk bias is (near-)unobservable through a
channel that also carries the dynamic tie flow.

---

## 4. Which changes helped

**None of the six changes increased the coupling gap.** The baseline `4.66%` stands;
two changes reduced it, three had essentially no effect, and one broke numerically.
The reasons are instructive and locate where A's coupling actually lives:

- **Tie-channel changes (1, 4, 5) did nothing** (gap and costs identical to baseline
  to 3 sig figs; change 5 even augments the state to `r_x=48` yet is unchanged).
  These all act *only on the tie-flow measurement*, and it turns out that channel is
  **nearly control-irrelevant** here: the control-induced weight `M(θ)=Kᵀ(R+BᵀPB)K`
  keys on the frequency states (what AGC prioritizes), and frequency is cheaply and
  well measured (fixed `v_f`), so tie-flow uncertainty barely gates control. De-noising,
  co-locating, or delaying the tie channel therefore moves nothing. This overturns the
  premise that "AGC acts on the estimated ACE" makes tie metering the critical axis —
  with frequency well measured, the ACE is effectively known regardless of tie noise.
- **Aggressive AGC (2)** *reduced* the gap (4.66→3.15%): shrinking `R` rescales the
  problem so control dominates proportionally, compressing the relative
  estimation-coupling.
- **Tight budget (3)** *collapsed* the gap (→0.03%): `B=5` overshot — the pool is so
  scarce that both axes sit pinned near the baseline with no room to reallocate, so
  joint ≈ sequential (the "too-tight" regime flagged earlier). A milder tightening
  might have helped; `B=5` did not.
- **Bias (6)** is infeasible as posed (unobservable bias → detectability guard).

**Conclusion (per the plan to keep only gap-improving changes): there is nothing to
keep — baseline A remains the best at ~4.66%.** The diagnosis says why: A's gap is
*limited by its sensing axis not being control-critical*. The genuinely effective
lever is therefore not any of the six above but a **reformulation that makes the
tie-flow/ACE the control-critical, poorly-observed quantity** — e.g. weight the ACE
component (`ΔP_tie`) far above frequency in `Q`, and/or degrade the frequency channel
too, so estimation quality actually gates achievable control. That is the natural next
step if we want a larger gap; it changes A's objective rather than tweaking a
parameter, so it is left as a proposal here rather than applied.
