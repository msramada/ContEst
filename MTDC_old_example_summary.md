# Old MTDC example — summary (paper + code)

Snapshot of the multi-terminal HVDC (MTDC) study **before** it is replaced by the
ℓ₁-SDP sensor-allocation-with-caps example. Captured so the original is documented
and recoverable.

---

## 1. Role in the paper

MTDC is the paper's **large-scale, linear, robust** study — the one case that
uses the **steady-state H∞ (game-Riccati) output-feedback** inner problem
(`src/Hinf.jl`) instead of the eKF–MPC path used by the three nonlinear studies
(ADCS, distillation, PLL). It is the scalability + robustness headline:
`r_x = 50` states, solved by an `O(r_x^3)` game Riccati.

Location in `ContEst_TeX/main.tex`:
- **Subsection** `\subsection{Multi-terminal HVDC droop coordination}`
  (`\label{sec:sim:mtdc}`), ~lines 1592–1739.
- **Allocation table** `\label{tab:mtdc-alloc}` (ℓ₁ sparsity/performance frontier).
- **Consolidated table** `\label{tab:formulations}` — the MTDC row
  `(50,25,25)`, `25+0` design (droop only), `H∞ output-fb.`,
  `J_est↓42.6% / J_det↓39.7% / J_tot↓43.4%`.
- **Summary text** `sec:sim:summary` — MTDC reported separately (worst-case
  guaranteed cost, different numerical scale).
- **Complexity section** `sec:sim:complexity` — MTDC is the exemplar of the
  `O(r_x^3)` game-Riccati route (≈2 ms at `r_x=50`), and the reason the H∞ SDP
  form is avoided (rank-deficient / ill-conditioned dual at the minimal-γ optimum).
- **Cross-references** in `sec:hinf` (Riccati twin) and `sec:sensor`
  (allocation exercised on the robust path).

## 2. Problem posed (paper)

Linearized device–network model, explicit-Euler discretized:

- `x_{k+1} = A(θ) x_k + B u_k + E w_k`, `y_k = C x_k + D_v v_k`,
  `z_k = [Q^{1/2}x; R^{1/2}u]`, with `A = I + Δt A_c`, `E = W^{1/2}`.
- **Only DC-bus voltages measured** (`C` selects `ΔV_i`, `r_y = 25`,
  `V = 10^{-3} I`); converter **currents estimated** ⇒ genuine output feedback.
- Fixed attenuation `γ² = 8`. Inner value splits (as in LQG separation):
  `J_c(θ) = tr(XW) [=J_det] + tr(MΣe) [=J_est]`, `M = Q`, `X` from the control
  game GARE and `Σe` from the **dual filter GARE**.
- Design gradient = sum of the two closed-form envelope contractions
  (`∂/∂A = 2XA_cl S` and its dual), chained with `∂A(θ)/∂θ` — no differentiation
  through the solver. Fixed-θ game separates into two Riccatis, but θ (droop)
  enters `A(θ)` and moves **both** X and Σe ⇒ co-design does not separate.

**Why H∞:** the failure mode is a *resonant peak* — converter lag makes excessive
droop provoke a lightly damped inter-terminal DC-voltage oscillation; a
renewable-infeed step near that resonance is amplified far beyond RMS. H∞ prices
the worst-case amplification; both `J_det` and `J_est` are sharply **U-shaped** in
the droop. A naive tight-regulation design over-provisions droop onto the
resonance wall (the baseline); ContEst lowers the gains to the interior optimum.

**State-space (per terminal `i`, ring of `n_t = 25`, `r_x = 50`, `Δt = 5 ms`):**
- `C_i ΔV̇_i = -Σ_{j~i} G_ij(ΔV_i-ΔV_j) + ΔI_i + infeed`
- `τ_i Δİ_i = -ΔI_i - (k_0 θ_i) ΔV_i + u_i`, `k_0 = 12`
- `y_{k,i} = ΔV_i + v_{k,i}`, `v ~ N(0, 10^{-3} I_25)`
- Heterogeneous renewable-infeed `σ_inf` per node; `Q` penalizes `ΔV`, `R = 2I`,
  `M = Q`; converter-stress design cost `J_des = c_k Σ_i θ_i`, `c_k = 5e-3`.
- Box `θ ∈ [0.5, 4]^25`; naive over-provisioned baseline `θ_nom = 3`.
  Interior robust optimum near `θ ≈ 1`; `J_c(θ=3) ≈ 1.7 J_c(θ=1)`.

## 3. Results reported (paper)

**Droop-only co-design (25 design vars, all 25 sensors on):**
- All 5 BFGS starts agree (empirical global optimum over Θ).
- Optimum droop profile `θ ∈ [0.66, 0.92]` (interior, none at a bound; noisier
  terminals keep higher droop).
- `J_est` ↓ **42.6%** (1.08 → 0.62), `J_det` ↓ **39.7%** (2.78 → 1.68),
  `J_c` ↓ **40.5%** (3.86 → 2.30), `J_tot` ↓ **43.4%** (4.24 → 2.40).
- Single control-side variable improves the estimation half purely by reshaping
  the shared dynamics `A(θ)` — the design-level coupling a place-then-control
  pipeline cannot exploit.

**Sparse sensing axis (ℓ₁ allocation), Table `tab:mtdc-alloc`** — sensor gains
`α_i ≥ 0`, `C(α) = diag(α) C`, penalty `λ_s ‖α‖₁` in `J_des`; joint droop+sensing:

| λ_s | sensors | J_est | J_det | J_c |
|-----|---------|-------|-------|-----|
| baseline | 25 | 1.08 | 2.78 | 3.86 |
| 0.00 | 25 | 0.35 | 1.67 | 2.02 |
| 0.03 | 15 | 0.80 | 1.67 | 2.47 |
| 0.06 | 12 | 1.01 | 1.67 | 2.68 |
| 0.10 |  6 | 1.27 | 1.67 | 2.94 |
| 0.20 |  2 | 1.53 | 1.66 | 3.20 |

- `J_det` holds at its robust optimum (~1.67) across rows ⇒ frontier isolates the
  price of sparsity; `J_est` degrades **gracefully** (network coupling makes most
  telemetry redundant): 12 sensors ≈ 25-sensor naive-baseline `J_est`; 6 sensors
  cost only ~17%.
- Headline: at `λ_s = 0.06`, thresholded to a **10-sensor** integer menu and
  re-solved, `J_est = 0.72`, `J_det = 1.67` — `J_c` down **38.2%** vs full
  telemetry while **removing 60% of the voltage sensors**, estimation in fact
  *better* than with all 25.

## 4. Code

**`benchmarks/example_mtdc.jl`** — builds the model and runs the report:
- Constants: `NT = 25` (⇒ `2NT = 50` states), `dt = 0.005`, ring neighbours
  `nb(i)`, `Cdc`, `τcnv`, `Gdc = 5`, `k0 = 12`, `γ² = 8`, heterogeneous `σinf`.
- `Qmat = diag(60 on ΔV, 0.2 on ΔI)`, `Rmat = 2I`, voltage-only measurement
  `Cmat` (`meas_rows = 2i-1`), `v_meas = 1e-3`.
- `model(θ) -> (Ad,Bd,W,Qmat,Rmat,Cmat,Vmeas)`; droop `k0·θ_i` enters `Ac[c,v]`.
- `J_des(θ) = c_droop Σθ`, `c_droop = 5e-3`; box `[0.5,4]^25`; `θ_nom = 3`.
- Evaluator `of_eval = Hinf_of_θ(2NT, NT, NT; γ²)`; report via `hinf_of_report`
  (gradient check + 5-start BFGS + cost/θ tables).

**`src/Hinf.jl`** — steady-state H∞ (game-Riccati) inner solver:
- `hinf_gare(A,B,E,Q,R,γ²)` — control GARE as an indefinite-weight DARE
  (`B̃=[B E]`, `R̃=diag(R,−γ²I)`), value `tr(XW)`; reuses `dare`/`dlyap` (LQR.jl).
- `hinf_filter_gare(A,C,W,V,Le,γ²)` — **dual** filter GARE, error covariance `Σe`.
- `Hinf_θ` (state-fb), `Hinf_of_θ` (output-fb: both GAREs, returns
  `J_c, ∇J_c, Ku, J_det=tr(XW), J_est=tr(QΣe)`), envelope gradients via the H₂
  identities on augmented data; `∂/∂θ` model Jacobians by ForwardDiff.
- Outer wrappers `hinf_contest_objective`, `hinf_report`, `hinf_of_report`.

**Run:** `julia --project benchmarks/example_mtdc.jl`

## 5. Key properties (what the old example demonstrated)

- The only **robust / H∞** and only **large-scale (50-state)** study.
- **Two-axis co-design**: control droop `θ_f = k` (in `A`) + sensing gains
  `θ_h = α` (in `C`), on one robust inner value.
- **Scalability headline**: `O(r_x^3)` game Riccati (~2 ms at `r_x=50`), the
  contrast against the `O(r_x^6)` interior-point path in `sec:sim:complexity`.
- Sensor allocation via ℓ₁ relaxation, pruning 60% of sensors at low cost.
