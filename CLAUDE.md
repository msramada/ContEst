# ContEst Co-Design

Research project on **ContEst**: bringing state/parameter estimation into control
co-design by jointly optimizing actuation (`θ_f`, enters `f`) and sensing
(`θ_h`, enters `h`) parameters against a single stochastic objective.

## Key locations

- **`ContEst_TeX/`** — the paper I am working on. `main.tex` is the main
  manuscript (Automatica format, `autart.cls`); references in
  `ContEst_TeX/References.bib`.
- **`References/`** — the background material required for the paper. Consult it
  for prior work, definitions, and framing before writing or citing.
- **`GUIDE.md`** — the guideline that must be followed when creating the
  simulations/examples. Follow it for problem structure, cost decomposition
  (`J_est`, `J_det`, `J_c`, `J_des`, `J_tot`), and reporting conventions.

## Simulations

- Examples are the `benchmarks/example_*.jl` files. Four studies are **included
  in the paper**: two nonlinear (eKF–MPC) studies — `example_distillation.jl`
  (feed-stage/sensor-tray placement) and `example_pll.jl` (PLL/DSE estimator
  co-design, bias–variance `V(θ)`) — a linear steady-state SDP study,
  `example_adcs_hinfest.jl` (spacecraft ADCS: `H₂` control under a hard
  control-effort-covariance cap `+` worst-case `H∞` robust estimation, shared
  power/mass budget across wheel authority `e_rw` and sensor precisions
  `α_st,α_gyro`; reported costs are Monte-Carlo realizations, M=100/T=30, of
  the fixed linear gains on the true nonlinear plant — see below), and
  `example_mtdc_sdp.jl` (multi-terminal HVDC, 30 states / `NT=15` ring
  terminals, **last** in Section 9): a MERGED co-design of the **droop
  `θ_f=k`** (enters `A`) and **sparse sensor gains `θ_h=α∈[0,1]`** (enter `C`),
  **SDP path only — no Riccati equation anywhere**. Control is a
  fixed-`γ²=20` `H∞` bounded-real-lemma SDP restricted to a **block-diagonal**
  (fully decentralized) pattern — the tighter **banded** (ring-neighbour)
  restriction is numerically **dual-degenerate** (confirmed: >100% disagreement
  vs. a central finite difference, sign included) so it cannot drive BFGS, a
  finding worth checking before reusing this file's control solver at a
  different `γ²` or pattern. Estimation is the `H₂` SDP in the
  observability-gramian (`P^ε`,`F`) form, restricted to the **banded**
  pattern (no degeneracy there), with an `ℓ1` penalty `λ_s‖α‖₁` for sparse
  sensing. The gradient check uses **central**, not forward, finite
  differences (`h=1e-3`) and tightened Clarabel tolerances
  (`tol_gap_abs=tol_gap_rel=tol_feas=1e-10`) on the control SDP — the default
  tolerance leaves ~1e-5 absolute solver noise that swamps the true local
  derivative at practical step sizes (confirmed empirically this session).
  `mtdc_sdp_report()` prints: SDP validity (structured vs. full-dense, no
  Riccati/Kalman reference), the gradient check (0.55% relative error, 30
  params), a multi-start joint co-design (droop `θ*∈[0.89,1.67]`, 12/15
  sensors kept, `J_est −31.8%`, `J_cont −33.4%`, `J_tot −21.6%`), and the
  headline **sensor-count frontier**: fixed-droop pruning vs. co-designed
  droop at each retained-sensor count, showing the co-designed droop with 6
  sensors matches the pruned baseline's 15-sensor estimation cost (60% fewer
  sensors, no loss in `J_est`) — droop coordination "buys back" voltage
  sensors, the paper's headline for this study. Older MTDC files
  (`example_mtdc.jl` base library, `example_mtdc_sparse.jl`,
  `example_mtdc_alloc.jl`) remain as historical/dev exploration — same pattern
  as ADCS's superseded `example_adcs.jl`/`_hinf.jl`/`_h2cap.jl` — and describe
  an earlier Riccati-based or differently-patterned MTDC design no longer used
  by the paper.
  Every study is posed so its optimum is a genuine **interior / non-monotone**
  trade (budgeted allocation, saturating benefit, or U-shape), not a parameter
  pinned at a bound. (`example_cstr.jl` was removed as trivial; the earlier LQG
  suite `example_inertia.jl`/`example_pss.jl`/`example_bess.jl`/`example_sensing.jl`
  was removed. The two-area WAMS study
  (`example_wams.jl`) and all WAMS experiments were removed: on both the LQG and
  the constrained eKF–MPC paths the co-design gate was ~0–1% over a competent
  sequential design and the joint optimum degraded inter-area damping — WAMS is
  sensor-dominated / near-separable, so it is not a defensible co-design headline.)
- Shared machinery is in `src/` (eKF, MPC, BFGS, `Simulate.jl`) and
  `report_contest.jl` (gradient check + multi-start BFGS + cost/parameter tables).
- Design optimization uses **multi-start BFGS with 5 initial points**: the
  baseline `θ_nom` plus 4 uniform random restarts (`n_starts = 5`), keeping the
  best minimizer. The optimization is still on the open-loop MPC value with the
  exact envelope gradient — **unchanged**.
- **Reported** eKF–MPC costs (`J_est`/`J_det`) come from **sampled closed-loop
  trajectories of the original cost** (`simulate_mc` in `src/Simulate.jl`), both
  components on the true state; the **baseline is a certainty-equivalence
  controller (eKF mean only)** via `certainty_equivalence_mpc`, the ContEst optimum
  the information-state MPC. Only reporting changed, not `best_θ`. MTDC and ADCS
  are on the SDP path, not eKF-MPC: MTDC reports exact steady-state SDP/Lyapunov
  values (no Monte Carlo — its plant is linear by construction, so there is no
  surrogate-vs-truth gap to check), while ADCS reports Monte-Carlo realizations
  of its fixed linear gains on the true (mildly nonlinear) plant.

## Commands

- Run an example's report: `julia --project benchmarks/example_<name>.jl`
- Regenerate a paper figure: `julia --project make_paper_figs.jl <case>`
  (`adcs` / `distillation` / `pll`)
- Build the paper: from `ContEst_TeX/`, run `pdflatex main.tex` twice
  (the `.bbl` is committed; `latexmk` may fail on pre-existing bibliography
  issues, so prefer direct `pdflatex` passes).

## Conventions

- Keep the paper's numbers, tables, and figures **consistent with the code** —
  after changing an example, re-run it and update the corresponding LaTeX
  table/figure/text.
- Match the existing style of each example and the paper subsections when adding
  new ones.
