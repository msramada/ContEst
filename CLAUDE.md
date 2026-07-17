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

- Examples are the `benchmarks/example_*.jl` files. The ones **included in the
  paper** are three nonlinear (eKF–MPC) studies —
  `example_adcs.jl` (spacecraft: shared power/mass budget split across
  sensors+wheels via `V(θ)`), `example_distillation.jl` (feed-stage/sensor-tray
  placement), and `example_pll.jl` (PLL/DSE estimator co-design, bias–variance
  `V(θ)`) — plus one **linear, robust (H∞) two-axis** study,
  `example_mtdc.jl` (multi-terminal HVDC, 30 states / `NT=15` terminals): a MERGED
  co-design of the **control droop `θ_f=k`** (enters `A`) via the **H∞ game Riccati**
  (`src/Hinf.jl` `hinf_gare`, fixed `γ²=8`, `J_det=tr(XW)`, U-shaped worst-case cost)
  AND the **sensing gains `θ_h=α`** (enter `C`) via an **ℓ1 covariance-SDP**
  (JuMP/Clarabel): posterior error-cov bound `X ⪰ (Σ_pred(θ)⁻¹ + Σⱼ αⱼ CⱼᵀCⱼ/v)⁻¹`
  with `Σ_pred(θ)=dlyap(A(θ),W)`, sparse ℓ1 sensor selection, and **per-converter-current
  estimation caps `[X]_{2i,2i} ≤ tol`** (`J_est=tr(QX)`). The caps + discrete selection
  make the estimation half **SDP-only** (no Riccati expresses per-state covariance caps
  with an active set); the cap duals are **shadow prices** localizing essential sensors
  (terminal 11), with an **infeasibility floor** below the all-on covariance. The droop
  enters `A`, so it reshapes BOTH `X` (game Riccati) and `Σ_pred(θ)` (estimation prior),
  coupling the two axes. Report `mtdc_report()` prints: (A) droop co-design (all 15
  sensors: `J_det −38.2%`, `J_est −60.8%`, `J_tot −49.7%`), (B) ℓ1 frontier, (C) capped
  ℓ1-SDP + shadow prices, (D) merged optimum (4 of 15 sensors `{4,7,11,12}`, `J_c −40.9%`).
  Droop step: `Optim` `Fminbox(LBFGS())` on `Jdet+Jest+c_k Σθ`; sensor step: `sensor_sdp`
  (Clarabel). (The earlier 50-state droop-only `Hinf_of_θ`/`hinf_of_report` version — a
  dual filter GARE `J_est=tr(QΣe)` with no sensor allocation — was replaced; see
  `MTDC_old_example_summary.md`. `Hinf_of_θ`/`hinf_filter_gare` remain in `src/Hinf.jl`
  but are no longer used by the example.)
  Every study is posed so its optimum is a genuine **interior / non-monotone**
  trade (budgeted allocation, saturating benefit, or U-shape), not a parameter
  pinned at a bound. (`example_cstr.jl` was removed as trivial; the earlier LQG
  suite `example_inertia.jl`/`example_pss.jl`/`example_bess.jl`/`example_sensing.jl`
  was removed, keeping only the MTDC study. The two-area WAMS study
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
  the information-state MPC. Only reporting changed, not `best_θ`. MTDC keeps its
  exact steady-state values (H∞ game-Riccati control + covariance-SDP estimation).

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
