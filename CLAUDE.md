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
  `V(θ)`) — plus one **large-scale linear, robust output-feedback H∞** study,
  `example_mtdc.jl` (multi-terminal HVDC droop coordination, control-side droop `θ_f`,
  U-shaped worst-case costs), which uses the steady-state **H∞ game Riccati** lower
  level rather than eKF–MPC (`src/Hinf.jl`, fixed attenuation `γ²=8`). It solves TWO
  indefinite-weight DAREs: a control GARE (`B̃=[B E]`, `R̃=diag(R,−γ²I)`, `E=√W`) giving
  `J_det=tr(XW)`, and its DUAL filter GARE (`B̃=[Cᵀ Leᵀ]`, `R̃=diag(V,−γ²I)`, `Le=√Q`)
  giving `J_est=tr(QΣe)` — DC-bus **voltages measured**, converter currents estimated.
  Both reuse `dare`/`dlyap` from `src/LQR.jl` with the same envelope gradients
  (`∂/∂A=2XA_cl S` and its dual `2Σe Ã_cl S̃`); the droop `θ_f` enters `A`, coupling
  both halves at the design level. As `γ→∞` these collapse to the H₂ `tr(PW)`/Kalman.
  (The earlier H₂/LQG version was replaced.) Evaluator/report: `Hinf_of_θ`,
  `hinf_of_report`.
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
  the information-state MPC. Only reporting changed, not `best_θ`. MTDC (LQG path)
  keeps its exact steady-state value.

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
