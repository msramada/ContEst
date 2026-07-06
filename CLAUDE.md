# ContEst Co-Design

Research project on **ContEst**: bringing state/parameter estimation into control
co-design by jointly optimizing actuation (`θ_f`, enters `f`) and sensing
(`θ_h`, enters `h`) parameters against a single stochastic objective.

## Key locations

- **`Paper_LaTeX/`** — the paper I am working on. `automatica.tex` is the main
  manuscript (Automatica format, `autart.cls`); figures live in
  `Paper_LaTeX/figs/`, references in `Paper_LaTeX/References.bib`.
- **`References/`** — the background material required for the paper. Consult it
  for prior work, definitions, and framing before writing or citing.
- **`GUIDE.md`** — the guideline that must be followed when creating the
  simulations/examples. Follow it for problem structure, cost decomposition
  (`J_est`, `J_det`, `J_c`, `J_des`, `J_tot`), and reporting conventions.

## Simulations

- Examples are the `benchmarks/example_*.jl` files. The ones **included in the
  paper** are four nonlinear (eKF–MPC) studies — `example_wams.jl` (grid:
  comms-budgeted PMU rate via `V(θ)` + damping provisioning),
  `example_adcs.jl` (spacecraft: shared power/mass budget split across
  sensors+wheels via `V(θ)`), `example_distillation.jl` (feed-stage/sensor-tray
  placement), and `example_pll.jl` (PLL/DSE estimator co-design, bias–variance
  `V(θ)`) — plus one **large-scale linear LQG** study, `example_mtdc.jl`
  (multi-terminal HVDC droop coordination, control-only `θ_f`, U-shaped `J_c`),
  which uses the steady-state LQG lower level (control Riccati + Lyapunov) rather
  than eKF–MPC.
  Every study is posed so its optimum is a genuine **interior / non-monotone**
  trade (budgeted allocation, saturating benefit, or U-shape), not a parameter
  pinned at a bound. (`example_cstr.jl` was removed as trivial; the earlier LQG
  suite `example_inertia.jl`/`example_pss.jl`/`example_bess.jl`/`example_sensing.jl`
  was removed, keeping only the MTDC study.)
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
  (`wams` / `adcs` / `distillation` / `pll`)
- Build the paper: from `Paper_LaTeX/`, run `pdflatex automatica.tex` twice
  (the `.bbl` is committed; `latexmk` may fail on pre-existing bibliography
  issues, so prefer direct `pdflatex` passes).

## Conventions

- Keep the paper's numbers, tables, and figures **consistent with the code** —
  after changing an example, re-run it and update the corresponding LaTeX
  table/figure/text.
- Match the existing style of each example and the paper subsections when adding
  new ones.
