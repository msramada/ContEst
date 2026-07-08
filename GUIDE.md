# ContEst — running the solver on a new example

This guide shows how to (1) define a new co-design problem, (2) wire it into the
information-state MPC, (3) optimise the design parameters `θ` with BFGS, (4) use
**multiple initial points** to escape local minima / certify the optimum, and
(5) read off the cost improvement. The module internals are documented in
[`src/README.md`](src/README.md).

The recipe is system-agnostic: every term below is something *you* supply for
your own plant, sensors, and design parameters. Copy the structure verbatim and
fill in the six ingredients.

> **Paper draft.** The [`ContEst_TeX/`](ContEst_TeX/) folder holds the draft of
> the accompanying paper — the manuscript I am preparing for submission to the
> *Automatica* journal. It is the authoritative reference for the theory behind
> this code (the bilevel co-design formulation, the envelope-theorem gradients,
> and the SDP / MPC-eKF inner solvers); consult it for derivations and notation.

---

## Pick the inner solver first — SDP vs MPC-eKF

ContEst ships **two** inner solvers, and the choice is dictated by your problem,
not by preference. Decide this before writing any code:

> - **Linear dynamics *and* no constraints → use the SDP path (`SDPs.jl`).**
>   The infinite-horizon steady-state LMI is exact for a linear, unconstrained
>   plant, returns the stabilising gain `K = L Σ⁻¹` in one solve, and has no
>   horizon to tune. Do **not** reach for the MPC here — it would only be an
>   approximate, finite-horizon, re-linearised version of a problem the SDP
>   solves exactly.
> - **Nonlinear dynamics *and/or* any input/state constraints → use the
>   MPC-eKF path (`MPCs.jl`).** This is the only path that can (a) re-linearise
>   the eKF/dynamics about the operating point each step and (b) enforce
>   inequality constraints `g(x,u) ≤ 0`. The SDP path supports neither.

| your problem is…                    | inner solver            | entry point             |
|-------------------------------------|-------------------------|-------------------------|
| linear **and** unconstrained        | **SDP / LMI**           | `SDP_θ` (`SDPs.jl`)     |
| nonlinear, **or** has constraints   | **MPC + eKF**           | `nonlinear_mpc_θ` (`MPCs.jl`) |

Rule of thumb: **the moment you have a single nonlinearity or a single
constraint, you are on the MPC-eKF path.** The rest of this guide (§0–§6) is
written for that path, since it is the general case. If your plant is genuinely
linear and unconstrained, follow the SDP usage in `SDPs.jl` /
[`src/README.md`](src/README.md) instead — the outer BFGS / multi-start loop
(§4–§5) and the cost-reporting tables (§6) apply unchanged; only the inner
`mpc_eval` call is replaced by the SDP `eval!`.

---

## 0. Setup (once)

```bash
julia --project -e 'using Pkg; Pkg.instantiate()'   # first run only
julia --project your_example.jl                      # run your script
```

## 1. The mental model

ContEst minimises a **total** design cost over physical/sensor parameters `θ`:

```
min_{θ ∈ [θ_lb, θ_ub]}   J(θ) = J_c(θ) + J_des(θ)
```

- `J_c(θ)` — closed-loop control+estimation cost, computed by the MPC
  (`MPCs.nonlinear_mpc_θ`). You never code this; it comes out of the solver.
- `J_des(θ)` — your **design/hardware cost**: what each upgrade *costs*
  (bigger actuator, finer sensor, plant retuning …). You define it and its
  gradient.

`θ` typically splits in two: parameters that enter the **dynamics** `f`
(actuator authority, plant tuning → better *control*) and parameters that enter
the **measurement** `h` (sensor sensitivities → better *estimation*, via the
`tr(Q Σ)` term).

## 2. Define the problem (six ingredients)

```julia
include("src/eKF.jl")
include("src/MPCs.jl")
include("src/BFGS.jl")    # cached objective + gradient check + BFGS / multi-start
using LinearAlgebra, Optim, Printf

# (a) nonlinear dynamics  x_{k+1} = f(x, u, θ)
function f_sys(x, u, θ) ...; return x_next end

# (b) parametric measurement  y = h(x, u, θ)   (θ here tunes sensors)
h_sys(x, ::AbstractVector, θ) = [θ[3]*x[1], θ[4]*x[2], ...]

# (c) noise covariances  (process W, measurement V)
const Covars = (W, V)

# (d) costs: deterministic ℓ(x,u) and stochastic ℓ̄(ζ,u) = ‖x‖²_Q+‖u‖²_R+tr(QΣ)
det_cost(x, u)   = (x-s_x)'*Q*(x-s_x) + (u-s_u)'*R*(u-s_u)
function stoch_cost(ζ, u)
    n = <state dim>
    xk = ζ[1:n]; Σk = vec_to_mat(ζ[n+1:end])
    return (xk-s_x)'*Q*(xk-s_x) + (u-s_u)'*R*(u-s_u) + tr(Q*Σk)
end

# (e) input/state inequality constraints  g(x,u) ≤ 0
constraint_fn(::AbstractVector, u) = [u[1]-umax; -u[1]-umax; ...]

# (f) assemble and compile the MPC evaluator
prob = ControlModel(f_sys, h_sys, Covars, n, m, o, N, det_cost, stoch_cost, constraint_fn)
const mpc_eval = nonlinear_mpc_θ(prob)
```

> Tips that save debugging time:
> - Build the dynamics so that `x = 0, u = 0` is an equilibrium for **every**
>   `θ` (e.g. write force/moment terms about the trim point). This keeps the
>   regulation target consistent as `θ` varies.
> - Penalise small-magnitude channels more heavily in `Q` than large ones so
>   all states are regulated comparably regardless of their physical units.

## 3. The upper-level objective + design cost

```julia
const θ_nom = ones(length(θ_init))          # baseline design
const λ_des = [...]                          # cost weight per parameter
J_des(θ)  = sum(λ_des .* (θ .- θ_nom).^2)    # quadratic upgrade cost
∇J_des(θ) = 2 .* λ_des .* (θ .- θ_nom)
```

`contest_objective` (from `src/BFGS.jl`) builds the **cached** `f`/`g!` pair —
`J(θ) = J_c(θ) + J_des(θ)` and its gradient — so BFGS does not re-solve the QP
when it asks for the value and gradient at the same `θ`:

```julia
contest_f, contest_g! = contest_objective(mpc_eval, x_ic, Σ_ic, u_lin, J_des, ∇J_des)
```

**Always verify the analytic gradient** against finite differences before
trusting the optimiser. `verify_gradient` returns `(∇_analytic, ∇_fd, rel_err)`;
the relative error should be `< 5e-2` (typically `~1e-5`):

```julia
∇_env, ∇_fd, rel_err = verify_gradient(contest_f, contest_g!, θ_init)
@show rel_err
```

## 4. Optimise with BFGS (single start)

`bfgs_design` wraps `Optim.Fminbox(BFGS())` — it handles the box bounds and uses
the analytic gradient. Extra keywords pass through to `Optim.Options`:

```julia
result = bfgs_design(contest_f, contest_g!, θ_lb, θ_ub, θ_init;
                     iterations = 80, g_tol = 1e-4)   # show_trace = true to watch it
θ_opt = result.minimizer
```

`result` is a standard `Optim` result, so `Optim.minimum(result)`,
`Optim.converged(result)`, and `Optim.iterations(result)` all apply.

## 5. Multiple initial points (multi-start) — do not skip this

`J(θ)` is **non-convex** in `θ` (the dynamics enter `f`/`h` nonlinearly), so
BFGS only ever returns a **local** minimum — the one in whose basin `θ_init`
happens to fall. A single run can therefore report a design that is *not* the
best available, with no warning. The remedy is to **start from several different
initial parameter points** and let each descend independently: the more of the
box you seed, the lower the chance that every run is trapped in the same
sub-optimal basin.

`multistart_design` does this for you — it runs `bfgs_design` from `θ_nom` plus
`n_starts-1` random points drawn uniformly inside `[θ_lb, θ_ub]`, then returns
the single best result:

```julia
best_J, best_θ, results = multistart_design(
    contest_f, contest_g!, θ_lb, θ_ub;
    n_starts = 8, θ_nom = θ_nom, seed = 20240624,   # seed → reproducible starts
    g_tol = 1e-5, iterations = 80,
)
```

Use **enough** starts to cover the box (8–16 is typical for a 5-parameter
design; scale up with `length(θ)`).

**Reading the outcome.** `results` holds one `(; x0, θ, J, iters)` per start, so
you can see how many distinct minima appeared:
- If all starts return the **same** `θ*`/`J*` → strong evidence the optimum is
  global.
- If starts return **different** minima → the problem is genuinely multimodal;
  `best_θ`/`best_J` already hold the lowest-cost one.

> **Report only the best minimum.** The optimal design you present in §6 is
> `best_θ` — the single local minimum with the lowest total cost. The other
> starts exist only to *certify* that choice (they show nothing better was
> found); their individual minimisers are intermediate diagnostics and are
> **omitted** from the final results. It is good practice to state how many
> starts were run and how many distinct minima they found.

## 6. Report the cost improvement — break the cost into its components

Reporting only the single number `J_tot` hides *where* a design change acts. The
stochastic running cost the MPC minimises is itself a sum of two physically
distinct pieces, and the design cost sits on top:

```
            ┌ deterministic control cost  J_det = Σ_t ‖x_t‖²_Q + ‖u_t‖²_R   (regulate the mean state + effort)
J_c  ───────┤                                                                ← "stochastic control cost (both)"
(stochastic)└ estimation cost            J_est = Σ_t tr(Q Σ_t)              (penalty on estimation uncertainty)

J_tot = J_c + J_des ,   J_des = hardware/design cost you pay for the upgrade
```

So there are **five** quantities worth tabulating, and the interesting result is
how much each one *moves* between the baseline design `θ_nom` and the optimum
`best_θ`:

| component | symbol | what it measures | typically… |
|-----------|--------|------------------|------------|
| Estimation cost            | `J_est` | `Σ tr(Q Σ_t)` — expected uncertainty | **decreases** (better sensors / observability) |
| Deterministic control cost | `J_det` | `Σ ‖x‖²_Q + ‖u‖²_R` — mean tracking + effort | decreases (more control authority) |
| Stochastic control cost (both) | `J_c` | `J_det + J_est` — the MPC's objective | decreases |
| Design cost                | `J_des` | hardware paid for the upgrade | **increases** (you bought something) |
| **Total**                  | `J_tot` | `J_c + J_des` — the bilevel objective | net **decrease** |

### Computing the breakdown

The **reported** comparison comes from **sampled closed-loop trajectories of the
original cost**, not the open-loop MPC prediction. `simulate_mc` (`src/Simulate.jl`)
runs the true nonlinear plant forward under a receding-horizon controller with the
eKF on the noisy measurements, and averages the realised cost over `n_samples`
noise/IC draws. **Both** components use the true state:

```
J_det = Σ_t ‖x_true_t − s_x‖²_Q + ‖u_t − s_u‖²_R     (prob.runningcost, true state)
J_est = Σ_t (x_true_t − x_est_t)' Q (x_true_t − x_est_t)   (realised estimation error)
J_c   = J_det + J_est
```

The **baseline** (`:baseline`) is driven by the certainty-equivalence controller
`certainty_equivalence_mpc` — it uses **only the eKF mean** (mean dynamics,
deterministic cost, blind to the covariance). The **optimum** (`:optimum`) is
driven by the information-state MPC (`mpc_eval(...; grad=false)`):

```julia
ctrl_info(xe, Σe, θ) = mpc_eval(xe, Σe, u_lin, θ; grad = false)[1]  # info-state MPC
ctrl_mean(xe, Σe, θ) = ce_mpc(xe, u_lin, θ)                          # mean-only MPC

function cost_breakdown(θ, role::Symbol)
    ctrl = role == :baseline ? ctrl_mean : ctrl_info
    r  = simulate_mc(prob, ctrl, θ, Matrix(Q), x_ic, Σ_ic)   # sampled closed loop
    Jc = r.J_det + r.J_est
    return (; J_est = r.J_est, J_det = r.J_det, J_c = Jc,
              J_des = J_des(θ), J_tot = Jc + J_des(θ))
end
```

`run_contest_report` calls `cost_breakdown(θ_nom, :baseline)` and
`cost_breakdown(best_θ, :optimum)` (it detects the 2-arg form; legacy 1-arg
breakdowns still work). The **design optimisation is unchanged** — multi-start BFGS
still minimises the open-loop MPC objective with the exact envelope gradient
(`contest_objective`), so `best_θ` and the parameter table do not depend on the
sampling; only the reported cost split does.

> A common seed is reused across the baseline and optimum runs, so the comparison
> is paired (common random numbers) and the *difference* is low-variance even at a
> few hundred samples. Bump `n_samples` if a reported percentage sits within a
> point or two of the Monte-Carlo noise.

### Presenting it

Tabulate every component with its **change** (`Δ = optimal − baseline`) and
percentage, plus a short design-parameter table. Print with, e.g.:

```julia
# % change is undefined when the baseline is 0 (e.g. J_des(θ_nom) = 0 by
# construction) — show a dash there and report the absolute Δ instead.
pct(b, o_) = abs(b) < 1e-9 ? "     —  " : @sprintf("%+6.2f%%", 100*(o_-b)/b)
for (name, b, o_) in (("J_est", base.J_est, opt.J_est),
                      ("J_det", base.J_det, opt.J_det),
                      ("J_c  ", base.J_c,   opt.J_c),
                      ("J_des", base.J_des, opt.J_des),
                      ("J_tot", base.J_tot, opt.J_tot))
    @printf("%s  %10.4f → %10.4f   Δ = %+9.4f   %s\n", name, b, o_, o_-b, pct(b, o_))
end
```

**Cost-component table** (the headline result):

| cost component                 | baseline `θ_nom` | optimal `θ*` | Δ (increase +) | % change |
|--------------------------------|------------------|--------------|----------------|----------|
| Estimation cost `J_est`        | …                | …            | …              | …        |
| Deterministic control `J_det`  | …                | …            | …              | …        |
| Stochastic control (both) `J_c`| …                | …            | …              | …        |
| Design cost `J_des`            | …                | …            | …              | …        |
| **Total `J_tot`**              | …                | …            | …              | …        |

**Design-parameter table** (one row per `θ`, best minimum only):

| parameter | baseline `θ_nom` | optimal `θ*` | enters | role |
|-----------|------------------|--------------|--------|------|
| `θ[1]` …  | …                | …            | `f`/`h`| control authority / sensor |

**How to read it:**
- The **net result** is the `J_tot` row — a reduction in the bilevel objective.
- Read it as a *trade*: `J_des` **goes up** (you bought hardware) while `J_c`
  **comes down**; the net is what the co-design actually buys.
- The `J_est` vs `J_det` rows show *which mechanism* paid off — a falling
  `J_est` means the win came from better **estimation** (sensor `θ` in `h`); a
  falling `J_det` means it came from better **control** (plant/actuator `θ` in
  `f`). Sensor channels only help when they are heavily weighted in `Q`.
- These results are for the **best minimum only** (§5); the other multi-start
  minimisers are not reported.
- If `J_tot` barely moves, the cost is usually dominated by an unavoidable
  transient-recovery / estimation-uncertainty floor — the design knobs only trim
  the *controllable* margin on top of it.

---

## Checklist for a new example

1. Implement `f`, `h`, `W`/`V`, `Q`/`R`, the running costs, and `constraint_fn`.
2. Assemble `ControlModel` → `nonlinear_mpc_θ` → `mpc_eval`.
3. Define `J_des`/`∇J_des`, then build the objective with `contest_objective`.
4. **Verify the gradient** with `verify_gradient` (`rel_err < 5e-2`).
5. `bfgs_design` from `θ_init`; then `multistart_design` (§5).
6. Tabulate the baseline-vs-optimal costs (§6) and report the reduction.

All four optimisation helpers — `contest_objective`, `verify_gradient`,
`bfgs_design`, `multistart_design` — live in [`src/BFGS.jl`](src/BFGS.jl); the
example only calls them.
