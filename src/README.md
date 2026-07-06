# `src/` — ContEst source modules

ContEst ("**Cont**rol–**Est**imation co-design") solves a **bilevel** design
problem: an *outer* loop picks physical/sensor **design parameters** `θ` while
an *inner* loop solves the closed-loop control+estimation problem and returns
both its optimal cost and the analytic gradient `∂J/∂θ`. The outer loop is an
ordinary gradient-based optimiser (BFGS); the inner loop is what these modules
implement.

```
upper level:   min_{θ ∈ [θ_lb, θ_ub]}   J(θ) = J_c(θ) + J_des(θ)
                                          │            │
                                          │            └─ design/hardware cost (you define)
                                          └─ inner control cost from the MPC (MPCs.jl)
```

There are **two independent inner solvers**:

| file        | role                                  | used by              |
|-------------|---------------------------------------|----------------------|
| `MPCs.jl`   | finite-horizon **information-state MPC** (a QP) — inner solver | all `example*.jl` |
| `SDPs.jl`   | infinite-horizon **steady-state SDP/LMI** — alternative inner solver | standalone (`test_SDP.ipynb`) |
| `eKF.jl`    | the extended Kalman filter both inner solvers build on | `MPCs.jl`           |
| `BFGS.jl`   | the **outer loop**: cached objective, gradient check, BFGS / multi-start | all `example*.jl` |

`ContEst.jl` simply `include`s all four.

The key idea shared by both solvers: the optimal value `J(θ)` is a function of
matrices that depend on `θ` (the linearised dynamics `A`, `B`, the weights
`Q`, `R`, the process noise `W`). The **envelope theorem** says the gradient of
the optimal value w.r.t. those matrices equals the optimiser's **dual
variables**, so `∂J/∂θ` is obtained by contracting the duals with `∂(matrix)/∂θ`
(computed by `ForwardDiff`) — no need to differentiate through the solver.

---

## `eKF.jl` — extended Kalman filter / information-state propagation

The plant is nonlinear:

```
x_{k+1} = f(x_k, u_k, θ) + w_k,   w ~ N(0, W)
y_k     = h(x_k, u_k, θ) + v_k,   v ~ N(0, V)
```

At design time we do not have real measurements, so we propagate the *expected*
estimate and its **covariance** `Σ` forward — the pair `(x, Σ)` is the
**information state**. The covariance is what couples sensor quality into the
control cost (via the `tr(Q Σ)` term in the running cost).

**Functions**

- `∇ₓf, ∇ᵤf, ∇ₓh, ∇ᵤh` — `ForwardDiff` Jacobians `∂f/∂x`, `∂f/∂u`, `∂h/∂x`,
  `∂h/∂u` evaluated at an operating point. These are the local linearisations
  `F = ∂f/∂x` and `H = ∂h/∂x` the EKF needs.

- `time_update((x,Σ), u, θ, f, Covars)` — the **prediction** step:
  ```
  F  = ∂f/∂x|(x,u,θ)
  x⁺ = f(x, u, θ)
  Σ⁺ = F Σ Fᵀ + W
  ```

- `measurement_update((x,Σ), u, y, θ, h, Covars; mode)` — the **correction** step:
  ```
  H = ∂h/∂x|(x,u,θ)
  L = Σ Hᵀ (H Σ Hᵀ + V + 1e-4·I)⁻¹      # Kalman gain  (ε·I conditions the inverse)
  Σ⁺ = (I − L H) Σ
  ```
  With `mode = "predict"` (the design-time default) the **mean is left
  unchanged** — there is no real innovation `y − h(x)` at design time, we only
  need how the covariance *would* shrink. With any other `mode` the usual mean
  correction `x⁺ = x + L(y − h(x))` is applied.

- `update(...)` — one full EKF step: `time_update` then `measurement_update`.
  This is the single building block `MPCs.jl` differentiates to get the
  information-state dynamics.

**Why it matters for co-design:** a *better* sensor (larger output sensitivity
in `h`, encoded in `θ`) makes `H` larger, which makes `L` larger and `Σ⁺`
smaller — directly lowering the `tr(Q Σ)` estimation-uncertainty cost.

---

## `MPCs.jl` — information-state MPC (the inner solver used by the examples)

This is the workhorse. It lifts the EKF into a **finite-horizon QP** whose
optimal value is the inner control cost `J_c(θ)`, and returns `∂J_c/∂θ`.

### `ControlModel` (the problem definition you fill in)

```julia
ControlModel(
    f, h,            # nonlinear dynamics & measurement, signatures f(x,u,θ), h(x,u,θ)
    (W, V),          # process- and measurement-noise covariances
    n, m, o,         # state / input / output dimensions
    N,               # prediction horizon (number of steps)
    runningcost,     # deterministic ℓ(x,u)   (reference cost)
    stoch_runningcost, # stochastic ℓ̄(ζ,u) = ‖x‖²_Q + ‖u‖²_R + tr(Q Σ)
    constraint_function # g(x,u) ≤ 0  (input/state inequality constraints)
)
```

### Information-state packing

The covariance is symmetric, so only its upper triangle is stored. The
information state is the flat vector `ζ = [x ; vech(Σ)]` of length
`n_info = n + n(n+1)/2`.

- `mat_uptriang_to_vec(Σ)` → `vech(Σ)` (column-stacked upper triangle).
- `vec_to_mat(v)` → the symmetric matrix (used inside cost terms like `tr(Q Σ)`).
- `info_state_dynamics(ζ, u, θ, f, h, Covars)` → one EKF `update` written as a
  map `ζ_{k+1} = Φ(ζ_k, u_k, θ)` on the flat info-state.

### `nonlinear_mpc_θ(prob) → eval`

Builds the JuMP/Clarabel model **once** and returns a closure `eval` that is
called repeatedly by the outer optimiser. The QP is

```
min_{ζ, u}   Σ_{t=1}^{N-1} ℓ̄(ζ_t, u_t)  +  ℓ̄(ζ_N, 0)
s.t.         ζ_1 = ζ₀                              (initial info-state)
             Σ_1 ⪰ 1e-6·I                          (PSD on initial covariance)
             ζ_{t+1} = A_info ζ_t + B_info u_t      (linearised info dynamics)
             g(x_{t+1}, u_t) ≤ 0                    (inequality constraints)
```

**`eval(x₀₀, Σ₀₀, u₀, θ)`** does, on each call:

1. **Linearise** the information-state dynamics about the operating point
   `(ζ₀, u₀)` using `ForwardDiff`:
   `A_info = ∂Φ/∂ζ`, `B_info = ∂Φ/∂u`. Both depend on `θ`.
2. **Rebuild** the dynamics + initial-state constraints with these matrices and
   `optimize!`.
3. Read off:
   - `J` — the optimal QP value = `J_c(θ)`.
   - `u_first = value.(u)[:,1]` — the first receding-horizon action.
   - `feasibility` — solver status.
4. **Gradient** via the envelope theorem. Let `λ_t` be the dual of the
   `t`-th dynamics constraint. Aggregate
   ```
   Λ_A = Σ_t λ_t ζ_tᵀ ,   Λ_B = Σ_t λ_t u_tᵀ
   ```
   then, with `∂A_info/∂θ` and `∂B_info/∂θ` from (nested) `ForwardDiff`,
   ```
   ∂J_c/∂θ_k = tr(Λ_Aᵀ ∂A_info/∂θ_k) + tr(Λ_Bᵀ ∂B_info/∂θ_k).
   ```

Returns `(u_first, J, ∇θJ, feasibility)`.

> **Important subtlety (fixed):** the parameter gradient in step 4 must
> re-linearise about the *same* operating input `u₀` that the QP in step 1 used.
> An earlier version overwrote `u₀` with the optimised first input before
> computing `∂A/∂θ, ∂B/∂θ`, which silently biased the gradient whenever the
> optimal input was non-zero (e.g. an active flight controller). The first
> receding-horizon action is now returned in a separate variable `u_first`,
> leaving `u₀` intact. You can always verify the gradient against finite
> differences — every `example*.jl` does this and prints the relative error
> (it should be `< 5e-2`, typically `~1e-5`).

---

## `SDPs.jl` — steady-state SDP/LMI design (standalone alternative)

An infinite-horizon alternative to the MPC: instead of a finite-horizon QP it
solves a single **semidefinite program** for the stationary closed loop, and
returns the optimal cost, its `θ`-gradient, and the stabilising gain
`K = L Σ⁻¹`. Demonstrated in `test_SDP.ipynb`; not used by `example*.jl`.

`SDP_θ(rₓ, rᵤ) → eval!` builds the model once and returns a closure
`eval!(myLinearModel, θ)`, where `myLinearModel(θ) → (A, B, W, Q, R)` returns
the linear model matrices at `θ`.

**The SDP** (decision variables `Σ` ⪰ 0 state covariance, `Z₀` input-cost
auxiliary, `L = K Σ`):

```
min_{Σ, Z₀, L}   tr(Q Σ) + tr(R Z₀)
s.t.   Σ ⪰ 1e-5·I
       [ Z₀   L ;  Lᵀ  Σ ] ⪰ 0                 # Schur:  Z₀ ⪰ L Σ⁻¹ Lᵀ  (input cost)
       [ Σ−W   AΣ+BL ; (AΣ+BL)ᵀ  Σ ] ⪰ 0       # Lyapunov LMI (closed-loop stability)
```

**Gradient** from the dual `S` of the Lyapunov LMI (partitioned into blocks
`S11, S12`):

```
∂J/∂A = −2 S12 Σᵀ ,  ∂J/∂B = −2 S12 Lᵀ ,  ∂J/∂W = S11 ,  ∂J/∂Q = Σ ,  ∂J/∂R = Z₀
```

chained with the `ForwardDiff` Jacobians of `(A,B,W,Q,R)` w.r.t. `θ` to form
`∂J/∂θ`. The closed-loop gain is recovered as `K = L Σ⁻¹` (`Lstar / Σstar`).

### Known issues (flagged during review — fix before relying on this module)

1. **Wrong function differentiated.** The Jacobian calls reference a global
   `myModel(θ)` instead of the closure argument `myLinearModel`. It happens to
   work in `test_SDP.ipynb` only because a global named `myModel` exists there;
   pass any other model and the gradient is wrong or errors. Replace `myModel`
   with `myLinearModel` in the five `ForwardDiff.jacobian` calls.
2. **`Σ` and `Z₀` are not symmetric variables.** They are declared as full
   square matrices and only their upper triangles enter the `Symmetric(...)`
   PSD constraints, leaving the lower-triangle entries **free**. They still
   appear in `tr(Q Σ)`, so for a non-diagonal `Q` the problem is unbounded.
   It is safe only because the notebook uses `Q = I`. Declare them
   `@variable(model, Σ[1:rₓ,1:rₓ], Symmetric)` (and likewise `Z₀`).

---

## `BFGS.jl` — the outer optimisation loop

The inner solver returns `J_c(θ)` and `∂J_c/∂θ`; this module turns that into an
optimal design. Examples **call** these helpers instead of inlining the
cached-objective / gradient-check / multi-start boilerplate.

- `cached_objective(eval_J) → (f, g!)` — wrap any `eval_J(θ) → (J, ∇J)` into the
  `(f, g!)` pair `Optim` expects, caching the last `θ` so the underlying solve
  runs **once** per point (Optim queries value and gradient separately).

- `contest_objective(mpc_eval, x0, Σ0, u_lin, J_des, ∇J_des) → (f, g!)` — the
  standard ContEst objective `J = J_c + J_des` (and its gradient), already
  cached. This is what an example normally calls.

- `verify_gradient(f, g!, θ) → (∇_analytic, ∇_fd, rel_err)` — check the analytic
  gradient against finite differences. Run this before trusting the optimiser;
  `rel_err` should be `< 5e-2` (typically `~1e-5`).

- `bfgs_design(f, g!, θ_lb, θ_ub, θ_init; …) → Optim result` — one
  box-constrained `Fminbox(BFGS())` solve using the analytic gradient. Extra
  keywords pass through to `Optim.Options`.

- `multistart_design(f, g!, θ_lb, θ_ub; n_starts, θ_nom, seed, …)
  → (best_J, best_θ, results)` — run `bfgs_design` from `θ_nom` plus random
  starts inside the box and keep the best. Because `J(θ)` is generally
  non-convex, this escapes local minima; if every start returns the same point
  it certifies the optimum as global.

---

## How the pieces connect

```
 θ ──► example*.jl: ControlModel(f, h, …)  ─────────────► MPCs.nonlinear_mpc_θ
                                                                │
   BFGS.jl outer loop  (bfgs_design / multistart_design)       │ each call:
   J(θ)=J_c+J_des,  ∇J=∇J_c+∇J_des  ◄────── (J_c, ∇J_c) ───────┤  linearise via
   via contest_objective                                        │  eKF.update → QP
                                                                ▼
                                                          Clarabel solve
                                                          + dual-based ∇θ
```

See **[`../GUIDE.md`](../GUIDE.md)** for a step-by-step recipe to build a new
example, run BFGS with multiple initial points, and read off the cost
improvement.
