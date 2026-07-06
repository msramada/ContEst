# `example_wams.jl` — construction summary

ContEst flagship: **wide-area damping co-design with budgeted, binary PMU placement.**

## Problem
- **Network:** 4-machine, two-area Kundur-style grid. States `x=[Δδ₁..₄, Δω₁..₄]` (`n=8`), one damping actuator per generator (`m=4`), one candidate PMU per bus streaming (angle, speed) (`o=8`). Nonlinear swing dynamics (`sin` synchronizing terms), deviation form so `f(0,0,θ)=0` for all `θ`; forward Euler, `dt=0.05`, horizon `N=12`.
- **Design `θ` (8 params):** `θ_f = e₁..₄ ∈ [0.3,3]` actuator effectivenesses (enter `f`); `θ_h = α₁..₄ ∈ {0,1}` **binary PMU deployment gates** (enter `h`, `h = diag(α)·h_base`).
- **Disturbance concentrated in area 2** (machines 3,4): large initial state kick and large prior covariance `Σ_ic` there, area 1 quiet; low uniform process noise. ⇒ PMUs in area 2 are far more valuable, so *placement* drives the estimation cost.

## Method
- Inner solver: information-state eKF + convex MPC (`nonlinear_mpc_θ`); design gradients from the dynamics costates (envelope theorem). Outer loop: projected BFGS (`bfgs_design`).
- **Sparse sensor allocation:** `J_des = λ_e·Σ(eᵢ−1)² + λ_s·‖α‖₁`. Sweep `λ_s` from 0 → 30; lock the binary selection at the first `λ_s` meeting the budget `K_max=2` (top-`K_max` gates → 1, rest → 0), then re-optimize actuators with the gates fixed.
- Estimation cost computed consistently with the MPC via linearized `A_info` covariance propagation (so `J_est ≤ J_c`, `J_det ≥ 0`).
- **Comparison:** ContEst (ℓ1-placed PMUs + co-tuned actuators) vs an equal-budget **sequential baseline** (PMUs at the convenient buses {1,2}, nominal actuators).

## Key iterations (what was tried and why)
1. Started from a 2-machine continuous-precision PMU example → enlarged to 4 machines with **binary 0/1** gates + ℓ1 budget (gain >1 disallowed: can't boost SNR beyond hardware).
2. **6 machines was too slow** (nested-AD over a 90-dim info-state) → settled on 4 machines.
3. Symmetric network ⇒ all gates fell together / `J_est` insensitive to placement → fixed by **concentrating the disturbance + prior covariance in one area** (the minimal, decisive change).
4. Fixed a negative-`J_det` artifact by switching `J_est` to the MPC-consistent linearized rollout.
5. Robustness/speed: capped BFGS iterations, reduced `Σ_ic` contrast to avoid QP ill-conditioning, and lock-selection-then-continue so the sweep reaches 30 without stalling.

## Result (gradient check rel.err ≈ 6e-6)
ℓ1 sweep deploys PMUs at the disturbed buses **{3,4}**; quiet-area gates → 0 by `λ_s≈12`, disturbed-area gates persist near the threshold to `λ_s=30`.

| cost | sequential (PMUs {1,2}) | ContEst (PMUs {3,4}) | change |
|---|---|---|---|
| `J_est` | 56.7 | 19.0 | **−66%** |
| `J_det` | 38.5 | 29.4 | −24% |
| `J_c`   | 95.1 | 48.4 | −49% |
| `J_tot` | 111.1 | 68.4 | **−38%** |

Reported in the paper as §9.1 + Fig. (ℓ1 selection + cost breakdown) + summary table.
