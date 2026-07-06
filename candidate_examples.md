# ContEst — Candidate Real-World Systems & Formulation Guide

**Purpose.** A working reference for selecting and formulating the application
examples in the *Control and Estimation Co-Design (ContEst)* paper. Each
candidate is mapped onto the bilevel problem (Problem 1) of the draft, with
explicit instructions for instantiating `f`, `h`, `θ`, `J_des`, and the
lower-level dual-control/LQG problem.

---

## 0. The ContEst discriminator (read first)

Standard **Control Co-Design (CCD)** only requires a design parameter `θ` that
enters the plant dynamics `f`. **ContEst additionally requires `θ` to enter the
output map `h`** — sensor placement, measurement configuration, sensor
precision — so the *same* design choice trades off against **both** closed-loop
control performance **and** estimator accuracy / observability, coupled through
the dual-control lower level.

A system is a good ContEst example only if **all four** hold:

1. **`θ` enters `h`** — sensor/measurement placement is a real, costly design decision.
2. **Real-time optimal estimation** — a Bayesian/Kalman filter produces the information state.
3. **Real-time optimal control** whose performance *depends on estimation quality*.
4. **Large-scale & significant** — important infrastructure or mission, many coupled loops.

### Generic instantiation template

For every example below, fill in the same slots:

```
State        x_k  ∈ R^{r_x}
Control      u_k  ∈ R^{r_u}
Output       y_k  ∈ R^{r_y}
Dynamics     x_{k+1} = f(x_k, u_k, w_k; θ)        # θ_f enters here  (CCD axis)
Output       y_k     = h(x_k, v_k; θ)             # θ_h enters here  (ESTIMATION axis)
Disturbance  w_k (process), v_k (measurement), p_0 (initial-state prior)

Design param θ = (θ_f, θ_h)                       # θ_h is what makes it ContEst
Design cost  J_des(θ)                             # must be nontrivial so trade-off is real
Constraints  θ ∈ Θ

Lower level  J_stoch(p_0^i; θ) = min_u (1/N) E { Σ ℓ(x_k,u_k) + ℓ^N(x_N) }
             policy is causal:  u_k = u_k(y_{0:k}, u_{0:k-1}, p_0)
             optimal control is a function of the filter density
                 p(x_k | y_{0:k}, u_{0:k-1}; θ)
```

**Upper-level gradient.** Use the **envelope theorem**: `dJ/dθ` needs the
sensitivity of the lower-level optimal value `J_stoch(θ)`, which equals the
*explicit* partial of the lower-level objective in `θ` at the optimal policy —
so you never differentiate through the controller/filter optimizer. In the LQG
case `J_stoch` is closed-form in the control and filter Riccati solutions,
giving **analytic gradients in `θ`** through the Riccati and Kalman-covariance
equations.

---

## 1. Wide-Area Monitoring & Control (WAMS) in bulk power grids — **recommended flagship**

**What / significance / scale.** Damping of low-frequency inter-area
oscillations in large interconnected transmission grids. Critical
infrastructure; real deployments exist (BPA / Sandia / Montana Tech wide-area
damping proof-of-concept with real-time PMU feedback). Scales from the 4-machine
Kundur testbed to the 16-machine, 68-bus NETS–NYPS benchmark.

**State / control / output.**
- `x_k`: generator rotor angles & speeds (linearized swing dynamics), plus
  exciter/PSS states.
- `u_k`: damping-control signals (PSS supplementary inputs and/or FACTS device
  set-points — STATCOM / TCSC).
- `y_k`: PMU measurements (bus voltage phasors, frequencies) at instrumented buses.

**`θ` split (the ContEst core).**
- `θ_h` (enters `h`): **PMU placement and measurement weighting** — a vector
  (binary relaxed to continuous) selecting instrumented buses and their
  precision. This sets inter-area mode observability and the achievable filter
  covariance.
- `θ_f` (enters `f`): **damping-actuator siting and tuning** — PSS locations/gains,
  FACTS device location and size.

**Design cost.** `J_des(θ)` = (per-PMU cost × number of PMUs) + actuator
hardware/siting cost. Constraints `Θ`: budget on PMU count, feasible actuator buses.

**Lower-level dual control.** Linearized swing dynamics with `w_k` = load
fluctuation (process noise), `v_k` = PMU noise. A Kalman filter forms
`p(x_k | y_{0:k}, …; θ)`; an LQG damping controller acts on the estimate.
`ℓ` penalizes inter-area oscillation energy (tie-line power / speed-difference
deviations) + control effort. Poor PMU placement → weak inter-area observability
→ degraded estimate → weak damping. ContEst jointly optimizes placement,
actuator tuning, and the implied LQG.

**Why ContEst, not CCD.** PMU placement is purely an `h`/observability decision;
CCD has no native slot for it. The damping controller can only act on modes the
chosen PMUs render observable, so estimation design and control performance are
inseparable.

**Julia tooling / benchmark.** `PowerSystems.jl` (network/machine data),
`PowerSimulationsDynamics.jl` (swing dynamics), `ControlSystems.jl` +
`MatrixEquations.jl` (LQG Riccati/Kalman solves), `Optim.jl` (continuous-relaxed
`θ`) or `JuMP` (combinatorial PMU placement). Start on the Kundur 2-area model
(ships as a standard case), then scale to 68-bus NETS–NYPS.

**Key references.** Garcia-Sanz 2019 (CCD); Aminifar et al. (optimal PMU
placement for observability); Zenelis & Wang (PMU-based wide-area damping,
arXiv:2108.01193); CERTS / BPA–Sandia wide-area damping proof-of-concept.

---

## 2. Spacecraft Attitude Determination & Control (ADCS / GNC)

**What / significance / scale.** Orientation estimation and control for every
space mission (JWST, Gaia, GRACE-FO). Tight, physically real coupling between
sensing and control; sub-arcsecond pointing missions (e.g. IRASSI) make the
trade-offs explicit. Smaller state dimension than a grid but the cleanest
control–estimation coupling.

**State / control / output.**
- `x_k`: attitude (quaternion / error-quaternion), angular rate, gyro/sensor biases.
- `u_k`: reaction-wheel / thruster torques.
- `y_k`: star-tracker, sun-sensor, gyro, magnetometer measurements.

**`θ` split.**
- `θ_h` (enters `h`): **sensor suite selection & placement** — which sensors,
  their mounting geometry, precision class, and misalignment budget. Sets the
  achievable EKF/MEKF covariance.
- `θ_f` (enters `f`): **actuator configuration** — reaction-wheel pyramid
  geometry, wheel sizing, thruster placement.

**Design cost.** `J_des(θ)` = sensor + actuator mass/power/$ budget. Constraints
`Θ`: redundancy requirements, mass/power caps, field-of-view geometry.

**Lower-level dual control.** Nonlinear attitude kinematics/dynamics; `w_k` =
disturbance torques (gravity-gradient, solar pressure, drag), `v_k` = sensor
noise. A **Multiplicative EKF (MEKF)** forms the information state; an optimal
(LQR/robust) attitude controller acts on it. `ℓ` penalizes pointing error +
control effort. Fine-pointing modes exhibit a genuine dual/"probing" effect.

**Why ContEst, not CCD.** Pointing accuracy is limited by the *estimator*, which
is set by sensor choice/placement (`θ_h`) — a decision CCD cannot represent.
Maps directly onto your paper-3 EKF/LQR notation.

**Julia tooling / benchmark.** `SatelliteToolbox.jl` (dynamics/environment),
custom MEKF, `ControlSystems.jl` for LQR. Benchmark on an IRASSI-like or
CubeSat-class ADCS spec.

**Key references.** Crassidis & Junkins, *Optimal Estimation of Dynamic Systems*;
Markley & Crassidis (attitude estimation, MEKF); IRASSI sub-arcsecond ADCS
(J. Guidance, Control, and Dynamics).

---

## 3. Chemical Process Plants — integrated design & control

**What / significance / scale.** Simultaneous plant design and control of
reactor/separation networks (CSTRs in series, distillation, the Tennessee
Eastman benchmark). Economically significant; large MINLP/dynamic-optimization
instances (thousands of variables). **Best research-gap story:** the field is
overwhelmingly design+control — estimation is the underdeveloped axis ContEst
adds. Strong tie to your CMU contacts (Biegler, Grossmann).

**State / control / output.**
- `x_k`: concentrations, temperatures, holdups, pressures.
- `u_k`: flow rates, heat duties, valve positions.
- `y_k`: available process measurements (temperatures, selected compositions /
  soft-sensor inputs).

**`θ` split.**
- `θ_h` (enters `h`): **measurement/sensor placement & selection** — which
  states are measured, online analyzer vs. soft sensor, sampling precision.
  (This is the dimension classic integrated design-and-control omits.)
- `θ_f` (enters `f`): **equipment sizing** — reactor volume, heat-exchanger area,
  holdup, recycle structure (the established CCD variables here).

**Design cost.** `J_des(θ)` = capital cost (equipment) + sensor/analyzer cost.
Constraints `Θ`: feasibility, safety margins, controllability constraints.

**Lower-level dual control.** Nonlinear process dynamics with `w_k` = feed/load
disturbances, `v_k` = measurement noise; EKF/moving-horizon estimator forms the
information state; (N)MPC or LQG controller acts on it. `ℓ` penalizes
setpoint/economic deviation + control effort, with probabilistic constraints for
purity/safety.

**Why ContEst, not CCD.** Adds sensor placement and estimator performance to the
established simultaneous design-and-control problem — directly extends the
Biegler/Grossmann line of work.

**Julia tooling / benchmark.** `ModelingToolkit.jl` / `DifferentialEquations.jl`
(process model), `JuMP` + `Ipopt` (dynamic optimization), custom MHE/EKF.
Benchmark on CSTR-in-series or a reduced Tennessee Eastman model.

**Key references.** Biegler & Grossmann (retrospective on optimization);
Flores-Tlacuahuac & Biegler (simultaneous design & control of CSTRs);
Ricardez-Sandoval et al. (Tennessee Eastman simultaneous design & control).

---

## 4. Floating Offshore Wind Turbines / Wind Farms

**What / significance / scale.** The native home of the CCD lineage (Garcia-Sanz;
NREL's WEIS toolset does CCD of floating offshore wind turbines). Significant and
well-tooled. **Caveat:** estimation is a thinner part of the loop than in #1–#2,
so this reads as "CCD with an estimation extension" rather than a native ContEst
showcase — good as a secondary/validation example.

**State / control / output.**
- `x_k`: rotor speed, blade-flap/structural modes, platform/tower motion.
- `u_k`: blade pitch, generator torque, trailing-edge flap commands.
- `y_k`: available turbine sensors + (optionally) LIDAR wind preview.

**`θ` split.**
- `θ_h` (enters `h`): **LIDAR / load-sensor configuration** — preview range,
  measurement set for wind/load estimation.
- `θ_f` (enters `f`): **blade & flap geometry, actuator configuration, rotor
  diameter, floater parameters** (the established WEIS CCD variables).

**Design cost.** `J_des(θ)` ≈ levelized cost of energy (LCOE) components +
sensor cost. Constraints `Θ`: structural loads (damage-equivalent loads),
deflection limits, stability.

**Lower-level dual control.** Aero-servo-elastic dynamics with `w_k` = turbulent
wind, `v_k` = sensor noise; a wind/load estimator (Kalman-type) feeds a
controller (e.g. ROSCO-style). `ℓ` penalizes load/deflection + power-tracking +
actuator effort.

**Why ContEst, not CCD.** Adds the wind/load-estimation and LIDAR-placement axis
on top of an established CCD pipeline.

**Julia tooling / benchmark.** Couple to OpenFAST or build a reduced
aero-servo-elastic model in `ModelingToolkit.jl`; `Optim.jl`/`JuMP` for the
upper level. Benchmark against an NREL reference turbine.

**Key references.** Garcia-Sanz 2019 (CCD); NREL WEIS toolset; Abbas et al.
(aero-servo-elastic co-optimization, *Wind Energy*, 2023); Feil et al. (control
co-design of a floating offshore wind turbine, *Applied Energy*, 2024).

---

## 5. Comparison & recommendation

| System | `θ_h` (estimation axis) | Estimation–control coupling | Scale / significance | Julia readiness | Fit to ContEst |
|---|---|---|---|---|---|
| **1. WAMS / power grid** | PMU placement & weighting | Very strong | Critical infrastructure, large | High | **Best** |
| **2. Spacecraft ADCS** | Sensor suite selection/placement | Very strong | Every space mission | Medium–High | Excellent |
| **3. Chemical plant** | Measurement/sensor placement | Strong | High economic value | Medium | Good (best research gap) |
| **4. Floating wind** | LIDAR / load-sensor config | Moderate | High, well-tooled | Medium | Good (CCD-leaning) |

**Recommendation.**
- **Flagship:** WAMS / PMU + wide-area damping. Strongest `θ`-in-`h` coupling,
  large-scale, and directly Julia-able. Develop it at two scales — the Kundur
  2-area model as the transparent illustrative case, then the 68-bus NETS–NYPS
  as the benchmark redesign showing ContEst beating sequential
  (PMU-then-controller) design.
- **Second example:** Spacecraft ADCS, to demonstrate generality and exercise
  your paper-3 EKF/LQR notation.
- **Framing example for related work:** Chemical-plant integrated design &
  control, to position ContEst as the estimation-aware extension of the
  Biegler/Grossmann program.

---

## 6. Framework references (your draft)

- Garcia-Sanz (2019), *Control co-design: an engineering game changer.*
- Feldbaum (1960), *Dual control theory.*
- Bar-Shalom & Tse (1974), *Dual effect, certainty equivalence, and separation.*
- Kumar & Varaiya, *Stochastic Systems: Estimation, Identification, and Adaptive Control.*
- Anderson & Moore, *Optimal Filtering.*
- Envelope theorem (convex optimization) — for the upper-level gradient `dJ/dθ`.
