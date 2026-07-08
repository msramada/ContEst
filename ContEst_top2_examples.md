# Two Strongest Candidate ContEst Examples — Formulation, Results, and Outlook

This document records, in full detail, the two cross-field examples that came out
strongest in the ContEst simulation sweep: **cybergenetics** (synthetic biology)
and **seismic protection of a building** (structural / civil engineering). For
each it gives the problem formulation *as actually implemented and solved*, the
state-space equations, the cost functions, the ContEst inner method used, the
numerical results (all with envelope/analytic gradients verified against finite
differences), and concrete development directions.

Common conventions follow the paper: continuous-time dynamics
$\dot x = A_c x + B_c u + E w$ are discretized to
$x_{k+1}=A_\theta x_k + B_\theta u_k + w_k$, $y_k = C_\theta x_k + v_k$, with
$w_k\sim\mathcal N(0,W)$, $v_k\sim\mathcal N(0,V(\theta))$, $Q\succeq0$, $R\succ0$.
The design vector splits as $\theta=(\theta_f,\theta_h)$: $\theta_f$ enters the
actuation/plant ($B_\theta$, and possibly $A_\theta$), $\theta_h$ enters the
sensing/estimation ($C_\theta$, the measurement covariance $V(\theta)$, or the
filter state itself). The reported costs are
$J_{\mathrm{det}}=\operatorname{tr}(P W)$ (control),
$J_{\mathrm{est}}=\operatorname{tr}(Q\Sigma_e)$ (estimation),
$J_{\mathrm{des}}$ (design/budget), and $J_{\mathrm{tot}}=J_{\mathrm{det}}+\alpha
J_{\mathrm{est}}+J_{\mathrm{des}}$ with $\alpha=1$. The **co-design gate** compares
the *joint* optimum against a *competent sequential* design (tune sensing
$\theta_h$ for estimation, then actuation $\theta_f$ for control).

---

## 1. Cybergenetics — optogenetic control of gene expression (synthetic biology)

### 1.1 Field and significance

Cybergenetics closes a real-time feedback loop around a living cell: light drives
a light-inducible promoter, a fluorescent reporter is measured, and a controller
regulates protein expression. The feature that makes this a *ContEst* problem
rather than a pure control problem is that **the reporter is itself a dynamic,
noisy sensor with a maturation lag**: newly translated fluorophore is dark and
matures at rate $k_{\mathrm{mat}}$, so fluorescence is a *delayed, filtered*
readout of the true protein — not a direct measurement. The reporter design
(fusion strength, maturation rate, brightness) therefore sets both the
measurement noise $V(\theta)$ *and the sensor dynamics themselves*. Designing the
circuit (promoter/RBS strengths — actuation) together with the reporter (sensing)
is a genuine co-design, and it is the cleanest demonstration of the paper's
distinctive claim that **the estimator is a first-class design variable**.

### 1.2 State-space model

State $x=[m,\,p,\,r_1,\,r_2]^\top$: target mRNA $m$, protein $p$, immature
(dark) reporter $r_1$, and mature (fluorescent) reporter $r_2$, all in deviation
form about the operating set-point. Light input $u$ induces transcription; the
Hill induction is linearized about the set-point (a standard cybergenetics
small-signal step), so the implemented model is linear:

$$
\begin{aligned}
\dot m   &= \theta_{f,\alpha}\,\alpha\,u \;-\; \delta_m\,m,\\
\dot p   &= \theta_{f,\beta}\,\beta\,m \;-\; \delta_p\,p,\\
\dot r_1 &= \theta_{h,\beta_r}\,\beta_r\,m \;-\; \big(\delta_r+\theta_{h,\mathrm{mat}}\,k_{\mathrm{mat}}\big)\,r_1,\\
\dot r_2 &= \theta_{h,\mathrm{mat}}\,k_{\mathrm{mat}}\,r_1 \;-\; \delta_r\,r_2.
\end{aligned}
$$

In matrix form $\dot x = A(\theta)x + B(\theta)u$ with

$$
A(\theta)=
\begin{bmatrix}
-\delta_m & 0 & 0 & 0\\
\theta_{f,\beta}\beta & -\delta_p & 0 & 0\\
\theta_{h,\beta_r}\beta_r & 0 & -(\delta_r+\theta_{h,\mathrm{mat}}k_{\mathrm{mat}}) & 0\\
0 & 0 & \theta_{h,\mathrm{mat}}k_{\mathrm{mat}} & -\delta_r
\end{bmatrix},
\qquad
B(\theta)=\begin{bmatrix}\theta_{f,\alpha}\alpha\\0\\0\\0\end{bmatrix}.
$$

The key structural feature: **$\theta_h=(\theta_{h,\beta_r},\theta_{h,\mathrm{mat}})$
enters $A(\theta)$** (the reporter sub-dynamics), i.e. the *sensing design augments
the filter state* — a coupling none of the paper's current studies exhibit.
$\theta_{f,\beta}$ (translation) also enters $A$, so control and estimation share
the dynamics: choosing $\theta_{f,\beta}$ changes both what the controller can do
and what the estimator can infer. Discretization is forward Euler at
$\Delta t=0.25$ ($A_\theta=I+\Delta t\,A(\theta)$, $B_\theta=\Delta t\,B(\theta)$;
the gene time-constants are slow, $\Delta t\,\delta\ll1$).

Constants used: $\delta_m=0.2,\ \delta_p=0.1,\ \delta_r=0.1,\ k_{\mathrm{mat}}=0.3,\
\alpha=\beta=\beta_r=\phi=1$.

### 1.3 Output equation

Fluorescence reads the **mature** reporter only, scaled by brightness $\phi$, with
photon/shot noise whose variance falls with brightness:

$$
y_k = \theta_{h,\phi}\,\phi\,r_{2,k} + v_k,
\qquad C(\theta)=[\,0\ \ 0\ \ 0\ \ \theta_{h,\phi}\phi\,],
\qquad v_k\sim\mathcal N\!\big(0,\,v_0/\theta_{h,\phi}\big),\ v_0=0.05 .
$$

Because only $r_2$ is measured, the true protein $p$ must be reconstructed through
the reporter's maturation lag — a genuine deconvolution.

### 1.4 Cost functions

Process noise $W=\operatorname{diag}(2{\times}10^{-3},2{\times}10^{-3},10^{-4},10^{-4})$
(expression/biological noise). Control weight $Q=\operatorname{diag}(0,10,0,0)$
(regulate protein $p$), $R=0.2$.

$$
J_{\mathrm{det}}=\operatorname{tr}(P W),\qquad
J_{\mathrm{est}}=\operatorname{tr}(Q\,\Sigma_e),\qquad
J_{\mathrm{des}}(\theta)=c_b\Big(\textstyle\sum_i\theta_i-B\Big)^2 .
$$

$\theta=(\theta_{f,\alpha},\theta_{f,\beta},\theta_{h,\beta_r},\theta_{h,\mathrm{mat}},\theta_{h,\phi})$,
box $[0.3,4]^5$, shared "part-library / burden" budget $B=5$, $c_b=3$
(baseline $\theta_{\mathrm{nom}}=\mathbf 1$). The maturation choice
$\theta_{h,\mathrm{mat}}$ is a **bias–variance knob** (fast-dim vs. slow-bright) —
the biological analog of the PLL-bandwidth study.

### 1.5 ContEst method

**H2 / LQG-SDP path** (Sec. 5): the inner value is the stationary LQG cost,
$J_{\mathrm{det}}=\operatorname{tr}(PW)$ from the control DARE and
$J_{\mathrm{est}}=\operatorname{tr}(Q\Sigma_e)$ from the dual Kalman filter DARE,
with the exact envelope gradient. Because $\theta$ enters $A(\theta)$, the filter
evaluator differentiates through $A(\theta)$ as well as $C(\theta),V(\theta)$.
Optimization is multi-start box-constrained BFGS; the co-design gate is joint vs.
the sequential (estimator-then-controller) pipeline. Gradient verified vs. finite
differences to relative error $4.5\times10^{-5}$.

### 1.6 Results

| Quantity | Baseline | ContEst joint | Change |
|---|---:|---:|---:|
| Estimation cost $J_{\mathrm{est}}=\operatorname{tr}(Q\Sigma_e)$ | 1.137 | 0.451 | **−60.4%** |
| Control cost $J_{\mathrm{det}}=\operatorname{tr}(PW)$ | 0.070 | 0.085 | +21.6% (worse) |
| Total $J_{\mathrm{tot}}$ | 1.208 | 0.536 | **−55.6%** |

- **Co-design gate:** sequential $J_{\mathrm{tot}}=0.785$ → joint $0.536$ = **+31.7%** (joint beats a competent sequential design by ~32%).
- **Design:** $\theta_f=(\theta_\alpha,\theta_\beta):(1,1)\to(1.25,0.30)$; $\theta_h=(\theta_{\beta_r},\theta_{\mathrm{mat}},\theta_\phi):(1,1,1)\to(1.55,1.21,0.69)$.
- Optimum is at a **corner** (translation strength $\theta_\beta$ floored).

Interpretation: this is an *estimation-dominated* regime — the protein is only
observed through the lagged reporter, so the binding cost is inference. ContEst
invests heavily in the reporter (raise fusion $\beta_r$ and maturation, trim
brightness), cutting the estimation cost 60% and the total 56%, and it does so in
a way a place-then-control pipeline cannot (the reporter dynamics couple into the
achievable control), hence the large +32% co-design gate.

### 1.7 Possible developments

1. **Restore the Hill nonlinearity on the eKF–MPC path.** The implemented model
   linearizes the Hill induction $u^n/(K^n+u^n)$. Solving the full nonlinear
   circuit with the eKF–MPC inner problem would (a) match the paper's other
   nonlinear studies, (b) capture the genuinely nonlinear light-dose response and
   phototoxicity bound $0\le u\le u_{\max}$, and (c) let the bias–variance
   maturation trade play out over a transient set-point change rather than a
   stationary operating point.
2. **Address the corner.** $\theta_\beta$ flooring means translation is being
   traded entirely for reporter quality. Re-scaling the metabolic-burden prices
   (separate $c_{\mathrm{burden}}$ on $\theta_\alpha,\theta_\beta$ vs.
   $c_{\mathrm{tox}}$ on brightness) or an $\ell_1$ part-selection term would move
   the optimum interior and answer the discrete "which part from the library"
   question in one convex sweep.
3. **Multi-gene / multi-reporter circuits.** Extending to two target genes with
   distinct reporters would create a genuine *reporter-allocation* problem
   (which species to fuse the bright/slow vs. dim/fast reporter to), a richer
   co-design than the single-reporter case.
4. **Validated biological model.** Replace the normalized constants with
   characterized part parameters (promoter/RBS strengths from a real library,
   measured maturation rates for common fluorophores) so the numbers are
   experimentally traceable — the prerequisite for a synthetic-biology venue.

---

## 2. Seismic protection of a building (structural / civil engineering)

### 2.1 Field and significance

Active/semi-active control of civil structures under earthquake loading is a
mature field and the natural home for the paper's LQG-SDP machinery: the plant is
linear, the cost is quadratic (interstory drift and acceleration), and both
actuator (brace/damper) authority and sensor (accelerometer/displacement)
precision are classically posed as placement problems. The estimation axis is
essential because a building is instrumented on only a few floors, the full
drift/velocity state is needed for control, and the ground motion is unmeasured.
It provides a hard, interpretable *physical* metric — peak interstory-drift RMS,
the standard proxy for structural damage.

### 2.2 State-space model

A 3-storey shear building, mass-normalized (floor mass $=1$), with story
stiffness $\kappa=280$, mass $M=I_3$, stiffness and damping matrices

$$
K=\kappa\begin{bmatrix}2&-1&0\\-1&2&-1\\0&-1&1\end{bmatrix},\qquad
C_d=0.005\,K+0.02\,M\ (\approx 2\%\ \text{modal damping}).
$$

With floor-displacement vector $q\in\mathbb R^3$ and state
$x=[q;\dot q]\in\mathbb R^6$, ground acceleration $\ddot x_g$ (disturbance) and
per-story active-brace forces $u\in\mathbb R^3$,

$$
\dot x=
\underbrace{\begin{bmatrix}0&I\\-M^{-1}K&-M^{-1}C_d\end{bmatrix}}_{A_c}x
+\underbrace{\begin{bmatrix}0\\M^{-1}\Gamma\,\theta_{f}\end{bmatrix}}_{B_c(\theta_f)}u
+\underbrace{\begin{bmatrix}0\\-\iota\end{bmatrix}}_{E}\ddot x_g,
\qquad
\Gamma=\begin{bmatrix}1&-1&0\\0&1&-1\\0&0&1\end{bmatrix},
$$

where $\Gamma$ maps interstory brace forces to floors and $\theta_f$ is the
(scalar) brace authority scaling all three actuators. The model is discretized
**exactly (zero-order hold)** at $\Delta t=0.01$ (forward Euler is unstable for
the lightly-damped structural modes). Ground motion enters as rank-1 process
noise on the floor velocities, $W=\sigma^2\,e\,e^\top+10^{-8}I$ with
$e=[0,0,0,1,1,1]^\top$, $\sigma^2=1$.

Interstory drifts are $d=T_d\,q$ with
$T_d=\begin{bsmallmatrix}1&0&0\\-1&1&0\\0&-1&1\end{bsmallmatrix}$, extracted from
the state by $T_{\mathrm{drift}}=[\,T_d\ \ 0\,]$.

### 2.3 Output equation

Floor-displacement sensors on all three floors, with per-floor precision
$\theta_h=(\alpha_1,\alpha_2,\alpha_3)$:

$$
y_k=C\,x_k+v_k,\qquad C=[\,I_3\ \ 0\,],\qquad
v_k\sim\mathcal N\!\big(0,\,V(\theta_h)\big),\quad
V(\theta_h)=\operatorname{diag}\!\big(v_0/\alpha_j^2\big),\ v_0=10^{-2}.
$$

### 2.4 Cost functions

Control weight penalizes interstory drift (damage) heavily over velocity:

$$
Q=\begin{bmatrix}T_d^\top(10^{3} I_3)T_d & 0\\ 0 & I_3\end{bmatrix},\qquad R=10^{-2}I_3 .
$$

$$
J_{\mathrm{det}}=\operatorname{tr}(P W),\quad
J_{\mathrm{est}}=\operatorname{tr}(Q\Sigma_e),\quad
J_{\mathrm{des}}(\theta)=c_b\Big(\textstyle\sum_i\theta_i-B\Big)^2 .
$$

Design vector $\theta=(\theta_f,\alpha_1,\alpha_2,\alpha_3)$, box $[0.3,5]^4$,
shared retrofit budget $B=4$ ($\theta_{\mathrm{nom}}=\mathbf 1$), $c_b=6$. A
physical **damage metric** is also reported: the closed-loop RMS interstory drift
$\sqrt{\operatorname{diag}(T_{\mathrm{drift}}\,S\,T_{\mathrm{drift}}^\top)}$, with
$S=\mathrm{dlyap}(A_\theta+B_\theta K,\,W)$ the stationary covariance under the LQR
gain $K$.

### 2.5 ContEst method

**H2 / LQG-SDP path** (Sec. 5), the ideal vehicle for the paper's exact-gradient
regularity result on a linear plant: $J_{\mathrm{det}}=\operatorname{tr}(PW)$ from
the control DARE, $J_{\mathrm{est}}=\operatorname{tr}(Q\Sigma_e)$ from the dual
Kalman DARE, exact envelope gradient, multi-start BFGS, and the joint-vs-sequential
co-design gate. Here $\theta_f$ enters $B_\theta$ and $\theta_h$ enters
$C,V$; the coupling that makes co-design non-trivial is the shared retrofit
budget. Gradient verified vs. finite differences to relative error $10^{-5}$.

### 2.6 Results

| Quantity | Baseline | ContEst joint | Change |
|---|---:|---:|---:|
| Estimation cost $J_{\mathrm{est}}$ | 15.877 | 21.100 | +32.9% (worse) |
| Control cost $J_{\mathrm{det}}$ | 121.967 | 31.410 | **−74.2%** |
| Total $J_{\mathrm{tot}}$ | 137.844 | 57.410 | **−58.4%** |
| **Peak interstory-drift RMS** (damage) | 0.152 | 0.060 | **−60.3%** |

- **Co-design gate:** sequential $J_{\mathrm{tot}}=73.93$ → joint $57.41$ = **+22.4%**.
- **Design:** $\theta_f:1\to3.31$ (brace authority raised); $\theta_h=(\alpha_1,\alpha_2,\alpha_3):(1,1,1)\to(0.30,0.99,0.30)$.
- Optimum is at a **corner** (two floor sensors floored — the co-design keeps only the informative mid-floor sensor).

Interpretation: seismic drift control is *actuation-dominated* — reducing drift
is fundamentally an actuation task — so ContEst pours the shared budget into brace
authority and keeps only the sensor that most improves the closed-loop drift
estimate. The headline is the **physical −60% peak-drift reduction**; the +22%
co-design gate shows the separation pipeline mis-allocates the budget (it
over-provisions estimation-optimal sensing that the closed loop does not need).

### 2.7 Possible developments

1. **Drift chance-constraint (breaks separation cleanly).** Replace the
   quadratic drift cost with a probabilistic damage limit
   $\Pr(|d_i|\le d_{\max})\ge1-\varepsilon$, i.e.
   $\kappa\,\sigma_{d_i}(\theta)\le d_{\max}$. Because the true closed-loop drift
   covariance depends on *both* the LQR gain and the Kalman error covariance
   (via an augmented Lyapunov equation), the back-off couples sensing into control
   feasibility: starving the sensors inflates drift uncertainty and violates the
   code limit, so sensing can no longer be floored. This is the principled route
   to an **interior** optimum in which both halves genuinely matter (a prototype
   confirmed the mechanism is correct; the plant must be made less
   actuation-dominated for it to un-floor the sensors).
2. **Better/limited actuation and richer sensing.** A single scalar brace
   authority corners immediately; per-story independent authorities, a
   rate/stroke-limited actuator (so good estimation is *necessary* to use scarce
   authority), and multiple binding drifts across floors would produce a genuine
   multi-axis interior allocation.
3. **Documented benchmark structure.** Replace the normalized shear building with
   a standard seismic-control benchmark (e.g. the ASCE/Spencer 3- or 20-storey
   building) so masses, stiffnesses, damping, and ground-motion records are
   traceable — the prerequisite for a structural-engineering venue.
4. **Recorded ground motions and time-domain validation.** Drive the optimized
   design with real earthquake records (not just stationary white excitation) and
   report peak drift, floor acceleration, and actuator force/stroke usage — the
   metrics a structural reviewer expects.
5. **$\ell_1$ sparse joint placement.** Use the sensor-allocation $\ell_1$
   relaxation (Sec. 7) to answer *where* to install a small set of dampers and
   sensors on a taller building in one convex sweep, rather than tuning a fixed
   instrumentation.

---

## 3. Summary and honest caveats

| Example | Field | ContEst path | $J_{\mathrm{tot}}$ | co-design gate | distinctive capability |
|---|---|---|---:|---:|---|
| Cybergenetics | Synthetic biology | H2/LQG-SDP ($\theta$ in $A$) | −55.6% | +31.7% | estimator (reporter) is a designed sensor; $\theta_h$ augments the filter state |
| Structural | Civil engineering | H2/LQG-SDP | −58.4% (drift −60.3%) | +22.4% | exact LQG-SDP gradient on a linear plant; hard physical damage metric |

Both are strong candidates, but two caveats apply to both as currently
formulated (and are the natural targets of the development sections above):

- **Both optima are corners** (one axis floored), whereas the paper emphasizes
  *interior* trades. Cybergenetics floors translation; structural floors two of
  three sensors. Re-scaling, richer sensing/actuation, or (for structural) the
  drift chance-constraint are the routes to interior optima.
- **One cost improves while the other worsens** — the signature of a *fixed*
  shared budget (estimation-dominated for cybergenetics, actuation-dominated for
  structural). This is expected and not a bug, but a reviewer will want at least
  one example where both halves improve simultaneously (achievable by loosening
  the budget, at the cost of a smaller co-design gate).

All gradients were verified against finite differences ($\sim10^{-5}$), confirming
the ContEst machinery operates correctly across the H2/LQG-SDP path with $\theta$
entering $B$, $C$, $V$, and $A$.
