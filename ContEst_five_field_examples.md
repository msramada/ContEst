# Five Cross-Disciplinary Case Studies for ContEst

Five examples from five distinct scientific fields, each a genuine Control-and-Estimation Co-Design problem: a design variable that enters **actuation** ($\theta_f$, through $\bar f$) *and* a design variable that enters **sensing/estimation** ($\theta_h$, through $\bar h$, the measurement covariance $V(\theta)$, or the filter state itself), coupled so tightly that the sequential *place-then-control* pipeline sits at a non-stationary point of the joint problem.

They are chosen so that, collectively, they exercise **every** inner regime in the paper:

| # | Field | Inner regime showcased | Sensing axis $\theta_h$ acts through | Discrete allocation |
|---|-------|------------------------|--------------------------------------|---------------------|
| 1 | Neuroscience (closed-loop DBS) | eKF–MPC (Sec. 6) | $\bar h$ + co-located electrodes | ℓ1 contact selection (Sec. 7) |
| 2 | Synthetic biology (cybergenetics) | eKF–MPC (Sec. 6) | augmented filter state $+\,V(\theta)$ | ℓ1 over a part/reporter library |
| 3 | Structural engineering (seismic control) | **LQG SDP** (Sec. 5) | $C_\theta$ + $V(\theta)$ | ℓ1 joint sensor/actuator placement |
| 4 | Fluid dynamics (active flow control) | **H∞ / game-Riccati** (Sec. 5.2) | $C_\theta$ + $V(\theta)$ | ℓ1 sparse sensing/actuation |
| 5 | Epidemiology (surveillance + NPI) | eKF–MPC (Sec. 6) | $V(\theta)=v_0/\theta_h$ (rate budget) | shared testing budget |

Throughout we follow the paper's conventions: continuous-time dynamics discretized by an explicit Euler step at sample time $\Delta t$ to the form $x_{k+1}=\bar f(x_k,u_k;\theta)+w_k$, $y_k=\bar h(x_k;\theta)+v_k$, with $w_k\sim\mathcal N(0,\Sigma_w)$, $v_k\sim\mathcal N(0,\Sigma_v(\theta))$; $Q\succeq 0$, $R\succ 0$; design split $\theta=(\theta_f,\theta_h)$; total cost $J_{tot}=J_{des}+J_c$, $J_c=J_{det}+J_{est}$, $J_{est}=\sum_t\mathrm{tr}(Q\Sigma_{t|t})$.

---

## 1. Neuroscience — Closed-loop deep brain stimulation (DBS)

**Field & significance.** Deep brain stimulation for Parkinson's disease and essential tremor suppresses pathological beta-band ($13$–$30$ Hz) synchrony in the basal ganglia–thalamocortical loop. Adaptive (closed-loop) DBS is the frontier: it modulates stimulation in response to a recorded biomarker rather than delivering fixed high-frequency pulses. The defining feature — and what makes this a ContEst problem rather than a pure control problem — is that **the DBS lead's contacts are simultaneously the sensors (local field potential, LFP) and the actuators (stimulation current)**. A contact that best *records* the pathological oscillation need not be the one that best *disrupts* it, and only a subset of contacts can be active under a charge-safety budget. Sequential design (choose a recording contact by anatomy, choose a stim contact by trial-and-error, then tune a fixed-gain controller) cannot reconcile this trade.

**Dynamics.** Model $n$ neural populations, each a damped oscillator (LFP proxy $v_i$ and its rate $\dot v_i$) coupled synaptically; the pathology is a large-amplitude beta limit cycle. Let $x_i=[v_i,\dot v_i]^\top$, $x=[x_1;\dots;x_n]\in\mathbb R^{2n}$:

$$
\ddot v_i + 2\gamma_i\dot v_i + \Omega_i^2 v_i \;=\; \sum_{j\ne i}\kappa_{ij}\,S(v_j)\;+\;\theta_{f,i}\,b_i\,u_i\;+\;\xi_i,
\qquad S(v)=\tanh(\eta v),
$$

with natural frequencies $\Omega_i$ in the beta band, damping $\gamma_i$, synaptic coupling $\kappa_{ij}$, stimulation gain $b_i$ scaled by the **design effectiveness $\theta_{f,i}$** of contact $i$, and process noise $\xi_i$. The saturating coupling $S(\cdot)$ is the genuine nonlinearity (it generates the limit cycle); in deviation form about the desynchronized equilibrium $\bar f(0,0;\theta)=0$. Discretize at $\Delta t$ to $x_{k+1}=\bar f(x_k,u_k;\theta)+w_k$, with the eKF Jacobian $F_k=\partial\bar f/\partial x$ carrying the $\tanh$ slopes.

**Output equation.** Each contact $i$ can record its local LFP with a design gain/precision $\theta_{h,i}\ge 0$ ($\theta_{h,i}=0$ ⇒ contact not used for sensing):

$$
y_k \;=\; \operatorname{diag}(\theta_h)\,C_{\mathrm{LFP}}\,x_k \;+\; v_k,
\qquad C_{\mathrm{LFP}} = I_n\otimes[\,1\;\;0\,],
\qquad v_k\sim\mathcal N\!\big(0,\;V(\theta_h)\big),
$$

with $V(\theta_h)=\operatorname{diag}(v_0/(\varepsilon+\theta_{h,i}^2))$: a higher recording gain lowers readout noise with saturating benefit (amplifier/impedance limit $\varepsilon$).

**Estimation objective.** Reconstruct the *distributed* oscillation state (phase and amplitude across all $n$ populations) from the few active contacts — far richer than the single-contact band-power biomarker used in deployed adaptive DBS:
$$
J_{est}=\sum_t\mathrm{tr}\!\big(Q\,\Sigma_{t|t}\big).
$$

**Control objective.** Suppress oscillation power with minimal, charge-balanced stimulation:
$$
J_{det}=\sum_{t}\Big(\|x_t\|_Q^2+\|u_t\|_R^2\Big),\qquad
Q=\operatorname{diag}(q_v,0)\otimes I_n\ (\text{penalize }v_i),\quad
|u_{i,t}|\le u_{\max}\ (\text{safety}).
$$

**Design variables.** $\theta_f=(\theta_{f,1},\dots,\theta_{f,n})\in[0,1]^n$ (per-contact stimulation effectiveness, in $\bar f$) and $\theta_h=(\theta_{h,1},\dots,\theta_{h,n})\in[0,\theta_{\max}]^n$ (per-contact recording precision, in $\bar h$).

**Design objective.** Tissue-safety/energy penalty plus a hardware budget forcing a sparse *active-contact* set (a lead has $\le 8$ contacts, and simultaneous sense+stim on the same contact incurs artifact):
$$
J_{des}(\theta)=c_e\textstyle\sum_i\theta_{f,i}^2 \;+\; \lambda_s\big(\|\theta_f\|_1+\|\theta_h\|_1\big),
$$
the ℓ1 term (Sec. 7) selecting which contacts stimulate and which record.

**How it beats sequential design.** A control-only co-design has no slot for $\theta_h$ and inherits the clinician's default recording contact; a pure sensor-placement method optimizes $\theta_h$ for observability alone and ignores that the informative contact competes with stimulation for the same hardware. ContEst resolves the sense-vs-stim contact allocation jointly, and — because the eKF exposes the distributed state — enables a covariance-aware controller that outperforms band-power triggering.

---

## 2. Synthetic biology — Optogenetic control of gene expression (cybergenetics)

**Field & significance.** Cybergenetics closes a feedback loop around a living cell: light actuates a light-inducible promoter, a fluorescent reporter is measured, and a controller regulates protein expression in real time. The subtlety that makes this ContEst — and that the cybergenetics literature routinely under-models — is that **the reporter is a dynamic, noisy sensor with its own maturation lag**: newly translated fluorophore is dark and matures with rate $k_{\mathrm{mat}}$, so fluorescence is a *delayed, filtered* readout of the true protein, not a direct measurement. The reporter choice (brightness, maturation rate, which species it is fused to) is a design variable that sets both $V(\theta)$ *and the sensor dynamics*. Designing the circuit's promoter/RBS strengths (actuation) together with the reporter (sensing) is a genuine co-design.

**Dynamics.** Target mRNA $m$, protein $p$, and a two-stage reporter (immature $r_1$, mature/fluorescent $r_2$). Light input $u$ enters through a Hill induction; growth dilution and degradation are linear:

$$
\begin{aligned}
\dot m &= \alpha_0 + \theta_{f,\alpha}\,\alpha\,\frac{u^{n}}{K^{n}+u^{n}} - \delta_m m,\\
\dot p &= \theta_{f,\beta}\,\beta\,m - \delta_p p,\\
\dot r_1 &= \theta_{h,\beta_r}\,\beta_r\,m - (\delta_r+\theta_{h,\mathrm{mat}}\,k_{\mathrm{mat}})\,r_1,\\
\dot r_2 &= \theta_{h,\mathrm{mat}}\,k_{\mathrm{mat}}\,r_1 - \delta_r r_2,
\end{aligned}
$$

state $x=[m,p,r_1,r_2]^\top$ (extended per additional target gene). The Hill term is the nonlinearity (linearize about the set-point operating concentration for the eKF; the deviation form preserves positivity locally). Note the **sensing design $\theta_h$ augments the state** through the reporter subsystem — a filtering feature none of the paper's current studies exhibit.

**Output equation.** Fluorescence reads the *mature* reporter only, with brightness $\phi$ and shot/photobleaching noise whose variance falls with brightness:
$$
y_k=\theta_{h,\phi}\,\phi\,r_{2,k}+v_k,\qquad v_k\sim\mathcal N\!\big(0,\,v_0/\theta_{h,\phi}\big).
$$

**Estimation objective.** Infer the true, unmeasured protein $p$ from the delayed noisy fluorescence — a nontrivial deconvolution of the maturation lag:
$$
J_{est}=\sum_t \mathrm{tr}\!\big(Q\,\Sigma_{t|t}\big),\qquad Q=\operatorname{diag}(0,q_p,0,0).
$$

**Control objective.** Track a reference protein level $p^\star$ with bounded light (phototoxicity):
$$
J_{det}=\sum_t\Big(q_p\,(p_t-p^\star)^2+R\,u_t^2\Big),\qquad 0\le u_t\le u_{\max}.
$$

**Design variables.** $\theta_f=(\theta_{f,\alpha},\theta_{f,\beta})$ — promoter and ribosome-binding-site (translation) strengths, chosen from characterized part libraries, in $\bar f$; $\theta_h=(\theta_{h,\beta_r},\theta_{h,\mathrm{mat}},\theta_{h,\phi})$ — reporter fusion strength, maturation rate, and brightness, in the sensor dynamics and $V(\theta)$. The maturation choice is a **bias–variance knob** (fast-dim vs. slow-bright), the exact analog of the PLL-bandwidth study but in a biological sensor.

**Design objective.** Metabolic burden of strong expression, reporter phototoxicity, and discreteness of the part library:
$$
J_{des}(\theta)=c_{\mathrm{burden}}\big(\theta_{f,\alpha}+\theta_{f,\beta}\big) + c_{\mathrm{tox}}\,\theta_{h,\phi} + \lambda_s\|\theta_h\|_1\ (\text{part selection}).
$$

**How it beats sequential design.** Standard practice fixes a convenient reporter (e.g., a slow-maturing GFP), builds the circuit, then designs a controller *assuming fluorescence equals protein*. ContEst accounts for the reporter's lag and noise in the achievable closed-loop cost and co-selects a reporter whose dynamics the controller can actually invert — a coupling invisible to any place-then-control workflow. It is also the cleanest demonstration that *estimator hardware* (the reporter) is a first-class design variable.

---

## 3. Structural engineering — Seismic protection of a building (LQG-SDP showcase)

**Field & significance.** Active and semi-active control of civil structures under earthquake and wind loading is a mature field (Skelton, Spencer, Dyke), and it is *the* natural home for the paper's LQG-SDP machinery: the plant is linear, the cost is quadratic (drift and acceleration), and both actuator (damper) placement and sensor (accelerometer) placement are classically posed as LMIs. The estimation axis is essential because a tall building is instrumented on only a few floors, yet the controller needs the full drift/velocity state — and the ground motion itself is unmeasured. This example is the ideal vehicle for the revision's **exact LQG regularity result**: under stabilizability/detectability and $Q,R,W,V\succ0$, the sensor-and-actuator-placement SDPs (Eqs. (9)–(10)) have unique duals equal to the Riccati/Kalman solutions, so the design gradient is exact.

**Dynamics.** An $n$-story shear building with mass, damping, stiffness matrices $M,\,C_d,\,K$, inter-story drift vector $q\in\mathbb R^n$, ground acceleration $\ddot x_g$ (disturbance), and control forces $u$:
$$
M\ddot q + C_d\dot q + K q \;=\; -M\,\iota\,\ddot x_g \;+\; \Gamma(\theta_f)\,u,
$$
where $\iota$ is the influence vector and $\Gamma(\theta_f)=\Gamma_0\operatorname{diag}(\theta_f)$ places dampers with per-floor authority $\theta_{f,i}\ge0$. In first-order form $x=[q;\dot q]\in\mathbb R^{2n}$,
$$
\dot x=\underbrace{\begin{bmatrix}0&I\\-M^{-1}K&-M^{-1}C_d\end{bmatrix}}_{A_c}x
+\underbrace{\begin{bmatrix}0\\M^{-1}\Gamma(\theta_f)\end{bmatrix}}_{B_c(\theta_f)}u
+\underbrace{\begin{bmatrix}0\\-\iota\end{bmatrix}}_{E}\ddot x_g,
$$
discretized to $x_{k+1}=A_\theta x_k+B_\theta u_k+w_k$ with $W=\Sigma_w$ the (whitened) ground-excitation covariance. **Linear ⇒ the inner problem is exactly the LQG SDP pair (9)–(10).**

**Output equation.** Accelerometers/interstory-drift sensors on a chosen subset of floors, with placement gates and precision in $\theta_h$:
$$
y_k=C(\theta_h)x_k+v_k,\qquad C(\theta_h)=\operatorname{diag}(\theta_h)\,C_{\mathrm{base}},\qquad v_k\sim\mathcal N(0,V(\theta_h)),
$$
$V(\theta_h)=\operatorname{diag}(v_0/\theta_{h,j}^2)$; $\theta_{h,j}=0$ removes sensor $j$. (Augment $x$ with a ground-motion shaping filter to estimate $\ddot x_g$; its states are unmeasured and enter $J_{est}$.)

**Estimation objective.** $J_{est}=\mathrm{tr}(Q\Sigma^\varepsilon)$ — reconstruct all drifts, velocities, and the ground excitation from sparse floor sensors.

**Control objective.** Minimize the stationary drift/acceleration variance (damage + occupant comfort) and control force:
$$
J_{det}=\mathrm{tr}(Q\Sigma)+\mathrm{tr}(RK\Sigma K^\top),\qquad
Q=\operatorname{blkdiag}(Q_{\mathrm{drift}},Q_{\mathrm{vel}}),\quad |u_j|\le u_{\max}.
$$

**Design variables.** $\theta_f\in[0,\bar\theta]^{n}$ (damper placement/authority, in $B_\theta$), $\theta_h\in[0,\bar\theta]^{m}$ (accelerometer placement/precision, in $C_\theta,V$).

**Design objective.** Hardware cost of dampers and sensors, with ℓ1 promoting a small installed set of each:
$$
J_{des}(\theta)=c_a\|\theta_f\|_1+c_s\|\theta_h\|_1 \quad(\text{joint sparse actuator + sensor placement}).
$$

**How it beats sequential design.** The textbook workflow is three separate steps: place sensors to maximize an observability Gramian, place actuators to maximize a controllability Gramian, then design LQG on the frozen instrumentation. ContEst shows the observability-optimal sensor set is *not* the closed-loop-optimal set once actuator locations and the spatial disturbance profile are accounted for jointly — the sensors should watch the floors the available actuators can actually protect against the actual ground-motion mode shape. Because the plant is linear, the improvement is certified by the exact LQG gradient with unique duals, and the ℓ1 relaxation answers the discrete placement question in one convex sweep.

---

## 4. Fluid dynamics — Active flow control / drag reduction (H∞ showcase)

**Field & significance.** Feedback control of transitional and turbulent flows (skin-friction drag reduction, separation control, cavity-tone suppression) is a flagship control-theory application. The design decisions are *where to place wall sensors* (pressure taps, shear sensors) and *where to place actuators* (synthetic jets, plasma actuators) on a body, to control a flow whose state (the vorticity/modal-amplitude field) is only sparsely observed. Crucially, the relevant disturbance — an upstream gust or free-stream turbulence — is an **adversarial, bounded-energy** signal, not white noise: what a secure design must bound is the *worst-case* amplification from disturbance to perturbation energy. This makes flow control the ideal showcase for the revision's **robust ($H_\infty$) inner value and its game-Riccati route (Sec. 5.2)**.

**Dynamics.** Reduced-order model of the flow by POD–Galerkin or balanced truncation about a base flow, modal amplitudes $a\in\mathbb R^{r}$ ($r$ modest, e.g. $10$–$50$):
$$
\dot a = A\,a + B(\theta_f)\,u + E\,w,\qquad z=\big[\,Q^{1/2}a\,;\,R^{1/2}u\,\big],
$$
where $A$ is the linearized (mean-flow) operator, $B(\theta_f)=B_0\operatorname{diag}(\theta_f)$ places actuators with authority $\theta_f$, $E$ is the gust-input map, $w\in\ell_2$ is the **bounded-energy** disturbance, and $z$ is the performance output stacking perturbation kinetic energy and actuation effort. (The neglected quadratic Galerkin term $a^\top\!\mathcal Q\,a$ is the nonlinearity; linearizing about the base flow yields the linear ROM used for $H_\infty$ synthesis.) Discretize to $x_{k+1}=A_\theta x_k+B_\theta u_k+E_\theta w_k$.

**Output equation.** Sparse wall sensors with placement/precision $\theta_h$:
$$
y_k=C(\theta_h)a_k+D_v v_k,\qquad C(\theta_h)=\operatorname{diag}(\theta_h)C_{\mathrm{base}},\qquad D_v=V(\theta_h)^{1/2}.
$$

**Estimation objective (worst-case).** Reconstruct the modal state from wall sensors with a guaranteed error-energy bound — the $H_\infty$ filtering value $\gamma^2_{\mathrm{est}}(\theta)$ from SDP (17), or its game-Riccati twin. This is the flow-estimation problem, notoriously ill-conditioned under sparse sensing.

**Control objective (worst-case).** Minimize the squared worst-case $\ell_2$ gain from gust $w$ to performance output $z$:
$$
\gamma^2_{\mathrm{cont}}(\theta)=\min_{K}\ \sup_{w\in\ell_2,\,w\ne0}\ \frac{\|z\|_{\ell_2}^2}{\|w\|_{\ell_2}^2},
$$
obtained from the bounded-real SDP (16) or the fixed-$\gamma$ game Riccati (20)–(22). Actuator authority is bounded by an energy budget.

**Design variables.** $\theta_f$ (actuator locations/authority, in $B_\theta$) and $\theta_h$ (sensor locations/precision, in $C_\theta,V$).

**Design objective.** Actuator power and manufacturability, with ℓ1 promoting a sparse, realizable sensor/actuator layout:
$$
J_{des}(\theta)=c_p\|\theta_f\|_1+c_s\|\theta_h\|_1.
$$

**How it beats sequential design.** Flow-control practice places sensors by an *open-loop* estimation criterion (maximize the ROM observability Gramian) and actuators by a controllability criterion, then designs the compensator — a fully sequential pipeline. Against a worst-case gust this is doubly fragile: the open-loop-optimal sensor set need not observe the modes the closed loop most needs, and the white-noise ($H_2$) criterion under-weights the resonant amplification a gust actually excites. ContEst co-places sensors and actuators against the *closed-loop worst-case gain*, and the game-Riccati route delivers the exact design gradient at $O(r^3)$ without the rank-deficient LMI dual that plagues minimal-$\gamma$ synthesis — precisely the numerical advantage argued in Sec. 5.2/9.6.

---

## 5. Epidemiology — Co-design of surveillance and intervention (sensing = testing)

**Field & significance.** Epidemic response allocates two scarce resources across regions: **intervention capacity** (isolation, treatment, targeted vaccination — the actuation) and **surveillance capacity** (testing, genomic sampling — the sensing). The pandemic made vivid that under-surveilled regions harbor hidden outbreaks: testing rate directly sets how noisily the true prevalence is observed. Public-health practice allocates testing by population or past incidence (a pure monitoring criterion) and interventions separately — a textbook sequential design. ContEst poses the question correctly: surveillance should be concentrated where it most improves *control outcomes*, which depends on mobility coupling and where intervention capacity exists, not on current case counts alone. The sensing axis maps to $V(\theta)=v_0/\theta_h$ exactly as the WAMS PMU-rate study, but here the latent compartment (exposed-but-not-yet-infectious) makes estimation genuinely hard.

**Dynamics.** Networked SEIR over $n$ regions with mobility coupling $\{m_{ij}\}$; state per region $[S_i,E_i,I_i]$ (with $R_i=N_i-S_i-E_i-I_i$):
$$
\begin{aligned}
\dot S_i &= -\beta_i\frac{S_iI_i}{N_i} - \sum_{j}m_{ij}\big(S_i/N_i-S_j/N_j\big)N_i,\\
\dot E_i &= \beta_i\frac{S_iI_i}{N_i} - \sigma E_i - \sum_j m_{ij}(\cdots),\\
\dot I_i &= \sigma E_i - \gamma I_i - \theta_{f,i}\,u_i\,I_i - \sum_j m_{ij}(\cdots),
\end{aligned}
$$
where $u_i\ge0$ is the intervention intensity in region $i$ (isolation/treatment that shortens the infectious period), scaled by the **provisioned capacity $\theta_{f,i}$**. The bilinear incidence $\beta_i S_iI_i/N_i$ and the control term $\theta_{f,i}u_iI_i$ are the nonlinearities; linearize about the current operating point for the eKF. Discretize to $x_{k+1}=\bar f(x_k,u_k;\theta)+w_k$, $x=[S;E;I]\in\mathbb R^{3n}$.

**Output equation.** Region $i$ reports noisy case counts whose observation gain and variance are set by its **testing rate $\theta_{h,i}$** (a fixed total testing budget is allocated across regions):
$$
y_{i,k}=I_{i,k}+v_{i,k},\qquad v_{i,k}\sim\mathcal N\!\big(0,\;v_0/\theta_{h,i}\big),\qquad \textstyle\sum_i\theta_{h,i}\le B_{\mathrm{test}}.
$$
Only $I$ is (partially) reported; $E$ and $S$ are unobserved and must be inferred.

**Estimation objective.** Infer true prevalence and the latent exposed pool from under-reported, noisy counts:
$$
J_{est}=\sum_t\mathrm{tr}\!\big(Q\,\Sigma_{t|t}\big),\qquad Q\ \text{weighting } I \text{ (and } E\text{) heavily}.
$$

**Control objective.** Minimize cumulative/peak infections at bounded socioeconomic cost of intervention:
$$
J_{det}=\sum_t\Big(\|I_t\|_Q^2+\|u_t\|_R^2\Big),\qquad 0\le u_{i,t}\le u_{\max},\ \textstyle\sum_i \theta_{f,i}\le B_{\mathrm{cap}}.
$$

**Design variables.** $\theta_f\in[0,\bar\theta]^n$ (intervention-capacity allocation, in $\bar f$) and $\theta_h\in[0,\bar\theta]^n$ (testing-rate allocation, in $V(\theta)$), each under a resource budget.

**Design objective.** Provisioning cost of capacity and testing, with the budgets entering as soft penalties:
$$
J_{des}(\theta)=c_b\Big(\textstyle\sum_i\theta_{f,i}-B_{\mathrm{cap}}\Big)^2+c_t\Big(\textstyle\sum_i\theta_{h,i}-B_{\mathrm{test}}\Big)^2.
$$

**How it beats sequential design.** The naive baseline tests in proportion to population and intervenes in proportion to reported incidence — two decoupled monitoring/response rules. ContEst reallocates surveillance toward regions where a mobility-driven hidden outbreak would most degrade the achievable control, *given* where intervention capacity sits, and provisions capacity in step. Because the informative regions for estimation are set by the network coupling (not local incidence), and the effective regions for control are set by capacity, the joint optimum is unreachable by any place-then-control rule — the epidemiological instance of the paper's central thesis, and a timely, high-impact demonstration.

---

## Summary: what each example proves about ContEst

| Example | Distinctive framework capability demonstrated | Sequential baseline it defeats |
|---------|-----------------------------------------------|-------------------------------|
| **DBS** | Co-located sensor/actuator hardware forces a genuine sense-vs-actuate allocation; ℓ1 contact selection on a nonlinear oscillator network | Anatomy-chosen recording + trial-and-error stim + fixed-gain trigger |
| **Cybergenetics** | Sensing design *augments the filter state* (reporter maturation lag) and shapes $V(\theta)$; bias–variance in a biological sensor | Convenient reporter assumed to equal the true state |
| **Structural** | Exact LQG-SDP path with unique duals (the revision's regularity result); joint sparse sensor+actuator LMI placement | Gramian-based sensor placement → separate actuator placement → LQG |
| **Flow control** | Robust $H_\infty$ inner value + game-Riccati gradient (Sec. 5.2), the right tool for adversarial gusts | Open-loop-Gramian sensor/actuator placement + $H_2$ compensator |
| **Epidemiology** | $V(\theta)=v_0/\theta_h$ rate-budget sensing with a genuinely latent state; network-coupled informativeness | Test ∝ population, intervene ∝ incidence (decoupled) |

Collectively the five span both inner regimes (LQG and $H_\infty$ linear; eKF-MPC nonlinear), all three sensing channels ($\bar h$, $V(\theta)$, and the filter state itself), and both continuous tuning and ℓ1 discrete allocation — so adopting any subset strengthens the paper's claim that ContEst is a *general* framework, not a power-systems method. For a first strong addition I would build out **Structural** (it certifies the new exact-gradient theory on a linear plant) and **Epidemiology** (highest novelty/impact and a clean surveillance-as-sensing narrative), keeping the other three as a "further applications" discussion or a companion paper.
