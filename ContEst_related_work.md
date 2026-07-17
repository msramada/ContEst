# Related Work Map: Riccati/Envelope Design Gradients and Co-Design

**Purpose.** Audit of prior art bearing on ContEst's claim that the design gradient of an LQ/H∞ inner value is obtained exactly and cheaply from inner duals. Organized by threat to the current framing.

**Verification key.**
- ✅ = retrieved and read (abstract or full text) during this pass
- 📄 = seen only in the reference list of a retrieved paper; bibliographic details need confirming before citing
- 📁 = already in the project folder

**Threat key.** 🔴 = must cite and reframe against · 🟠 = must cite · 🟡 = cite for completeness

---

## Family A — The gradient formula is classical: "parametric LQ"

> **Bottom line:** the two-gramian formula (one Riccati for P, one Lyapunov for Σ, contract) is Levine–Athans 1970. ContEst's §5 Riccati route rederives it. This family alone forces a reframing of what §5 contributes.

### A1. Levine & Athans (1970) — 🔴 ✅
*"On the determination of the optimal constant output feedback gains for linear multivariable systems," IEEE TAC 15(1):44–48.*
**Method.** Static output-feedback LQ posed as parameter optimization over the gain. The cost is tr(P·Σ₀)-type; its gradient requires solving two coupled Lyapunov equations (cost gramian P, state-covariance gramian Σ), then contracting. The dK terms cancel by stationarity.
**vs. ContEst.** This *is* the ContEst §3 derivation, with the design parameter being K instead of θ. The cancellation ContEst calls "the envelope theorem in its rawest form" is the Levine–Athans cancellation. Structural identity — cite in §5, not §2.

### A2. Mäkilä & Toivonen (1987) — 🔴 ✅
*"Computational methods for parametric LQ problems — A survey," IEEE TAC 32(8):658–671.*
**Method.** Survey of the whole 1970s–80s machinery for **LQ values parametrized by design data**, including gradient/Hessian formulas and quasi-Newton outer loops.
**vs. ContEst.** The title names ContEst's inner problem. This is the single most dangerous "we already knew this" citation. Companion: Toivonen & Mäkilä, *"On Newton's method for solving parametric linear quadratic control problems,"* Int. J. Control 46:897–911 (1987) — a Newton outer loop on a parametric LQ, i.e. ContEst §8 with BFGS→Newton. ✅

### A3. Toivonen (1985); Toivonen & Mäkilä (1985); Moerder & Calise (1985) — 🟡 ✅
*"A globally convergent algorithm for the optimal constant output feedback problem,"* Int. J. Control 41(6):1589–1599; *"A descent Anderson–Moore algorithm for optimal decentralized control,"* Automatica 21:743–744; *"Convergence of a numerical algorithm for calculating optimal output feedback gains,"* IEEE TAC 30(9):900–903.
**Method.** Convergence theory for descent on the Levine–Athans gradient.
**vs. ContEst.** Prior art for §8's stationarity argument in the LQ case. Weaker threat (they descend on gains, not design), but they establish that gradient descent on a Riccati-derived LQ value is a solved genre.

### A4. Rautert & Sachs (1997) — 🟡 ✅
*"Computational design of optimal output feedback controllers," SIAM J. Optimization 7(3):837–852.*
**Method.** Structure-exploiting Newton/quasi-Newton on the output-feedback LQ cost using the same gramian gradients.
**vs. ContEst.** Closest classical antecedent to §8's "exact gradients ⇒ good secant pairs ⇒ superlinear BFGS" argument.

### A5. Fazel, Ge, Kakade & Mesbahi (2018) and the policy-gradient line — 🟠 ✅
*"Global convergence of policy gradient methods for the linear quadratic regulator," ICML.*
**Method.** Rediscovery of ∇J(K) = 2((R+B⊤PB)K + B⊤PA)Σ_K — Levine–Athans in modern dress — plus gradient dominance / global convergence analysis. Successors: Mohammadi–Zare–Soltanolkotabi–Jovanović; Bu–Mesbahi *"LQR through the lens of first order methods"* (arXiv:1907.08921) ✅; Polyak et al. *"Optimizing static linear feedback: gradient method,"* SIAM J. Control Optim. (2021) ✅.
**vs. ContEst.** Same formula, different variable (policy vs. design). Relevant because reviewers from this community will recognize ContEst's §3 immediately. Survey to cite: Hu, Zheng, Mesbahi, Fazel & Başar, *"Toward a theoretical foundation of policy optimization for learning control policies,"* Annual Reviews in Control (2023) ✅ — it cites Levine–Athans and Mäkilä–Toivonen as the ancestors, confirming the lineage is openly acknowledged in that community.

### A6. Observer-based dynamic LQ landscape (2026) — 🟠 ✅
*"On the Optimization Landscape of Observer-based Dynamic Linear Quadratic Control," arXiv:2604.10635.*
**Method.** Closed-form gradients of the OD-LQR cost w.r.t. **both** controller gain K and observer gain L, via two Lyapunov solutions (S_{K,L}, Ω_{K,L}); characterizes when the stationary point collapses to the standard separation pair.
**vs. ContEst.** Uncomfortably close to §5's Jcont/Jest split and Appendix A's separation-within-the-inner-solve. Differs in that the variables are *gains*, not physical design θ, and there is no sensing/plant axis. Worth reading in full before finalizing §5.

---

## Family B — Riccati sensitivity and automatic differentiation

### B1. Kao & Hennequin (2020) — 🔴 ✅
*"Automatic differentiation of Sylvester, Lyapunov, and algebraic Riccati equations," arXiv:2011.11430.* Code: `tachukao/autodiff-inverse-lqr`.
**Method.** Forward- and reverse-mode derivatives of the solutions of all three matrix equations. The reverse mode of the ARE reduces to a Lyapunov solve against the closed-loop matrix. Demonstrated on an inverse-control problem.
**vs. ContEst.** This is the general version of ContEst's Appendix-B(c) implicit-function step, packaged as an AD primitive. **Direct challenge to the "no implicit differentiation" pitch:** they make ARE differentiation a one-line library call. ContEst's genuine edge is narrower and must be stated as such — when the value is *exactly* tr(PW), the ARE adjoint has the closed form Σ and no Lyapunov solve with a general RHS is needed. That is a property of the objective, not a defect of AD.

### B2. Kenney & Hewer; Konstantinov et al.; Byers; Sun — 🟡 📄
*"The sensitivity of the algebraic and differential Riccati equations"* and successors (perturbation theory for AREs).
**Method.** Condition numbers and perturbation bounds for ARE solutions; the linearization is the closed-loop Stein/Lyapunov operator, invertible iff the closed loop is stable.
**vs. ContEst.** Exactly Proposition 2's argument, done in the 1980s–90s for the H₂ DARE. ContEst's contribution there is the *fixed-γ game* version. Cite to anchor Prop. 2; verify the exact refs.

### B3. Mania, Tu & Recht (2019) — 🟡 ✅
*"Certainty equivalence is efficient for linear quadratic control," arXiv:1902.07826.*
**Method.** Explicit upper bounds on the sensitivity of the Riccati solution to (A,B,Q,R) perturbations, used to bound CE controller suboptimality; extends to LQG.
**vs. ContEst.** Not co-design, but it is the quantitative theory of ∂P/∂θ. Relevant to Remark 9(iii)'s open question — "error bounds relating the surrogate stationary point to a stationary point of (3)" — this line is where such bounds would come from.

---

## Family C — Co-design with a Riccati inner problem

### C1. He, Kaneko, Howell, Li & Martins (Feb 2026) — 🔴 ✅ **← the collision**
*"Efficient Adjoint-based Design Optimization with Optimal Control," arXiv:2602.15242.*
**Method.** Nested closed-loop CCD: outer plant design **d**, inner LQR. Three coupled residuals (equilibrium solve, ARE, closed-loop time march) differentiated by a **coupled adjoint** exploiting feed-forward structure → three smaller adjoint solves by block back-substitution. The ARE adjoint follows Kao & Hennequin and reduces to a Lyapunov equation. Cost independent of the number of design variables. Verified vs. FD to 5–7 digits. Cases: cart-pole (76.8% cost cut), quadrotor blade with 20 BEM design variables (10% LQR cost cut for 3% hover-power penalty, Pareto front vs. sequential design).
**vs. ContEst.**
- **Same thesis, same target.** Sequential design is suboptimal → nested CCD with an analytic (non-FD) design gradient through the Riccati equation.
- **Kills Table 1's "gradient w.r.t. design: often FD" cell** for the CCD/MDO column. They explicitly position *against* Herber–Allison and Sundarrajan–Herber for using FD.
- **Names ContEst §5.2 and §6 as their future work:** "extends naturally to other closed-loop control formulations beyond LQR, such as MPC and H∞."
- **What they lack — ContEst's actual moat:** no sensing/output map, no estimator, no stochastic/information-state objective, no sensor allocation, no H∞, no eKF. Their objective is a deterministic trajectory cost from one initial condition, not tr(QΣ)+tr(QΣᵋ).
- **Action:** cite in §1 and §2; rewrite Table 1's gradient row; add them to the Positioning paragraph as line (ii).

### C2. Herber & Allison (2019) — 🟠 📁 (your [22])
*"Nested and simultaneous solution strategies for general combined plant and control design problems," J. Mech. Design 141(1):011402.*
**Method.** Taxonomy of nested vs. simultaneous CCD; LQR-inner examples; design derivatives by finite differences.
**vs. ContEst.** Already cited. Note that C1 supersedes its gradient story — do not let Table 1 imply that FD is the state of the art in CCD.

### C3. Fathy, Reyer, Papalambros & Ulsoy (2001) — 🟠 📁 (your [16])
**Method.** Coupling conditions between plant and controller optimization; nested formulation with LQR inner.
**vs. ContEst.** Already the anchor for "sequential is not stationary." Companion worth adding: Reyer, Fathy, Papalambros & Ulsoy (2001), *"Comparison of combined embodiment design and control optimization strategies using optimality conditions"* 📄.

### C4. Aerospace structure/control simultaneous optimization (1985–1990) — 🟠 📄
Hale, Lisowski & Dahl (1985) [your ref 21]; Eastep, Khot & Grandhi (1987); Onoda & Haftka (1987) [your ref 30]; Rao (1988) *"Combined structural and control optimization of flexible structures"*; Haftka (1990); Belvin & Park (1990); **D. F. Miller (1990), "Combined structural and control optimization: a steepest descents approach"** (Control and Dynamic Systems, vol. 32) — this last one is reference [4] of your own `2.pdf`.
**Method.** Nested structural design over an LQ/eigenvalue inner problem, with design sensitivities of the LQ cost obtained analytically or by forward differentiation.
**vs. ContEst.** This is the 1980s ancestor of C1 and of ContEst's §5. Miller's title alone ("steepest descents") indicates gradient-based structure/control co-design. **Get and read Miller (1990) before submission** — it is cited inside a paper already in your project folder, so "we didn't see it" won't fly.

### C5. Grigoriadis, Zhu & Skelton (1996); Shi & Skelton (1996) — 🟠 📁 (your [20]; `2.pdf`)
*"Optimal redesign of linear systems," JDSMC 118(3):598–605; "An algorithm for integrated structure and control design with variance bounds," CDC.*
**Method.** Alternate between a convex control LMI (fixed plant) and a convex plant-parameter LMI (fixed Lyapunov matrix). Convergence to a local optimum guaranteed; explicitly non-convex jointly.
**vs. ContEst.** The LMI-alternation alternative to gradient descent. ContEst's advantage over this is real and worth sharpening: alternation freezes P or θ at each step and (per the follow-on "Convexifying LMI methods" paper ✅) converges slowly with no rate guarantee, whereas ContEst descends θ directly with exact curvature. This is a better contrast to feature than the FD contrast, which C1 has taken away.

### C6. Cunis, Kolmanovsky & Cesnik (2023) — 🟡 📄
*"Integrating nonlinear controllability into a multidisciplinary design process," J. Guidance, Control, and Dynamics 46(6):1026–1037.*
**Method.** Analytic derivatives of open-loop optimal control for MDO; nonlinear controllability as a design criterion.
**vs. ContEst.** The "design for controllability" thread you asked me to track, on the aerospace side. Open-loop, no estimation.

### C7. Kaneko & Martins (2025); Sundarrajan & Herber (2021) — 🟡 📄
*"Simultaneous design and trajectory optimization strategies for computationally expensive models," AIAA J. 63(2):420–438; "Towards a fair comparison between the nested and simultaneous CCD methods using an active suspension case study," ACC.*
**vs. ContEst.** Background for §8's nested-vs-simultaneous choice.

---

## Family D — The estimation/sensing axis

### D1. Belabbas (2016) — 🔴 ✅
*"Geometric methods for optimal sensor design," Proc. Royal Society A* (arXiv:1503.05968).
**Method.** Minimize J = tr(LK) over the sensing matrix C, K the filter Riccati solution. Differentiating the ARE gives a Lyapunov equation in K̇; contracting yields ∇J(C). Cast on the Grassmannian with the "normal metric"; proves the objective is Morse with a unique minimum and the gradient flow converges globally for small SNR and stable A.
**vs. ContEst.** **The estimation-half Riccati gradient already exists.** Structurally this is ContEst's ∂J_est/∂C_θ = −2G⋆⊤P^ε A_f Σ^ε. Differences: (i) no control axis — the plant and controller are fixed, so no co-design; (ii) continuous-time, orthonormal C only; (iii) he gets *global* convergence, which ContEst does not claim. **Must cite.** Its existence is also good news: it means the estimation gradient is well-posed and someone has done the hard geometry.

### D2. Tzoumas, Carlone, Pappas & Jadbabaie (2021) — 🟠 📁 (your [38])
*"LQG control and sensing co-design," IEEE TAC 66(4):1468–1483.*
**vs. ContEst.** Already correctly identified as nearest single reference. Unaffected by this pass.

### D3. Unified MIO for sensor scheduling/selection (2023) — 🟠 ✅
*"A unified approach to optimally solving sensor scheduling and sensor selection problems in Kalman filtering," arXiv:2304.02692.*
**Method.** Mixed-integer optimization exploiting Kalman optimality; solves to **global optimality** for 35–50 states in seconds; explicitly states that its general form "captures ... Linear Quadratic Gaussian (LQG) control and sensing co-design."
**vs. ContEst.** Threatens §7's framing that exact sensor allocation is hopeless (NP-hard, 2^ry) and must be relaxed. At MTDC's ry = 25 this method may solve the allocation exactly. Either cite and scope §7's claim to the *joint* continuous-θ + discrete-b problem (which is the honest distinction — they fix the plant), or be prepared for a reviewer to ask why the ℓ₁ relaxation is needed.

### D4. Optimal PMU placement for DAE power systems (2025) — 🟡 ✅
*arXiv:2502.03338.*
**Method.** Exact MISDP reformulation minimizing tr of the steady-state DARE covariance for descriptor/DAE power models; assimilated sensing precision S(γ) = Σ γᵢCᵢ⊤Rᵢ⁻¹Cᵢ.
**vs. ContEst.** Directly comparable to §9.4's sparse voltage-sensor allocation, and in the same application domain. A natural head-to-head baseline for the MTDC study: they get a certificate, ContEst gets the joint droop+sensing coupling. Also relevant to the planned Julia power-systems work.

### D5. Sensor placement complexity line — 🟡 ✅
Zhang, Ayoub & Sundaram [your 41]; Tzoumas, Jadbabaie & Pappas [your 39]; resilient/graph-based placement (arXiv:2006.14036).
**vs. ContEst.** Already covered by §2; no change.

---

## Consolidated verdict

| Claim as currently written | Status after this pass |
|---|---|
| "Design gradient exact from inner duals, no differentiation through the optimizer" (LQG/Riccati) | **Not new.** Levine–Athans 1970; Mäkilä–Toivonen 1987. Reframe as recalled. |
| "No implicit differentiation / no KKT differentiation needed" | **Narrow it.** True *because* the inner value is tr(PW); Kao–Hennequin + He et al. differentiate the ARE routinely and cheaply. |
| Table 1: CCD/MDO gradient = "often FD" | **False as of Feb 2026.** He et al. (C1). Rewrite the cell. |
| Sensor/estimation design via Riccati gradients is absent | **False.** Belabbas 2016 (D1). |
| Exact sensor allocation is intractable, hence ℓ₁ | **Over-claimed.** D3 solves it globally at your scales for fixed plant. Scope to the joint problem. |
| **Plant + controller + estimator in one stochastic objective, gradients from inner duals, across LQG / H∞ / eKF-MPC** | **Survives.** This is the paper. |

**The defensible triangle:** Belabbas has estimation-Riccati gradients without control; He et al. has control-Riccati gradients without estimation; Tzoumas et al. has sensing+LQG co-design but by menu selection through separation. ContEst is the intersection — and the envelope theorem is the *unifying lens*, not the novel mechanism.

## Immediate to-do

1. Obtain **Miller (1990)** and **Mäkilä–Toivonen (1987)** — the two most likely "this is old" reviewer strikes.
2. Read **He et al. (arXiv:2602.15242)** in full; add to §1/§2; rewrite Table 1's gradient row.
3. Read **Belabbas (2016)** in full; cite in §5 and §2 (sensor/estimator design line).
4. Rewrite the §4 opening: the envelope theorem is the *unification*, not the discovery.
5. Re-scope §7's NP-hardness framing against D3.
6. Check the process-systems side (Grossmann / Biegler / Pistikopoulos "design for controllability") — not covered in this pass; my expectation is DAE-NLP + back-off rather than Riccati gradients, i.e. no collision, but unverified.
