# Literature Review — Control and Estimation Co-Design (ContEst)

**Prepared for:** M. S. Ramadan & M. Anitescu, *Control and Estimation Co-Design*
**Scope:** Prior work relevant to extending Control Co-Design (CCD) to jointly optimize plant, controller, **and** estimator under a stochastic (dual-control) objective, solved as a bilevel program with envelope-theorem sensitivities, and demonstrated on a power-systems example in Julia.

---

## How this review is organized

The draft sits at the intersection of several literatures that have historically developed separately. Each section below is one of those threads, ordered roughly to follow the logic of your paper: the co-design idea itself (§1–§5), the stochastic/dual-control machinery that makes *estimation* part of the objective (§6), the estimation-side tools that let sensor/actuator choices enter the cost (§7–§8), the mathematical sensitivity results that support the bilevel solve (§9), the estimation/covariance foundations behind your notation (§10), and the power-systems + Julia tooling for the simulation (§11).

Entries already cited in your current draft are marked **[in draft]**. The central novelty claim — that the *estimator* is co-designed alongside plant and controller, not just the controller — appears genuinely under-addressed in every thread below, which is a good sign for the contribution.

---

## 1. Control Co-Design (CCD): the framing you extend

**Garcia-Sanz, M. (2019). "Control Co-Design: An engineering game changer." *Advanced Control for Applications*, 1(1), e18. [in draft — Paper 1]**
The conceptual anchor for your abstract. Argues that treating control as the last, sequential stage of design is suboptimal, and that plant and controller should be conceived together from the start. Introduces the three-pillar view of CCD (control-inspired paradigms, formal co-optimization, co-simulation). Useful as the paper you position ContEst against: it never treats the *estimator/sensing* configuration as a first-class design variable.

**Allison, J. T. & Herber, D. R. (2014). "Multidisciplinary Design Optimization of Dynamic Engineering Systems." *AIAA Journal*, 52(4), 691–710.**
The formal backbone of modern CCD. Extends MDO, which was built around static models, to systems with dynamics, and lays out the taxonomy of solution strategies (sequential, iterative, nested, simultaneous). Your bilevel formulation is essentially the "nested" strategy with the estimator folded into the inner objective — this paper is the natural reference for that structure.

**Herber, D. R. & Allison, J. T. (2019). "Nested and Simultaneous Solution Strategies for General Combined Plant and Controller Design Problems." *ASME J. Mechanical Design*, 141(1), 011402.**
Removes many restrictions earlier co-design theory imposed and gives optimality conditions for both nested and simultaneous formulations. Directly relevant: your upper-level (θ) / lower-level (u) split is exactly the "general combined plant and controller" problem, and their optimality conditions are the deterministic counterpart to what you derive stochastically.

**Fathy, H. K., Reyer, J. A., Papalambros, P. Y. & Ulsoy, A. G. (2001/2003). "Nested Plant/Controller Optimization with Application to Combined Passive/Active Automotive Suspensions" and "On the Coupling between the Plant and Controller Optimization Problems."**
The papers that formalized *why* sequential design is suboptimal — they characterize the bidirectional coupling between plant and control problems and show when the sequential solution is provably wrong. Strong support for the "not optimal" claim in your abstract; you extend the coupling argument to a third block (the estimator).

**Sundarrajan, A. K. & Herber, D. R. (2021). "Towards a Fair Comparison between the Nested and Simultaneous Control Co-Design Methods using an Active Suspension Case Study." *ACC 2021*, 358–365.**
A careful head-to-head of the two solution strategies. Relevant because you will likely need to justify a nested/bilevel solve over a monolithic one, and this is the reference that discusses the trade-offs concretely.

**Allison, J. T., Guo, T. & Han, Z. (2014). "Co-Design of an Active Suspension Using Simultaneous Dynamic Optimization." *ASME J. Mechanical Design*, 136(8), 081003.**
Canonical simultaneous-CCD example via direct transcription. A useful template for how a co-design objective is discretized and solved, and a candidate comparison point.

**Deshmukh, A. P. & Allison, J. T. (2016). "Multidisciplinary Dynamic Optimization of Horizontal Axis Wind Turbine Design." *Structural and Multidisciplinary Optimization.***
Wind-turbine CCD showing measurable energy gains over sequential design. Relevant to your power-systems application and to motivating CCD in energy systems.

**Nash, A. L. & Jain, N. (2019). "Combined Plant and Control Co-Design for Robust Disturbance Rejection in Thermal-Fluid Systems." *IEEE Trans. Control Systems Technology*, 28(6), 2532–2539.**
CCD explicitly organized around disturbance rejection — closer in spirit to your stochastic objective than the deterministic-performance CCD papers, since disturbances are where estimation starts to matter.

**Silva, J. G. D. et al. / Gadelha et al. (2024). "Engineering and Technology Applications of Control Co-Design: A Survey." *IEEE Access*, DOI 10.1109/ACCESS.2024.3412416.**
Recent applications survey. Good for a one-paragraph "state of CCD" in your intro and for confirming, by omission, that estimator co-design is not standard practice.

**Herber, D. R. (2020). "What is Control Co-Design? Overview of CCD Research, Open Questions and Impact."**
A talk/whitepaper cataloguing open problems in CCD (e.g., appropriate models for CCD vs. control-only). Handy for framing ContEst as addressing one of the named gaps.

---

## 2. Integrated Structure-and-Control Design (the classical co-design lineage)

**Shi, G. & Skelton, R. E. (1996). "An Algorithm for Integrated Structure and Control Design with Variance Bounds." *Proc. 35th IEEE CDC*, Kobe, 167–172. [in draft — Paper 2]**
Your Paper 2. Alternating LMI/Lyapunov algorithm that jointly tunes structural parameters (stiffness, damping, actuator location) and a state-feedback gain to minimize control effort under output-variance (RMS) constraints. The variance-bound formulation is a direct ancestor of your `tr(QΣ)` / `tr(QΣε)` cost decomposition — worth citing when you introduce the covariance-based costs.

**Hale, A. L., Lisowski, R. J. & Dahl, W. E. (1985). "Optimal Simultaneous Structural and Control Design of Maneuvering Flexible Spacecraft." *J. Guidance, Control, and Dynamics*, 8(1), 86–93.**
One of the earliest simultaneous structure/control optimizations. Establishes the historical pedigree of co-design well before the term existed; useful in a "prior art" paragraph.

**Lim, K. B. & Junkins, J. L. (1989). "Robustness Optimization of Structural and Controller Parameters." *J. Guidance, Control, and Dynamics*, 12(1), 89–93.**
Jointly optimizes structure and controller for robustness. Complements Shi–Skelton and reinforces that the plant/control coupling was recognized early in aerospace.

**Grigoriadis, K. M., Zhu, G. & Skelton, R. E. (1996). "Optimal Redesign of Linear Systems." *ASME J. Dynamic Systems, Measurement, and Control.***
The LMI/redesign machinery underlying Paper 2. Relevant if you want to connect your convex (LMI) inner-problem formulation to this lineage.

**Skelton, R. E. (1989). "Model Error Concepts in Control Design." *Int. J. Control*, 49(5), 1725–1753.**
Argues modeling and control are not separable. A conceptual precursor to co-design and, indirectly, to the estimation-matters argument, since model error is what estimation must handle.

---

## 3. Integrated Design and Control in Process Systems Engineering (the CMU / PSE thread — Biegler, Grossmann)

**Grossmann, I. E. & Morari, M. (1984). "Operability, Resiliency, and Flexibility: Process Design Objectives for a Changing World."**
The foundational statement that dynamic operability/controllability is a *property of the design itself*, not something control can fully fix afterward. This is the process-systems version of your thesis and a strong citation for the "designed-in" argument.

**Biegler, L. T. (2010). *Nonlinear Programming: Concepts, Algorithms, and Applications to Chemical Processes.* SIAM.** and **Kameswaran, S. & Biegler, L. T. (2006). "Simultaneous Dynamic Optimization Strategies: Recent Advances and Challenges." *Computers & Chemical Engineering*, 30, 1560–1575.**
Biegler's simultaneous (full-discretization / direct-transcription) approach to dynamic optimization is the numerical engine that made large integrated design-and-control problems solvable. Directly relevant to how you will actually solve the ContEst bilevel/NLP, and to citing the CMU dynamic-optimization tradition you flagged.

**Flores-Tlacuahuac, A. & Grossmann, I. E. (2007–2011). Integrated design, control, and scheduling via Mixed-Integer Dynamic Optimization (MIDO); e.g., polymer grade-transition work.**
Shows how design + control (and later scheduling) are posed as a single MIDO and why solving them separately is suboptimal. The MIDO framing is a close cousin of your bilevel program; useful if any of your design variables θ are discrete (e.g., sensor *placement*).

**Pistikopoulos, E. N. & Diangelakis, N. A. (2015); Diangelakis, Burnak, Katz & Pistikopoulos (2017). "Process Design and Control Optimization: A Simultaneous Approach by Multi-Parametric Programming." *AIChE Journal*, 63(11), 4827–4846.**
Multiparametric programming derives an explicit closed-form controller that is then embedded in the design optimization — an elegant way to collapse the bilevel structure. Worth contrasting with your envelope-theorem approach to handling the inner problem.

**Sakizlis, V., Perkins, J. D. & Pistikopoulos, E. N. (2004). "Recent Advances in Optimization-Based Simultaneous Process and Control Design." *Computers & Chemical Engineering.***
A standard review that categorizes integrated design/control formulations. Good for a compact survey citation.

**Ricardez-Sandoval, L. A., Budman, H. M. & Douglas, P. L. (2009). "Integration of Design and Control: A Review and New Perspectives." and Yuan, Z. et al. (2012) review.**
Two widely cited reviews of the integrated design-and-control field. Between them they cover the controllability-measure vs. dynamic-optimization schools. Use to situate ContEst and to point out that estimator design is essentially absent from both schools.

---

## 4. Design for Controllability / Input–Output Controllability

**Skogestad, S. & Postlethwaite, I. (2005). *Multivariable Feedback Control: Analysis and Design*, 2nd ed. Wiley.**
Contains the crisp statement that input–output controllability — the best achievable performance — is fixed by sensor/actuator locations and the plant, and cannot be recovered by the controller ("even the best control system cannot make a Ferrari out of a Volkswagen"). This is almost a slogan for your paper and a natural framing citation, especially because it explicitly ties controllability to *sensor and actuator placement* — your θ.

**Skogestad, S. (1996). "A Procedure for SISO Controllability Analysis — with Application to Design of pH Neutralization Processes." *Computers & Chemical Engineering*, 20, 373–386.**
The quantitative "design for controllability" procedure. Relevant background for the controllability-analysis school and a concrete example of design choices limiting achievable control.

**Morari, M. (1983). "Design of Resilient Processing Plants III: A General Framework for the Assessment of Dynamic Resilience."**
Early formalization of assessing achievable dynamic performance at the design stage. Historical anchor for design-for-controllability.

---

## 5. Multidisciplinary Design Optimization (MDO) — the NASA / Alexandrov thread

**Alexandrov, N. M. & Lewis, R. M. (2000/2004). "Analytical and Computational Aspects of Collaborative Optimization for Multidisciplinary Design" (*AIAA Journal*, 38(2)) and "Reconfigurability in MDO Problem Synthesis, Parts 1–2" (AIAA 2004-4307/4308).**
Alexandrov's core MDO contributions: how you *formulate* a coupled multidisciplinary problem (collaborative optimization, disciplinary autonomy) strongly affects tractability. Directly relevant to how you decompose ContEst across the plant/control/estimation "disciplines," and to the bilevel-vs-monolithic choice.

**Alexandrov, N. M., Dennis, J. E., Lewis, R. M. & Torczon, V. (1998). "A Trust-Region Framework for Managing the Use of Approximation Models in Optimization." *Structural Optimization*, 15, 16–23** (and the related Approximation/Model Management Framework papers).
Provenly convergent management of variable-fidelity models inside optimization. Highly relevant if the inner stochastic-control cost is expensive to evaluate and you use surrogate/approximate `Jstoch` evaluations in the θ-loop.

**Alexandrov, N. M. & Hussaini, M. Y., eds. (1997). *Multidisciplinary Design Optimization: State of the Art.* SIAM.**
The reference MDO volume. Good single citation for the MDO field and for the "disciplinary autonomy" perspective you may lean on.

**Sobieszczanski-Sobieski, J. & Haftka, R. T. (1997). "Multidisciplinary Aerospace Design Optimization: Survey of Recent Developments." *Structural Optimization*, 14, 1–23.** and **Analytical Target Cascading (Kim, Michelena, Papalambros & Jiang, 2003).**
Broader MDO survey plus ATC, a hierarchical decomposition method. Relevant if you frame ContEst as a hierarchical (upper-level design / lower-level operation) decomposition — ATC is the MDO analog of your bilevel split.

---

## 6. Dual Control & Stochastic Optimal Control (the conceptual core — why estimation is *in* the cost)

**Feldbaum, A. A. (1960–61). "Dual Control Theory, I–IV." *Avtomatika i Telemekhanika.* [in draft]**
The origin of the dual effect: control simultaneously regulates and probes/learns the state. This is exactly the mechanism by which your objective depends on the filtering density, so it belongs at the head of your stochastic-control section (as it already does).

**Bar-Shalom, Y. & Tse, E. (1974). "Dual Effect, Certainty Equivalence, and Separation in Stochastic Control." *IEEE TAC*, 19(5), 494–500. [in draft]** and **Tse, E. & Bar-Shalom, Y. (1973). "An Actively Adaptive Control for Linear Systems with Random Parameters via the Dual Control Approach."**
Defines precisely when the dual effect vanishes (certainty equivalence / separation) and when it does not. Central: your claim that estimator design affects closed-loop cost is precisely a statement that the dual effect is present, so this is the reference that makes your point rigorous.

**Kumar, P. R. & Varaiya, P. (1986/2015). *Stochastic Systems: Estimation, Identification, and Adaptive Control.* SIAM. [in draft]**
The information-state / Bayesian-filter formulation you invoke in Lemma-based derivations. Also the standard reference for why the optimal policy is a function of the conditional density, not the state.

**Bertsekas, D. P. *Dynamic Programming and Optimal Control* (Vol. I–II).**
Standard reference for stochastic DP with imperfect state information and the information-state reformulation. Useful to cite alongside Kumar–Varaiya for the Bellman/DP backbone.

**Mesbah, A. (2018). "Stochastic Model Predictive Control with Active Uncertainty Learning: A Survey on Dual Control." *Annual Reviews in Control*, 45, 107–117.**
The modern survey. Cleanly separates *explicit* dual control (exploration bonus / persistent excitation) from *implicit* dual control (approximating the information-state Bellman equation). This is the single best citation to connect ContEst to contemporary work and to explain how "estimation quality" enters a control objective in practice.

**Filatov, N. M. & Unbehauen, H. (2000). "Survey of Adaptive Dual Control Methods." *IEE Proc. Control Theory & Applications*, 147(1).**
The pre-MPC survey; complements Mesbah with the classical adaptive-dual-control literature.

**Heirung, T. A. N., Ydstie, B. E. & Foss, B. (2017). "Dual Adaptive Model Predictive Control." *Automatica*, 80, 340–348.**
A concrete, implementable dual MPC. Useful as a comparison point for how anticipated future covariance shapes the control action — the same covariance objects that appear in your `Jest(θ)`.

**Klenske, E. D. & Hennig, P. (2016). "Dual Control for Approximate Bayesian Reinforcement Learning." *JMLR.*** and **Bayard, D. S. & Eslami, M. (1985). "Implicit Dual Control..."**
Two representative approximate-dual-control approaches (Bayesian and implicit). Relevant to positioning your envelope-theorem treatment as an alternative to scenario/ADP approximations of the information-state problem.

---

## 7. Sensor / Actuator Placement, Observability, and Estimator Accuracy

**Summers, T. H., Cortesi, F. L. & Lygeros, J. (2016). "On Submodularity and Controllability in Complex Dynamical Networks." *IEEE Trans. Control of Network Systems*, 3(1), 91–101.**
Shows several controllability/observability-Gramian metrics are submodular, so greedy selection has constant-factor guarantees. Relevant if any component of θ is a *combinatorial* sensor/actuator placement, and for grounding "information gathering quality" in Gramian metrics.

**Tzoumas, V., Jadbabaie, A. & Pappas, G. J. (2016). "Sensor Placement for Optimal Kalman Filtering: Fundamental Limits, Submodularity, and Algorithms." *ACC 2016*.**
The most direct precedent for the estimation half of ContEst: places sensors to bound the *Kalman filter error covariance*, proving submodularity of a log-det criterion. Your `Jest(θ) = tr(QΣε)` is a close relative (A-optimal flavor); cite this to show the estimator-design problem is well-posed and to contrast the co-design (joint with control/plant) vs. estimation-only viewpoint.

**Joshi, S. & Boyd, S. (2009). "Sensor Selection via Convex Optimization." *IEEE Trans. Signal Processing*, 57(2), 451–462.**
Convex-relaxation approach to choosing k of m sensors to minimize estimation-error volume, tied explicitly to D-optimal experiment design. The cleanest bridge between "sensor configuration as a design variable" and convex optimization — very useful for making θ-selection tractable.

**Krener, A. J. & Ide, K. (2009). "Measures of Unobservability." *CDC 2009*.**
Local-observability-Gramian measures for nonlinear systems. Relevant because your general (nonlinear) formulation `h(x,v;θ)` needs a nonlinear notion of how θ affects observability.

**Qi, J., Sun, K. & Kang, W. (2015). "Optimal PMU Placement for Power System Dynamic State Estimation by Using Empirical Observability Gramian." *IEEE Trans. Power Systems*, 30(4), 2041–2054.**
Directly connects §7 to your power-systems application: chooses PMU locations to maximize empirical observability for dynamic state estimation, validated with an unscented KF on WSCC/NPCC systems. A natural template for the estimation-design variable in your Julia example.

**Zhang, H., Ayoub, R. & Sundaram, S. (2017). "Sensor Selection for Kalman Filtering of Linear Dynamical Systems: Complexity, Limitations and Greedy Algorithms." *Automatica*, 78, 202–210.**
Establishes NP-hardness and characterizes when greedy works for the a priori/a posteriori error-covariance objectives. Good for honestly stating the computational difficulty of the estimator-design sub-problem.

---

## 8. Optimal Experimental Design & Information Metrics (formalizing "information-gathering quality")

**Pukelsheim, F. (2006). *Optimal Design of Experiments.* SIAM (Classics).**
The reference text on A-, D-, E-optimality as convex criteria on the Fisher information matrix. Since your estimation cost is a trace of an error covariance (an A-optimal-type criterion), this is the principled home for the "information quality" language in your abstract.

**Uciński, D. (2005). *Optimal Measurement Methods for Distributed-Parameter System Identification.* CRC Press** (and the D-optimal monitoring-network work).
Extends optimal experimental design to sensor *placement* for dynamical/distributed systems, using convex functions of the FIM. Bridges classical OED and the spatial sensor-placement problem you face.

*(See also Joshi & Boyd, §7, which explicitly links sensor selection to D-optimal design — the trace/log-det/min-eigenvalue tri{A,D,E}-optimality maps onto the covariance-based costs `tr(QΣ)`, `tr(QΣε)` in your linear/LQG section.)*

---

## 9. Parametric Sensitivity & the Envelope Theorem (machinery for the bilevel upper-level gradient)

**Milgrom, P. & Segal, I. (2002). "Envelope Theorems for Arbitrary Choice Sets." *Econometrica*, 70(2), 583–601.**
The modern, general envelope theorem: the derivative of a value function w.r.t. a parameter equals the partial derivative of the objective at the optimizer, without differentiating the optimizer. This is precisely the result that lets you differentiate `Jstoch(p0;θ)` w.r.t. θ at the upper level while treating the optimal `u` as fixed — the mathematical justification for your bilevel gradient. Cite it where you invoke the envelope theorem.

**Bonnans, J. F. & Shapiro, A. (2000). *Perturbation Analysis of Optimization Problems.* Springer.**
The authoritative treatment of value-function differentiability, directional derivatives, and sensitivity under constraints. The rigorous backstop for the envelope argument when your lower-level problem is constrained (`u∈U`, `x∈X`) and the minimizer may be non-unique.

**Fiacco, A. V. (1983). *Introduction to Sensitivity and Stability Analysis in Nonlinear Programming.* Academic Press.**
Classical NLP sensitivity (how optimal value and solution move with parameters). The go-to citation for differentiating a parametric NLP value function — useful if you compute θ-gradients through KKT conditions rather than the pure envelope form.

**Danskin, J. M. (1967). *The Theory of Max-Min and Its Applications.* Springer.**
Danskin's theorem — the min/max special case of the envelope theorem. Often the cleanest citation when the inner problem is a minimization and you want the derivative of the min.

**Dempe, S., Mordukhovich, B. S. & Zemkoho, A. B. (2012). "Sensitivity Analysis for Two-Level Value Functions with Applications to Bilevel Programming." *SIAM J. Optimization*, 22(4), 1309–1343.** and **Colson, B., Marcotte, P. & Savard, G. (2007). "An Overview of Bilevel Optimization." *Annals of OR*, 153, 235–256.**
Sensitivity of bilevel value functions and a broad survey of bilevel methods (KKT reformulation, MPEC). These place your θ/u problem in the formal bilevel-programming literature and provide the tools (and caveats) for computing sensitivities when the envelope conditions are delicate.

---

## 10. Estimation Foundations & Covariance-Based Control (behind your notation, Paper 3)

**Anderson, B. D. O. & Moore, J. B. (1979/2012). *Optimal Filtering.* [in draft]**
Your Kalman-filter reference and the source of the steady-state covariance/Riccati machinery your `Σ`, `Σε`, `Pε` equations rest on.

**Kalman, R. E. (1960); Kalman, R. E. & Bucy, R. S. (1961).**
The primary filtering references. Worth an explicit citation when you first write the observer/error dynamics (eqs. for `x̂`, `ε`).

**Ramadan, M. S. (2024). [L-CSS paper — Paper 3]. *IEEE Control Systems Letters*, Vol. 8.**
Your own notation reference for the covariance/moment reformulation of the LQ cost (`E{xᵀQx} = tr(QΣ) + xᵀQx` form), EKF, and control-estimation coupling. It is the notational and technical bridge to ContEst; cite it for the moment-based cost rewriting and for the varying-observability/SNR numerical setup you can reuse.

**Okamoto, K., Goldshtein, M. & Tsiotras, P. (2018). "Optimal Covariance Control for Stochastic Systems Under Chance Constraints." *IEEE L-CSS*, 2(2), 266–271** (and Bakolas; Rapakoulias & Tsiotras, 2023, SDP covariance steering).
Covariance steering treats the state *covariance* as the controlled object under chance constraints — the same covariance quantities that ContEst designs through θ. A modern, convex/SDP-friendly viewpoint that could inform your inner-problem formulation and chance-constraint handling.

---

## 11. Power-Systems Application & Julia Tooling (for the simulation)

**Lara, J. D., Henriquez-Auba, R., Bossart, M., Callaway, D. S. & Barrows, C. (2023). "PowerSimulationsDynamics.jl — An Open-Source Modeling Package for Modern Power Systems with Inverter-Based Resources." arXiv:2308.02921 / NREL Sienna.**
The most likely simulation backbone for your example: a Julia package for time-domain power-system dynamics with rich synchronous-machine and inverter libraries, quasi-static-phasor and EMT-dq models, and (importantly) compatibility with Julia's ML/AD ecosystem for parameter fitting and optimization. This is where θ (actuator/sensor/controller tuning) can be swept and `Jstoch` evaluated.

**Henriquez-Auba, R. et al. (2020). "LITS.jl — An Open-Source Julia-Based Simulation Toolbox for Low-Inertia Power Systems." arXiv:2003.02957.**
The predecessor toolbox; its modular inverter meta-model (filter, converter, inner/outer control, PLL frequency estimator) is a clean place to expose *estimator* design variables (e.g., PLL tuning) — a nice, concrete instance of co-designing an estimator.

**PowerSystems.jl (NREL Sienna).**
Data layer used by the simulators; relevant for building the test network and managing time series. Cite for reproducibility of the example.

**JuliaGrid (2025). "An Open-Source Julia-Based Framework for Power System State Estimation." arXiv:2502.18229.**
Julia framework covering state estimation, observability analysis, and optimal PMU placement at scale. Directly supports the *estimation* side of your simulation and the sensor-placement θ.

**Zhao, J. et al. (2019). "Power System Dynamic State Estimation: Motivations, Definitions, Methodologies, and Future Work." *IEEE Trans. Power Systems*, 34(4), 3188–3198.**
The authoritative DSE survey. Good for motivating why dynamic state estimation (and thus estimator design) matters in modern grids with high inverter penetration — the application-side justification for ContEst.

**Ghahremani, E. & Kamwa, I. (2016). "Local and Wide-Area PMU-Based Decentralized Dynamic State Estimation in Multi-Machine Power Systems." *IEEE Trans. Power Systems*, 31(1), 547–562.**
Representative EKF/UKF-based DSE with real PMU data. Useful as a realistic estimator model to co-design against.

---

## Suggested "closest prior work" to address head-on

If a reviewer asks "what is genuinely new," these are the works nearest your contribution, and the sentence-level distinction you can draw:

- **Tzoumas–Jadbabaie–Pappas / Joshi–Boyd** design the *estimator's* sensing but hold the plant and controller fixed.
- **Herber–Allison / Fathy / Shi–Skelton** co-design *plant and controller* but treat sensing/estimation as given (full-state feedback or a fixed observer).
- **Mesbah / Heirung dual control** couple control and *learning/estimation* online, but for a *fixed* physical design — they don't optimize θ (the plant/sensor hardware).

ContEst appears to be the first to close all three loops at once — plant, controller, **and** estimator — inside a single stochastic (dual-control) objective, using the envelope theorem to make the bilevel solve tractable. That triangulation is your defensible novelty.

---

## Notes on gaps / next searches worth running

1. **Economic/robust CCD under uncertainty** (e.g., robust MDSDO, chance-constrained CCD) — a few results exist and would strengthen the "stochastic CCD" positioning.
2. **POMDP / belief-space planning** (Platt, Kaelbling, Todorov) — the robotics analog of information-state control; relevant if you want to connect to that community.
3. **Simultaneous perturbation / differentiable optimization layers** (e.g., differentiating through Riccati/LQR solutions) — an implementation alternative to the envelope theorem for computing θ-gradients.
4. Confirm exact bibliographic details (year/volume/pages) for the classics (Feldbaum parts, Kalman 1960, Bertsekas edition, Milgrom–Segal) before the final bib.
