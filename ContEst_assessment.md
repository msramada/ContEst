# ContEst — Automatica Readiness Assessment

*Control and Estimation Co-Design (Ramadan & Anitescu). Reviewer's-eye assessment of the current draft (`main.tex`) and the intended contribution.*

## Bottom line

The **idea** has a real shot at Automatica. The **manuscript** is not close yet, and even the idea, in its current framing, is more likely to land as a CDC/ACC paper than an Automatica paper unless it is sharpened into a hard theoretical result plus a serious case study.

These are two separate judgments and should not be conflated.

## The current draft is a skeleton, not a submission

What exists today would be desk-rejected anywhere, simply because it is unfinished:

- The Introduction is an empty header.
- The Conclusion says "Later."
- There is a placeholder subsection "A more serios example" [sic] with no content.
- Working notes remain in the body ("seems like the cutest thing," "LPV-MPC discussion stuff," the to-do list around line 393).
- There are no numerical results.
- It is typeset in `ieeeconf` (a conference class), not Automatica format — suggesting the paper was mentally aimed at ACC/CDC first.

So the literal answer to "is it good enough" is *no*, but that is trivially true of any ~30%-drafted paper and is not the useful answer.

## Is the contribution Automatica-worthy? Conditionally, and not yet

### What works in its favor

- The unifying thread is clean and, as far as I know, not stated this cleanly before: moving from control co-design (plant + controller) to co-design of **plant + controller + estimator**, via the observation that the stochastic-OC cost depends on the *filtering density* and therefore on the output/sensing design.
- Lemmas 1–3 (smoothing theorem → cost depends on the conditional density → hence on the output equation) make that thread rigorous rather than hand-wavy.
- The LQG special case, where the LQR cost splits cleanly into `J_est(θ) = tr(Q Σᵉ)` and `J_cont(θ)`, each posed as an SDP/LMI, is concrete and tractable. Good backbone.

### What a good reviewer will press on

**1. Novelty must be provable, not architectural.**
Automatica rarely accepts "framework" or "useful way to think about it" papers. Much of the current content assembles known pieces — dual control and the information state (Feldbaum, Bar-Shalom), the separation structure, LMI/Gramian reformulations of Lyapunov problems. The two literatures being bridged already exist:
- Integrated plant-and-control design — Grossmann and Biegler's design-for-controllability and simultaneous design-under-uncertainty work.
- Sensor/estimator design — Joshi–Boyd sensor selection, Summers/Tzoumas submodular observability, Zare–Jovanović covariance/actuator–sensor design.

The reviewer's question will be: *what is the non-obvious theorem here, beyond "combine co-design with sensor design"?* One crisp result that is not immediate from the setup is needed.

**2. The envelope theorem is the most promising route to that theorem — and is currently the least developed part.**
The natural contributive result: use the envelope theorem to differentiate the lower-level value function `J_stoch(θ)` with respect to the design parameters `θ` *without* solving the full bilevel sensitivity — with rigorous conditions for validity (existence/uniqueness of the lower-level optimizer, differentiability, SDP non-degeneracy / strict complementarity so the dual `S` in the gradient section is unique). The current gradient sketch (`∂J_cont/∂A_θ = −2 S₁₂ Σ⋆`) is this idea in embryo, but with mismatched notation and no stated hypotheses. Turning this into a theorem with conditions and an algorithm would be an Automatica-scale methodological contribution. Not doing so leaves a survey-flavored framework paper.

**3. Rigor gaps reviewers will catch.**
- The separation principle is invoked to minimize `J_est` and `J_cont` "separately." That holds *for fixed θ*, but θ couples `A_θ, B_θ, C_θ`, the noise covariances, and hence both `K` and `G` — so at the design level the subproblems are **not** separable. The paper must state that separation is used only within the inner solve; as written it reads like a claim of design-level separability.
- The `J_est` derivation (≈ lines 259–262) has a stray `Σᵉ` inside the trace on line 260 that should not be there.
- The rewritten lower-level Problem (≈ line 388) is malformed: a doubled `+ +`, and the objective/constraints look truncated/unbalanced.
- `\label{eq:stateSpace}` is defined twice (lines 133 and 316).
- Chance constraints (`x_k ∈ X`, `u_k ∈ U` "with some probability") are asserted but never handled.
- In the eKF/nonlinear case, separation fails and `J_stoch(θ)` is generally nonconvex and possibly nonsmooth — the envelope-theorem gradients need far more care there than in the LQG case. The nonlinear section should not write checks the theory cannot cash.

**4. No validation.**
Automatica expects either deep theory + an illustrative example, or a compelling application study. The planned power-systems simulation is not optional — it is load-bearing for acceptance.

## Recommendations

**Scope down.** The draft currently tries to hold general nonlinear stochastic OC + eKF + LQG + SDPs + gradients + examples at once; Automatica papers are tightly scoped. The strongest paper hiding in here is:

> The linear/LQG control-and-estimation co-design problem, with an envelope-theorem-based gradient result (stated with conditions) that makes the bilevel problem efficiently solvable, validated on a real power-systems case study.

Leave the general nonlinear/eKF treatment as clearly-flagged future work or a follow-on.

**Suggested publication trajectory.** Publish the framework + LQG result + a compact example as a CDC or ACC paper first (the `ieeeconf` formatting suggests this was already the lean), then extend to the full Automatica version with the rigorous envelope-theorem theorem, the conditions, and a substantial power-systems study. This is the standard way this kind of contribution matures and it de-risks the Automatica submission.

## Suggested next steps

- **(a) Positioning:** draft the Introduction and a sharp "contributions + relation to prior work" paragraph that stakes the novelty claim explicitly against Grossmann/Biegler and the sensor-selection literature.
- **(b) Core theory:** develop the envelope-theorem gradient result properly — state assumptions, the theorem, and identify where it breaks in the eKF case. *This is the step that decides whether this is an Automatica paper.*
