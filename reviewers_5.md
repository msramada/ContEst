# Simulated Referee Panel — *Control and Estimation Co-Design (ContEst)*

> **What this is.** A simulated peer-review panel for a top control journal
> (*Automatica* / *IEEE TAC*). Five referees with deliberately different
> backgrounds and priorities read the current manuscript (`ContEst_TeX/main.tex`)
> and its supplements. Each review is written in the referee's voice; the reviews
> intentionally disagree with one another where real panels would. A
> cross-cutting synthesis and a suggested revision-priority list close the file.
>
> Recommendations use the usual scale: *Accept*, *Minor revision*,
> *Major revision*, *Reject & resubmit*, *Reject*.

## Editor's snapshot

| # | Referee background & primary lens | Recommendation |
|---|---|---|
| R1 | Convex/variational analysis; sensitivity & bilevel theory | **Major revision** |
| R2 | Control co-design / MDO / mechatronics; novelty & baselines | **Reject & resubmit** |
| R3 | Estimation, sensor selection, optimal experiment design | **Major revision** |
| R4 | Power-systems & applications; modeling fidelity | **Major revision** |
| R5 | Numerical optimization / scientific computing; scalability & rigor | **Minor revision** |

**One-line consensus.** The envelope-gradient idea is elegant and the theory
(LQG/H∞ regularity) is solid, but the *empirical case for co-design itself* is not
yet made: the headline reductions are measured against naive baselines, the one
competent-sequential comparison is essentially a null result, the nonlinear
studies optimize an unquantified surrogate, and the paper has **no figures** and
**no statistical error bars**.

---

## Reviewer 1 — Convex-analysis / sensitivity theorist

**Background & priorities.** I work on perturbation analysis of conic programs and
bilevel optimization. I care about whether the central theorem is correct, whether
its hypotheses are checkable, and whether the "exactness" claims survive contact
with the nonlinear regime and the nonsmooth penalty.

**Summary assessment.** The paper's mathematical spine — reading the design
gradient of an inner value from the returned dual via the envelope theorem — is
correct and, frankly, the nicest thing in the paper. Promoting Assumption 1 to
verifiable sufficient conditions (Prop. on LQG regularity; the fixed-γ H∞
analog) is exactly the right move and the appendices are sound. But the paper
repeatedly blurs *exact gradient of the surrogate* with *exact gradient of the
problem*, and its convergence/optimality language outruns what is proved.

**Strengths.**
- Theorem 1 + the identification of the LMI dual with the Riccati solution
  (`S₁₁★ = P`) is clean and genuinely useful; the appendix proofs are correct.
- The "exactness ledger" table is the right instinct for separating the three
  cost objects.

**Major concerns.**
1. **The nonlinear headline rests on a surrogate whose gap is never bounded.** In
   the eKF regime you optimize \(\widehat J_c\) — eKF moment closure + linearized
   information dynamics + certainty-equivalent policy. Theorem 1 then gives the
   *exact gradient of \(\widehat J_c\)*, not of \(J_\mathrm{stoch}\). Three of the
   four studies live here. You are candid about this (Remark on where exactness
   stops), but candor is not a bound. At minimum I want either (a) an a-posteriori
   certificate — e.g. an estimate of \(\|\nabla \widehat J_c - \nabla
   J_\mathrm{stoch}\|\) via a higher-fidelity rollout — or (b) a theorem relating a
   stationary point of \(\widehat J_c\) to an approximate-stationary point of
   \(J_\mathrm{stoch}\). Without one, "exact design gradients" in the abstract is
   misleading for the nonlinear case and should be qualified.
2. **Assumption 1(A3) is still assumed for the nonlinear MPC.** Your regularity
   propositions cover the LQG and fixed-γ H∞ SDPs. For the convex information-state
   MPC (Eq. for \(\widehat J_c\)) you again fall back to "verified a posteriori."
   State conditions on the lifted \((A^\mathrm{info}_\theta,B^\mathrm{info}_\theta)\)
   under which strict complementarity / unique dual holds, or explicitly restrict
   the exactness claim.
3. **Outer-loop optimality claims exceed the proof.** "Every accumulation point is
   stationary" is fine *for the smooth surrogate*; but you then call 5-start minima
   "the global optimum over Θ to numerical tolerance" (ADCS, PLL). For a nonconvex
   program, five starts certify nothing global — drop "global" or justify with a
   coverage/estimation argument.
4. **Epigraph ⇔ proximal equivalence is asserted, not shown.** The claim that the
   ℓ1 epigraph lift and the ISTA-type proximal step "either removes the kink and
   keeps the convergence theory" deserves one precise statement: for the
   *bound-constrained* problem with `α ≥ 0`, note that the lift is exact and the
   penalty is linear (so BFGS applies verbatim); for the sign-free case the
   equivalence is different. As written it reads as folklore.
5. **The best theorem is hidden in the supplement.** The dual-invariance result
   (returned-dual gradient equals \(\nabla_\theta V\) even without strict
   complementarity) is the strongest statement you have and it makes Remark
   (degenerate case) actionable. Put at least its statement in the body.

**Minor.** Define "nondegenerate" where Assumption 1 is stated, not only by
citation. The directional-derivative formula and the subgradient fallback should
share notation with the ledger.

**Recommendation: Major revision.** The core is right; the exactness rhetoric must
be brought into line with what is actually proved for each regime.

---

## Reviewer 2 — Control co-design / MDO

**Background & priorities.** Plant–controller co-design, nested vs. simultaneous
architectures, and the empirical question that co-design papers live or die on:
*does the joint design actually beat a competent sequential/alternating pipeline,
or just a strawman?*

**Summary assessment.** The paper's thesis is stated forcefully — sequential
design sits at a non-stationary point of the joint problem, so co-design wins. I
went looking for the evidence and did not find it. The one study that compares
against a *competent* sequential (and alternating) baseline reports a **0.3%**
improvement and concedes the problem is "in practice close to separable." Every
other study compares only against a naive/uniform baseline. So the marquee
6–43% reductions largely measure *the distance from a deliberately poor starting
design*, not the value of co-design over good engineering practice. Until that is
fixed, the paper does not support its own headline.

**Strengths.**
- Bringing sensing/estimation into the CCD objective is a genuinely worthwhile
  scope extension, and the positioning table is a good start.
- The MTDC droop result and its new sparse-sensing extension are the most
  convincing pieces.

**Major concerns.**
1. **The central claim is under-tested.** The thesis is "the sequential solution
   need not be stationary." Demonstrate it: for *every* study report (i) naive
   baseline, (ii) estimator-first-then-control sequential, (iii) alternating, and
   (iv) joint. Right now only PLL has (ii)–(iii), and there co-design ≈ sequential.
   I need at least one study where joint co-design *decisively* beats a competent
   sequential/alternating pipeline; otherwise the honest conclusion is "co-design
   rarely helps beyond a good sequential pass," which is a different (and still
   publishable, but much weaker) paper.
2. **Baselines are strawmen.** "Naive over-provisioned droop \(k=3\)", "uniform
   budget", "naive default placement" — no competent engineer would ship these. A
   reduction against them is not evidence of co-design value. Re-baseline against
   best sequential practice.
3. **Where is the coupling actually strong?** Your own LQG analysis shows the axes
   couple through \(\theta\); but coupling *strength* is what determines whether
   joint beats sequential. Add a diagnostic (e.g. the off-diagonal cross-sensitivity
   of Fathy et al.) per study, predicting when co-design pays — that would turn the
   negative PLL result into a feature ("ContEst tells you when you don't need it").
4. **MDO positioning is prose, not placement.** You say ContEst "sits within" MDO.
   Which architecture does it *instantiate*? It looks like a nested/MDF scheme whose
   discipline solve is convex and whose coupling derivative is supplied by duality;
   say that, and contrast with collaborative optimization's consistency constraints.
5. **Speculative deployment claims.** The PLL "slow, event-triggered design layer"
   paragraph is interesting but unsupported (no triggered loop, no switching
   stability). Either demonstrate it or cut it to one sentence of future work.

**Minor.** The abstract's "6–43%" should be qualified by baseline type. "Dual
control" in the framing invites the estimation community's objections (see R3).

**Recommendation: Reject & resubmit.** The idea and machinery are good, but the
empirical argument for the paper's own thesis has to be rebuilt around competent
baselines before the contribution can be judged.

---

## Reviewer 3 — Estimation / sensor selection / OED

**Background & priorities.** Kalman filtering, sensor selection, optimal
experiment design. I judge the estimation half: the filter approximation, the
sensing-design formulation, and the comparison to the sensor-selection literature
the paper positions against.

**Summary assessment.** The estimator-as-design-variable idea (tuning a PLL
bandwidth that enters both \(\bar f\) and \(V(\theta)\)) is the freshest
contribution and I'd like to see it foregrounded. But the paper leans on "dual
control" branding it doesn't deliver, the eKF approximation is used without any
fidelity check, and the sensor-allocation section — now finally exercised on MTDC
— has no optimality guarantee and no head-to-head against the methods it cites.

**Strengths.**
- Co-designing the *estimator itself*, not just sensor placement, is novel and
  well-motivated by the PLL example.
- The A-optimal reading \(\mathrm{tr}(Q\Sigma^\varepsilon)\) correctly connects to
  OED, and the dual filter GARE for the robust case is elegant.

**Major concerns.**
1. **"Dual control" is oversold.** The title and abstract invoke dual control, but
   the method is (your words) caution-aware and certainty-equivalent-in-mean — it
   prices uncertainty but does not probe. That's *cautious* control, not dual
   control. Either retitle/reframe around "estimation-aware" co-design, or add one
   example with genuine probing. As written this will draw fire from anyone who
   works on the dual effect.
2. **eKF moment-closure gap is never checked.** Every nonlinear result depends on
   the eKF being an adequate information state. Provide a sanity check on at least
   one study against a UKF or a particle filter (bias in \(\Sigma_{t|t}\) vs. true
   conditional covariance). Otherwise the estimation costs \(J_\mathrm{est}\) you
   report and optimize may be systematically off.
3. **Sensor allocation lacks guarantees and a head-to-head.** The ℓ1 relaxation +
   threshold + re-solve is standard, but you provide no recovery/optimality bound
   and no comparison. On the MTDC sparse-sensing study, run exhaustive selection
   (feasible for small active sets) and greedy/submodular selection
   (Tzoumas 2016; Zhang 2017) and report the optimality gap of your ℓ1 solution.
   The multimodality you note makes this comparison essential.
4. **The Tzoumas 2021 "head-to-head" is not one.** You compare against "the
   separation-style sequential pipeline that such menu methods induce," not against
   their actual co-selection algorithm. Implement the real baseline on a common
   linear/LQG instance.
5. **Criterion choice unjustified.** Why A-optimal (\(\mathrm{tr}\)) rather than
   D-optimal (\(\log\det\)) or a task-relevant criterion? The choice interacts with
   the sensor-selection guarantees (log-det is submodular; trace is not always).

**Minor.** State the prior \(\Sigma_0\) sensitivity. Clarify whether \(V(\theta)\)
scaling in ADCS/PLL keeps the channel noise physically consistent as gains change.

**Recommendation: Major revision.** Reframe the dual-control claim, validate the
eKF, and give the sensor allocation a real comparison; the estimator-design angle
can carry the paper if these are addressed.

---

## Reviewer 4 — Power systems & applications

**Background & priorities.** I evaluate whether the case studies are faithful
enough that the conclusions would survive on a realistic system, and whether the
reported wins are physics or artifacts of tuned constants.

**Summary assessment.** The breadth (spacecraft, distillation, grid) is
attractive, but each model is a low-order caricature run over very short horizons,
and several results depend on constants chosen precisely to manufacture the
interior optimum the framework wants to show. I'm not yet convinced a practitioner
would see these gains on a real asset, and — decisively for an applications
reviewer — there is not a single time-domain plot.

**Strengths.**
- The MTDC study asks a real question (droop vs. converter-lag resonance) and the
  H∞ worst-case framing is the right tool for it.
- The PLL bias–variance bandwidth trade is physically meaningful.

**Major concerns.**
1. **Models are toys and horizons are tiny.** \(N=12\text{–}15\) steps,
   linearization about the origin, and "genuine nonlinearity that vanishes at the
   origin" (ADCS, distillation, PLL) — this is a regime where the eKF/MPC-from-IC
   is almost linear. Do the claimed nonlinear benefits persist for large
   excursions and realistic horizons? Show it.
2. **Tuned-to-order constants.** MTDC uses a uniform ring, droop leverage
   \(k_0=12\), and a fixed \(\gamma^2=8\) "chosen with margin"; the U-shape that
   drives the result is a consequence of these. Report sensitivity to \(\gamma\),
   \(k_0\), and topology (e.g. a non-uniform / meshed real DC grid). If the
   interior optimum disappears for reasonable alternatives, the study is fragile.
3. **Strawman baselines (again).** Over-provisioned uniform droop \(k=3\) is not
   what any operator runs. Compare to a tuned droop.
4. **No realistic benchmark.** For a power-systems claim I'd expect at least one
   standard test system (an IEEE bus system, or a published MTDC benchmark) rather
   than a synthetic ring.
5. **No figures whatsoever.** For an applications audience this is disqualifying on
   its own: I need closed-loop time responses (baseline vs. ContEst under the same
   disturbance), the cost breakdown, and the sensor-sparsity frontier plotted.
6. **Reproducibility is a placeholder.** The repository link is anonymized/pending;
   given the modeling concerns, the code and exact constants must be available to
   check.

**Minor.** State integration method/validity, per-unit conventions, and the
Monte-Carlo sample size and seeds. Clarify why distillation reports the open-loop
surrogate while the others report closed-loop Monte-Carlo — this inconsistency
undermines cross-study comparison.

**Recommendation: Major revision.** Add one realistic benchmark, sensitivity
studies to the manufactured constants, competent baselines, and — non-negotiable —
figures.

---

## Reviewer 5 — Numerical optimization / scientific computing

**Background & priorities.** Complexity, solver reproducibility, statistical
soundness of reported numbers, and honest scope for scalability.

**Summary assessment.** This is the most careful of the ContEst claims and I'm
broadly positive — the envelope route really does give the gradient at
\(O(r_x^2 r_\theta)\) on top of one inner solve, and the Riccati path is genuinely
scalable. My concerns are about *honesty of scope* and *statistical hygiene*
rather than correctness: the word "scalable" is doing more work than the
experiments support, and the Monte-Carlo numbers are reported without any
uncertainty.

**Strengths.**
- The \(O(r_x^6)\) (interior-point covariance/MPC) vs. \(O(r_x^3)\) (Riccati)
  analysis is clear and the crossover measurements are convincing.
- Reusing one inner solve's primal/dual for the gradient is the right design and is
  demonstrated to match finite differences tightly.

**Major concerns.**
1. **"Scalable" is oversold for the nonlinear path.** The eKF–MPC studies run at
   \(r_x = 6, 9, 12\); the information state is \(O(r_x^2)\) and the QP hits the
   same \(O(r_x^6)\) wall you diagnose. Only the *linear* Riccati path reaches
   hundreds of states. Say this in the abstract/intro, not only in the complexity
   subsection, so readers don't infer large-scale *nonlinear* co-design.
2. **Monte-Carlo results have no error bars.** Every reported \(J_\mathrm{est}\),
   \(J_\mathrm{det}\), reduction % is a point estimate with no sample size, no
   confidence interval, no seed disclosure. For a stochastic evaluation this is not
   acceptable — report \(N\), CIs (or shaded quantile bands in the figures you don't
   yet have), and confirm the reductions exceed Monte-Carlo noise. The PLL text even
   says a difference is "within Monte-Carlo noise" — quantify that noise everywhere.
3. **"Global to numerical tolerance" from 5 starts is not a certification.** Either
   report a coverage argument (e.g. many more starts, basin statistics) or downgrade
   the language to "best of five starts; all agreed."
4. **Timing methodology is thin.** One solver (Clarabel/Riccati), one machine, no
   package versions, no repetitions/variance on the timings. Give the standard
   reproducibility block.
5. **No figures.** At minimum plot the MTDC sparsity/performance frontier and the
   per-study cost breakdowns; these are cheap and would materially strengthen the
   numerical section.

**Minor.** The one-time gradient-check statement is good; add the FD step size and
that central differences were used. Multistart seed is given for one study — give
it for all.

**Recommendation: Minor revision.** The numerics are sound; fix the scope wording,
add statistical uncertainty and figures, and soften the global-optimum language.

---

## Cross-cutting synthesis (for the Editor)

### Consensus concerns (ranked by how many referees raised them)

1. **No figures (R4, R5; implied R1).** Unanimous among the empirically-minded
   referees; disqualifying for an applications audience. *Highest-value, lowest-cost
   fix.*
2. **Baselines are naive; the co-design thesis is under-tested (R2, R4).** The one
   competent-sequential comparison (PLL) is a ~0.3% null result; large reductions
   elsewhere are measured against strawmen. This is the deepest problem — it
   touches the paper's central claim.
3. **Nonlinear results optimize an unquantified surrogate (R1, R3).** "Exact
   gradient" is exact for \(\widehat J_c\), not \(J_\mathrm{stoch}\); the gap is
   never bounded and the eKF is never validated against a higher-fidelity filter.
4. **No statistical uncertainty on Monte-Carlo numbers; over-strong optimality
   language (R1, R5).** "Global to tolerance" from 5 starts; no CIs/seeds/\(N\).
5. **"Dual control" branding vs. delivered capability (R3, echoed by R2).** The
   method is caution-aware, not probing; retitle/reframe or demonstrate probing.
6. **Weak head-to-heads with the nearest methods (R2, R3).** Tzoumas 2021 and
   greedy/submodular sensor selection are cited as the closest work but never run.

### Points where the referees disagree

- **More theory vs. more application.** R1 wants the dual-invariance theorem *in
  the body*; R4 finds the paper already too theory-heavy relative to its toy
  models and wants space spent on realistic benchmarks and plots. The editor will
  have to arbitrate the page budget.
- **Severity.** R5 sees a minor revision (numerics are fine, presentation gaps);
  R2 sees a reject-and-resubmit (the thesis isn't demonstrated). The split is
  real: whether the paper's *machinery* or its *claim* is being judged.
- **What the headline should be.** R3 would foreground the estimator-design (PLL)
  angle; R2 would foreground a decisive joint-vs-sequential win (which doesn't yet
  exist); R5 would foreground the envelope-gradient complexity result.

### Suggested revision priorities (editor's digest)

1. Add **figures** (per-study cost breakdown with error bars; ≥1 closed-loop
   time response; the MTDC sparsity frontier). *(all)*
2. Re-baseline **every** study against a competent sequential (and alternating)
   pipeline; find or clearly concede where joint co-design decisively wins. *(R2, R4)*
3. Add **statistical rigor**: sample sizes, seeds, confidence intervals; soften
   "global optimum" to "best of N starts." *(R5, R1)*
4. **Bound or empirically probe the surrogate gap** and validate the eKF against a
   higher-fidelity filter on one study. *(R1, R3)*
5. **Reframe "dual control"** to "estimation-/caution-aware co-design," or add a
   probing example. *(R3)*
6. Run **real head-to-heads** (Tzoumas 2021; greedy sensor selection) on a common
   instance; report optimality gaps. *(R2, R3)*
7. State the **scalability scope** (nonlinear path is small-\(r_x\); only the
   Riccati path scales) up front, and add sensitivity to the tuned constants
   (\(\gamma\), \(k_0\), topology, horizon). *(R4, R5)*
8. Bring the **dual-invariance statement into the body** and tighten the
   surrogate/exactness language throughout. *(R1)*

**Net.** A strong methodological core (envelope-gradient co-design + LQG/H∞
regularity) wrapped in an empirical case that does not yet prove the paper's own
thesis. With competent baselines, figures, statistical rigor, and a bounded/validated
surrogate, this becomes a clear accept; without them, the referees split between
major revision and resubmission.
