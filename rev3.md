# Review of "Control and Estimation Co-Design (ContEst)" — Reviewer 3

**Venue:** *Automatica* (regular paper)
**Recommendation:** Major revision.

The paper is well-motivated, the core idea (read the design gradient of an
inner control/estimation value off the duals the solver already returns, and
drive a first-order outer loop with it) is clean and, in the LQG/H∞ regimes,
correct and elegant. The unification of the control and estimation halves under
one envelope identity, and the honest "exactness ledger" separating the true
value / convex surrogate / realized cost, are genuine strengths. However, before
it is publishable I believe several issues of **claim calibration, experimental
rigor, and technical precision** must be addressed. None appears fatal; most are
about tightening claims and making the empirical case defensible to a skeptical
control audience.

---

## 1. Summary of the contribution

ContEst poses joint design of actuation (`θ_f`, entering `f`) and
sensing/estimation (`θ_h`, entering `h` or the noise level `V`) as a two-stage
stochastic program whose second-stage value is a stochastic optimal-control
(dual-control) cost expressed through the information state. The design gradient
is obtained from the inner solver's duals via the envelope theorem, so the outer
BFGS loop needs no differentiation through the optimizer. The framework is
instantiated in LQG (a pair of covariance SDPs / Riccati twins), robust H∞
(bounded-real-lemma SDPs and a fixed-γ game-Riccati route), and a nonlinear
eKF–MPC surrogate. Sensor allocation enters as a mixed-integer or ℓ1-relaxed
term. Four studies (spacecraft ADCS, distillation, PLL/DSE grid, multi-terminal
HVDC) report 6–50% total-cost reductions.

---

## 2. Major comments

### M1. The headline improvements conflate two effects — add controller-matched ablations.
For the three eKF–MPC studies, the **baseline** is a certainty-equivalence (CE)
controller on the eKF mean at `θ_nom`, while the **ContEst** point is the
covariance-aware information-state MPC at `θ*`. The reported reduction therefore
mixes (a) the re-designed parameters — the actual co-design contribution — with
(b) a controller upgrade (CE → info-state MPC). The text acknowledges this
(Sec. 7, "reflects both the re-designed parameters and the estimation-aware
control they enable"), but for a co-design paper the co-design gain must be
isolated. **Please add an ablation** holding the controller fixed:
`(θ_nom, info-state MPC)` vs `(θ*, info-state MPC)`, and ideally also
`(θ_nom, CE)` vs `(θ*, CE)`. Report the decomposition. Without this, a reviewer
cannot tell how much of the 7.7%/6.5%/34.8% is co-design versus a better online
controller.

### M2. The "naive baseline" is a straw man — add a competent sequential baseline to every study.
Only the PLL study (Sec. 7.3, "Beyond the naive baseline") compares against a
sequential estimator-then-control pipeline, and there the joint optimum beats it
by only ~0.3% on the surrogate. This is the honest and important number. The
other three studies compare only against `θ_nom` (uniform budget / naive
placement / over-provisioned droop), which inflates the apparent co-design
benefit. **Please report the sequential-design baseline (tune `θ_h` for an
estimation criterion, then `θ_f` with `θ_h` frozen, optionally alternating) for
ADCS, distillation, and MTDC as well.** If, as for PLL, the gap is small, say so
plainly — that is a scientifically valuable finding about *when* co-design pays,
and it strengthens rather than weakens the paper. As written, three of four
studies risk the criticism "you only beat a design nobody would deploy."

### M3. No uncertainty quantification on Monte-Carlo results.
The eKF–MPC costs are Monte-Carlo averages of the realized closed-loop cost, but
the number of trajectories, the sampling protocol, and **confidence intervals /
standard errors are not reported.** For ADCS (7.7% total) and distillation
(6.5%), the improvement is small enough that statistical significance is not
self-evident. Please report `N_MC`, seed protocol, and 95% CIs (or standard
errors) on every reported cost and reduction. A 6.5% reduction without an error
bar is not a defensible claim at *Automatica*.

### M4. Inconsistent reporting metric across studies.
Distillation reports the **open-loop surrogate** `Ĵ_c` (with
`J_est = Σ_t tr(QΣ_{t|t})` from a covariance rollout) rather than the sampled
closed-loop cost used for ADCS/PLL, "because the many-stage information-state MPC
is expensive to simulate." The surrogate is admittedly optimistic. This makes the
distillation number not directly comparable to the others and potentially
unrealized in closed loop. Please either (i) simulate the distillation study in
closed loop like the others (even at reduced `N_MC`), or (ii) explicitly mark the
distillation row as a surrogate-only number in Tables 6–7 and the abstract range,
and verify at least that the *ranking* (baseline vs `θ*`) is preserved in a small
closed-loop check.

### M5. Global-optimality language is too strong.
Two passages state that agreement of 5 multi-start runs is "strong empirical
evidence that it is globally optimal over Θ." Five starts is weak evidence for
global optimality of a nonconvex problem; the distillation study itself finds
three minima with five starts, which undercuts the claim elsewhere. The paper
already contains the correct hedge ("multi-start alone certifies only that it is
the best of the starts"). Please **delete the "strong evidence of global
optimality" phrasing** (Sec. 7 intro and Sec. 7.3) and keep only the hedged
statement, or substantially increase the number of restarts and report the
basin-of-attraction statistics.

### M6. Regularity hypotheses `Q,R,W,V ≻ 0` are violated by two studies.
Proposition 1 (LQG regularity) and its strict-complementarity argument
(App. A, "`P ⪰ Q ≻ 0`") assume `Q ≻ 0`. But the PLL study uses
`Q = blkdiag(15I, 120I, 0_3)` (PLL states unweighted) and the MTDC study
penalizes only `ΔV`, so `Q` is **singular** in both. This is exactly a case where
(A3) can fail and the exact-gradient guarantee degrades to the subgradient regime
of Remark 2. Two requests: (i) In Remark 2, add "singular `Q`" to the list of
degeneracy sources (currently only singular `W` / lost detectability are named).
(ii) State explicitly, per study, whether the exactness hypotheses hold, or add a
small Tikhonov term to `Q` and confirm the results are unchanged. As written,
there is a gap between the theorems (which need `Q ≻ 0`) and the experiments
(which use `Q ⪰ 0`).

### M7. The "second stage = stochastic optimal control" is only ever a tractable surrogate — make this crisper up front.
The abstract and introduction repeatedly call the second stage "the stochastic
optimal control problem." The true dual-control problem
(Eqs. 3.1a–b via the Bayesian filter) is intractable and is *never* solved: it is
replaced by exact convex programs (LQG/H∞) or a linearized eKF surrogate. The
exactness ledger (Table 5) handles this beautifully in Sec. 6, but the abstract
still reads as if the stochastic control problem itself is solved. Please add one
clause to the abstract/intro clarifying that `J_stoch` is instantiated by a
tractable convex inner value (exact in LQG/H∞, an eKF surrogate otherwise).

### M8. Cost-of-gradient overclaims.
Contribution 2 says gradient computation "comes with no (or trivial) cost" and
the abstract calls it "essentially free," yet Sec. 7.6 correctly gives the
contraction cost as `O(r_x^2 r_θ)`. Please reconcile: replace "no cost" with
"negligible relative to the inner solve," consistent with Sec. 7.6. This is a
one-line fix but the current phrasing invites an easy reviewer objection.

---

## 3. Technical / correctness issues

### T1. Duplicate equation label `eq:Jstoch`.
Eqs. (3.1) finite-horizon and infinite-horizon **both carry `\label{eq:Jstoch}`**
(lines ~310 and ~315). This produces a "multiply-defined label" warning and makes
every `\eqref{eq:Jstoch}` ambiguous (Problem 1, Lemma 2 proof, etc.). Give them
distinct labels (e.g. `eq:Jstoch-fh`, `eq:Jstoch-ih`) and point each reference to
the intended one.

### T2. The infinite-horizon cost equation is malformed.
Eq. (3.1, infinite horizon):
`J = min_κ lim_{N→∞} E{ ℓ(x_k, κ(x_k) + ℓ^N(x_N) }` has (i) a missing closing
parenthesis after `κ(x_k)`; (ii) a free index `k` under a limit in `N` with no
sum; (iii) a terminal cost `ℓ^N(x_N)` inside an infinite-horizon limit, which is
not meaningful. Presumably the intended object is an average cost
`min_κ lim_{N→∞} (1/N) E Σ_{k=0}^{N-1} ℓ(x_k, κ(x_k))`. Please rewrite.

### T3. Finite-horizon sum double-counts the terminal stage.
Eq. (3.1, finite horizon) sums `Σ_{k=0}^{N} ℓ(x_k,u_k) + ℓ^N(x_N)`, but `u_k` is
only defined for `k = 0,…,N-1`, so the `k=N` stage term `ℓ(x_N,u_N)` is undefined
and the terminal cost appears twice. Should be `Σ_{k=0}^{N-1} ℓ + ℓ^N`.

### T4. Open-loop vs feedback inconsistency in the second stage.
The "bilevel → two-stage" narrative (Sec. 3) states the control is the *recourse*
reacting to information `Z_k`, i.e. a feedback policy, but the finite-horizon
`J_stoch` minimizes over the **open-loop** sequence `u_{0:N-1}`, while the
infinite-horizon one minimizes over a policy `κ`. Please make the second stage a
policy `u_k = κ_k(Z_k)` consistently (or clarify that the eKF–MPC receding-horizon
implementation supplies the feedback and the written program is the per-step
open-loop value). As stated, the formulation and the recourse story disagree.

### T5. Subgradient regime attributed to the wrong non-uniqueness.
The abstract ("subgradient in the worst cases where the **minimizer** is
non-unique") and intro (line ~146) attribute the subgradient case to a
non-unique *primal minimizer* (A2). But the theory (Assumption (A3), Remark 2,
App. C) makes clear the exact-gradient formula fails when the **dual** is
non-unique, i.e. (A3) — strict complementarity/nondegeneracy — fails. Non-unique
primal (A2 failing) is a separate issue. Please correct the abstract/intro to
"non-unique dual" (or "when strict complementarity fails").

### T6. Symbol collision `X` in the MTDC study.
In Sec. 7.4, `X` denotes both the **H∞ game-Riccati solution** (`J_det = tr(XW)`,
Eq. via game GARE) and the **estimation posterior covariance decision variable**
in Eq. (7.x) `min_{X⪰0} tr(QX) + λ_s 1ᵀα`. Same symbol, same subsection, two
different objects — this is genuinely confusing. Rename the estimation covariance
(e.g. `Π` or `Σ_post`).

### T7. Cap notation inconsistent within the MTDC study.
The estimation SDP writes the cap as `[X]_{2i,2i} ≤ σ̄_i²`; the prose (item (C))
and the Table 4 caption write `[Σ]_{2i,2i} ≤ tol`. Unify the variable (`X` vs
`Σ`) and the bound symbol (`σ̄_i²` vs `tol`).

### T8. Envelope "Theorem" is a known result — consider status/attribution.
Theorem 1 (`prop:envelope`) is essentially the Danskin/Bonnans–Shapiro value
sensitivity for convex conic programs (correctly cited). Presenting a classical
result as a numbered "Theorem" of the paper may draw fire. Either downgrade to
"Proposition (restated from [Bonnans–Shapiro])" or add a sentence stating the
theorem is classical and the contribution is its *uniform application* to the
control+estimation halves. (The Related Work already says this; mirror it at the
theorem.)

### T9. Feasibility of the design box Θ.
The inner SDPs/Riccati require stabilizability/detectability at the evaluated `θ`.
If the box `Θ` contains points where the plant loses stabilizability (e.g. a
droop or damping driven to a degenerate value), the inner problem is infeasible
and the outer gradient undefined. Please state how infeasible `θ` are handled
(projection, barrier, or a guarantee that `Θ` is contained in the
stabilizable/detectable set for each study).

### T10. γ² = 8 is fixed without justification (MTDC / H∞).
The robust study fixes `γ² = 8`. Since `γ` shifts both the worst-case control cost
and (through `A(θ)`) the estimation prior, the reported optimum and the 38–61%
reductions depend on it. Please justify the choice and show a short sensitivity
sweep of the co-design outcome vs `γ` (or at least argue the qualitative
conclusions are `γ`-robust).

---

## 4. Presentation and organization

### P1. Flow: probability → abstract convexity → concrete problems.
Section 3.1 (Bayesian filter + smoothing lemma) sits between the deterministic
problem statement and the abstract envelope section, and its payoff (cost depends
on the filtering density → estimation enters) is only used later in the LQG
separation and Lemma 4. The register switches probability → abstract conic
optimization → concrete filters in three consecutive moves. Consider (i) moving
§3.1 to lead into the LQG/eKF instantiations where it is consumed, (ii) adding a
one-paragraph roadmap at the end of §3 and a bridge sentence opening §4, and
(iii) demoting the smoothing lemma (Lemma 1, the tower property) to an inline
citation. This is presentational, but it would materially improve readability.

### P2. Lemma 1 (Smoothing) is textbook.
The tower property stated as a formal lemma interrupts the development. Fold it
into the proof of Lemma 2 as a cited step.

### P3. Abstract omits scope items.
The abstract no longer mentions sensor allocation or the number/breadth of case
studies. A one-clause addition ("...and sensor selection via an ℓ1 relaxation,
validated on four studies across aerospace, process, and power systems") would
better set expectations.

### P4. Transparency about modest gains.
Two of four studies show <8% total improvement (ADCS 7.7%, distillation 6.5%),
and the headline "6–50%" is driven by PLL and MTDC. The paper is commendably
honest, but the discussion (or conclusion) should state plainly *when* ContEst
delivers large gains (strong `θ`-coupling, shared budgets, estimator tuning) vs
small ones (near-separable problems), so readers can predict applicability. This
also pre-empts the "sometimes co-design barely helps" objection.

### P5. Mixed metrics in the abstract range.
The "6–50%" range mixes sampled closed-loop reductions (eKF studies) with the
steady-state H∞/covariance-SDP number (MTDC 49.7%). Please note these are on
different scales, or give the ranges separately (e.g. "6–35% closed-loop across
the three nonlinear studies; ~50% on the robust steady-state HVDC study").

---

## 5. Minor / typographical

- **Contribution bullets** lack consistent terminal punctuation (bullets 1, 3
  end without a period). Standardize.
- Line ~128: "the sequential solution (if exists)" → "(if it exists)".
- Line ~154: "the gradients computation comes with no (or trivial) cost." →
  rephrase (see M8) and fix grammar ("gradient computation").
- Line ~165: "as a pair of or more semidefinite programs" → "as a pair (or more)
  of SDPs".
- Line ~302: "The state dynamics `f_θ` is locally bounded" → subject/verb; and
  clarify *what* is bounded (Lipschitz? bounded on compact sets?). "Locally
  bounded" for a dynamics map is unusual phrasing — did you mean locally
  Lipschitz, or that trajectories stay bounded?
- Line ~308: "some mechanical, thermal, structural,...etc considerations" —
  informal "...etc"; rewrite as "(mechanical, thermal, structural, etc.)".
- Sec. 7.6 / App. A cross-check the year of `grigoriadis1998integrate`: the
  citation key says 1998 but the commented BibTeX at the end of the file (lines
  ~2109–2114) says 1996 and is marked "VERIFY". Resolve the year and remove the
  placeholder comment block before submission.
- The GitHub URL is an anonymized placeholder — fine for review, but ensure a
  data/code-availability statement is present for the final version.
- Table 3 (positioning) `✓/✗` entries: a few are debatable (e.g. dual/adaptive
  control arguably does "design the estimator online"). Consider a footnote
  defining each column crisply to avoid disputes; "estimator itself designed" is
  defined, but "closed-loop control in objective" and "nonlinear regime" are not.
- Eq. (3.1) expectation is "over `(x_0, w_{0:N-1}, v_{0:N})`" — confirm the `v`
  index range matches the measurement indices used (`y_{0:N}`).

---

## 6. Questions to the authors

1. In the eKF–MPC studies, what fraction of each reported reduction is due to the
   controller change (CE → info-state MPC) vs the parameter redesign? (See M1.)
2. How many Monte-Carlo trajectories underlie each cost, and what are the
   confidence intervals? (See M3.)
3. Against a competent sequential (estimation-then-control) baseline, what are the
   co-design gains for ADCS, distillation, and MTDC? (See M2.)
4. Do the exactness hypotheses (`Q,R,W,V ≻ 0`) hold in the PLL and MTDC studies
   given their singular `Q`? If not, are those runs in the subgradient regime, and
   does that affect convergence? (See M6.)
5. How sensitive is the MTDC co-design to the fixed `γ² = 8`? (See T10.)
6. How is an infeasible `θ` (loss of stabilizability/detectability within `Θ`)
   handled during the outer search? (See T9.)
7. The gradient check reports relative error "below `10⁻⁴` (best `1.6×10⁻⁶`)."
   Which regime gives `10⁻⁴`? For the exact SDP/Riccati cases one expects
   agreement near the finite-difference floor (`~10⁻⁶`–`10⁻⁸`); a `10⁻⁴`
   discrepancy there would suggest an implementation issue rather than the eKF
   surrogate. Please report the check per regime and state the step sizes.

---

## 7. Overall

The methodology is sound and the writing is (mostly) clear; the envelope-from-duals
idea is a nice, reusable observation and the LQG/H∞ derivations and their
appendix proofs are careful. The weaknesses are almost entirely in **empirical
calibration**: isolate the co-design effect (M1), use fair baselines (M2), quantify
uncertainty (M3), and reconcile the theory's positivity hypotheses with the
experiments (M6). Fixing the labeling/formula bugs (T1–T3, T5–T7) and softening
the global-optimality and "free-gradient" language (M5, M8) are quick but
important. With these addressed, this would be a solid *Automatica* contribution.
