# ContEst Revision Plan: Addressing the Four Review Points

This document expands each review point into concrete actions, with specific suggestions for new text, results, and structure, keyed to the current draft (`automatica.pdf`, sections/equations as numbered there).

---

## Review Point 1 — Sharpen the theory contribution: assumptions, exactness limits, exact vs. surrogate-exact

**What the reviewer is really asking.** Right now, Assumption 1 is stated abstractly (Slater, uniqueness, nondegeneracy, strict complementarity), verified only *a posteriori* by finite differences, and the exact/surrogate boundary lives in a remark (Remark 8). The reviewer wants (i) hypotheses a reader can *check from the model data*, (ii) an unambiguous statement of which quantities are exact gradients of what, and (iii) the limits of that exactness stated as results, not asides.

### 1.1 Replace "assumed regularity" with verifiable sufficient conditions in the LQG case

This is the highest-value theory addition, and it is provable with classical tools. Add a lemma of the following form (as a new **Lemma/Proposition in Section 5**, full proof in an appendix — see Point 3):

> **Proposition (LQG regularity).** Suppose $(A_\theta, B_\theta)$ is stabilizable, $(A_\theta, C_\theta)$ is detectable, and $Q, R, W, V \succ 0$, with all data $C^1$ in $\theta$. Then for the control SDP (9) and the estimation SDP (10): (a) Slater's condition holds; (b) the primal minimizer is unique and given in closed form by the stabilizing DARE solution ($\Sigma^\star$ the closed-loop covariance, $K^\star$ the Riccati gain, and dually for the filter); (c) the dual of the Lyapunov LMI is unique and strictly complementary — in fact the $(1,1)$ dual block of the control LMI **is** the control Riccati solution $P$, and the $(1,1)$ dual block of the estimation LMI **is** the steady-state Kalman error covariance $\Sigma^\varepsilon$. Consequently Assumption 1 holds on a neighborhood $\Theta_0$ and Corollary 1 applies.

Proof ingredients (all standard, so the proof is compact): a strictly feasible point built from any stabilizing gain plus an $\varepsilon I$ inflation (Slater); a Lyapunov eigenvector argument showing every feasible point encodes a stabilizing gain ($W \succ 0$); the completion-of-squares identity $J(K) - \mathrm{tr}(PW) = \mathrm{tr}((K-K^\star)^\top(R+B^\top P B)(K-K^\star)\Sigma_K)$ for primal uniqueness; a null-space parametrization of the dual forced by complementary slackness, whose stationarity conditions reduce to the closed-loop Stein equation $\Lambda = Q + K^{\star\top}RK^\star + A_{cl}^\top \Lambda A_{cl}$, i.e. $\Lambda = P$ uniquely. Strict complementarity follows from $\mathrm{rank}$ counting with $P \succeq Q \succ 0$ and $R \succ 0$.

**Why this lands well.** It (a) removes the "verified a posteriori" weakness; (b) explains structurally *why* the paper's Section 9.6 Riccati gradient formulas ($\partial J_c/\partial W = P$, etc.) coincide with the LMI-dual formulas (11)–(14) — they are literally the same object; (c) gives a clean corollary: controllability + observability + positive-definite weights ⇒ everything in Sections 4–5 is exact and differentiable. Note in the text that controllability/observability imply the weaker stated hypotheses.

**Reflect in paper:** new Proposition + short proof sketch in Section 5 (3–5 lines), full proof in Appendix A; rewrite Corollary 1 to cite the proposition instead of "amounts to strict feasibility and nondegeneracy"; add one sentence to Remark 5 noting degeneracy can only occur when these hypotheses fail (e.g., singular $W$ or loss of detectability at isolated $\theta$).

### 1.2 State the H∞ regularity honestly and separately

Section 5.2 currently gestures at the problem (the minimal-γ BRL LMI is rank-deficient and its dual ill-conditioned) but Corollary 2 still invokes Assumption 1. Split the claim:

- **Fixed-γ game-Riccati route:** state a mirror regularity result — for fixed $\gamma$ above the optimal attenuation, $R + B^\top X B \succ 0$ and Schur $A_{cl}$ give a locally unique, $C^1$ stabilizing GARE solution (implicit function theorem on the Riccati residual; its linearization is the closed-loop Stein operator, invertible when $\rho(A_{cl}) < 1$), hence $V_{wc}(\theta)$ is $C^1$ with the closed-form gradient (22). This is what Section 9.5 actually uses, so make it the theorem.
- **Minimal-γ SDP route:** state explicitly, as a short remark, that at $\gamma = \gamma^\star$ strict complementarity generically *fails* and only Remark 5's subgradient statement is available. This turns an implementation footnote into a correct scoping of the theory — reviewers reward that.

### 1.3 Promote Remark 8 to a formal "exactness ledger"

Define notation that separates the three objects the paper currently overloads: the true stochastic value $J_{stoch}(\theta)$ of (2), the information-state surrogate $\widehat{J}_c(\theta)$ of (23), and the Monte-Carlo closed-loop cost reported in Section 9. Then state as a boxed proposition (or a small table in Section 6):

| Regime | Inner problem | Gradient (7)/(24) is exact for | Gap to true value |
|---|---|---|---|
| LQG (Sec. 5) | SDPs (9)–(10) / Riccati | $J_{cont}, J_{est}$ = true stationary LQG cost | none |
| H∞ fixed γ (Sec. 5.2) | GARE (20) | $V_{wc}(\theta)$ | none (given the ℓ2 disturbance model) |
| Nonlinear (Sec. 6) | convex program (23) | the surrogate $\widehat{J}_c(\theta)$ | eKF moment closure + linearization of the filter map + certainty-equivalent policy class |

Then make the Monte-Carlo evaluation of Section 9 do theoretical work: it *empirically bounds* the surrogate gap (optimized designs improve the true closed-loop cost, not just the surrogate). Add one sentence per study reporting surrogate value vs. Monte-Carlo value at $\theta^\star$ so the reader can see the gap size directly — this is cheap (you already have both numbers) and preempts the "you optimized a surrogate" objection.

### 1.4 Tighten the dual-control claim

The phrase "genuine stochastic optimal control (dual control) cost" is exposed: after linearizing the information-state dynamics about a nominal trajectory, the covariance trajectory in (23) may no longer respond to $u$, so the surrogate captures **caution** (the $\mathrm{tr}(Q\Sigma_{t|t})$ term) but not **probing**. Either (a) soften: reserve "dual control" for Problem 1/(2) and say the surrogate is a *caution-aware, certainty-equivalent-in-mean* approximation (cite Bar-Shalom–Tse's caution/probing taxonomy [13], Mesbah [38] for passive vs. active learning); or (b) keep the claim and demonstrate probing numerically in one small example. Option (a) costs one paragraph and closes the flank.

### 1.5 One sentence on the outer loop

Section 8 currently makes no convergence claim, and L-BFGS on the nonsmooth ℓ1 term is a known objection. Add: under the new Proposition (so $J$ is $C^1$ in the LQG/GARE regimes), projected gradient/quasi-Newton accumulation points are stationary; with the ℓ1 term, replace "adds the subgradient" by an explicit **proximal-gradient (ISTA-type) step on $\lambda_s\|\alpha\|_1$** with the standard stationarity guarantee, keeping L-BFGS on the smooth part. Two sentences plus one citation (e.g., Beck's *First-Order Methods*) fix this.

---

## Review Point 2 — Strengthen novelty positioning

**What the reviewer is really asking.** Section 2's "Positioning" paragraph argues in prose; the reviewer wants the deltas to be *falsifiable* — stated per literature, per capability, ideally with a head-to-head.

### 2.1 Add a positioning table

A compact table early in Section 2 (or replacing the "Positioning" paragraph) makes the claim auditable at a glance:

| Capability | CCD / MDO [25,33,8,37] | Sensor selection [34,56,59] | LQG sensing co-design [55] | Info-architecture LMIs [36,47] | Dual/adaptive MPC [26,38,32] | Diff. opt. layers [1,9,15] | **ContEst** |
|---|---|---|---|---|---|---|---|
| Plant/actuation designed | ✓ | ✗ | ✗ | partial | ✗ | ✗ | ✓ |
| Sensing/output map designed | ✗ | ✓ | ✓ (menu) | ✓ | ✗ | ✗ | ✓ (continuous + ℓ1 menu) |
| Estimator itself designed (e.g. PLL bandwidth, $V(\theta)$) | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✓ |
| Closed-loop control in the objective | ✓ | ✗ | ✓ | partial | ✓ (online) | ✓ | ✓ |
| Nonlinear regime | ✓ | ✗ | ✗ | ✗ | ✓ | ✓ | ✓ (eKF surrogate) |
| Gradient w.r.t. design | often FD | greedy/relaxation | menu enumeration | alternating LMIs | n/a | implicit diff. of KKT | exact from returned duals |

### 2.2 Sharpen the per-literature delta sentences

- **MDO:** don't just say ContEst "sits within" the MDO architecture literature — say *which* architecture it instantiates. ContEst is a nested (MDF-like) architecture whose discipline solve is convex and whose coupling derivative is supplied by duality instead of a sensitivity sweep; contrast with collaborative optimization [5], where consistency constraints replace the exact value gradient. One paragraph, and it converts the MDO citation block from décor into positioning.
- **Differentiable optimization:** make the complexity delta quantitative. Implicit differentiation of the KKT/cone map [1,2,9,15] solves an $m \times m$ linear system per design direction ($m = O(r_x^2)$ for the covariance SDPs), i.e. $O(r_x^6)$ per directional derivative, and differentiates the *minimizer*, which ContEst never needs; the envelope route reads the already-computed dual and contracts at $O(r_x^2 r_\theta)$. You already have the finite-difference comparison paragraph (Section 3); add one sentence extending it to implicit differentiation and cross-reference Section 9.6.
- **LQG sensing co-design [55] (Tzoumas et al.):** this is the nearest neighbor and deserves a *head-to-head experiment*, not only prose. On a linear benchmark with a discrete sensor menu, run (a) [55]'s separation-based co-selection, (b) ContEst's ℓ1 relaxation + threshold + re-solve, (c) exhaustive enumeration where feasible. Report cost and runtime. Even parity with [55] plus the added axes (continuous precisions, actuation, estimator tuning, nonlinear regime) is a strong positioning result; a win is stronger.
- **Dual/adaptive MPC:** the delta is offline hardware vs. online policy — [32,38] adapt the *inputs* for a fixed design; ContEst optimizes the *design* against the (approximated) closed-loop information dynamics. State it in exactly those words, and tie to the caution/probing precision of Point 1.4 so the two claims are consistent.
- **Sensor placement / OED:** you already note the A-optimal connection [42]; add the "least costly identification" line (Gevers–Bombois–Hjalmarsson) as the identification-community sibling of designing $V(\theta)$ against closed-loop cost — reviewers from that community will look for it.

### 2.3 Add a "what would break each baseline" sentence per study

Each Section 9 study already has a *Comments* paragraph; give each one a single crisp sentence of the form "a method of class X cannot represent this design variable because Y" (e.g., WAMS: menu-selection methods have no slot for a continuous shared-rate budget; PLL: no sensor-placement method can tune $V(\theta)$ through an estimator state in $\bar f$). This distributes the novelty argument onto evidence.

### 2.4 Consider one stronger empirical comparator

The current baselines are naive/uniform $\theta_{nom}$. Add, at least for WAMS, a **smart-sequential baseline** (optimize $\theta_h$ for a pure estimation criterion first, then $\theta_f$ for control with sensing frozen) and, ideally, an **alternating** baseline ($\theta_h$-step, $\theta_f$-step, iterate). Beating uniform is expected; beating a competent sequential pipeline *is the thesis of the paper* ("the sequential solution need not be a stationary point"). This is the single most persuasive novelty evidence you can add for the cost of a few runs.

---

## Review Point 3 — Move heavy derivations to appendices; focus the body

**What the reviewer is really asking.** The body currently interleaves framework, two full inner-problem derivations (H2 and H∞), five study formulations with complete state-space models and constants, and a complexity analysis. Automatica readers want the contribution arc in ~10 body pages.

### 3.1 Proposed reorganization

**Keep in body (trimmed):**
- Sections 1–4 nearly as-is (Section 2 compressed by the positioning table of Point 2.1; the finite-difference paragraph in Section 3 can shrink to 3 sentences with the details in the complexity appendix).
- Section 5: the two SDPs, the new regularity Proposition (statement + 3-line proof sketch), gradient formulas (11)–(14), and the weights/noise subsection 5.1 (it is short and load-bearing for the PLL study).
- Section 6: Lemma 6, program (23), costate gradient (24), and the exactness table of Point 1.3.
- Section 7 and 8: essentially as-is (Section 8 gains the two proximal-step sentences of Point 1.5).
- Section 9: per study keep *Significance* (trimmed to ~5 lines), the design-variable definition, the results table, one figure, and *Comments*; everything else moves out.

**Move to appendices:**
- **Appendix A** — Proof of the LQG regularity Proposition (Point 1.1) and derivation of the dual identities $S_{11}^\star = P$, $S_{11}'^\star = \Sigma^\varepsilon$; this also *replaces* the currently unproven pointer to [6,48] in Corollary 1.
- **Appendix B** — The H∞ development: derivation of (16)–(17) from the bounded-real lemma, the game-Riccati route (20)–(22), and the fixed-γ regularity result (Point 1.2). Keep in the body only a half-page Section 5.2 stating: robust inner value exists, same conic form, same envelope gradient, Riccati preferred for numerics — with forward references.
- **Appendix C** — Study details: the full state-space models (28)–(34), all constants, noise levels, priors, and design-cost coefficients, one subsection per study. The body keeps only the equations that define the design coupling (e.g., the PLL line (31c) and $V(\theta)$ (32), the WAMS rate map $V_j = v_0/r_j$).
- **Appendix D** — The complexity study (current Section 9.6) minus its "Takeaway" paragraph, which stays in the body (it is a genuine contribution message). Timing methodology, solver versions, and the $O(r_x^6)$ vs. $O(r_x^3)$ measurements go here.

### 3.2 Practical notes

- Automatica permits appendices in the same manuscript; keep each appendix opened by one sentence saying what the body statement it supports is, so appendices remain skimmable.
- After the move, re-check every `\ref`/`\eqref`: equation numbers (11)–(34) will shift; do the renumber last.
- The word budget freed (~3 pages) is exactly what Points 1.1, 2.1, and 2.4 consume — net length roughly constant, density of contribution much higher.

---

## Review Point 4 — Placeholders, build issues, publication-ready tables/figures

A concrete checklist against the current build:

### 4.1 Build artifacts
- **"1.0.0.1 Contributions"** in Section 1 is a sectioning bug — a `\paragraph`/`\subsubsection` numbered because `secnumdepth` is too deep, or a stray `\subsubsection` under an unnumbered structure. Fix with `\paragraph*{Contributions.}` or a bold run-in macro (`\noindent\textbf{Contributions.}`), and grep the source for any other auto-numbered run-in heads.
- Compile with `-halt-on-error` and clear **every** warning class: undefined references, multiply-defined labels, overfull `\hbox` (the wide 4×4 LMIs (16)–(17) are prime suspects in Automatica's double-column format — consider `\resizebox` on the LMI or the appendix move of Point 3.1, which solves it for free), and missing citations.
- Verify the Automatica class (`autart`/elsarticle-style) requirements: running title, corresponding-author footnote format, keyword count, and the reference style (Automatica uses author–year `elsarticle-harv`; the current numeric bracketed citations suggest the wrong bibliography style is loaded — this alone can trigger a desk revision).

### 4.2 Figures
- The current text references no figures at all — earlier drafts had per-study closed-loop responses and cost breakdowns. Reinstate **one figure per study** (two panels: (a) closed-loop response baseline vs. ContEst under the *same noise realization*, (b) cost breakdown $J_{est}/J_{det}/J_c/J_{tot}$), plus the ℓ1/budget sweep where applicable. For the Monte-Carlo results (Point 1.3), add error bars or shaded quantile bands to the cost-breakdown bars — this makes the stochastic evaluation visibly rigorous.
- Uniform figure hygiene: identical fonts (match body font size after scaling), labeled axes with units (s, rad, K, p.u.), consistent baseline/ContEst colors and line styles across all studies, vector PDF export, no titles inside panels (captions carry the message).
- Captions should state the takeaway, not just the content ("ContEst reallocates the fixed rate budget to area 2, cutting $J_{est}$ 10.5%"), Automatica-style.

### 4.3 Tables
- Consistent significant figures within each table (Table 1 mixes 0.19/0.93 with 104.9/319.8 — fine, but be uniform per column) and consistent percent formatting in Table 5.
- Table 1's $J_c$ baseline (319.3) vs. $J_{tot}$ baseline (319.8) with $J_{des} = 0.5$: correct, but a reader will pause — add $J_c$ rows explicitly to Tables 1–4 or a footnote $J_{tot} = J_c + J_{des}$ so every column is self-checking.
- Every table gets: units where applicable, a caption that names the baseline precisely, `booktabs` rules (no vertical lines), and a pointer to the appendix subsection holding the constants.
- The text promises "specific [FD] values are reported per case" — verify each study actually reports its gradient-check number (they do now; keep it that way after the appendix move, one number per study, e.g. in the results table footer).

### 4.4 Consistency sweep
- Terminology: pick one of "lower level"/"inner problem"/"second stage" and use it throughout (all three currently appear).
- Notation: $\ell$ vs. $\ell_s$ (fixed in this draft — keep it), $J_c$ vs. $\widehat{J}_c$ after Point 1.3's notation split, $\theta_h$ entering $\bar h$ vs. $V(\theta)$ per study (Section 9 preamble already summarizes this — turn that paragraph into a small "design-axis map" table for instant readability).
- Cross-references: Section 9.5 refers to "Section 9.6" complexity results and vice versa; re-verify after reorganization. Abstract still says "duals of Lyapunov constraints … costates" — extend it with one clause mentioning the robust (H∞) instantiation, which is new in this revision and currently invisible from the abstract.

---

## Suggested priority order

| Priority | Item | Effort | Review points served |
|---|---|---|---|
| 1 | LQG regularity Proposition + Appendix A proof | medium (math is classical; we have the proof skeleton) | 1 |
| 2 | Appendix reorganization + renumber | low | 3, 4 |
| 3 | Positioning table + sharpened deltas + abstract clause | low | 2 |
| 4 | Smart-sequential/alternating baseline (at least WAMS) + [55] head-to-head | medium | 2 |
| 5 | Exactness ledger + dual-control precision + surrogate-vs-MC numbers | low | 1 |
| 6 | Figures reinstated with MC error bars; table/build checklist | medium | 4 |
| 7 | H∞ fixed-γ regularity statement; proximal outer-step sentences | low | 1 |

Each numbered action above maps one-to-one onto a bullet in the response-to-reviewers letter, which makes the rebuttal easy to write: quote the review point, cite the new Proposition/appendix/table by label, done.