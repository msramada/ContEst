# Referee report — "ContEst: Control and Estimation Co-Design …"

## Recommendation: **Major revision**

The paper makes a genuinely useful observation — that for the *value* of a convex
inner program the design gradient is available for free from the dual the solver
already returns, and that this unifies the LQG/H∞/eKF-MPC co-design variants under
one envelope-theorem identity. The framing (bringing the estimator/sensing axis
into control co-design) is well-motivated and the writing is mostly clear.
However, the central theorem's proof is only a sketch, one newly-featured
exactness result is asserted rather than proved, and the numerical section has
methodological confounds — including one example that contradicts the paper's own
stated evaluation protocol. These are fixable, but they are substantive.

---

## Major comments — mathematics and proofs

**M1. The proof of Theorem 1 (envelope gradient) is a heuristic, not a proof.**
The proof (l.491–502) differentiates
$V(\theta)=L(z^\star(\theta),\mu^\star(\theta),S^\star(\theta),\theta)$ via the
chain rule and asserts the terms in $\partial_\theta z^\star$ and
$\partial_\theta(\mu^\star,S^\star)$ vanish. Two problems:

- It invokes **differentiability of the solution map**
  $\theta\mapsto(z^\star,\mu^\star,S^\star)$. That is *stronger* than what you need
  and, more awkwardly, it is exactly the object your own selling point (l.533–540,
  "no implicit differentiation") says you avoid — establishing $C^1$-ness of the
  solution map is normally done through the IFT on the KKT system. This reads as
  internally inconsistent.
- For the **conic** term, $\nabla_S L=G(z^\star,\theta)\in\mathcal K$ is *not*
  zero; the claimed cancellation of $\langle \partial_\theta S^\star,
  G(z^\star,\theta)\rangle$ is the crux and is precisely what needs strict
  complementarity + the perturbation theory you cite. The one-line
  "$=\partial_\theta L$" glosses it.

Fix: either invoke Bonnans–Shapiro (Thm 4.24-type) directly and drop the
chain-rule hand-wave, or give the clean Danskin/saddle argument
$V(\theta)=\min_z\max_{S\in\mathcal K^\ast}L$ with a supporting-minorant
identification of the gradient. The latter needs only the *unique* primal–dual,
not a differentiable solution map.

**M2. Corollary "H∞ fixed-γ design: exact gradient" rests on an unproven
regularity claim.** You now (correctly) restructure the H∞ section around the
fixed-γ SDP being the exact case, but the corollary just states
"$\eqref{eq:Hinf_fixed}$ is strictly feasible with a nondegenerate, strictly
complementary dual, so Assumption 1 holds." There is no proof that the fixed-γ BRL
SDP has a unique/nondegenerate dual — in contrast to the H₂ case, where Appendix B
proves dual uniqueness constructively ($S_{11}^\star=P$). Since the ADCS control
gradient depends on this, it needs support. The cleanest route is already in your
hands: the fixed-γ GARE has a stabilizing solution $X(\theta)$ that is locally
$C^1$ for $\gamma$ above the optimal attenuation (Remark "Scalable game-Riccati
gradient"), so $V_{\rm wc}=\operatorname{tr}(XW_\theta)$ is differentiable; then
argue the SDP dual recovers that gradient. As written, the exactness is asserted.

**M3. Appendix C (subgradient existence) proves the right thing but by a detour.**
Your weak-duality bound gives a smooth *minorant*
$g(\theta')=\min_z L(z,\lambda^\star,\theta')$ with $g\le V$ and $g(\theta)=V(\theta)$;
that alone yields $\partial_\theta L(z^\star,\lambda^\star,\theta)=\nabla
g(\theta)\in\partial^- V(\theta)\subseteq\partial_{\rm Clarke}V(\theta)$ directly —
cleaner, and it *identifies the returned dual as the subgradient*, which is what
you want. The current "Lipschitz ⇒ Rademacher ⇒ pass to the limit" ending
(l.1601–1608) actually needs outer-semicontinuity of the dual map to conclude the
*specific returned* dual is a limit of nearby gradients, which you do not
establish. Recommend recasting Prop. 2's proof around the minorant.

**M4. Two of the four regularity assumptions are conflated in the abstract/intro.**
Line 142 says the gradient is exact "under a unique inner minimizer"; exactness
needs *both* a unique minimizer (A2) *and* a unique dual / strict complementarity
(A3), as your contribution bullet (l.148) and Assumption 1 correctly state. Align
the intro.

**M5. Info-state formulation vs. state feedback (l.275).** The infinite-horizon
inner problem is written $\min_\kappa \lim_k \mathbb E\,\ell(x_k,\kappa(x_k))$ —
full-state feedback — which contradicts the paper's thesis that the control is a
function of the *information state*, not $x_k$. Should be $\kappa(\mathcal Z_k)$ /
a function of the filtering density. Same care in eq. (mpc): the linearized
information-state dynamics are *affine* (a constant offset enters the covariance
recursion via $+W$); the offset (and its $\theta$-dependence if $W=W_\theta$) is
dropped from eq:gradmpc without comment.

---

## Major comments — numerical studies

**N1. The eKF-MPC baseline is confounded (PLL; and the eKF-ADCS variant).** For the
nonlinear studies the *baseline* runs a **certainty-equivalence controller** (eKF
mean only) at $\theta_{\rm nom}$, while the *optimum* runs the **information-state
MPC** at $\theta^\star$ (`example_pll.jl:147–153`, `example_adcs.jl:128–131`;
stated in passing at l.1371). So the reported reduction mixes two changes: the
design move $\theta_{\rm nom}\!\to\!\theta^\star$ **and** a controller-class upgrade
CE→info-state. To attribute the gain to *co-design*, the baseline must hold the
controller fixed (info-state MPC at $\theta_{\rm nom}$). As reported, the numbers
overstate the design contribution and are not an apples-to-apples design
comparison.

**N2. The distillation study contradicts the paper's stated protocol.** The paper
says (l.1004–1006, 1028–1033) that for the nonlinear studies it reports
**Monte-Carlo realized closed-loop** cost "to confirm each optimized design
improves the true cost, not merely the surrogate." But `example_distillation.jl`
(l.171–186) **skips MC entirely** ("MC sampling skipped for this study") and
reports the **open-loop MPC surrogate value** $J_c$ plus a design-time covariance
rollout, with $J_{\rm det}$ *defined as the residual* $J_c-J_{\rm est}$. Thus the
distillation reductions (22.8% / 1.9% / 6.5%) are exactly the surrogate/open-loop
numbers the paper claims not to rely on. Either run the MC or correct the
methodology text; as it stands this is a direct paper-vs-code contradiction, and
the 1.9% control figure in particular is unverified against the true cost.

**N3. No comparison against a competent *sequential* design.** The paper's
structural argument is that co-design captures an interior trade "a control-only
design cannot represent" (l.1122–1126). But every reduction is measured against a
naive/uniform $\theta_{\rm nom}$, and you honestly caveat (l.1374–1380) that the
baselines are "plausible rather than carefully hand-tuned." The decisive co-design
experiment is *joint vs. sequential*: place sensors first by a standard criterion
(A-optimal / greedy submodular) and design control second, then show the joint
optimum beats that pipeline. Without it, the experiments show only that
"optimizing beats not optimizing," not that *jointness* pays — which is the paper's
actual claim. This is the single most important addition, and it is directly on
point for an Automatica co-design paper.

**N4. Monte-Carlo reporting is incomplete.** The sample count $M$, horizon $T$, and
(critically) **confidence intervals** are not given in the text; a single seed
(20240624) is used. Please report $M$, standard errors/CIs, and preferably use
**common random numbers** (paired baseline/optimum) for a variance-reduced
comparison. Absent CIs, small effects (distillation control 1.9%) may be within
Monte-Carlo noise.

**N5. Public artifact does not match the manuscript.** The repository's "run all
paper examples" driver runs a **four-study** set — `example_adcs.jl` (eKF-MPC)
*and* `example_mtdc.jl` (HVDC) — whereas the paper describes **three** studies with
ADCS as the H∞+H₂ *SDP* study (`example_adcs_hinf.jl`). The reproducibility
artifact should reproduce the manuscript's studies.

---

## Minor / consistency

- **"Four" vs "three" studies:** contribution bullet (l.153) claims "four
  realistic co-design studies"; the body, Table 6, and summary all say three.
  (Consistent with N5 — a leftover from when MTDC was included.)
- **ADCS arithmetic (l.1119–1120):** $J_c$ is reported as $-31.9\%$ and
  $J_{\rm tot}$ as $-31.7\%$, both from an identical $0.224\!\to\!0.153$; they must
  be equal ($31.7\%$). Also the displayed 3-digit costs don't quite reproduce the
  headline percentages ($0.150\!\to\!0.091$ is $39.3\%$, not $39.5\%$;
  $0.074\!\to\!0.062$ is $16.2\%$, not $16.4\%$) — give one more digit or
  reconcile.
- **Notation:** eq. (mpc) uses $B_\theta$ but eq:gradmpc uses $\mathcal B_\theta$ —
  unify. `$J_{\mathrm{sc}^\star}$` brace typo (l.864). Typos: "accuratel"
  (l.294), "similarily", "the the".
- **H∞ filter "exact transpose dual"** (l.640–648): fine for what you do, but flag
  that you only combine *state-feedback* H∞ with H₂ estimation — H∞
  *output-feedback* is not separable, and the current wording could mislead.
- **Assumption A4** is called "a modeling convenience," but it is load-bearing for
  the $O(r_x^2 r_\theta)$ cost claim; say so.

---

## What is solid

The H₂/LQG regularity proof (Appendix B) is the strongest part — the constructive
dual-uniqueness argument giving $S_{11}^\star=P$ and $S_{11}'^\star=\Sigma^\varepsilon$
is correct and elegant, and it rightly grounds the "same object" claim linking the
LMI-dual and Riccati gradients. The complexity analysis ($O(r_x^6)$ IP scaling,
with the $2.6\text{s}\to30\text{s}$ timing matching $(30/20)^6$) is convincing, and
the honesty of Remark "Where exactness stops" about the eKF-MPC surrogate is
commendable — the issue in N2 is that the distillation *code* doesn't live up to
it.
