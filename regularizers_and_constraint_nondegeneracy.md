# Regularizers, LMI constraints, and constraint (non)degeneracy

This note expands on which regularizers / LMI constraints can be added to the inner
SDP while keeping the envelope‑theorem assumptions **A1–A4**, and then explains in
detail **what constraint (non)degeneracy is** and **why degeneracy breaks the exact
gradient.**

---

## 1. Recap — the admissible class

Split additions by *where* they live.

- **Design‑cost side** $J_{\mathrm{des}}(\theta)$: any convex penalty on $\theta$
  (sparsity, group‑sparsity, budgets). These never touch the inner SDP's A1–A4;
  nonsmooth ones are handled by the outer prox/epigraph. Always compatible.

- **Inner SDP side:** A1 (convex, $C^1$ data after conic lifting, Slater, attainment)
  and A4 ($\theta$ enters through finitely many $C^1$ data) admit essentially **any
  convex, conic‑representable term whose data are $C^1$ in $\theta$ and that keeps a
  strictly feasible point**: trace terms, gain norms, nuclear norm, log‑det/Tikhonov,
  variance/output/ellipsoid caps, D‑stability (pole‑region) LMIs, polytopic quadratic
  stability, gain‑scheduling/LPV LMIs, passivity/IQC, structural gain patterns.

**A2 and A3 are the gate.** They are *not* automatic under convex additions. Two
universal recipes make them hold for the whole family at once:

1. **A strictly convex regularizer** in the variables — a small Tikhonov
   $\varepsilon\lVert\cdot\rVert^2$ or $-\varepsilon\log\det\Sigma$ — makes the inner
   objective strictly convex, forcing **A2** (unique primal) and generically restoring
   A3.
2. **Nondegenerate active constraints** — the subject of the rest of this note — give
   **A3** (unique dual, differentiable value).

The augmented envelope gradient simply gains one term per constraint,
$\nabla_\theta V=\partial_\theta\phi+\sum_i\langle\lambda_i^\star,\partial_\theta(\text{data}_i)\rangle$.

---

## 2. What the assumptions really need

Restated in the language that matters for perturbation:

- **A2** = the primal minimizer is unique.
- **A3** = the dual multiplier is unique **and** strict complementarity holds, so the
  optimal value $V(\theta)$ is differentiable and $\nabla V=\partial_\theta\mathcal L$
  at that single multiplier.

The classical SDP sensitivity results (Alizadeh–Haeberly–Overton 1997; Shapiro 1997;
Bonnans–Shapiro 2000) tie these to three geometric properties of the primal–dual pair:

$$
\textbf{dual nondegeneracy}\ \Rightarrow\ \text{unique primal (A2)},\qquad
\textbf{primal nondegeneracy}\ \Rightarrow\ \text{unique dual (A3)},
$$
$$
\textbf{strict complementarity}\ +\ \text{nondegeneracy}\ \Rightarrow\ V\ \text{differentiable}.
$$

So "constraint nondegeneracy" is exactly the property that buys A3 (and its dual
buys A2).

---

## 3. What "nondegenerate" means

### 3a. Scalar inequality constraints — LICQ

For a program $\min_z f(z)$ s.t. $g_i(z)\le 0$, let $\mathcal A=\{i:g_i(z^\star)=0\}$
be the **active set** at the optimum. The **linear‑independence constraint
qualification (LICQ)** is:

$$
\{\nabla g_i(z^\star)\}_{i\in\mathcal A}\quad\text{are linearly independent.}
$$

Under LICQ the KKT multiplier $\lambda^\star$ is **unique**. LICQ *failing* — active
constraint gradients linearly dependent — is exactly (scalar) degeneracy.

### 3b. LMI / conic constraints — primal & dual nondegeneracy

For an LMI $M(z)\succeq 0$ that is active at $z^\star$ (i.e. $M(z^\star)$ is singular),
degeneracy is more subtle because the "active set" is a *face* of the PSD cone, not a
finite index set. Let the columns of $Q$ span $\ker M(z^\star)$ (dimension = the rank
drop). Following Alizadeh–Haeberly–Overton / Shapiro:

- **Primal nondegeneracy:** the linear map $z\mapsto Q^\top\, \big(D_zM(z^\star)\big)\, Q$
  is **surjective** onto the symmetric matrices on $\ker M$. Intuitively: the
  constraint's directions "reach" every way the kernel could move — the active face is
  hit *transversally*. This is the LMI analogue of LICQ, and it implies the **dual is
  unique**.
- **Dual nondegeneracy:** the analogous surjectivity on the dual side; it implies the
  **primal is unique**.
- **Strict complementarity:** the primal slack $M(z^\star)$ and its multiplier $S^\star$
  have complementary full rank,
  $$
  \operatorname{rank}M(z^\star)+\operatorname{rank}S^\star=\dim .
  $$
  Their ranges are orthogonal complements (no "shared" kernel direction).

A constraint is **degenerate** if any of these fails: dependent active directions
(primal), or a rank overlap $\operatorname{rank}M^\star+\operatorname{rank}S^\star<\dim$
(loss of strict complementarity).

> In the LQG appendix this is exactly why A3 goes through: at the optimum
> $\operatorname{rank}M_2^\star=n$, its multiplier has $\operatorname{rank}S^\star=n$
> (core $\Xi=P\succ0$), and $n+n=2n$ = full — strict complementarity holds, and the
> lifted form makes the dual a single point.

---

## 4. What a degenerate constraint looks like

Concrete, minimal examples.

**(i) Redundant / duplicated constraints (dependent gradients).** Add the same
variance cap twice, or include a polytope vertex dominated by the others in a
quadratic‑stability set. At the optimum both are "active" but their normals are
linearly dependent, so the multiplier splits arbitrarily — the dual is a line, not a
point.

**(ii) A "tie" between two active constraints.** The canonical parametric example:
$$
V(\theta)=\min_x\ x\quad\text{s.t.}\quad x\ge\theta_1,\ x\ge\theta_2,
\qquad x^\star=\max(\theta_1,\theta_2).
$$
When $\theta_1=\theta_2$ both are active; the multiplier set is the whole simplex
$\{\lambda\ge0:\lambda_1+\lambda_2=1\}$ — **non‑unique** — and $V(\theta)=\max(\theta_1,\theta_2)$
has a **kink** there. Away from the tie exactly one constraint is active,
nondegenerate, $V$ smooth. Degeneracy = the tie.

**(iii) SDP rank overlap (loss of strict complementarity).** The minimal‑$\gamma$
$H_\infty$ synthesis LMI at $\gamma=\gamma^\star$: the certificate makes the
bounded‑real matrix drop rank *and* the dual drops rank so that
$\operatorname{rank}M^\star+\operatorname{rank}S^\star<\dim$. The optimal face is not a
vertex; the dual is non‑unique and $\gamma^\star(\theta)$ is nonsmooth. This is
example (ii) lifted to matrices (a "tie" between worst‑case frequencies).

**(iv) Structural gain equalities.** Coupled decentralized‑$K$ patterns can make the
active equality constraints linearly dependent (a fixed‑pattern gain that is already
implied by others), again splitting the multiplier.

---

## 5. Why degeneracy is a problem

Everything the envelope theorem promises rests on a **unique** multiplier and a
**smooth** value. Degeneracy destroys both:

1. **The dual is not a singleton.** With multiplier set $\mathcal D^\star(\theta)$ not a
   point, "read the multiplier and contract it against $\partial_\theta(\text{data})$"
   is **ambiguous** — different solvers (or the same solver on a different run) return
   different $\lambda^\star\in\mathcal D^\star$, giving different gradients.

2. **The value is only directionally differentiable.** In place of a gradient one has
   $$
   V'(\theta;d)=\min_{\lambda\in\mathcal D^\star(\theta)}\big\langle\partial_\theta\mathcal L(z^\star,\lambda,\theta),\,d\big\rangle,
   $$
   which is generally **not linear in $d$** (it's a min over a set), so $V$ has
   **kinks** — exactly the $\max(\theta_1,\theta_2)$ / $H_\infty$‑norm nonsmoothness.
   $\partial_\theta\mathcal L$ at any returned $\lambda$ is then a **subgradient**, not
   the gradient.

3. **The solution map can be ill‑posed.** Degeneracy usually coincides with a
   non‑unique or non‑Lipschitz primal (A2 also fails), so the sensitivity $\partial
   z^\star/\partial\theta$ need not even exist — the outer optimizer's secant
   information becomes unreliable and BFGS curvature degrades.

In short: **degeneracy turns the exact‑gradient guarantee into a subgradient at best,
and can make the very quantity you differentiate nonsmooth.** That is why A3 is the
delicate assumption and why nondegeneracy is required to keep the framework exact.

---

## 6. Subgradients: what they are, when you get one, and why they help

### What a subgradient is

For a convex function $V:\mathbb R^p\to\mathbb R$, a vector $g$ is a **subgradient** at
$\theta$ if it defines a globally under‑estimating affine minorant:

$$
V(\theta')\ \ge\ V(\theta)+\langle g,\ \theta'-\theta\rangle\qquad\forall\,\theta'.
$$

The set of all such $g$ is the **subdifferential** $\partial V(\theta)$. Geometrically,
each $g$ is the slope of a hyperplane that touches the graph of $V$ at $\theta$ and lies
below it everywhere. Where $V$ is differentiable, $\partial V(\theta)=\{\nabla
V(\theta)\}$ — a single element, the ordinary gradient. Where $V$ has a **kink**,
$\partial V(\theta)$ is a whole set (an interval/polytope of admissible slopes).

The connection to the directional derivatives of §5 is exact: the subdifferential is
the set whose support function is the directional derivative,

$$
V'(\theta;d)=\max_{g\in\partial V(\theta)}\langle g,d\rangle
\qquad(\text{convex case}),
$$

so "$V$ has a kink" and "$\partial V$ is set‑valued" are the same statement. (For a
nonconvex but locally Lipschitz $V$ — the general co‑design value — the analogue is the
**Clarke subdifferential**, the convex hull of limits of nearby gradients; the same
picture holds locally.)

### When you get a subgradient instead of a gradient

Precisely when $V$ is **not differentiable** at $\theta$ — i.e. exactly the degeneracy
of §3–§5:

- **non‑unique dual** (primal degeneracy / LICQ failure): the envelope map returns
  $\partial_\theta\mathcal L(z^\star,\lambda,\theta)$ for *some* $\lambda\in\mathcal
  D^\star$, and each choice is a subgradient —
  $$
  \partial V(\theta)=\big\{\partial_\theta\mathcal L(z^\star,\lambda,\theta):\ \lambda\in\mathcal D^\star(\theta)\big\};
  $$
- **loss of strict complementarity** (rank overlap): a kink along the offending
  direction;
- **active‑set ties** (example (ii)), the **minimal‑$\gamma$ $H_\infty$** optimum
  (example (iii)), and nonsmooth regularizers used *without* the epigraph lift.

When A1–A3 hold (nondegenerate, strictly complementary), $\mathcal D^\star$ is a
singleton, so $\partial V(\theta)=\{\nabla V(\theta)\}$: **the subgradient *is* the
gradient.** Thus "subgradient vs. gradient" is exactly "degenerate vs. nondegenerate."

### Why a subgradient is still useful

- **Optimization still works.** Projected/proximal **subgradient** and **bundle**
  methods converge on nonsmooth convex problems using *any* $g\in\partial V$; the outer
  quasi‑Newton loop degrades gracefully to such a step. You lose the superlinear rate,
  not correctness.
- **Optimality certificate.** $\theta^\star$ is stationary iff
  $0\in\partial V(\theta^\star)$ (plus the box's normal cone) — a clean test that
  remains valid at kinks, where $\nabla V$ does not exist.
- **A guaranteed descent direction exists.** The steepest‑descent direction is $-g^\star$
  with $g^\star=\arg\min_{g\in\partial V(\theta)}\lVert g\rVert$ (the minimum‑norm
  subgradient); it points downhill whenever $0\notin\partial V$.
- **It is what nonsmooth $H_\infty$ synthesis already uses.** The minimal‑$\gamma$ value
  is nonsmooth, and subgradient/bundle methods (HIFOO, Apkarian–Noll) are the standard,
  effective tools there — so the fallback is a well‑trodden path, not a dead end.

So a subgradient is a principled continuation with a weaker (linear rather than
superlinear) rate; it keeps the method usable even when a design sits on a degenerate
face.

---

## 7. How to detect and ensure nondegeneracy

- **Rank test.** Compute the primal slack $M(z^\star)$ and its multiplier $S^\star$;
  strict complementarity ⇔ $\operatorname{rank}M^\star+\operatorname{rank}S^\star=\dim$.
  A rank overlap flags degeneracy.
- **Remove redundancy.** Drop dominated polytope vertices, duplicate caps, and
  linearly dependent structural equalities before solving — these are the commonest
  sources of dependent active constraints.
- **Add a strictly convex regularizer** (recipe 1: Tikhonov / log‑det). It pins the
  primal (A2) and, generically, perturbs the optimum off degenerate faces, restoring
  strict complementarity (A3). This is the cheapest universal safeguard.
- **Keep caps strictly feasible.** A variance/output cap set *above* the achievable
  floor is inactive or nondegenerately active; a cap pushed *to* the floor is the
  matrix analogue of a tie (example (ii)/(iii)) and should be avoided.
- **The finite‑difference gradient check** is a practical degeneracy detector: if the
  analytic envelope gradient matches central differences (as in the studies,
  rel. err $\sim10^{-6}$), the point is nondegenerate; a mismatch localizes a
  degenerate design.

---

## 8. The graceful fallback

Where degeneracy is unavoidable at isolated designs, the method does **not** fail: the
returned multiplier gives a valid **subgradient** (§6), and the projected/quasi‑Newton
outer loop degrades to a subgradient step (Remark on the degenerate case). Recipe 1
(a strictly convex regularizer) removes almost all such points, so in practice the
exact‑gradient regime holds everywhere the studies operate.

---

### References

- F. Alizadeh, J.-P. Haeberly, M. Overton, *Complementarity and nondegeneracy in
  semidefinite programming*, Math. Programming 77 (1997).
- A. Shapiro, *First and second order analysis of nonlinear semidefinite programs*
  / differentiability of optimal values (1997).
- J.F. Bonnans, A. Shapiro, *Perturbation Analysis of Optimization Problems*,
  Springer (2000).
- S. Boyd, L. El Ghaoui, E. Feron, V. Balakrishnan, *Linear Matrix Inequalities in
  System and Control Theory*, SIAM (1994).
