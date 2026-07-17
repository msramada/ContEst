# Concise A1–A4 verification for the fixed‑γ $H_\infty$ case

This mirrors the H₂ verification. The important scoping point first:

> **A1–A4 hold for the _fixed‑γ_ $H_\infty$ (guaranteed‑cost) problem, with $\gamma$ held
> strictly above the optimal attenuation $\gamma^\star(\theta)$.** They do **not** hold
> for the _minimal‑γ_ synthesis $\min\gamma^2$: there the bounded‑real LMI is rank
> deficient, strict complementarity (A3) fails, and the value $\gamma^{\star2}(\theta)$
> is only directionally differentiable (subgradient only).

So "the $H_\infty$ case" below means the fixed‑γ guaranteed cost, which is exactly what
the MTDC study and Proposition `prop:hinfreg` use.

---

## Setup

Linear plant with disturbance channel and performance output,

$$
x_{k+1}=A_\theta x_k+B_\theta u_k+E_\theta w_k,\qquad
z_k=C_z x_k+D_z u_k,\qquad E_\theta=W^{1/2},
$$

with $C_z=[Q^{1/2};0]$, $D_z=[0;R^{1/2}]$. Fix an attenuation level $\gamma$. The
worst‑case (guaranteed) control cost is $\operatorname{tr}(XW)$, where $X\succeq0$ is
the stabilizing solution of the discrete $H_\infty$ **game** Riccati equation — a DARE
in the **augmented data**

$$
\tilde B_\theta=[\,B_\theta\ \ E_\theta\,],\qquad
\tilde R=\operatorname{diag}(R,\,-\gamma^2 I),\qquad
A_{\mathrm{cl}}=A_\theta+\tilde B_\theta\tilde K,\ \
\tilde K=-(\tilde R+\tilde B_\theta^\top X\tilde B_\theta)^{-1}\tilde B_\theta^\top XA_\theta .
$$

Equivalently, as a **convex** program (the conic form the assumptions apply to), with
$Y=P^{-1}$, $L=KY$,

$$
V(\theta)=\min_{Y\succ0,\,L,\,Z}\ \operatorname{tr}(WZ)
\quad\text{s.t.}\quad
\underbrace{\begin{bmatrix}-Y&0&(A_\theta Y+B_\theta L)^\top&(C_zY+D_zL)^\top\\
0&-\gamma^2 I&E_\theta^\top&0\\
A_\theta Y+B_\theta L&E_\theta&-Y&0\\
C_zY+D_zL&0&0&-I\end{bmatrix}}_{\text{fixed‑}\gamma\text{ BRL}}\preceq0,\quad
\begin{bmatrix}Y&I\\ I&Z\end{bmatrix}\succeq0,
$$

where the Schur slack gives $Z\succeq Y^{-1}=P$, so minimizing $\operatorname{tr}(WZ)$
drives $P$ to the minimal storage $X$ and $V(\theta)=\operatorname{tr}(XW)$.

**Standing hypotheses.** $(A_\theta,B_\theta)$ stabilizable, $Q,R,W\succ0$, and
$\gamma>\gamma^\star(\theta)$ — i.e. $\gamma$ is a *feasible, non‑minimal* attenuation.
This last condition is the $H_\infty$ replacement for "$R\succ0$" in H₂: it guarantees
$R+B_\theta^\top XB_\theta\succ0$ (control‑block positivity) and $A_{\mathrm{cl}}$ Schur.

---

## The key fact: it is the H₂ proof with augmented data

The game Riccati is a DARE in $(A_\theta,\tilde B_\theta,Q,\tilde R)$ — the **only**
change from H₂ is $B_\theta\mapsto\tilde B_\theta$, $R\mapsto\tilde R$ (indefinite),
$P\mapsto X$. Every algebraic step of the H₂ appendix used the DARE structure and
$A_{\mathrm{cl}}$ Schur, **never** $R\succ0$; it used only the admissibility
$R+B_\theta^\top XB_\theta\succ0$, which $\gamma>\gamma^\star$ supplies. Hence the whole
A1–A4 argument transfers verbatim.

---

## A1 — convexity, smooth data, strong duality

$\operatorname{tr}(WZ)$ is linear; the BRL and the Schur block are LMIs (affine in
$(Y,L,Z)$) — a convex program, with data affine in
$(A_\theta,B_\theta,E_\theta,C_z,D_z,W)$, hence jointly $C^1$ in $(z,\theta)$. **Slater
holds because $\gamma>\gamma^\star$:** a strictly suboptimal certificate satisfies the
BRL *strictly* (that is the meaning of $\gamma$ being above the optimum), and
$Z\succ Y^{-1}$ is strict; so the interior is nonempty and strong duality with attained
primal/dual solutions follows. (At $\gamma=\gamma^\star$ this interior collapses — the
degenerate case.)

## A2 — unique primal minimizer

Minimizing $\operatorname{tr}(WZ)$ with $W\succ0$ drives $Z=Y^{-1}=P$ to the **minimal
storage** $P=X$, the unique stabilizing game‑Riccati solution. The game completion of
squares (valid because it needs only $R+B_\theta^\top XB_\theta\succ0$, not $R\succ0$)

$$
\operatorname{tr}(P_KW)-\operatorname{tr}(XW)
=\operatorname{tr}\!\big((K-K^\star)^\top(R+B_\theta^\top XB_\theta)(K-K^\star)\,\Sigma_K\big)\ \ge\ 0
$$

vanishes iff $K=K^\star$ (the central controller), so $(Y^\star,L^\star,Z^\star)=
(X^{-1},K^\star X^{-1},X)$ is unique.

## A4 — smooth parameter entry

The variables and cones are $\theta$‑independent; $\theta$ enters only through
$A_\theta,B_\theta,E_\theta$ (and, if co‑designed, $C_z,D_z,W$), each $C^1$ in $\theta$.
So $\nabla_\theta V=\langle\lambda^\star,\partial_\theta(\text{data})\rangle$ once
$\lambda^\star$ is unique.

## A3 — unique dual, strict complementarity

Two equivalent routes, both from $\gamma>\gamma^\star$:

- **Riccati/IFT (Proposition `prop:hinfreg`).** The stabilizing game‑Riccati solution
  $X$ is unique, and the GARE residual linearizes to the closed‑loop Stein operator
  $\delta X\mapsto\delta X-A_{\mathrm{cl}}^\top\delta X A_{\mathrm{cl}}$, invertible
  since $A_{\mathrm{cl}}$ is Schur. So $X(\theta)\in C^1$ and $V=\operatorname{tr}(XW)$
  is differentiable with a **unique** sensitivity — the A3 conclusion.
- **Dual lifting (H₂ appendix, augmented data).** Complementary slackness lifts the BRL
  multiplier onto the kernel of the active certificate, collapsing it to a core that the
  stationarity conditions pin to $X$ (the same Stein equation as H₂, now with
  $\tilde R,\tilde B$). Strict complementarity holds because $\gamma>\gamma^\star$ keeps
  $A_{\mathrm{cl}}$ Schur and $R+B_\theta^\top XB_\theta\succ0$, so the primal‑slack and
  dual ranks add to full — exactly the H₂ rank count with $P\mapsto X$.

## Conclusion

Under stabilizability, $Q,R,W\succ0$, and $\gamma>\gamma^\star(\theta)$, the fixed‑γ
$H_\infty$ guaranteed‑cost problem satisfies A1–A4, so $V(\theta)=\operatorname{tr}(XW)$
is differentiable and its gradient is read from the inner dual — equivalently, from the
game‑Riccati envelope identity

$$
\frac{\partial\,\operatorname{tr}(XW)}{\partial A_\theta}=2\,X A_{\mathrm{cl}} S,
\qquad S=A_{\mathrm{cl}} S A_{\mathrm{cl}}^\top+W,
$$

chained with $\partial_\theta A_\theta$ (and the analogous $\tilde B_\theta,W$ terms).
This is the H₂ formula with $(B,R,P)\mapsto(\tilde B,\tilde R,X)$, and as
$\gamma\to\infty$ it collapses back to the H₂ result $\operatorname{tr}(PW)$.

**Estimation half.** The dual filter game GARE is regular under the analogous conditions
(swap $A_\theta\!\leftrightarrow A_\theta^\top$, $B_\theta\!\leftrightarrow C_\theta^\top$,
etc.), giving $V_{\mathrm{est}}(\theta)=\operatorname{tr}(Q\Sigma_e)\in C^1$.

---

## The boundary case (why "fixed‑γ" is essential)

At $\gamma=\gamma^\star(\theta)$ (minimal‑γ synthesis): the worst‑case gain is *attained*,
the Hamiltonian has imaginary‑axis eigenvalues, the BRL certificate is rank deficient,
and $\operatorname{rank}(\text{primal slack})+\operatorname{rank}(\text{dual})<\dim$ —
**strict complementarity fails**. The optimal face is not a vertex, the multiplier is
non‑unique, and $\gamma^{\star2}(\theta)$ is nonsmooth. There A2/A3 break and one only
has a subgradient. Fixing $\gamma$ strictly above $\gamma^\star$ is precisely what moves
the problem into the regular A1–A4 regime.
