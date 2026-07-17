# Why the SDP dual is a rank‑$n$ lifting of the Riccati $P$

Context: for the H₂ **control** covariance SDP, the $2n\times2n$ Lagrange multiplier $T$
of the Lyapunov LMI is not the $n\times n$ Riccati matrix $P$ — it is a rank‑$n$
*lifting* of $P$. This note explains the mechanism. It rests on two facts:

1. at the optimum the primal LMI collapses onto a rank‑$n$ subspace, and
2. complementary slackness forces the dual onto the orthogonal complement of that
   subspace.

---

## Fact 1 — the active primal constraint is a rank‑$n$ outer product

The Lyapunov LMI is

$$
M_2=\begin{bmatrix}\Sigma-W & A\Sigma+BL\\[2pt] \star & \Sigma\end{bmatrix}\succeq 0 .
$$

At the optimum, minimizing $\operatorname{tr}(Q\Sigma)$ with $Q\succ0$ drives it to the
**equality**

$$
\Sigma^\star = A_{\mathrm{cl}}\,\Sigma^\star A_{\mathrm{cl}}^\top + W,
\qquad
A\Sigma^\star + BL^\star = A_{\mathrm{cl}}\,\Sigma^\star,
\qquad A_{\mathrm{cl}} = A+BK^\star .
$$

Substituting,

$$
M_2^\star
=\begin{bmatrix}A_{\mathrm{cl}}\Sigma^\star A_{\mathrm{cl}}^\top & A_{\mathrm{cl}}\Sigma^\star\\[2pt]
\Sigma^\star A_{\mathrm{cl}}^\top & \Sigma^\star\end{bmatrix}
=\underbrace{\begin{bmatrix}A_{\mathrm{cl}}\\ I\end{bmatrix}}_{G}\ \Sigma^\star\
\begin{bmatrix}A_{\mathrm{cl}}\\ I\end{bmatrix}^{\!\top}.
$$

Here $G$ is $2n\times n$ of full column rank and $\Sigma^\star\succ0$, so
$M_2^\star = G\,\Sigma^\star G^\top$ has **rank $n$** (not $2n$): the constraint sits on
a face of the PSD cone. Its range is the **graph of $A_{\mathrm{cl}}$**,

$$
\operatorname{range}(M_2^\star)=\operatorname{range}(G)
=\operatorname{range}\!\begin{bmatrix}A_{\mathrm{cl}}\\ I\end{bmatrix},
$$

and its kernel is

$$
\ker M_2^\star=\ker G^\top
=\Big\{[v_1;v_2]:A_{\mathrm{cl}}^\top v_1+v_2=0\Big\}
=\operatorname{range}\!\begin{bmatrix}I\\ -A_{\mathrm{cl}}^\top\end{bmatrix}.
$$

(Because $M_2^\star v=G\Sigma^\star G^\top v=0 \iff G^\top v=0$, since $G$ is injective
and $\Sigma^\star\succ0$.)

---

## Fact 2 — complementary slackness pushes the dual onto that kernel

For the PSD cone, complementary slackness is $\langle T,M_2^\star\rangle=0$ with
$T,M_2^\star\succeq0$. Two PSD matrices with zero inner product satisfy
$T\,M_2^\star=0$, i.e.

$$
\operatorname{range}(T)\subseteq\ker(M_2^\star)
=\operatorname{range}\!\begin{bmatrix}I\\ -A_{\mathrm{cl}}^\top\end{bmatrix}.
$$

A PSD matrix whose range lies in a subspace $\operatorname{range}(J)$ **must** factor as
$T=J\,\Xi\,J^\top$ for some $\Xi\succeq0$ on that subspace. With
$J=\begin{bmatrix}I\\ -A_{\mathrm{cl}}^\top\end{bmatrix}$ ($2n\times n$),

$$
\boxed{\,T=\begin{bmatrix}I\\ -A_{\mathrm{cl}}^\top\end{bmatrix}\Xi
\begin{bmatrix}I\\ -A_{\mathrm{cl}}^\top\end{bmatrix}^{\!\top},\qquad \Xi\in\mathbb S^n_+ .\,}
$$

**That is the lifting.** The $2n\times2n$ dual is not free — it has only the
$n(n{+}1)/2$ degrees of freedom of a single $n\times n$ matrix $\Xi$, sitting in the
$n$‑dimensional kernel coordinates.

---

## The geometry (why it is clean)

$$
G^\top J=[A_{\mathrm{cl}}^\top\ \ I]\begin{bmatrix}I\\ -A_{\mathrm{cl}}^\top\end{bmatrix}
=A_{\mathrm{cl}}^\top-A_{\mathrm{cl}}^\top=0,
$$

so the **primal range (graph of $A_{\mathrm{cl}}$) and the dual range
(graph of $-A_{\mathrm{cl}}^\top$) are orthogonal complements** in $\mathbb R^{2n}$.
That is complementary slackness made visible: the primal solution lives on one
$n$‑dimensional subspace, the dual on the perpendicular one. Neither can be full rank;
each carries exactly one $n\times n$ matrix — $\Sigma^\star$ for the primal, $\Xi$ for
the dual.

---

## Why the lifted core equals $P$

Here is the full derivation. It has two multipliers, both lifted, and the crucial move
is that the $L$‑stationarity condition cancels the $B$‑cross terms in the
$\Sigma$‑stationarity condition, leaving a pure Stein equation.

### Both multipliers are lifted

The same complementary‑slackness argument (Facts 1–2) applied to the *other* constraint
$M_1=\begin{bmatrix}Z_0&L\\ L^\top&\Sigma\end{bmatrix}$ gives, at the optimum,
$M_1^\star=\begin{bmatrix}K^\star\\ I\end{bmatrix}\Sigma^\star\begin{bmatrix}K^\star\\ I\end{bmatrix}^{\!\top}$
(rank $n$, kernel $\operatorname{range}\begin{bmatrix}I\\ -K^{\star\top}\end{bmatrix}$),
so its multiplier $S$ is lifted through an $m\times m$ core $\Theta\succeq0$. Writing both
multipliers in block form:

$$
S=\begin{bmatrix}I\\ -K^{\star\top}\end{bmatrix}\Theta\begin{bmatrix}I\\ -K^{\star\top}\end{bmatrix}^{\!\top}
=\begin{bmatrix}\Theta & -\Theta K^\star\\ -K^{\star\top}\Theta & K^{\star\top}\Theta K^\star\end{bmatrix},
\qquad
T=\begin{bmatrix}I\\ -A_{\mathrm{cl}}^\top\end{bmatrix}\Xi\begin{bmatrix}I\\ -A_{\mathrm{cl}}^\top\end{bmatrix}^{\!\top}
=\begin{bmatrix}\Xi & -\Xi A_{\mathrm{cl}}\\ -A_{\mathrm{cl}}^\top\Xi & A_{\mathrm{cl}}^\top\Xi A_{\mathrm{cl}}\end{bmatrix}.
$$

So the block reads: $S_{11}=\Theta,\ S_{12}=-\Theta K^\star,\ S_{22}=K^{\star\top}\Theta K^\star$
and $T_{11}=\Xi,\ T_{12}=-\Xi A_{\mathrm{cl}},\ T_{22}=A_{\mathrm{cl}}^\top\Xi A_{\mathrm{cl}}$.

### The Lagrangian and its three stationarity conditions

With $S,T\succeq0$ the multipliers of $M_1,M_2$,

$$
\mathcal L=\operatorname{tr}(Q\Sigma)+\operatorname{tr}(RZ_0)-\langle S,M_1\rangle-\langle T,M_2\rangle,
$$

and, expanding the inner products,

$$
\langle S,M_1\rangle=\operatorname{tr}(S_{11}Z_0)+2\langle S_{12},L\rangle+\operatorname{tr}(S_{22}\Sigma),
$$
$$
\langle T,M_2\rangle=\operatorname{tr}\!\big(T_{11}(\Sigma-W)\big)+2\operatorname{tr}\!\big(T_{12}^\top(A\Sigma+BL)\big)+\operatorname{tr}(T_{22}\Sigma).
$$

Setting the derivative of $\mathcal L$ in each primal variable to zero:

$$
\text{(I)}\quad \partial_{Z_0}:\ R-S_{11}=0,
$$
$$
\text{(II)}\quad \partial_{L}:\ S_{12}+B^\top T_{12}=0,
$$
$$
\text{(III)}\quad \partial_{\Sigma}:\ Q-S_{22}-T_{11}-T_{22}-\big(A^\top T_{12}+T_{12}^\top A\big)=0.
$$

(The $\partial_\Sigma$ term $A^\top T_{12}+T_{12}^\top A$ is the symmetrization of
$2\operatorname{tr}(T_{12}^\top A\Sigma)$; note it uses the **open‑loop** $A$, because the
constraint block is $A\Sigma+BL$.)

### Substitute the lifted blocks

**From (I):** $\Theta=R$. Hence immediately $S_{12}=-RK^\star$ and $S_{22}=K^{\star\top}RK^\star$.

**From (II):** $S_{12}=-B^\top T_{12}=-B^\top(-\Xi A_{\mathrm{cl}})=B^\top\Xi A_{\mathrm{cl}}$.
Combined with $S_{12}=-RK^\star$,

$$
\boxed{\,B^\top\Xi A_{\mathrm{cl}}=-RK^\star.\,}\tag{$\star$}
$$

This is the key relation. (It is consistent: for $\Xi=P$ one checks
$B^\top P A_{\mathrm{cl}}=R(R+B^\top PB)^{-1}B^\top PA=-RK^\star$.)

**From (III):** insert $S_{22}=K^{\star\top}RK^\star$, $T_{11}=\Xi$,
$T_{22}=A_{\mathrm{cl}}^\top\Xi A_{\mathrm{cl}}$, $T_{12}=-\Xi A_{\mathrm{cl}}$:

$$
K^{\star\top}RK^\star=Q-\Xi-A_{\mathrm{cl}}^\top\Xi A_{\mathrm{cl}}
+A^\top\Xi A_{\mathrm{cl}}+A_{\mathrm{cl}}^\top\Xi A.\tag{$\star\star$}
$$

### The cross terms cancel — this is where ($\star$) enters

The two terms $A^\top\Xi A_{\mathrm{cl}}+A_{\mathrm{cl}}^\top\Xi A$ still carry the
open‑loop $A$. Write $A=A_{\mathrm{cl}}-BK^\star$:

$$
A^\top\Xi A_{\mathrm{cl}}=A_{\mathrm{cl}}^\top\Xi A_{\mathrm{cl}}-K^{\star\top}\underbrace{B^\top\Xi A_{\mathrm{cl}}}_{=-RK^\star},\qquad
A_{\mathrm{cl}}^\top\Xi A=A_{\mathrm{cl}}^\top\Xi A_{\mathrm{cl}}-\underbrace{A_{\mathrm{cl}}^\top\Xi B}_{=(-RK^\star)^\top}K^\star,
$$

where both underbraces use ($\star$) (and $A_{\mathrm{cl}}^\top\Xi B=(B^\top\Xi A_{\mathrm{cl}})^\top=(-RK^\star)^\top=-K^{\star\top}R$). Hence

$$
A^\top\Xi A_{\mathrm{cl}}+A_{\mathrm{cl}}^\top\Xi A
=2A_{\mathrm{cl}}^\top\Xi A_{\mathrm{cl}}+K^{\star\top}RK^\star+K^{\star\top}RK^\star
=2A_{\mathrm{cl}}^\top\Xi A_{\mathrm{cl}}+2K^{\star\top}RK^\star.
$$

The $B$‑dependence is gone. Substituting this into ($\star\star$):

$$
K^{\star\top}RK^\star=Q-\Xi-A_{\mathrm{cl}}^\top\Xi A_{\mathrm{cl}}
+2A_{\mathrm{cl}}^\top\Xi A_{\mathrm{cl}}+2K^{\star\top}RK^\star
=Q-\Xi+A_{\mathrm{cl}}^\top\Xi A_{\mathrm{cl}}+2K^{\star\top}RK^\star.
$$

Cancel $2K^{\star\top}RK^\star$ from both sides and rearrange for $\Xi$:

$$
\boxed{\,\Xi=A_{\mathrm{cl}}^\top\,\Xi\,A_{\mathrm{cl}}+\big(Q+K^{\star\top}RK^\star\big).\,}
$$

### Conclusion

This is the discrete **Stein (Lyapunov) equation**. Because $A_{\mathrm{cl}}$ is Schur it
has a **unique** solution, and that solution is exactly the Lyapunov form of the control
Riccati matrix, $P=A_{\mathrm{cl}}^\top PA_{\mathrm{cl}}+Q+K^{\star\top}RK^\star$. Hence
$\Xi=P$, the core is unique, and the full dual is

$$
T=\begin{bmatrix}I\\ -A_{\mathrm{cl}}^\top\end{bmatrix}P
\begin{bmatrix}I\\ -A_{\mathrm{cl}}^\top\end{bmatrix}^{\!\top}.
$$

The interpretation matches the algebra: $\Xi=P$ is the **shadow price of the Lyapunov
(dynamics) constraint** — the sensitivity of the optimal cost to perturbing it — which
is precisely the cost‑to‑go / value matrix of the LQR problem.

---

## One‑line summary

The optimal covariance forces the constraint onto the graph of $A_{\mathrm{cl}}$
(rank $n$); complementary slackness forces the multiplier onto the orthogonal graph of
$-A_{\mathrm{cl}}^\top$ (also rank $n$); so the $2n\times2n$ dual is the rank‑$n$
lifting of a single $n\times n$ matrix, and stationarity makes that matrix the Riccati
$P$.
