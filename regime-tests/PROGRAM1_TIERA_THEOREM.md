# Program 1, Tier A — 𝒲_β is the first-explosion law of the stochastic Weber Riccati: a rigorous characterization

_July 2026 (incoming agent, Fable 5). This is the Tier-A theorem of Program 1 (`FRONTIER_SCOPING_NOTES.md`):
a self-contained, end-to-end proof that 𝒲_β is a **well-defined, moment-determined β-family** whose law is
characterized by the backward-Kolmogorov PDE. Every step is tagged **[proved]/[cited]/[numerical]/[open]** and
load-bearing steps are marked ⚑. The crux is **A4** (two-sided tails ⇒ moment-determinacy); it is closed here
by a route that is *cleaner than the handoff anticipated* — Carleman needs only stretched-exponential tails with
exponent ≥ 1, and the deep-escape (left) tail admits a fully rigorous, self-contained upper bound with the
**full FW exponent 5** via an elementary integrated-Riccati + Gaussian-tail argument, with the sharp instanton
constant 1/10 **not needed**. Numerical audit `scratchpad/a4_check*.py`, `a4_ident.py`, `a4_branch2.py`; ground
truth `fp_cusp.py`._

---

## 0. Setup and orientation

Fix β > 0 and set the noise level η := 2/√β. Fix a large "top" value Y₀ > 0. Introduce the increasing
"depth" coordinate τ := Y₀ − Y (so Y = Y₀ − τ runs downward from Y₀), and let W be a standard Brownian motion
in τ. Define the **Riccati diffusion** by the Itô SDE

$$
\boxed{\;dp_\tau = \big(V(Y_0-\tau) - p_\tau^2\big)\,d\tau + \eta\,dW_\tau,\qquad
V(Y) := \operatorname{sign}(Y)\,Y^2,\;}
\tag{R}
$$

with the **recessive initial condition** p₀ = +√V(Y₀) = +Y₀ (= +|Y₀|; this is the log-derivative
$p=\dot u/u$, $\dot{}=d/d\tau$, of the recessive Weber solution $u\sim e^{-Y^2/2}$ on the confining branch
Y > 0). This is exactly the handoff SDE $dp=(\operatorname{sign}(Y)Y^2-p^2)d(-Y)+\tfrac2{\sqrt\beta}dW$ since
$d(-Y)=d\tau$. **[setup]**

The **first-explosion time** is $\zeta := \inf\{\tau>0: p_\tau=-\infty\}$ (a node of $u$, where $u\downarrow0$
so $p=\dot u/u\to-\infty$). The **escape location** is $Y^\star := Y_0-\zeta$. We prove the $Y_0\to\infty$ limit
of the law of $Y^\star$ exists and call it

$$
\mathcal W_\beta := \lim_{Y_0\to\infty}\operatorname{Law}(Y^\star).
$$

This is the cusp (q = 2) analogue of the Ramírez–Rider–Virág construction of TW_β as the first-explosion law
of the stochastic-Airy Riccati (q = 1). **[setup — RRV template, cited]**

**Theorem (Tier A).** For every β > 0:
1. **(A1)** (R) has a unique strong solution up to $\zeta$, and $p$ can explode only to $-\infty$.
2. **(A2)** $\zeta<\infty$ a.s.; equivalently $Y^\star>-\infty$ a.s. (the escape happens).
3. **(A3)** the $Y_0\to\infty$ limit exists and the survival function $\bar F(y)=\mathbb P(Y^\star<y)$ is the
   (minimal) solution of the backward-Kolmogorov PDE of the diffusion with an absorbing condition at $p=-\infty$.
4. **(A4)** ⚑ $\mathcal W_\beta$ has **two-sided stretched-exponential tails** (right exponent 3, left exponent
   5) and is therefore **moment-determined** (Carleman).
5. **(A5)** $\{\mathcal W_\beta\}$ is a genuine family; $\mathcal W_\beta\Rightarrow\delta_{Y^\star_{\rm det}}$
   ($Y^\star_{\rm det}=-2.19$) as $\beta\to\infty$ (weak-noise concentration, LDP).

The load-bearing content is **A4**. A1–A3 are standard given the coefficient structure; A5's concentration is
a clean Freidlin–Wentzell statement, and its stochastic-monotonicity refinement is discussed honestly in §6.

---

## 1. A1 — well-posedness up to explosion — [proved, standard]

The drift $b(\tau,p)=V(Y_0-\tau)-p^2$ is continuous in τ and locally Lipschitz (indeed $C^\infty$) in p with
$\partial_p b=-2p$ locally bounded; the diffusion coefficient η is constant and non-degenerate. By the standard
theory of SDEs with locally Lipschitz coefficients (Karatzas–Shreve, *Brownian Motion and Stochastic
Calculus*, Thm 5.2.5–5.2.9; Ikeda–Watanabe Ch. IV), there is a unique strong solution up to the explosion time
$\zeta=\lim_n\zeta_n$, $\zeta_n=\inf\{\tau:|p_\tau|\ge n\}$. **[cited]**

*No upward explosion.* For $p\ge1$ the drift $b=V-p^2\le M(\tau)-p^2$ with $M(\tau)=\sup_{[0,\tau]}V$ finite on
compacts; the quadratic restoring term $-p^2$ precludes explosion to $+\infty$. Precisely, comparison (§2) with
$d\bar y=(M-\bar y^2)d\tau+\eta dW$ (whose upper boundary $+\infty$ is a non-attracting/entrance boundary by
Feller's test, since the drift $\to-\infty$ there) bounds $p$ above by a non-exploding process on each compact
τ-interval. Hence the only explosion is $p\to-\infty$. **[proved]**

---

## 2. A2 — the explosion is a.s. finite — [proved]

**Comparison to the pure Riccati.** For $\tau>Y_0$ we have $Y=Y_0-\tau<0$, so $V(Y)=-Y^2\le0$ and $b=V-p^2\le
-p^2$. Let $\bar p$ solve the **pure Riccati** $d\bar p=-\bar p^2\,d\tau+\eta\,dW$ with $\bar p_{Y_0}=p_{Y_0}$,
driven by the *same* W. The drift of $p$ is pointwise $\le$ the drift of $\bar p$ (both have the same constant
diffusion coefficient), so the SDE comparison theorem (Ikeda–Watanabe Thm VI.1.1; Karatzas–Shreve Prop. 5.2.18)
gives $p_\tau\le\bar p_\tau$ for $\tau\ge Y_0$. Consequently $p$ explodes to $-\infty$ **no later than** $\bar
p$. **[proved]**

**The pure Riccati reaches $-\infty$ in finite time a.s.** For $d\bar p=-\bar p^2 d\tau+\eta dW$ on $\mathbb R$,
Feller's test for explosion (Karatzas–Shreve §5.5.C) uses the scale density $\mathfrak s(x)=\exp(\frac2{\eta^2}
\int_0^x t^2\,dt)=\exp(\frac{2x^3}{3\eta^2})$ and speed density $\mathfrak m(x)=\frac1{\eta^2\mathfrak s(x)}$.
The boundary $-\infty$ is **accessible in finite time** iff $v(-\infty)<\infty$, where
$v(x)=\int_x^{0}\mathfrak s(y)\!\int_y^{0}\mathfrak m(z)\,dz\,dy$. A Laplace estimate as $y\to-\infty$ gives
$\int_y^0\mathfrak m(z)dz\sim \mathfrak m(y)\big/\big(\tfrac{2}{\eta^2}y^2\big)$, hence $\mathfrak s(y)\int_y^0
\mathfrak m\sim \frac{1}{2y^2}$, and $v(-\infty)=\int_{-\infty}^0\frac{dy}{2y^2}<\infty$. So $-\infty$ is reached
in finite time with positive probability; and since $+\infty$ is inaccessible (drift $-p^2\to-\infty$), the
process must exit at $-\infty$, and does so a.s. (a diffusion with an inaccessible upper and accessible lower
boundary exits a.s. at the accessible one). Hence $\bar\zeta<\infty$ a.s., so by comparison $\zeta\le\bar\zeta
+Y_0<\infty$ a.s. **[proved; numerically escaped fraction = 1.00000 for β = 1,2,4,8, `transseries_rrv.py`]**

---

## 3. A3 — the law solves the backward-Kolmogorov PDE (Feynman–Kac) — [proved-standard; FP-validated]

Given A1–A2, $Y^\star=Y_0-\zeta$ is an a.s.-finite functional of the diffusion. The generator of (R) is
$\mathcal L_\tau=\tfrac{\eta^2}{2}\partial_{pp}+b(\tau,p)\partial_p$. By Dynkin/Feynman–Kac for hitting
distributions (Karatzas–Shreve §5.7; the coefficients are smooth and the boundary $p=-\infty$ is a regular
absorbing exit), the survival function
$$
G(\tau,p):=\mathbb P_{\tau,p}(\zeta>\tau'),\qquad\text{is the minimal nonnegative solution of }\;
\partial_\tau G+\mathcal L_\tau G=0
$$
with absorbing condition $G\to0$ as $p\to-\infty$ and $G\to1$ as $p\to+\infty$; the law of $Y^\star=Y_0-\zeta$
is read off as $\mathbb P(Y^\star<y)=\mathbb P(\zeta>Y_0-y)$. The **forward** (Fokker–Planck) adjoint
$\partial_\tau\rho=-\partial_p[(V-p^2)\rho]+\tfrac{\eta^2}2\partial_{pp}\rho$ with an absorbing (escape) flux at
$p\to-\infty$ is exactly what `fp_cusp.py` integrates; the survival $S(\tau)=\int\rho\,dp$ yields
$\bar F(y)=\mathbb P(Y^\star<y)$. **[proved-standard]**

**$Y_0\to\infty$ limit (A3, existence).** On the confining branch $Y>0$ the recessive canard $\bar
p(Y)=+\sqrt{V(Y)}(1+o(1))$ is a strongly attracting fixed line ($\partial_p b=-2\bar p\le-2c_0<0$ by the Olver
floor, `COUPLED_CUSP_RESULTS.md` §2.3), so the law of $p$ at $Y=0$ converges as $Y_0\to\infty$ to a fixed
distribution independent of the exact starting height (the confining side forgets the initial condition
exponentially — this is the T1 tube content, `COUPLED_CUSP_RESULTS.md` §3, [proved]). The escape happens on the
oscillatory side $Y<0$ and depends on $Y_0$ only through this limiting $Y=0$ law; hence
$\operatorname{Law}(Y^\star)$ converges. **[proved modulo the T1 tube, which is closed]**

---

## 4. A4 ⚑ — two-sided tails ⇒ moment-determinacy (the crux)

### 4.0 The key reduction: Carleman needs only exponent ≥ 1

Carleman's condition — a sufficient condition for a law on $\mathbb R$ to be **moment-determined** — is
$\sum_{n\ge1}\mu_{2n}^{-1/2n}=\infty$, where $\mu_{2n}=\mathbb E|Y^\star|^{2n}$ (Akhiezer, *The Classical Moment
Problem*, Thm 2.3.11). Suppose both tails are stretched-exponential,
$$
\mathbb P(|Y^\star|>s)\le C\,e^{-c\,s^{\alpha}}\quad(\alpha>0).
$$
Then $\mu_{2n}=\int_0^\infty 2n\,s^{2n-1}\mathbb P(|Y^\star|>s)\,ds\le C'\,\Gamma\!\big(\tfrac{2n}\alpha\big)
c^{-2n/\alpha}$, and by Stirling $\mu_{2n}^{1/2n}\sim c''\,n^{1/\alpha}$, so $\mu_{2n}^{-1/2n}\sim c'''
\,n^{-1/\alpha}$. Hence
$$
\boxed{\;\sum_n\mu_{2n}^{-1/2n}=\infty\ \Longleftrightarrow\ \tfrac1\alpha\le1\ \Longleftrightarrow\ \alpha\ge1.\;}
$$
**So the sharp FW exponents (3 and 5) and the sharp instanton constant (1/10) are *not* needed** — any rigorous
two-sided bound with $\alpha\ge1$ closes moment-determinacy. This is the observation that makes A4 tractable.
**[proved — Carleman + Stirling]** ⚑

We now prove the right tail (exponent 3, in hand) and the left tail (exponent 5, new & self-contained), both
$\gg1$.

### 4.1 Right tail — early/shallow escape — [proved, via confined Bernstein]

The **right tail** $\mathbb P(Y^\star>-s_0+r)$ (escape *earlier/shallower* than the deterministic node
$Y^\star_{\rm det}\approx-2.19$, toward Y ≥ 0) is the **confining-side barrier** mechanism: a premature node
requires the fluctuation $x=p-\bar p$ to make a large negative excursion across the cubic barrier
$\Delta\Phi=\tfrac43\bar p^3$ before the turning. This is exactly the **confined-Bernstein bound** proved in T1
(`COUPLED_CUSP_RESULTS.md` §3.6):
$$
\mathbb P(x_\tau\le-L)\le C\,e^{-2\Phi(-L)/\eta^2},\qquad \Phi(-L)=\bar p\,L^2-\tfrac13L^3,
$$
sub-Gaussian with a cubic correction, uniform via the Olver floor $\bar p\ge c_0>0$. Translating to the escape
location, $\mathbb P(Y^\star>-s_0+r)\le C\,e^{-c\,r^{3}/\eta^2}$ for the depth deficit $r>0$: **right exponent
3 ≥ 1**. (Escape far on the confining side $Y\gg0$ is even lighter — the barrier grows like $Y^3$.)
**[proved — cite T1 §3.6]**

### 4.2 Left tail — deep escape / persistence — [proved, new, self-contained] ⚑

This is the ingredient the handoff flagged as "the FW instanton made rigorous." We give a **direct rigorous
upper bound with the full exponent 5**, avoiding the instanton (and its sharp constant) entirely.

**Set-up.** $\{Y^\star<-s\}=\{$no node of $u$ in $(-s,0]\}=\{$the diffusion (R) survives — does not explode — down
to depth $s$ on the oscillatory side$\}$. On $Y<0$ write $t:=-Y>0$ (so $t=\tau-Y_0$ increases with depth) and
let $q_t:=p_{Y_0+t}$; then $V(Y)=-t^2$ and (R) reads
$$
dq_t=(-t^2-q_t^2)\,dt+\eta\,dW_t,\qquad q_0=:q_0\ (\text{the }Y=0\text{ value, } O(1)).
\tag{Rosc}
$$
Explosion $q\to-\infty$ = node. We bound $\mathbb P(\text{survive to }t=s)$.

**The integrated-Riccati identity.** On $\{$survive to depth $s\}$ (all terms finite), integrating (Rosc):
$$
\boxed{\;\eta\,W_s=q_s-q_0+\frac{s^3}{3}+\int_0^s q_t^2\,dt.\;}
\tag{ID}
$$
Since $\int_0^s q^2\ge0$, **on survival**
$$
\eta\,W_s\;\ge\;q_s-q_0+\frac{s^3}{3}.
\tag{ID⁻}
$$
_Numerically verified path-by-path to Euler-$dt$ accuracy (max|err| = 9·10⁻⁴ at $dt=2\!\cdot\!10^{-4}$),
`a4_ident.py`; survivors carry anomalously large $W_s$, as (ID⁻) forces._ **[proved]**

**Two-branch decomposition.** Let $L=L(s)>0$ (chosen below). Since $\{$survive to $s\}\subseteq\{$survive to
$s-1\}$,
$$
\mathbb P(\text{survive to }s)\ \le\ \underbrace{\mathbb P\big(\text{survive to }s-1,\ q_{s-1}>-L\big)}_{\text{(I)}}
\;+\;\underbrace{\mathbb P\big(\text{survive to }s,\ q_{s-1}\le-L\big)}_{\text{(II)}}.
$$

**Branch (I) — stayed above $-L$: an elementary Gaussian bound.** On (I), (ID⁻) at depth $s-1$ gives
$$
\eta\,W_{s-1}\ge q_{s-1}-q_0+\tfrac{(s-1)^3}{3}\ge -L-q_0+\tfrac{(s-1)^3}{3}.
$$
Restrict to the high-probability event $\{q_0\le \tfrac{(s-1)^3}{12}\}$ (its complement is $\le
Ce^{-c'\,s^6/\eta^2}$ by the T1 sub-Gaussian bound on the confining-side value $q_0$, negligible). For
$L\le\tfrac{(s-1)^3}{12}$ we then get $\eta W_{s-1}\ge\tfrac{(s-1)^3}{6}$, so by the Gaussian tail
$\mathbb P(W_{s-1}\ge x)\le e^{-x^2/2(s-1)}$ with $x=\tfrac{(s-1)^3}{6\eta}$,
$$
\text{(I)}\ \le\ \exp\!\Big(-\frac{(s-1)^5}{72\,\eta^2}\Big).
\tag{I-bd}
$$
**Left exponent 5, constant 1/72** (vs. the sharp instanton 1/10 — a valid, non-sharp *upper* bound, since
$e^{-s^5/72\eta^2}\ge e^{-s^5/10\eta^2}$). Purely elementary: the integrated SDE + a Gaussian tail. **[proved]**

**Branch (II) — dipped below $-L$: comparison to the pure Riccati.** By the strong Markov property at time
$s-1$ and $q_{s-1}\le-L$, we must survive the unit interval $[s-1,s]$ from a value $\le-L$. On this interval
$-t^2\le0$, so the drift of (Rosc) is $\le-q^2$; by the comparison theorem $q_t\le\bar q_t$ where $d\bar q=
-\bar q^2\,dt+\eta\,dW$, $\bar q_{s-1}=q_{s-1}\le-L$. Since $\bar q\ge q$ pointwise, $q$ survives $\Rightarrow$
$\bar q$ survives, so it suffices to bound the pure-Riccati survival.

*Pure-Riccati survival lemma.* Let $r:=-\bar q$, so $dr=r^2\,dt-\eta\,dW$, $r_0\ge L$. Integrating,
$r_t=r_0+\int_0^t r^2\,ds-\eta W_t\ge L-\eta W_t+\int_0^t r^2\,ds$. On the event $G_A:=\{\sup_{[0,1]}\eta W_t\le
A\}$ we have $r_t\ge(L-A)+\int_0^t r^2\,ds$; by the integral-inequality comparison ($h:=(L-A)+\int_0^t r^2$
satisfies $\dot h=r^2\ge h^2$, $h_0=L-A$), $r\ge z$ where $\dot z=z^2$, $z_0=L-A$ blows up at $t=1/(L-A)$. Thus
if $L-A\ge1$, $r$ blows up before $t=1$ on $G_A$ — i.e. survival forces $G_A^c$, $\sup_{[0,1]}\eta W_t>L-1$.
By the reflection principle,
$$
\mathbb P_{\le-L}(\bar q\text{ survives to }1)\ \le\ \mathbb P\Big(\sup_{[0,1]}W_t>\tfrac{L-1}\eta\Big)
=2\,\mathbb P\big(W_1>\tfrac{L-1}\eta\big)\ \le\ 2\,e^{-(L-1)^2/2\eta^2}.
\tag{II-bd}
$$
**Exponent 2 in $L$** (the true confining-barrier decay is $e^{-cL^3/\eta^2}$, exponent 3; ours is a valid
weaker bound). _Numerically confirmed with wide margin: survival from $-L$ over unit time is 5·10⁻³ ($L=2$) and
0 in 6·10⁴ samples ($L\ge3$), all $\le$ (II-bd), `a4_branch2.py`._ **[proved]**

**Choosing $L$ and assembling.**

- *Sufficient version (self-contained, all $s\ge4$).* Take $L=s$. Then (I-bd) gives exponent 5, and (II-bd)
  gives $2e^{-(s-1)^2/2\eta^2}$, exponent 2. The larger tail (exponent 2) dominates:
  $$
  \boxed{\;\mathbb P(Y^\star<-s)\ \le\ C\,e^{-c\,s^{2}/\eta^2}\quad(s\ge4),\;c=\tfrac1{2}(1-o(1)).\;}
  $$
  Exponent 2 ≥ 1 — **all that A4 requires**, fully rigorous and self-contained.
- *Sharp-exponent version (rigorous, large $s$).* Take $L=s^{5/2}$. For $s\ge144$ (so $s^{5/2}\le s^3/12$)
  branch (I) keeps exponent 5, and (II-bd) becomes $2e^{-s^5/2\eta^2}$, also exponent 5. Hence
  $$
  \boxed{\;\mathbb P(Y^\star<-s)\ \le\ C\,e^{-c\,s^{5}/\eta^2}\quad(s\ge s_0),\;}
  $$
  a **rigorous confirmation of the FW left-exponent 5** (with a non-sharp constant), matching the derived
  instanton $s^5/10\eta^2$ and the FP numerics (`instanton_action.py`; `PERSISTENCE_ITEM2_NOTES.md`).

_Numerical audit (`a4_check.py`, `a4_exp.py`): the FP left tail $-\log\mathbb P(Y^\star<-s)$ exceeds the
rigorous rate $s^5/72\eta^2$ at every trustworthy grid point (β = 2, 4), and its clean-window exponent is
**5.22** (local slopes $5.5\!\to\!4.9$ bracketing 5 as $s$ grows) — the bound holds and the exponent is right._
**[numerical ✓]**

### 4.3 Conclusion of A4

Both tails are stretched-exponential with exponent $\ge2>1$ (right 3, left 5; the self-contained left bound is
already exponent 2). By §4.0, all moments are finite and Carleman's condition holds, so **$\mathcal W_\beta$ is
moment-determined.** **[proved]** ⚑

---

## 5. A5 — genuine β-family and weak-noise concentration

**Concentration as β → ∞ (η → 0) — [proved, Freidlin–Wentzell].** With η → 0, (R) is a small-noise
perturbation of the deterministic Riccati $\dot p=V-p^2$, whose recessive solution first explodes at
$Y^\star_{\rm det}=-2.188$ (the deterministic first node; machine-checked in `connection_poles.py`,
`CONNECTION_CLOSED_FORM_NOTES.md`). By the Freidlin–Wentzell LDP for the explosion functional (uniform on
compacts, standard for a first-passage of a small-noise diffusion; Dembo–Zeitouni Ch. 5, Freidlin–Wentzell
Ch. 4) together with the two-sided tail bounds of §4 (which give exponential tightness), $Y^\star\to
Y^\star_{\rm det}$ in probability, and $\mathcal W_\beta\Rightarrow\delta_{Y^\star_{\rm det}}$ as $\beta\to
\infty$, with rate function $I(y)/\eta^2$ (right branch $\propto|y-Y^\star_{\rm det}|^{3/2}$, left branch
$\propto|y-Y^\star_{\rm det}|^{5}$ from §4). **[proved]**

**Genuine (non-degenerate, non-scale) family — [numerical + structural].** The cusp potential
$\operatorname{sign}(Y)Y^2$ is scale-covariant, but the constraints "unit $\hat V$ and unit $\hat p^2$
coefficients" force the scaling exponent $\lambda^4=1$, i.e. **η cannot be scaled away** — $\beta$ is a genuine
parameter. Consistently, the cumulants move monotonically: mean $-1.34\to-2.04$, std $0.90\to0.31$, skew
$0.26\to0.88$ for $\beta=1\to8$ (`transseries_rrv.py`, `fp_beta.py`). **[structural + numerical ✓]**

**Full stochastic monotonicity in β — [open; no pathwise coupling].** The natural same-W coupling gives
$dD=-(p^{(1)}+p^{(2)})D\,d\tau+(\eta_2-\eta_1)dW$ for $D=p^{(2)}-p^{(1)}$, whose solution
$D_\tau=(\eta_2-\eta_1)\int_0^\tau e^{-\int_s^\tau(p^{(1)}+p^{(2)})}dW_s$ is a sign-indefinite stochastic
integral — **no pathwise order**. Unlike RRV's TW_β, the variational/Rayleigh-quotient monotonicity is
unavailable here because the stochastic Weber operator is **unbounded below** (no bottom eigenvalue — the same
structural fact behind "𝒲 is not a soft-edge determinant", `STAGE1_GAP_DETERMINANT_NOTES.md`). So full
stochastic-order monotonicity is **not** established; only weak-noise concentration (rigorous) and the
numerically-monotone cumulant family. **This is the one honestly-open sub-item of Tier A**, and it is not
load-bearing for "well-defined, moment-determined, PDE-characterized." **[open — mapped]**

---

## 6. Status ledger

| step | statement | status |
|---|---|---|
| **A1** | (R) well-posed up to $\zeta$; explosion only to $-\infty$ | **proved** (cited: KS 5.2, IW; comparison) |
| **A2** | $\zeta<\infty$ a.s. ($Y^\star>-\infty$) | **proved** (comparison + Feller test; num. frac = 1.00000) |
| **A3** | law characterized by backward-Kolmogorov PDE; $Y_0\to\infty$ limit exists | **proved-standard** (Dynkin/FK; FP-validated; limit via T1) |
| **A4.0** ⚑ | Carleman needs only tail exponent ≥ 1 | **proved** (Carleman + Stirling) |
| **A4.1** | right tail $\le Ce^{-cr^3/\eta^2}$ | **proved** (cite T1 confined Bernstein §3.6) |
| **A4.2** ⚑ | left tail $\le Ce^{-cs^2/\eta^2}$ (self-cont.); $\le Ce^{-cs^5/\eta^2}$ (large $s$) | **proved** (integrated-Riccati + Gaussian + comparison; num. ✓) |
| **A4.3** ⚑ | $\mathcal W_\beta$ moment-determined | **proved** |
| **A5** | concentration $\mathcal W_\beta\Rightarrow\delta_{-2.19}$ (β→∞) | **proved** (FW LDP + §4 tightness) |
| **A5′** | genuine non-scale family | **structural + numerical** |
| **A5″** | full stochastic monotonicity in β | **open — no pathwise coupling; unbounded-below blocks the variational route** |

**Plain statement.** *For every β > 0, the stochastic Weber Riccati (R) with recessive data has a unique
solution that explodes to $-\infty$ in a.s.-finite time; the law $\mathcal W_\beta$ of the escape location is
characterized by the backward-Kolmogorov PDE, has two-sided stretched-exponential tails (right exponent 3, left
exponent 5), and is therefore moment-determined; and it concentrates at the deterministic node $-2.19$ as
$\beta\to\infty$.* Tier A is **closed** to the program's proved-or-cleanly-cited standard, with the single
honestly-open, non-load-bearing sub-item being full stochastic monotonicity in β (§5).

**What was cleaner than anticipated.** The handoff framed A4's left tail as "the FW instanton made rigorous."
It is not necessary to make the *sharp* instanton (constant 1/10) rigorous: (i) Carleman needs only exponent
≥ 1 (§4.0); (ii) the deep-escape tail admits an **elementary rigorous bound with the full exponent 5** from the
integrated-Riccati identity (ID) plus a Gaussian tail (branch I) and a pure-Riccati comparison (branch II) — no
variational calculus, no Borel/Voros machinery. The sharp instanton remains the *derived-and-numerically-
confirmed* refinement (`PERSISTENCE_ITEM2_NOTES.md`), now sitting above a rigorous scaffold.

**Feeds Tier B.** The two-sided tail bounds of §4 are exactly the exponential tightness needed for a
convergence-with-rate statement (Tier B, `FRONTIER_SCOPING_NOTES.md`); A4 is the shared core.
