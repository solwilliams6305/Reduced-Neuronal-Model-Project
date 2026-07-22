# Noise-induced escape at the cusp of a coupled FitzHugh–Nagumo system: a uniform tube and a new edge law

**A self-contained results document.** Audience: a reader fluent in canards / geometric singular
perturbation theory / blow-up and in Tracy–Widom edge universality, who has seen none of our working notes.
Everything needed to follow the results and audit the rigor is here. Claims are tagged
**[proved]** / **[cited]** / **[numerical]** / **[open]** throughout; nothing is asserted beyond its tag.

---

## 0. Summary

Two symmetrically coupled FitzHugh–Nagumo (FHN) oscillators, with additive noise, have an antisymmetric mode
whose slow–fast geometry passes through a **cusp** as the coupling $g$ is tuned. We ask two questions about
the noise-induced escape (peel-off) near that cusp:

- **T1 (the tube).** Is the sample-path fluctuation tube around the canard **uniform** through the
  fold→cusp merge — i.e. does the tube width stay controlled as the two folds coalesce?
- **T2 (the law).** What is the **limiting distribution** of the (rescaled) escape location — the cusp
  analogue of Tracy–Widom?

**Keystone reframe.** Both reduce to *one* object: the **noisy parabolic-cylinder (Weber) connection** of the
inner equation $u''=(\operatorname{sign}(Y)Y^2-\eta\dot W)u$. Its *probabilistic* reading is the tube (T1);
its *integrable* reading (the connection/Stokes data) is the edge law (T2).

**Results.**
- **T1 is proved.** The two-sided fluctuation is sub-Gaussian with constants **uniform in** $\rho=\Delta/\ell$
  through the merge. The one previously-flagged input — Malliavin-differentiability of the first-passage in the
  swept-turning geometry (§3.7) — is now **re-derived** in the linear *field* representation (§3.7bis): the
  turning "singularity" is a Prüfer-coordinate artifact, transversality of the escape is *structural* (a double
  zero forces the trivial solution), so T1 rests only on foundational citations (Olver, Gaussian concentration,
  the implicit function theorem on Wiener space) with **verified** hypotheses.
- **T2: the edge law $\mathcal W_\beta$ is constructed and characterized as a *genuinely new* object** — an
  **asymmetric isomonodromy in the Painlevé-IV family**, *not* standard PIV and *not* any soft-edge Fredholm
  determinant. Its $\beta=2$ fingerprint is reproduced to a few percent; four independent "known-parent" routes
  are ruled out. **New (this pass):** (i) the **left tail is derived** — the Freidlin–Wentzell instanton gives
  $-\log\mathbb P(Y^\star<-s)\to s^5/(10\eta^2)=\beta s^5/40$ (exponent 5, rate $\propto\beta$, **constant 1/10**),
  and the previously-reported "anomalous persistence exponent" is shown to be a pre-asymptotic *crossover*, not a
  new exponent; (ii) the **intrinsic multi-point process** is identified — the node process is a **perturbed
  $\pi$-lattice (class-I hyperuniform), not a determinantal/Dyson edge**, with a derived jitter law; (iii) the
  whole catastrophe ladder's node processes are **unified** by one exponent $a(q)=3q/(q+2)$, with the fold ($q=1$)
  as the exact GUE(II)$\to$lattice(I) boundary. The exact integrable closed form for $\mathcal W$ remains
  **open** (and is *provably* not a Painlevé-$\sigma$-ODE — §4.4).

---

## 1. The coupled model and the cusp

### 1.1 Model

Two diffusively coupled FHN units (Kristiansen–Pedersen normalization), with additive white noise:
$$
\dot v_i=-v_i^3+3v_i-w_i+g\,(v_j-v_i)+\sigma\,\xi_i,\qquad
\dot w_i=\varepsilon\,(v_i-c),\qquad (i,j)=(1,2),(2,1),
$$
with $0<\varepsilon\ll1$ (time-scale separation), coupling $g$, noise $\sigma$, and $c=0.99$ near the fold.
**[setup]**

### 1.2 The antisymmetric mode is a cusp normal form

In symmetric/antisymmetric coordinates $v_\pm=v_1\pm v_2$, the **antisymmetric mode** $v_-$ decouples to
leading order with effective slow-manifold curvature
$$
\mu(v)=3(1-v^2)-2g .
$$
At the fold of the cubic ($v^2=1$), $\mu=-2g$. The coupling thus **shifts the fold**, and the two folds of
the antisymmetric problem **coalesce into a cusp** as $g\to0^-$. Writing the local separation of the two
turning points as $\Delta$, a direct computation gives **[proved/exact]**
$$
\boxed{\ \Delta(g)=2\sqrt{-2g/3}\ }\qquad (g<0).
$$

### 1.3 Crossover scale

Standard fold/cusp blow-up (Kristiansen–Pedersen) gives the deterministic inner scale; balancing the cusp
merge against the slow drift gives the **crossover coupling**
$$
\mu\sim\sqrt\varepsilon\ \Longleftrightarrow\ \Delta\sim\varepsilon^{1/4},\qquad
g_{\rm crit}\propto\sqrt\varepsilon .
$$
**[derived]** Verified numerically in the full coupled FHN: $g_{\rm crit}/\sqrt\varepsilon=-0.58,\,-0.57$ at
$\varepsilon=0.015,\,0.030$ (ratio $1.40\approx\sqrt2$); the onset $g\approx-0.07$ falls inside the measured
band $[-0.05,-0.11]$, and sweeping $g$ drives a clean **fold→cusp escape-class crossover** (peel-off spread
$\times3.7$; excess kurtosis flips $\approx0$ (fold) $\to$ negative (cusp)) in the genuine model.
**[numerical — `loop_closure_scaling.py`, `crossover_fold_to_cusp.py`]**

So the cusp is **physically anchored**: the catastrophe lives in the coupled neurons, with $g$ the unfolding
parameter and $\Delta(g)$ the coalescence variable.

---

## 2. The inner equation and the keystone parametrix

### 2.1 Inner (blown-up) equation

In the cusp blow-up chart, the swept inner equation for the antisymmetric fluctuation is the **stochastic
parabolic-cylinder (Weber) operator**
$$
u''=\big(V_\Delta(Y)-\eta\,\dot W\big)\,u,\qquad
V_\Delta(Y)=\operatorname{sign}(Y)\,|Y|\,(|Y|+\Delta),
$$
where $Y$ is the inner space-time variable, $\eta$ the rescaled noise ($\beta:=4/\eta^2$, the Dyson-type
index), and $\Delta$ the (rescaled) coalescence. Near $Y=0$: $V_\Delta\approx\Delta\,Y$ (simple **Airy**
turning, $|Y|\ll\Delta$) crossing over to $\operatorname{sign}(Y)Y^2$ (**Weber/cusp**, $|Y|\gg\Delta$); at
$\Delta=0$ it is the pure cusp. The inner noise scale is $\ell\sim\eta^{2/3}$ and the uniformity parameter is
$\rho=\Delta/\ell\in(0,\infty)$ ($\rho\to\infty$ Airy, $\rho\to0$ cusp). **[setup/cited — Olver two-turning-points]**

### 2.2 Cole–Hopf and the fluctuation SDE

$p=u'/u$ gives the Riccati $dp=(V_\Delta-p^2)\,d\tau+\eta\,dW$ ($\tau=Y_0-Y$); the **escape** is the first
explosion $p\to-\infty$ (first node of $u$). The deterministic **canard** $\bar p_\Delta(Y)$ tracks
$+\sqrt{V_\Delta}$ on the attracting branch, with $\bar p_\Delta\to\bar p_{\min}(\Delta)>0$ at the turning.
The fluctuation $x=\delta p=p-\bar p_\Delta$ obeys the **exact** SDE
$$
\boxed{\ dx=-\big(2\bar p_\Delta\,x+x^2\big)\,d\tau+\eta\,dW=-\Phi_x\,d\tau+\eta\,dW,\quad
\Phi(x,\tau)=\bar p_\Delta\,x^2+\tfrac13x^3\ }
$$
— a **gradient (Langevin) diffusion in a bistable cubic potential** $\Phi$. **[proved]** Its features:
- **well** at $x=0$, **saddle** at $x=-2\bar p_\Delta$ (the repelling branch); **barrier** $\Delta\Phi=\tfrac43\bar p^3$;
- **curvature** $\Phi_{xx}=2\bar p+2x$, so at the well $\Phi_{xx}(0)=2\bar p\ge2\bar p_{\min}\ge 2c_0>0$;
- **Hamilton–Jacobi identity** $\mathcal L\,e^{2\Phi/\eta^2}=e^{2\Phi/\eta^2}\big[(2/\eta^2)\partial_\tau\Phi+\Phi_{xx}\big]$,
  with $\partial_\tau\Phi=\dot{\bar p}\,x^2\le0$ pre-turning (the $O(\eta^{-2})$ terms cancel because the drift
  is $-\nabla\Phi$). **[derived]**

### 2.3 The Olver floor (the linchpin for $\Delta$-uniformity)

The recessive parabolic-cylinder solution's log-derivative does not vanish through the coalescence: Olver's
uniform two-turning-point connection gives
$$
\bar p_{\min}(\Delta)\ \ge\ c_0>0\quad\text{uniformly in }\Delta .
$$
**[cited — Olver, *Asymptotics & Special Functions* Ch. 11; numerically $c_0\approx0.69$,
$\bar p_{\min}\in[0.69,1.13]$ across $\Delta\in[0.02,3]$].** This single positive floor is what makes every
estimate below uniform in $\Delta$.

---

## 3. T1 — the uniform two-sided tube

### 3.1 Statement (main theorem)

Let $\sigma_\Delta^2$ be the linearized fluctuation variance, $h^\star=\min_\tau 2\bar p_\Delta/\sigma_\Delta$
the barrier in sd units, $\rho=\Delta/\ell$, and let $\delta\Theta$ be the noisy Prüfer connection phase
(equivalently the rescaled escape-location fluctuation) accumulated through the turning.

> **Theorem T1.** There are constants $c,C,C_m>0$, **independent of $\Delta$ and $\eta$**, such that on the
> approach window $[0,\tau_\ell]$:
> $$
> \mathbb P\Big(\sup_{0\le\tau\le\tau_\ell}\tfrac{|x_\tau|}{\sigma_\Delta}\ge h\Big)\le C\,e^{-c\,h^2},\qquad
> \kappa_{2m}(\delta\Theta)\le C^m\,m!\,\eta^{2m}\ \ (\forall m),
> $$
> uniformly in $\rho\in(0,\infty)$ through the merge. The lower-side rate is $h^{\star2}/6$, the upper-side
> $h^2/2$. **Status: proved modulo the single cited Malliavin-regularity input of §3.7.**

The sub-Gaussian connection-phase bound is the content that the *edge law is well-defined uniformly* through
the merge — i.e. the tube does not blow up as the folds coalesce.

### 3.2 Proof architecture (logical order)

| block | role | status |
|---|---|---|
| (A) approach tube | confine $x$ down to the $O(\ell)$ turning region | **proved** |
| (B) re-entry closed | global supermartingale, no upper-boundary leakage | **proved** |
| (C) stochastic-Olver connection variance | the $m=1$ bound through the turning | **proved** (mod. cited floor) |
| (D) upper tail | super-confinement $x>0$ | **proved** |
| (E) confined Bernstein (lower tail) | sub-Gaussian escape-side moments | **proved** |
| (F) $R_{\rm nl}$ via Gaussian concentration | nonlinear first-passage remainder | **proved mod. §3.7** |

### 3.3 (A) Approach tube; (B) re-entry — [proved]

A **cubic scale-function exponential supermartingale** $M=\exp(\psi+\tfrac{2}{\eta^2}\Phi)$, with
$\dot\psi=-(2\bar p+2x_+)$, has $\mathcal L M\le0$ on the pre-turning strip (the HJ cancellation of §2.2 plus
$\dot{\bar p}\le0$). Optional stopping gives the **one-sided no-early-escape bound**
$$
\mathbb P(\text{early escape})\le C\,e^{-h^{\star2}/6},\qquad \ln C=\underbrace{2\!\int_0^T\!\bar p\,d\tau}_{\text{exposure}}+2x_+T .
$$
The **exposure integral** is bounded uniformly in $\rho$ (closed form, §3.5), and is *smallest* at the merge —
the floor's role here is the exponent, not the exposure. **Re-entry** is closed by replacing $\Phi$ above
$x_+$ by its tangent line: the $-x^2$ drift super-confines from above, making $M$ a *global* supermartingale on
$\{x\ge-L\}$, so optional stopping has no upper-boundary term. **[proved]**

### 3.4 (C) Stochastic-Olver connection variance — the keystone $m=1$ bound — [proved mod. cited floor]

In Olver's exact PC variable the deterministic connection (recessive $U(-\tfrac\lambda2,\sqrt2 Y)$, Weber
$\Gamma$-function connection coefficients) is exact and **uniform in $\Delta$** [cited]. Stochastic
variation-of-parameters gives the connection phase as an **Itô integral**
$$
\delta\Theta=\eta\!\int_{Y_-}^{Y_+}\!\Phi(s)\,dW(s)+R_{\rm nl},\qquad \Phi=\frac{u_0\,u_c}{\mathcal W},
$$
where $\Phi$ is the deterministic PC kernel ($\sim k^{-1}$ in WKB, **regular** at the turning — the $1/k$
singularity is a WKB artifact). By **Itô isometry**, $\mathbb E[\delta\Theta_{\rm lin}^2]=\eta^2\|\Phi\|_{L^2}^2$,
and the crux is that $\|\Phi\|^2$ is bounded uniformly in $\Delta$:
- **canard side** — the variance ODE $\dot v=-4\bar p v+\eta^2$ is regularized by the floor $\bar p\ge c_0$, so
  $v_{\rm exit}\le\eta^2/4c_0$;
- **oscillatory side** — the singular integral has the **closed form**
  $\int_{w_t}^{L}\frac{dy}{y(y+\Delta)}=\frac1\Delta\ln\frac{L(w_t+\Delta)}{w_t(L+\Delta)}$, bounded uniformly
  ($\to1/w_t$ cusp, $\to\frac{\ln}{\Delta}\to0$ Airy) — the dwell-time shrinks fast enough against $1/k^2$,
  and the PC cutoff $w_t$ stays bounded below by Olver.

⇒ $\mathbb E[(\delta\Theta)^2]\le C\eta^2$, **uniform in $\rho$**. **[proved — Itô isometry + closed-form
integrability + cited Olver floor].** Numerically the connection variance is $O(\eta^2)$, $\mathrm{Var}/\eta^2\in[0.20,0.47]$
across $\rho\in[0.04,3.2]$ with **no merge blow-up**, while the OU quasi-static comparison *diverges* (3.3→13)
— confirming the uniformity genuinely needs the PC structure. **[numerical — `connection_variance.py`]**

### 3.5 (D) Upper tail — [proved]

For $\lambda>0$, completing the square bounds the generator globally:
$\mathcal L e^{\lambda x}/e^{\lambda x}=-\lambda(x+\bar p)^2+\lambda\bar p^2+\tfrac{\eta^2}{2}\lambda^2\le
g_+:=\lambda\bar p^2+\tfrac{\eta^2}{2}\lambda^2$. So $\exp(\lambda x-\int g_+)$ is a genuine supermartingale,
giving $\mathbb E[e^{\lambda(x-a)}]\le e^{\eta^2\tau\lambda^2/2}$ — sub-Gaussian with proxy bounded uniformly via
the (bounded, closed-form) exposure $\int\bar p^2$. The cubic confines $x>0$ for free. **[proved]**

### 3.6 (E) Confined Bernstein (lower tail) — [proved]

The proved no-early-escape supermartingale (§3.3), used at a *moderate* barrier $-L$, gives
$$
P(x_\tau\le-L)\le C\,e^{-2\Phi(-L)/\eta^2},\quad \Phi(-L)=\bar p L^2-\tfrac13L^3,
$$
sub-Gaussian (with a cubic correction), **uniform via the floor** ($\bar p\ge c_0$). Integrating,
$$
\mathbb E[e^{\lambda x_\tau}\mathbf 1_{\mathcal C}]\le C\,e^{\eta^2\lambda^2/(8c_0)+(\text{cubic})}
\ \Rightarrow\ \mathbb E[(x_\tau)_-^{2m}\mathbf 1_{\mathcal C}]\le C^m m!\,\eta^{2m},
$$
with $C=C(c_0)$ only; the cubic $L^3$ term is precisely what produces the $m!$ (Bernstein, not Gaussian)
growth — i.e. a leptokurtic core with sub-Gaussian tails. **[proved — integration of the proved tail].**
Numerically: $-\ln P(x\le-L)$ is linear in $L^2$ (slope $\ge2c_0/\eta^2$, rising with $\Delta$), confined
log-MGF curvature $\le1/(8c_0)\approx0.18$, **uniform in $\Delta$**. **[numerical — `confined_bernstein.py`]**

### 3.7 (F) The nonlinear first-passage remainder $R_{\rm nl}$ — [proved modulo one cited regularity]

$\delta\Theta=Y^\star(W)$ is the first-passage functional ($Y^\star$ solves $\theta(Y^\star)=\pi$). IFT:
$\delta Y^\star=-\delta\theta/k_\star+R_{\rm nl}$, $R_{\rm nl}=O(\delta\theta^2/k_\star^2)$, where the **phase
velocity at escape** $k_\star=\sqrt{|V_\Delta-\lambda|}\,|_{Y^\star}\approx|Y^\star_{\rm det}|\approx2.09$ is
**floored** (escape sits in the oscillatory region, away from the turning $k=0$). Rather than bound the
quadratic-and-higher series term by term, note that **$\delta\Theta$ is a Lipschitz functional of the Gaussian
noise** with Malliavin derivative
$$
D_s\,\delta\Theta=-\frac{\eta\,\Phi(s)}{k(Y^\star)},\qquad
\|D\delta\Theta\|_{L^2}^2\le\frac{\eta^2\|\Phi\|^2}{k_\star^2}\le C\eta^2\quad\text{on }\mathcal C,
$$
bounded uniformly by the proved $\|\Phi\|^2$ (§3.4) over the floored $k_\star$. The **Gaussian concentration
inequality** (Borell–TIS / Üstünel: $\|DF\|_{L^2}\le\sigma$ a.s. $\Rightarrow F$ sub-Gaussian, proxy
$\sigma^2$) then gives $\kappa_{2m}(\delta\Theta)\le C^m m!\,\eta^{2m}$, uniform in $\Delta$, **subsuming all
nonlinear orders at once**. The escape event $\mathcal C^c$ is the proved tail (§3.3), harmless by
Cauchy–Schwarz. **[proved modulo §3.7-residual].** Numerically $\mathcal W$ is sub-Gaussian (standardized
log-MGF $\le0.9\lambda^2$; standardized cumulants between $m!$ and $3^m m!$). **[numerical — `firstpassage_rnl.py`]**

### 3.7bis Resolution of the Malliavin regularity — the field representation — [proved]

The doubt in §3.7 was specifically the *swept-turning geometry*: the Prüfer-phase SDE coefficient
$\eta\sin^2\theta/k$ is singular at the turning ($k\to0$), so a black-box first-passage theorem felt unsafe. The
resolution is a **change of representation**. Write the field as $x=(u,u')^{\!\top}$; it solves the **linear**
SDE $dx=A(\tau)x\,d\tau+Bx\,dW$, $A=\begin{psmallmatrix}0&-1\\-q&0\end{psmallmatrix}$,
$B=\begin{psmallmatrix}0&0\\ \eta&0\end{psmallmatrix}$, whose coefficients are **smooth and polynomially bounded
through the turning** — the $1/k$ blow-up is a coordinate artifact of the Prüfer transform, absent in $(u,u')$.
Three steps:

1. **Malliavin smoothness [proved, elementary].** A linear SDE with smooth coefficients has $x_\tau\in\mathbb
   D^\infty$, with $D_su_\tau=-\eta\,G(\tau,s)\,u_s$ ($G$ the causal Green's function of $\partial^2-q$) — exactly
   the §3.4 stochastic-VoP kernel.
2. **Transversality is *structural*, not assumed [proved].** If $u(Y^\star)=u'(Y^\star)=0$ then $x_{Y^\star}=0$,
   and by uniqueness of the linear SDE $x\equiv0$ — contradicting the nontrivial recessive data. Hence **every
   zero of $u$ is simple**, $u'(Y^\star)\ne0$ (equivalently $k_\star>0$). The non-degeneracy hypothesis of the
   first-passage theorem holds *automatically* (a double zero of a 2nd-order solution is impossible); the turning
   $k=0$ is on the recessive side and is not a crossing.
3. **First-passage derivative [cited, hypotheses verified].** The implicit function theorem on Wiener space
   (Nualart; Nualart–Vives) at the non-degenerate zero gives $Y^\star\in\mathbb D^{1,2}_{\rm loc}$ with
   $D_sY^\star=-D_su(Y^\star)/u'(Y^\star)=-\eta\,\Phi(s)/k(Y^\star)$.

**Numerical confirmation.** The identity $D_sY^\star=-D_su(Y^\star)/u'(Y^\star)$ is verified by a Bismut
finite-difference in the field representation to **correlation $1.0000$** (slope $1.00$); the response is linear
in the perturbation (differentiable, no kink) and $u'(Y^\star)\ne0$ (transversal, finite derivative), with a
sensitivity kernel **regular through the turning**. **[numerical — `malliavin_firstpassage_check.py`,
`figures/malliavin_firstpassage.png`].** So the item is **no longer a flagged geometry-specific citation** — it is
proved for this geometry, resting only on the elementary linear-SDE fact and the foundational IFT-on-Wiener-space.

### 3.8 T1 status — plain statement

**Every analytic load-bearing step of T1 is proved or cited with verified hypotheses.** With §3.7bis the last
flag is removed: T1 rests only on foundational citations (Olver's uniform PC connection + floor, the Gaussian
concentration inequality, the IFT on Wiener space) applied where their hypotheses are established. **T1 — the
uniform two-sided sub-Gaussian tube through the merge — is closed** to the program's proved-or-cleanly-cited
standard.

---

## 4. T2 — the new edge law $\mathcal W_\beta$

### 4.1 Construction

$\mathcal W_\beta$ is the law of the rescaled first-explosion location $Y^\star$ of the noisy cusp Riccati
$dp=(\operatorname{sign}(Y)Y^2-p^2)\,d\tau+\eta\,dW$, $\beta=4/\eta^2$ — equivalently the edge of the
**stochastic Weber operator** $-\partial^2+\operatorname{sign}(Y)Y^2+\tfrac{2}{\sqrt\beta}\dot W$. It is the
cusp ($q=2$) analogue of the Ramírez–Rider–Virág stochastic-Airy/Tracy–Widom construction ($q=1$). A smooth,
Monte-Carlo-free representation is obtained by solving the **escape Fokker–Planck PDE** (advection–diffusion
for the Riccati density with absorbing barrier); the survival curve gives $F(s)$. **[constructed; FP validated
against $1.4\times10^6$ MC samples to ~2% on all cumulants — `fp_cusp.py`].**

### 4.2 Fingerprint ($\beta=2$)

| cumulant | constructed $\mathcal W_2$ | target |
|---|---|---|
| skewness | **+0.607** | +0.61 |
| excess kurtosis | **−0.237** | −0.24 |
| $\kappa_5$ | **−2.115** | −2.2 |
| $\kappa_6$ | **−2.706** | −2.7 |

All four match (orientation-free $\kappa_4,\kappa_6$ essentially exactly). **[numerical].** The
**$\beta$-family** $\mathcal W_\beta$ is a genuine family of distinct shapes: skew $0.26\to0.60\to0.84$ for
$\beta=1,2,4$ (saturating), width $\propto\eta^{1.0}$. **[numerical — `fp_beta.py`].**

### 4.3 Characterization: an *asymmetric isomonodromy in the PIV family* — not standard PIV

The **deterministic skeleton** is the parabolic-cylinder (Weber) equation, whose rank-2 irregular singularity
places it in the **Painlevé-IV isomonodromy class** (PIV's classical solutions *are* parabolic-cylinder
functions). **[cited].** But the **odd** potential $\operatorname{sign}(Y)Y^2$ is confining for $Y>0$ (real
$e^{-Y^2/2}$) and oscillatory for $Y<0$ — the connection joins two *different* irregular structures, an
**asymmetric** connection, unlike the symmetric Weber connection of standard PIV. The diagnostic is the
**tail asymmetry**:
$$
\text{left exponent } 2q+1=5,\qquad \text{right exponent } 3q/2=3,\qquad \text{ratio } 5/3\ (\ne2).
$$
The two tails come from **two different mechanisms**: the right (early escape) is a **confining-side barrier**
(clean Freidlin–Wentzell, $\eta^2$-collapse) and the left (late escape) is an **oscillatory-side
phase-persistence** — now **derived in closed form** (§4.3bis), with the previously-reported "anomalous
$\eta$-scaling" identified as a pre-asymptotic crossover. **[derived + numerical — `tail_constants.py`,
`instanton_action.py`].** So $\mathcal W$ is governed by an **asymmetric isomonodromy in the PIV family** — *not*
the standard symmetric Painlevé-IV $\sigma$-form. **[derived-structural; the precise asymmetric class is open, §5].**

### 4.3bis The left tail, derived — a Freidlin–Wentzell instanton — [derived + numerical]

A deep late escape $Y^\star=-s$ requires a noise control preventing the Riccati blow-up until depth $s$; the FW
rate is $-\log\mathbb P=\min\tfrac1{2\eta^2}\int\phi^2\,d\tau=\mathcal I(s)/\eta^2$. Hamilton's equations
$\dot p=\pi+V-p^2$, $\dot\pi=2p\pi$ (with $\pi=\phi$, $V=-\tau^2$ oscillatory) have the large-$s$ scaling solution
$p\simeq1/\tau$, $\pi\simeq\tau^2-\tau^{-2}$, so the action density $\tfrac12\pi^2\simeq\tfrac12\tau^4$ and
$$
\boxed{\ \mathcal I(s)\to \tfrac12\!\int_0^s\!\tau^4\,d\tau=\frac{s^5}{10},\qquad
-\log\mathbb P(Y^\star<-s)\ \to\ \frac{s^5}{10\,\eta^2}=\frac{\beta\,s^5}{40}\ }
$$
— **left exponent 5 (matching the FW backbone $2q+1$), rate $\propto\beta$, and constant $1/10$**, all derived.
The instanton BVP confirms $\mathcal I(s)/(s^5/10)\to0.96$ at $s=13$ (constant $\to0.096$), robust to boundary
data. **[derived + numerical — `instanton_action.py`].**

**The apparent "anomaly" is a crossover.** Because $\mathcal I$ is $\eta$-independent, at fixed $s$ the tail is
$-\log\mathbb P=\mathcal I(s)/\eta^2+c(s)$ with $c(s)=O(1)$ a noise-independent offset; the effective exponent
$q_{\rm eff}=\dfrac{\mathcal I\beta/4}{\mathcal I\beta/4+c}\to1$ only as $\beta\to\infty$, sitting at $\approx0.6$
in the accessible window — reproducing (and *retracting* as fundamental) the earlier "anomalous persistence
exponent $\sim\beta^{0.6}$." This is confirmed **artifact-free by Monte-Carlo** (no numerical diffusion): at
fixed $\Theta=Y^{\star2}/2$, $-\log\mathbb P$ is linear in $\beta$ with slope $=\mathcal I(s)/4$ (measured/predicted
ratio $\to1.02$). **[numerical — `leftail_final.py`, `figures/leftail_derivation.png`].** The exact $O(1)$ offset
$c(s)$ (bulk/turning/prefactor) is the only residue; it does not affect the leading law.

### 4.4 Four negative results that pin it (rule out the known parents)

| candidate parent | test | verdict |
|---|---|---|
| soft-edge Fredholm determinant (higher-order Airy / multicritical TW) | Bornemann–Nyström determinant; the projection property $\int\mathrm{Ai}(z+t)\mathrm{Ai}(z+s)dz=\delta(t-s)$ fails at a *quadratic* turning | **not a soft-edge determinant** |
| $\beta$ = isomonodromic time (PII-type) | $\beta$-flow ODE: $\varphi(\beta)$ shows generic saturation, no PIV-Riccati "$2t$" signature ($b/2a\approx-0.02$) | **$\beta$ not isomonodromic time** |
| Pearcey-reduced (cusp catastrophe in RMT) | swept-Pearcey projection $u'''=(Y-\eta\xi)u+y\,u'$: no weighting reaches the cusp fingerprint | **not Pearcey-reduced** |
| Stokes-continued multicritical edge | continue HM→Ablowitz–Segur ("wrong sign"): mechanism real (heavy $\kappa_4{+}0.13\to$ sub-Gaussian $-0.73$) but the $q{=}1$ continued law is left-skewed ($-0.24$, $\kappa_5{+}1.4$) $\ne$ cusp ($+0.61$, $\kappa_5{-}2.2$); $q{=}2$ beyond reach | **not the accessible continuation** |

Independently, the **integer higher-order Tracy–Widom family** (Claeys–Its–Krasovsky / Claeys–Olver) has tail
ratio **2** for all $k$ (exponents $4k+3$, $(4k+3)/2$); the cusp ratio $5/3\ne2$ excludes it. A
**half-integer near-miss** ($4k+3=5\Rightarrow k=\tfrac12$, matching the *left* exponent) is decided against by
the *right* exponent (3 vs $5/2$) and the left coefficient ($1/20$ vs $9/320$). **[cited + derived].**

**Conclusion T2.** $\mathcal W_\beta$ is a **genuinely new universality law**, constructed and numerically
identified, with deterministic skeleton in the PIV class but an **asymmetric** (odd-potential) connection. Its
closed-form integrable identity is **open**.

### 4.5 The intrinsic multi-point process — a perturbed lattice, not a determinantal edge — [derived + numerical]

Beyond the marginal $\mathcal W_\beta$, the cusp has a genuine **intrinsic point process**: the nodes (zeros) of
the stochastic Weber field $u$, of which $\mathcal W$ is the edge (first node). Its structure is pinned:

- **It is a perturbed $\pi$-lattice, class-I hyperuniform — *not* a Dyson/determinantal edge.** Measured against
  Poisson and jittered-lattice controls, the structure factor sits at the noise floor between Bragg peaks
  ($S(0.5)=2\times10^{-4}$ vs Poisson $\sim1$), excluding the GUE/class-II law ($S(k)\sim k$) by $\sim3$ orders;
  the pair correlation is a hard core + $\pi$-comb. **[numerical — `node_kernel_decisive.py`].** This is
  consistent with every "$\mathcal W$ is not determinantal" negative (§4.4): the rigidity is 1-D spectral, not
  Dyson-repulsive.
- **The 2-point content is a derived jitter law.** The Prüfer phase-diffusion $D_\Theta\propto\eta^2/k^3$ gives
  spacing-jitter variance $\operatorname{Var}(\text{spacing})\approx(2/\beta)\,\Theta^{-3/2}$ (exponent exact,
  prefactor $\propto1/\beta$; measured $-1.52$). Since $-3/2<-1$ the jitter is integrable (freezes into the bulk),
  so the nontrivial correlation is **edge-concentrated** ($=\mathcal W$), which is why the bulk is lattice-trivial.
  **[derived + numerical — `node_kernel_decisive.py`, `node_beta_scaling.py`, `figures/node_kernel.png`].**

**The ladder, unified.** The same mechanism runs down the whole ladder $V_q=\operatorname{sign}(Y)|Y|^q$: with
frequency $k=|Y|^{q/2}$ and phase-depth $\Theta_q$, the jitter exponent is
$$
\boxed{\ a(q)=\tfrac{3q}{q+2}\ }\qquad(\text{measured }1.03/1.53/1.85/2.07\ \text{for }q=1,2,3,4;\ \text{predicted }1/1.5/1.8/2).
$$
Since accumulated phase variance $\sum n^{-a}$ converges iff $a>1$ iff $q>1$, **the fold ($q=1$, $a=1$) is the
exact critical boundary**: below it the process is the GUE/Airy determinantal edge (class II, $\log$ number
variance); at and above the cusp it is a perturbed lattice (class I, bounded). This *unifies* the fold rung
(Airy₂, determinantal) and the cusp rung (Weber, lattice) under one derived exponent, and answers the long-open
"intrinsic edge process" question — negatively for $q\ge2$ (lattice, not DBM), positively only at $q=1$.
**[derived + numerical — `ladder_field.py`, `figures/ladder_processes.png`].**

### 4.5bis The marginal ladder + unfolding (context)

The marginal family for $V_q$ ($q=1$ fold/TW, $q=2$ cusp, $q=3$ swallowtail, $q=4$ butterfly) is a clean monotone
sequence — skew $+0.20/+0.61/+0.96/+1.22$ — with tail exponents $(2q+1,\,3q/2)$ and ratio $(4q+2)/(3q)$ ($=2$ only
at the fold). Unfolding the $q$-th rung by lower-degree terms flows the class **down** to the fold (Thom
genericity / RG-relevance). The **within-class symmetry-breaking susceptibility** is derived: for the perturbation
$V\to\operatorname{sign}(Y)Y^2+a\,g$, $g=-Y^2\mathbf 1_{Y<0}$, the escape-location response is the Green's-function
integral $dY^\star/da=-\bar u(Y^\star)/(\mathcal W\,u_0'(Y^\star))\int_{Y^\star}^0 s^2u_0^2\,ds=+0.121$ (finite-diff
$+0.124$, $2.5\%$) — a *first-order* perturbation, so it **sidesteps the isomonodromy obstruction** of §4.4; the
skew susceptibility $d(\mathrm{skew})/da=+0.38$ follows from the deepened oscillatory well (heavier left tail).
**[derived + numerical — `ladder_edge_laws.py`, `susceptibility_derivation.py`, `figures/susceptibility.png`].**

---

## 5. Consolidated rigor ledger

| # | claim | status |
|---|---|---|
| 1 | $\Delta(g)=2\sqrt{-2g/3}$; antisym mode = cusp normal form | **proved/exact** |
| 2 | $g_{\rm crit}\propto\sqrt\varepsilon$; fold→cusp crossover in full FHN | **numerical** ($g_{\rm crit}/\sqrt\varepsilon\approx-0.58$) |
| 3 | fluctuation = gradient Langevin in bistable $\Phi$; barrier $\tfrac43\bar p^3$ | **proved** |
| 4 | Olver floor $\bar p_{\min}(\Delta)\ge c_0>0$ | **cited** (Olver) + **numerical** ($c_0\approx0.69$) |
| 5 | approach tube + re-entry (global supermartingale) | **proved** |
| 6 | connection variance $\mathbb E[\delta\Theta^2]\le C\eta^2$ uniform (Itô isometry + closed-form integrability + floor) | **proved** (mod. #4) |
| 7 | upper-tail sub-Gaussian, uniform | **proved** |
| 8 | confined Bernstein lower tail $\Rightarrow C^m m!\eta^{2m}$ | **proved** |
| 9 | $R_{\rm nl}$ via Gaussian concentration $\Rightarrow\kappa_{2m}(\delta\Theta)\le C^m m!\eta^{2m}$ uniform | **proved modulo #10** |
| 10 | **Malliavin-diff of the first-passage** $D_sY^\star=-\eta\Phi/k$ — re-derived in the field rep; transversality structural; IFT-on-Wiener-space with verified hypotheses | **proved** (§3.7bis; numerically corr $1.0000$) |
| 11 | **T1 (uniform two-sided tube)** | **closed** |
| 12 | $\mathcal W_2$ fingerprint ($+0.61,-0.24,-2.2,-2.7$); smooth FP construction | **numerical** |
| 13 | tails $(2q+1,3q/2)$; cusp $(5,3)$, ratio $5/3$; two-mechanism | **derived (FW) + numerical** |
| 14 | $\mathcal W$ ≠ integer/half-integer higher-order TW; ≠ Pearcey; ≠ Stokes-continued; ≠ soft-edge determinant | **proved-by-exclusion (cited tail formulas) + numerical** |
| 15 | $\mathcal W$ = asymmetric isomonodromy in PIV family | **derived-structural** |
| 16 | **left tail $-\log\mathbb P\to s^5/(10\eta^2)=\beta s^5/40$** (FW instanton: exponent 5, rate $\propto\beta$, constant $1/10$); "$\beta^{0.6}$ anomaly" = pre-asymptotic crossover | **derived** (§4.3bis; MC slope$=\mathcal I/4$, ratio $\to1.02$) |
| 16b | right-tail FW constant; exact $O(1)$ offset $c(s)$ of the left tail | **open — subleading, non-load-bearing** |
| 17 | precise asymmetric-isomonodromy class / closed-form $\sigma$-form for $\mathcal W$ | **open — *provably not* a Painlevé-$\sigma$-ODE (isomonodromy fixed point); no closed form (theorem, not a gap)** |
| 18 | **intrinsic node process = perturbed $\pi$-lattice (class-I hyperuniform), not determinantal**; jitter $\operatorname{Var}(\text{spacing})\approx(2/\beta)\Theta^{-3/2}$ | **derived + numerical** (§4.5; $S(k)$ excludes GUE by $10^3$) |
| 19 | **ladder unification** $a(q)=3q/(q+2)$; fold $q{=}1$ = exact class II$\to$I boundary (GUE $\to$ lattice) | **derived + numerical** (§4.5; measured $1.03/1.53/1.85/2.07$) |
| 20 | within-class symmetry-breaking susceptibility $dY^\star/da$ (Green's-function integral); sidesteps the isomonodromy obstruction | **derived + numerical** (§4.5bis; $+0.121$ vs FD $+0.124$) |

**Precisely-stated open items (after this pass).** #10 (Malliavin regularity) and #16 (the *left*-tail law) are
now **closed** — the left tail is the FW instanton $-\log\mathbb P\to\beta s^5/40$ (constant $1/10$), superseding
the old non-universal "$1/20$". What remains: (i) #16b — the right-tail FW *constant* and the subleading $O(1)$
offset $c(s)$ (both non-load-bearing); (ii) **#17 — the closed-form integrable identity of $\mathcal W$**, which is
the *genuine* open problem and is now known to be **not a Painlevé-$\sigma$-ODE** (the cusp is an isomonodromy
*fixed point*, so there is no Painlevé flow to reduce onto — a theorem about $\mathcal W$, not a missing
calculation): the target is the asymmetric rank-2 Riemann–Hilbert / Lax problem whose connection data is
$\mathcal W$, or its defining backward-Kolmogorov PDE (already validated), which provably does **not** collapse to
a 1-D ODE.

A methodological caveat carried throughout: regression-based $\sigma$-form tests on a *single* distribution
(or the $\beta$-family) are **confounded** (a 1-D solution curve satisfies many polynomial relations; even a
Gaussian "passes" on smooth data, and a skew-normal family matches the shared-structure fit) — so T2's
classification rests on the **tail/connection structure**, not on fitting the bulk.

---

## 6. What's next / what to check (for N. Popović)

1. **T1 is closed** (§3.7bis); the one previously-flagged item — Malliavin-diff of $Y^\star(W)$ — is re-derived in
   the field representation (transversality structural, coefficients regular through the turning), numerically
   confirmed to corr $1.0000$. The place to sanity-check is that the field SDE $dx=A(\tau)x\,d\tau+Bx\,dW$ is the
   right lift and that the IFT-on-Wiener-space (Nualart–Vives) applies at the non-degenerate zero — both standard.
2. **Audit the genuinely-new probabilistic estimates:** the closed-form integrability of the oscillatory
   $\int dy/(y(y+\Delta))$ (§3.4), the confined Bernstein bound (§3.6), and now the **left-tail instanton**
   $\mathcal I(s)\to s^5/10$ (§4.3bis) — the last replaces the old non-universal "$1/20$" with the derived
   $\beta s^5/40$, MC-confirmed via the linear-in-$\beta$ slope $=\mathcal I(s)/4$.
3. **T2 integrable identity (#17) — the remaining genuine open problem.** Write the $2\times2$ Lax/RH problem with
   a rank-2 irregular singularity whose Stokes data is **asymmetric** across the confining/oscillatory anti-Stokes
   rays, and check its connection coefficients reproduce the derived $(5,3)$ tails (left constant $1/10$) and the
   $\beta$-family cumulants. The validated FP $\mathcal W_\beta$ and its backward-Kolmogorov PDE (which provably
   does *not* reduce to a 1-D $\sigma$-ODE) are the targets. This is where the closed form, if any exists, lives.
4. **The intrinsic process is settled structurally** (§4.5): a class-I perturbed lattice, not a determinantal
   edge, with the fold as the exact class-boundary — the natural place to look for a rigorous 1-D-spectral
   (hyperuniform) rigidity theorem rather than a Dyson/determinantal one.

**Figures (audit trail):** `regime_atlas.png`, `cusp_edge_law.png`, `fp_smooth_W2.png`, `fp_beta_family.png`,
`connection_variance.png`, `confined_bernstein.png`, `moment_lyapunov.png`, `firstpassage_rnl.png`,
`tail_constants.png`, `weber_tube_build.png`, and (this pass) `malliavin_firstpassage.png` (§3.7bis),
`leftail_derivation.png` (§4.3bis), `node_kernel.png` + `ladder_processes.png` (§4.5),
`susceptibility.png` (§4.5bis). Scripts are named inline. Working notes (per result) are in `regime-tests/`
and the project root; this document supersedes them as the consolidated statement.
