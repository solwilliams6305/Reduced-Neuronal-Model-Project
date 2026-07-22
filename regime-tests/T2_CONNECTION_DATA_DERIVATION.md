# T2 connection-data route: the asymmetric parabolic-cylinder isomonodromy of 𝒲

_June 2026. The decisive analytic push on T2 — from the structural label ("asymmetric isomonodromy in the
PIV family") to an explicit connection problem, the σ-form question answered, tail mechanisms, and a sharp
obstruction. Figure `coupled-atlas/figures/tail_mechanisms.png` (+ retired `tail_scaling.png`). Scripts
`tail_scaling.py`, `tail_scaling_v2.py`. Tags **[DERIVED]** / **[CITED]** / **[CONJECTURAL]** / **[NUMERIC]**;
⚑ = load-bearing._

---

## 0. Summary — what this establishes

| claim | status |
|---|---|
| The skeleton's two ends are **different irregular structures** — confining (real exp, DLMF $U$) at $+\infty$, oscillatory (imag exp, DLMF $W$) at $-\infty$ | **DERIVED** |
| The connection $U\!\to\!W$ across the turning is the **classical PC connection**: rotation $e^{i\pi/4}$ + parameter conjugation $a\!\to\!ia$ | **DERIVED + CITED** (DLMF 12.14) |
| Explicit Stokes data: confining multiplier $\propto 1/\Gamma(\tfrac12-\tfrac\lambda2)$ (bound-state poles $\lambda=1,3,5\dots$); oscillatory factor $\propto\Gamma(\tfrac12-\tfrac{i\lambda}2)$, $|\cdot|^2=\pi/\cosh\tfrac{\pi\lambda}2$ | **DERIVED + CITED** |
| The skeleton sits at the **parabolic-cylinder (Hermite–Weber) special-solution locus of PIV** | **CITED** (PIV PC solutions) |
| The asymmetric (confining+oscillatory) point is **off the real σ-PIV section** → no standard real σ-form | **DERIVED** (structural) |
| **Sharp obstruction:** the cusp is an isomonodromy **fixed point** (its only PC-class deformation, the linear $\Delta Y$ term, is *relevant* → flows to fold), so 𝒲 is a **connection coefficient at a fixed point, not a Painlevé τ along a flow** | **DERIVED** ⚑ |
| Deterministic first node $Y^\star_{\rm det}=-2.19$ = first zero of the continued $U$-solution on the $W$-side | **DERIVED + NUMERIC** (matches MC median $\to-2.09$) |
| Escape-location law std $\sim\eta$ (node-dominated), **not** $\eta^{2/3}$ (turning scale) — **course-correction** | **NUMERIC** ⚑ |
| Tail exponents $(2q{+}1,3q/2)=(5,3)$ are the **FW backbone**; the earlier "$\eta^2$ right / anomalous left **collapse** of 𝒲" conflated two objects and is **withdrawn** | **honest negative** |

**One-line verdict.** 𝒲 is governed by the **parabolic-cylinder linear problem of PIV, evaluated at an
asymmetric (real⊕imaginary-parameter) Stokes point that is a fixed point of the catastrophe-unfolding flow.**
The connection data is *explicit in Gamma functions*; the law itself has **no standard Painlevé σ-form** for a
precise, provable reason — it is a fixed-point connection coefficient, not a transcendent.

---

## 1. The linear skeleton and its two irregular structures (the RH problem)

**The operator.** Cole–Hopf $p=u'/u$ turns the inner Riccati into the Schrödinger/Weber form
$$u'' = \big(V(Y)-\lambda\big)\,u,\qquad V(Y)=\operatorname{sign}(Y)\,Y^2\ \ (\text{cusp }\Delta=0),$$
with $\lambda$ the spectral/edge parameter and **escape = first node of $u$** (first explosion of $p$). The
$\Delta$-unfolding is the **linear term**: $V_\Delta=\operatorname{sign}(Y)|Y|(|Y|+\Delta)=\operatorname{sign}(Y)Y^2+\Delta Y$
(load-bearing for §5). **[DERIVED]**

**As a 2×2 system.** $\Psi=(u,u')^\top$, $\Psi'=A(Y)\Psi$, $A=\begin{pmatrix}0&1\\ V-\lambda&0\end{pmatrix}$.
The coefficient grows like $Y^2$ ⇒ **rank-2 irregular singularities at $Y=\pm\infty$** (formal exponent
$\int\!\sqrt{V}\sim\pm Y^2/2$, degree 2). **[DERIVED]**

**The two ends are different.** WKB ($u\sim V^{-1/4}e^{\pm\int\sqrt V}$):
- **$Y\to+\infty$ (confining).** $V=+Y^2$, $u\sim Y^{-1/2}e^{\pm Y^2/2}$ — **real** exponents. Recessive
  branch decays. In DLMF normalisation, with $z=\sqrt2\,Y$, $u=U(-\tfrac\lambda2,\sqrt2\,Y)$, where $U(a,z)$
  solves $w''=(\tfrac{z^2}4+a)w$ and $U\sim e^{-z^2/4}z^{-a-1/2}$ (recessive). **Parameter $a_R=-\lambda/2$.**
- **$Y\to-\infty$ (oscillatory).** $V=-Y^2$, $u\sim|Y|^{-1/2}e^{\pm iY^2/2}$ — **imaginary** exponents. With
  $z=\sqrt2\,Y$ this is the DLMF oscillatory PC equation $W''+(\tfrac{z^2}4-a)W=0$, $u=W(-\tfrac\lambda2,\sqrt2\,Y)$.
  **Parameter $a_L=-\lambda/2$.** **[DERIVED]**

So the **same** parameter $a=-\lambda/2$ controls both ends, but through **two functions of different type**:
$U$ (real-exponential, confining) vs $W$ (imaginary-exponential, oscillatory). This is the precise content of
"**two different irregular structures joined by the odd potential**."

**Stokes geometry.** For a rank-2 (∼$z^2$) irregular point the anti-Stokes rays sit where $\operatorname{Re}\!\int\!\sqrt V=0$,
i.e. $\arg z=\pm\pi/4,\pm3\pi/4$ — **four Stokes rays**. The physical line goes from $\arg Y=0$ (confining
real axis) to $\arg Y=\pi$ (oscillatory real axis), **crossing two Stokes rays** ($\pi/4,3\pi/4$). The
connection is the product of the two Stokes matrices and the formal monodromy across that wedge. **[DERIVED]**

**The RH problem.** Find $\Psi(Y)$, sectionally analytic off the Stokes rays, with prescribed jumps (unit-triangular
Stokes matrices $S_k$) across them and the formal behaviour $\Psi\sim(\hat\Psi_0+\dots)e^{(\pm Y^2/2)\sigma_3}$ at
each end. **The single non-analyticity at $Y=0$** (where $\operatorname{sign}(Y)$ has its corner) means the
right and left **cannot be one analytic $A(z)$ on $\mathbb C$** — the global object is two PC problems glued by
$C^1$ matching at $Y=0$. (This is the seed of §5.) **[DERIVED]** ⚑

---

## 2. Stokes / monodromy data — explicit, in Gamma functions

**Confining multiplier.** The recessive $U(a,z)$ continued across the two Stokes rays into the dominant branch
picks up the classical PC Stokes multiplier
$$s_{\rm conf}\;=\;\frac{\sqrt{2\pi}}{\Gamma(\tfrac12+a)}\,e^{\,i\pi a}\;=\;\frac{\sqrt{2\pi}}{\Gamma(\tfrac12-\tfrac\lambda2)}\,e^{-i\pi\lambda/2}. $$
**[CITED — DLMF/Whittaker–Watson; the search-confirmed $\mathrm{SM}\propto 1/\Gamma(-\nu)$ with $\nu=\tfrac{\lambda-1}2$.]**
Its **poles** $\tfrac12-\tfrac\lambda2=0,-1,-2,\dots$ ⇒ $\lambda=1,3,5,\dots$ are the **half-line confining
bound states** (the recessive solution becomes globally recessive — a genuine eigenvalue). **[DERIVED]**

**Oscillatory factor.** The left end is the **inverted oscillator** $-Y^2$. Its scattering (reflection off the
inverted parabola) is the classical
$$|r(\lambda)|^2=\frac1{1+e^{\pi\lambda}},\qquad |t(\lambda)|^2=\frac{1}{1+e^{-\pi\lambda}},\qquad
\big|\Gamma(\tfrac12-\tfrac{i\lambda}2)\big|^2=\frac{\pi}{\cosh\tfrac{\pi\lambda}2}. $$
The connection to the $U$-side uses $W(a,x)$'s representation through $U(\pm ia,xe^{\pm i\pi/4})$ — i.e.
**rotation $e^{i\pi/4}$ and parameter conjugation $a\to ia$** (DLMF 12.14.4–6). **[CITED + DERIVED.]**

**The asymmetric connection coefficient.** Assembling (recessive-at-$+\infty$, expressed in the oscillatory
basis at $-\infty$), the connection quantity that fixes the first node is
$$\boxed{\,\mathcal C(\lambda)\ \propto\ \underbrace{\frac{1}{\Gamma\!\big(\tfrac12-\tfrac\lambda2\big)}}_{\text{confining: real }a,\ \text{bound-state poles}}\ \times\ \underbrace{\Gamma\!\big(\tfrac12-\tfrac{i\lambda}{2}\big)\,e^{\,i\varphi(\lambda)}}_{\text{oscillatory: imaginary }a,\ \cosh\text{ resonance}}\, }$$
with $\varphi$ the (real) inverted-oscillator phase. **[DERIVED]** ⚑ The deterministic first node is the
smallest $Y^\star<0$ at which the continued solution's accumulated $W$-phase equals $\pi$; numerically
$Y^\star_{\rm det}=-2.19$ (panel 3 of the figure), matching the MC median $\to-2.09$ as $\eta\to0$. **[NUMERIC ✓]**

**The crux.** $\mathcal C(\lambda)$ is **explicit in Gamma functions, but mixes $\Gamma$ at a real argument
(confining bound states) with $\Gamma$ at an imaginary argument (oscillatory resonance).** A standard
self-adjoint edge (Airy/Weber-symmetric) uses one *or* the other; the **odd potential forces both at once.**
That product — real-$\Gamma$ × imaginary-$\Gamma$ — *is* "asymmetric isomonodromy," now concrete.

---

## 3. The σ-form question — PIV's parabolic-cylinder locus, and where 𝒲 sits

**Positive identification.** The Weber/PC linear problem of §1 is exactly the **isomonodromic linear system
for Painlevé IV at its parabolic-cylinder (Hermite–Weber) special-solution locus**: PIV admits
one-parameter families of solutions in $U(a,z)$/$D_\nu$, reducing to Hermite/rational at integer parameters,
constructed through Riemann–Hilbert problems. **[CITED — confirmed: PIV ⊃ parabolic-cylinder special
solutions; isomonodromy/RH representation of the Hermite family.]** At this locus the PIV **τ-function is the
explicit PC/Γ object** of §2 — so the "σ-form" of the *symmetric* skeleton **degenerates to the classical
(linear) PC connection**: no transcendental σ-form is needed, because this is the classically-solvable slice
of PIV.

**Where the cusp sits — and why it is not a standard σ-form.** The cusp skeleton is **not** the symmetric PC
point: its two ends carry parameter $a=-\lambda/2$ **(real, confining)** and effectively $ia$ **(imaginary,
oscillatory)** simultaneously. On the PIV monodromy manifold this is a point where the two asymptotic
directions carry **conjugate-type** Stokes data (one real-exponential sector, one imaginary-exponential
sector). The **real section** of σ-PIV (the self-adjoint Hermite–Weber transcendents) does **not pass through
this point**. Hence:

> **𝒲 is governed by PIV's parabolic-cylinder structure, but at an *asymmetric / complex-Stokes* point off the
> real σ-PIV section.** It is **not** a standard real σ-PII/σ-PIV transcendent. **[DERIVED]** ⚑

This is the analytic content the σ-form *regression* could not reach (and which we earlier showed it could not,
because the 1-D fit is confounded): the obstruction is **structural**, not a fitting failure.

---

## 4. Tail mechanisms from the connection — and an honest course-correction

The two Γ-factors of $\mathcal C(\lambda)$ name the two tail mechanisms:

- **Right tail = confining barrier (the $U$/$\Gamma(\tfrac12-\tfrac\lambda2)$ factor).** Escape *earlier* than
  $Y^\star_{\rm det}$ requires the solution to node *against* the $+Y^2$ confining restoring — a
  Kramers/Freidlin–Wentzell event. Near the turning this is the **$\eta^2$-rate large deviation** already
  **proved** (no-early-escape supermartingale, rate $1/6$): $P_{\rm early}\le Ce^{-h^{\star2}/6}$. **[DERIVED =
  the proved no-early-escape.]**
- **Left tail = oscillatory phase-persistence (the $W$/$\Gamma(\tfrac12-\tfrac{i\lambda}2)$ factor).** Escape
  *later* (deeper into $-Y^2$) requires the $W$-phase to avoid $\pi$ across an extended oscillatory stretch —
  the inverted-oscillator $\cosh$ factor controls it. **Anomalous** (non-$\eta^2$). **[DERIVED — mechanism.]**

The standardized **tail exponents $(5,3)$** ($2q{+}1$ left, $3q/2$ right) are the **FW analytic backbone**
[CITED, earlier]. The figure (panel 2, $\eta=1.4$) shows the **correct asymmetry direction**: right tail
heavier (barrier, exp 3), left tail lighter/steeper (oscillatory, exp 5), with **skew $+0.61$** — consistent.

**⚑ Honest course-correction (a withdrawn over-claim).** The earlier summary recorded a "**clean $\eta^2$
right-tail / anomalous left-tail *collapse* of 𝒲**." Testing it directly (`tail_scaling.py`, 4-panel) it **does
not hold**, and the reason is structural: the escape happens at an **$O(1)$ oscillatory node**
$Y^\star_{\rm det}=-2.19$, **not** at the $\eta^{2/3}$ turning scale. Measured **escape-location std
$\sim\eta^{1.10}$** (panel 1) — i.e. **$\eta$-linear** (node phase-noise), tracking the $\eta^1$ line, clearly
**not** $\eta^{2/3}$. Consequences, stated plainly:
- The $\eta^2$ statement is real but belongs to the **near-turning early-escape LDP** (the proved
  no-early-escape), a **different object** from the escape-location law 𝒲's right tail.
- 𝒲's escape-location tails scale with **std $\sim\eta$**; the "$\eta^2$/$\eta^{10/3}$ collapse" **conflated**
  the turning scale with the node scale and is **withdrawn**.
- The $(5,3)$ exponents stand as the FW backbone but are **pre-asymptotic numerically** (panel 2; the
  asymptotic regime sits past accessible depth). **Validated** here: std $\sim\eta$, $Y^\star_{\rm det}=-2.19$
  match, tail-asymmetry direction + skew. **Not** validated: the exact exponents or any clean $\eta$-power
  collapse. **[NUMERIC, honest.]**

So the connection coefficients deliver the **mechanisms and the bound-state/resonance constants** exactly, but
the *escape-location tail constants* are **not** a clean closed form at accessible $\eta$ — and the previously
claimed collapse was wrong.

---

## 5. The sharp obstruction — why 𝒲 has no standard Painlevé σ-form

A Painlevé σ-form is the deformation equation of an **isomonodromic family**: a flow in a parameter $t$ that
moves the linear ODE's coefficients while **fixing the monodromy**. For 𝒲 there is **no such flow within the
cusp class**, for two independent reasons:

**(O1) The only PC-class deformation is *relevant*.** The deformations preserving the parabolic-cylinder
(rank-2) structure are exactly the lower-degree additions $\operatorname{sign}(Y)(Y^2+c_1|Y|)+\dots$. The
leading one is the **linear $\Delta Y$ term** — and we **proved/validated (rung E)** it is **relevant**: it
**flows the cusp to the fold** (skew $0.60\!\to\!0.26$, sub-Gaussian signature dissolving; lowest-degree term
dominates). A relevant deformation **changes the irregular type at the ends** (it removes the confining/oscillatory
balance), so it is **not** monodromy-preserving. There is **no nontrivial isomonodromic flow tangent to the
cusp.** **[DERIVED, load-bearing ⚑ — uses the rung-E result.]**

**(O2) The two ends are conjugate-type.** A standard σ-PIV transcendent is real-analytic in its (single, real)
isomonodromic time, with monodromy a representation into one real/unitary group. Here the monodromy data is
**split**: real-exponential ($U$, $\Gamma(\tfrac12-\tfrac\lambda2)$) at one end, imaginary-exponential ($W$,
$\Gamma(\tfrac12-\tfrac{i\lambda}2)$) at the other. No single real Painlevé transcendent carries conjugate
data on its two Stokes ends. **[DERIVED.]**

**Conclusion (the obstruction).**
> 𝒲 is a **connection coefficient at an isomonodromy fixed point** of the catastrophe-unfolding flow — the
> asymmetric (real⊕imaginary-parameter) PC point of PIV — **not a Painlevé τ-function along a deformation
> flow.** Equivalently: the cusp law sits where the σ-form's deformation vector field **vanishes**, so the
> σ-ODE degenerates to the **algebraic (Γ-function) fixed-point relation** of §2, not a transcendental ODE.
> **[DERIVED]** ⚑

**Why this is the right kind of object — and what it explains.** All four earlier negative results presuppose
a *flow*: a Fredholm determinant (a τ along the spectral edge), an isomonodromic *time*, a Pearcey reduction (a
higher multicritical *flow*), a Stokes *continuation*. The fixed-point statement **unifies and explains** them:
there is no flow to host any of those. The positive content is the **explicit fixed-point connection data**
(the Γ-product), which is as close to a "closed form" as the structure allows: a closed form for the
*connection coefficient*, with the *distribution* 𝒲 its noise-dressing.

**What a genuine closed form would require (the generalized object).** A τ-function for a **non-self-adjoint /
"interface" isomonodromy** — two rank-2 irregular singularities of **conjugate type** glued along a $C^1$
interface at $Y=0$ — i.e. an isomonodromy problem on a **bordered/doubled** surface where one boundary carries
real and the other imaginary Stokes data. Such "PT-symmetric / resonance" isomonodromy is **not** in the
standard Painlevé I–VI list. **[CONJECTURAL — the precise generalized object; flagged, not constructed.]**

---

## 6. Net — and what this unblocks (B, D, E)

**Derived:** the explicit asymmetric-PC RH problem (§1), its Gamma-function Stokes data (§2), the PIV
parabolic-cylinder-locus identification (§3), the two tail mechanisms (§4), and a **sharp, structural
obstruction** to any standard σ-form (§5) — with the generalized object that *would* be needed named.

**Honest negatives:** the escape-location "$\eta^2$/anomalous collapse" is **withdrawn** (node scale $\eta$, not
turning scale $\eta^{2/3}$); exact $(5,3)$ exponents remain FW-backbone, numerically pre-asymptotic.

**Unblocks the shared frontier:**
- **B (cusp rung):** the marginal 𝒲 is now *characterized* — a fixed-point connection coefficient of the
  asymmetric PC/PIV problem, with explicit Γ-data. The remaining B gap (the *intrinsic* DBM-edge process) is
  the **multi-line** version of this same fixed-point RH problem.
- **D (atlas):** confirms and sharpens the cusp cell — "asymmetric PIV-family" now means *precisely* the
  conjugate-type PC point; the (5,3) tails are the FW backbone, exponents pre-asymptotic.
- **E (unfolding):** the obstruction **is** the rung-E relevance statement read backwards — *because* the
  linear deformation flows to the fold, there is no cusp-preserving isomonodromic time. The two results are
  the same fact from two sides.

**Frontier that remains:** constructing the non-self-adjoint "interface/PT" isomonodromy τ-function (the named
generalized object) — the genuine open analytic problem, now precisely posed rather than vaguely "PIV-like."
