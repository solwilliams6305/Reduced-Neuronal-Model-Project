# Headline — the two-sided uniform tube and the construction of the new cusp edge law 𝒲_β

_June 2026. The capstone of the Tier-3 program. Figures `coupled-atlas/figures/cusp_edge_law.png` (the new
law) and the prior `no_early_escape_test.png`, `prefactor_uniformity_test.png`. Scripts
`cusp_edge_law.py`, `no_early_escape_test.py`. Builds on `NO_EARLY_ESCAPE_NOTES.md`,
`G1_PREFACTOR_UNIFORMITY_NOTES.md`, `UNIFORM_WEBER_TUBE_THEOREM_SCOPING.md`. Tags
[PROVED]/[NUMERIC]/[CITED]/[HEURISTIC]/[CONJECTURAL] used precisely; a numerical match is **never** called a
proof._

---

## Stage 1 — re-entry closed: the one-sided bound is now fully ρ-uniform, no gaps

The one remaining hole in G1 was the upper-boundary re-entry (a path that exits the strip upward and returns
to escape later). It closes cleanly using a structural fact:

**The nonlinearity confines from above.** The fluctuation drift $b(x)=-x(2\bar p+x)$ has $b(x)\le-2\bar p x$
for $x>0$ and $b(x)\to-\infty$ like $-x^2$ — so the process **cannot run away upward**; large $x>0$ is
strongly mean-reverting. There is therefore no genuine "upper escape," only the downward one.

**Globally-valid supermartingale (linearly-truncated potential).** Replace the bistable potential $\Phi$ by
$$U(x,\tau)=\begin{cases}\Phi(x,\tau)=\bar p x^2+\tfrac13x^3,& x\le x_+\\[2pt]\Phi(x_+,\tau)+\Phi_x(x_+,\tau)(x-x_+),& x>x_+\end{cases}$$
($C^1$, with a **concave** kink at $x_+$). Set $M_\tau=\exp(\psi+\tfrac{2}{\eta^2}U)$, $\dot\psi=-(2\bar p+2x_+)$.

- For $x\le x_+$: the original computation gives $\mathcal LM/M=\dot\psi+2\bar p+2x+\tfrac{2}{\eta^2}\dot{\bar p}x^2
  =2(x-x_+)+\tfrac{2}{\eta^2}\dot{\bar p}x^2\le0$ (using $x\le x_+$ and $\dot{\bar p}\le0$ pre-turning). **[PROVED]**
- For $x>x_+$ (linear tail, $U_{xx}=0$): the Hamilton–Jacobi term $U_\tau+bU_x+U_x^2=
  \dot{\bar p}x_+^2+2\dot{\bar p}x_+(x-x_+)-c_+(x-x_+)(2\bar p+x+x_+)$ with $c_+=\Phi_x(x_+)>0$ — **every term
  $\le0$** for $x>x_+>0$, $\dot{\bar p}\le0$. So $\mathcal LM/M\le\dot\psi\le0$. **[PROVED]**
- The concave kink at $x_+$ contributes a non-positive local-time term (Itô–Tanaka), preserving the
  supermartingale. **[PROVED]**

Hence $M$ is a supermartingale on **all** of $\{x\ge-L\}$ with no upper boundary; optional stopping at
$\tau_{-L}\wedge T$ (only the escape stop) gives the no-early-escape bound **with no re-entry term**:
$$\boxed{\ \mathbb P(\text{early escape})\le C\,e^{-h^{\star2}/6},\quad \ln C=\underbrace{2\!\int_0^T\!\bar p\,d\tau}_{\mathcal E}+2x_+T\ }$$
both pieces bounded uniformly in $\rho$ (G1, §2 there). **The one-sided pre-turning bound is now fully
ρ-uniform with no remaining gaps.** [PROVED, modulo the cited deterministic floor G2 / $\dot{\bar p}\le0$.]

---

## Stage 2A — the two-sided uniform tube (C1)

**Upper side is the easy direction.** For $x>0$, $b(x)\le-2\bar p x$, so by the SDE comparison theorem the
fluctuation is dominated by the linear OU process $dx^{\rm OU}=-2\bar p\,x^{\rm OU}d\tau+\eta dW$ (same noise).
Thus $\mathbb P(\sup x\ge h\sigma)\le\mathbb P(\sup x^{\rm OU}\ge h\sigma)\le C'e^{-h^2/2}$ — the standard
Berglund–Gentz linear tube, **valid with no smallness assumption** because the nonlinearity is favorable
upward. **[PROVED]**

**Two-sided.** Combining with Stage 1 (lower/escape side, rate $h^{\star2}/6$):
$$\mathbb P\Big(\sup_{0\le\tau\le\tau_\ell}\frac{|x(\tau)|}{\sigma_\Delta}\ge h\Big)\le C\,e^{-c\,h^2},$$
constants uniform in $\rho$; the binding (smaller) rate is the lower/escape side. This is **(C1)** of the
scoped theorem, on the approach window $[0,\tau_\ell]$ up to the entry of the $O(\ell)$ turning neighborhood.
**[PROVED, uniform in ρ.]**

**"Through and past the turning."** The tube controls the *approach*; inside the $O(\ell)$ turning
neighborhood the canard $\bar p$ ceases (the attracting branch ends at $Y=0$) and **all** trajectories peel
off — there the tube estimate hands over to the *connection problem*, i.e. to (C2). So "past the turning" is
precisely where (C1) ends and the **edge law** begins; it is not a confinement statement but a
distributional one. [DERIVED — the scope is set by where $\bar p$ exists.]

---

## Stage 2B — construction of the new cusp edge law 𝒲_β (the headline)

**Definition.** At the cusp ($\Delta=0$) the inner equation is the **stochastic Weber operator**
$u''=(\operatorname{sign}(Y)Y^2-\eta\dot W)u$. Cole–Hopf $p=u'/u$ gives the swept Riccati
$$dp=(\operatorname{sign}(Y)Y^2-p^2)\,d\tau+\eta\,dW,\qquad \tau=Y_0-Y,$$
started recessive on the canard $p(Y_0)=+\sqrt{V}$. Its **first explosion** $p\to-\infty$ at $Y=Y^\star$ is
the first node of the recessive solution — i.e. where the noisy recessive→continuation **connection first
vanishes**. Define
$$\boxed{\ \mathcal W_\beta:=\text{law of the standardized }Y^\star,\qquad \beta=4/\eta^2.\ }$$
This is the **noisy parabolic-cylinder connection data**, and the **cusp ($q=2$) analogue of the
Ramírez–Rider–Virág stochastic-Airy/Tracy–Widom construction ($q=1$)**. [DEFINITION — exact.]

**Computed at β=2 ($\eta=\sqrt2$), 1.4×10⁶ Monte-Carlo realizations** (orientation $s=+(Y^\star-\text{mean})/\text{std}$,
the eigenvalue convention):

| cumulant | constructed 𝒲₂ | cusp fingerprint (target) | match |
|---|---|---|---|
| skewness | **+0.607** (SE 0.002) | +0.61 | ✓ |
| excess kurtosis | **−0.237** (SE 0.004) | −0.24 | ✓✓ (orientation-free) |
| κ₅ | **−2.115** | −2.2 | ✓ (within 4%) |
| κ₆ | **−2.706** | −2.7 | ✓✓ (orientation-free) |

**All four cumulants match** — the two orientation-independent ones (κ₄, κ₆) essentially exactly, the two
odd ones (skew, κ₅) in magnitude to within sampling/discretization (~4%). **[NUMERIC — strong.]**

**Identification.** The object defined purely as the *noisy parabolic-cylinder connection data* reproduces
the independently-established cusp universality fingerprint. This is the **construction and numerical
identification of the genuinely new T2 edge law** $\mathcal W_\beta$ — the law the falsifier campaign showed
descends from no known parent. **[NUMERIC identification; not an analytic proof — see ledger.]**

**Tails (honest).** Left tail (late escape, $s\to-\infty$) effective exponent **2.5**; right tail (early
escape) **1.84** — the **left is steeper than the right**, the correct Freidlin–Wentzell asymmetry
(predicted exponents 5 vs 3). The *absolute* exponents are **pre-asymptotic**: at the accessible depth
($t\lesssim3.4$, $P\gtrsim10^{-5}$) the law has not reached the $e^{-|s|^5/20}$ asymptotic regime (the
crossover sits near $t\approx3.5$). So the tail asymmetry is **confirmed**; the exponent-5 value itself is
**consistent but not numerically confirmed** (would need ~$10^9$ rare-event samples). [NUMERIC + honest caveat.]

---

## The full theorem (C1)+(C2)+(C3), at the strongest honestly-supportable level

> **Theorem (uniform parabolic-cylinder tube + edge law), current status.**
> For the inner Riccati SDE through the cusp merge, with $\rho=\Delta/\ell$ and barrier scale
> $h^\star=4\bar p_{\min}^{3/2}/\eta$:
>
> **(C1)** [PROVED, uniform in ρ on the approach window] two-sided sub-Gaussian confinement
> $\mathbb P(\sup|x|/\sigma_\Delta\ge h)\le C e^{-c h^2}$, $c,C$ independent of $\rho$ (lower-side rate
> $h^{\star2}/6$ via the cubic-scale-function supermartingale, re-entry closed; upper-side $h^2/2$ via OU
> comparison; prefactor bounded via the exposure integral).
>
> **(C2)** endpoints: $\rho\to\infty$ → $\mathrm{TW}_\beta$ [CITED, RRV]; $\rho\to0$ → $\mathcal W_\beta$,
> **constructed** as the noisy parabolic-cylinder connection data and **numerically identified** with the
> cusp fingerprint (skew +0.61, exk −0.24, κ₅ −2.2, κ₆ −2.7) [NUMERIC].
>
> **(C3)** uniform crossover: the **tube is uniform in ρ** [PROVED]; the **law interpolation**
> $\mathrm{TW}_\beta\!\leftrightarrow\!\mathcal W_\beta$ is [NUMERIC] (the crossover fingerprint work).

### Ledger — proved / numerically validated / conjectural

| statement | status |
|---|---|
| fluctuation = bistable cubic SDE; barrier $\tfrac43\bar p^3$ | **PROVED** |
| one-sided no-early-escape $P\le Ce^{-h^{\star2}/6}$, rate $c=\tfrac16$ ρ-independent | **PROVED** (re-entry now closed) |
| prefactor $C(\rho)$ bounded uniformly (exposure, monotone, min at merge) | **PROVED** |
| two-sided tube (C1), uniform in ρ, on the approach window | **PROVED** |
| floor $\bar p_{\min}(\rho)\ge c_0$, $\dot{\bar p}\le0$ pre-turning, turning curvature | **CITED** (Olver uniform PC) |
| $\rho\to\infty$ endpoint = $\mathrm{TW}_\beta$ | **CITED** (RRV) |
| **$\mathcal W_\beta$ defined as noisy PC connection data** | **PROVED** (definition) |
| **$\mathcal W_\beta$ reproduces the cusp fingerprint (4/4 cumulants)** | **NUMERIC — strong** |
| tail asymmetry left-steeper-than-right (FW 5 vs 3) | **NUMERIC** |
| exponent-5 left tail $e^{-|s|^5/20}$ exactly | **NUMERIC, pre-asymptotic — not confirmed** |
| law interpolation $\mathrm{TW}_\beta\leftrightarrow\mathcal W_\beta$ uniform in ρ | **NUMERIC** |
| $\mathcal W_\beta$ = a Painlevé-IV / closed-form connection law (analytic) | **CONJECTURAL** |
| (C1) inside the $O(\ell)$ turning neighborhood (post-turning two-sided) | **OPEN** |
| stochastic uniform PC connection bounds (scoping §3(ii)(a)) | **OPEN** |

### What remains conjectural / open (no dressing-up)
1. **Analytic identity for $\mathcal W_\beta$.** We have a *construction* (the stochastic Weber explosion law)
   and a *numerical* fingerprint match — **not** a closed-form/Painlevé-IV expression nor an analytic proof
   that the connection data equals the cusp law. That analytic identification is the deepest part of T2 and
   stays conjectural.
2. **The exponent-5 tail** is consistent but pre-asymptotic numerically.
3. **(C1) through the $O(\ell)$ turning core** (the two-sided estimate *inside* the connection region) and the
   **stochastic uniform PC connection bounds** remain open (scoping §3(ii)(a) / §4 post-turning).

**Net.** (C1) — the uniform two-sided tube — is **proved** (uniform in ρ, both sides, prefactor bounded,
re-entry closed) on the approach. (C2) — the prize — is **achieved at the level of construction + numerical
identification**: the new universality law $\mathcal W_\beta$ is *defined* as the noisy parabolic-cylinder
connection data and *reproduces the cusp fingerprint to within a few percent on all four cumulants*. What is
honestly **not** done is the closed-form/analytic identification of $\mathcal W_\beta$ (the Painlevé-IV
content) and the connection estimate inside the turning core. The headline stands: **the genuinely new cusp
edge law is constructed and numerically identified; the tube that delivers it is uniform and proved.**
