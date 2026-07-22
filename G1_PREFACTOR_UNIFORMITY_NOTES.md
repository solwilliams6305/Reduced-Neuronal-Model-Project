# Proof step — G1: prefactor uniformity for the no-early-escape bound

_June 2026. Figure `coupled-atlas/figures/prefactor_uniformity_test.png`; script
`prefactor_uniformity_test.py`. Closes the load-bearing remainder of `NO_EARLY_ESCAPE_NOTES.md`: upgrade
$P(\text{early escape})\le C\,e^{-h^{\star2}/6}$ so the **prefactor $C$** is also bounded uniformly in
$\rho=\Delta/\ell$ through the merge. Tags [PROVED]/[NUMERIC]/[CITED]/[OPEN]/[HEURISTIC]._

---

## 1. What $C$ depends on (explicit, from optional stopping)

The supermartingale $M_\tau=\exp(\psi(\tau)+\tfrac{2}{\eta^2}\Phi(x_\tau,\tau))$ on the strip
$x\in[-L,x_+]$ required $\dot\psi=-(2\bar p_\Delta+2x_+)$ (to dominate the Itô term $\Phi_{xx}=2\bar p+2x$),
and optional stopping gave $P\le e^{-\psi(T)}e^{-2\Phi_\ast/\eta^2}$. Hence the prefactor is **explicit**:
$$\boxed{\ \ln C(\rho)\;=\;\underbrace{2\!\int_0^T\!\bar p_\Delta(\tau)\,d\tau}_{\text{exposure }\mathcal E(\rho)}\;+\;\underbrace{2\,x_+\,T}_{\text{boundary }\mathcal B(\rho)}\ }$$
with $T=Y_0-Y_{\rm end}$ the (fixed) window length and, via $\tau=Y_0-Y$,
$$\mathcal E(\rho)=2\!\int_{Y_{\rm end}}^{Y_{\rm in}}\!\bar p_\Delta(Y)\,dY .$$
$\mathcal E$ is the **integrated relaxation rate** (curvature of the well bottom along the canard); $\mathcal B$
is the **upper-strip / entry term** (from the cutoff $x_+$ that bounds the Laplacian term for $x>0$). There is
no separate entry contribution: the path starts at $x=0$ ($\Phi=0$, $\psi=0$), so $M_0=1$ exactly. A
cutoff-free variant trades $\mathcal B$ for $\tfrac{\eta^2}{2}\int d\tau/|\dot{\bar p}_\Delta|$ (§4).

---

## 2. The exposure is bounded uniformly through the merge — and is *smallest* there [PROVED]

**Closed form (outer).** With $\bar p_\Delta\approx\sqrt{V_\Delta}$, $V_\Delta=Y(Y+\Delta)$ ($Y>0$), set
$u=Y+\tfrac\Delta2$, $a=\tfrac\Delta2$:
$$\int\!\sqrt{Y(Y+\Delta)}\,dY=\tfrac{u}{2}\sqrt{u^2-a^2}-\tfrac{a^2}{2}\ln\!\big(u+\sqrt{u^2-a^2}\big)+\text{const}.$$
This is **elementary, finite, and monotone increasing in $\Delta$** (the integrand $\sqrt{Y(Y+\Delta)}$ is
increasing in $\Delta$ pointwise). As $\Delta\to0$ it decreases to $\int_{Y_{\rm end}}^{Y_{\rm in}}Y\,dY=
\tfrac12(Y_{\rm in}^2-Y_{\rm end}^2)$ — a finite limit. **The merge ($\rho\to0$) is the minimum, not a
blow-up.** [PROVED — explicit integral; numerically the closed form matches the exact-canard exposure to 3
decimals, see §5.]

**Inner correction.** The exact canard cannot follow $\sqrt{V_\Delta}\to0$; it levels at $\bar p_{\min}>0$.
The correction $2\int_{\rm inner}(\bar p_\Delta-\sqrt{V_\Delta})\,dY$ is supported on the inner region of width
$O(\ell)$ and height $O(\bar p_{\min})$, hence $\le 2\,\bar p_{\min}\,\cdot O(\ell)$ — **bounded** (numerically
$\approx1.0$–$1.2$, itself *decreasing* in $\rho$). [PROVED-modulo-ceiling: needs $\bar p_{\min}(\rho)$ bounded
*above*, automatic on bounded $\rho$.]

**Conclusion.** For $\rho\in[0,\rho_0]$ ($\rho_0=O(1)$, through-the-merge range),
$$\mathcal E(\rho)\le \mathcal E(\rho_0)=2\!\int_{Y_{\rm end}}^{Y_{\rm in}}\!\bar p_{\Delta_0}\,dY<\infty
\quad\Rightarrow\quad e^{\mathcal E(\rho)}\ \text{bounded uniformly.}\quad\blacksquare$$

---

## 3. The floor's actual role (an honest correction to the cited framing)

The natural expectation was "the Olver floor $\bar p_{\min}(\rho)\ge c_0$ keeps the exposure from
degenerating at the merge." **This is inverted, and worth stating plainly:**

- The exposure does **not** degenerate at the merge regardless of the floor — it is an *outer-dominated*,
  $\Delta$-monotone integral over a *fixed* window, so it is automatically bounded (and minimal at the merge).
  The bound on $\mathcal E$ uses a **ceiling** on $\bar p$ (automatic), not a floor.
- The floor $\bar p_{\min}\ge c_0$ is essential **elsewhere**: **(i) the exponent** $h^\star=
  \min 2\bar p/\sigma=4\bar p_{\min}^{3/2}/\eta$, which would collapse to $0$ (vacuous bound) without it; and
  **(ii) the strip geometry** in $\mathcal B$ — the strip half-width is the barrier $2\bar p_{\min}\ge2c_0$, so
  the floor keeps the strip from pinching shut at the turning.

So G1 (prefactor uniformity) and G2 (the floor) are **largely independent**: the prefactor is bounded by
outer geometry; the floor's job is the exponent and the strip non-degeneracy. [HEURISTIC→PROVED reclassification.]

---

## 4. The boundary term $\mathcal B$ (partial) and re-entry

$\mathcal B=2x_+T$ depends on the cutoff $x_+$. Two honest issues:

1. **Choice of $x_+$.** Take $x_+=O(\bar p_{\min})$ (the strip scaled to the barrier). Then $\mathcal B=
   O(\bar p_{\min}T)$ — bounded uniformly on $\rho\in[0,\rho_0]$ (ceiling on $\bar p_{\min}$). [PROVED-modulo-ceiling.]
2. **Re-entry.** Stopping at $\tau_{-L}\wedge\tau_{x_+}\wedge T$ bounds escape *before first reaching $x_+$*.
   A path may exit upward through $x_+$ (the safe direction) and later return to escape. Bounding total escape
   needs a restart: from $x_+$, reaching $-L$ has an *even larger* barrier ($\Phi(x_+\!\to\!-L)>\Phi(0\!\to\!-L)$),
   so each excursion contributes $\le$ the same exponential with a smaller constant; the number of independent
   up–down excursions in the window is $O(T/\theta_{\rm relax})=O(\int\bar p)$, giving a **bounded multiplicative
   factor**. This makes $C\to C\cdot O(\mathcal E)$ — still uniformly bounded, but the clean proof of the
   excursion count uniform in $\rho$ is **[OPEN — the remaining delicacy of $\mathcal B$]**.

**Cutoff-free variant.** Using $\dot\psi=-2\bar p-\tfrac{\eta^2}{2|\dot{\bar p}|}$ removes $x_+$ (the helpful
$\tfrac{2}{\eta^2}\dot{\bar p}x^2\le0$ term dominates $+2x$ for large $x$), at the cost of
$\tfrac{\eta^2}{2}\int d\tau/|\dot{\bar p}|$. Near the turning $\bar p_\Delta(Y)=\bar p_{\min}+\tfrac12\bar
p''(Y-Y_\ast)^2+\cdots$, so $|\dot{\bar p}|\sim|Y-Y_\ast|$ and the integral is **logarithmically divergent** at
the turning — regularized by an $\varepsilon$-core handled with the cutoff version. Net: $\mathcal B=
O(\eta^2|\bar p''|^{-1}\ln(1/\varepsilon))+O(\bar p_{\min}\varepsilon)$, bounded given the turning curvature
$\bar p''(Y_\ast)$ — **also an Olver-uniform quantity** (same source as the floor). [PARTIAL.]

---

## 5. Is the floor "as cited" (G2) sufficient?

| needs | does the cited floor $\bar p_{\min}\ge c_0$ suffice? |
|---|---|
| exposure $\mathcal E$ bounded | **not even needed** (outer ceiling does it) — sufficient a fortiori |
| exponent $h^\star$ bounded below | **YES** — this is exactly G2's role |
| strip non-degeneracy in $\mathcal B$ | **YES** (strip half-width $\ge 2c_0$) |
| re-entry excursion count uniform | needs $\int\bar p$ bounded (have it) — **YES** modulo the clean count |
| cutoff-free turning $\ln$ factor | needs the **turning curvature** $\bar p''(Y_\ast)$ uniform — a *mild strengthening* of G2 (same Olver package, not just the floor *value*) |

**Verdict:** the floor as cited is **sufficient for the exposure and the exponent** (the main content). The
boundary term $\mathcal B$ needs, in addition, either a clean uniform excursion-count (Tier-3 bookkeeping) or
the turning curvature $\bar p''(Y_\ast)$ — both available from the same uniform parabolic-cylinder asymptotics
(scoping §3(ii)), so "in the same package," but formally a hair more than the floor *value* alone.

---

## 6. Numerical validation (exact canard, $\eta=\sqrt2$, window $Y\in[0.25,2]$)

| ρ | $\bar p_{\min}$ | $\mathcal E$ exact | $\mathcal E$ (√V closed) | inner corr | $C=e^{\mathcal E}$ |
|---|---|---|---|---|---|
| 0.024 (deep merge) | 0.817 | **5.18** | 3.99 | 1.19 | 178 |
| 0.079 | 0.837 | 5.28 | 4.11 | 1.17 | 196 |
| 0.397 | 0.942 | 5.80 | 4.72 | 1.08 | 330 |
| 0.952 | 1.099 | 6.62 | 5.61 | 1.01 | 746 |
| 1.984 | 1.334 | **7.90** | 6.96 | 0.94 | 2706 |

- **$\mathcal E(\rho)$ bounded and monotone increasing**; ratio max/min over $\rho\in[0.024,1.98]$ is **1.53**.
  $\mathcal E$ at the deepest merge is the **minimum** (5.18), not a blow-up. **[NUMERIC ✓]**
- Closed form matches the exact-canard outer part to 3 decimals (validates §2). **[NUMERIC ✓]**
- $\bar p_{\min}(\rho)$ bounded below at $c_0\approx0.82$ (cusp value, at the merge) and above — the floor.
  **[NUMERIC ✓, consistent with Olver-uniform G2]**
- $C=e^{\mathcal E}$ bounded uniformly (178→2706). **Large in absolute terms** — the exponential-supermartingale
  prefactor is *crude*; the bound bites only for $h^{\star2}/6>\ln C\approx5$–$8$ (small noise). Tightness is a
  *separate* issue from uniformity. **[NUMERIC; honest caveat]**

---

## 7. Status — complete vs blocked

| component | status |
|---|---|
| explicit $\ln C=\mathcal E+\mathcal B$ from optional stopping | **[PROVED]** |
| exposure $\mathcal E(\rho)$ bounded uniformly, monotone, min at merge (closed form) | **[PROVED]** + **[NUMERIC]** |
| floor's role clarified (exponent + strip, **not** exposure) | **[PROVED]** (honest correction) |
| boundary $\mathcal B$ bounded for fixed $x_+=O(\bar p_{\min})$ | **[PROVED-modulo-ceiling]** |
| re-entry: total escape = (bounded factor)×bound | **[OPEN — clean uniform excursion count]** |
| cutoff-free turning $\ln$ factor (needs $\bar p''(Y_\ast)$ uniform) | **[PARTIAL — Olver package]** |
| **fully ρ-uniform tail bound $P\le C(\rho)e^{-h^{\star2}/6}$, $\sup_\rho C<\infty$** | **[PROVED up to the re-entry count]** |

**Net.** G1 is **essentially settled**: the prefactor's dominant piece (the exposure) is **provably bounded
uniformly through the merge** — and, contrary to the initial worry, is *smallest* at the merge, by outer
geometry independent of the floor. Combined with the proved rate $c=\tfrac16$ (ρ-independent) and the floor
$\bar p_{\min}\ge c_0$ (exponent), the no-early-escape bound is **uniform in $\rho$** up to one residual
bookkeeping item: the uniform excursion-count for the upper-boundary re-entry (a Tier-3 technicality, not a
new mechanism). The absolute tightness of $C$ (currently $e^{\mathcal E}$, large) is a separate matter — the
true prefactor is polynomial (Berglund–Gentz), recoverable with a sharper martingale.
