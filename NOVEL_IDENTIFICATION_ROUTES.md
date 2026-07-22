# Identifying 𝒲₂ — literature scan + novel routes

_June 2026. A literature-grounded reset on what the cusp escape law 𝒲₂ can and cannot be, and a ranked set
of concrete identification strategies. Tags [ESTABLISHED]/[LITERATURE]/[CONJECTURAL]/[TESTABLE]._

## The decisive literature finding: 𝒲₂ is NOT an (integer) higher-order Tracy–Widom

The multicritical / higher-order TW family (Claeys–Its–Krasovsky; Claeys–Olver, arXiv:1111.3527) — the
$\det(I-K^{(k)}_s)$ from critical unitary ensembles whose edge density vanishes like $(b-x)^{2k+1/2}$ — has
**exactly** these tail exponents (their eqs. 1.22–1.23):
$$\text{left } |s|^{4k+3},\qquad \text{right } s^{(4k+3)/2}\ \Rightarrow\ \text{left/right ratio} = 2\ \text{for every }k.$$
So the family is left-tails 3, 7, 11, … (k=0,1,2). The **cusp ladder** has (Freidlin–Wentzell, this project)
$$\text{left } 2q+1,\quad \text{right } 3q/2\ \Rightarrow\ \text{ratio }(4q+2)/(3q),\ \text{which is 2 only at }q=1.$$
For the cusp (q=2): tails **(5, 3)**, ratio **5/3 ≠ 2**. Since $5\notin\{3,7,11,\dots\}$ **and** the ratio is
wrong, **𝒲₂ is provably not a member of the integer higher-order-TW family** [ESTABLISHED, via literature
tail formulas]. This upgrades the earlier inconclusive numerical attempt to a clean refutation, and explains
*why*: the RMT multicritical edge is a **symmetric soft edge** degenerating (both tails tied to one
exponent), whereas the cusp is an **odd swept potential** $\operatorname{sign}(Y)Y^2$ — intrinsically
asymmetric, so its two tails carry **independent** exponents. That asymmetry is the structural fingerprint.

## The one live "near-miss" — a half-integer multicritical edge (k=½)

The cusp's **left** exponent 5 = 4k+3 at **k=½**. So 𝒲₂ *could* be the analytic continuation of the CIK
family to half-integer k. This is decided by two precise numbers:

| quantity | project / FW value | CIK at k=½ | verdict |
|---|---|---|---|
| left-tail exponent | 5 | 5 | **match** |
| right-tail exponent | 3 (=3q/2) | 5/2 (=(4k+3)/2) | **differ** (3 vs 2.5) |
| left-tail coefficient | 1/20 = 0.050 | $\frac{1}{20}\frac{\Gamma(5/2)^2}{\Gamma(3/2)^2\Gamma(3)^2}=\frac{9}{320}=0.028$ | **differ** |

If the cusp's **right exponent is really 3** and the **left coefficient really 1/20**, 𝒲₂ is genuinely new.
If instead the right exponent is **5/2** and the coefficient **9/320**, then 𝒲₂ **is** the half-integer
(k=½) multicritical edge — a beautiful, citable identification. The current numerics (right-tail fit ≈1.8,
pre-asymptotic) cannot tell 2.5 from 3. **This is the single question to settle**, and it is computable.

## Ranked novel identification routes

**R1 — Instanton / Freidlin–Wentzell optimal-fluctuation tails [TESTABLE NOW, decisive].**
The two tails are the two optimal-escape actions. Solve the FW/instanton (Euler–Lagrange for the
weak-noise action of the swept Riccati) for *early* and *late* escape, extracting **both exponents and
their coefficients exactly**. This settles R-vs-new outright: it tells us whether (right exp, left coeff) =
(3, 1/20) [new] or (5/2, 9/320) [= k=½ multicritical]. Cheap in the sandbox (a deterministic 2-point BVP /
shooting), unlike the rare-event tails that Monte-Carlo cannot reach. **Do this first.**

**R2 — Dynamic-pitchfork weak-noise scaling function [physical home].**
The cusp escape *is* slow passage through a pitchfork (= the cusp normal form). Berglund–Gentz
(arXiv:math/0008208) give the σ/ε^{1/4} concentration and the escape window $[\sqrt\varepsilon,\,
c\sqrt{\varepsilon|\log\sigma|}]$ — the exact scales of this project — but not the law. Maier–Stein
(arXiv:cond-mat/9506097) give the **universal scaling function** for weak-noise escape *at* a bifurcation
(non-Arrhenius, WKB boundary layer at the separatrix). **Identify 𝒲₂ with that scaling function**, extended
from the mean escape time to the full distribution. This is the canonical home for 𝒲₂ and likely yields a
closed scaling form.

**R3 — Painlevé-IV σ-form test [best shot at a closed integrable form].**
The integer PII-hierarchy is now excluded, which *revives* the project's earlier Painlevé-**IV** guess —
and the mechanism supports it: the cusp inner operator is the **parabolic-cylinder (Weber)** operator, and
Weber/Hermite weights are exactly where PIV lives (Freud/Hermite edges, band-touching). **Test the σ-form:**
does $H(s)=\frac{d}{ds}\log\det$ of the cusp law satisfy the PIV σ-equation
$(H'')^2 + 4(H')^3 + \cdots$ (Jimbo–Miwa form)? The project's β-flow test checked a PII-type Riccati; redo it
against the **PIV** σ-form. If it fits, 𝒲₂ has a Painlevé-IV closed form.

**R4 — Asymmetric Riemann–Hilbert problem [the structural requirement].**
The independent left/right exponents (5 ≠ 2×3) mean any RH characterization of 𝒲₂ must have **asymmetric
controlling exponents** — unlike the CIK RH (eq. 2.8) whose single $\theta(\zeta)\sim\zeta^{(4k+3)/2}$ ties
both tails. So 𝒲₂'s parametrix needs a $\theta$ with **different growth on the two rays** (an odd-potential
/ non-self-adjoint Lax pair). This is a concrete structural constraint that narrows the search and is itself
a publishable observation about the catastrophe-ladder family.

**R5 — Half-integer CIK continuation [if R1 says (5/2, 9/320)].**
If R1 finds the right exponent 5/2 and coefficient 9/320, solve the CIK RH problem (Φ, eqs. 2.6–2.8) at
**4k+3 = 5** (contour angle π/5, θ ∼ ζ^{5/2}) numerically (Olver's RH collocation, or our validated
Fredholm machinery with the corrected k=½ kernel) and confirm the full distribution matches 𝒲₂.

**R6 — Integrated-Gaussian first-passage / persistence [long shot, conceptual].**
The left exponent 2q+1 and the swept structure echo first-passage of **higher-order/integrated Gaussian
processes** (random acceleration: persistence θ=1/4; integrated fBm). Relate the late-escape (left) tail to
the persistence of the linearized fluctuation; may explain the 2q+1 exponent combinatorially.

## Recommended program

1. **R1 now** — instanton tails (exponents + coefficients), the decisive 30-line BVP. Settles "genuinely new
   (3, 1/20)" vs "half-integer multicritical (5/2, 9/320)."
2. Depending on R1: either **R3** (PIV σ-form, for a closed form) or **R5** (k=½ RH, confirm the multicritical
   identity).
3. **R2** in parallel — the dynamic-pitchfork scaling function as the physical identification.

The honest headline: the elegant "higher-order TW" reframe is **refuted in its integer form** by the tail
ratio (2 vs 5/3), but it leaves a sharp, testable residue — a half-integer (k=½) near-miss decided by one
coefficient — and it redirects the search toward **Painlevé-IV / the dynamic-pitchfork scaling function**,
both of which are concretely testable next.

## Sources
- [Claeys & Olver, *Numerical study of higher order analogues of the Tracy–Widom distribution*, arXiv:1111.3527](https://arxiv.org/abs/1111.3527)
- [Claeys, Its, Krasovsky, *Higher order analogues of the TW distribution and the Painlevé II hierarchy*, arXiv:0901.2473](https://arxiv.org/abs/0901.2473)
- [Berglund & Gentz, *Pathwise description of dynamic pitchfork bifurcations with additive noise*, arXiv:math/0008208](https://arxiv.org/abs/math/0008208)
- [Maier & Stein, *A scaling theory of bifurcations in the symmetric weak-noise escape problem*, arXiv:cond-mat/9506097](https://arxiv.org/abs/cond-mat/9506097)
- [*Reaction rates and the noisy saddle-node bifurcation: RG for barrier crossing*, arXiv:1902.07382](https://arxiv.org/abs/1902.07382)
- [*Record statistics of integrated random walks and the random acceleration process*, arXiv:2109.05582](https://arxiv.org/abs/2109.05582)
- [Betea, Bouttier, Walsh, *Multicritical Schur measures and higher-order analogues of the TW distribution*, arXiv:2307.05303](https://arxiv.org/abs/2307.05303)
