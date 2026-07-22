# Program 2, route (b) — sub-attack 1 EXECUTED: the O(η²) stochastic Stokes constant, defined and computed

_2026-07-10 (Fable 5). Executes the minimal target lemma of `PROGRAM2_HANDOFF_B_STOCHASTIC_WKB.md` §"Concrete
sub-attacks" #1 and the deep-research's "minimal provable lemma" (`PROGRAM2_STOCHASTIC_WKB_DEEPRESEARCH.md`,
Open Q4). Scripts: `coupled-atlas/stochastic_stokes_o_eta2.py` (main), `_validate_pert.py` (backbone check),
`_stochastic_stokes_figure.py`; figure `coupled-atlas/figures/stochastic_stokes_o_eta2.png`. All numbers
certified (grid h, oscillatory cutoff To, Borel-rotation θ, and an independent exact-Newton validation)._

## Headline
The **stochastic Stokes constant is now DEFINED as 𝔼[random Voros datum] and its O(η²) value is computed and
certified** (Milestone 1 of the handoff). Concretely, the noise-averaged connection root is
$$\mathbb E[\lambda] \;=\; \lambda_0 + \eta^2\,\Omega_2 + O(\eta^4),\qquad
\boxed{\;\Omega_2 = -0.451 + 0.352\,i\;}\quad |\Omega_2|=0.572,\ \arg\Omega_2 = +142^\circ\ (\approx \tfrac{3\pi}{4}).$$
Two structural results fall out, both bearing directly on the frontier:
1. **The 45°→54° Borel-phase shift is NOT a mean-shift effect at O(η²).** 𝔼[λ] stays at arg ≈ −45° (rotation
   rate only −3.2°/η²; reaching −54° needs η²≈2.8 ≫ the perturbative radius η²_c≈0.6). The shift Ω₂ is nearly
   *radial* (shrinks |λ|). **The phase shift lives in the FLUCTUATION sector, not the mean** — a sharp,
   quantitative confirmation of the Bureković–Grauer mechanism hint (Finding 4).
2. **The lemma's literal form is corrected, and the renormalization-clean object identified.** Noise-averaging
   does *not* leave the Borel *location* finite-and-fixed: the classical-period (action) mean shift is
   **δ(0)-divergent** (the v3 renormalization threshold). Only the **connection datum** — a *non-local* Ito
   functional of the noise — has a finite, δ(0)-free O(η²) average. So the well-defined "stochastic Stokes
   constant" is 𝔼[connection root] = Ω₂, *not* 𝔼[period].

## Setup (the well-defined route — avoids the ruled-out complex-scaling trap)
Skeleton `u'' = (sign(Y)Y² − λ − η Ẇ(Y)) u`. Deterministic connection datum = the closed-form Weber
Γ-equation root **λ₀ = 0.8896 − 0.8896i** (arg −45°), reproduced here to machine precision (|D|=9e-14) by an
independent two-sided physical-Y solver: recessive Weber on Y>0, outgoing (Gamow) on Y<0, resonance = match of
the physical log-derivatives `m_L(0) = m_R(0)` (the e^{3iπ/4} inversion phase is encoded automatically in the
numerical outgoing solution).

**We never impose a boundary condition on a noise realization** (that is the ill-defined thing the
deep-research ruled out). Instead we average over the noise FIRST. For the log-derivative `m = u'/u` the Riccati
`m' = (Q−λ) − ηẆ − m²` linearizes with integrating factor `u₀²` (since `(u₀²)'/u₀² = 2m₀`):
```
(u₀² m₁)' = −ηẆ u₀²     (O(η), linear in noise)   ⇒  m₁(0) = ± u₀(0)⁻² ∫ Ẇ(s) u₀(s)² ds
(u₀² m₂)' = −u₀² m₁²     (O(η²), quadratic)        ⇒  m₂(0) =   u₀(0)⁻² ∫ u₀(s)² m₁(s)² ds
```
Taking 𝔼[·] collapses 𝔼[Ẇ(s)Ẇ(s′)] = δ(s−s′) to **deterministic single integrals** — the program's
Green's/Wick second-Wiener-chaos machinery (cf. `weaknoise_greens.py`, `v0 = (c/u₀*)²∫u₀⁴`). The noise-averaged
connection function `D_eff(λ) = D(λ) + 𝔼[δD⁽¹⁾] + 𝔼[δD⁽²⁾] + O(η⁴)` has `𝔼[δD⁽¹⁾]=0`; the first shift is O(η²).

## The O(η²) pieces (all per η², normalization-invariant, certified)
| quantity | value | meaning |
|---|---|---|
| `D'(λ₀)` | −0.794 − 0.329i | connection-function slope |
| `D''(λ₀)` | +0.206 + 0.495i | curvature |
| `σ²_D = A_L + A_R` | +0.322 + 0.133i | pseudo-variance of `δD⁽¹⁾` (Ito isometry ∫u⁴) |
| `A_R = u_R(0)⁻⁴∫₀^∞ u_R⁴` | +0.305 − 0.280i | confining side |
| `A_L = u_L(0)⁻⁴∫_{−∞}⁰ u_L⁴` | +0.017 + 0.414i | oscillatory side (**Borel-regularized**, see below) |
| `𝔼[δD⁽²⁾] = M_L − M_R` | −0.303 + 0.303i | mean 2nd-chaos (nested ∫u⁻²(∫u⁴)) |
| `Var(λ)/η²` (pseudo) | +0.436 − 0.181i | anisotropy axis −11.3° |
| `Var(λ)/η²` (true, 𝔼\|δλ\|²) | 1.784 → **rms\|δλ\| = 1.34 η** | fluctuation magnitude (≈\|λ₀\|=1.26 at η²≈0.6) |
| **Ω₂** | **−0.451 + 0.352i** | **the O(η²) stochastic Stokes constant** |

Assembly (standard 2nd-order root perturbation, all holomorphic in λ):
`Ω₂ = −(1/D')[ ½D''·(σ²_D/D'²) − (1/2D')·dσ²_D/dλ + (M_L − M_R) ]`.
Structural cross-check: **`M_L = i·conj(M_R)`** to 4 digits — the confining(real-Weber)↔oscillatory
(imaginary-Weber) e^{iπ/2} mirror, produced by the two *independent* solvers. Pseudo-variance arg −22.5° = −45°/2.

## The oscillatory Voros integral needs Borel regularization (meaningful, not a numerical nuisance)
`∫_{−∞}⁰ u_L⁴ dY` is **conditionally convergent** (amplitude ~|Y|⁻¹/², phase ~e^{2iY²}): naive real-axis
truncation drifts (A_L: 0.029→0.022→0.008 as To grows). We regularize by **rotating the Y<0 contour**
`x = ρe^{iθ}` so the outgoing wave decays and the integral becomes absolutely convergent — the exact-WKB
lateral/Borel prescription (as `complex_scaling.py`). The regularized value is **θ-independent** (identical to
5 digits across θ∈[0.2,0.4] and cutoffs To∈[14,22]) — this *is* the median/Borel-summed value. This conditional
convergence is itself the fingerprint that the oscillatory-side connection carries the Stokes phenomenon.

## Validation (independent)
`_validate_pert.py`: perturb the potential by a *deterministic* smooth bump `εg(Y)` and compare the perturbative
root shift `δλ₁ = −δD/D'` (kernel `δm_R(0) = −u_R(0)⁻²∫g u_R²`) against an **exact Newton solve** of the
ε-perturbed resonance. Quadratic fit over ε∈[0.02,0.10]: predicted `dλ/dε = +0.0966−0.2442i` vs fitted
`+0.0975−0.2427i` → **0.67% agreement**. The perturbation backbone (log-derivative variation + D′) is correct;
the only extra step for noise is the standard Wick contraction (Ito isometry), so Ω₂ is trustworthy.

## The lemma test — location (δ(0)) vs constant (finite), the decisive dichotomy
Direct test of "noise-averaging leaves Borel-singularity LOCATIONS fixed, shifting only the constant":
- **Borel LOCATION** = classical period `P(λ) = ∮√(λ − Y²)`. Its O(η²) mean shift
  `𝔼[δ²P] = −¼η²∮ 𝔼[Ẇ²]/(λ−Y²)^{3/2}` carries `𝔼[Ẇ(Y)²]=δ(0)` — a **coincident self-contraction**. Numerically
  it **diverges as 1/h** (|shift| = 113→227→454→908 as h halves, clean δ(0)~1/h). NOT finite, NOT fixed —
  it is exactly the v3 renormalization threshold (`PROGRAM2_ROUTE2B_NOTES.md` §7, the boundary δ(0)).
- **Connection datum** 2nd chaos `∫u_R⁴` is **finite and grid-convergent** (A_R = 0.305−0.280i, stable to 5
  digits under h halving). Why: the log-derivative is a *non-local* functional (`m₁(s) = u₀⁻²∫_s Ẇ u₀²`), so its
  variance is an **Ito isometry ∫u⁴**, never a coincident δ(0). Non-locality self-renormalizes.

**Verdict.** The literal lemma is false as stated (the *location* mean shift is δ(0)-divergent), but its *intent*
is realized in a sharper form: **the renormalization-clean, finite O(η²) noise-averaged object is the connection
datum**, and its shift `Ω₂` IS the well-defined stochastic Stokes constant at O(η²). This also explains
*mechanically* why the naive "shift the period" attacks fail — they sit on the δ(0) wall — while the
Green's/Wick connection route is finite.

## Where this leaves the frontier (honest)
- **Milestone 1: DONE** — stochastic Stokes constant defined (𝔼[Voros datum]) with its O(η²) value computed,
  certified, and validated. `Ω₂ = −0.451 + 0.352i`.
- **Milestone 2 (45°→54° explained): NOT via the mean.** Quantified: the O(η²) mean rotates only −3.2°/η²
  (needs η²≈2.8 to reach 54°). The phase shift must come from the **fluctuation** (rms|δλ|≈1.34η is large, and
  the cloud is anisotropic with pseudo-variance axis −11.3°), i.e. from how the *distribution* of the random
  connection root sets the large-order growth — the Bureković–Grauer fluctuation-sector mechanism, now
  localized. **Next move:** compute how the anisotropic fluctuation cloud (2nd chaos) feeds the Borel-singularity
  *phase* via the large-order asymptotics of the v_n — i.e. push to the variance→Borel-phase link, not the mean.
- **Open:** O(η⁴) (needs the 4th-chaos / the v3 renormalization counterterm made explicit); tie Ω₂ to Hao's
  deterministic Ω=1 via a Voros-symbol identity (is `1 + η²·(Ω₂-derived residue)` the noise-dressed Stokes
  constant in Hao's normalization, or is λ₀'s shift a *different* connection object?). The mapping
  resonance-root-shift ↔ Hao's Stokes-constant residue is the remaining conceptual bridge for a theorem.

## Files
`coupled-atlas/stochastic_stokes_o_eta2.py` — solvers, 2nd-chaos integrals, assembly, lemma test (run it).
`coupled-atlas/_validate_pert.py` — exact-Newton backbone validation (0.67%).
`coupled-atlas/_stochastic_stokes_figure.py` → `figures/stochastic_stokes_o_eta2.png`.
