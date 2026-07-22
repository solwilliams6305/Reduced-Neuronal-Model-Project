# R3 — Painlevé-IV σ-form test for 𝒲₂: inconclusive (MC-noise-limited), plus structural context

_June 2026. Script `sigma_form_test.py`. Goal: test whether the cusp law's
$\sigma(s)=\frac{d}{ds}\log F(s)$ satisfies a Painlevé σ-form (the PIV candidate), since a Painlevé σ-form
is a 2nd-order 2nd-degree ODE $(\sigma'')^2=\text{poly}(s,\sigma,\sigma')$. Tags
[VALIDATED]/[NUMERIC]/[INCONCLUSIVE]/[STRUCTURAL]._

## Method and controls

$\sigma$-form is covariant under affine $s\mapsto as+b$, so all laws are standardized. Fit $\log F$ by a
polynomial over the bulk, take $\sigma,\sigma',\sigma''$ analytically, then **sparse + cross-validated**
regression of $(\sigma'')^2$ onto a monomial basis spanning the PII and PIV σ-forms (PIV's distinguishing
block is $(s\sigma'-\sigma)^2$). Sparse (best size-5 subset) defeats overfitting; CV (fit left half of the
range, predict right) tests whether the relation is *global* (a true σ-form) rather than an accidental fit.

| law | CV-R² | meaning |
|---|---|---|
| **TW₂ smooth** (Fredholm, no noise) | **+0.994** | calibration ✓ — method recovers the PII σ-form |
| **TW₂ sampled** (1.4M draws, same pipeline) | **−3.6** | a *true Painlevé* + matched MC noise — **FAILS** |
| Gaussian (null) | −36 | non-Painlevé null — fails (as expected) |
| **cusp 𝒲₂** (1.4M MC) | −0.06 | **not interpretable** (baseline fails) |

## Verdict: inconclusive — and why

The method **works on smooth data** (TW₂-smooth CV-R² 0.99). But the matched-noise control **TW₂-sampled
fails** (−3.6): **1.4×10⁶ Monte-Carlo samples are not enough to detect a σ-form** through the 3rd-derivative
($\sigma''$) pipeline — the sampling noise in $\sigma''$ destroys the signal. Since a *known* Painlevé law
fails the test at this data quality, the cusp's score cannot be interpreted (its −0.06 happening to beat
TW₂-sampled's −3.6 is not meaningful when the baseline itself fails).

**A σ-form test on 𝒲₂ requires a *smooth* representation of $F(s)$** — from the escape **Fokker–Planck PDE**
or from the **R2 dynamic-pitchfork scaling function** — not Monte-Carlo. That is the prerequisite, and it is
close to the open problem itself (a smooth/integrable representation of 𝒲₂ is most of what we are after).

## Structural context (what we *can* say)

- The standard **soft-edge Painlevé family** (Airy/PII and the multicritical PII-hierarchy) all have tail
  ratio **left:right = 2** (exponents $4k{+}3$, $(4k{+}3)/2$). The cusp's tails are **(5, 3), ratio 5/3** —
  **not in that family** (consistent with R1).
- **But** "escape-driven" does **not** preclude Painlevé: the fold law (TW) is itself a Riccati-explosion
  escape law (RRV) *and* a PII σ-form. So the cusp could still secretly be a Painlevé — just not a
  standard soft-edge one. A Painlevé with ratio-5/3 tails would be non-standard (possibly an **asymmetric
  isomonodromy**, R4 — different exponents on the two Stokes rays, matching the cusp's independent left/right
  tails).
- The **Weber/PIV resemblance** is genuine at the level of the **deterministic inner operator** (the cusp
  skeleton is the parabolic-cylinder operator, whose Painlevé is PIV). Whether that lifts to a PIV identity
  for the *noisy* edge law 𝒲₂ is exactly what the σ-form test would decide — and cannot, yet.

## Status

| item | status |
|---|---|
| σ-form sparse+CV pipeline; PII recovered from smooth TW₂ | **VALIDATED** (CV-R² 0.99) |
| σ-form test of 𝒲₂ from 1.4M MC | **INCONCLUSIVE** (noise-limited; Painlevé baseline fails) |
| 𝒲₂ ∉ standard soft-edge Painlevé family (ratio 5/3≠2) | **STRUCTURAL** (from R1) |
| 𝒲₂ = PIV (or any Painlevé) specifically | **UNDECIDED** — needs smooth $F$ |
| prerequisite: smooth $F(s)$ via escape Fokker–Planck PDE / R2 | **NEXT STEP** |

**Net.** R3 did not decide PIV — honestly, because detecting a σ-form needs a smooth $F(s)$ that
Monte-Carlo can't provide (demonstrated by the TW₂-sampled control failing). The pipeline is validated and
ready. The productive order is now: **(i)** get a smooth $\mathcal W_2$ from the escape Fokker–Planck PDE
(or R2's dynamic-pitchfork scaling function), **(ii)** then the σ-form / PIV test becomes decisive. The
structural facts meanwhile place 𝒲₂ outside the standard soft-edge Painlevé family while leaving a
non-standard (asymmetric) isomonodromy open.
