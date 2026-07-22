# β-family consistency test — verdict: the σ-form regression is confounded; PIV undecided by numerics

_June 2026. Figure `coupled-atlas/figures/fp_beta_family.png`; scripts `fp_beta.py`, `fp_beta_figure.py`.
The decisive non-confounded test attempt for Painlevé-IV. Tags [NUMERIC]/[CONFOUND]/[ROBUST]/[VERDICT]._

## What was built

Solved $\mathcal W_\beta$ via the escape Fokker–Planck PDE for **β = 1, 2, 4, 8** (η = 2/√β). A genuine
family of **distinct shapes**:

| β | η | mean | std | skew | exk |
|---|---|---|---|---|---|
| 1 | 2.00 | −1.243 | 0.904 | +0.258 | −0.78 |
| 2 | 1.41 | −1.593 | 0.690 | +0.600 | −0.24 |
| 4 | 1.00 | −1.817 | 0.483 | +0.843 | +0.64 |
| 8 | 0.71 | −1.943 | 0.321 | +0.843 | +1.18 |

Skew rises 0.26 → 0.84 (saturating toward the bare law); excess kurtosis −0.78 → +1.18. **[NUMERIC]**

## Natural variable — reasoned, then corrected by the data

The σ-form's $(s\sigma'-\sigma)^2$ block is **not shift-covariant**, so the right **origin** and **scale**
matter. Origin: the deterministic edge $Y^\star_{\rm det}=-2.09$ (the η→0 spectral edge). Scale: the
operator balance $u''=(\operatorname{sign}(Y)Y^2-\eta\dot W)u$ with $Y=\eta^a\hat Y$ gives potential
$\sim\eta^{4a}$, noise $\sim\eta^{1+3a/2}$, equal at $a=2/5$. **But the data corrects this:** the
escape-location **width scales as $\eta^{1.00}$** (measured), not $\eta^{2/5}$. The $\eta^{2/5}$ is the
*inner spectral* scale; the *escape-location* (bulk) width is linear in η — the fluctuation around the
canard is OU-tube-like ($\propto\eta$). **[NUMERIC — honest correction to the a-priori reasoning.]**

## The test is confounded — and why

Joint shared-structure fit $(\sigma'')^2 = A(s\sigma'-\sigma)^2 + B\sigma'^3 + [D_\beta\sigma'^2+E_\beta\sigma']$,
with **A,B shared** across β (the rigid PIV signature) and $D_\beta,E_\beta$ the per-β monodromy $\theta(\beta)$:

| scale $\eta^p$ | cusp joint $R^2$ | skew-normal NULL $R^2$ |
|---|---|---|
| p=0.30 | 0.9998 | 0.9999 |
| p=0.40 | 0.9997 | 0.9998 |
| p=0.50 | 0.9996 | 0.9998 |
| p=0.667 | 0.9995 | 0.9997 |

**The null fits as well as the cusp ($R^2\approx0.9997$ both) — the test has no discriminating power.**
[CONFOUND] The reason is now clear and instructive:
1. PIV genuinely *needs* per-β $\theta(\beta)$ freedom (the shapes vary). But that per-β freedom (2 params
   each) is exactly enough to refit *any* smooth family — re-opening the 1-curve confound at the family level.
2. **Deeper:** the σ-form regression operates on the **bulk** ($F\in[0.05,0.95]$), and there $\mathcal W_\beta$
   is a **generic skewed-Gaussian-like core** (width $\propto\eta$ — that is *why* a skew-normal null matches
   it). The Painlevé-distinguishing structure lives in the **tails** (the (5,3) exponents) and the
   **connection/Stokes data** — which a bulk regression cannot see, and which MC/FP cannot resolve far enough.

So both the single-law test (1-curve confound) and the β-family test (per-β-θ confound, bulk-generic) are
**confounded**. Regression on the distribution cannot decide PIV.

## Verdict

- **PIV is UNDECIDED by numerics.** The σ-form regression — single-law or β-family — is confounded; it
  cannot distinguish $\mathcal W_\beta$ from a generic skewed family. I will **not** claim PIV confirmed on a
  test that a skew-normal also passes. **[VERDICT]**
- **Robust facts stand:** $\mathcal W_2$ has tail exponents **(5, 3)** (R1, analytic), ratio 5/3 — **outside
  the standard soft-edge Painlevé family** (all ratio 2). So it is *not* the standard PII-hierarchy; whether
  it carries a PIV or a new (asymmetric) isomonodromy is **not decidable by bulk regression**. **[ROBUST]**
- **What would decide it:** (i) an **analytic** derivation of the σ-form from the stochastic-Weber operator;
  or (ii) the **tail constants** (the e^{−|s|⁵/20} coefficient and the right-tail constant) matched against a
  candidate transcendent; or (iii) the **connection/Stokes data** of the parabolic-cylinder parametrix
  (the §0 keystone object). None is a bulk-regression task.

## Rigor ledger (the tube program)

| item | status |
|---|---|
| **re-entry count** (boundary term of the no-early-escape bound) | **DONE [PROVED]** — the linearly-truncated potential gives a *global* supermartingale on $\{x\ge -L\}$ (the $-x^2$ confines from above, no upper boundary), so the one-sided bound is **fully ρ-uniform pre-turning, no gaps** (`TWO_SIDED_TUBE_AND_CUSP_LAW.md`) |
| **two-sided tube (C1)** | **PROVED on the approach window** $[0,\tau_\ell]$ — lower side rate $h^{\star2}/6$ (cubic supermartingale), upper side $h^2/2$ (OU comparison), constants uniform in ρ |
| **post-turning $O(\ell)$ core** (where C1 hands to the edge law C2) | **OPEN / deferred** — the connection region; needs the uniform parabolic-cylinder-with-noise parametrix |

## Net

The Fokker–Planck β-family delivered a clean family of smooth $\mathcal W_\beta$ (skew 0.26→0.84, width
$\propto\eta$), but **settled that the σ-form regression cannot decide PIV** — the skew-normal null passes
the same test (R²≈0.9997 both), because the bulk is generically skewed and the Painlevé signature lives in
the tails/connection data. Honest verdict: **PIV undecided by numerics; $\mathcal W_2$ is established as
genuinely new (tails 5,3, outside the soft-edge family); its integrable identity (PIV or a new asymmetric
isomonodromy) requires analytic derivation or the connection/Stokes data, not distribution-fitting.** The
tube program's bookkeeping is closed where claimed: re-entry **done**, C1 **proved on the approach**,
post-turning core **open**.
