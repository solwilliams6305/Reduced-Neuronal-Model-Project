# Smooth 𝒲₂ via the escape Fokker–Planck PDE — and why the single-law σ-form test is confounded

_June 2026. Figure `coupled-atlas/figures/fp_smooth_W2.png`; scripts `fp_cusp.py`, `fp_figure.py`.
The follow-up to R3: get a smooth 𝒲₂ (no MC noise) to "unlock" the σ-form test. The FP succeeded; the
σ-form test turned out to be confounded for a single distribution. Tags
[VALIDATED]/[NUMERIC]/[CONFOUND]/[NEXT]._

## Win: a smooth, validated 𝒲₂

Solved the escape Fokker–Planck PDE for the Riccati density
$\partial_\tau\rho=-\partial_p[(V-p^2)\rho]+\tfrac{\eta^2}{2}\partial_{pp}\rho$, $V=\operatorname{sign}(Y)Y^2$,
conservative finite-volume (upwind advection + central diffusion), **absorbing** at $p=p_{\min}$ (escape),
**reflecting** at $p=p_{\max}$ (no rightward escape). Survival $S(\tau)=\int\rho\,dp$ gives the
escape-location CDF $F(y)=S(Y_0-y)$ — **smooth, Monte-Carlo-free**.

**Validation against MC (1.4M):** cumulants match to ~2%:

| | mean | std | skew | exk | κ₅ | κ₆ |
|---|---|---|---|---|---|---|
| FP (smooth) | −1.594 | 0.689 | **+0.601** | **−0.244** | **−2.10** | **−2.60** |
| MC target | −1.628 | 0.684 | +0.61 | −0.24 | −2.2 | −2.7 |

The smooth density overlays the MC histogram in both the bulk and the (log-scale) tails. **[VALIDATED]**
A clean representation of the new law is now in hand (saved `fp_cusp_F.npy`).

## The σ-form test is confounded for a single law [important methodological finding]

The plan was: with smooth $F$, test whether $\sigma=(\log F)'$ obeys a Painlevé σ-form
$(\sigma'')^2=\text{poly}(s,\sigma,\sigma')$. **It does — but so does everything**, vacuously:

| law (smooth) | "1st-order Riccati" $R^2$ | "2nd-order σ-form" $R^2$ |
|---|---|---|
| cusp 𝒲₂ (FP) | 1.0000 | 1.0000 |
| TW₂ (Painlevé control) | 1.0000 | 1.0000 |
| Gaussian (null) | 1.0000 | 1.0000 |

A single distribution is a **1-parameter solution curve** $s\mapsto(\sigma,\sigma',\sigma'')$. Any three
smooth functions of one variable satisfy *many* polynomial relations locally, so the regression fits with
$R^2=1$ for the cusp, for TW₂, **and for a Gaussian** alike. The earlier "CV-R² 0.94 passes" was really
detecting **smoothness vs MC-noise**, not Painlevé structure (the Gaussian only "failed" before because its
*Monte-Carlo* $\sigma''$ was noisy; a *smooth* Gaussian passes too — it even obeys a 1st-order Riccati
$\sigma'=-\sigma^2-z\sigma$). **So the single-law σ-form regression has no power to identify PIV.** [CONFOUND]

The σ-form's real content is that **one fixed equation holds across a whole family** (with only the
monodromy parameters $\theta(\beta)$ varying), or that the **exact coefficients** match in the natural
variable. Neither is testable from one law.

## What still stands, and the right next test

- **Robust identity facts (R1, unchanged):** 𝒲₂ has tail exponents **(5, 3)** (ratio 5/3), placing it
  **outside the standard soft-edge Painlevé family** (all ratio 2). Genuinely new.
- **The correct, non-confounded σ-form test is the β-family consistency test:** compute $\mathcal W_\beta$
  for several $\beta$ (cheap now — the FP solves in seconds per $\beta$), and check whether a **single**
  σ-form — fixed $(s\sigma'-\sigma)^2$ and $\sigma'^3$ structure, only $\theta(\beta)$ varying in the
  $\sigma'$-slots — fits **all** $\beta$ jointly (multiple curves pinning one surface). The open
  difficulty there is identifying the **natural variable** (the β-dependent rescaling of $Y^\star$), since
  the σ-form is not the same in the standardized variable across $\beta$.

## Status

| item | status |
|---|---|
| smooth 𝒲₂ via escape FP PDE, validated vs MC (~2%) | **[VALIDATED]** — `fp_cusp_F.npy` |
| single-law σ-form regression decides PIV | **[CONFOUND]** — vacuous ($R^2$=1 for TW₂/Gaussian/cusp) |
| 𝒲₂ ∉ standard soft-edge Painlevé family (tails 5,3) | **[ROBUST]** (R1) |
| 𝒲₂ = PIV (or any Painlevé) | **[OPEN]** — needs the β-family / exact-coefficient test |

**Net.** The Fokker–Planck route delivered its concrete promise — a **smooth, validated $\mathcal W_2$** — but
also exposed that a single distribution cannot be tested for a σ-form by regression (the 1-D-curve
confound; even TW₂ and a Gaussian "pass"). The honest path to deciding PIV is the **β-family consistency
test** (a fixed σ-form across $\mathcal W_\beta$), which the cheap FP solver now makes feasible, modulo
identifying the natural variable. Meanwhile the rigorous tail facts keep 𝒲₂ outside the standard
soft-edge Painlevé family.
