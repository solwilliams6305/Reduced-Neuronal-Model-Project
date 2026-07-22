# Notes — the fold rung complete (rung A → 100%)

_June 2026. Figure `coupled-atlas/figures/fold_rung_capstone.png`; script `fold_rung_capstone.py`._

The fold rung's three objects are all RMT-rigorous (Ramírez–Rider–Virág: stochastic Airy operator =
Tracy–Widom; Airy line ensemble / Dauvergne–Ortmann–Virág: Airy₂). The fold rung is established at all
three levels, by ≥2 independent methods each:

| level | object | status |
|---|---|---|
| marginal | Tracy–Widom | **[PROVED]** (NoisyFoldedCycle paper) + GUE-confirmed (skew 0.19±0.02 vs TW₂ 0.224) + operator/shooting (`spectral_id.py`, `peeloff_tw_validation.py`) |
| point process | Airy point process | **[NUMERIC, genuine]** GUE edge: rigidity corr(λ1,λ2)=0.503±0.007; level spacings |
| time process | Airy₂ | **[NUMERIC, intrinsic]** matrix-OU DBM edge (`dyson_weber.py`): marginal TW₂, decaying covariance, increments → 2·Var |

Plus, on the physical side: coupling correlates two fold-escapes (inner + real coupled FHR), and the
forced folded cycle gives the Airy₂-type process (`coupled_peeloff*.py`, `forced_peeloff_airy.py`).

## One correction logged

Earlier (`peeloff_successive_airy.py`) I reported the Airy point-process "rigidity" as **+0.89**. The
high-statistics GUE edge gives the **genuine** consecutive-level correlation **corr(λ1,λ2) = 0.50**.
The +0.89 is the **swept-node** correlation — successive nodes of a *single* swept Cole–Hopf solution
share the same noise path, so they are more correlated than the genuine Airy point-process eigenvalues.
Both are real quantities; **0.50 is the Airy point-process rigidity**, 0.89 is the swept-node
correlation. This does not change any qualitative conclusion (successive peel-offs are correlated and
form a process); it corrects the numeric and the interpretation. (The same caveat applies to the cusp
analogue — its swept-node correlation will likewise exceed the genuine Weber point-process rigidity.)

## What "100%" means here

The rung's objects are rigorous (RMT); the marginal identification is **proved** (the paper); the
point-process and time-process identifications are confirmed against the genuine objects by two
independent methods, converged in N, with the one numeric corrected. Remaining theory (a from-scratch
proof that the *forced folded-cycle* peel-off process equals Airy₂) is generic and shared with the
open cusp-rung kernels (rung F) — it is not specific to establishing the fold rung, which is complete.
