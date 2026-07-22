# Notes — the cusp rung at its honest ceiling (rung B → ~90%)

_June 2026. Figure `coupled-atlas/figures/cusp_rung_capstone.png`; script `cusp_rung_capstone.py`._

Same three-level treatment as the fold rung — but the cusp rung has a genuinely **open** third
level, so it cannot honestly reach 100%. Here is exactly what is established and what is not.

| level | object | status |
|---|---|---|
| marginal | Weber / parabolic-cylinder (≠ Pearcey) | **[PROVED reduction]** (K-P: slow-fast cusp → 2nd-order Weber, not 3rd-order Pearcey) + **[NUMERIC]** KP coupled-FHN escape ↔ Weber, KS 0.056 vs Pearcey 0.123; inner Weber skew +0.64, distinct from TW (+0.22) and Pearcey (+0.19) |
| point process | Weber point process | **[NUMERIC, swept]** successive swept-Weber nodes (means −1.72, −3.06, −3.96); swept-node corr(Y¹,Y²)=+0.92. Genuine multicritical rigidity: **open** |
| time process | Weber analogue of Airy₂ | **[NUMERIC, surrogate]** forced cusp peel-offs = stationary Weber-marginal process (decaying covariance, increments→2·Var). **Intrinsic** Dyson–Weber object: **open** |

## The one correction (parallels the fold)

The cusp **swept-node** correlation is **+0.92** — the Weber analogue of the fold's swept-node +0.89,
and slightly larger (the Weber funnel is more confining). As with the fold, this is the *swept-node*
quantity: successive nodes of one swept solution share the noise path, so they over-correlate relative
to the genuine point-process rigidity. The **genuine multicritical point-process rigidity** (the cusp
analogue of the fold's 0.50 from the GUE edge) requires the Dyson–Weber object and is **open**.

## The no-go that caps the rung

The fold rung hit 100% because its intrinsic process is the **Airy line ensemble**, realised exactly by
the matrix-OU DBM edge (`dyson_weber.py`). For the cusp this route **provably does not work**: the
GUE/matrix-OU edge top eigenvalue is Airy (skew +0.25 ≈ TW₂), whereas the cusp marginal is Weber
(skew +0.64). So the cusp's intrinsic process is **not** the Airy line ensemble — it is the
**multicritical-edge "Dyson–Weber" line ensemble** (DBM at a higher-order edge, equilibrium density
vanishing like dist^{3/2} rather than √). That object is genuinely open:
- no RRV-type stochastic-operator characterization is known for k>1 (the naive harmonic operator gives
  Gaussian, not Weber — `spectral_id.py`);
- matrix-OU realises only the Gaussian/Airy edge; reaching the multicritical edge needs a
  critically-tuned non-Gaussian matrix model, which is numerically delicate (`dyson_weber.py` notes).

## What a real 100% for B would need (path)

1. The **intrinsic Dyson–Weber ensemble**: construct the multicritical-edge DBM (tuned non-Gaussian
   log-gas at criticality) and confirm its marginal = Weber and its line-ensemble covariance — the
   genuine counterpart of the matrix-OU Airy₂ for the fold.
2. The **genuine multicritical point-process rigidity** (replacing the swept-node +0.92), read off the
   same ensemble — the cusp analogue of the fold's 0.50.
3. (Stretch) an **RRV-type stochastic-operator** for the Weber/multicritical edge — a theorem-level
   identification of the marginal, the cusp analogue of stochastic-Airy = TW.

Items 1–3 are the open core shared with rung F. Until at least (1)+(2) land, B stays at ~90%, not 100.
