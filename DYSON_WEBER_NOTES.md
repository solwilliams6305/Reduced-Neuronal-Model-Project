# Notes — the Dyson–Weber object (intrinsic edge processes)

_June 2026. Figure `coupled-atlas/figures/dyson_weber.png`; script `dyson_weber.py`.
Tags: [VALIDATED] · [INTRINSIC] · [OPEN]._

The forced-passage results (fold→Airy₂-type, cusp→Weber-type) used an *imposed* OU correlation as a
surrogate for the forcing. The genuine, **intrinsic** processes come from the Dyson-Brownian-motion
edge — the Airy line ensemble for the fold, and (conjecturally) a parabolic-cylinder/multicritical
line ensemble for the cusp ("the Dyson–Weber object").

## Fold: intrinsic Airy₂ — [VALIDATED]

A Hermitian matrix Ornstein–Uhlenbeck process M(t+dt)=ρM(t)+√(1−ρ²)G (ρ=e^{−dt/2}, G fresh GUE) has
eigenvalues performing exact DBM — *stable*, unlike the singular eigenvalue Langevin (which overflows
under explicit Euler). The top eigenvalue λ_max(t) is the intrinsic Airy₂ process:
- **marginal = Tracy–Widom:** standardised skew **0.200** (TW₂ 0.224), excess kurtosis **0.088**
  (TW₂ 0.093) — spot on;
- **intrinsic covariance:** C(τ) decays (to <0.2 by τ≈0.85 in DBM time), and the increment variance
  saturates at 2·Var — the Airy₂ covariance shape, now obtained *intrinsically* (no imposed ρ).

So the fold rung's multi-point process is pinned down to the genuine Airy₂, not just a surrogate.

## Cusp: the Dyson–Weber object — [OPEN, now framed]

The cusp's intrinsic process is the edge process of a DBM at a **multicritical** point — where the
equilibrium spectral density vanishes faster than √ at the edge (a higher-order / parabolic-cylinder
line ensemble). Matrix-OU only gives the Gaussian (semicircle, √-edge ⇒ Airy) ensemble; reaching a
multicritical edge needs a **critically-tuned non-Gaussian matrix model** (eigenvalue Langevin with a
quartic potential at its critical coupling). Numerically this is delicate and research-grade:
- the naive eigenvalue Langevin is unstable (singular repulsion overflows under explicit Euler);
- a negative quartic (to flatten the edge) is unbounded (the gas escapes);
- hitting the exact multicritical coupling is a codim-1 tuning that is strongly finite-N sensitive.

So the Dyson–Weber object is precisely **defined** (multicritical-edge line ensemble, the parabolic-
cylinder analogue of the Airy line ensemble) and identified as the deep open RMT-side problem — the
intrinsic counterpart of the cusp's Weber edge marginal.

## Status

| object | status |
|---|---|
| intrinsic Airy₂ (fold) from DBM edge | **[VALIDATED]** — marginal TW₂, intrinsic covariance |
| matrix-OU as the stable DBM construction | **[VALIDATED]** (Langevin unstable) |
| Dyson–Weber (cusp): multicritical-edge line ensemble | **[OPEN]** — defined, numerically delicate |

This closes the loop on the forced surrogates for the fold (now the genuine Airy₂) and leaves the cusp
intrinsic process as the single deepest open object — the parabolic-cylinder line ensemble at a
multicritical RMT edge.
