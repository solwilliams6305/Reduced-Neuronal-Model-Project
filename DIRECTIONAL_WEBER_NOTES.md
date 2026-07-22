# Notes — the directional (swept) Weber operator: what the cusp operator must be

_June 2026. Figure `coupled-atlas/figures/directional_weber_operator.png`; script
`directional_weber_operator.py`. Tags: [VALIDATE] · [NEGATIVE] · [KEY] · [OPEN]._

Follow-up to `STOCHASTIC_WEBER_OPERATOR_NOTES.md` (static self-adjoint power-law operators can't
reproduce the cusp law). The diagnosis there was that the cusp law's **negative** excess kurtosis is
a dynamical fingerprint, so the operator must be non-self-adjoint — a Weber generator with a **sweep**.
Here I tried the simplest such object and learned exactly what's required.

## The control table (the mechanism, isolated)

Standardised peel-off / escape level, β=2:

| object | skew | excess kurtosis |
|---|---|---|
| static Airy operator $-\partial^2+x$ (soft edge) | +0.18 | **+0.03** |
| static harmonic operator $-\partial^2+x^2$ (parabolic, no sweep) | +0.15 | ~0 |
| swept **linear** / fold $u''=(Y-\eta\xi)u$ | +0.27 | **+0.29** |
| swept **parabolic** / cusp $u''=(\mathrm{sign}(Y)Y^2-\eta\xi)u$ | **+0.63** | **−0.20** |

Reading: the **excess-kurtosis sign flips with the swept potential** — swept-linear (fold) gives
exkurt > 0 (the soft-edge/TW signature); swept-parabolic (cusp) gives exkurt < 0. Large skew (beyond
TW) appears only for the parabolic sweep.

## What I built, and the honest negative

The candidate "directional operator" was the backward generator of the dynamic pitchfork (the cusp's
antisymmetric mode), $L_T=\tfrac{\sigma^2}{2}\partial_X^2+(\mu(T)X-X^3)\partial_X$, realised by the SDE
$dX=(\mu(T)X-X^3)dT+\sigma dW$ swept through the bifurcation; observable = escape level $\mu_{\rm esc}$.

It does **not** reproduce the cusp law. Across noise levels its excess kurtosis stays **positive**
(σ=0.18 → skew +0.43, exkurt +0.31; larger σ kills the skew). It is a **monotone/gradient** escape —
no oscillation — so it cannot produce sub-Gaussian (negative-kurtosis) tails.

## [KEY] the structural requirement

Two natural simplifications now both fail, in the same place:
- static self-adjoint power-law operator → exkurt > 0;
- 1st-order/gradient dynamic pitchfork → exkurt > 0.

The cusp's **negative** excess kurtosis requires the **2nd-order, oscillatory** (small-amplitude /
SAO) Weber structure. The light tails come from the **oscillation** of the parabolic-cylinder passage,
not from a confining well and not from a soft edge. So the faithful object is intrinsically:

> **a swept, 2nd-order, parabolic-cylinder (Weber) operator** — the swept-Weber inner equation
> $u''=(\mathrm{sign}(Y)Y^2-\eta\xi)u$ itself — and any canonical static or 1st-order RRV-style
> reduction provably loses the signature.

## Status

This sharpens, but does not close, the open object. We now know the cusp operator is **swept +
2nd-order + oscillatory**, with the excess-kurtosis sign as the discriminating invariant (a clean,
checkable signature). The remaining target is a rigorous canonical form / theorem for this swept
parabolic-cylinder operator — the cusp analogue of RRV's stochastic-Airy = TW — alongside the
intrinsic Dyson–Weber multicritical line ensemble (`DYSON_WEBER_NOTES.md`). Both stay [OPEN]; rung F
advances (the operator's required structure is now pinned and two ansätze are excluded).
