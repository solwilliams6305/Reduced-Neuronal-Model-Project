# Notes — trying to construct an RRV-type stochastic Weber operator

_June 2026. Figure `coupled-atlas/figures/stochastic_weber_operator.png`; script
`stochastic_weber_operator.py`. Tags: [VALIDATE] · [NEGATIVE] · [HEURISTIC] · [OPEN]._

**Question.** Is there an RRV-type stochastic operator whose ground state is the cusp (Weber) escape
law — the cusp analogue of "stochastic Airy operator = Tracy–Widom"?

**Setup.** Build $H = -\partial_x^2 + V(x) + \frac{2}{\sqrt\beta}b'(x)$ on $[0,\infty)$, Dirichlet,
by the same tridiagonal + Sturm-bisection method that reproduces TW for the fold. RRV's fold operator
has the **linear** $V(x)=x$ (the fingerprint of a √-edge: Weyl count $N(E)\sim E^{1/2+1/p}$ gives
$E^{3/2}$ at $p=1$). The earlier no-go (`spectral_id.py`) showed $V=x^2$ (super-linear/harmonic)
confines like a bulk → Gaussian. So the multicritical edge should need a **sub-linear** $V=x^p,\ p<1$
(Weyl: the $m$-th multicritical edge sits at $p=1/m$).

## What worked

- **[VALIDATE]** $p=1$ reproduces TW (operator skew +0.27 ≈ swept fold +0.27 ≈ TW₂ +0.22); $p=2$
  reproduces the harmonic no-go (skew +0.05, Gaussian). The method is sound.
- **[HEURISTIC ✓]** Sub-linear $p<1$ does **raise** the ground-state skew above TW (plateau ≈ +0.32 for
  $p\in[0.4,0.85]$) — the Weyl direction is correct.

## What did not (the honest negative)

The static power-law family **cannot** reproduce the swept Weber cusp law (skew +0.62, exkurt **−0.20**):

| | skew | excess kurtosis |
|---|---|---|
| swept Weber (cusp target) | **+0.62** | **−0.20** (light, sub-Gaussian tails) |
| operator family $V=x^p$, any $p$ | saturates ≈ +0.32 | **positive** (heavy, soft-edge tails) |

Two independent failures:
1. **Skew saturates** ≈ +0.32, far short of +0.62 — no power reaches it.
2. **Excess-kurtosis sign is wrong.** Every soft-edge operator has positive excess kurtosis; the swept
   Weber law has **negative**. A β-probe at $p=0.45$ confirms the trap: stronger noise ($\beta=0.5$)
   finally drives exkurt negative (−0.24) but then the skew **collapses** to +0.21. No $(p,\beta)$
   delivers skew +0.62 **and** exkurt < 0 together.

## Reading — where the obstruction lives

The mismatch is diagnostic, not just a failure:
- The cusp law's **negative** excess kurtosis (sub-Gaussian, light tails) is the signature of
  **parabolic-cylinder (Weber) confinement** — the opposite of a soft edge's one heavy tail.
- Its large **+0.62 skew** is largely a **swept-passage** (dynamical) effect: the system is driven
  through $V(Y)=\mathrm{sign}(Y)|Y|^2$ directionally, and the first node inherits the asymmetry.

So the faithful cusp operator is **not** a softened soft-edge power law; it must encode **sweep +
confinement** (e.g. a parabolic-cylinder operator with a directional/non-self-adjoint term), which is
why the static self-adjoint power-law ansatz cannot work. This **rules out the simplest ansatz** and
**localises** the open problem — genuine progress, but the operator itself remains **[OPEN]**.

## Status update

- The RRV-route to the cusp marginal is now characterised: simplest ansatz excluded, obstruction pinned
  (kurtosis sign ⇒ confinement; skew ⇒ sweep). Rung F advances; rung B's intrinsic object is unchanged
  (still open). The remaining constructive target is a **swept/confined parabolic-cylinder operator**,
  alongside the **Dyson–Weber** multicritical line ensemble (`DYSON_WEBER_NOTES.md`).
