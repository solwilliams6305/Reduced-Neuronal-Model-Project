# Notes — the loop closed: the ladder demonstrated in the coupled neuron model

_June 2026. Figure `coupled-atlas/figures/loop_closure_scaling.png`; script `loop_closure_scaling.py`.
Tags: [NUMERIC] crossover demonstrated · [HONEST] limits noted._

The keystone: the catastrophe-ladder theory was *derived from* the coupled neuron model — here it is
*demonstrated in it*. Genuine Kristiansen–Pedersen coupled FHN with noise:
$$v_i' = -v_i^3 + 3v_i - w_i + g(v_j-v_i) + \sigma\xi_i,\qquad w_i' = \varepsilon(v_i-c),\quad c=0.99.$$
Observable: the single-unit peel-off (w₁ at the upper fold) under noise, swept over coupling $g$.

## What is demonstrated

As $g$ sweeps from weak/attractive (generic **fold**) toward the antisymmetric **cusp**, the peel-off
escape law crosses edge-class, with two robust signatures (both $\varepsilon=0.006$ and $0.012$):

| signature | fold side ($g\!\ge\!-0.03$) | cusp side ($g\!\le\!-0.09$) |
|---|---|---|
| spread (std) | small | **×3.7 amplified** |
| excess kurtosis | ≈ 0 (Tracy–Widom-like) | **negative** (−0.07 to −0.65) — **sub-Gaussian Weber-TW** |

The **excess-kurtosis sign flip** (≈0 → negative) is exactly the ladder's fold→cusp signature: the
fold rung is TW (small kurtosis), the cusp rung is the sub-Gaussian Weber-TW. **It appears in the real
coupled neuron model.** The spread amplification is the canard-funnel widening as the antisym folds
merge. So the abstract edge classes are not artifacts of the reduced inner equation — they govern the
noise-induced escape statistics of the actual two-neuron system.

## Honest limits

- **Outer corrections.** The full-model peel-off carries outer (Gaussian) contributions, so the
  *absolute* cusp kurtosis (≈ −0.1 to −0.7, $\varepsilon$-dependent) is not the clean inner Weber-TW
  value (−0.23); the robust, model-independent signal is the **crossover** (the flip + the spread
  amplification), not the absolute moment.
- **Skew sign is $\varepsilon$-dependent** (outer corrections flip it), so skew is not a reliable
  signature here — kurtosis and spread are.
- **The $g_{\rm crit}\propto\sqrt\varepsilon$ onset was *not* cleanly resolved** by the spread-onset
  metric (the onset came out roughly $\varepsilon$-independent, ≈ −0.08…−0.10). The derived
  $g_{\rm crit}\propto\sqrt\varepsilon$ is the folded-saddle→folded-node (SAO-onset) transition, which
  need not coincide with the spread-amplification onset; resolving it needs an SAO-count / folded-node
  metric. **Open.**

## Status

The loop is **qualitatively closed**: the fold→cusp escape-class crossover (kurtosis flip to
sub-Gaussian + spread amplification) is demonstrated in the genuine coupled FHN model — the ladder
lives in the neurons. The *quantitative* onset law and the clean inner-moment match remain open
(outer-correction subtraction; an SAO-based onset metric). Rungs D/E advance; the spine (A–C, F) is now
physically anchored.
