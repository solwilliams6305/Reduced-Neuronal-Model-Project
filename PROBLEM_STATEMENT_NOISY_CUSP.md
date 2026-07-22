# Problem statement — The Noisy Folded Cusp (Pearcey-class canard escape)

_June 2026. A scoping document for the strongest creative direction (see CREATIVE_DIRECTIONS.md §I.1),
with first numerical evidence already in hand._

---

## Conjecture

> The noise-induced peel-off (canard escape) law of a folded **limit cycle at a cusp** of the
> cycle manifold is a **non-Tracy–Widom edge universality class** — the cusp analogue of the
> fold's TW law. More generally, **each catastrophe of the cycle manifold carries its own
> noise-induced escape universality**, with fold → Tracy–Widom as rung one of a
> "diffraction-catastrophe ladder."

The fold rung is the NoisyFoldedCycle paper: a fold + white noise → the stochastic Airy
operator → TW_β. The cusp rung is the natural — and reachable — sequel.

## Why this is the right next problem

1. **Coupling manufactures the cusp.** Two symmetrically coupled folded cycles do not fold
   as "two folds"; at the synchrony locus they degenerate into a **cusp**. This is exactly
   what Kristiansen–Pedersen (arXiv:2202.12027, SIADS 2023) found for coupled FitzHugh–Nagumo
   MMOs — "singularities at a cusp, not a fold," with the small-oscillation count governed by
   the **Weber (parabolic-cylinder) equation**. So the cusp is not exotic: it is what our
   coupled project already sits on.
2. **The deterministic correspondence is classical.** Berry–Upstill diffraction catastrophes:
   fold ↔ Airy function, cusp ↔ Pearcey function. The paper made the fold↔Airy correspondence
   *stochastic*; the cusp↔Pearcey correspondence is the obvious next stochastic statement.
3. **The deterministic ladder exists in RMT.** Multicritical random-matrix edges are governed
   by the higher-order Airy / higher-order Tracy–Widom hierarchy (Le Doussal–Majumdar–Schehr;
   Cafasso–Claeys–Girotti — ⚠ verify exact attributions). The cusp is the first multicritical rung.
4. **It resolves our hardest gap.** The coupled noisy blow-up (the analytic gap we flagged as
   having no current idea) *is* the cusp blow-up. Solving the cusp solves that.

## First numerical evidence (this session)

`coupled-atlas/peeloff_cusp_ladder.py` isolates the universality question with a controlled
model: the inner Cole–Hopf equation with a tunable **turning order** k,

> u″ = ( sign(Y)·|Y|^k − η ξ ) u ,   Y = Y₀ − T ,   (k = 1 is the fold → Airy → TW),

which is precisely the higher-order-Airy family. Result (η = √2, β = 2 at the fold):

| k (turning order) | peel-off skew | excess kurt |
|---|---|---|
| 0.7 | 0.34 | 0.99 |
| **1.0 (fold)** | **0.18 ≈ TW₂** | 0.09 |
| 1.5 | 0.43 | −0.25 |
| 2.0 | 0.60 | −0.26 |
| 3.0 | 0.96 | 0.19 |

The escape law **moves systematically with the turning order**, and the fold (TW) sits at a
distinguished, least-skewed point. This is direct evidence that the ladder is real — different
singularity order ⇒ different edge law (figure: `figures/peeloff_cusp_ladder.png`).
**Caveat:** k is a model knob for singularity order, not the literal cusp; absolute moments are
setup-dependent (the trend at fixed setup is the result); k = 1 reproduces TW (code validated).

## The precise problem (four steps to a theorem)

1. **Inner equation at the coupled cusp.** Carry out the blow-up at the synchrony (cusp) point
   of the coupled folded cycle and identify the inner equation and its turning order k_cusp.
   Expected: a 2nd-order Weber-type equation (per Kristiansen–Pedersen) or a higher-order-Airy
   equation. *This is the hard GSPT step — the coupled noisy blow-up.*
2. **Add noise.** Cole–Hopf-linearise the cusp inner Riccati; identify the resulting
   "stochastic Weber / stochastic Pearcey" operator (the cusp analogue of RRV's stochastic Airy).
3. **Edge law.** Characterise its ground-state / first-node distribution: is it a named
   higher-order Tracy–Widom, the Pearcey-edge law, or genuinely new? Match moments to a
   multicritical-edge reference.
4. **Blow down.** Transport to physical (σ, ε) variables; derive the cusp critical-noise law
   σ\*_cusp (the analogue of the paper's σ\* = 2√π·√ε₂·G) and cross-check in the coupled FHR.

## Line of attack

- ✅ **Step 0 (model):** turning-order probe — law moves with k (done, above).
- ⏭ **Step 1 (numeric cusp):** simulate the actual coupled *symmetric-subspace* inner equation
  through the cusp; confirm the peel-off law is non-TW and read off k_cusp. (Cheap; next.)
- ⏭ **Step 2 (reference):** build a multicritical-edge / higher-order-Airy sampler; test whether
  the cusp law matches a known higher-order-TW or is new.
- ⏭ **Step 3 (analysis):** the coupled cusp blow-up (the rigorous core; nearest tools — RRW
  averaging for the deterministic coupled folded singularity, K–P cusp blow-up, Berglund–Gentz
  tubes for the noise).
- ⏭ **Step 4:** blow-down + physical critical-noise law + FHR cross-check.

## Risks and honest unknowns

- **Pearcey vs Weber naming.** The pure diffraction cusp (Pearcey) is a *3rd-order* problem,
  which breaks the clean 2nd-order Cole–Hopf → Schrödinger route; but the *slow-fast* cusp
  reduces to the 2nd-order Weber equation (K–P). So the right object may be a "stochastic Weber"
  edge law rather than a literal Pearcey. Determining which is part of Step 1–2.
- **Known or new?** The cusp law may coincide with an established higher-order TW, or be a new
  distribution. Either is publishable; the test is Step 2.
- **The blow-up (Step 3)** remains the genuine analytic frontier; Steps 0–1 give strong evidence
  without it, but the theorem needs it.

## Payoff

A second rung of a noise-induced escape universality ladder — a new edge law beyond Tracy–Widom,
grounded in catastrophe theory, reachable from the coupled folded cycle we already built, and
resolving the coupled-blow-up gap in the process. Working title: **"The Noisy Folded Cusp:
Pearcey-class statistics of coupled canard escape."**

## References

Williams 2026 (NoisyFoldedCycle); Kristiansen–Pedersen 2023 (arXiv:2202.12027, coupled-FHN cusp);
Jelbart–Kuehn–Kuntz 2024 (arXiv:2208.01361, folded-limit-cycle blow-up); Roberts–Rubin–Wechselberger
2015 (coupling → folded singularities of the cycle manifold); Wechselberger (folded node / Weber);
Ramírez–Rider–Virág 2011 (stochastic Airy operator); Berry–Upstill (diffraction catastrophes);
Le Doussal–Majumdar–Schehr, Cafasso–Claeys–Girotti (higher-order / multicritical edge — ⚠ verify).
