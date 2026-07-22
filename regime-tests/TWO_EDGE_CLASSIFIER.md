# Direction D: a noise-only classifier for excitability class — the two-edge fingerprint

*The noisy folded cycle has two edges; they map onto the two classical excitability classes, with
**opposite, measurable noise signatures**. A Type-I (SNIC) neuron carries the noise on the **phase**
edge (spike timing, `σ^{2/3}`); a Type-II (Hopf / fold-of-cycles) neuron carries it on the
**amplitude** edge (oscillation level, `σ^{1/2}`). The threshold scaling exponent alone classifies
the cell. Companion: `two_edge_classifier.py`, figure `figures/two_edge_classifier.png`. Tags
**[R]/[N]**.*

---

## The classifier [R]

| | **Type-I (SNIC / QIF)** | **Type-II (Hopf / fold-of-cycles)** |
|---|---|---|
| edge | **phase / time** | **amplitude / level** |
| noisy order parameter | spike **timing** (ISI) | oscillation **amplitude** |
| threshold noise scaling | **`r ∼ σ^{2/3}`** | **`A ∼ σ^{1/2}`** |
| f–I onset | **continuous** (`f → 0`) | **discontinuous** (`f → f_H ≠ 0`) |
| fluctuation law | **quartic FPT** (skewed; not TW) | Rayleigh — **TW for fold-of-cycles**\* |

\*A *generic* Hopf gives **Rayleigh** amplitude fluctuations. The **Tracy–Widom** fine structure of
the amplitude edge is specific to the **fold-of-cycles canard peel-off** (the burster) — the
project's main result — where the order parameter is the canard-escape *level* and the inner operator
is the stochastic Airy operator. The classifier's robust, model-independent content is the **exponent
and onset** dichotomy; TW is the extra fingerprint carried by the bursting member.

## The numerics [N]

`two_edge_classifier.py` simulates both members at threshold and fits the order-parameter exponent:

```
   Type-I  (noisy QIF, I=0)            threshold rate  r ~ sigma^0.668     (predict 2/3 = 0.667)
   Type-II (noisy Stuart-Landau, mu=0) threshold amp   A ~ sigma^0.494     (predict 1/2 = 0.500)
```

Both exponents hit their predicted values to <1.5%. The figure shows the three discriminators:

* **Onset dichotomy (panel A).** Type-I firing frequency rises continuously from zero (noise-smoothed
  `√I`); Type-II frequency jumps from rest to `f_H = ω/2π ≠ 0` at threshold — Hodgkin's original
  Type-I/Type-II distinction, recovered.
* **Exponent classifier (panel B).** On log–log axes the Type-I rate sits on the `σ^{2/3}` line and
  the Type-II amplitude on the `σ^{1/2}` line — two distinct, cleanly separated slopes. This is the
  core: **measure the threshold order-parameter exponent → read off the class.**
* **Fingerprint table (panel C).** The full two-channel signature.

## Why the exponents are what they are [R]

* **Type-I, `σ^{2/3}`.** At the SNIC the ISI is the first-passage time of `dv=(v²+I)dt+σdW`; the inner
  rescaling `v=σ^{2/3}ρ, t=σ^{-2/3}s` is noise-free, so the rate carries one factor of `σ^{2/3}`
  (the phase edge; see `QIF_PHASE_EDGE.md`).
* **Type-II, `σ^{1/2}`.** At the Hopf point the radial normal form `dr=(−r³+σ²/2r)dt+σdW_r` balances
  `r³ ∼ σ²/r`, giving `r ∼ σ^{1/2}` (the amplitude edge). The fold-of-cycles refinement replaces the
  smooth `−r³` relaxation by the canard peel-off, turning the Rayleigh amplitude into `TW_β`.

## Reproduce

```
python3 two_edge_classifier.py     # -> figures/two_edge_classifier.png  (prints both exponents)
```

## Verdict & how to use it

The two excitability classes are **noise-distinguishable without resolving the deterministic
bifurcation**: a Type-I/SNIC cell shows continuous f–I onset and a `σ^{2/3}` timing edge with a
quartic ISI law; a Type-II/Hopf cell shows a discontinuous onset frequency and a `σ^{1/2}` amplitude
edge; and the **bursting / fold-of-cycles** member additionally carries the **`TW_β`** peel-off-level
law. Experimentally: classify with (i) the f–I onset, (ii) the threshold scaling exponent of the
relevant order parameter (timing vs amplitude) under controlled noise, and — for the headline — (iii)
intracellular canard-escape amplitude statistics in a burster to test `TW_β`. This closes the
Type-I half (done end-to-end, Direction #1 / A) against the Type-II/amplitude half, and is the
experimental face of `TWO_CHANNEL_SYNTHESIS.md`.
