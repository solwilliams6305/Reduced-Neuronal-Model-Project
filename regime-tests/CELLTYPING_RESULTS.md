# Excitability cell-typing from extracellular spikes — a secondary result (honest summary)

*A consolidated, caveated results note for the excitability arm of the inversion project. Full working
detail: `EXCITABILITY_CHARACTERIZER.md`. Code: `excitability_characterizer.py`, `morris_lecar_check.py`,
`ml_atlas_separability.py`. Figures: `figures/excitability_characterizer.png`, `…/morris_lecar_check.png`,
`…/ml_atlas_separability.png`.*

## Claim

Read a neuron's **excitability class** (Type-I/SNIC vs Type-II/Hopf) and its **proximity to spiking
threshold** from **extracellular spike statistics across a few drive levels** — a dynamical invariant
that normally needs intracellular access — with a goodness-of-fit that **refuses cells it cannot
characterise**.

## Method

Three fused features: (1) near-threshold-weighted ISI-**shape** class (skewed quartic-FPT = Type-I vs
near-symmetric Hopf jitter = Type-II); (2) **f–I onset** (Hodgkin Class 1 vs 2: Type-I reaches low rate
near rheobase, Type-II is bounded) — the **adaptation-robust arbiter**, since adaptation regularises the
shape but does not move the onset bifurcation; (3) the **W1 residual** as a reject gate. Plus a
confound-corrected proximity estimate across drive.

## What we found, in order of honesty

1. **In-family (normal-form cohort, 18 cells): 100%** three-way accuracy (Type-I / Type-II /
   out-of-family) after using f–I onset to neutralise the adaptation confound; out-of-family rejected
   with a ~10× residual margin. **This is an in-family number** — the cohort and the atlas share the
   reduced-model family.

2. **Conductance-model reality check (Morris–Lecar): all cells auto-rejected.** Residual ≈ 0.091, at the
   garbage level (0.108), far above in-family (0.007). Morris–Lecar leaves the noise-driven first-passage
   regime at moderate drive (regular firing, skew flips negative) — outside the normal-form atlas. The
   self-diagnosing residual makes the mismatch **visible** (refuses rather than mislabels).

3. **Kill-or-confirm (train the atlas on Morris–Lecar itself): PARTIAL.** Shape-only accuracy **69%**
   all-drive (vs 89–100% for normal forms), **78% near rheobase**; the classes are distinct near
   threshold (skew 0.66 vs 1.36) but separation collapses at higher drive (Type-II → 25–38%). f–I onset
   still separates directionally (0.32 vs 0.56), gap narrowed from the idealised 0.47/0.87.

## Honest bottom line

The excitability-from-spikes premise is **real but modest**. It lives in a **narrow near-threshold
window** and is **weak**: the realistic ceiling is ≈ **80%** with a biophysically-built atlas,
near-rheobase sweeps, and shape combined with f–I onset. That is **population-level cell-type
characterisation, not per-cell certainty** — a sober claim, but a defensible and useful one, and the
reject mechanism is what keeps it honest. The idealised-normal-form 100% was, in hindsight, flattering.

## What would make it a primary result

A biophysically-trained atlas (Morris–Lecar / Hodgkin–Huxley ISIs across drive and noise), a combined
shape + f–I-onset classifier restricted near rheobase, then validation against **Allen Cell Types**
intracellular ground truth (`allen_excitability.py`, local). Until then this is a **secondary result**,
honestly bounded.
