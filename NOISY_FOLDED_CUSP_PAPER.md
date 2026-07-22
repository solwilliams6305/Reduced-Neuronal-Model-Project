# The Noisy Folded Cusp: Weber-class canard escape and a catastrophe ladder of noise-induced escape

**S. Williams — working draft, June 2026.** Companion to *The Noisy Folded Limit Cycle:
Tracy–Widom Statistics of Canard Escape*. Figures in `coupled-atlas/figures/`; scripts named inline.
Status tags: [PROVED] established/derived exactly · [DERIVED] asymptotic balance · [NUMERIC] simulated ·
[VALIDATED] reproduces a known result · [OPEN] conjecture/theorem-to-be.

---

## Abstract

The noisy folded limit cycle places single-unit canard escape in the Tracy–Widom (TW) / KPZ edge
class: a fold of the cycle manifold plus white noise gives the stochastic Airy operator, whose ground
state is TW. We ask what happens when two such units are **electrically coupled**. Using the
Kristiansen–Pedersen coupled FitzHugh–Nagumo model, we show that coupling does not merely correlate
two copies of the fold law — it changes the **singularity**: the antisymmetric mode is exactly a
**cusp**, and its noise-induced escape is a distinct, **Weber-class** edge law, *not* the third-order
Pearcey one might expect from the diffraction-catastrophe correspondence. We derive the fold
separation Δ(g)=2√(−2g/3) from the slow manifold, show a coupling-driven **fold→cusp crossover** with
onset g_crit ∝ √ε (derived and verified), and construct the cusp's **multi-point process** — a
stationary Weber-marginal process, the cusp analogue of the Airy₂ process. Together these populate a
*catastrophe ladder of noise-induced escape*: fold→Airy/TW→Airy₂; cusp→Weber→Weber-process; coupling,
forcing and heterogeneity are the unfolding parameters that move between rungs.

---

## 1. Introduction

The deterministic backbone is geometric blow-up of folded singularities (Krupa–Szmolyan; Jelbart–
Kuehn–Kuntz for the folded limit cycle). The noisy folded *cycle* paper proved: the inner Riccati
Cole–Hopf-linearises to the stochastic Airy operator (Ramírez–Rider–Virág), so the peel-off level is
TW_β, β=4/η². A closing problem (§9.4 there) asks for the joint law of successive peel-offs in a
**forced or coupled** folded cycle.

The diffraction-catastrophe hierarchy (Berry–Upstill) attaches a special function to each elementary
catastrophe: fold↔Airy, cusp↔Pearcey, …. Since the fold rung of a *noisy* version of this hierarchy
is exactly the paper's result, it is natural to ask whether the higher rungs exist — and coupling is
the physical knob that reaches them.

## 2. Model and framing

Kristiansen–Pedersen coupled FHN (gap-junction coupling on the fast variable):

  v_i' = −v_i³ + 3v_i − w_i + g(v_j − v_i),  w_i' = ε(v_i − c),  i=1,2.

Single units (g=0) are relaxation oscillators (|c|<1); KP prove that symmetric repulsive coupling
produces MMOs organised by a **cusp** of the critical manifold, with the small-oscillation count
governed by the **Weber (parabolic-cylinder)** equation. We locate this regime numerically at
c≈0.99, g≈−0.12 (`kp_cusp_explore.py`).

Organising question: **which singularity → which edge class → which multi-point process**, with
g, forcing, heterogeneity as unfolding parameters.

## 3. Δ(g): the antisymmetric mode is a cusp  [PROVED — exact]

In symmetric/antisymmetric coordinates v_{1,2}=v_s±δ, the antisymmetric mode obeys

  δ' = −δ³ + μ δ − δw,  **μ(v_s,g) = 3(1−v_s²) − 2g**,

i.e. exactly the cusp normal form. At the symmetric fold v_s=1, μ=−2g, and the antisymmetric folds
sit at δ=±√(μ/3), giving the **fold separation**

  **Δ(g) = 2√(−2g/3)**  (repulsive g<0; for g≥0, μ≤0, synchrony is stable → bare fold).

The cusp is μ=0 (g=0); the Weber funnel is the folded-node regime g<0. This replaces the earlier
geometry-argued map with an exact derivation (`DELTA_G_BLOWUP_NOTES.md`).

## 4. The cusp edge class is Weber, not Pearcey

[PROVED] (KP): the slow-fast cusp reduces to the **2nd-order Weber** equation, not the 3rd-order
Pearcey ODE — the slow flow collapses Pearcey to Weber.
[NUMERIC] (`weber_vs_pearcey.py`): the genuine KP cusp escape law matches Weber-class and rejects
Pearcey — skew 0.61 vs Weber 0.64 (Pearcey 0.12, fold 0.21); excess kurtosis negative (Weber sign,
Pearcey wrong sign); KS 0.056 to Weber vs 0.123 to Pearcey; the compressed lower tail tracks Weber.
Figure `weber_vs_pearcey.png`.

## 5. The fold↔cusp crossover  [DERIVED scaling + NUMERIC]

Coupling drives a continuous crossover of the peel-off edge class (`crossover_fold_to_cusp.py`): as
repulsive coupling merges the two folds into the cusp, the peel-off spread amplifies ×3.6 and the
shape swings from fold-class to cusp-class. The deterministic backbone is **Olver's uniform
connection** (isolated turning point → Airy; two coalescing → Weber), realised in the A₃ blow-up
chart. The operative onset is the folded-saddle→folded-node transition at μ~√ε, i.e.

  **g_crit ∝ √ε**  (verified: g_crit/√ε = −0.58, −0.57 at ε=0.015, 0.030; ratio 1.40 ≈ √2),

with onset inside the measured band −0.05…−0.11 (`delta_g_blowup.py`, figures `crossover_*`,
`delta_g_blowup.png`).

## 6. The cusp successive-peel-off process — Weber analogue of Airy₂  [NUMERIC]

Successive forced passages of the cusp inner equation (`forced_cusp_process.py`) give a **stationary
Weber-marginal process** (skew 0.57, distinct from the fold's TW skew 0.17), with decaying
autocovariance under forcing (C: 0.76→0.20 vs iid ≈0) and 2·Var-saturating increments — the same
process signatures as the fold's Airy₂-type result, now on the Weber marginal (figure
`forced_cusp_process.png`). This completes the spine symmetrically:

| rung | marginal | process |
|---|---|---|
| fold | Tracy–Widom | Airy₂-type |
| cusp | Weber | Weber-type |

## 7. The stochastic step  [DERIVED / NUMERIC / OPEN]

- **Noise scale [DERIVED]:** pushing σ dW through the A₃ blow-up gives η_cusp = σ/(√2 ε^{2/5}); the
  same Itô-on-blow-up reproduces the fold's η = σ/√ε₂ (validation).
- **Tube [NUMERIC]:** the deviation from the canard is an OU process with rate a(s)=−3D̄²+M; the
  tube variance follows η²/(2|a|) (ratio 0.85) and the confinement is Gaussian (κ≈1.85). Matching
  height d\*~η^{2/3}/M^{1/6}. (`cusp_tubes.py`.)
- **Spectral id:** the tridiagonal Sturm-operator method **reproduces TW at the fold [VALIDATED]**;
  for the cusp the naive harmonic operator is **provably wrong** (Gaussian, skew 0.05 ≠ 0.62) — the
  correct operator is a **soft-edge multicritical / higher-order-Airy** one [OPEN], obstruction
  pinned (`spectral_id.py`).

## 8. Discussion and open problems

The picture is a **catastrophe ladder of noise-induced escape**, with coupling as the unfolding
parameter. The fold rung is the prior paper; the cusp rung (this work) is a genuinely distinct edge
class (Weber) with its own multi-point process, reached continuously from the fold by a √ε-scaled
crossover. Open:

1. **Intrinsic processes (both rungs):** the genuine Airy₂ (fold) and "Dyson–Weber" (cusp) edge
   processes from the Dyson-Brownian-motion edge, with intrinsic covariance (our forced results use
   an imposed-OU surrogate). [the next push]
2. **The soft-edge multicritical operator** whose ground state is the cusp peel-off (the cusp
   analogue of RRV's stochastic-Airy = TW).
3. **The rigorous coupled cusp blow-up + tube theorem** (constants, uniform-in-M; the stochastic
   uniform Airy↔Weber connection).
4. Higher rungs (swallowtail, …) and the genuine coupled-FHN *forced* cusp.

## 9. Reproducibility

`kp_cusp_explore.py` · `kp_cusp_noise.py` · `weber_vs_pearcey.py` · `crossover_fold_to_cusp.py` ·
`delta_g_blowup.py` · `crossover_airy_to_weber.py` · `cusp_tubes.py` · `spectral_id.py` ·
`stochastic_cusp.py` · `forced_cusp_process.py` (all numpy/matplotlib, self-contained).
Notes: `DELTA_G_BLOWUP_NOTES.md`, `NOISY_CUSP_CROSSOVER_NOTES.md`, `STOCHASTIC_CUSP_NOTES.md`,
`CUSP_TUBES_SPECTRAL_NOTES.md`, `CUSP_PROCESS_NOTES.md`. Lit map: `COUPLED_REDUCED_NEURON_LIT_MAP.md`.

## References (abbreviated)

Williams 2026 (NoisyFoldedCycle); Kristiansen–Pedersen 2023 (arXiv:2202.12027); Jelbart–Kuehn–Kuntz
2024; Ramírez–Rider–Virág 2011; Wechselberger (folded node / Weber); Berry–Upstill (diffraction
catastrophes); Olver (uniform asymptotics, coalescing turning points); Le Doussal–Majumdar–Schehr,
Cafasso–Claeys–Girotti (multicritical / higher-order edges); Dauvergne–Ortmann–Virág (directed
landscape); Roberts–Rubin–Wechselberger 2015 (coupling → folded singularities).
