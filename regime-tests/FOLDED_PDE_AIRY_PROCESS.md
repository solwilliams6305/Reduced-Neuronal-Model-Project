# Airy process from a folded PDE — what coupling actually delivers it (banked #4)

*Question: can a spatially-extended / many-mode folded system produce the Airy₂
**process** (the KPZ-class curve), not just the TW_β **marginal** (rung 1) or the
single-operator Airy point **process** (rung 2)? Your C-arc found nearest-neighbour
arrays give only a local exponential shadow. This note isolates **which coupling
structure supplies the genuine log-gas repulsion** the Airy₂ process needs.*

Scripts: `folded_pde_airy_process.py` · figure: `figures/folded_pde_airy_process.png`.

## Setup

N coupled Cole–Hopf fields linearised about the canard give
`U''(x) = [ x·I + C − η·Ξ(x) ] U(x)`, so the **coupling matrix `C`** is the
random-matrix ensemble whose **edge** governs the joint peel-off law. The direction
reduces to: *for which `C` is the swept top-eigenvalue the Airy₂ process?*

## Results

**[N] Reference object — Dyson Brownian motion edge IS Airy₂.** A Hermitian OU process
(stationary law = GUE), edge-rescaled top eigenvalue `ξ(t)=N^{1/6}(λ_max−2√N)`:
marginal matches TW₂ (mean −1.78, var 0.81, skew ✓), and the increment variance
`Var[ξ(t+Δ)−ξ(t)]` is **linear at small lag (locally Brownian) and saturates at
1.62 ≈ 2·Var(TW₂)=1.626**. That is the Airy₂ fingerprint, reproduced parameter-free.
This pins the target the folded system must hit.

**[N] Mechanism — only disordered all-to-all coupling repels like a log-gas.**
Level-spacing ratio `⟨r⟩` (Poisson 0.386, GOE 0.531, GUE 0.603) and the standardised
edge law, N=160, 1200 samples:

| coupling `C`                    | `⟨r⟩` | edge skew | verdict |
|---------------------------------|-------|-----------|---------|
| **random all-to-all (GUE)**     | 0.599 | +0.16     | TW₂ edge **+ Dyson repulsion → Airy₂** ✓ |
| nearest-neighbour (Anderson)    | 0.391 | +0.71     | Poisson, **localised** → not TW ✗ |
| uniform mean-field (rank-1)     | 0.601*| +1.12     | single **BBP outlier**, not TW ✗ |

`*` the mean-field bulk `⟨r⟩` is high only because rank-1 interlacing between sorted
on-site values fakes spacing rigidity; its **edge** is the detached collective mode
(BBP), not a TW edge — so it does **not** give the process. Real-symmetric coupling
gives `⟨r⟩=0.533` (GOE) → **β=1**, i.e. the **Airy₁** sibling (the symmetry-class
selection flagged earlier: real coupling → TW₁/Airy₁, complex → TW₂/Airy₂).

**[refuted, methodological] The increment-variance *shape* is NOT a discriminator.**
My first instinct — "linear-then-saturating increments ⇒ Airy₂" — is wrong: *any*
bounded stationary process has that shape, and a standardised nearest-neighbour top
eigenvalue mimics it perfectly. The real discriminators are (i) genuine repulsion
`⟨r⟩` and (ii) the **edge marginal** (TW vs Gumbel vs BBP). Keep this; it would have
been an easy false positive.

## Refined verdict (sharpens banked #4)

Airy₂ needs the transverse operator to be a **long-range AND disordered** random
matrix. The two natural "long-range" fixes fail in *opposite* ways:

- **nearest-neighbour / any fixed-range** spatial coupling → Anderson **localisation**
  → Poisson edge (your exponential shadow, now explained by localisation);
- **uniform mean-field** (homogeneous all-to-all) → rank-1 → **self-averages to one
  BBP collective mode** (propagation-of-chaos), not a log-gas edge;
- and a **homogeneous (translation-invariant) folded PDE** is Fourier-diagonal → a
  smooth dispersion spectrum (regular "picket-fence" spacings), also no Airy.

So the needle to thread is a **disordered, genuinely long-range** folded medium:
quenched-random or sparse-random-network coefficients, **not** a clean homogeneous
PDE. "Airy₂ from a *disordered* folded medium" is the corrected target.

Direct tie-ins: the mean-field failure **is** banked **#7 (BBP transition)** — the
collective mode is the spike; tuning `κ` through the edge is the detachment. And the
result upgrades #4's "all-to-all" to "**disordered** all-to-all".

## Next (if pursued)

1. Build the faithful model: a folded reaction–diffusion / neural-field PDE with
   **quenched spatial disorder** (random coefficients) or a **random sparse network**
   of folded units; confirm its canard-linearisation edge → GOE/GUE `⟨r⟩` and a TW
   edge, then sweep the slow variable and measure the top-line increment variance
   against the 2·Var(TW) plateau.
2. Symmetry class as a prediction: real heterogeneity → Airy₁; add a directed/flux
   term (broken reciprocity) → Airy₂. Recordable distinction.
3. Caveat: all of the above is the **edge/operator** caricature `D + C`; the genuine
   dynamical folded-PDE simulation (escape front in space–time) is the next, heavier
   confirmation, and where any non-universal corrections would show up.
