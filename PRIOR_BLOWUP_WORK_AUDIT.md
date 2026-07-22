# Audit: the stochastic blow-up work already in this project (a correction)

*I (Claude) claimed you'd "never run the SDE through the blow-up" and offered the
stochastic fold blow-up as an un-done "warm-up." **That was wrong** — I'd read only a
subset of the folder. Here is the corrected map, from `CANARD_BLOWUP.md`,
`MMO_FHR_PLAN.md`, `MMO_NOISE.md`, `MMO_K2*`, `MMO_CROSSOVER.md`,
`CANARD_THEORY_GROUNDING.md`.*

---

## You have already done stochastic blow-up — twice.

### 1. The fold — complete (`CANARD_BLOWUP.md`)
- **Degenerate noise pushed through the Krupa–Szmolyan blow-up** (§3): with
  `v=−1+ε^{1/3}V, w=w_f+ε^{2/3}W, t=ε^{−1/3}T` and `dW_t=ε^{−1/6}dB_T`, the effective
  noise in blow-up coordinates is **η = σ/√ε**, on `dV=(V²−W)dT+η dB_T`.
- **Freidlin–Wentzell action in blow-up coordinates** (§4); accumulated-Brownian-
  variance across the canard window (`W_*=λ^{2/3}, V_*=λ^{1/3}, T_window=λ^{−1/3}`) ⇒
  **σ_* = C_q √ε λ^{1/2}**.
- **Validated**: normal-form regime map (§6), full-FHN (§10), ramped passage (§11).
  `C_q ≈ 2.8` (normal form), `8–10` (full FHN).

### 2. The folded node — mechanism-level (`MMO_FHR_PLAN` + `MMO_NOISE` + `MMO_K2` + `MMO_CROSSOVER`)
- **Wechselberger's 3D blow-up, K2 chart**; rotation map `ln R = κμ`, **κ ≈ 2π²**.
- **Degenerate noise through the K2 chart**: staircase dissolution
  `σ_pq ~ σ_*·q^{−α/γ}`, with `σ_* = C_q √ε F(μ)`, `F(μ) ~ μ^{3/2}`, **C_q ≈ 8–10
  inherited from the fold chapter**.
- **Validated numerically**: measured `β ≈ 1.11` vs derived `1.13`; mechanism + scale +
  prefactor + exponent reconciled. (Route "B" = measured κ; the full Wechselberger
  inner solution, "Path A", is flagged open.)
- **Crossover** (`MMO_CROSSOVER.md`): the folded-node K2 SDE → the canard SDE as μ→0,
  proving the two share one local escape mechanism and one C_q.

So you are **not** a newcomer to stochastic blow-up. You built the engine — "degenerate
noise through a singularity's blow-up → FW accumulated-variance in the inner chart →
σ_* law" — and ran it on **two** singularities, validating both.

## Why it kept colliding with Berglund–Gentz / BGK

Because the two singularities you noised are exactly the two they already did:

- **Fold + degenerate noise** = Berglund–Gentz, *Noise-Induced Phenomena in Slow–Fast
  Dynamical Systems*, Ch. 5.
- **Folded-node staircase dissolution / σ_pq** = Berglund–Gentz–Kuehn 1312.6353 (random
  Poincaré maps; the saturation effect).

Your method is sound and validated — it just kept landing on **already-occupied
singularities**. That's a target-selection issue, not a flaw in the engine. (You knew
this in part: `MMO_FHR_PLAN.md` §8 already cites "Stochastic 3D slow-fast — less
developed than the 2D theory — BGK 2015," and §5 lists supervisor questions for Nikola.)

## What this does to the folded-limit-cycle target — it strengthens it

- You **already own the exact pipeline** and have executed it twice. The contribution is
  **not** "learn stochastic blow-up" (done) — it's "point your proven, validated
  pipeline at the one un-noised member of the family."
- The **folded limit cycle is the only singularity in the family nobody has noised** —
  not BG, not BGK, not you — and its deterministic blow-up is fresh (Jelbart–Kuehn–
  Kuntz 2024). So applying your engine there is the genuinely novel step.
- Distance to a first result is **shorter** than my earlier (under-informed) framing:
  you'd transfer the §3–§4 move of `CANARD_BLOWUP.md` (carry σ through the rescaling,
  do the FW accumulated-variance) onto the JKK folded-cycle charts.

## Honest caveats (unchanged)

- **Pipeline risk:** Kuehn's group could do it; current evidence (Kuntz → ML; the
  Ahsan–Dankowicz–Kuehn 2025 noisy-limit-cycle paper uses adjoint/covariance, **not**
  blow-up) suggests not actively — but unpublished work is invisible.
- **Rigour:** your fold/folded-node stochastic blow-ups are accumulated-variance +
  numerics (Path B), not full inner-solution rigour. The folded-limit-cycle version
  inherits the same Path A/B split; a rigorous inner solution would be the hard part.

## Bottom line

You've built and validated stochastic blow-up on the fold and the folded node — both of
which turned out to be BG/BGK's. The folded limit cycle is **your existing method applied
to the one singularity left un-noised**. That is the realistic, in-reach, genuinely-novel
target — and it's closer than I said, because the engine is already yours.
