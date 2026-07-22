# MMO keystone — the canard/folded-node crossover (Phases A+B, normal form)

**Status:** inheritance **confirmed at the normal-form level**, with a refinement.
The folded-node escape continuously becomes the canard escape as the folded node
degenerates to a folded saddle-node (μ → 0): η_*(μ) → the canard plateau with a
**shared prefactor C_q** (flat ≈ 2.0 for μ ≲ 0.1 in these units). The refinement:
there is **no μ^{3/2} in the local normal-form escape** — η_* is set by the
constant drift, so the μ^{3/2} geometric factor of `MMO_NOISE.md` is a
**global-return (funnel-filling) effect**, not the local escape. That is exactly
the universal-local + model-specific-global split the rest of the chapter found.

Phases A+B only (normal-form analysis + SDE numerics); full-FHR Phases C–D are a
separate follow-on. Companion: `mmo_crossover_nf.py`, `figures/mmo_crossover_nf.png`.

---

## 1. Question

The canard chapter's σ_* = C_q√ε·λ^{1/2} (C_q ≈ 2.8 normal form / 8–10 full FHN) and
the MMO chapter's σ_*(μ) = C_q√ε·μ^{3/2} (C_q ≈ 8–10, the full-FHN canard value)
share a prefactor. Is that a coincidence, or is the canard law the **μ → 0 limit**
of the folded-node law? Phases A+B test it at the universal (normal-form) level.

## 2. A1 — the normal-form bridge (analytic)

Folded-node K2 normal form (Wechselberger): `dV = (V²−W)dT + η dB`,
`dW = (μV − (1+μ)/2)dT`. As **μ → 0** the rotation term `μV` vanishes and
`dW/dT → −1/2` — exactly the canard Krupa–Szmolyan form `dV=(V²−W)dT+η dB,
dW=−λ dT` with **λ = 1/2**. The noise enters identically as η = σ/√ε in both. So
C_q is *the same object in the same coordinates*; only the geometric factor
(λ^{1/2} ↔ μ^{3/2}) can differ. (Convention caveat: the `−(1+μ)/2` constant is from
project notes — the *structural* μ→0 ⇒ constant-drift fact is convention-robust and
is all A1 uses; the μ_c prefactor below is not, so it is not relied on
quantitatively. Verify the constant against Wechselberger 2005 §4 before C–D.)

## 2′. A2 — crossover scale

Balancing the rotation term against the constant drift over the canard window,
`μ·V_* ~ drift` with `V_* ~ λ^{1/3} ≈ 0.79`, `drift ≈ 1/2`, gives **μ_c ≈ 0.6**.
Below μ_c the dynamics is a single canard passage (canard-like escape); above it
the funnel forms. The swept μ ≤ 0.5 therefore sit at/below μ_c — in the
canard-like regime — exactly as the numerics show.

## 3. B — normal-form SDE numerics

Reusing the canard escape criterion verbatim (start V = −√W on the attracting
branch, record W_hit at the first V-crossing of V_cross = 1, R_hit = W_hit/W_*,
W_* = λ_eff^{2/3}, λ_eff = (1+μ)/2; escape threshold = where median R_hit crosses 0;
Θ = η/√λ_eff, C_q ≡ Θ_crit), sweeping μ from 0 (pure canard) to 0.5:

| μ | λ_eff | C_q = Θ_crit | η_* | μ^{3/2}·C_q (ref) |
|---|---|---|---|---|
| 0.00 | 0.500 | 2.01 | 1.42 | — |
| 0.02 | 0.510 | 1.99 | 1.42 | 0.006 |
| 0.05 | 0.525 | 2.02 | 1.47 | 0.022 |
| 0.10 | 0.550 | 2.13 | 1.58 | 0.064 |
| 0.20 | 0.600 | 2.29 | 1.77 | 0.180 |
| 0.50 | 0.750 | 2.51 | 2.17 | 0.711 |

**Read-off:**
- **η_*(μ) → the canard plateau** (1.42) smoothly as μ → 0 — the two laws meet.
- **C_q is flat ≈ 2.0 for μ ≲ 0.1** (the deep canard regime) and rises mildly to
  2.5 by μ = 0.5 (onset of funnel physics as μ → μ_c ≈ 0.6). Total drift 26%, all
  at the large-μ end; at small μ it is constant. **Shared prefactor confirmed.**
- **η_*(μ) does NOT follow μ^{3/2}** (the dotted reference sits 3–200× below the
  data and far steeper). η_* is set by the constant-drift canard scale
  `≈ C_q√λ_eff`, not the rotation term.

(Absolute-value note: with these faster settings the canard anchor gives C_q ≈ 2.0
vs the canard chapter's 2.8 at finer dt/longer T — the *absolute* C_q is
settings-dependent; the **constancy across μ and the μ → 0 continuity** are the
convention-robust results, and the μ-sweep is anchored to its *own* μ = 0 canard.)

## 4. Verdict

**The keystone inheritance holds, as a continuous limit:** the folded-node escape
*is* the canard escape with a shared prefactor C_q, recovered exactly as μ → 0.
This upgrades the "C_q ≈ 8–10 in both chapters" from a measured coincidence to a
continuous μ → 0 limit at the normal-form level.

**Refinement (a real result):** the μ^{3/2} geometric factor is **not** in the
local normal-form escape — that escape is canard-like (set by the constant drift)
across the whole μ-sweep. So μ^{3/2} is a **global-return / funnel-filling**
property (the f(c) of `MMO_K2.md`), not local escape physics. This is the same
universal-local (canard escape, shared C_q) + model-specific-global (funnel-filling
f(c)) decomposition the deterministic side found — now mirrored on the noise side.
The chapter's cross-chapter unification is therefore: **one local escape mechanism
(canard, C_q), two global geometries (single fold ↔ folded-node funnel).**

## 5. What's established / open

| Piece | Status |
|---|---|
| μ→0 ⇒ canard normal form (λ=1/2) | **derived** (A1, structural, convention-robust) |
| crossover scale μ_c ≈ 0.6 | derived (A2; prefactor convention-dependent) |
| η_*(μ) → canard plateau, C_q shared (≈const, μ ≲ 0.1) | **validated** (B) |
| μ^{3/2} is global, not local escape | **shown** (B; η_* canard-like, no μ^{3/2}) |
| absolute C_q = 2.8 (vs 2.0 here) | settings-dependent; finer dt/longer T to match |
| full-FHR crossover σ_*(control) through the FSN, C_q ≈ 8–10 | **open — Phases C–D** (separate handoff; Krupa–Wechselberger 2010 FSN scaling) |

## 6. Reproduce

```
python3 regime-tests/mmo_crossover_nf.py
```
Writes `data/mmo_crossover.npz`, `results/mmo/mmo_crossover.txt`,
`figures/mmo_crossover_nf.png`.
