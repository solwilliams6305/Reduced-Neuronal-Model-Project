# Notes — the cusp rung is Weber, and the Airy→Weber crossover

_June 2026. Companion to `PROBLEM_STATEMENT_NOISY_CUSP.md` and `CREATIVE_DIRECTIONS.md`.
Figures in `coupled-atlas/figures/`: `weber_vs_pearcey.png`, `crossover_airy_to_weber.png`.
Tags: [PROVED] established/citable · [NUMERIC] shown here by simulation · [HEURISTIC] argued, not proved · [OPEN] needed for a theorem._

---

## Part 1 — the cusp rung's identity: Weber, not Pearcey

**Result.** The genuine Kristiansen–Pedersen coupled-FHN cusp escape is **Weber-class
(2nd-order, quadratic potential), not Pearcey (3rd-order).**

- **[PROVED]** (K-P, arXiv:2202.12027): the slow-fast cusp's reduced inner equation is the
  **Weber / parabolic-cylinder** equation, SAOs counted by parabolic-cylinder zeros — a
  *2nd-order* object. The literal diffraction cusp (Pearcey) is *3rd-order*; the slow flow
  collapses it to Weber. So at the level of the reduced equation, "Weber not Pearcey" is a theorem.
- **[NUMERIC]** (`weber_vs_pearcey.py`): the standardised KP cusp escape matches the
  Weber-class inner law and rejects Pearcey:

  | law | skew | exkurt | KS to KP cusp |
  |---|---|---|---|
  | fold / Airy (TW) | +0.21 | +0.03 | 0.108 |
  | **Weber-class (2nd, quadratic)** | **+0.64** | **−0.19** | **0.056** |
  | Pearcey (3rd order) | +0.12 | +0.08 | 0.123 |
  | KP cusp escape (orient.) | +0.61 | −0.63 | — |

  Skew matches Weber almost exactly (0.61 vs 0.64) and is double the Pearcey/fold value;
  excess kurtosis has the Weber sign (negative); KS to Weber (0.056) is >2× closer than to
  Pearcey (0.123); and the KP cusp's *compressed lower tail* (1%-ile −1.6) tracks Weber (−1.7),
  not fold/Pearcey (−2.2). **Deviations:** kurtosis *magnitude* (−0.63 vs −0.19) and the far
  right tail — the full-model escape carries outer (Gaussian) corrections, and the quadratic-
  turning proxy is not the exact even-well parabolic cylinder.
- **[HEURISTIC]** orientation: the full-model peel-off (w at fold release) is oriented opposite
  to the inner peel-off level (early escape = small w = left tail, vs inner early escape = large
  Y_node = right tail), so the comparison negates the KP sample. The parabolic-funnel mechanism
  is what supplies the extra spread and skew.

---

## Part 2 — the Airy → Weber crossover (two folds merging)

**Result.** As coupling drives the two symmetric-mode folds together, the peel-off edge class
crosses over **continuously** from fold/Airy/TW to cusp/Weber. This is the stochastic shadow of
a classical uniform-asymptotics fact.

- **[PROVED]** (Olver, *Asymptotics and Special Functions*, coalescing-turning-points chapter):
  an **isolated** turning point of u″=(V−·)u gives the **Airy** function; **two coalescing**
  turning points give the **parabolic-cylinder / Weber** function. With V≈Y²−Δ² (turning points
  at ±Δ; `crossover_airy_to_weber.py` panel A), Δ≫ℓ is "isolated" (Airy) and Δ≲ℓ is "coalesced"
  (Weber), where ℓ is the inner blow-up scale. The effective turning order k_eff runs 1→2 as Δ/ℓ:∞→0.
- **[NUMERIC]** the inner crossover: sweeping the turning order, the peel-off skew rises smoothly
  **0.25 (k=1, Airy/TW) → 0.61 (k=2, Weber)** — and that endpoint equals the KP cusp |skew| 0.61
  from Part 1 (panel B). The **physical** crossover (`crossover_fold_to_cusp.py`, panel C):
  sweeping coupling from two folds into the symmetric cusp swings the peel-off skew and amplifies
  the spread **×3.6**, with the transition at the cusp onset **g ≈ −0.05 … −0.11** — i.e. where
  Δ(g) ~ ℓ. Inner and physical crossovers agree in direction and endpoint.
- **[HEURISTIC]** the identification k_eff(Δ/ℓ) and the map g ↦ Δ(g) are argued from the geometry
  (the symmetric-mode fold separation shrinks toward the cusp), not derived.

### What would make it a theorem [OPEN]

1. **Δ(g):** derive the symmetric-mode critical-manifold fold separation as a function of coupling
   from the coupled-FHN geometry; show Δ(g_cusp)=0.
2. **Coupled cusp blow-up:** the Krupa–Szmolyan / JKK-style blow-up at the codim-2 cusp, with the
   rescaling that resolves Δ against the inner scale ℓ — i.e. the *uniform* Airy↔Weber connection
   realised inside the blow-up charts.
3. **Noise through the charts:** carry additive noise through both regimes (Berglund–Gentz tubes);
   show the inner exit measure → TW for Δ≫ℓ and → the Weber/parabolic-cylinder edge law for Δ≲ℓ.
4. **The crossover as a uniform stochastic connection:** a stochastic analogue of Olver's
   coalescing-turning-point asymptotics — the actual new theorem.

The deterministic uniform connection (1–2 without noise) is classical and should be the tractable
first rigorous step; the stochastic version (3–4) is the genuine contribution and rests on the
coupled blow-up that is the project's standing analytic gap.

---

## One-line status

Cusp rung = **Weber** (proved by K-P; numerically matched here, KS 0.056 vs Pearcey 0.123).
Fold↔cusp = a **continuous Airy→Weber crossover** (classical uniform asymptotics; reproduced
inner *and* physically; theorem pending the coupled cusp blow-up).
