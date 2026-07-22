# MMO Phase 1.5 — folded-node μ, parametrisation fix, and an honest α

**Status:** Phase 1.5 (Tasks A′–C′ of `MMO_FHR_PLAN.md`) executed. Three results,
two of which sharpen the plan and one of which corrects it.

Companion scripts: `mmo_fhr_foldednode.py` (Task B′), `mmo_fhr_staircase.py`
(Tasks A′/C′, from Phase 1). Figure: `figures/mmo_fhr_foldednode.png`.

---

## 1. Parametrisation correction — the plan's §0 FHR does not make MMOs

`MMO_FHR_PLAN.md` §0 writes a *canonical* FHR in the Rinzel form

```
v' = v − v³/3 − w − y + I,   w' = ε(v + a − b w),   y' = δ(c v + d − y),
a=0.7, b=0.8, ε=0.08, c=1.0, d=0.8, δ∈[0.001,0.01].
```

Tested directly (`mmo_fhr_foldednode.py` §1): this system produces **pure
spiking — zero SAOs — at every δ ∈ {0.005, 0.01, 0.02, 0.04} and I ∈ {0.5…1.3}**.
The folded node exists analytically there (see §2), but the **global return
re-injects the trajectory outside its funnel**, so it never executes the small
rotations. (At larger I the folded singularity is even a folded *focus* — no real
canards at all.) So the plan's §0 parameter set is **not** an MMO regime; it is an
elliptic burster / pure spiker.

The model that **does** produce the MMO staircase (Phase 1b) is the variant

```
v' = v − v³/3 − w + y + I,   w' = ε(v + a − b w),   y' = ε δ (c − v),
a=0.7, b=0.8, ε=0.08, δ=0.2, I=0.30,  c the bifurcation control (≈ −0.7…−0.96).
```

**Recommendation:** replace `MMO_FHR_PLAN.md` §0's equations/parameters with this
working variant (or re-derive a Rinzel-form parameter set that actually feeds the
funnel — a bounded search, but the §0 numbers as written do not). Everything
downstream (folded-node analysis, σ_pq) is unaffected by which FHR variant is
used, as long as it is one with genuine MMOs.

## 2. Folded-node μ — verified, but the SAO formula is a bound, not an equality

**The reduction (Task B′).** Treating (w, y) as the two slow variables, the
desingularised reduced flow on the critical manifold `w = v−v³/3+y+I` has, at the
fold `v=−1`, a folded singularity whose Jacobian is

```
J = [[1+δ, −b], [2 b δ (c+1), 0]],   tr = 1+δ,   det = 2 b δ (c+1).
```

For the working FHR this is a **folded node** (real eigenvalues, det>0) across the
whole MMO band, with eigenvalue ratio `μ = λ_weak/λ_strong` falling toward 0 as
`c → −1` (det → 0).

**Wechselberger test.** The plan (§2.1) asserts `#SAOs ≈ (1−μ)/(2μ)`. Computing
μ(c) and the *observed* maximum SAOs per spike:

| c | μ | (1−μ)/(2μ) (bound) | observed max SAOs |
|---|---|---|---|
| −0.72 | 0.071 | 6.5 | 1 |
| −0.78 | 0.054 | 8.7 | 1 |
| −0.84 | 0.038 | 12.5 | 3 |
| −0.90 | 0.023 | 21.0 | 7 |
| −0.94 | 0.014 | 36.0 | 16 |

The observed SAO count **grows in lockstep with the bound** as μ → 0 — confirming
that the folded-node μ genuinely governs the SAO structure — but it sits
**strictly below** `(1−μ)/(2μ)` (≈15 % of the bound at large μ, rising to ≈45 % at
the deepest point). This is exactly right: `(1−μ)/(2μ)` is the **maximal**
rotation number (the number of secondary canards), realised only as the global
return injects deep into the funnel. **The plan's "≈" should be "≤"** — the SAO
count per spike is bounded by `(1−μ)/(2μ)` and fills up toward it as the
re-injection moves toward the strong canard. This is a sharper, correct statement
of the deterministic structure than the plan's equality.

## 3. Honest α — Δ_pq ~ 1/q^α with α ≈ 1.0, not the circle-map [2,3]

> **⚠ SUPERSEDED by `MMO_ALPHA_DERIVATION.md`.** The α ≈ 1.0 below was measured
> from `ρ=p/q` samples within a tolerance on a coarse grid — an artifact. The
> clean SAO-count scaling (`s_obs ~ (c+1)^{−1.74}`) gives **realised α ≈ 1.6**,
> and the folded-node geometry (`μ ∝ (c+1)` + Wechselberger `s_max=(1−μ)/(2μ)`)
> gives a **derived ceiling α = 2**. So **α ≠ 1** either way. Retain the text
> below only as the record of the coarse first pass.

Phase 1's coarse fit gave α ≈ 1.1. A finer sweep (Δc = 0.0025, plateaus resolved
at q = 2,3,4,5: 1/2, 2/3, 1/3, 3/4, 1/4, 2/5, 3/5, 1/5) gives **α ≈ 0.96**. So
across both resolutions `α ≈ 1.0`, robustly **below** the standard critical
circle-map value `α ∈ [2,3]` the plan inherited from the (now-retired) 2D
Boyland/Glass–Mackey framing.

This is the "substantive result" Task C′ anticipated, and it is *consistent with
the redirect*: the FHR staircase is **not** a critical circle map — its plateaus
are organised by the folded node's **secondary canards** (`2k+1` canards for
`μ ∈ (1/(2k+3), 1/(2k+1))`, Wechselberger 2005), a different mechanism with
different width scaling. Pinning α from first principles is precisely Phase 2's
return-map deliverable; the numerics say it should target ≈ 1, not [2,3].
(Caveat: still resolution-limited — the high-q widths are a few grid spacings —
so α ≈ 1.0 is a Phase-1.5 measurement, not a 2-figure result.)

## 4. What this does to the chapter

| Plan claim | Phase-1.5 status |
|---|---|
| §0 canonical FHR is an MMO regime | **wrong as written** — pure spiking; use the working variant (§1) |
| folded node is the right engine | **confirmed** — μ governs SAO structure across the band |
| `#SAOs ≈ (1−μ)/(2μ)` | **corrected to `≤`** — it is the upper bound (secondary-canard count) |
| `Δ_pq ~ 1/q^α`, α ∈ [2,3] | **α ≈ 1.0**, below [2,3]; folded-node not circle-map |

None of this derails the chapter — it grounds it. The folded-node framework is
the right one; the deterministic deliverables (Phase 2) are now sharper: derive
the SAO-count bound and the α ≈ 1 plateau scaling from Wechselberger's return map,
not from a circle-map analogy.

**Sharper supervisor question (refines `MMO_FHR_PLAN.md` §5):** the FHR MMO
staircase has α ≈ 1, not the circle-map [2,3], because it is folded-node-organised
(secondary canards) rather than a critical circle map. Does Wechselberger's
return map predict α ≈ 1 directly, and is `(1−μ)/(2μ)` the right object to anchor
the noise σ_pq dissolution on (the bound) or should σ_pq be indexed by the
realised rotation number instead?

## 5. Reproduce

```
python3 regime-tests/mmo_fhr_foldednode.py   # parametrisation check + μ test (Task B')
python3 regime-tests/mmo_fhr_staircase.py    # staircase + α (Tasks A'/C')
```
Outputs in `regime-tests/results/mmo/` and `figures/`.
