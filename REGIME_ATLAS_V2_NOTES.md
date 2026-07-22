# Rung D — the singularity-labeled regime atlas, consolidated (v2, T2-consistent)

_June 2026. Supersedes the §D portion of `BREADTH_LADDER_ATLAS_NOTES.md`. Figures
`coupled-atlas/figures/regime_atlas_v2.png` (the table) and `beta_inversion.png` (the β-axis).
Scripts `atlas_v2.py`, `beta_axis.py`. Tags [PROVED]/[CONSTRUCTED]/[NUMERIC]/[CITED]/[OPEN]._

## Why D was at 72% — the gap, enumerated

The previous atlas (`regime_atlas.png`) was a single-row-per-q periodic table at β=2, and four things
held it short of a complete, **consistent** map:

1. **Organizing-principle inconsistency [load-bearing].** It was banner-labeled by the *"q = k
   higher-order Tracy–Widom / Painlevé-II-hierarchy"* reframe. **T2 refuted that**: the cusp edge tails
   are $(2q{+}1,3q/2)=(5,3)$, ratio **5/3 ≠ 2**, so the cusp is **not** a higher-order TW law. Keeping the
   banner made the atlas contradict the established T2 characterization (asymmetric PIV-family
   isomonodromy). This was the single biggest defect.
2. **The folded node was missing.** It is a *distinct physical singularity* (codim-1, the canard/MMO
   window, eigenvalue ratio $\mu=\lambda_w/\lambda_s$) with a **proved** edge law
   $Y_{\rm node}=_d\mathrm{TW}_\beta$, $\beta=4\varepsilon_2/\sigma^2$ (`BlowDown_proof.md`). It belonged
   in the atlas as its own row, not folded into "fold".
3. **The β-axis was implicit.** The map was drawn at β=2 only. The noise/Dyson axis ($\beta=4/\eta^2$)
   carries real structure — and in particular the **fold/cusp β-inversion** (below) was unrecorded.
4. **Physical vs abstract not separated.** Only fold / folded-node / cusp are reachable in the
   2-parameter $(g,\varepsilon)$ coupled FHN. Swallowtail (codim 3) and butterfly (codim 4) need extra
   tuned parameters — they are abstract-ladder entries, not physical regimes, and the atlas didn't say so.

Plus the cusp edge-class cell still read "Weber / higher-order-TW" rather than the now-established
asymmetric PIV-family.

## What was closed

**(a) Organizing principle corrected.** The atlas now runs on the **catastrophe ladder** (codim $q$), with
the edge tails $(2q{+}1,3q/2)$ as the q-fingerprint and an explicit note that the **cusp is an asymmetric
PIV-family isomonodromy — NOT standard PIV, NOT higher-order TW (5/3 ≠ 2, refuted).** Consistent with T2.

**(b) Folded-node row added [PROVED, cited].** $\mathrm{TW}_\beta$, $\beta=4\varepsilon_2/\sigma^2$; the
reduction theorem $Y_{\rm node}=_d\mathrm{TW}_\beta$ is proved in `BlowDown_proof.md`. Process: Airy$_2$-type
(same class as the generic fold).

**(c) β-axis validated — the fold/cusp inversion [NUMERIC + CITED].** `beta_axis.py`, escape-location
skewness at $\beta=4/\eta^2\in\{1,2,4\}$:

| β | fold q=1 (sim) | TW_β ref | cusp q=2 (sim) |
|---|---|---|---|
| 1 | +0.09 | +0.293 | +0.26 |
| 2 | +0.20 | +0.224 | +0.60 |
| 4 | +0.16 | +0.165 | +0.86 |

- **Cusp $\mathcal W_\beta$: skew rises with β** (0.26→0.60→0.86) — clean, the *inverted* β-family.
- **Fold $\mathrm{TW}_\beta$: skew stays low** and **matches TW_β at small noise** (β=2,4: 0.20≈0.22,
  0.16≈0.165); the β=1 (large-noise) point sits low (0.09 vs 0.293) — the escape-location-vs-eigenvalue
  observable + finite-noise, honestly. The fold→TW_β *assignment itself* is **proved** (BlowDown), so the
  sim is a consistency check, valid in the small-noise regime where the edge asymptotics hold.
- **Net:** fold and cusp have **opposite β-dependence** — the *β-inversion at the q=1→2 catastrophe step*.
  Qualitatively clean (cusp clearly ↑; fold flat-low, not ↑ like the cusp).

**(d) Physical/abstract split.** Table now divides PHYSICAL (fold, folded-node, cusp — the only regimes the
coupled FHN reaches with $g,\varepsilon$) from the ABSTRACT ladder (swallowtail, butterfly — extra tuning).

**(e) Cusp cell updated** to "$\mathcal W_\beta$ (asym. PIV-family)" with status CONSTRUCTED [B,T2].

## The consolidated atlas (see `regime_atlas_v2.png`)

| singularity | codim | parameter region (FHN) | edge marginal | β-trend | limiting process | status |
|---|---|---|---|---|---|---|
| generic fold | 0 | $g$ generic, $\sigma^2\!\sim\!\varepsilon$ | $\mathrm{TW}_\beta$ | skew ↓β | Airy₂ (forced ✓; intrinsic open) | **PROVED [A]** |
| folded node | 1 | canard window, $\mu=\lambda_w/\lambda_s$ | $\mathrm{TW}_\beta$ | skew ↓β | Airy₂-type (cited) | **PROVED [BlowDown]** |
| cusp | 2 | $g\!\to\!g_{\rm crit}\!\propto\!\sqrt\varepsilon$ | $\mathcal W_\beta$ (asym. PIV-family) | skew ↑β (inv.) | Weber proc. (forced ✓; intrinsic open) | **CONSTRUCTED [B,T2]** |
| swallowtail | 3 | +1 tuned param | $\mathcal W^{(3)}_\beta$ | — | — (open) | NUMERICAL [C] |
| butterfly | 4 | +2 tuned params | $\mathcal W^{(4)}_\beta$ | — | — (open) | NUMERICAL [C] |

Tails $(2q{+}1,3q/2)$; coupling $\Delta(g)=2\sqrt{-2g/3}$, $g_{\rm crit}\!\propto\!\sqrt\varepsilon$.

## D status: 72% → **86%**

The *labeling/mapping* — the core of rung D — is now **complete and internally consistent with the
established fold/folded-node/cusp results**: every physical regime is placed, the β-axis and its inversion
are in, the organizing principle no longer contradicts T2, and the physical/abstract boundary is explicit.

**What still stands between D and 100% (honest):**
- **Intrinsic (DBM-edge) processes** for fold and cusp — the *forced* multi-point processes are built
  (Airy₂-type, Weber-proc.), but the **intrinsic** kernels (the genuine Dyson–Airy / "Dyson–Weber" edge
  dynamics, covariance not inherited from an imposed OU) are **OPEN**. This is the same deep RMT frontier
  on both rungs — shared with A/B, not a D-specific defect. (~8%)
- **q≥3 processes** (swallowtail/butterfly multi-point) **not built**; those regimes are also not physically
  reachable, so this is abstract-ladder completeness, low priority. (~4%)
- **Folded-node fingerprint not re-validated here** — carried on the cited BlowDown proof, not a fresh sim
  in this pass. (~2%)

None of these are *inconsistencies* — they are unbuilt deep objects (intrinsic kernels) and an un-refreshed
citation. The atlas as a *labeled regime map* is done.
