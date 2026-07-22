# Inverter hardening — off-grid fitting, near-threshold class weighting, colored-noise axis

*Closes the three "still open" items in `INVERSION_MVP_ROADMAP.md` (Phase 1+/2). Code:
`harden_inverter.py` (imports the existing atlas machinery in `bifurcation_classifier`; nothing
downstream modified). Figure: `figures/harden_inverter.png`. Tags: **[R]** proved / **[N]** numerical
/ **[H]** heuristic.*

```
python3 harden_inverter.py coloratlas   # build + cache the (ν, τ_s) colored atlas (~10 s, one-time)
python3 harden_inverter.py harden       # all three numerical tests (≈18 s)
python3 harden_inverter.py fig          # -> figures/harden_inverter.png
```

## (A) Off-grid (interpolated) fitting — grid-snap removed

The atlas stores a mean-normalised quantile function `Q(p)` on a **coarse** parameter grid; the old
`classify3`/`fit_adapt` returned `argmin` over those nodes, so `param̂` snapped to a node. Because
`Q(p)` is smooth and monotone in the class parameter, we **linearly interpolate the quantile curves
between nodes** and minimise `W1` on a dense continuous axis (`fit_interp_1d`; bilinear for the 2D
`(ν,b)`/`(ν,τ_s)` grids, `_dense_2d`). `param̂` is now continuous.

- **[N]** On SNIC data drawn at `ν` deliberately *between* atlas nodes (node spacing ≈ 0.28), mean
  `|ν̂−ν|` falls **0.079 → 0.031 (≈60% reduction)**; recovered `ν̂` lands on the true value instead
  of the nearest node (e.g. true `ν=2.05`: grid `1.892` → off-grid `2.008`).
- **[H]** Validity rests on `Q(p)` being smooth/monotone in `p` between adjacent nodes (true for the
  one-parameter edge families); linear interpolation of the quantile curves is the cheapest estimator
  consistent with that and is what removes the snap. A spline would refine it but is not needed at MVP
  atlas resolution.

## (B) Near-threshold weighting for the class call

**Justification — separability is a near-threshold resource. [N]** Pooling large samples of SNIC
(Type-I) spikes and measuring the `W1` distance to the **Hopf** atlas (the "room for error" before a
sweep is misread as Type-II), that distance **collapses with drive**:

| drive `ν` | 0.2 | 0.5 | 0.9 | 1.4 | 2.0 | 2.8 |
|---|---|---|---|---|---|---|
| ISI CV | 0.48 | 0.42 | 0.32 | 0.25 | 0.20 | 0.17 |
| `W1` to Hopf | 0.066 | 0.050 | 0.052 | 0.045 | 0.029 | **0.017** |

As drive rises the Type-I ISI becomes regular (CV→0) and its shape converges onto Hopf's near-symmetric
jitter, so a **high-drive sweep sits on the SNIC↔Hopf boundary and flips under sampling noise**. The
informative class evidence lives **near threshold**.

**Method.** `aggregate_class` pools per-sweep soft class evidence across a cell's sweep set, weighting
each sweep by a near-threshold factor `1/(1+(ν̂/ν₀)²)` (`ν₀=0.8`) times a sample-size weight.

**What it does — and doesn't — buy (Monte Carlo, 80 Type-I cells, ~26 ISIs/sweep). [N]**

| sampling | metric | uniform | near-threshold |
|---|---|---|---|
| balanced (drives spread across edge) | correct Type-I call | 99% | 98% |
| | mean P(Type-I) | 0.71 | **0.74** |
| high-drive (mostly above rheobase) | correct Type-I call | 76% | 76% |
| | mean P(Type-I) | 0.56 | **0.63** |

The honest read: cell-level pooling is **already robust on balanced sampling** (both ≈100%), so the
weighting's effect there is **better-calibrated confidence**, not higher accuracy. Its value shows in
the realistic **high-drive-dominated** regime, where uniform voting sits at a near coin-flip posterior
(P=0.56) and weighting pulls it to **P=0.63** by leaning on the few informative near-threshold sweeps.
It does **not** rescue the binary call rate when a cell is sampled *only* far above threshold (there is
simply little class signal left) — that is a sampling-design limit, not a fixable estimator flaw.
*Takeaway for real data: record and weight near-rheobase sweeps for the class call.*

## (C) Colored-noise confound axis

Real input noise is correlated. `colored_snic_intervals(ν, τ_s)` drives the SNIC QIF with OU noise
(the `qif_validation` convention, → white as `τ_s→0`); a 2D **(ν, τ_s)** atlas + joint `fit_colored`
sits alongside the existing `(ν, b)` adaptation axis.

- **[N] ν is recovered; the white-only atlas is biased.** On colored data (`τ_s=0.6`), fitting against
  the white atlas gives mean `|ν error| = 0.092`; the **joint (ν, τ_s) fit cuts it to 0.052 (43%)**.
- **[N] Coloring is detected.** `τ̂_s` rises monotonically with true `τ_s` (`τ_s = 0.05/0.3/0.6/1.0
  → τ̂_s ≈ 0.37/0.20/0.44/0.68`), so the axis flags colored noise.
- **[N] But `τ_s` is only weakly identified — a real degeneracy.** The `W1` valley is **≈8× steeper
  along ν than along τ_s** (sensitivity 0.168 vs 0.020). So `ν` is well-constrained while `τ̂_s` is
  biased and noisy — the colored-noise analogue of the documented low-`ν` `(ν,b)` adaptation
  degeneracy. The method **reports that coloring is present and corrects `ν` for it**, rather than
  precisely measuring `τ_s` — which is the behaviour an early-warning tool actually needs.

## Net effect on the inverter

`ν̂` is now continuous (no grid-snap) and more accurate under both confounds in the atlas; the
class call is computed from a cell's most informative (near-threshold) sweeps with a calibrated
posterior; and colored noise joins adaptation as a fitted/flagged confound axis. Remaining honest
limits, unchanged in character: the class signature is intrinsically weak far above threshold, and the
second confound parameter (`b` at low `ν`, `τ_s` generally) is weakly identified — both **surfaced**
by the fit rather than silently absorbed.
