# Phase-2 prerequisite: does Wechselberger's 2-slow folded node apply at finite δ?

**Status:** checked. **Verdict: YES — not a blocker.** The working FHR at δ=0.2
(and at δ=1.0) is a genuine 1-fast/2-slow system, the folded singularity is a
real folded node, and MMOs exist at the cleanest 2-slow limit δ=1.0. The
α-ceiling=2 result is δ-robust. The plan's aspiration toward δ→0 is the
*three-timescale* regime and is unnecessary — it would move **out** of
Wechselberger's reach.

Companion script: `mmo_fhr_timescale.py`. Figure: `figures/mmo_fhr_timescale.png`.

---

## The concern

`MMO_PHASE1_5.md` flagged that the μ-reduction treats w and y as co-equal O(ε)
slow variables, but the working model has `y' = εδ(c−v)` with δ=0.2 — "not
infinitesimal." Does Wechselberger's 2-slow folded-node theorem apply, or is the
system sliding into the three-timescale regime where Krupa–Popovic–Kopell (2008)
is needed instead?

## The resolution

**Timescale bookkeeping.** A 2-slow folded node requires both slow variables at
O(ε). Here w-rate = ε and y-rate = δε. **δ is an O(1) *ratio*** — so both rates
are O(ε) and the system is genuinely 2-slow for any δ = O(1). The concern had the
direction backwards: 2-slow does **not** want δ infinitesimal; δ = O(1) is exactly
the 2-slow setting. Genuine three-timescale is the **δ → 0** limit (y
asymptotically slower than w), which is what the plan's §0 aspired to with
δ ∈ [0.001, 0.01] — and *that* is the regime that would require Krupa–Popovic–
Kopell, not the working δ = 0.2.

**Numerical confirmation (`mmo_fhr_timescale.py`).**

1. The folded singularity is a **real folded node** (real eigenvalues, μ∈(0,1))
   for every δ ∈ [0.05, 1.0] and c in the MMO band. μ ∝ (c+1) at fixed δ (so the
   α-ceiling = 2 derivation is δ-robust), and μ rises with δ then saturates at
   ≈ 0.16 (c=−0.7) — never reaching O(1), i.e. the node is always in the
   many-rotation (small-μ) part of its range, but always a valid node.

2. **MMOs exist at δ = 1.0** — the cleanest co-equal 2-slow limit (c ≈ −0.93:
   48L/144S). They persist across δ ∈ [0.1, 1.0] at a δ-dependent c. (An earlier
   too-narrow c-scan wrongly suggested δ=1.0 had no MMOs; a broader scan finds
   them at more negative c.) So a *clean* 2-slow folded node with MMOs is
   available at δ=1.0, where Wechselberger 2005 applies with no three-timescale
   caveat at all.

## Recommendation for Phase 2

Launch the K2-chart return-map calculation at **δ = 1.0** (cleanest 2-slow
folded node, μ ≈ 0.03–0.04 in the MMO band, MMOs present), with δ = 0.2 as a
cross-check. Both are legitimate 2-slow folded nodes. **Do not** push δ → 0 (the
plan's §0 value): that exits the 2-slow regime into three-timescale territory and
would require different machinery for no benefit — MMOs and the folded node are
already present and clean at δ = O(1).

Caveat retained: μ is small (≈0.03–0.16) throughout the MMO band, so the node sits
in the many-secondary-canard part of its range (near, but not at, the
folded-saddle-node μ→0 boundary). Wechselberger's theorem holds for all μ∈(0,1),
but the small-μ regime is numerically delicate (many accumulating secondary
canards) — expect the K2-chart matching to need care there, as the plan's
hazard §6.4 anticipated.

## Reproduce

```
python3 regime-tests/mmo_fhr_timescale.py
```
Writes `results/mmo/mmo_fhr_timescale.txt` and `figures/mmo_fhr_timescale.png`.
