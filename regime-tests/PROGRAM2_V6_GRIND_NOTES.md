# PROGRAM2 session note — v₆ computed exactly & grid-certified; Borel phase re-pinned to θ ≈ 48° ± 3° (2026-07-15)

_Executes `PROGRAM2_HANDOFF_A_V6_GRIND.md` (route (a) of the endgame). Engine work + compute grind + Borel
re-analysis, 2026-07-10 → 07-15. Primary results below; all runs quad-hybrid transfer engine unless noted._

## Headline results

1. **v₆ = −1.90 ± 0.15** (grid-certified, see ladder below). **The predicted sign flip to v₆ ≈ +2 is REFUTED**
   — the 6-coefficient Darboux extrapolation (§3quinque of `PROGRAM2_CONSOLIDATION.md`) was wrong, exactly the
   kind of unidentifiability artifact suspected in the Milestone-2 post-mortem.
2. **The Borel phase moves DOWN, not up: θ = 48° ± 3°** (single-pair Darboux fit on v₀…v₆ with sharpened
   ladder; shallow residual basin 45–50°, hard exclusion of θ ≥ 53°). The old "54–63°, NOT 45°" verdict is
   **retracted**: with v₆ in hand, θ = 45° (λ₀ inheritance) sits INSIDE the fit basin (resid 0.153 vs min 0.122
   at 48–50°). The λ₀-inheritance hypothesis is revived but not yet pinned — "≈48°, 45° allowed" is the honest
   statement. Borel–Padé [2/4] independently gives the pair at 48.2° (its other orders scatter 57–60° — Padé is
   the noisier tool at 7 coefficients, as expected).
3. **Sharpened ladder anchors:** v₄(∞) = −0.451 ± 0.003 (plateau n=24–32), v₅(∞) = −1.19 ± 0.01 (geometric
   n=16–32; the old −1.1 was ~8% low). Ladder now: v = (0.134, 0.111, 0.104, −0.030, −0.451, −1.19, −1.90).
   Signs +,+,+,−,−,−,− : the negative run is FOUR long — consistent with period-8 (θ=45°) and with θ≈48°;
   inconsistent with period-6 (θ=60°).
4. **First-ever v₇ datum:** v₇(n=10) = −178.6 (pre-asymptotic; pipeline proof + cost calibration only).
   **v₇ is the decisive next discriminator:** single-pair fits predict v₇ ≈ +1.3 (θ=45°) vs +8.2 (θ=50°) vs
   +24 (θ=56°) — far outside error bars of any estimate. Also needed for the two-singularity (pair + real
   instanton) 8-parameter fit.

## The v₆ grid ladder

| n  | v₆(n)      | Var(Y₇) term | note |
|----|------------|--------------|------|
| 10 | +80.942    | +50.855      | pre-asymptotic (dual-engine exact match) |
| 12 | −1.737     | +43.590      | pre-asymptotic (dual-engine exact match; error budget run) |
| 14 | +4.663     | +42.108      | pre-asymptotic oscillation |
| 16 | −0.427     | +40.804      | smooth branch begins |
| 18 | −1.743     | +40.038      | parity-wobble point (2-step) |
| 20 | −1.906     | +39.560      | 4-step ladder |
| 24 | −1.896     | +39.018      | plateau: 20→24 move = +0.011 |

v₆(∞) = −1.90 ± 0.15 (plateau read; wobble bound from n=18–24 spread). Var(Y₇)(∞) ≈ +38.6 ± 0.2
(clean geometric tail, ratio ≈ 0.44 per 4 steps).

## Certification (per program workflow)

- **Independent dual-engine validation:** full v₆ assembly reproduced EXACTLY (every term, every printed digit)
  by two independent moment engines at n=10 and n=12 (generic transfer DP vs quad-hybrid).
- **Boundary/renormalization-convention systematic:** `_v6_breakdown.py` at n=12 — the boundary+order>5
  convention-caveat moment class contributes −0.00101, i.e. **0.1% of v₆**. Negligible.
- **Sanity gates (handoff §"plan" item 3):** Var(Y₇) dominant & positive ✓; (5,9) cross large & negative
  (−21 to −26 across grids) ✓; all 7 terms converge smoothly even where the assembled residue oscillates ✓.
- **Extrapolation methodology validated on v₅** (known answer): the 4-step-grid plateau/geometric fit on
  n=16,20,24(,28,32) recovers v₅ to ~5%; naive 1/n Richardson FAILS (pre-asymptotic n≤14 is sign-oscillating);
  **term-wise extrapolation FAILS** (assembles −10.6 vs −1.19 — per-term tails not asymptotic by n=20). v₆ was
  therefore read from its 16/20/24 plateau, which is v₄-like (fast), not v₅-like (slow).

## Engine work (reusable)

- **Gaussian quad fast path** (`chaos_transfer._moment_quad`, routed in `_moment_bulk_gen`): pure-bulk moments
  where every atom has order ≤ 2 are Gaussian expectations of quadratic+linear forms → connected-diagram
  expansion (chain-cycles = traces, leg-paths = vector sandwiches; weights: 2-cycle 2·tr, k≥3-cycle 2^k per
  distinct arrangement, path 2^m) assembled by counted matching recursion. **Machine-precision match with the
  DP on all classes; ~2×10⁴–3×10⁴× on the v₆ bottleneck class** (3 chains + 8 legs: 12.5 s → 1 ms at n=14).
  This is what made n=20/24 affordable.
- **Rank-2 G rewrite: RULED OUT (first-class negative result).** G = U⊗Ψ − Ψ⊗U is rank 2, and the separable
  DP (state = slot-count + pending factor, n-independent) is implemented in `_moment_sep` — but it is
  numerically DEAD in double precision: u₀ spans ~16 orders of magnitude across the interval, the branch
  products U·Ψ exceed |G| by ~10⁴, and the per-bond cancellation compounds (measured rel. error 10⁵ even for
  single-bond chains, 10²⁰ for the slow class). No constant SL(2) basis change fixes both ends of a Wronskian
  pair. Kept in `chaos_transfer.py` behind `SEP_ENABLE=False` with a post-mortem comment. Do not re-attempt
  without ≥ quad precision (which erases the speed win).
- **Order-15 ladder built** (v₇-ready): `_gen15.py` (term-by-term substitution+expand — the OOM-safe pattern;
  note the handoff's "git history" pointer is stale, there is no git repo), `_yexprs_15.txt` (Y₁₅ = 91,159
  U-terms), `_ybase_14/15.pkl` (4,574 / 7,470 base monomials). Driver `_yfile` extended; `chaos_diagram`
  symbol ranges extended to order 15 (backward-compatible).
- New tools: `_borel_analysis.py` (Borel–Padé + Darboux/Dingle + two-singularity fit + θ-scan),
  `_extrap_terms.py` (term-wise extrapolation — kept as the documented failed method), `_v6_sweep.py`,
  `_vk_smallgrid*.py`, `_vk_ladder*.py`, `_v5_terms.py`, `_test_quad.py`, `_test_sep.py`.

## Cost calibration (for planning v₇)

Single core, quad-hybrid engine: v₆ per grid ≈ 2.1 h (n=10), 2.4 h (n=12), 4.7 h (n=14), 8 h (n=16),
11 h (n=18), 14.8 h (n=20), 29 h (n=24). v₇(n=10) = 28.8 h (Var(Y₈) alone 22.5 h) — **~14× v₆ per grid**.
A v₇ ladder to n=24 at this cost is ~2–3 weeks single-core: NOT viable as-is. The residual bottleneck is
≥3-chain moments containing an order-≥3 chain (not quad-eligible, generic DP state ~n^#chains). Recommended
route: numba/C rewrite of the DP inner loop (pure numeric dict-loop; 50–100× plausible) → v₇ ladder in days.
v₇(n=12) left running as a cross-check point for that rewrite.

## What this changes upstream

- `PROGRAM2_CONSOLIDATION.md` §3quater/3quinque: the "54–63°, NOT 45°" phase claim is RETRACTED (6-coefficient
  artifact). New: **θ = 48° ± 3°, 45° allowed, ≥53° excluded**; radius |ζ| ≈ 1.8–2.6 (single-pair fit; real
  instanton at ~1.9 unresolved from the pair — three-singularity competition confirmed, needs v₇).
- The Thread-4 λ₀→Borel bridge is REOPENED as plausible (phase consistent with inheritance).
- Deliverable (ii) of the handoff (real-instanton separation) remains open pending v₇ — scoped as follow-up
  with the DP rewrite as its first step.

## ADDENDUM (2026-07-15, same day) — real-pole subtraction executed: 45° DISFAVORED, θ ≈ 50°

Milestone-2's prescribed real-instanton subtraction (never previously executed) was run as a z_r-profiled
two-singularity fit: pin the FW real pole across z_r ∈ {1.4…2.6} × a_r ∈ {−½,0,+½}, fit pair + residue,
scan θ. Dual independent implementations, cell-by-cell agreement, reconciler re-verified. Findings:

1. **Subtraction does NOT move the phase.** In the overfit-protected α=0 variant (2 dof): θ = 50° (band
   49–51) in all 24 (z_r,a_r) cells and all 6 v₆/v₃ robustness variants — identical to the unsubtracted
   baseline. The stiff fit's phase is immovable at ~50°.
2. **The fitted real-pole residue is NEGATIVE everywhere trustworthy** (C_r = −0.002…−0.034). The
   45°-masquerade mechanism (bias experiment) requires POSITIVE C_r; positive C_r appears only in α-free
   fits with α railed at its +1.4 ceiling (20/24 cells; 6 params on 7 data = 1 dof) — a pure overfit
   signature. 45° enters the band only through those railed fits.
3. **Verdict: exact-45° (radial dressing / λ₀ phase inheritance) is DISFAVORED at 7 coefficients** — not
   excluded (caveats: α=0 stiffness assumption; pair+one-real-pole model; z_r profiled not theory-pinned),
   but the earlier "45° allowed" reading is demoted: it survives only via overfitting artifacts. The
   robust statement is now **θ ≈ 50° ± 2°, a genuine small (~5°) rotation off λ₀'s 45° as the leading
   hypothesis**, with ≥53° still hard-excluded.
4. **Consequence for theory:** the fluctuation-rotation question (Milestone-2 Mechanism C) is partially
   re-activated but must deliver BOTH θ≈50° (needs untheorized c≈0.5, natural c=1 excluded) AND
   |ζ|≈1.8–2.0 (Mechanism C predicts ≈1.24) — no current mechanism gives both. The v₇ gate sharpens:
   pre-registered v₇ = +1.3 (45°) vs +5.5 (48.6°) vs +8.2 (50°); a converged v₇ near +6–8 confirms
   genuine rotation, near +1 reopens 45°.

## ADDENDUM (2026-07-19) — v₇ engine campaign: 5 iterations, NOT reached; scoped as future work

v₇ needs Var(Y₈) (8th-order Wiener chaos). Its moment set (6450 distinct in ⟨Y₈Y₈⟩ at n=10) contains
classes that broke the engine in FIVE successive ways; each fix exposed the next scaling wall. Full record so
the next attempt starts from the frontier, not the beginning:

1. **numba DP** (`chaos_dp_numba.py`, gated 2026-07-15): rewrote the pure-python dict-DP inner loop as a
   numba mixed-radix kernel. Battery 4.3e-15, v₆(n=10) anchor in 390 s vs 2.1 h (~19×). CORRECT + in
   production. But cost is still ~n^(#chains) → walls for many-chain moments.
   *(Note 2026-07-18: a `chaos_dp_numpy.py` vectorized-DP fallback was added — same DP, numpy-vectorized
   per-timestep, for environments without numba. Battery-clean; same n^(#chains) scaling.)*
2. **Held–Karp quad** (`_moment_quad` subset-DP): the all-order≤2 path enumerated chain permutations
   (factorial). Replaced with a subset-DP → 64× on the 8-chain class. CORRECT + in production.
3. **Routing fix** (2026-07-17): `moment_hybrid` sent pure-bulk order≤2 many-chain moments (e.g.
   (2,2,2,2)+8 legs) to the DENSE engine (35 s) instead of quad (2 ms) — a 17,000× miss. Now all pure-bulk
   routes through `moment_transfer`. Battery-clean; CORRECT + in production.
4. **Hybrid soft-folding** (`chaos_hybrid.py`): fold order-2 chains + legs into the quad connected-diagram
   machinery, DP only over the order-≥3 hard core. **Math PROVEN correct** (`_tmp_hybrid_proof.py`: 200+
   mixed moments, max rel 3.95e-14 at n=6,8). Soft sector genuinely grid-independent (1 hard chain: flat
   0.000 s at n=10…24). NOT wired to production.
5. **The wall that stopped it:** the hybrid's hard-core contraction (`_contract_hard`) memoizes the grid-TIME
   of every half-open "bridge" (an order-2 soft chain linking two hard ports). For the real v₇ killer moments
   with an order-5 hard chain (e.g. (2,2,2,2,3,5), 8 hard ports, ~4 simultaneous bridges) this reintroduces
   n^(#open-bridges) → measured **12.7 s (n=10) → 5203 s (n=16)**, ~n^13, WORSE than the numba DP it replaces.

**Exact fix for a future attempt (well-scoped):** rework `_contract_hard` as a grid-time-ORDERED sweep that
carries open bridges as forward vector-messages (accumulate Σ[t,·] as a running message) instead of
memoizing their scalar open-times — targets ~n²–n³ (n^(#hard chains)). The soft-folding derivation is
correct and reusable; only the hard-core port/bridge contraction needs the transfer rewrite. No numerical
obstruction (bridges are full n×n S-matrices, not the rank-2 factorization that was ruled out for
cancellation). Then: gate on the real 5/6/7-chain Var(Y₈) battery + v₇(n=10)=−178.6084 anchor, run the
ladder.

**v₇(n=10) anchor status:** −178.60840 confirmed by the ORIGINAL pure-python engine (its (8,8)=+367.82382 and
(1,15)=−119.56156 terms reproduced across two independent runs); no faster engine has yet reproduced the full
value end-to-end (all fast engines wall on the same many-chain Var(Y₈) moments).

## DECISION (2026-07-19): FINALIZED ON v₆

v₆ = −1.90 ± 0.15 (θ ≈ 50°, exact-45° disfavored) is the deliverable. v₇ — deliverable (ii), the
real-instanton separation and the sharp 45-vs-50 discrimination — is deferred as scoped future work (the
`_contract_hard` transfer rewrite above is its single first step). v₇ was always the STRETCH goal; v₆ already
establishes the headline (exact-45° disfavored, θ≈50° leading), so the science is complete without it.

## Next steps (priority order) — for a future session

1. **`_contract_hard` transfer rewrite** — ✅ DONE 2026-07-20, but see the addendum below: it is NOT the v₇
   gate. Do not re-attempt as "the gate."
2. **v₇ grid ladder + plateau read** → discriminates θ=45° vs 50° sharply (pre-registered +1.3 vs +8.2) and
   feeds the 8-parameter two-singularity fit (deliverable (ii)). **BLOCKED by assembly-scale cost, not by any
   contractor — see addendum. Not feasible on the current single machine (~days/grid at n=20–24).**
3. Optional v₆ hardening: n=28 (~cheap now with numba, was ~2.5 days) to shrink ±0.15 → ±0.05.
4. Meijer-G approximant (deferred — with 7 coefficients Darboux+Padé already agree).

## ADDENDUM (2026-07-20) — `_contract_hard` rewrite DONE + CORRECT, but it is NOT the v₇ gate (premise retracted)

The scoped `_contract_hard` transfer rewrite (2026-07-19 addendum, item "Exact fix") was implemented and
validated. **But benchmarking it against the ACTUAL production engine overturns the "it is the one gate"
premise.** Honest record so no future session re-runs this:

**What was built (banked, reusable):** `chaos_hybrid._contract_hard_sweep` (gated by `HYBRID_SWEEP=True`).
Same recursion as the old `_contract_hard`, but `solve` now RETURNS A TENSOR over the open-bridges' open-times
and memoizes only on the DISCRETE state (ci, si, lastt, remain, open-bridge indices) — open-times became numpy
axes summed in C, killing the n^(#open bridges) memo blowup.
- CORRECT: machine-precision (rel ≤ 1.5e-14) vs the reference DP on the full `_tmp_hybrid_proof` battery
  (200+ mixed classes) AND on the real killer classes ((2,2,2,2,3,5) etc.).
- FAST: on (2,2,2,2,3,5), old hybrid `_contract_hard` = the wall (5203 s @ n=16, per the 07-19 addendum);
  the rewrite does it in **2.85 s** (~1800×), scaling ~n^3.2 instead of ~n^6. The "no numerical obstruction"
  claim held.

**Why it is NOT the gate (the decisive finding):** the prod path is the **numba DP** (`chaos_dp_numba`, in
`_moment_bulk_gen`), added AFTER Handoff A was written. Measured, prod numba DP vs the (fixed) hybrid, on the
shapes that appear in Var(Y₈):
| shape | chains/hard | numba | hybrid | winner |
|---|---|---|---|---|
| (2,2,2,2,3,5) | 6/2 @ n=20 | 1.3 s | 6.8 s | **numba** |
| (2,2,2,2,2,2,3,3) | 8/2 @ n=20 | 3.5 s | 15 s | **numba** |
| (2,2,2,2,2,2,2,3,3) | 9/2 @ n=16 | 33.9 s | 36.0 s | ~tie |

The hybrid (n^(#hard chains)) only overtakes numba at **≥9–10 total chains** — but those need 18–20 ports and
⟨Y₈Y₈⟩ tops out at combined order 16 (≤8 chains), so **those shapes never appear.** The numba DP also prunes
aggressively (not naive n^(#chains)). Net: **swapping in the hybrid does not speed up v₇.**

**The REAL v₇ bottleneck = assembly scale, not any contractor.** Var(Y₈)+cross = **27,033 distinct moments**
(152k pre-memo calls; vk_terms(7) = {(8,8),(1,15)…(7,9)}; |Y₁₅|=7470 base monomials). Multi-core driver built
(`coupled-atlas/_v7_parallel.py`, Phase-1 enumerate → Phase-2 parallel compute → Phase-3 assemble; runs clean
on 10 cores, n=8 end-to-end = 508 s). But the aggregate wall scales ~n^5.3 (n=8: 8.5 min → n=10: ~30 min on
9 cores), extrapolating to n=16 ~6 h, n=20 ~20 h, n=24 ~2 days **per grid** — and n≤14 is pre-asymptotic, so a
usable v₇ needs the n=16–24 plateau. **Not feasible on one machine.** (The old 28.8 h single-core v₇(n=10) is
now ~30 min on 9 cores — a real speedup, but nowhere near enough for the ladder.)

**DECISION UNCHANGED: finalized on v₆.** v₇ would need either a cluster (the assembly is embarrassingly
parallel — `_v7_parallel.py` is ready) or a fundamentally cheaper assembly (dedup sub-moments across workers /
a smarter moment-generating-function contraction), not a new contractor. Deliverable (ii) stays open. The
`_contract_hard` rewrite is a correct engine improvement kept for the record; it is not the lever.
Investigation scripts: `_v7_baseline.py`, `_v7_timing.py`, `_v7_route_probe{,2}.py`, `_v7_crossover{,8}.py`,
`_v7_parallel.py`.
