# HANDOFF A — compute v₆ (and v₇) exactly with the transfer engine, then pin the Borel phase

> **STATUS 2026-07-15: EXECUTED (main deliverable done)** — see `PROGRAM2_V6_GRIND_NOTES.md`.
> v₆ = **−1.90 ± 0.15** grid-certified (n=10…24 ladder, dual-engine validated, boundary budget 0.1%).
> **Sign-flip prediction (+2) REFUTED; θ re-pinned to 48° ± 3° — the "54–63°, NOT 45°" claim is retracted,
> λ₀'s 45° is allowed again.** Quad fast path added to the engine (~10⁴×); rank-2 G RULED OUT (catastrophic
> cancellation — see note). v₇: order-15 ladder + first datum done (n=10: −178.6, pre-asymptotic); full v₇
> ladder needs the DP-loop numba/C rewrite (~14× v₆ cost per grid as-is) — scoped as follow-up.

_Self-contained handoff for a fresh Claude Code session. Route (a) of Program 2's endgame: the transfer-operator
engine is BUILT and VALIDATED (2026-07-10); what remains is the compute grind + the Borel-plane analysis it feeds.
Read `PROGRAM2_CHAOS_ENGINE.md` §"Transfer-operator rewrite" and `PROGRAM2_CONSOLIDATION.md` §0 for background._

## The goal (one sentence)
Compute the weak-noise variance coefficient **v₆** (and ideally **v₇**) of 𝒲 exactly and grid-converged, then
redo the Borel/Darboux analysis on the extended ladder to pin the singularity phase θ (currently ~54–63°, vs the
deterministic λ₀'s 45°) and separate the complex-conjugate Borel pair from the competing real instanton.

## Where things stand (what is DONE, so don't rebuild it)
- **The transfer engine works and is validated.** `coupled-atlas/chaos_transfer.py` computes any Wiener-chaos
  moment `<∏ atoms>` via a grid-time sweep that NEVER materializes the n^k symmetric tensor (the wall that blocked
  v₆). Validated to **machine precision at orders 6 and 7** against the dense engine (`chaos_diagram.py`) built at
  MAXORD=6/7 on small grids. It reproduces the known ladder: **v₃=−0.030, v₄=−0.449/−0.460** (match), and its
  boundary convention is self-consistent (v₃ agrees with dense to 3 digits at every grid — the δ(0) boundary
  divergences cancel correctly in the assembly).
- **The order-13 functionals exist.** `_yexprs_13.txt` (Y₁…Y₁₃ in U-symbols) and `_ybase_1.pkl … _ybase_13.pkl`
  (their base-atom decompositions) are built and cached. Y₁₃ = 2764 base monomials, Y₁₂ = 1641.
- **Known ladder (targets to reproduce):** v₀=0.134, v₁=0.111, v₂≈0.10, v₃≈−0.03, v₄≈−0.45, v₅≈−1.1.
  **Prediction to test: v₆ ≈ +2** (a sign flip back to positive — the Darboux fit predicts this; confirming it is
  the scientific payoff).

## The formula
v_k = Var(Y_{k+1}) + 2·Σ_{j=1..k} Cov(Y_j, Y_{2k+2−j}). For **v₆** the (a,b) pairs are:
`(7,7)` [the Var], and `(1,13),(2,12),(3,11),(4,10),(5,9),(6,8)` [each ×2]. Y-indices needed: 1…13 (all cached).

## How to run it
```
cd coupled-atlas
python3 _v6_driver.py v6 12          # computes v6 at grid n=12, prints per-term progress
```
- `_v6_driver.py`: `compute_vk(k, n, mom, verbose)` assembles v_k; `vk_terms(k)` gives the pairs; `_yfile(order)`
  picks the Y-expression file. Default engine is `chaos_transfer.moment_hybrid` and `CT.BND_ENGINE='transfer'`
  (the validated consistent boundary convention — KEEP THIS).
- `_v6_breakdown.py`: same computation but separately reports the contribution of boundary+order>5 moments (the
  small "convention-caveat" set — 107 of 3155 distinct moments), for an honest error budget.
- The engine auto-memoizes (`_HMEMO`, `_RMEMO`, `CD._MCACHE`); `CT.clear_all()` is called per setup.

## The plan
1. **Run v₆ across grids** n = 10, 12, 14 (optionally 16). Each grid is ~tens of minutes (see performance below).
   Launch in the background; the driver prints each (a,b) term as it finishes so you can watch progress.
2. **Richardson-extrapolate** v₆(n) → n→∞. v_n carries an O(1/n) grid error from the boundary evaluations
   (`beval`/`reduce_head`), so extrapolate in 1/n. Cross-check the trend against v₃,v₄,v₅ (whose n-dependence is
   documented: e.g. v₃ = −0.060/−0.030/−0.021 at n=20/30/40).
3. **Sanity gates:** the dominant term is Var(Y₇), itself dominated by `<U₇₀ U₇₀>` (the 2-chain self-contraction,
   fast). Confirm `2·Cov(Y₃,Y₉)`-type terms are large and negative (they drove v₄,v₅). If v₆ comes out large and
   **positive** (~+1 to +3), that is the predicted sign flip — the headline result.
4. **(Stretch) v₇:** needs Y₁₅ → extend `_gen13.py`→ order 15 (`sp.solve` version; ~10–15 min, watch memory — it
   OOM'd at order 12 once; the term-by-term-substitution trick in the git history of `_gen13.py` fixes that) and
   build `_ybase_14.pkl,_ybase_15.pkl`. Then `vk_terms(7)` needs Y up to 15.
5. **Borel/Darboux on the extended ladder** (the actual scientific deliverable): with v₀…v₆ (or v₇), redo the
   large-order analysis. Recommended tools (from `PROGRAM2_CONSOLIDATION.md` §3ter): **Meijer-G / hypergeometric
   approximants** (converge at 3–5 orders, native branch cut — beat Borel–Padé for few terms; refs Mera–Pedersen–
   Nikolic arXiv:1802.06034, Crew–Trinh arXiv:2208.07290) + Darboux/Dingle late-terms inversion, Borel–Padé as
   cross-check. Deliverables: (i) pin θ (is it really ~54°, i.e. NOT λ₀'s 45°?); (ii) separate the complex pair
   from the real s⁵/10 instanton (needs v₆,v₇ — 8 params > 6 data before that); (iii) radius |ζ|.

## Performance — the one thing that needs work
v₆ has ~3155 distinct moments (memoized). Most are fast. The bottleneck is a class of **order-14 "few-chains +
many-legs" moments** (e.g. 3 order-2 chains + 8 order-1 legs, from Var(Y₇)): ~3.6 s each, dozens of them, and
NEITHER dense (14-atom multigraph blowup) NOR transfer (n^{#chains} state) is fast for them. This makes v₆
~30–60 min/grid. Two concrete optimizations if you want it faster (both are clean, bounded work):
- **A fast all-order-≤2 Gaussian-moment path.** Moments where every atom is order ≤2 (I₂ quadratic + I₁ linear
  forms) are a Gaussian expectation computable as a hafnian/Pfaffian with 2×2 blocks — O(#atoms³) instead of the
  DP. This kills the dominant slow class. Add it as a routing branch in `moment_hybrid`.
- **Exploit rank-2 G.** The bond `G[i,j]=U[i]Ψ[j]−Ψ[i]U[j]` is rank-2, so a chain's last-time axis collapses to a
  2-vector message — this can vectorize/accelerate the multi-chain sweep. Bigger job.
Routing lives in `chaos_transfer.moment_hybrid`; the current rule (documented inline) is: all-degree-1 → hafnian;
boundary → transfer (consistent convention) if order>5 else dense; pure bulk → transfer for order>5 or (≤2 chains
& ≥6 legs), else dense.

## Gotchas (learned the hard way — do not re-discover)
- **Boundary moments are individually grid-DIVERGENT** (δ(0) self-contractions); they only become finite after
  cancellation in the assembled v_k. This REQUIRES one consistent boundary convention across all of them — keep
  `BND_ENGINE='transfer'`. Mixing dense (order≤5) and transfer (order>5) boundary moments is acceptable ONLY
  because both converge to the same continuum (verified: the two conventions give v₃ agreeing to 3 digits at every
  grid); do not "fix" one to match the other at finite n.
- `reduce_head` (transfer's boundary reduction) differs from dense's `beval` at finite n for j≥2 boundary
  derivatives (a discretization branch), but both converge — this is fine and expected, not a bug.
- Don't try to reach v₆ by building dense symmetric tensors at order 6/7 on the real grids — that is the memory
  wall (order-7 at n=14 is ~8 GB). The whole point of the transfer engine is to avoid it.
- Two cheap shortcuts to the phase were RULED OUT this session (see `PROGRAM2_STOCHASTIC_WKB_DEEPRESEARCH.md`):
  the FP-generator spectrum (wrong object) and PDE small-η re-extraction (resolution wall). Don't re-run them.

## Success criteria
- v₆ grid-converged to ±~0.2, sign and rough magnitude robust across grids and across the boundary-error budget.
- The Borel analysis on v₀…v₆ gives a stable θ (settling the 54°-vs-45° question) and, with v₇, disentangles the
  real instanton. That closes `PROGRAM2_CONSOLIDATION.md`'s "one remaining large build (deferred)".

## Key files
`coupled-atlas/`: `chaos_transfer.py` (engine), `chaos_diagram.py` (dense ref + Ybase/setup), `chaos_engine.py`
(atom kernels), `_v6_driver.py`, `_v6_breakdown.py`, `_gen13.py` (Y-expression generator), `_yexprs_13.txt`,
`_ybase_*.pkl`. Notes: `PROGRAM2_CHAOS_ENGINE.md`, `PROGRAM2_CONSOLIDATION.md`, `PROGRAM2_ROUTE2B_NOTES.md`.
