# MMO chapter — Phase 3.1 (exponent tightening): handoff to Opus

This is a research handoff. Phase 3 validated the MMO noise-dissolution
**mechanism** but left the **exponent** not precision-tested. Your job is to
pin it — and, critically, to first decide whether the apparent exponent gap
is a bookkeeping artefact or genuine missing physics, before anyone invests
in the hard Path-A K2 inner solution. Two sub-tasks, **cheapest first**:
(α) recompute the 3B prediction with the *exact* μ at each plateau instead of
the idealised μ ∝ 1/q — nearly free, reuses existing code, and is the decisive
test; (β) a dt-/H-converged SAO-counting-under-noise study that resolves q = 2
and q = 5–6 and measures β with real error bars.

## 1. Brief context

The project is a regime-by-regime theory of how degenerate noise (σ on the
fast variable v only) affects slow–fast neuronal dynamics. Five chapters:
excitable, canard (σ_* = C_q√ε·λ^{1/2}), tonic, resonator, and **MMO** —
mode-locked L^a S^b patterns in 3D FitzHugh–Rinzel, σ_pq the noise that
dissolves each Farey plateau p/q.

MMO state: deterministic side closed (Phases 1.5–2: folded node, realised
α ≈ 1.55, κ ≈ 2π², funnel-filling f(c)). **Phase 3 validated the noise
mechanism** (`MMO_NOISE.md`): degenerate noise dissolves the staircase
high-q-first (dt-converged), with σ_* = C_q√ε·F(μ), F(μ) ~ μ^{3/2}, and
**C_q ≈ 10 — the full-FHN canard chapter's prefactor** (the cross-chapter
unification the chapter was built to deliver). The one open item is the
**exponent**: measured β ≈ 1.15 (a 2-point q = 3, 4 slope) vs the 3B-derived
β ≈ 1.5, γ ≈ 1. The chapter currently sits at "bin 1–2: mechanism + order +
scale + prefactor validated; exponent consistent but not precision-tested."

## 2. Read first (in order)

1. `regime-tests/MMO_NOISE.md` — the Phase-3 result you are tightening
   (σ_pq table, the 3B Path-B derivation, the honest-limits §3).
2. `regime-tests/MMO_K2.md` — the deterministic source of the three objects
   Task α needs: the funnel-filling **f(c)**, the smallest SAO amplitude
   **a_min(c)** (per plateau), and the folded-node μ.
3. `regime-tests/MMO_ALPHA_DERIVATION.md` — μ(c) = 2bδ(c+1)/(1+δ)², the
   ceiling s_max = (1−μ)/(2μ), and the measured f(c) (0.21 → 0.44 across the
   band).
4. `regime-tests/CANARD_BLOWUP.md` §4 — the σ_* accumulated-variance argument
   that 3B generalises (so you can extend it cleanly if Task α points to real
   physics rather than bookkeeping).
5. `regime-tests/mmo_noise_staircase.py` — the Phase-3A script you will
   extend (EM integrator, the w-peak counter, `sigma_pq`). Note `mu(c)` lives
   in `mmo_k2_return_map.py` (and `mmo_fhr_foldednode.py`) — reuse it.

## 3. What's already established (rely on, don't re-derive)

- **Working FHR, δ = 0.2** (the operating point where α, f(c), κ, and the
  Phase-3 σ_pq were all pinned — keep it, for self-consistency; the δ-guidance
  in `MMO_PHASE3_PROMPT.md` §3 governs any switch).
  `dv = (v−v³/3−w+y+I)dt + σ dW`, `dw = ε(v+a−bw)dt`, `dy = εδ(c−v)dt`,
  `(a,b,ε,δ,I) = (0.7,0.8,0.08,0.2,0.30)`, c the bifurcation control.
- **The law:** σ_pq ~ σ_*·q^{−α/γ}, σ_* = C_q√ε·F(μ), F(μ) ~ μ^{3/2}, γ ≈ 1
  ⇒ β = α/γ ≈ 1.5 with α = 1.55. Prefactor C_q ≈ 10 (= canard).
- **The counter:** spikes by v-hysteresis (cross +0.5, reset −0.5); SAO loops
  as prominence-H peaks of **w** (smooth — noise enters only v, w
  low-pass-filters it). v-based SAO counting is broken under noise (≈10³
  spurious maxima). ρ = nL / (w-peaks).
- **Measured Phase-3 σ_pq** (H = 0.008, N = 70, dt = 0.02): q = 3 → 0.022,
  q = 4 → 0.016; q = 2 unresolved (> 0.045, σ-grid too small for the wide
  plateau); q ≥ 5 below the SAO-countability floor. β ≈ 1.15 from the two
  resolved points. σ_pq is H-sensitive at the factor-2 level (q = 3: 0.022 at
  H = 0.008 vs 0.010 at H = 0.03).

## 4. The decisive reframing — exact μ(q), not μ ∝ 1/q (Task α)

3B's β = 1.5 used the *leading* idealisation μ ∝ 1/q. The exact relation is
not 1/q. At the ρ = 1/q plateau the SAO count is s = q − 1, and
s = f(c)·s_max = f(c)·(1−μ)/(2μ), so the μ that sets σ_pq is **fixed by the
plateau's actual c**:

```
μ(c_pq) = 2 b δ (c_pq + 1) / (1+δ)²          [closed form, = mu(c)]
σ_pq^pred(q) = C_q · √ε · μ(c_pq)^{3/2}      [exact 3B, no μ∝1/q idealisation]
```

Because the measured **f(c) rises** 0.21 → 0.44 across the band, the
self-consistent μ(q) = f(q)/(2(q−1)) falls **slower** than 1/q, which *lowers*
the predicted local slope below 1.5 — toward the measured 1.15. **But do not
assume it closes the gap:** there is a competing correction — at small q the
(q−1)-vs-q distinction *raises* the local slope (d ln(q−1)/d ln q > 1 for
q = 3, 4). The net is an empirical question. Compute the exact-μ(q) prediction
curve over the measured q-range and read off its slope; that number — not 1.5
— is the honest 3B target the measurement should be compared against.

**Task α steps (≈1 session, no new simulation):**
1. For each plateau, take its **measured centre c_pq** (from the deterministic
   staircase, `mmo_fhr_staircase.py`, or the Phase-3 table) and compute μ(c_pq)
   via `mu(c)`.
2. Form σ_pq^pred(q) = C_q√ε·μ(c_pq)^{3/2} (C_q ≈ 10 fixed from Phase 3).
3. Fit its log-log slope vs q over (i) the measured q = 3, 4 and (ii) the
   extended q = 2…6. Report β_pred^exact for both. Compare to β = 1.5
   (idealised) and the measured β.
4. Verdict: how much of the 1.15-vs-1.5 gap is the μ∝1/q idealisation? If
   β_pred^exact ≈ measured β, the gap was bookkeeping and **Path A is not
   needed for the exponent**. If a residual remains, Task β + Path A are
   justified.

## 5. Task β — the dt-/H-converged counting study

The exponent measurement is gated by SAO-counting-under-noise. Tighten it.

1. **dt convergence as a logged sweep, not a spot-check.** Run σ_pq at
   dt ∈ {0.04, 0.02, 0.01, 0.005} for q = 3 (and one high-q plateau) and show
   σ_pq has plateaued by dt = 0.01. Lock dt = 0.01 for production.
2. **Per-plateau, geometry-grounded H.** The single global H = 0.008 is the
   limiting systematic (factor-2 sensitivity; high-q loops fall below it). Set
   H **per plateau** as a fraction of the deterministic smallest-loop amplitude
   **a_min(c_pq)** from K2 — large enough to reject v-noise leakage into w,
   small enough to catch the genuine deep-funnel loops. Report σ_pq vs H around
   that choice to show it is now stable, not factor-2.
3. **Resolve q = 2.** Extend the σ-grid above 0.045 (its σ_pq > 0.045) and
   refine near the knee so q = 2 yields a finite σ_pq, not a lower bound.
4. **Resolve q = 5, 6.** Locate their plateau centres from the deterministic
   staircase (q grows fast as c → −1: q = 5 ≈ c ≈ −0.86…−0.87). Use longer T
   and/or more seeds (N ≥ 200) for statistics on the small loops; the
   geometry-grounded H (step 2) is what makes them countable.
5. **Measure γ directly.** At a fixed in-plateau c, fit the rotation-number
   spread std(ρ)(σ) ~ σ^γ in the small-σ regime (below dissolution). This is
   an independent check on γ ≈ 1, separate from γ = α/β.
6. **Re-fit β** over q = 2…6 with genuine error bars (bootstrap over seeds and
   over the H-choice). Report β_meas ± its honest systematic spread.

## 6. Falsifiability / decision tree

- **Outcome 1 — the gap was bookkeeping.** β_pred^exact (Task α) ≈ β_meas
  (Task β), both ≈ 1.1–1.3. The μ∝1/q idealisation was the whole story; the
  exponent is reconciled at the mechanism level and **Path A is unnecessary
  for it**. Flip the README MMO row to a flat "derived + validated."
- **Outcome 2 — the 2-point slope was the artefact.** With q = 2, 5, 6 and
  converged counting, β_meas moves up to ≈ 1.4–1.5 and matches β_pred^exact.
  Again reconciled; the original 1.15 was the under-resolved measurement.
- **Outcome 3 — a real residual.** β_pred^exact and a now-trustworthy β_meas
  disagree beyond error (e.g. exact-μ predicts ≈ 1.4, converged measurement is
  solidly ≈ 1.15 ± 0.1). That is genuine physics the bound + measured-κ Path B
  misses → **Path A (Wechselberger K2 inner solution) is justified**; frame it
  precisely (what the residual is, which observable it lives in) as the next
  handoff. This is the same Path-A frontier as the deterministic side
  (closed-form a_min(c), first-principles κ = 2π²).

State which outcome you land in explicitly.

## 7. Hazards

1. **SAO-countability floor (the limiting issue).** High-q deep-funnel loops
   are tiny; the geometry-grounded per-plateau H (§5.2) is the fix, but q ≥ 7
   may remain genuinely uncountable under noise — if so, say so and cap the fit
   at the resolvable q.
2. **Don't manufacture a clean β.** If after converged counting β_meas still
   carries a ±0.2 systematic, report the band — a tightened-but-honest
   β = 1.3 ± 0.2 is a real result; a forced two-figure number is not.
3. **q = 2 wide plateau.** Its dissolution knee is at larger σ where the
   ⟨ρ⟩-drift is gentle; make sure the locked-fraction crossing is real and not
   grid-edge extrapolation.
4. **Cost / δ fallback.** If converged counting at δ = 0.2 (×N seeds × σ-grid
   × q = 2…6) is too slow, do **not** mix regimes — re-pin the whole chain
   (α, f(c)) at δ = 1.0 and run everything there, per the δ-guidance in
   `MMO_PHASE3_PROMPT.md` §3. The folded-node structure survives the move
   (timescale check).

## 8. Deliverables

1. Revise **`MMO_NOISE.md`** §3–§5 with the converged σ_pq(q), β ± error, the
   direct γ, and a new short subsection "Exponent tightening (Phase 3.1)"
   documenting Task α (exact-μ(q) vs idealised) and the Outcome 1/2/3 verdict.
2. Extend **`mmo_noise_staircase.py`**: the dt-sweep, per-plateau a_min-scaled
   H, q = 2 extended σ-grid, q = 5–6, the exact-μ(q) 3B curve, and the direct-γ
   measurement. **While there, fix the stale stdout string** (the summary line
   asserting `β∈[~1.5(3B), ~2–4(meas)]` contradicts the computed/ documented
   β ≈ 1.15 printed just above — make it read the actual fitted β).
3. Update **`figures/mmo_noise_staircase.png`**: σ_pq(q) over q = 2…6 with both
   the idealised q^{−3/2} line and the **exact-μ(q) 3B curve** overlaid, plus a
   panel for the direct γ fit (std(ρ) vs σ).
4. Add **banner line (8)** to `MMO_FHR_PLAN.md`: "(8) Phase 3.1 exponent
   tightening: [outcome] — β_meas = …, β_pred^exact = …, γ = …; [bookkeeping /
   real residual → Path A]."
5. If Outcome 1 or 2: update the README **§1, §7, §9** MMO rows from
   "exponent ballpark / not precision-tested" to the reconciled statement.
   If Outcome 3: keep the honest "exponent open" and add the precise Path-A
   framing to §9.2.

## 9. Reading priority if time-constrained

1. `MMO_NOISE.md` (what you're tightening).
2. `MMO_K2.md` — f(c), a_min(c), and `mu(c)` (Task α's inputs).
3. `CANARD_BLOWUP.md` §4 (the σ_* argument, for Outcome 3 / Path A framing).

Do **Task α first** — it is nearly free and may settle the exponent without
any new simulation. Only run the full Task β counting study if Task α leaves a
residual worth resolving.
