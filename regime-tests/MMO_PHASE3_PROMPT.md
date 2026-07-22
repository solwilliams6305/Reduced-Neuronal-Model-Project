# MMO chapter — Phase 3 (noise dissolution): handoff to Opus

This is a research handoff. You are picking up the FHN-MMO chapter at the
point where the **deterministic side has closed** (Phases 1→2) and the
remaining work is the chapter's distinctive contribution: the **noise**
that dissolves the Farey staircase. The task has a numerical half (measure
the q-dependent dissolution threshold σ_pq) and an analytical half (derive
the noise scale σ_*(μ,ε,δ) and the rotation-spread exponent γ, then
predict σ_pq). Methodology is **measure-first, then derive** — exactly how
the chapter handled the deterministic exponent α (measured 1.55, then K2
reproduced it).

## 1. Brief context

The broader project is a regime-by-regime theory of how **degenerate
noise** (amplitude σ on the fast variable v only) affects slow–fast
neuronal dynamics. Five chapters:

- excitable — σ_crit = √(ε/C);
- canard — σ_* = C_q · √ε · λ^{1/2} (validated normal-form + full FHN + Van der Pol);
- tonic — CV = σ · A_mid, matched-asymptotic adjoint Floquet;
- resonator — subcritical-Hopf B ∼ δ² corner law + phase-gated instanton;
- **MMO** (current) — mode-locked L^a S^b patterns in 3D FitzHugh–Rinzel,
  with σ_pq the noise that dissolves each Farey plateau p/q.

The deterministic MMO mechanism is folded-node canard theory
(Wechselberger 2005) — established mathematics that Phases 1–2 applied to
FHR. **The noise extension is the chapter's new contribution.** Closing
Phase 3 brings MMO to the same level as canard/tonic, and delivers the
cross-chapter unification σ_pq ∼ σ_* · q^{−α/γ} linking back to the canard
chapter's σ_*.

## 2. Read first (in order, before computing anything)

1. `README.md` §1, §7 — master σ_crit table (the MMO row now reflects the
   3D FHR state), and §4 (the canard chapter you will generalise).
2. `regime-tests/MMO_FHR_PLAN.md` — full plan. Read the **banner**
   (six-correction trail, ending at the Phase-2 K2 result) and §3 (noise
   framework) + §4 Phase 3 (Tasks G′/H′/I′) + §6 (3D hazards).
3. `regime-tests/MMO_K2.md` — the deterministic closure you are adding
   noise to: the K2 return map, the folded-node rotation constant κ ≈ 2π²,
   the funnel-filling f(c), and α ≈ 1.45–1.55.
4. `regime-tests/MMO_ALPHA_DERIVATION.md` — especially **§6**, which states
   the Phase-3 target law σ_pq ∼ σ_* · q^{−α/γ} and explains why the plan's
   original σ_pq ∼ σ_*/√q (which assumed α = 1, γ = 2) is **retired**.
5. `regime-tests/CANARD_BLOWUP.md` **§4** — the 2D σ_* = C_q √ε λ^{1/2}
   derivation by accumulated Brownian variance (η = σ/√ε, W_* = λ^{2/3},
   V_* = λ^{1/3}, T_window = λ^{−1/3}, η² · T_window ∼ V_*² ⇒ η_* = C_q λ^{1/2}).
   **Task 3B is the 3D folded-node analogue of this calculation.**
6. `regime-tests/mmo_fhr_staircase.py` and `regime-tests/mmo_k2_return_map.py`
   — the deterministic FHR integrator, the `rotation_number` peak-counter,
   and the per-episode SAO machinery. **You will reuse these wholesale**;
   the only change for 3A is the integrator (RK4 → Euler–Maruyama with a
   σ·dW term on v).

## 3. What's already established (rely on without re-deriving)

**Working FHR model** (the one that actually produces MMOs — the canonical
Rinzel form does not; see `MMO_PHASE1_5.md` §1):

```
dv = ( v − v³/3 − w + y + I ) dt  +  σ dW       (fast, noisy)
dw = ε ( v + a − b w ) dt                        (slow)
dy = ε δ ( c − v ) dt                            (slow modulator)
```

with `(a, b, ε, δ, I) = (0.7, 0.8, 0.08, 0.2, 0.30)`, and **c the
bifurcation control**: sweeping `c ∈ [−0.96, −0.70]` traverses the Farey
staircase. Noise is **degenerate** — only on v, the same convention as
every other chapter.

**Keep δ = 0.2 as the primary operating point — do not switch to δ = 1.0.**
Phase 3 validates σ_pq ∼ σ_*·q^{−α/γ}, and **every input to that law must
come from one operating point**: α = 1.55, the K2 funnel-filling f(c), and
κ ≈ 2π² (which F(μ) and the dissolution geometry are built on) were all
pinned at δ = 0.2. Running the noise sweep at δ = 1.0 while importing α and
f(c) from δ = 0.2 would quietly mix regimes — exactly the inconsistency this
project exists to catch. δ = 0.2 is already genuinely 2-slow (O(1) slow-rate
ratio, real folded node — `MMO_TIMESCALE_CHECK.md`); δ = 1.0 is not *more*
defensible, just the other edge of the same regime, so there is no
cleanliness to buy by switching.

**δ = 1.0 has one real job: the δ-robustness cross-check** — testing whether
F(μ) and γ (hence the σ_pq law itself) are δ-robust, the same way μ ∝ (c+1)
was shown δ-robust. Run it *after* the δ = 0.2 chain, as a check, not a
substitute. **Fallback:** if the δ = 0.2 noise sweeps prove too expensive to
resolve σ_pq per plateau (at δ = 1.0 the y-drift is 5× faster, MMO periods
~5× shorter, sweeps correspondingly cheaper), the clean move is *not* to mix
— re-pin α and f(c) at δ = 1.0 (re-run the fine α sweep; the MMO band sits
nearer c ≈ −0.93 there) and run the **whole chain** self-consistently at
δ = 1.0. The timescale check guarantees the folded-node structure survives
the move.

**Folded node.** Desingularised reduced flow at the fold v = −1 has
Jacobian `J = [[1+δ, −b],[2bδ(c+1), 0]]`, so `μ = λ_weak/λ_strong ≈
2bδ(c+1)/(1+δ)² ∝ (c+1)`. Over the MMO band `μ ∈ [0.03, 0.16]` — small,
the **many-secondary-canard** regime.

**Deterministic structure.** Each `L¹Sˢ` episode = one global return (L
spike + re-injection toward the fold) + s small oscillations (SAOs)
through the funnel. So `ρ = 1/(s+1)`, `q = s+1`, `p = 1` for the principal
sequence. Wechselberger ceiling `s_max = (1−μ)/(2μ)` (an upper bound, α_ceiling = 2);
realised `α = 1.55 ± 0.06` (bracket [1.48, 1.60]); funnel-filling
`f(c) = (2/κ)·ln(a_max/a_min)/(1−μ)` with rotation constant κ ≈ 2π²
(measured 19.6 ± 1.8, constant across c; SAO amplitudes grow by a per-turn
factor `R = e^{κμ}`).

**The Phase-3 target law** (from `MMO_ALPHA_DERIVATION.md` §6):

```
plateau (p,q) dissolves when  δρ_noise(σ) ∼ Δ_pq ∼ q^{−α},
with  δρ_noise ∼ σ^γ   ⇒   σ_pq ∼ σ_* · q^{−α/γ},
and   σ_*(μ, ε, δ) = C_q · √ε · F(μ)   (3D analogue of canard's √ε λ^{1/2}).
```

α ≈ 1.55 is **known**. The two unknowns Phase 3 supplies are **F(μ)** (the
folded-node geometry factor) and **γ** (the rotation-spread exponent).
Then β = α/γ and the staircase dissolution is fully predicted. Order-of-magnitude:
ε = 0.08 ⇒ √ε ≈ 0.28; canard C_q ≈ 8–10 (full FHN); μ ∈ [0.03, 0.16];
expect σ_* somewhere in O(0.01–0.1), so the σ-sweep should span
σ ∈ [1e-4, ~0.2] log-spaced.

## 4. The task

Two halves, in this order.

**3A — numerical σ_pq sweep (measure-first).** Add degenerate noise to the
FHR integrator and measure, for each low-q plateau, the σ at which it
dissolves. Extract: σ_pq(q), the dissolution **order** (predicted: high-q
first), and γ **directly** from the rotation-number spread δρ(σ) at a fixed
in-plateau c.

**3B — analytical σ_*(μ,ε,δ) and γ.** Extend the canard chapter's
accumulated-Brownian-variance argument into the K2 chart of the
folded-node blow-up to get F(μ) and γ in closed form. Predict σ_pq ∼ σ_* ·
q^{−α/γ} and check against 3A.

The chapter's deliverable is the **validated** law: measured (3A) and
derived (3B), agreeing within error, with the same universal-vs-model-specific
split as the canard/tonic chapters (universal exponent law + model-specific
prefactor).

## 5. Specific steps — 3A (numerical, do this first)

1. Read the listed files (§2). Re-run `mmo_fhr_staircase.py` to confirm the
   deterministic staircase reproduces (ρ locking at 1/2, 1/3, 2/3, …) and
   record the deterministic plateau windows [c₁(p,q), c₂(p,q)] and widths
   Δ_pq for q = 2, 3, 4, 5.
2. Build `mmo_noise_staircase.py`: copy the FHR RHS and `rotation_number`
   counter from `mmo_fhr_staircase.py`, replace the RK4 step with
   **Euler–Maruyama**:
   ```
   v += dt*f_v + σ*sqrt(dt)*randn();   w += dt*f_w;   y += dt*f_y
   ```
   (noise on v only). Additive noise ⇒ EM is strong order 1.0, adequate;
   the real constraint is **dt small enough that the drift near spikes is
   resolved and EM discretisation noise ≪ σ**. Use dt ≤ 0.01 and verify
   dt-convergence (halve dt; σ_pq must not move — this is the chapter's
   chief numerical hazard, see §7).
3. **Robust event counter under noise.** The deterministic counter splits
   peaks at v = 0 (spike vs subthreshold). Under noise, near-threshold loops
   flip class spuriously. Replace with an amplitude classifier **with
   hysteresis**: a spike requires v to exceed a high threshold (e.g. +1.0)
   and reset below a low one; everything else with a local v-max in the
   subthreshold band counts as an SAO. Validate it reproduces the σ = 0
   counts before trusting it at σ > 0.
4. **σ_pq, primary definition (locked-fraction).** Fix c at the plateau
   centre c*(p,q). For a grid of σ (log-spaced), run N ≥ 200 independent
   seeds; classify each realisation's ρ. Define σ_pq as the σ where the
   locked fraction (realisations still at ρ = p/q within tolerance) falls
   below 50%.
5. **σ_pq, cross-check definition (variance crossing).** σ_pq where the
   ensemble std(ρ) first exceeds Δ_pq/2 (half the deterministic plateau
   width). The two definitions should agree to O(1); report both.
6. **γ directly.** At a fixed in-plateau c*, measure δρ_noise(σ) = std(ρ)
   over seeds for σ **below** dissolution, and fit δρ ∼ σ^γ. (Canard
   intuition suggests γ near 1 or 1/2 — measure, don't assume.)
7. **Fit σ_pq vs q.** Log-log fit σ_pq(q) over the resolved plateaus →
   exponent β. With α = 1.55 known, infer γ = α/β and check consistency
   with the directly-measured γ from step 6.
8. **Dissolution order.** Sweep the full c-staircase at increasing σ and
   record which plateaus vanish first. Prediction: high-q (deep c → −1,
   small μ, narrow plateaus) dissolve first; the staircase collapses toward
   the 1/2 and 1/1 steps. Produce this as the qualitative signature figure.

## 6. Specific steps — 3B (analytical, after 3A gives the target)

The goal is F(μ) and γ. **Do not reconstruct Wechselberger (2005) §4's
parabolic-cylinder inner solution from memory** — same intellectual-honesty
rule as the K2 prompt. Use the two-path split.

### Path A (preferred if tractable)
Full K2-chart calculation. In the folded-node blow-up the noise picks up
the same effective amplitude η = σ/√ε as the canard chart. Identify the K2
analogues of the canard window: the trajectory threads s ≈ s_max secondary
canards separated by sectors whose width shrinks with μ. Set up the
accumulated-variance criterion for noise to push the trajectory **across one
secondary-canard sector** (i.e. change s by ±1, hence change ρ): η² ·
T_funnel ∼ (sector width)². Read off η_*(μ) = F(μ) and the dependence of
δs (hence δρ) on σ above threshold, giving γ. If you can write the sector
width and T_funnel in μ explicitly via the parabolic-cylinder count — do
it; otherwise switch to Path B.

### Path B (fallback — likely the operationally correct path)
Route around the inner solution using the **measured** rotation map from
`MMO_K2.md`: SAO amplitudes grow geometrically by `R = e^{κμ}` per turn,
κ ≈ 2π². So the funnel's secondary-canard boundaries sit at
geometrically-spaced amplitudes, gap ∼ μ in log-amplitude. A noise-induced
log-amplitude jitter that accumulates over the funnel passage (canard-style
accumulated Brownian variance, η = σ/√ε) crosses one gap when the
jitter ∼ μ; this sets η_*(μ) = F(μ) and δs ∼ (jitter/gap), giving δρ ∼ σ^γ.
This reuses the deterministic global return numerically (already in
`mmo_k2_return_map.py`) and Wechselberger's count symbolically — the same
Path-B logic that closed the deterministic side.

Either path: **compose F(μ) and γ with the known α = 1.55** to predict
σ_pq ∼ σ_* · q^{−α/γ}, and overlay on the 3A measurement.

## 7. Falsifiability criteria

- **σ_pq(q) is a clean power law, dissolution is high-q-first, and the
  3B-derived (F(μ), γ) reproduce the measured β = α/γ within error** →
  the noise law is validated; the chapter closes. Write up §10.
- **σ_pq(q) is a power law and dissolution order is right, but 3B misses
  the prefactor / γ** → the mechanism (Farey dissolution) holds but the
  K2-chart FW action needs Path A; frame as the next step (parallels the
  deterministic side, where Path B closed the mechanism and Path A — the
  closed-form a_min(c) and first-principles κ — was left open).
- **σ_pq(q) is NOT a power law, OR dissolution order is not monotone in q**
  → the Farey-staircase-dissolution picture itself needs revisiting (e.g.
  noise-induced mixing between non-adjacent plateaus, or the small-μ
  many-canard regime breaking the q-indexing). Document precisely what the
  data shows.
- **Watch for the numerical-noise false positive:** if σ_pq does not survive
  halving dt, the "dissolution" is discretisation noise, not physics —
  re-run before claiming any threshold (§9.1).

## 8. Hazards

1. **Numerical noise vs real noise (the chief hazard).** The v-equation is
   stiff near spikes; a too-coarse dt injects discretisation noise that
   mimics dissolution and biases σ_pq downward. Every σ_pq must be shown
   dt-converged (halve dt, σ_pq stable). This is `MMO_FHR_PLAN.md` §6
   hazard 2 — treat it as the limiting error, as SAO-count dt-sensitivity
   was for the deterministic α.
2. **SAO/spike misclassification under noise.** Near-threshold loops flip
   class when noise jiggles a peak across the deterministic v = 0 split.
   The hysteresis classifier (§5 step 3) is mandatory; validate it against
   σ = 0 counts.
3. **Narrow high-q plateaus.** Deep in the band (c → −1, μ small) plateaus
   are a few c-grid spacings wide and the SAO trains are long. Resolving
   their dissolution needs a fine c-grid AND many seeds; these are the
   points where both 3A statistics and the 3B Path-B gap estimate are
   least accurate (mirror of the deterministic f(c) being worst at c = −0.94).
4. **σ_* is an ε^{1/2} scale, not O(1).** With √ε ≈ 0.28 and small μ, σ_*
   may be ≲ 0.05. Make sure the σ-grid resolves the threshold rather than
   stepping over it; pilot a coarse sweep first to locate the knee.
5. **Don't over-reach Path A.** If reconstructing the parabolic-cylinder
   inner solution runs past a page, stop and use Path B — reconstructed
   math has the same status as numerics here; claims must come from proofs
   that exist or computations that ran.

## 9. Deliverables

1. `MMO_NOISE.md` — Phase-3 writeup. Structure parallel to `CANARD_BLOWUP.md`
   / `MMO_K2.md`: §1 question, §2 setup (noisy FHR + σ_pq definitions),
   §3 measured σ_pq(q) + γ + dissolution order (3A), §4 derived F(μ) + γ and
   the σ_pq ∼ σ_*·q^{−α/γ} comparison (3B), §5 what's established / open,
   §6 hazards. Be explicit about which 3B path (A or B) was used and the
   dt-convergence evidence behind each σ_pq.
2. `mmo_noise_staircase.py` — self-contained EM noise sweep; reads the
   deterministic plateau windows (or recomputes them), writes outputs to
   `data/` and `figures/`. Reuses the FHR RHS + counter from
   `mmo_fhr_staircase.py`.
3. `figures/mmo_noise_staircase.png` — (top) σ_pq vs q log-log with the
   3B-predicted line overlaid; (middle) ρ-distribution broadening with σ at
   fixed c; (bottom) the staircase at increasing σ showing high-q-first
   collapse.
4. Update the `MMO_FHR_PLAN.md` banner with a new line:
   "(7) Phase 3 noise dissolution: [result] — σ_pq ∼ σ_*·q^{−α/γ},
   measured β = …, γ = …, F(μ) = …; [closed / what's missing]." Be honest.

## 10. If the calculation closes

The chapter is then complete. Final integration:

- `MMO_CHAPTER.md` — consolidate Phases 1.5 / 2 / 3 into one chapter doc
  parallel to `CANARD_BLOWUP.md` and `TONIC_PHASE.md`.
- Flip the **README §1 and §7 MMO rows** from "deterministic mechanism
  derived + confirmed; σ_pq pending" to **"derived + validated"**, with the
  σ_pq ∼ σ_*·q^{−α/γ} law and the measured (F(μ), γ).
- Add the MMO row's universal-vs-model-specific split to the §7 structural
  finding and §9.3 methodological insight: universal ceiling α = 2 and
  universal noise mechanism σ_pq ∼ σ_*·q^{−α/γ}; model-specific funnel-filling
  f(c) and prefactor C_q via FHR's global return. This is the fifth instance
  of "noise destroys deterministic structure where transverse stability is
  weakest, laws universal with model-specific prefactors."

## 11. If the calculation doesn't close

Document precisely what fails. Likely scenarios:

- **Path B's gap heuristic is wrong** — the secondary-canard sector width
  isn't simply ∼ μ in log-amplitude; identify the right scaling, flag Path A.
- **Path A genuinely required** — the K2 inner solution carries
  rotation-spread information the bound + measured κ don't capture; frame as
  the next session, parallel to the deterministic Path-A open items
  (closed-form a_min(c), first-principles κ = 2π²).
- **σ_pq not a power law in q** — the Farey-dissolution framing needs
  revisiting in the small-μ many-canard regime; this is itself a result and
  a sharp supervisor question.

In all cases the framework is sharper and the supervisor meeting becomes
"which direction next," not "what's the result."

## 12. Reading priority if time-constrained

Three things before starting:
1. `MMO_K2.md` (the deterministic object you add noise to; the κ, f(c)).
2. `MMO_ALPHA_DERIVATION.md` §6 (the σ_pq ∼ σ_*·q^{−α/γ} target law).
3. `CANARD_BLOWUP.md` §4 (the 2D σ_* accumulated-variance derivation that
   3B generalises).

Everything else is supporting context to read as needed. Start 3A
immediately after — it reuses existing code and hands 3B its target.
