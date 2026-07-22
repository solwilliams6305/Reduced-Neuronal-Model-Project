# MMO chapter — Phase 2 K2 calculation: handoff to Opus

This is a research handoff. You are picking up the FHN-MMO chapter at
the point where Phase 1.5 has cleared all prerequisites and the
remaining analytical step is one well-defined calculation. The task is
to compute the FHR funnel-filling f(c) from the K2-chart return map
matched to the global return, and check whether the result reproduces
the measured exponent α = 1.55 ± 0.06.

## 1. Brief context

The broader project is a regime-by-regime theory of how noise affects
stochastic FitzHugh–Nagumo dynamics. Four chapters have closed:
excitable (σ_crit = √(ε/C)), canard (σ_* = C_q · √ε · λ^{1/2}, validated
in both normal form and full FHN, plus Van der Pol cross-model), tonic
(CV = σ · A_mid with matched-asymptotic adjoint Floquet at A_mid ≈ 0.39),
and resonator (subcritical-Hopf B ∼ δ² corner law + phase-gated
commitment). The current chapter is **MMO** — the mode-locked
L^a S^b patterns in 3D autonomous FHR (FitzHugh–Rinzel), with σ_pq the
noise that dissolves each plateau p/q in the Farey staircase.

Two earlier framings of the chapter were retracted in Phase 1 / 1.5:
the 2D FHN MMO scheme failed at Jordan-curve topology (autonomous 2D
flows can't support L^a S^b deterministically), and an initial α ≈ 1.0
fit was a coarse-grid artefact. The chapter now sits on a five-correction
trail (Jordan obstruction → bound-not-equality on Wechselberger's SAO
count → α=1 retracted → timescale-direction inverted → α tightened to
1.55 ± 0.06). All five are documented in `MMO_FHR_PLAN.md` banner.

## 2. Read first (in order, before computing anything)

1. `README.md` §1, §6 — master σ_crit table + project framing.
2. `regime-tests/CANARD_BLOWUP.md` §1–§6 — the 2D matched-asymptotic
   methodology this generalises. The Krupa–Szmolyan blow-up there is
   the K2 chart in the folded-saddle-node limit μ → 0.
3. `regime-tests/TONIC_PHASE.md` §14, §19 — the corrected −2V adjoint
   inner equation. Demonstrates the "same blow-up applied to the dual
   equation" pattern that this chapter extends.
4. `regime-tests/MMO_FHR_PLAN.md` — full plan with Phase 1–1.5 trail
   in the banner.
5. `regime-tests/MMO_PHASE1_5.md`, `MMO_ALPHA_DERIVATION.md`,
   `MMO_TIMESCALE_CHECK.md` — diagnostic Phase 1.5 documents. Read
   especially `MMO_ALPHA_DERIVATION.md` for what's already analytically
   established about μ ∝ (c+1) and the structural α ∈ [1.5, 2] bracket.

## 3. What's already established

The Phase 1.5 outputs you can rely on without re-deriving:

**Working FHR model** (the one that actually produces MMOs):

```
dv/dt = v − v³/3 − w + y + I,
dw/dt = ε ( v + a − b w ),
dy/dt = ε δ ( c − v ),
```

with (a, b, ε, δ) = (0.7, 0.8, 0.08, 1.0). I held at a value that puts
the system in the MMO band; c is the bifurcation control, sweeping
c ∈ [−0.96, −0.7] covers the deterministic Farey staircase. δ = 1.0
is the cleanest 2-slow limit; do NOT use the canonical Rinzel form
with δ → 0 (that's three-timescale, Krupa–Popovic–Kopell territory,
unnecessary for this task).

**Folded node geometry.** The desingularised reduced flow at the v = −1
fold gives a folded node for all δ ∈ [0.05, 1.0]. The eigenvalue ratio
μ = λ_weak / λ_strong satisfies μ ∝ (c + 1) at fixed δ. Over the MMO
band μ ∈ [0.03, 0.16] — small, in the many-secondary-canard regime.

**Wechselberger ceiling.** s_max = (1 − μ)/(2μ) ⇒ since μ ∝ (c+1),
the ceiling SAO count s_max ∝ 1/(c+1), giving the universal
ceiling α = 2.

**Measured realised α.** Tightened on a 22-point c-grid:
α = 1.55 ± 0.06 (full bracket [1.48, 1.60], deep-subrange 1.58 ± 0.06).
Realised s_obs(c) ∼ (c+1)^{−1.74} vs ceiling s_max(c) ∼ (c+1)^{−1.0}.

**Funnel-filling.** f(c) := s_obs(c) / s_max(c) rises from 0.21 at
c = −0.7 to 0.44 at c → −1. This is the quantity to derive.

**Noise corollary** (for after deterministic closes):
σ_pq ∼ σ_* · q^{−α/γ}, with σ_*(μ, ε, δ) = C_q · √ε · F(μ) the 3D
analogue of the canard chapter's σ_* law.

## 4. The task

Compute f(c) — the realised-to-ceiling SAO-count ratio — from
deterministic FHR via the K2-chart return map matched to the global
return. Predict α from f(c) · s_max(c), compare with measured 1.55.

## 5. Two-path strategy

You should not attempt to reconstruct Wechselberger (2005) §4's
parabolic-cylinder inner solution from memory. That's a 40-page proof
and reconstructing it without the paper risks fabrication. Use the
following principled split:

### Path A (preferred if tractable)

Full K2 matching. Use Wechselberger's K2 normal form dV/dT = V² − W,
dW/dT = μV − (μ+1)/2 (or equivalent; double-check the convention
against the canard chapter's normal form which is the μ → 0 limit).
Set up K1 (incoming attracting manifold) and K3 (outgoing fast jump)
matching boundary conditions. Compose with FHR's global return to
close the loop.

Realistically, the K2 inner solution at finite μ involves parabolic
cylinder functions and the secondary canard count. If you can write
down the K2 Poincaré map s(r_in, μ) explicitly using these special
functions — go for it. If you cannot, switch to Path B.

### Path B (fallback — likely the operationally correct path)

Treat K2 as a black box that produces s_max(μ) by Wechselberger's
theorem, then compute the realised injection distribution P(r_in | c)
from FHR's deterministic global return numerically. The funnel-filling
is then f(c) = ⟨s(r_in, μ(c))⟩_P / s_max(μ(c)).

The model for s(r_in, μ) inside the funnel can be heuristic:
trajectories landing in the inner funnel (close to the folded node)
hit the full s_max; trajectories landing on the outer rim do fewer
secondary canards. Parameterise this as s(r_in, μ) ≈ s_max · (1 −
(r_in / r_funnel)^p) or similar, fit p to the data, see if it
gives α = 1.55 ± 0.06.

This routes around the K2 inner solution by computing the global
return numerically (which is just FHR integration) and using
Wechselberger's bound symbolically. It's the operationally tractable
path and probably what closes the chapter.

## 6. Specific steps

1. Read the listed files (§2). Verify the project context.
2. Re-derive μ(c) from the desingularised reduced flow (one-page
   calculation; cross-check `MMO_ALPHA_DERIVATION.md`).
3. Pick a working point: c = −0.85, δ = 1.0, I such that the MMO band
   is crossed. Verify deterministic FHR simulation produces L^a S^b
   patterns at this point and the SAO count matches Phase-1.5
   measurements.
4. Define a Poincaré section transverse to the cycle just outside K2's
   boundary (in original (v, w, y) coordinates, this is a 2D plane
   intersecting the cycle once per period; pick something like
   v = 0 with v̇ > 0).
5. Trace deterministic trajectories from many initial conditions on
   the K3-exit side of K2 back to first re-entry into K2 (the K1 side).
   This is just forward FHR integration. Tabulate the K1-entry
   distribution P(r_in | c) for c ∈ [−0.95, −0.70].
6. Compute the realised SAO count s_obs(c, r_in) for each (c, r_in)
   on the grid. This is the count of small-amplitude loops between
   successive K3-exits.
7. Match to Wechselberger's bound: at each c, ⟨s_obs(c)⟩_P = realised
   count; s_max(c) = (1 − μ(c))/(2μ(c)) = ceiling.
8. Compute f(c) = ⟨s_obs⟩_P / s_max. Fit log(s_obs) vs log(c+1) and
   extract α.
9. Compare predicted α with measured 1.55 ± 0.06.

If you went via Path A: also derive the s(r_in, μ) closed form
analytically and compare with the numerics from step 6.

## 7. Falsifiability criteria

- α = 1.55 ± 0.1: confirmed. The mechanism is FHR's global return
  + Wechselberger's bound. Write up §10 of MMO_FHR_PLAN.md.
- α = 1.45 or 1.65: ambiguous. Re-check Phase-1.5 SAO counting at
  finer dt (factor 2 reduction) before claiming a miss. The counts
  at c → −0.95 involve long trajectories with small loops; if the
  counts are dt-sensitive, the measured 1.55 itself shifts.
- α = 1.3 or 1.8: real miss. The K2 inner solution is doing more work
  than the bound captures. Either Path B's heuristic for
  s(r_in, μ) is wrong, or the K2 inner calculation is genuinely
  necessary. Document what failed and what additional machinery
  is needed.

## 8. Hazards

1. **Many-secondary-canard regime.** μ ∈ [0.03, 0.16] means
   Wechselberger's secondary canards are 5–15 per fold passage. The K2
   matching is most delicate here; primary canard is not the
   load-bearing object. Wechselberger's theorem still applies but
   numerical evaluation of the inner solution may require careful
   asymptotics.

2. **δ = 1.0 boundary of "2-slow."** The Phase-1.5 timescale check
   verified Wechselberger applies, but δ = 1.0 is at the upper edge
   of the 2-slow regime. If you find the result is sensitive to δ
   near 1.0, the chapter may need a finite-δ correction to
   Wechselberger's theorem. Note this if it surfaces.

3. **SAO counting at deep c → −0.95.** Trajectories are long and
   loops are small-amplitude; counting can become dt-sensitive. The
   measured α = 1.55 rests on the counts staying accurate in this
   regime. If your global-return computation lands at α = 1.45, rerun
   SAO counting at half the current dt before concluding K2 is wrong.

4. **Path A's parabolic-cylinder reconstruction.** If you attempt
   Path A and find yourself reconstructing more than a page of
   Wechselberger (2005) §4 from memory, stop. Switch to Path B. The
   project's intellectual-honesty pattern is that reconstructed math
   has the same status as numerical evidence — claims need to come
   from proofs that exist or computations that ran, not from
   semi-remembered theorem statements.

## 9. Deliverables

1. `MMO_K2.md` — derivation/computation of f(c), predicted α with
   error bars, comparison with measured. Structure parallel to
   `CANARD_BLOWUP.md` (§1 question, §2 setup, §3 K2 / global return,
   §4 result, §5 hazards / what's open). Be honest about which path
   (A or B) was used.

2. `mmo_k2_return_map.py` — companion script. Self-contained, reads
   from `data/` if needed, writes outputs to `data/` and `figures/`.

3. `figures/mmo_k2_return_map.png` — the f(c) plot with measured
   α-fit overlaid. Top panel: f(c) vs c; bottom panel: log-log fit
   of s_obs(c) vs (c+1) with measured and K2-predicted lines.

4. Update banner in `MMO_FHR_PLAN.md` Phase 1.5 section with a new
   line:
   "(6) Phase 2 K2 calculation: [result]." Be honest — if the
   calculation closes (α matches), say so; if it doesn't, say what's
   missing.

## 10. If the calculation closes

The chapter has two pieces left:

- **σ_pq derivation.** Compute γ — the noise propagation factor through
  one pass of the K2-derived return map. The conjecture is
  σ_pq ∼ σ_* · q^{−α/γ} with α = 1.55 from above and γ from the FW
  action on R. This parallels the canard chapter's σ_* derivation
  exactly. ~1 session.

- **Noise-broadened staircase numerics.** Verify the σ_pq prediction
  by FHR + noise simulations: measure the σ at which each plateau
  dissolves. ~1 session.

- **Write-up.** `MMO_CHAPTER.md` with full derivation + numerical
  validation. Master README §1 table updated with the MMO row going
  from "planned (PhD scope)" to "derived + validated."

The chapter then closes to the same level as canard / tonic, with
the same universal-vs-model-specific split: universal ceiling α = 2
+ universal noise mechanism σ_pq ∼ σ_* · q^{−α/γ}; model-specific
f(c) funnel-filling + γ via FHR's global return geometry.

## 11. If the calculation doesn't close

Document precisely what fails. The likely scenarios:

- **Path B's heuristic s(r_in, μ) is wrong.** The funnel-filling
  parametrisation needs more structure. Identify what.
- **Path A is genuinely required.** Wechselberger's K2 inner solution
  carries non-trivial information beyond the bound. Frame as the
  next session's task.
- **Three-timescale matters after all.** If δ = 1.0 is too "fast" for
  Wechselberger to apply cleanly, the chapter may need to investigate
  smaller δ in the Krupa–Popovic–Kopell three-timescale regime
  separately.

In all three cases, the deterministic side of the chapter isn't fully
closed but the framework is sharper. The supervisor meeting becomes
about which direction to push next rather than "what's the chapter's
result."

## 12. Reading priority if time-constrained

If you can only read three things before starting:
1. `MMO_FHR_PLAN.md` (full banner trail)
2. `MMO_ALPHA_DERIVATION.md` (the structural α ∈ [1.5, 2] derivation
   you'll build on)
3. `CANARD_BLOWUP.md` §1–§6 (the 2D analogue that this generalises)

Everything else is supporting context that you can read as needed.
