# MMO Phase 1 — findings and a load-bearing redirect

**Status:** Phase 1 (Tasks A–C) executed. The headline finding is a **correction
to the plan's premise**: the deterministic `L^a S^b` devil's staircase the plan
(MMO_PLAN.md §1–§2) attributes to autonomous 2D FHN **cannot exist there** — it
is precluded by planar topology. Genuine MMOs and their staircase do appear once
a third slow variable is added (3D FitzHugh–Rinzel), which is the
topologically-correct home for the chapter and is already on the project roadmap.

Companion scripts: `mmo_2d_check.py` (Phase 1a), `mmo_fhr_staircase.py` (Phase 1b).
Figures: `figures/mmo_2d_check.png`, `figures/mmo_fhr_staircase.png`.

---

## 1. The premise fails in 2D — and it has to

MMO_PLAN.md §1/§2.1 assumes that inside the canard explosion window the
**autonomous 2D** FHN limit cycle realises `L^a S^b` patterns (a large spikes +
b small canard loops per period) tracing a Farey devil's staircase. Two
independent obstructions:

**(a) Planar topology (the decisive one).** An autonomous planar flow has unique
trajectories, so a stable limit cycle is a simple closed (Jordan) curve. A curve
that executes "a large loops then b small loops" in one period (a+b > 1) must
self-intersect — forbidden. Hence a 2D FHN attractor can carry **at most one
loop per period**: it is a pure spike train (all L) or a pure subthreshold/canard
cycle (all S), never a mixture. There is no room for an `L^a S^b` symbol sequence
and therefore no rotation number to mode-lock. Genuine MMOs require ≥3 phase-space
dimensions (a folded node supplies the rotational SAOs) or non-autonomous forcing.

**(b) The window is exponentially thin.** Even the S↔L amplitude transition (the
canard explosion) occupies an I-interval of width `exp(−c/ε)` — the README notes
this repeatedly as below numerical resolution. So even the single-loop transition
the plan hoped to resolve as a staircase is not numerically accessible in I.

**Numerical confirmation (`mmo_2d_check.py`).** Sweeping I across `[I_H1−0.02,
I_H1+0.12]` at ε = 0.04 and 0.08, classifying every steady-state oscillation as L
(spike, v>0) or S (subthreshold): **0 of 18 I-values** show steady-state L and S
coexistence at either ε. The attractor is a fixed point below `I_H1` and a pure
relaxation cycle (all L, amplitude CV = 0) above it; `ptp(v)` jumps 0 → 3.9 in a
single step. ρ(I) is a 0→1 step with **no intermediate plateaus**. The plan's
deterministic staircase is simply not there.

## 2. The staircase is real in 3D FitzHugh–Rinzel

The minimal FHN-family carrier of `L^a S^b` MMOs adds one slow variable:

```
v' = v − v³/3 − w + y + I        (fast)
w' = ε (v + a − b w)             (recovery)
y' = ε δ (c − v)                 (slow drive; δ<1 ⇒ three timescales)
```

with the project's constants `a=0.7, b=0.8, ε=0.08`, here `δ=0.2, I=0.30`. The
slow `y` is a moving effective current `I_eff = I + y`; the global return drags
`I_eff` back and forth across the (v,w)-Hopf, so each pass produces a spike (L)
followed by a burst of subthreshold canard oscillations (S) — genuine MMOs.

**Numerical result (`mmo_fhr_staircase.py`).** Sweeping `c` (which sets where the
slow drift settles, hence the SAO count) gives a clean **devil's staircase** in
the rotation number `ρ = L/(L+S)`:

```
c > −0.70   : ρ = 1        (pure spiking, L)
c ≈ −0.72…−0.74 : ρ = 2/3  (L² S¹ plateau)
c ≈ −0.76…−0.78 : ρ = 1/2  (L¹ S¹ plateau)
c ≈ −0.82…−0.83 : ρ = 1/3  (L¹ S² plateau)
… descending through 1/4, 1/5 …
c < −0.96   : ρ = 0        (pure subthreshold, S)
```

The sample trajectory at the 1/2 plateau shows the textbook MMO signature: tall
spikes (v ≈ 1.7) interleaved with subthreshold S-loops (v ≈ −0.8). Plateau
widths **decrease with denominator**: Δ(1/2) ≈ 0.025, Δ(2/3) ≈ 0.020, Δ(1/3) ≈
0.015, with q = 4,5 plateaus at the grid-resolution floor (~0.010).

**On α.** A coarse fit gives `Δ_pq ~ 1/q^α` with `α ≈ 1.1`, **below** the
circle-map prediction `α ∈ [2,3]`. This is *not* a reliable measurement: only
~4 distinct q are resolved on a `Δc ≈ 0.005` grid, and the high-q widths are at
the resolution floor. The defensible Phase-1 statement is qualitative — **the
staircase exists and Δ_pq falls with q** — and pinning α is exactly the Phase-2
analytical job (return map via blow-up) plus a finer numerical sweep.

## 3. What this means for the chapter (the redirect)

The MMO chapter is viable and interesting, but its model must change from
autonomous 2D FHN to a 3-timescale carrier. Two options:

| Option | MMO mechanism | Pros | Cons |
|---|---|---|---|
| **A. 3D FitzHugh–Rinzel** (recommended) | folded node → autonomous `L^a S^b` | autonomous (keeps the project's framing); canard-chapter blow-up enters directly as the (v,w) fold passage during each S; on the roadmap (FitzHugh–Rinzel for bursting); matches the plan's L^a S^b language exactly | one more variable; the return map is 2D→1D reduced via the folded node |
| **B. Periodically-forced 2D FHN** | circle-map p:q mode-locking | stays 2D; matches the plan's *cited* circle-map literature (Boyland, Glass–Mackey); cleanest devil's staircase | non-autonomous; changes the question from `L^a S^b` to forcing:spike locking; less connected to the canard σ_* |

**Recommendation: Option A (FitzHugh–Rinzel).** It is the faithful carrier of the
plan's `L^a S^b` / canard-SAO picture, stays autonomous, and keeps the
cross-chapter unification intact: each S is a (v,w) fold passage governed by the
**same Krupa–Szmolyan blow-up** as the canard chapter, so `σ_pq` should still tie
back to `σ_* = C_q √ε λ^{1/2}` (MMO_PLAN.md §3). The plan's four analytical
engines survive the move — and the **entry–exit engine (Kaklamanos–Kuehn–Popovic–
Sensi 2025)** becomes *more* central, since the folded-node return map is exactly
an entry–exit relation across the canard fold.

## 4. Reframed task sequence

Phases 2–4 of MMO_PLAN.md carry over with the model swapped to FHR:

- **Phase 1 (done here):** detector + staircase + coarse Δ_pq, **in FHR** (Tasks
  A–C), plus the 2D-obstruction result that justifies the model choice.
- **Phase 1.5 (new, ~days):** finer `c`-sweep (Δc ≈ 0.001, longer records) to
  resolve higher-q plateaus and measure α properly; map the (δ, I) region where
  MMOs live and confirm the Farey/mediant ordering.
- **Phase 2:** folded-node return map via the blow-up; derive Δ_pq ~ 1/q^α. The
  entry–exit framework is the natural engine. (Unchanged in spirit; new model.)
- **Phase 3:** degenerate noise σ on v; noise-broadened staircase; σ_pq; the
  σ_pq ~ σ_* q^{−β} unification check. (Unchanged.)
- **Phase 4:** `MMO_CHAPTER.md` + README integration.

## 5. The supervisor question this sharpens

MMO_PLAN.md §7 asked three generic questions. Phase 1 turns the first into a
sharp, decision-ready one:

> **The deterministic `L^a S^b` staircase is topologically impossible in
> autonomous 2D FHN; it appears cleanly in 3D FitzHugh–Rinzel (folded node).
> Do we (A) take the chapter into FitzHugh–Rinzel — keeping it autonomous and
> letting the entry–exit/folded-node machinery carry the return map — or (B)
> stay 2D with periodic forcing and the circle-map machinery? And for (A), is
> the entry–exit framework the right engine for the folded-node return map, or
> is there a cleaner singular-limit reduction?**

This is a better-posed pre-MMO meeting than the original plan, because it is
backed by the obstruction result and a working staircase rather than an
untested premise.

## 6. Reproduce

```
python3 regime-tests/mmo_2d_check.py        # 2D: no staircase (obstruction)
python3 regime-tests/mmo_fhr_staircase.py   # 3D FHR: genuine MMO staircase
```
Outputs land in `regime-tests/results/mmo/` and `figures/`.
