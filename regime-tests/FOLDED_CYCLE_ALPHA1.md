# Folded limit cycle — Part 3: the α=1 coupled core (the averaged law)

**Status:** the α=1 rotating-phase Channel-A law is **derived and validated** at
leading order (Path B): fast rotation **averages the canard coefficients**
(Khasminskii), giving

```
σ_*(α=1)  =  C_q √ε₂ · √( ⟨a⟩ ⟨c⟩ / ⟨b⟩ ),      ⟨·⟩ = uniform phase average,
```

**not** any average of `G(θ)=√(ac/b)` itself, and **not** the minimum. A
discriminating test (modulate `b`) confirms it cleanly. A first swing at the A–B
fusion (§5) shows **Channel A preempts Channel B** for accessible parameters (both
flat in ω); the large-log corner and the geometry-decoupling under rotation, the
Fourier-hierarchy rigour, and the inner solution remain open.

This is the part that turns "two channels compared" into a *coupled* statement:
the phase the cycle carries rotates *through* the phase-modulated fold, and the
escape sees the rotation-averaged geometry.

Companion scripts: `alpha1_averaged_law.py` (the averaged law), `alpha1_fusion.py`
(the A–B fusion first swing). Figures: `figures/alpha1_averaged_law.png`,
`figures/alpha1_fusion.png`. Guardrail: `ALPHA1_METHODOLOGY.md`. Prereq:
`FOLDED_CYCLE_WORKED_MODEL.md` §6 (the θ-dependent fold the phase rotates through).

---

## 1. The guardrail (what was already settled)

`ALPHA1_METHODOLOGY.md` tested the natural conjecture — "fast rotation makes the
cycle escape at its *weakest* phase, `min_θ G`" — and **refuted it**: `η_*` is flat
in the rotation rate `ω` (ratio ≈ 1.01 out to ω=150; min predicted 0.50). The
mechanism: canard escape needs a finite local Riccati time `T_loc ~ 1/√(aY)`; a
fast phase dwells in the weak window only `δ/ω ≪ T_loc`, far too briefly to
complete an escape there. So **the weak phase is invisible to the rare escape**;
the result is an *average*, not the minimum. What was left open: *which* average.

## 2. Derivation — Khasminskii coefficient-averaging

The α=1 inner SDE (rescaling chart, rotating phase):

```
dr = ( b(θ) r² − a(θ) y ) dT + η dB,    θ̇ = ω,    dy = −c(θ) dT.
```

When `θ` rotates fast relative to the (r,y) escape dynamics, the slow variables
see the **time-averaged drift** (Khasminskii averaging). The drift averages
term-by-term, since `r, y` are ~constant over one fast phase-period:

```
⟨ b(θ) r² − a(θ) y ⟩_θ  =  ⟨b⟩ r² − ⟨a⟩ y,        ⟨c(θ)⟩_θ = ⟨c⟩.
```

So the effective system is an **autonomous canard with arithmetic-mean
coefficients**:

```
dr = ( ⟨b⟩ r² − ⟨a⟩ y ) dT + η dB,    dy = −⟨c⟩ dT,
```

whose escape threshold (`CANARD_BLOWUP.md`/`FOLDED_CYCLE_NOISE.md` §4) is

```
┌───────────────────────────────────────────────────────────┐
│   σ_*(α=1)  =  C_q √ε₂ · √( ⟨a⟩ ⟨c⟩ / ⟨b⟩ ).               │
└───────────────────────────────────────────────────────────┘
```

**This is "average the coefficients, then form G" — not "average G".** It is
distinct from every naive candidate:

| candidate | effective `G²` | when it would arise |
|---|---|---|
| **coefficient-avg (correct)** | `⟨a⟩⟨c⟩/⟨b⟩` | drift averages (Khasminskii) — escape can't resolve the fast phase |
| barrier-avg / RMS of G | `⟨a c / b⟩ = ⟨a⟩⟨c⟩·⟨1/b⟩` | if the *barrier* averaged (it does not) |
| arithmetic ⟨G⟩ | `⟨√(ac/b)⟩²` | — |
| min G (refuted) | `min(ac/b)` | Laplace/rare-event over phases (killed by `T_loc`) |

These coincide when only `a` or `c` is modulated (they enter `G²` linearly), so
the **clean discriminator is to modulate `b`** (it enters as `1/b`): then
`1/⟨b⟩ ≠ ⟨1/b⟩` by Jensen, and coefficient-avg vs barrier-avg separate.

## 3. The discriminating test

`alpha1_averaged_law.py`: inner SDE with `a=c=1`, `b(θ)=b₀(1+δcos2πθ)`, `δ=0.8`,
so `⟨b⟩=b₀`, `⟨1/b⟩=1/(b₀√(1−δ²))=1.667/b₀`, `b_max=1.8 b₀`. Predicted ratios of
`η_*(fast, modulated)` to `η_*(frozen, b=b₀)`:

```
coefficient-avg (THIS):  √(b₀/⟨b⟩)   = 1.000
barrier-avg (RMS of G):  √(b₀⟨1/b⟩)  = 1.291
min (refuted):           √(b₀/b_max) = 0.745
```

Measure `η_*` (where `⟨y_escape⟩` crosses 0, the canard convention) vs `ω`:

```
REFERENCE (frozen b₀):  η_* = 1.76
modulated b, ω = 0   :  ratio 1.06
modulated b, ω = 5   :  ratio 1.02
modulated b, ω = 20  :  ratio 1.01
modulated b, ω = 50  :  ratio 1.05
```

**Result: ratio ≈ 1.0–1.05, flat in ω** — sitting on coefficient-averaging
(1.00), and decisively away from barrier-averaging (1.29, ~29% off) and min
(0.745, ~26% off). The `⟨y_escape⟩(η)` curves for all ω overlay the frozen-`b₀`
reference (figure, right). So:

> **α=1 fast rotation selects coefficient-averaging:
> `σ_*(α=1) = C_q √ε₂ √(⟨a⟩⟨c⟩/⟨b⟩)`.** The escape sees the canard built from the
> phase-averaged coefficients, not the phase-averaged barrier.

This is the genuine α=1 contribution: a *derived* averaged law with the average
*pinned* by a test built to separate the candidates — the same discipline that
refuted the min route.

## 4. Why this is the coupled statement (not BG/BGK)

In the frozen-phase (α=2) chapter the threshold was `σ_*^A(θ)=C_q√ε₂√(ac/b)(θ)` —
phase-*resolved*. Here the phase *rotates through* that modulation during the
passage, and the coupling collapses the θ-resolved family into a single
rotation-averaged canard with `√(⟨a⟩⟨c⟩/⟨b⟩)`. Neither BG (no phase) nor BGK
(slow rotation, no sustained fast cycle) has this rotation-averaging of a
canard fold — it requires the sustained fast oscillation, which is exactly the
folded-limit-cycle's defining feature.

## 5. Honest status and what's open

Validated (Path B, leading order, fast ω):
- `σ_*(α=1)=C_q√ε₂√(⟨a⟩⟨c⟩/⟨b⟩)`, coefficient-averaging, flat in ω, with the
  competing averages (barrier/RMS, min) ruled out by a clean separation test.

Open (the rest of Part 3):
1. **The A–B fusion — first swing done (`alpha1_fusion.py`, `figures/alpha1_fusion.png`).**
   With BOTH channels on and the phase rotating, the key coupling is that Channel
   B's noise `η/(2πr)` *blows up as Channel A drives `r→0` at the fold*. Result
   (a=c=1, unmodulated, entry `y₀=5`): `σ_*^A ≈ 1.9` is **flat in ω** (CV 0.02),
   and at that threshold the phase is only `Var(θ)≈0.45 < 1` decohered —
   **Channel A PREEMPTS Channel B**: amplitude escape destroys the cycle before
   the phase fully randomizes, so `σ_*^B` is never reached (`Var(θ)` at escape
   saturates ~0.5 even at large η, because higher noise escapes *earlier*). So the
   coupled `σ_*` is set by Channel A, and rotation doesn't change that. The "race"
   is not two independent thresholds: the lower one (A here) fires and removes the
   substrate (the cycle) the other needs.
   **Geometry under rotation — tested, conjecture REFUTED (`folded_cycle_geometry_decoupling.py`).**
   I conjectured rotation would *decouple* the channels' geometry (A averages `⟨a⟩`,
   B's `1/r²` averages `⟨1/a⟩`), breaking the race's geometry-independence. The
   numerics kill it: modulating `a(θ)=a₀(1+δcos)`, `σ_*^A` is **flat** in δ (CV 0.00,
   A sees `⟨a⟩`) AND Channel B's variance coefficient is **flat** in δ at both ω=2
   and ω=20 — it does *not* track `⟨1/a⟩` (which would grow 1.67× at δ=0.8). Reason:
   the amplitude `r` is slaved to the *averaged* canard (`r²≈⟨a⟩y`), so Channel B's
   `1/r²` sees `1/⟨a⟩` — the **same** `⟨a⟩` functional as Channel A. **So the A-vs-B
   race stays GEOMETRY-INDEPENDENT under rotation** (robust; the regime diagram is
   defensible). **Still open:** the *large-log / small-ε₂* corner, where `σ_*^B`
   drops below `σ_*^A` — does B then preempt A (the deflated-headline "B wins"
   regime)? Untested.
2. **The `|ln ε₂|` under rotation.** `ALPHA1_METHODOLOGY.md` argues it should come
   from the ramp/turning-point (Airy) structure, not localization; to be derived.
3. **Fourier-hierarchy control (Path A backbone).** Make the averaging rigorous by
   bounding the n≠0 harmonics (damped at `n²η²/r²`) — the route to error control.
4. **Finite-ω correction / the weighting "where escape occurs".** The leading
   answer is uniform coefficient-averaging; the O(1/ω) correction (escape-location
   weighting) is open. The measured small excess (~1.05 vs 1.00) may be its onset.
5. **Rigorous inner solution (Path A, months).** Matched asymptotics through the
   JKK charts; the noisy Weber/Airy inner — the theorem-grade endpoint.

**Novelty gate (unchanged).** Everything remains standard ingredients recombined;
the deep Part-3 spend (items 1–5) should gate on the Popović reply
(`POPOVIC_NOVELTY_EMAIL.md`, `EVIDENCE_KUEHN_PIPELINE.md`).

## 6. References

- `ALPHA1_METHODOLOGY.md` (min refuted; averaging vindicated), `FOLDED_CYCLE_NOISE.md`
  (η, σ_*^A, Channel B), `FOLDED_CYCLE_WORKED_MODEL.md` (the θ-dependent fold),
  `CANARD_BLOWUP.md` (canard `C_q`), `TONIC_CMID_AIRY.md`/`TONIC_PHASE.md` (Airy/iPRC).
- Khasminskii, *Stochastic Stability of Differential Equations* (averaging);
  Kapitza effective-potential (high-frequency averaging).
- Jelbart–Kuehn–Kuntz 2024, arXiv:2208.01361 (deterministic scaffold).
