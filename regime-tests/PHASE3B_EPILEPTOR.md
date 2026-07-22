# Phase 3b — does the early-warning transfer to the Epileptor? A domain-of-validity result

*Ports the Phase-3 pipeline to the published **Epileptor** (Jirsa, Stacey, Quilichini, Ivanov,
Bernard, *Brain* 2014). Code: `epileptor_phase3b.py` (+ `epileptor_figure.py`). Figure:
`figures/epileptor_phase3b.png`. Tags: **[R]** proved / **[N]** numerical / **[H]** heuristic.*

```
python3 epileptor_phase3b.py sims     # integrate + cache the Epileptor run (~10 s)
python3 epileptor_phase3b.py analyze  # transfer diagnostics
python3 epileptor_phase3b.py fig      # -> figures/epileptor_phase3b.png
```

## Why the Epileptor is the right stress test

The Epileptor's seizure-like event is bracketed by the two **Type-I** bifurcations our classifier
knows — a SNIC at onset and a (nominal) homoclinic at offset — with the slow permittivity variable
`z` as an **internal, autonomous** drift carrying the fast subsystem between them. That is a stronger
test than Phase 3a's externally-imposed current ramp: the parameter drift is endogenous, and the
offset is a different Type-I bifurcation from the SNIC we already did.

## What we find (canonical `x0=-1.6`, 9 seizures, T=16000)

**[N] The discharge train is spike-and-wave, and is clean only at the complex level.** 20% of raw
discharge ISIs are `< 1.0` (the fast intra-complex spikes of the spike-and-wave structure, from the
`x2` subsystem) — a multimodal train that is **not** a first-passage sequence. Merging sub-discharges
into spike-and-wave **complexes** (the dominant rhythm) recovers a **unimodal Type-I** ISI shape: the
classifier calls SNIC with residual `W1 = 0.014`, comparable to genuine SNIC/homoclinic references
(`0.004–0.011`). So on the complex train the inversion is **well-posed and not producing garbage**.

**[N] But `ν̂` never warns.** Across the seizure (onset→offset) the windowed `ν̂` sits at
`2.86 / 2.99 / 2.98 / 2.76` — **high (~2.9, far from threshold) and flat**, never approaching the
Phase-3a alarm (`0.4`), let alone 0. There is **no `ν→0` early-warning of seizure offset.**

**[N] Because the offset is not a period-diverging bifurcation of the discharge cycle.** The
inter-complex period ratio `mean(last 5)/median = 0.86 < 1`: the discharges **accelerate** into an
abrupt, multi-variable termination rather than slowing toward a homoclinic (which would *diverge*,
ratio `≫1`). There is simply no diverging ISI precursor to detect — `z` sweeps the system out of the
ictal state faster than the discharge period can blow up.

## The domain-of-validity statement (the real Phase-3b contribution)

**[H] The edge-law early-warning requires an *adiabatic* crossing.** The universal first-passage law
(CV→0.57, period divergence, the skewed quartic shape) only develops if the system **dwells near the
bifurcation** long enough — i.e. the control parameter must cross slowly relative to the (diverging)
event timescale. Phase 3a's noise-driven ramped SNIC satisfies this and `ν̂→0` warns with 0% false
alarms. The Epileptor's fast relaxation-oscillator offset does **not**: its discharge period is set by
the fast `(x1,y1)` timescale, not by proximity to the offset, and the termination is non-adiabatic.

This is a clean boundary on the moonshot, established on a *real* published model:

| | Phase 3a (ramped SNIC) | Phase 3b (Epileptor offset) |
|---|---|---|
| approach to bifurcation | adiabatic (slow ramp) | non-adiabatic (fast `z`-driven exit) |
| discharge period near transition | diverges (`J(ν)→∞`) | ~constant, then accelerates |
| ISI in universal regime? | yes | no (complex-level Type-I, but far-from-threshold) |
| `ν̂→0` early-warning | **works** (96% detect, 0% FA) | **does not transfer** |

**[N] Constructively, the method does not misbehave.** Reading a developed seizure as a *robust Type-I
oscillation far from its bifurcation* (low `W1`, large `ν̂`) is the **correct** description — a
mid-seizure limit cycle is not near its birth/death point in the rescaled sense — and the inversion
**does not false-alarm**. The honest conclusion is about *applicability*, not correctness: the tool
targets the **noise-driven pre-ictal / interictal proximity-to-onset** regime, where individual events
are genuine first-passage times, rather than fully-developed seizure termination.

## Honest limits / open (Phase 3b-cont)

- **No clean noise-driven excitable Epileptor control yet.** The natural positive control on this model
  — freezing `z` just past onset and eliciting noise-induced first-passage events — was hard to access
  at small noise in the time available; left open. (Phase 3a *is* the noise-driven-SNIC positive
  control, just not inside the Epileptor.)
- **Single operating point** (`x0=-1.6`). The Saggio–Jirsa onset/offset taxonomy contains genuinely
  period-diverging homoclinic offsets; testing whether `ν̂→0` fires for *those* classes would sharpen
  the adiabaticity boundary into a per-bifurcation map.
- The pipeline itself (stream ISIs → complex-merge → window → invert → alarm) is model-agnostic and
  ran unchanged; only the verdict is model-specific.
