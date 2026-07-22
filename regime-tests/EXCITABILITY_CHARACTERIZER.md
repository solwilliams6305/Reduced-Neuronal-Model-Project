# Excitability class + proximity-to-threshold from extracellular spikes (the useful reframe)

*Points the edge-theory machinery at a question comp neuro actually asks. Code:
`excitability_characterizer.py` (+ `excitability_figure.py`, `allen_excitability.py`). Figure:
`figures/excitability_characterizer.png`. Tags: **[R]** proved / **[N]** numerical / **[H]** heuristic.*

```
python3 excitability_characterizer.py cohort     # build + cache synthetic cohort (~20 s)
python3 excitability_characterizer.py validate   # cohort metrics
python3 excitability_characterizer.py fig        # -> figures/excitability_characterizer.png
```

## The reframe

The moonshot was "early-warning of a tipping point" — dramatic, but oversold and hard to validate
(Phase 3/3b). The same machinery answers a question computational neuroscience genuinely has a hole
for: **is this neuron Type-I (SNIC) or Type-II (Hopf) excitable, and how far above its spiking
threshold is it?** Excitability class sets the phase-response curve, synchronization, resonance and
gain — it is a fundamental dynamical invariant — and measuring it normally needs **intracellular**
access (f–I curves, dynamic clamp). Reading it off **extracellular spike statistics** alone (a few
drive levels, e.g. Allen Long-Square sweeps) would extract that invariant from cheap, ubiquitous
recordings. That is a tool, not a curiosity.

## Method — three features, because no single one suffices

`characterize_cell(sweeps)` fuses:

1. **ISI-shape class** — the near-threshold-weighted off-grid 3-class classifier (Type-I = the skewed
   quartic-FPT shape, SNIC/homoclinic; Type-II = the near-symmetric Hopf-rotation jitter).
2. **f–I onset (Hodgkin Class 1 vs 2)** — Type-I fires at arbitrarily low rate near rheobase
   (min/max sweep rate → small); Type-II jumps to a bounded onset frequency (min/max → O(1)).
3. **Goodness-of-fit residual** — the min `W1` to the in-family atlas; large ⇒ **out-of-family**
   (reject). This is the self-diagnosis that keeps the tool honest.

(f–I onset alone conflates Type-II with Poisson — both have bounded min/max rate — so the residual is
what tells the two apart: Type-II = high onset + **low** residual; out-of-family = high residual.)
Plus a **confound-corrected `ν̂` across drive** (proximity-to-threshold) and an **adaptation flag**.

## Validation — synthetic cohort, known ground truth (mimics the Allen design)

18 cells: Type-I (varied rheobase/σ, half with adaptation), Type-II (Hopf), out-of-family (Poisson,
bursting). Truth known per cell.

- **[N] Three-way class accuracy = 100%** (Type-I / Type-II / out-of-family), after the adaptation fix
  below. (Before the fix it was 94%: the lone error was a Type-I-**with-adaptation** cell read as
  Type-II — adaptation regularises the ISIs toward the Hopf shape.)
- **[N] Adaptation confound fixed by making f–I onset the arbiter.** Adaptation flips the ISI-shape
  vote but does **not** move the onset bifurcation, so a Type-I cell keeps its low-rate Class-1 f–I
  onset (the misclassified cell had onset 0.49, firmly Type-I; all true Type-II are 0.86–0.89, a clean
  gap). The fused rule now calls Type-II only if f–I onset **and** the shape vote agree, defaulting to
  Type-I on conflict — neutralising the confound and taking the cohort to 100%.
- **[N] Out-of-family rejection works.** Residual: in-family ≤ 0.03 (Type-II 0.009, Type-I 0.029) vs
  out-of-family **0.304** — a clean ~10× margin; Poisson and bursting cells are both rejected.
- **[N] The f–I onset recovers the Hodgkin class:** Type-II 0.87 (bounded onset rate) vs Type-I 0.47
  (reaches low rate) — the Class-1/Class-2 distinction, from spike rates alone.
- **[N] Proximity-to-threshold:** Type-I `ν̂` recovery mean `|ν̂−ν| = 0.43` across the sweep set — the
  realistic small-sample/confounded regime (near-threshold sweeps are the accurate ones, cf. the
  hardening note); the absolute scale is coarse but the **trend** across drive is monotone and usable.

## The honest catch — it does NOT yet transfer to a conductance model

*Code: `morris_lecar_check.py` (+ `morris_lecar_figure.py`), `figures/morris_lecar_check.png`.*

The cohort above is built from the same reduced first-passage normal forms as the atlas — so 100% is an
**in-family** number. The real test is a conductance model. Running the identical `characterize_cell`
on **Morris–Lecar** cells (standard Type-I/SNIC and Type-II/Hopf parameter sets, Rinzel–Ermentrout,
with noise):

- **[N] All ML cells are auto-rejected as out-of-family.** Residual: in-family normal forms `0.007`,
  reject gate `0.07`, **ML `0.091`**, garbage (Poisson) `0.108`. The conductance model sits in the
  reject band, essentially at the garbage level — the normal-form atlas does **not** faithfully
  represent conductance-model ISIs.
- **[N] The mechanism:** near rheobase ML ISIs are first-passage-like (positive skew), but ML **leaves
  that regime** at moderate drive and fires regularly with a hard period ceiling (skew flips negative,
  CV → 0.2) — outside the noise-driven atlas. The cell-level residual is pulled up by those sweeps.
- **The silver lining:** the **self-diagnosing residual makes the gap visible** — the method *refuses*
  the cells rather than confidently mislabelling them (which it could easily have done, e.g. ML Type-II
  reads as SNIC where it does attempt a call). A method that knows it's out of its depth is exactly
  what you want before touching real (conductance-governed) neurons.

**Implication, stated plainly:** the 100% in-family accuracy does **not** imply real-neuron accuracy;
external validity to biophysical neurons is the open problem.

## The kill-or-confirm test — does the distinction survive in biophysics at all?

*Code: `ml_atlas_separability.py` (+ `ml_separability_figure.py`), `figures/ml_atlas_separability.png`.*
Before rebuilding the atlas, the load-bearing question: if we TRAIN the atlas on Morris–Lecar itself,
do held-out ML cells of known class separate by **shape alone** (scale-free quantile/`W1`, rate
normalised away)?

- **[N] VERDICT: PARTIAL.** Shape-only per-sweep accuracy = **69%** all-drive (vs 89–100% for the
  idealised normal forms), rising to **78% near rheobase**. The classes *do* separate near threshold
  (near-rheobase skew 0.66 vs 1.36), but separation **collapses at higher drive** — Type-II crashes to
  25–38% in the mid-drive *regular-firing* regime (it reads as Type-I). The f–I onset still separates
  directionally (0.32 vs 0.56) but the gap narrowed from the idealised 0.47/0.87.
- **Not killed, but re-spec'd.** The Type-I/II distinction is real in conductance-model ISIs, but it
  lives in a **narrow near-threshold window** (exactly where the universal first-passage shape is
  supposed to live) and is **weak** — realistic ceiling ≈ 80%, with a biophysically-built atlas,
  near-rheobase sweeps, and shape **+** f–I onset combined. That is **population-level cell-type
  characterisation, not per-cell certainty** — a more modest but honest claim than the in-family 100%.

## Why this is (still) the useful direction

It converts the project from "predict tipping" (a crowded, under-delivering field) into "**measure a
fundamental dynamical property cheaply, with a method that refuses cells it can't characterise**." The
class call is mechanistic (it *is* the bifurcation type), not a phenomenological ISI fit, and the
reject option is the feature clinicians/experimentalists actually want from an automated readout. And
it is **validatable now**: the Allen Cell Types data has both the spike trains and intracellular
ground truth (rheobase, f–I, transgenic line), so the synthetic accuracy here can be checked against
real labels.

## The real-data path (local)

`allen_excitability.py` (run locally, `pip install allensdk`) applies the identical `characterize_cell`
to real Allen Long-Square sweeps per cell and tabulates (class, `ν̂`, residual, accept/reject) against
the cell's metadata. Same code path as the synthetic validation — the synthetic cohort is the dry run.

## Honest limits

- **Synthetic-only so far** — every number above is on models where we know the answer; the Allen run
  is the first contact with data and is not yet done here (policy: no in-sandbox data pull).
- **Adaptation costs ~1 cell** in 18 — it mimics Type-II; the joint `(ν,b)` fit flags it but the fused
  call can still flip. A near-threshold-only class vote would help (adaptation bites least near rheobase).
- **Type-II model is Stuart–Landau rotation** — real Type-II neurons have richer subthreshold
  resonance; the ISI-shape signature should survive but wants a conductance-model check (Morris–Lecar
  Type-II).
- **f–I onset needs near-rheobase sweeps** to see the Class-1 low-rate tail; sparse high-drive-only
  sampling degrades it (same near-threshold dependence as the class call).
- **Reject threshold** (`W1 = 0.07`) is calibrated on the synthetic in-family spread; real data will
  need its own calibration set.
