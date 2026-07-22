# Validating the QIF phase-edge predictions against real neuroscience data

*Can the three phase-edge predictions of `QIF_PHASE_EDGE.md` be tested on real recordings? Yes for
two of three, with **open** data — but the test has to be designed around two confounds that the
pure QIF lacks (colored synaptic noise; spike-frequency adaptation). This note sets up the test,
shows on synthetic ground truth exactly how each confound distorts each signature (so a real result
can be read correctly), and ships a ready-to-run Allen Cell Types pull. Companions:
`qif_validation.py` (synthetic, runs anywhere), `allen_qif_validation.py` (real data, run locally),
figures `figures/qif_validation.png` and `figures/allen_qif_validation.png`. Tags **[N]/[R]/[H]**.*

---

## 0. The honest framing

The **f–I** prediction is, on its own, **not new**: the firing rate of the white-noise QIF was
solved exactly by **Brunel & Latham (2003)**. So "real Type-I neuron near rheobase ≈ noisy QIF" is
partly *expected*. What is new here is (i) reading it as a **spectral-edge / `σ^{2/3}`** law, (ii)
plotting it as an explicit data **collapse** `r·σ^{-2/3}` vs `ν=I/σ^{4/3}`, and (iii) the **two-edge
fingerprint** (SNIC/phase vs fold-of-cycles/amplitude). The *discriminating* content is the
fingerprint and the ISI-law shape, not the f–I curve alone.

Three predictions, restated as data tests:

| | observable | prediction | best data |
|---|---|---|---|
| **P1** | `r·σ^{-2/3}` vs `ν` | collapses onto `Φ(ν)=1/J(ν)`; threshold `∝σ^{2/3}` | dynamic clamp (independent `σ`) |
| **P2** | ISI CV vs drive | `1 → 0.57 → 0`; **0.57 at threshold** | Allen Long-Square (open) |
| **P3** | near-rheobase ISI shape | quartic FPT: skewed, **mode away from 0**, not exponential | Allen near-rheobase (open) |

---

## 1. Synthetic ground truth — what each confound does [N]

`qif_validation.py` runs the adaptive QIF with colored input noise
`dv=(v²+I−w)dt+ξ`, `dw=−w/τ_w dt` (`w→w+b` per spike), `ξ` either white or OU(`τ_s`), across two
noise levels and `ν∈[−0.4,2.4]`. Threshold (`ν=0`) summary:

```
   condition                 collapse r·σ^{-2/3}(0)   ISI CV(0)   skew(0)
   ideal QIF   (white)            0.166                0.57         1.83     <- the prediction
   colored     (OU τ_s=0.6)       0.14  (−16%)         0.62 (↑)     1.95
   adaptation  (sAHP τ_w=6)       0.13  (−20%)         0.51 (↓)     1.50
```

Read-outs (figure `qif_validation.png`):

* **P1.** The ideal QIF sits **on** `Φ(ν)` at *both* noise levels (collapse to ~2%) — the estimator
  works. **Colored noise** drops the whole collapse curve by a near-uniform ~15% (the Brunel–Latham
  `√τ_s` rate reduction). **Adaptation** bends it down **rate-dependently** (more suppression at
  high `ν`, where faster firing accumulates more `w`) — a *gain* distortion, not an offset.
* **P2.** Ideal QIF hits the **0.57** SNIC value at threshold. **The two confounds split it in
  opposite directions:** colored noise pushes CV(0) **up** (0.62), adaptation pulls it **down**
  (0.51, negative feedback regularizes). So **threshold CV is itself a confound diagnostic** —
  `>0.57` ⇒ colored-noise-dominated, `<0.57` ⇒ adaptation-dominated.
* **P3.** The threshold ISI is visibly **not exponential**: a clear **mode away from zero** and skew
  ≈ 1.8 (exponential: monotone, skew 2; Gaussian: skew 0). The single sharpest Poisson discriminator
  is the **short-ISI gap**: the quartic-FPT law puts essentially **zero** mass below `0.2⟨ISI⟩`
  (`frac=0.000`), versus ~18% for an exponential. Adaptation only mildly sharpens the shape (skew
  1.83→1.50); the ISI shape is the **robust** signature, least spoiled by the confounds.

> **Conclusion of the synthetic test.** All three estimators are sound, and the confounds leave
> **named, separable fingerprints**: P3 (ISI shape, esp. the short-ISI gap) is robust; P2 (threshold
> CV) is usable *and* diagnoses which confound dominates; P1 (collapse) is the most confound-sensitive
> and needs independent `σ` control.

## 2. Real data — the Allen Cell Types pull [R plan]

`allen_qif_validation.py` (run locally; `pip install allensdk`) pulls a Type-I-proxy cohort
(spiny / regular-spiking excitatory cells), reads precomputed spike times from the **Long Square**
sweeps, and emits the **same** 3-panel figure:

* **P2 — strong test.** ISI CV vs firing rate and vs `I−I_rheo`, pooled over cells; does CV cross
  `1 → 0` and cluster near **0.57** just above rheobase? Open, directly computable.
* **P3 — good test.** Pooled near-rheobase normalized-ISI density vs the exponential null and the
  white-QIF quartic-FPT reference; look for the mode away from 0 and skew ≈ 1.5–1.9.
* **P1 — partial only.** The Allen **Noise** protocol fixes the stimulus CV (0.2) and scales the
  *mean* to 0.75/1/1.5× rheobase, so `σ` co-varies with the mean — **no clean independent `σ`
  sweep**. A genuine `σ^{2/3}` collapse needs **dynamic-clamp** data (Chance–Abbott–Reyes 2002;
  Rauch 2003; Higgs–Slee–Spain 2006), where injected `σ` is controlled. Stated as a caveat in the
  script, by design.

**Caveats to honour.** (i) *Type-I/II by proxy:* spiny/regular-spiking is a stand-in for the
gold-standard phase-response-curve classification (all-positive PRC = Type-I); refine via
`CRE_LINE`/`DENDRITE`. (ii) *Adaptation:* separate the first few ISIs (transient) from the
steady-state tail, or pick weakly-adapting cells — and remember adaptation pulls CV **below** 0.57.
(iii) *Refractoriness:* an absolute refractory period also removes short ISIs, mimicking part of the
quartic-FPT gap — control by checking the gap scales with `⟨ISI⟩` (edge prediction) rather than
sitting at a fixed few-ms floor (refractory).

## 3. The decisive experiment (beyond single neurons) [H]

The headline `TW_β` lives on the **amplitude** edge, which spike times cannot see. The clean
two-edge test contrasts a **SNIC/Type-I** cell (phase edge: `σ^{2/3}`, quartic ISI) with a
**fold-of-cycles / bursting** cell, and resolves the **canard-escape amplitude** distribution
intracellularly in the burster — predicted to be `TW_β` (peel-off *level* law). That needs
voltage, not just spikes, in a bursting preparation — which is exactly why it is unexplored, and is
where the genuinely novel validation would land.

## 4. Reproduce

```
python3 qif_validation.py compute white   # synthetic ground truth (runs anywhere, ~10 s each)
python3 qif_validation.py compute ou
python3 qif_validation.py compute adapt
python3 qif_validation.py plot            # -> figures/qif_validation.png

pip install allensdk                      # real data, local only
python3 allen_qif_validation.py           # -> figures/allen_qif_validation.png
```

---

## Verdict

The phase-edge predictions **are** testable on real, partly open data **[R/N]**. **P2** (CV
crossover, `0.57` at threshold) and **P3** (quartic-FPT ISI shape, the short-ISI gap) are directly
checkable on the Allen Cell Types database; **P1** (`σ^{2/3}` collapse) needs dynamic-clamp `σ`
control and is partly a re-test of the known Brunel–Latham QIF rate. The synthetic ground truth shows
the estimators are sound and that the two confounds (colored noise, adaptation) leave **opposite,
separable** marks — so a real dataset can be read rather than guessed. The one test that would
validate the *headline* (`TW_β`) is intracellular canard-amplitude statistics in a burster, not a
spike-train measurement.
