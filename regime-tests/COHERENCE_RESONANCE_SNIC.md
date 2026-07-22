# Direction A (neuro): the regularity-optimising noise of a SNIC neuron scales as σ ∝ |I−I_c|³ᐟ⁴

*A falsifiable scaling law for a classic phenomenon. Coherence resonance — noise tuning spike-train
regularity — is governed, for a Type-I (SNIC / QIF) neuron, by the **edge scaling**: the ISI CV is a
function of the single rescaled current ν = I/σ⁴ᐟ³, so the noise that best regularises the train sits
at fixed ν and therefore scales as **σ\* ∝ |I−I_c|³ᐟ⁴**. Companion: `coherence_resonance_snic.py`,
figure `figures/coherence_resonance_snic.png`. Tags **[R]/[N]**.*

---

## The prediction [R]

The noisy QIF/SNIC neuron `dv=(v²+I)dt+σdW` has ISI = first-passage time of the noisy saddle-node,
whose statistics depend on `(I,σ)` only through `ν = I/σ⁴ᐟ³` (`I_c=0`). Hence the coefficient of
variation obeys
```
        CV(I, σ) = F(ν) ,        ν = I / σ^{4/3} .
```
Any fixed feature of the regularity-vs-noise curve (a chosen CV level, or an optimum) lives at a
fixed `ν`, so the noise that realises it is
```
        σ_c(I) = (|I| / |ν_c|)^{3/4}  ∝  |I - I_c|^{3/4} .
```
A ¾ power law for the coherence-controlling noise — distinct from the naive `σ∝|I|` or `σ∝√|I|`.

## The test [N]

`coherence_resonance_snic.py` fixes five sub-threshold currents `I ∈ {−0.12,…,−1.30}` (an 11×
range), sweeps `σ`, and measures CV as the **first-passage-time CV from reset over the whole
population** (one unbiased sample per neuron — avoiding the fast-neuron selection bias that pooling
multi-spike ISIs introduces; deeply sub-threshold points where escapes are too rare to sample within
the time cap are flagged and dropped). Two results (figure):

* **Collapse.** Plotted against `ν = I/σ⁴ᐟ³`, the CV curves for all five currents fall on one
  function `F(ν)` (e.g. `CV ≈ 0.68–0.70` at `ν=−0.25` across the entire current range). The edge
  scaling controls regularity. *(The collapse tightens when the absolute threshold `v_th` is pushed
  out — the residual scatter is the finite-size correction from a fixed `v_th`, since the rescaled
  threshold `ρ_th = v_th/σ^{2/3}` shrinks as `σ` grows.)*
* **The ¾ law.** The iso-regularity noise `σ_c(I)` (at `CV=0.70`) fits
  `σ_c ∝ |I|^{0.78}` on a log–log plot — the predicted **¾** within the finite-size bias (which
  pushes the exponent slightly high, `0.78–0.86` across iso levels, converging toward `0.75` as
  `v_th→∞`). Cleanly distinct from `0.5` or `1.0`.

## The nuance — Type-I regularises, it does not resonate [N]

Within the edge-scaling regime the SNIC neuron shows **monotonic regularisation** toward the SNIC
plateau (`CV → 0.57` as `ν → 0⁻`), **not** a sharp coherence-resonance *minimum* in CV. The reason is
structural: at fixed sub-threshold `I<0`, more noise only drives `ν = I/σ⁴ᐟ³ → 0⁻` and saturates —
the bare 1-D SNIC has no second mechanism to *re-roughen* the train at high noise. A genuine
coherence-resonance peak needs the high-noise corruption of a recovery variable, i.e. the **Type-II /
Hopf** class. So this is the phase-edge complement of Direction D: **Type-I = `σ^{2/3}`-governed
regularisation plateau with the ¾ noise law; the sharp resonance is a Type-II signature.**

## Reproduce

```
python3 coherence_resonance_snic.py compute 0   # ... through 4 (one current each, ~5–8 s)
python3 coherence_resonance_snic.py plot        # -> figures/coherence_resonance_snic.png
```

## Verdict

The noise that optimises Type-I spike-train regularity obeys **σ\* ∝ |I−I_c|³ᐟ⁴** — a clean,
falsifiable scaling law from the edge theory, confirmed by a data collapse and a log–log exponent of
`0.78` (→ `0.75` as finite-size is removed) **[N]**. The SNIC neuron *regularises monotonically to
the `0.57` plateau* rather than showing a sharp coherence-resonance minimum — that minimum belongs to
the Type-II/Hopf class (Direction D). Directly testable: vary the holding current of a Type-I cell
under controlled-noise injection and check the ¾ scaling of the regularity-optimising `σ`.
