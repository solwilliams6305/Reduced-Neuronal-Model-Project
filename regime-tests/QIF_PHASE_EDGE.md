# Direction #1: the phase edge **is** the noisy Type-I (QIF / θ) neuron

*The phase channel of the noisy folded cycle is not just analogous to a spiking neuron — it is
**literally** the canonical Type-I neuron model. The quadratic-integrate-and-fire (QIF) neuron, and
its smooth twin the θ-neuron, are the **SNIC normal form**; adding white current noise makes them
**exactly** the off-critical noisy saddle-node studied as the phase edge. So every phase-edge result
becomes a concrete, falsifiable statement about near-threshold spiking: the inter-spike-interval
(ISI) law, the noise-smoothed f–I curve, and the spike-train regularity. Companion code:
`qif_firing_statistics.py`, `offcritical_phase_edge.py`, figure `figures/qif_firing_statistics.png`.
Tags **[R]** rigorous/exact, **[N]** numerical, **[H]** heuristic.*

---

## 0. The dictionary [R]

The noisy QIF neuron and the off-critical noisy saddle-node are the **same SDE**:

```
   noisy QIF / θ-neuron            off-critical inner Riccati (phase edge)
   dv = (v² + I) dt + σ dW   <-->  dR = (μ + R²) dT + σ dB        (Y = −μ)
        v   spike variable              R   inner Riccati variable
        I   input current               μ   detuning  (= −Y, distance past the saddle-node)
        σ   current-noise amplitude     σ   inner-noise amplitude
   reset v→v_reset at v→+∞         <-->  rotation R: −∞ → +∞  (one first-passage = one spike)
```

The QIF is the normal form of the **saddle-node on an invariant circle (SNIC)** bifurcation — the
defining signature of **Type-I excitability**. Its smooth conjugate is the Ermentrout–Kopell
θ-neuron, `θ' = (1−cosθ) + (1+cosθ)I`, mapped by `v = tan(θ/2)`. So:

> **A noisy Type-I neuron near its spiking threshold is a noisy saddle-node, and one ISI is one
> first-passage of `dR = (μ+R²)dT + σdB`.** Everything proved for the phase edge transfers verbatim.

Rescaling `v = σ^{2/3} ρ`, `t = σ^{-2/3} s` removes `σ` and leaves a **one-parameter family** in the
single rescaled current
```
        dρ = (ν + ρ²) ds + dW ,         ν = I / σ^{4/3}            (rescaled detuning).
```
`ν` is the only knob: it interpolates sub-threshold (`ν<0`), threshold (`ν=0`), and supra-threshold
(`ν>0`) firing.

---

## 1. Three observables, three predictions

### 1a. The f–I curve is a **noise-smoothed Type-I onset**, `r(I,σ) = σ^{2/3} Φ(ν)` [R + N]

The firing rate is the inverse mean ISI, `r = 1/⟨ISI⟩ = σ^{2/3} Φ(ν)`, with the **universal scaling
function**
```
        Φ(ν) = 1 / J(ν) ,        J(ν) = canonical mean first-passage time of  dρ=(ν+ρ²)ds+dW .
```
Three limits, all confirmed numerically (figure, panels A–B):

* **supra-threshold** `ν ≫ 1`: `Φ(ν) → √ν/π`, i.e. `r → √I/π` — the **deterministic Type-I
  square-root onset** recovered as noise becomes irrelevant.
* **threshold** `ν = 0`: `r = σ^{2/3} / J(0)` — the **universal `σ^{2/3}` scaling** of the SNIC edge.
  No deterministic onset survives; the rate is set entirely by noise.
* **sub-threshold** `ν ≪ 0`: Kramers/Arrhenius **noise-induced firing**,
  `Φ(ν) → (√|ν|/π)·exp(−(8/3)|ν|^{3/2})`. The barrier `ΔU = (4/3)|I|^{3/2}` gives the escape weight
  `exp(−2ΔU/σ²) = exp(−(8/3)|ν|^{3/2})` — the same `|·|^{3/2}` Arrhenius exponent that governs the
  original `σ_crit` story, here as the **low end of the f–I curve**.

The noise **rounds the sharp Type-I knee** at `I=0` into a smooth onset of width `ΔI ∼ σ^{4/3}` — the
experimentally familiar smearing of the rheobase.

### 1b. The ISI law is the **quartic FPT law — and is NOT Tracy–Widom** [R]

This is the sharp structural statement, and the contrast that makes the two-edge picture worth
stating. The ISI is the first-passage time of the noisy saddle-node; its Laplace transform solves an
**anharmonic (quartic) Schrödinger operator**
```
        −χ'' + (ρ⁴ + 2νρ² + 2ρ + ν²) χ = 2λ χ .
```
This operator is **non-integrable** — there is no Painlevé/Tracy–Widom closed form for the ISI
distribution, unlike the **amplitude** edge (the canard peel-off **level**), whose operator is the
stochastic Airy operator and whose law **is** `TW_β`. The phase-edge ISI is a genuinely different,
skewed law: **not exponential** (so not a Poisson/renewal escape) and **not Gaussian** (so not a
diffusive jitter), but the quartic FPT law interpolating between them.

What *does* have a closed form is the **mean** ISI constant at threshold. The canonical mean FPT is
```
        J(0) = (√π/3) · 12^{1/6} · Γ(1/6) ≈ 4.976           [closed form, unit-diffusion convention]
```
(a `Γ(1/6)` constant — the **phase** analogue of the amplitude edge's `−2.338` first Airy zero). So
the **threshold firing rate has a closed-form prefactor**: `r(0,σ) = σ^{2/3}/J(0)`.¹

### 1c. The spike train **regularity crossover** in the ISI CV [N]

The ISI coefficient of variation `CV = std(ISI)/⟨ISI⟩` is a clean, dimensionless regularity readout
that crosses over with the *single* parameter `ν` (figure, panel C):
```
        ν ≲ 0   CV → 1     Poisson-like, sub-threshold (Kramers escapes are nearly memoryless)
        ν = 0   CV ≈ 0.57  the SNIC critical value (quartic FPT)
        ν ≫ 0   CV → 0     regular firing, supra-threshold (deterministic clock + small jitter)
```

---

## 2. Numerics — the confirmation [N]

Direct simulation of `dv = (v²+I)dt + σ dW` with reset, two noise levels `σ∈{0.5,0.3}`, four
currents per level (`qif_firing_statistics.py`):

```
 sigma     I       nu    r(QIF)   sigma^{2/3}·Phi(nu)   CV       r·sigma^{-2/3}
  0.50  -0.159   -0.40   0.0584      0.0556            0.74         0.093
  0.50   0.000    0.00   0.1032      0.1001            0.58         0.164
  0.50   0.278    0.70   0.1839      0.1761            0.37         0.292
  0.50   0.635    1.60   0.2666      0.2532            0.24         0.423
  0.30  -0.080   -0.40   0.0427      0.0395            0.72         0.096
  0.30   0.000    0.00   0.0742      0.0712            0.57         0.166
  0.30   0.141    0.70   0.1291      0.1252            0.37         0.289
  0.30   0.321    1.60   0.1872      0.1801            0.24         0.419
```

Two facts:

1. **Universal collapse.** The last column `r·σ^{-2/3}` is, at each `ν`, the **same for both noise
   levels** (`0.094, 0.165, 0.290, 0.42` to ~2%) — the `σ^{2/3}` SNIC scaling and the universal
   `Φ(ν)` curve, verified by data collapse across a 1.67× change in `σ`.
2. **The theory line is right.** `r(QIF)` tracks `σ^{2/3}Φ(ν)` to ~3–8% across the whole crossover
   (the small systematic over-shoot is the finite reset/threshold `±14` and finite `dt` — both push
   the simulated rate slightly above the idealised infinite-domain MFPT).

The CV column reproduces the predicted monotone crossover `0.73 → 0.57 → 0.37 → 0.24`.

---

## 3. Why it matters: the two-edge **neural fingerprint** [R framing / H prediction]

The noisy folded cycle has **two distinct edges**, and they map onto **two distinct kinds of
neuron** with **opposite noise signatures**:

| | **phase / time edge** | **amplitude / level edge** |
|---|---|---|
| neuron | **SNIC / Type-I** (QIF, θ) | **fold-of-cycles / Type-II burster** |
| what's noisy | *when* it spikes (ISI) | *how big* the canard escape is (peel-off level) |
| scaling | `σ^{2/3}` (threshold rate) | `σ_* ∼ √ε₂` (peel-off onset) |
| law | **quartic FPT** (skewed; **no TW**) | **Tracy–Widom `TW_β`** |
| constant | `J = Γ(1/6)`-closed-form mean ISI | `−2.338` (first Airy zero) |

This is a **classifier**: an experimenter who measures (i) the f–I scaling exponent near rheobase
(`σ^{2/3}` ⇒ SNIC phase edge) and (ii) the ISI CV crossover gets a noise-only discriminant between
Type-I (SNIC) and the fold-of-cycles/bursting class — without needing to resolve the deterministic
bifurcation structure. The headline `TW_β` belongs to the **amplitude** edge; the **phase** edge —
the ordinary Type-I spiking neuron — carries the **quartic** ISI law instead. Stating both, and that
they are *different*, is the contribution.

---

## 4. Reproduce

```
python3 qif_firing_statistics.py     # -> figures/qif_firing_statistics.png  (3 panels)
```
Panel A: noise-smoothed Type-I f–I curves for `σ∈{0.5,0.3}` with the deterministic `√I/π` asymptote
and the QIF points. Panel B: the universal collapse `r·σ^{-2/3}` vs `ν=I/σ^{4/3}` onto
`Φ(ν)=1/J(ν)`, with the `√ν/π` (supra-threshold) and `(√|ν|/π)e^{−(8/3)|ν|^{3/2}}` (Arrhenius)
asymptotes. Panel C: the ISI-CV regularity crossover. `Φ(ν)=1/J_of_nu(ν)` is imported from
`offcritical_phase_edge.py` (mean-FPT quadrature).

---

## Final status (Direction #1)

The phase edge of the noisy folded cycle **is** the near-threshold statistics of a noisy Type-I
(QIF / θ) neuron — an exact identification, not an analogy **[R]**. It delivers three falsifiable
predictions: a noise-smoothed `r=σ^{2/3}Φ(ν)` f–I curve with the deterministic `√I/π` and Arrhenius
`e^{−(8/3)|ν|^{3/2}}` limits **[R+N]**; an ISI law that is the **quartic** first-passage law — skewed,
**provably not Tracy–Widom**, with a `Γ(1/6)` closed-form mean-ISI constant **[R]**; and a
`CV: 1→0.57→0` spike-train regularity crossover in the single parameter `ν=I/σ^{4/3}` **[N]**. The
`σ^{2/3}` collapse and the f–I match are confirmed by data collapse across two noise levels to ~2%.
Together with the amplitude edge's `TW_β` peel-off law, this gives a **two-edge neural fingerprint**
that distinguishes SNIC (Type-I) from fold-of-cycles (bursting) neurons by their noise statistics
alone.

---

¹ Convention note: `qif_firing_statistics.py` imports `J_of_nu` from `offcritical_phase_edge.py`,
whose `ρ`-equation uses diffusion `D=1/2`; this gives `J_of_nu(0)=6.294 = 2^{1/3}·4.976`, i.e. the
closed-form constant above rescaled by the `2^{1/3}` that the `D=1/2` vs `D=1` choice carries through
the `σ^{2/3}` rescaling. The figure's theory column uses `J_of_nu(0)` consistently, so
`r(0,σ)=σ^{2/3}/J_of_nu(0)` matches the simulation; the `Γ(1/6)` statement is the canonical
(unit-diffusion) form of the same number.
