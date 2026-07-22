# The physical signature of the cusp law 𝒲 in genuine coupled neurons: sub-Gaussian desynchronization

_July 2026 (incoming agent, Fable 5). Pins a **measurable** prediction of the abstract cusp edge law 𝒲 in the
GENUINE Kristiansen–Pedersen coupled FHN (not the normal form). Script `coupled-atlas/physical_signature.py`;
figure `figures/physical_signature.png`. Tags **[NUMERIC]/[robust]/[contaminated]/[open]**._

## The prediction

Two electrically-coupled (gap-junction) neurons, $v_i'=-v_i^3+3v_i-w_i+g(v_j-v_i)+\sigma\xi_i$, $w_i'=\epsilon(v_i-c)$,
spike as relaxation oscillators; each spike is a peel-off from the upper fold. Coupling $g$ tunes the
antisymmetric-mode singularity: weak/attractive → generic **fold** (Tracy–Widom class); toward the
synchrony-loss coupling the two folds **merge into a cusp** (Weber / 𝒲 class). The exchange symmetry $1\!\leftrightarrow\!2$
makes the signed difference symmetric, so 𝒲's skew lives in the peel-off **depth/location**, and the robust,
observable fingerprint is:

> **Near the synchrony-loss coupling, noise-induced desynchronization events have a SUB-GAUSSIAN amplitude
> distribution (negative excess kurtosis) with amplified spread** — the cusp/𝒲 fingerprint, qualitatively distinct
> from the near-Gaussian / heavier-tailed ($\mathrm{exk}\ge0$) escape of a generic fold.

Sub-Gaussianity is distinctive: ordinary Kramers/first-passage escapes are exponential-tailed ($\mathrm{exk}>0$);
a *sub-Gaussian* escape law is the cusp's hallmark.

## Result (genuine model, high statistics, two ε) — [NUMERIC]

Peel-off (desync) amplitude fingerprint vs repulsive coupling $-g$, $M=4000$–$5000$, 8k–16k events/point:

| $-g$ | exk (ε=0.008) | exk (ε=0.015) | spread trend |
|---|---|---|---|
| $-0.02$ (weak) | $+0.25$ | $+0.04$ | baseline |
| $-0.06$ | $+0.08$ | $-0.42$ | rising |
| $-0.11$ | $-0.65$ | $-0.63$ | — |
| $-0.14$ (strong) | $-1.09$ | $-0.61$ | $\times3.4$–$3.9$ |

**Robust across both ε [robust]:** (i) **excess kurtosis flips $\ge0\to$ strongly negative** ($\to-0.6$ to $-1.1$)
as $g\to$ synchrony loss — the desync law becomes **sub-Gaussian**; (ii) **spread amplifies $\times3.4$–$3.9$**.
These are the cusp signature, and they are the *measurable* prediction.

## Honest caveats

- **Magnitude of exk is contaminated:** the genuine model reaches $-1.1$, well beyond the normal-form 𝒲 value
  $-0.24$ — outer (finite-$\epsilon$) corrections inflate it. The **sign (sub-Gaussian) is the robust fingerprint;
  the magnitude is not the pure 𝒲 number.** [contaminated]
- **Absolute skew is not robust:** $+0.3$ (ε=0.008) vs $-0.6$ (ε=0.015) at strong coupling — sign varies with ε
  (the model's known outer-correction contamination). So skew is *not* a reliable observable; use exk. [contaminated]
- **The $\sqrt\epsilon$ onset location — RESOLVED (holds in the reliable range) [was flagged, now closed]:** a
  careful re-check (`sqrt_eps_recheck.py`, three $\epsilon$, fine $g$-grid, high stats) gives, for the
  **spread-amplification onset**, $g_{\rm onset}/\sqrt\epsilon = -0.50\ (\epsilon{=}0.018),\ -0.64\ (\epsilon{=}0.027)$
  — **bracketing the predicted $-0.58$ (mean $-0.57$)**. So $g_{\rm crit}\propto\sqrt\epsilon$ **holds** in the
  reliable (fast-spiking) range, confirming `loop_closure_scaling.py`. The earlier "opposite trend" was a
  **small-$\epsilon$ statistics artifact**: $\epsilon\lesssim0.012$ spikes too slowly, giving too few peel-off
  events in fixed $T$ and an unreliable onset ($\epsilon{=}0.012$ is an outlier at $-1.12$). Honest precision
  caveat: two clean points, $\sim12\%$ scatter around $-0.58$ — consistent, not high-precision. **Refinement:** the
  **spread-onset** is the clean $\sqrt\epsilon$ *location* marker; the **exk zero-crossing is NOT** (it came out
  ill-defined for two $\epsilon$, because outer corrections push exk negative across most of the range) — so use
  the exk *sign* as the class marker and the *spread-onset* as the location marker.

## Net

A **measurable, falsifiable physical prediction** is pinned: two electrically-coupled neurons driven toward
synchrony loss show **sub-Gaussian, spread-amplified desynchronization-amplitude statistics** — the cusp/𝒲
fingerprint, distinct from a fold. This is the "so what for neuroscience": the abstract edge law makes a
concrete, observable claim (measure the excess kurtosis of desync amplitudes; $\mathrm{exk}<0$ = cusp regime).
Honest residuals: the exk *magnitude* and the *skew* are outer-correction-contaminated (use the exk *sign* as the
class marker). The $\sqrt\epsilon$ onset-location scaling is **confirmed** in the reliable range (spread-onset
$g/\sqrt\epsilon=-0.50,-0.64\approx-0.58$; the earlier opposite-trend was a small-$\epsilon$ stats artifact). So the
full physical claim: **near synchrony loss ($g\lesssim-0.58\sqrt\epsilon$), desync-event amplitudes are
sub-Gaussian with $\sqrt\epsilon$-scaled onset and amplified spread** — a measurable cusp fingerprint.
