# Taxonomy sweep — the early-warning observable splits along the two-edge line

*Phase 3b-cont. Code: `taxonomy_sweep.py` (+ `taxonomy_figure.py`). Figure:
`figures/taxonomy_sweep.png`. Tags: **[R]** proved / **[N]** numerical / **[H]** heuristic.*

```
python3 taxonomy_sweep.py sims     # integrate + cache all four classes (~4 s)
python3 taxonomy_sweep.py analyze  # the per-class early-warning map
python3 taxonomy_sweep.py fig      # -> figures/taxonomy_sweep.png
```

## The question

Phase 3 showed the edge-law `ν̂→0` warning needs an **adiabatic, period-diverging** approach. So *which*
bifurcations qualify? We drift the bifurcation parameter slowly through **four** canonical onset/offset
bifurcations and, per class, measure: does the inter-event period diverge? does windowed `ν̂` (the
timing inversion) warn? does the oscillation amplitude collapse — gradually (a precursor) or abruptly?

Models (each with a slow linear drift crossing its bifurcation):

- **SNIC** — saddle-node/QIF normal form `dv=(v²+I)dt+σdB` (the model the atlas is built on), `ν` drifts
  through 0;
- **Saddle-homoclinic** — the project's reduced log-FPT model (`ISI = T_ret + λ⁻¹ln(δ/|μ+ση|)`),
  `μ→0` (genuine dynamical SH needs Bogdanov–Takens path-following — out of scope; the reduced model
  carries the defining log-divergence);
- **Supercritical Hopf** — Stuart–Landau `ż=(b+iω)z−|z|²z`, `b` through 0 (amplitude `~√b`);
- **Fold-of-cycles** — `ṙ=br+r³−r⁵`, `b` through the fold at `−0.25` (amplitude jumps to 0).

## The map (one run each)

| bifurcation | edge | period change | timing `ν̂` warns | amplitude precursor | amplitude→0 |
|---|---|---|---|---|---|
| SNIC | **phase** | **2.2× (diverges)** | **yes** (`ν̂` min `−0.11`) | n/a | no |
| Saddle-homoclinic | **phase** | **2.2× (diverges)** | **yes** (via period; `ν̂`→0.64) | n/a | no |
| Supercrit Hopf | **amplitude** | 1.0× (finite) | no (`ν̂` blind) | **yes** (→0.34 before crossing) | yes |
| Fold-of-cycles | **amplitude** | 1.0× (finite) | no (`ν̂` blind) | **no** (held at 0.85, then jumps) | yes |

**[N] Phase edges (SNIC, saddle-homoclinic).** The inter-event **period diverges** as the parameter
approaches threshold (2.2× over the windows before crossing), so the **timing** channel warns. For the
SNIC the calibrated inversion drops `ν̂→0` (min `−0.11`, crossing the `0.4` alarm) because the SNIC ISIs
*are* the atlas's first-passage law. For the saddle-homoclinic the warning comes through the period
divergence directly; its ISIs are log-tailed (a different shape from the SNIC atlas, so `ν̂` reads
`0.64` — one would use the homoclinic-class distance, which the classifier already provides).

**[N] Amplitude edges (supercritical Hopf, fold-of-cycles).** The **period stays finite** (1.0×), so the
timing inversion is **blind** — correctly. The early-warning observable here is the **amplitude**, which
collapses to 0. And the two amplitude edges differ in a way that matters operationally: the
**Hopf** amplitude declines **gradually** (to 0.34 of baseline *before* the crossing) — a genuine
**precursor**; the **fold-of-cycles** holds its amplitude (0.85) and then **jumps** — **no precursor in
either channel**, the genuinely catastrophic, warning-free tipping.

## Why this matters — it *is* the paper's two-channel split

**[H]** The early-warning observable partitions along exactly the line of the project's headline
two-edge theory:

- **phase edge** = Type-I, `σ^{2/3}`, the quartic-FPT timing channel → **`ν̂→0` / period divergence**;
- **amplitude edge** = Type-II, `σ^{1/2}`, the canard peel-off / Tracy–Widom level channel →
  **amplitude collapse**.

So the inversion's **class call is not a detail — it selects the early-warning observable.** Identify a
Type-I (SNIC/homoclinic) approach and watch the timing (`ν̂→0`); identify a Type-II (Hopf/fold) approach
and watch the amplitude; and recognise a **fold-of-cycles** as *un-warnable* from precursors at all. This
turns the steady-state two-edge universality result into an **operational early-warning prescription**,
and it situates the earlier phases: Phase 3a is the phase-edge (adiabatic SNIC) success; the Epileptor
(Phase 3b) offset is *neither* period-diverging *nor* adiabatic, so it sits outside both channels.

## Honest limits

- **Mixed model fidelity.** SNIC and saddle-homoclinic use the noise-driven normal-form / reduced-FPT
  models (the regime where the universal ISI shape lives); Hopf and fold use integrated Stuart–Landau /
  amplitude oscillators. A fully dynamical SH oscillator (BT path-following) and a noisy fold-of-cycles
  with an amplitude-edge `ν̂`-analogue inversion are the natural next steps.
- **Threshold-based criteria.** "Warns" uses simple cutoffs (`ν̂<0.4`, period `>2×`, amplitude ratios);
  a per-class matched-FPR benchmark like Phase 3a would quantify lead/false-alarm trade-offs.
- **One drift rate / noise level per class.** The adiabaticity boundary (how slow the drift must be,
  from Phase 3b) is held fixed here, not swept — combining the two into a (drift-rate × class) map is
  the clean follow-up.
