# Survey brief — exposition/visualisation project
Compiled 2026-07-27. Inputs: two Explore agents + direct verification.

## User brief (from clarifying questions)
- **Format:** ONE core engine → both interactive web explorables AND rendered narrated video.
- **Scope:** a *sequenced course* (each module assumes the last, mirroring the argument chain).
- **Audience:** grad / mathematically literate (knows ODEs + complex analysis; has never seen
  resurgence or catastrophe theory). Kirsanov's band.
- **Hardest to convey (all four selected):** resurgence/trans-series, Stokes phenomenon,
  the catastrophe ladder, canards/blow-up.
- Style targets: 3Blue1Brown, Artem Kirsanov.
- Explicitly: do NOT feel obliged to build on existing assets.

## The headline story (what the course must ultimately deliver)
> Noise-induced escape at a slow-fast catastrophe singularity follows a universal pipeline whose
> output law depends only on the singularity *type* (fold, cusp, swallowtail, …), not the model.
> That law is always a resurgent trans-series; for the fold it is exactly Tracy–Widom.

Argument chain (each link is a candidate module):
1. FHN neuron → slow-fast dynamics with a fold. (`PROJECT_CONTEXT.md`)
2. Geometric blow-up (Krupa–Szmolyan) desingularises the fold → canonical form. (`THEORY_ROADMAP.md` L1.5–1.7)
3. Add noise → first-passage problem, Freidlin–Wentzell action.
4. Cole–Hopf: inner Riccati → Schrödinger with white-noise potential. (L3.1)
5. That operator IS the stochastic Airy operator. (L3.4)
6. Exit law = Tracy–Widom, β = 4/η². (L3.5, `NoisyFoldedCycle_paper.tex`)
7. Generalise → catastrophe ladder indexed by singularity order q. (`CATASTROPHE_LADDER.md`)
8. Higher rungs non-integrable → only closed form is a resurgent trans-series.
   (`W_ResurgentTransseries_paper.tex`, `SwallowtailTransseries_paper.tex`)
9. Wiener-chaos transfer engine computes high-order weak-noise coefficients exactly.
10. Borel–Padé resummation reproduces ground truth (~1% at cusp/β=2).

Cross-cutting themes: blow-up architecture; the ladder; perturbation→resurgence;
identification with RMT; integrable vs non-integrable (fold rigorous → cusp conjectural →
swallowtail provably non-integrable ⇒ "resurgence is the refuge").

## Concepts, by visual character (condensed from the inventory)
**Intrinsically geometric (animate these):** phase plane/nullclines; bifurcation; slow-fast flow &
relaxation oscillation; fold of equilibria & of limit cycles; canards; canard explosion; blow-up
charts K1/K2/K3 + matching; noise tube around a trajectory; Riccati blow-up R→∞; Borel plane
singularities & Stokes-line crossing; Airy zeros; the A_k unfolding/caustics.

**Intrinsically analytic (need invented visual metaphor — the hard part):** asymptotic/factorial
divergence; Borel transform & summation; trans-series; resurgence; Stokes phenomenon; Padé;
Darboux large-order fitting; stochastic Airy operator; Tracy–Widom; Painlevé II/IV;
Lax pair/isomonodromy; Wiener chaos; Cole–Hopf.

Difficulty tiers: T1 undergrad (phase plane, bifurcation, FHN, Euler–Maruyama, Monte Carlo);
T2 grad (GSPT/Fenichel, SDE, Fokker–Planck, FPT, quasipotential, Kramers, Airy/Weber/Bessel,
Padé); T3 specialist (blow-up, stochastic Airy operator, Borel/resurgence/trans-series, Stokes,
Painlevé, isomonodromy, Wiener-chaos transfer, instanton action).

Verification status matters for honesty in exposition: TW-at-fold **proved**; Bessel-1/(q+2)
skeleton + instanton action s^(2q+1)/[2(2q+1)] **verified q=1..5**; Borel phases **numeric**
(swallowtail θ pinned 33–42°, best 39.4; cusp 50°±2); cusp=Painlevé IV **conjecture**;
swallowtail **proved non-integrable**; Painlevé reduction q≥3 **open**.

## Existing assets — verified inventory
**Standalone HTML + vanilla JS/2D-canvas (no build step), all polished & working:**
- `two_channels_animation.html` — Morris–Lecar, two competing noise channels (phase vs
  amplitude), 5 canvases, ~1400 lines. Light theme.
- `morris_lecar_animator.html` — ML Type-II phase plane + voltage trace, Euler–Maruyama live,
  slow-passage protocol. Dark theme.
- `shooting_and_spectral_function.html` — shooting for node location vs noise; 4 panels;
  demonstrates the *bias* problem that G solves; burst sampling + histogram. Light theme.
- `coupled-atlas/coupled_explorer.html` — gap-junction coupled FitzHugh–Rinzel pair; 6 sliders,
  4 presets, synchrony metric. Dark theme.

**Build-required:** `fhn/` Vite 8 + React 19 app; `src/App.jsx` (20 KB) mirrors `kernel.py`
geometry. Dev-ready, not built. Also root `fhn_animator.jsx` (20 KB, standalone copy).

**Manim:** `archive/animate_fhn.py` — one archived scene (FHN excitability + stochastic loop,
4 annotated phases). Manim CE **0.20.1 confirmed installed and importable**. `media/` holds only
Pango text→SVG cache + one Tex render; **no rendered videos exist**.

**Static:** ~130 files in `figures/` (matplotlib, mixed exploratory/publication quality);
`W_make_figures.py` → publication PDFs (weak-noise coefficient ladder w/ sign-coloured stems +
resurgent envelope; Borel plane w/ conjugate pair + instanton root). Paper figs
`swtl_fig_{borel,ladder,resum}.pdf`, `cap_fig_{ladders,skeleton}.pdf`.

**CORRECTION to audit agent:** it reported `G_linear_vs_wiggle_contrast.svg`,
`G_projection_picture.svg`, `G_wiggle_to_tracywidom.svg` as empty placeholders (0 lines).
FALSE — they are 15 KB / 3.5 KB / 14 KB of real hand-authored SVG (46 rects, 8 texts,
8 polylines, 8 circles + embedded `<style>` in the first; similar in others), just minified to a
single line. These are a genuine existing diagram vocabulary. `folded_cycle_hazard_separatrix.svg`
(16 KB, 5 paths) is likewise real.

## Installed toolchain (verified)
Python: Manim CE 0.20.1, matplotlib 3.10.0, numpy 1.26.4, scipy 1.17.1, plotly 6.3.0,
ManimPango 0.6.1 — on **conda base**, no project venv.
Node v24.3.0, npm 11.4.2. React 19.2.6 + Vite 8.0.12 in `fhn/`. Root `package.json` has only
`lucide-react`. **Not present:** d3, three.js, p5, plotly.js.

## De-facto current architecture — and its core flaw
```
Python (kernel.py, phase*.py, chaos_transfer.py)  ← source of truth
  ├→ W_make_figures.py (matplotlib)   → PDF/PNG
  ├→ archive/animate_fhn.py (Manim)   → (no video rendered)
  └→ *.html / *.jsx  ← HAND-PORTED physics, no automatic sync
```
**The drift risk is the thing to fix.** JS re-implementations of the dynamics were hand-written
from the Python; nothing keeps them in step. This is exactly what the user's "single core engine"
requirement should solve — and is the strongest architectural argument in the whole survey.

Existing visual style: two palettes (warm light `#fbfaf8`/`#20242b`; deep-navy dark
`#0b1020`/`#e2e8f0`), system sans UI, monospace numerics, 8–12 px radii, colour-coded roles
(phase=green, amplitude=orange, error=red, SNIC=blue). Consistent enough to keep as a base.

## Gaps with no visual at all
No 3D/foliation view of the slow-fast reduction; no parameter atlas/heatmap over
(g, τ, Δc, σ); no bifurcation landscape that transitions smoothly through regimes; and —
most importantly for this project — **nothing whatsoever for the four hardest concepts**:
resurgence/trans-series, Stokes phenomenon, the A_k ladder, blow-up desingularisation.
The Borel plane exists only as a static matplotlib scatter in `W_make_figures.py`.

## Open question the deep research must settle
Realistic mechanism for "one kernel, two outputs": rewrite kernel in TS and call from Python;
compile Rust/C → WASM + Python extension; Pyodide in browser; or precompute to JSON/binary and
have both sides merely *render*. Plus: Manim CE vs Motion Canvas/Remotion for the video side,
given Manim is already installed but the web side is where interactivity lives.
