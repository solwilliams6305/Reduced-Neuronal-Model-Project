# Architecture + sequence — gate 1 proposal
Draft 2026-07-27. Prior-art sweep still in flight; the tooling call below is from first
principles and is the part most likely to be revised by it.

---

## 1. The kernel problem, and the actual fix

The brief is "one core engine". The temptation is to pick one *language* and force everything
through it — Pyodide in the browser, or WASM everywhere, or rewrite the numerics in TS. All three
are wrong here, because they answer the wrong question.

**Split by computational cost, not by language.** This program has two disjoint classes:

| | Class A — interactive-rate | Class B — offline |
|---|---|---|
| what | ODE/SDE integration, Airy/Weber evaluation, phase portraits, contour deformation, domain colouring | Wiener-chaos transfer engine → v_k, Borel–Padé resummation, Richardson extrapolation, high-order ladders |
| cost | microseconds–milliseconds | minutes–hours |
| must run live? | **yes** — it *is* the interaction | **no** — nobody re-derives v₇ by dragging a slider |
| where it lives | small TypeScript kernel | stays in Python, exactly as now |

Class B output is a **versioned data artifact** (JSON/binary) committed to the repo. Class A is a
few hundred lines of TS. Neither needs to know the other's language.

**Drift is then controlled by a contract, not by a shared runtime.** Python emits golden test
vectors — trajectories, special-function values, coefficient ladders at fixed inputs — and a TS
test suite asserts agreement to tolerance. That is the real fix for the hand-porting problem
identified in the survey: not one language, but one *set of reference numbers* that both
implementations must reproduce. It also catches the failure mode that actually bites, which is
silent numerical divergence rather than outright breakage.

Rejected, with reasons:
- **Pyodide** — ~10 MB payload, slow cold start, and buys nothing because Class B doesn't need to
  be live and Class A is trivial to write in TS.
- **Rust/C → WASM + Python extension** — the "correct" answer in the abstract; here it imposes a
  build toolchain and slow iteration on code that is mostly a Runge–Kutta and some special
  functions. Revisit only if a Class A computation turns out to need it.
- **Precompute everything** — kills the interactivity that is the whole point.

## 2. Video: derive it from the explorable, retire Manim for this course

The protocol already fixed that the interactive is primary and video derived. Architecturally
that means video should be *rendered from the same code*, not re-authored:

**Recommended:** build in TS/canvas, capture video by driving the interactive headlessly
(Remotion, since React 19 + Vite 8 are already installed, or Playwright frame capture → ffmpeg).
One implementation of every visual.

**Honest cost of retiring Manim.** Two real losses: (a) its LaTeX typography is better than
KaTeX in the browser; (b) its animation primitives (`Transform` and friends) are well-tuned and
would be partially reimplemented. Against that: keeping Manim means a *second independent
implementation of every visual*, which is the exact drift problem in a new dimension — and the
animations this course actually needs (contour deformation past a cut, blow-up charts, domain
colouring) are domain-specific things Manim doesn't provide either. Manim is already installed
and has rendered zero videos in this project; the sunk cost is near zero.

Keep Manim available for one-off set pieces. Do not put the course's spine through it.

## 3. Proposed sequence

Forced by `PREREQ_GRAPH.md`: the geometry and asymptotics trunks are independent until they meet
at the v_k ladder. So they run as **two interleaved threads that converge**, rather than one long
chain. Thread γ can start at episode 1.

**Ep 0 — the hook.** A neuron spikes; the timing is random; here is the histogram. Two questions
planted: *what geometry produces this* (α) and *why this shape* (γ).

**Thread α — geometry**
| | |
|---|---|
| α1 | Fast–slow, the fold, relaxation oscillation. Canard shown as *phenomenon* only — striking, no theory. |
| α2 | **Where the tools die.** Fenichel, normal hyperbolicity, and its failure at the fold. Ends on the breakdown. |
| α3 | **Blow-up.** Opens on that failure. The sphere, the charts, the Riccati falling out. |
| α4 | Noise, first passage, Cole–Hopf → Airy. Airy zeros as escape levels. |
| α5 | Tracy–Widom — **and the RMT reveal as the payoff**, not a prerequisite. |

**Thread γ — asymptotics** (runs from early, in parallel)
| | |
|---|---|
| γ1 | The series that diverges. Partial sums visibly getting *worse*. |
| γ2 | Borel: the transform, the plane, and **the cut** — visible from first appearance. |
| γ3 | **The linchpin.** The contour ambiguity across the cut *is* exactly e^{−A/x} ⇒ trans-series. |
| γ4 | Stokes — how the exponentially small term switches on. **After** γ3, never before. |

**Convergence**
| | |
|---|---|
| δ1 | **The junction.** The v_k ladder. Both threads land: the divergence of γ1 is *this* divergence; the cut of γ2 sits at *this* instanton action. |
| δ2 | The ladder: turning order q, higher rungs, non-integrability, resurgence as the refuge. |
| δ3 | Coda: the two-ladders split (swept vs genuine Berry–Upstill). |

Thirteen slots is a skeleton, not a commitment — several will merge.

## 4. Pilot recommendation: γ2+γ3

Build the Borel plane and the ambiguity identity first.

- It is the **hardest node** and the linchpin beat of the whole course.
- It has the **worst existing visualisation** (in this repo: one static matplotlib scatter).
- It **stresses the architecture hardest** — complex-plane rendering, domain colouring, live
  contour deformation at 60fps. A geometric pilot would leave all of that untested.
- Its concrete series is naturally **the v_k ladder**, so it pulls in the junction node δ1 —
  meaning the pilot exercises both trunks despite sitting in one.

The safer alternative is a cheap geometric pilot (α1) to validate the pipeline before attempting
anything hard. That is a legitimate call and costs perhaps a day less. I favour the harder pilot
because the architecture question we most need answered is whether the web stack can carry the
*analytic* modules — and α1 would tell us nothing about that.

## 5. Prior-art sweep — RETURNED 2026-07-27

Verdicts (GOOD = link it, don't rebuild · PARTIAL = exists but weak · GAP = nothing decent):

| topic | verdict | note |
|---|---|---|
| Borel plane / resurgence / trans-series | **GAP** | *Largest gap of the seven and the most visually tractable.* Zero public interactive Borel-plane explorers. Best that exists is static PDF figures (Sauzin arXiv:1405.0356; Aniceto–Başar–Schiappa; CERN 2024 lectures arXiv:2511.15528). |
| Stokes phenomenon / exponential asymptotics | **GAP** | Nothing animates a Stokes crossing or the error-function smoothing. Nearest thing to a tool is Feldbrugge's Picard–Lefschetz Integrator (jfeldbrugge.github.io/Picard_Lefschetz_Integrator), a code page not an explainer. Prose: plus.maths.org "Stokes phenomenon". |
| Canards / folded nodes / GSPT / **blow-up** | **GAP** | *"Arguably the single most animation-hungry idea on this list and completely unserved."* No public interactive or video shows the blow-up sphere construction. Canonical figures: Desroches et al., SIAM Review 54 (2012); Krupa–Szmolyan 2001. |
| Painlevé transcendents | **GAP** | Fornberg–Weideman pole-field images are striking and camera-ready; no interactive, no video. |
| Diffraction catastrophes (wave-dressed) | **GAP** | — |
| Catastrophe unfoldings / caustics (ray-level) | **PARTIAL** | Dujardin's *Catastrophe Teacher* applets are ingenious but **Java, i.e. bit-rotted**. AMS Feature Column (Tony Phillips); Berry's own gallery (~50 static images). Elementary geometry is covered several times over. |
| Tracy–Widom | **PARTIAL** | Quanta's "At the Far Ends of a New Universal Law" covers *why it matters* well; Tao's lecture series is chalk-talk. **Nothing visualises the edge-scaling mechanism** (semicircle edge → Airy kernel → TW), so "why Airy?" is unexplained everywhere. |
| **Airy / Pearcey functions** | **GOOD — do not rebuild** | **DLMF §36.3** (dlmf.nist.gov/36.3), authored by Berry and Howls: 21 figures covering Pearcey, swallowtail across 4 parameter slices, elliptic/hyperbolic umbilics, modulus *and* phase, with rotatable interactive 3D. Reference-grade. What it lacks is narrative, not pixels. |

**Consequences for the plan.**
1. **Pilot choice γ2+γ3 is confirmed** — the biggest gap coincides with the linchpin beat and with
   the most tractable animation. Nothing to borrow; everything to invent.
2. **δ2 (the ladder) must link DLMF §36.3 rather than replot special functions.** Our contribution
   there is narrative and the swept-vs-genuine distinction, not surfaces. This also pairs with the
   ledger entry warning against showing Thom surfaces as if they were the object.
3. **α5 (Tracy–Widom) has a specific unserved target**: the *mechanism*, not the significance.
   Quanta already does significance better than we would.
4. **Blow-up (α3) is the second-biggest opportunity** — a genuine void, and the one place 3D is
   justified rather than decorative.

**Tooling question — answered, and better than proposed above.** Remotion's dual output is a
first-class supported pattern, not a hack: `@remotion/player` embeds *the same React composition*
interactively in a page, and `renderMedia` renders that identical composition to MP4
(remotion.dev/docs/reusability). Existence proof: **GitHub Unwrapped** (githubunwrapped.com/about,
open source) drives an interactive site and per-user rendered video from one codebase at 10k+
scale. ~39k stars, ~900k installs/month as of 2026. So §2 above should be read as *Remotion Player
+ renderMedia*, not headless capture — one composition, two outputs, no second implementation.

**One honest caveat:** no public explorable-explanation project was found that pairs a TypeScript
kernel with Python-generated reference data. Cross-language golden-file validation is ordinary
engineering practice and nothing about it is technically risky, but we would be inventing the
workflow rather than following one.
