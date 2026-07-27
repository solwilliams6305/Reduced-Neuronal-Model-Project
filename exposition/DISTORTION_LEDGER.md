# Distortion ledger
Where each planned visual LIES. Per the agreed probe protocol, every module ships with a
completed entry here, **stated on-screen in the artifact itself** — not buried in notes.

Draft; entries are hypotheses about metaphors not yet built, to be corrected by the domain
expert as each module is designed.

| Module | Where the picture lies |
|---|---|
| **Borel plane** | Singularities drawn as **points**. They are generically **branch points with cuts**, and the cut is the thing the contour must avoid. Hide it and Stokes phenomenon is unmotivated — there is no visible reason rotating the contour should matter. Compounding: Borel–Padé renders a cut as a *string of accumulating poles*, so numerically-found poles must never be drawn as if each were a genuine isolated singularity. |
| **Stokes phenomenon** | A Stokes line drawn as a sharp boundary implies the coefficient jumps discontinuously. It does not — the jump is smooth on the scale of the exponentially small term (Berry's error-function smoothing). Most textbook pictures lie here. |
| **Blow-up** | Drawing the sphere suggests points were *added* to the space. The blow-down is a diffeomorphism away from the origin; the exceptional divisor is bookkeeping, not new geometry. Secondary: three charts side by side read as three different pictures rather than three coordinate views of one object. |
| **Catastrophe ladder** | Any swallowtail that can be drawn is a **slice** of a higher-dimensional unfolding — drawing "the swallowtail" invites believing the object has been seen. A ladder drawn as a linear sequence hides the two-ladders split (swept vs genuine), which is the interesting part. |
| **Catastrophe ladder — second entry, from PREREQ_GRAPH F3** | Showing Thom A_k surfaces at all risks the conflation the papers' positioning section exists to prevent: the swept multicritical family (second-order ODE, multicritical swept turning potential) is **not** the genuine n-th-order Berry–Upstill diffraction catastrophe. Any swallowtail surface shown early must be captioned as an analogy, not the object. |
| **Trans-series** | Written as a sum it looks convergent. It is not. The notation *is* the lie, so any visual showing terms "adding up to a total" faithfully reproduces it. |
| **Uncertainty regions generally** | Precedent set in `W_make_figures.py` (fixed, commit 65d55e7): a salient best-fit marker drawn on top of a faint admissible region inverts the visual hierarchy and asserts precision the region disclaims — worse when the interval is asymmetric about the best value, as |ζ| ∈ [1.8, 2.6] with best 1.9 is. Rule: the **region** is the primary mark; the point estimate is subordinate and sits visibly where it actually falls inside its own interval. |

---

## Module γ2+γ3 (pilot) — distortion statement
Completed 2026-07-27 against `pilot_borel/storyboard.py`. **Needs expert correction.**

1. ~~**The conjugate pair is a red herring in this arc, and it is the most visually salient thing
   on panel 6.**~~ **RESOLVED 2026-07-27 — expert confirmed this was the main defect.** The pair
   is now *absent* from the pilot entirely. Rationale: the ambiguity comes **only** from the real
   singularity sitting *on* the integration contour; the pair at 1.9e^{±i50°} is off-contour and
   causes no ambiguity for real positive x. It is a fact about *this series*, not about the
   *mechanism* — so showing it here is the classic error of displaying the full object before the
   simple case has landed. Panel 6 now carries one obstruction and nothing else.
   Two seeds planted instead of an explanation: panel 5 says "they decay, but not smoothly —
   hold that thought" (previously "the wobble is the conjugate pair", which introduced a term
   with no referent), and panel 6 says "this plane has more in it — later". The pair, and the
   coefficient wobble it causes, become the payoff of a later module on the plane's finer
   structure.
2. **The model series is not our series.** Panel 2–3 use m_k = k!/A^{k+1}, a prototype with a
   *single* real singularity. Labelled on-panel, but the viewer must not conclude our v_k behave
   this way — panel 4 exists precisely because they do not.
3. ~~**Δ = e^{−A/x} is stated too cleanly on panel 8.**~~ **RESOLVED 2026-07-27.** Panel 8 now
   states Δ ∼ C·x^b·e^{−A/x} and plots both curves; their near-parallelism on log axes shows the
   prefactor does not touch the exponential character, so the concession costs the punchline
   nothing. C is named as a Stokes constant, and C and b are both marked as refinements.
4. **The cut's placement is a choice, not a fact.** Drawing it along the positive real axis from A
   to ∞ is conventional; the branch point is canonical, the cut is not. A viewer may think the
   ray itself is intrinsic.
5. **Panels 7–8 are schematic** — no lateral Borel sum is actually computed. Marked on-panel.
6. **Seven coefficients cannot establish the least-term law**, which is why panel 4 says so
   outright. This is a *strength* of the sequence (it motivates Borel) but must not be quietly
   dropped when the module is animated.

## Standing rule
If a visual is compatible with every possible outcome, it is decoration, not evidence — the
falsification probe exists to catch this. See `PROBE_PROTOCOL.md`.

---

## Standing entry — unacknowledged slicing (all α3–α4 panels)
The inner chain's real object lives in more dimensions than can be drawn. Every panel must carry
a small inset showing **where the drawn slice sits** in the full space. Without it the viewer
silently concludes the 2D picture *is* the object — the same class of error as showing a
swallowtail surface as though it were the unfolding. Suspected to be part of why the chain was
hard to hold in the first place. See `DESIGN_inner_chain.md` §5.
