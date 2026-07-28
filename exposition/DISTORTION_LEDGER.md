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
| **Catastrophe ladder** | Any swallowtail that can be drawn is a **slice** of a higher-dimensional unfolding — drawing "the swallowtail" invites believing the object has been seen. A ladder drawn as a linear sequence hides the **three**-way split, which is the interesting part. Note this was recorded as a *two*-ladder split (swept vs genuine); the swallowtail paper's positioning adds a third — the random-matrix multicritical edge / "higher-order Tracy–Widom", which collides on both words and is **integrable by construction** where ours is provably not. That is the most dangerous of the three conflations and the one δ3 exists to prevent. |
| **Coda δ3 — three columns** | Placing the three families side by side implies a common construction. Only the **word** is common: families 1 and 2 are ODEs indexed by turning order / equation order, family 3 is a limit law of an ensemble. The third column therefore carries no coefficient ladder beside the others, and the asymmetry is deliberate — drawing one would manufacture a parallel that does not exist. |
| **Coda δ3 — the skew curves** | The separation is the claim, and the values are **measured, not derived**: Monte-Carlo (2.8×10⁶ escapes per rung) for the swept family, numerical integration for Berry–Upstill. Sampling error is ≈ the last digit shown — far below the separation, but not zero. Three points per curve is also too few to assert a *trend*; the honest claim is a shared first rung and opposite second differences. |
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

---

## Module α4a (R↔u twin panel) — distortion statement
Completed 2026-07-27 against `pilot_borel/twin_panel.py`. Rows 1–2 are exact; row 3 is honest
but weaker than it looks.

1. **Row 3 does NOT exhibit Tracy–Widom, and must never be captioned as though it did.**
   At η=0.55 on this grid we are nowhere near the scaling regime. Measured skew +0.11 against a
   standard error of 0.12 over 400 runs — *not significantly different from zero*, and well short
   of TW's positive skew. The panel now says "not the scaling limit — this is not TW" on its face.
   What row 3 legitimately shows is only that **the first zero fluctuates**; that its limit law is
   TW is the content of the theorem, not of this picture.
2. **Every Airy zero inside the window must be marked.** An early draft computed four zeros while
   five lay in range, leaving one asymptote in row 2 with no dot above it — which silently
   falsifies the panel's entire claim. Guarded now by filtering `ai_zeros(12)` to the window
   rather than taking a fixed count. *Any change to `Y_LO` must preserve this.*
3. **The rows must never drift out of vertical register.** The shared x-axis *is* the argument.
   If the module is ever animated or re-laid-out, alignment is load-bearing, not cosmetic.
4. **R is clipped at |R|>12 for legibility**, so the asymptotes are drawn as finite excursions.
   The clipping is a rendering choice; the divergence is real and unbounded.
5. ~~**Only the recessive (Ai) solution is shown.**~~ **RESOLVED 2026-07-27 — expert requested a
   caption.** Row 1 now states on its face that u = Ai is the solution which *decays* as Y→+∞,
   that the escape problem is what selects it, and that a different boundary condition (Bi, or a
   mixture) has different zeros — "so these are not *the* zeros, they are this condition's zeros."
   The choice is no longer invisible.
6. ~~**Row 3 broke the one-invariant discipline.**~~ **RESOLVED 2026-07-27.** Row 3 previously
   drew 14 noisy realisations in the same blue as rows 1–2 while the curve from rows 1–2 was
   absent — so the row introduced a fresh cast of objects at exactly the point where the design
   depends on one object persisting. Now the deterministic curve is plotted in row 3 in the same
   colour and weight as row 1, carrying the same red first-zero dot, and the ensemble is demoted
   to grey scatter *around* it. Colour now encodes role: **blue = the invariant, grey = scatter,
   red = zeros.** This is load-bearing for the whole α3–α4 chain — the fix was applied here
   before the pattern could propagate. See `DESIGN_inner_chain.md` §2.

---

## Module α4a-interactive (twin-panel scrubber) — distortion statement
Completed 2026-07-27 against `pilot_borel/twin_scrubber.html`. First web artifact.

1. **The blue curve is not Ai.** It is the first-order Euler solution of the same ODE, differing
   from the true Airy function by ~0.4% relative at N=30000. Stated in the page's own note. The
   R↔u correspondence is unaffected — it is exact for *any* u solving the equation — but the
   curve should not be read as a plot of a special function.
2. **The gap between where R leaves the frame and the dashed line is a rendering artifact.**
   R is clipped at |R|=12, and since R ≈ −1/(Y−z) near a zero, the curve exits about 1/12 ≈ 0.083
   before the actual pole. A viewer may read that offset as a real mismatch between the blow-up
   and the zero. It is not.
3. ~~**With noise on, one realisation is shown, not an ensemble.**~~ **RESOLVED 2026-07-27.**
   16 peers now draw in grey behind the scrubbed run, the η=0 solution stays visible as a dashed
   blue reference, the peers' first zeros appear as a red tick strip on the zero line, and the
   readout states the run's first zero against the peers' span, median, and its own percentile —
   so "is this typical?" is answered numerically rather than left to be assumed. Colour encodes
   role exactly as in the still: **blue = the trajectory in hand, grey = scatter, red = zeros.**
4. ~~**Only the recessive (Ai) boundary condition** … the caption from row 1 is *absent*.~~
   **RESOLVED 2026-07-27.** Carried over verbatim beneath the title.
   *General lesson filed below.*
5. **The noise discretisation is crude.** Euler–Maruyama at this step size is not trustworthy for
   η beyond roughly 1, which is why the slider stops at 1.2 rather than going further. The scheme
   is chosen to match the Python exactly, not to be the best available integrator.

---

## Standing entry — fixes do not propagate from a still to its interactive
The boundary-condition caption was added to `twin_panel.py` and then shipped missing from
`twin_scrubber.html`, reintroducing in the interactive a defect that had just been closed in the
figure. Nothing catches this: the two artifacts share a kernel contract (golden vectors) but no
*narrative* contract, so captions, caveats and on-screen distortion statements can silently
diverge. Until something better exists, **every ledger entry closed on a still must be re-checked
against its interactive before that interactive ships**, and vice versa.

---

## Module α3 (blow-up as anisotropic zoom) — distortion statement
Completed 2026-07-27 against `blowup/blowup.html`. **Wants expert correction.**

1. **The sphere never appears.** This shows a *chart* — a weighted rescaling of (x, y) — not the
   blow-up construction as usually presented. A viewer may conclude blow-up *is* this 2D
   rescaling. It is the part that carries the intuition, but the full object glues several such
   charts onto an exceptional divisor. Deliberate: DESIGN_inner_chain.md §4 says opening on the
   sphere is the standard reason this topic loses people.
2. ~~**ε is held fixed, so this is not the whole blow-up.**~~ **RESOLVED 2026-07-27 — expert said
   "clearly not the whole picture".** ε now has its own slider, and the widget carries all three
   weights. The trajectory is invariant only at (1, 2, 3), verified bit-for-bit: identical point
   count and endpoint (2.170, −0.542) at λ = 1, 32 and 1000.
3. ~~**Only the critical manifold is drawn — there are no dynamics.**~~ **RESOLVED 2026-07-27.**
   A trajectory of the view system is integrated by RK4 and drawn, showing the fold passage. The
   attracting/repelling labels are now earned by something on screen.
3b. **New, from that change: the payoff is now the EQUATION, not the picture.** The view system is
   dX/dT = X² − λ^(2−p)·Y, dY/dT = −λ^(1+p−q)·ε̄, displayed live with its exponents. They vanish
   together exactly at p=2, q=3 and λ leaves the equation. That — not a tidier picture — is what
   a blow-up buys, and it is the honest statement of the construction's purpose.
4. **y = x² is the normal form, not a generic fold.** A real fold is only locally, and after
   coordinate changes, of this shape. The widget silently works in already-normalised coordinates.
5. **Nothing is added to the space** — captioned on-screen, per the standing blow-up entry above.
6. **p is continuous, which is a pedagogical fiction worth keeping.** Letting the viewer slide
   through non-integer weights is what makes p=2 *discoverable*; but the weights of a blow-up are
   determined by the equation, not tuned. The widget is a search device, not a model of practice.

---

## Module α4c (shooting → the spectrum) — distortion statement
Completed 2026-07-28 against `shooting/shooting.html`. Closes the inner chain.

1. **The wall is imposed, not intrinsic.** u(0)=0 on a half-line is what makes the spectrum
   discrete; the Airy operator on the *whole* line has no eigenvalues at all. The escape problem
   is what puts the wall there. Captioned on-screen.
2. **u is normalised for display, so its height means nothing.** An eigenfunction is defined only
   up to scale. Only the zeros and the sign of the miss carry information. Captioned.
3. **The shot reaches ~5×10¹⁵ crossing the forbidden region.** Legitimate in double precision at
   these ranges, but it means any comparison of two implementations must be *relative*. An
   absolute tolerance of 1e−8 produced a false PORT MISMATCH badge on first run. Fixed; the check
   is now relative and reads 1.5e−14.
4. **"Interior nodes" excludes x=0 by construction.** At an eigenvalue u(0)=0 *is* the boundary
   condition, and its residual sign is numerical noise — counting that interval added a phantom
   node at some eigenvalues and not others. Count and markers both now start strictly inside.
   Verified: node count = n−1 for all eight eigenvalues below λ=14.
5. **The miss curve is precomputed and therefore fixed.** It cannot respond to any change of
   XMAX or scheme made later in the browser. If the kernel is edited, the curve must be rebuilt —
   the golden vectors will catch a divergence, but only at the one λ they sample.

---

## Module γ8 (Stokes phenomenon) — distortion statement
Completed 2026-07-28 against `stokes/stokes.html`. The fourth and last of the concepts the user
named as hardest to convey. **Wants expert correction — entries 2 and 5 especially.**

1. **Convention clash, and it is a real trap.** Some authors call Im ζ = 0 the *Stokes* lines
   (used here, after Berry); others reserve that name for Re ζ = 0. Reading across sources
   without checking which is meant gets the geometry exactly backwards. Captioned on-screen.
2. **The erf profile is least trustworthy exactly where it is most visible.** It is the leading
   term of Berry's smoothing, asymptotic in |ζ| — so at r ≈ 1, where the ramp is 45° wide and the
   demonstration is at its most persuasive, the formula is at its weakest. This tension is
   intrinsic to the module and is captioned, but a viewer will remember the wide ramp and not the
   caveat. **The strongest candidate for revision in this module.**
3. **The multiplier is plotted, not measured.** Berry's formula is drawn; nothing here extracts a
   Stokes multiplier from Ai(z) numerically and compares. The connection identity *is* verified
   (residual 1.5e−14, fixing the constant to exactly i), but that verifies the constant, not the
   switching profile.
4. **Only the Stokes line at arg z = 0 is analysed.** The other two rays are drawn and labelled
   but nothing is computed across them. Symmetry makes this harmless, and nothing says so.
5. **S₀ is drawn as a real ramp from 0 to 1, but the Airy Stokes constant is i.** The multiplier
   is *complex*; plotting |S|/S₀ as a real quantity hides that the switched-on term arrives with a
   phase. Defensible for a first pass — the magnitude is the story — but it is a real omission and
   a viewer would not guess it.
6. **The radius slider is logarithmic and stops at ~300.** Beyond that the ramp is narrower than a
   pixel, so the widget cannot show the limit it is describing; it can only show the trend toward
   it.

---

## Module δ1 (the junction) — distortion statement
Completed 2026-07-28 against `junction/junction.html`. The node where the two trunks meet.

1. **The closed form is a LEADING-ORDER action, and the page shows it as exact.** p = t^q is the
   leading behaviour of the optimal path, not the path. A boundary-value solve of the full
   Euler–Lagrange system approaches it *from below* — 0.839 of it at q=1 and s=10. The table
   reports that honestly, but the hero line still reads like an identity.
2. **This constant was wrong in this project's own notes.** An earlier draft carried
   s^(2q+1)/[4(2q+1)] — too small by a factor 2 — until the q=2 case was checked against a BVP
   solve and disagreed. Recorded on-page, because a viewer should know the formula has a history.
3. **The shaded region is drawn under p, not under p²/2.** It marks the *support* of the action
   integral, not its value; the height is p and the area is not the number quoted. Legible as a
   "cost accumulates here" cue, misleading if read quantitatively. **Weakest point of the module.**
4. ~~**The right-hand panel plots the tail exponent, not the coefficients.**~~ **RESOLVED
   2026-07-28.** A third panel now draws the Domb–Sykes ratios, and a static twin (`ratio_still.py`)
   renders the same numbers. But the resolution changed the claim, and the honest version is
   weaker and more interesting than the one it replaces:
   * The ratios **do not settle** to 1/A. They swing — sign flip at k=2, spike above 3.7 at k=3.
   * That swing is **diagnostic, not noise**: a conjugate pair at ρ=1.9, θ=50° with one fitted
     phase reproduces the whole sequence to **0.46%**.
   * But the **modulus is not recoverable from seven coefficients**. The root test |b_k|^(−1/k)
     plateaus at 2.5–2.7 against a true |ζ|≈1.9 and is not converging.
   So the junction's claim is now shown *and* bounded: the coefficients do encode the geometry, and
   seven of them cannot say so plainly. That is precisely why the programme needed matched-K
   calibration rather than a textbook Domb–Sykes — which the module now demonstrates instead of
   assuming.
4b. **The pair's ρ and θ are FIXED, not fitted.** Only the phase and an overall scale are free, so
   the 0.46% is a two-parameter check of a previously reported result, not evidence for it. A
   viewer could read the close overlay as independent confirmation of ρ=1.9, θ=50°. It is not.
5. **The ladder shown is the cusp only.** The swallowtail coefficients are actively being revised
   (v₇ changed on 2026-07-28) and are deliberately excluded rather than shown as settled.

---

## Module γ1 (the divergent series) — distortion statement
Completed 2026-07-28 against `divergence/divergence.html`. Opens the analysis thread.

1. **The series alternates, which makes it the EASY case.** Euler's series is Borel summable with
   no ambiguity. The series this course actually cares about does not alternate, and that
   difference is exactly where γ2–γ3 begin. Captioned, but a viewer meeting divergence for the
   first time may generalise "divergent but summable" too readily.
2. **Truncating at the least term is a rule of thumb, not a theorem.** Excellent here, standard
   practice, and asymptotic in nature. Captioned.
3. **The floor is e^(−1/x) only up to an algebraic prefactor.** Measured ratios of the minimum
   error to e^(−1/x) run 1.75 → 5.56 as x falls from 0.45 to 0.05; the prefactor is ~√(2π/x). The
   dashed line on the error panel is therefore *parallel to* the floor, not the floor itself.
4. **The left panel clips.** Partial sums past the optimal point leave the frame within a few
   terms; the vertical range is fixed at ±0.55 around the true value. Without clipping every
   other feature would be a flat line, but the divergence is visually understated as a result.

---

## Module α1+α2 (fast–slow, and where it breaks) — distortion statement
Completed 2026-07-28 against `fastslow/fastslow.html`. Merged deliberately: α3 opens on the
failure, so something must close on it, and "here is where the tools stop working" only lands once
the tools have been seen working.

1. **Planar FHN only.** In higher dimensions folded singularities acquire structure (folded nodes,
   folded saddles) and the story branches. Nothing here shows that, and a viewer could take the
   planar picture as the general case.
2. **λ = 1 − v² is a rate, glossed as "the strength of the pull".** It is the eigenvalue of the
   fast subsystem linearised at the manifold. The gloss is fair but it is a linearisation, and
   near the fold — precisely where the module points — linearisation is what stops being valid.
3. **The jump happens visibly PAST the fold.** At ε=0.08 the trajectory overshoots by an O(ε^(2/3))
   margin before departing. This is real and is exactly what α3 exists to resolve, but nothing
   on-screen says so, and a viewer may read it as sloppiness in the drawing.
4. **The critical manifold runs off the top and bottom of the frame.** The cubic exceeds the
   plotted w-range at |v| ≳ 2.2. Harmless, but the branches appear to terminate.
5. **"Fenichel does not apply" is stated, not demonstrated.** The module shows λ → 0, which is the
   hypothesis failing; it does not show a slow manifold failing to persist. That is the right
   scope — but the leap from "hypothesis fails" to "conclusion fails" is asserted.
