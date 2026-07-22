# HANDOFF B — the stochastic exact-WKB frontier: derive the "stochastic Stokes constant" of 𝒲

_Self-contained handoff for a fresh Claude Code session. Route (b) of Program 2's endgame: the genuine
multi-year frontier (`PROGRAM2_PROVED_CLOSED_FORM_ROUTE.md` Phase 3, the one load-bearing BLANK). Scoped by a
105-agent deep-research pass (2026-07-10, `PROGRAM2_STOCHASTIC_WKB_DEEPRESEARCH.md`). This is a research project,
not a compute grind — treat it as such._

## The goal (one sentence)
Build exact-WKB / resurgence for the noise-dressed Weber operator and derive a **stochastic Stokes constant** —
the noise-dressed analogue of the deterministic connection root λ₀ — that (i) links 𝒲's perturbative and
instanton sectors and (ii) explains why the stochastic Borel phase is ~54–63°, NOT λ₀'s exact 45°. The endgame
is a **proved partial resurgent representation** of 𝒲.

## The object
𝒲 = first-node (first-explosion) law of the stochastic Weber operator `u'' = (sign(Y)·Y² − η·Ẇ(Y))·u`, β = 4/η²
(Ẇ = white noise in the spatial variable Y). The cusp/q=2 analogue of Tracy–Widom (which is the q=1 stochastic-
Airy case). No 1-D closed form exists (PROVED: isomonodromy fixed point ⇒ no Painlevé reduction; unbounded-below
first-passage ⇒ no Fredholm determinant). The only viable "closed form" is a resurgent trans-series.

## What is ALREADY PROVED / in-hand — the two halves to JOIN (do not rebuild)
- **Deterministic Weber exact-WKB is fully closed.** Nikolaev 2024 (arXiv:2410.17224) proves WKB resurgence on
  Riemann surfaces (meromorphic potentials). **Hao 2025 (arXiv:2507.06922) gives the closed-form Weber Stokes
  datum**: disc = −½·log(1+e^{−2πi/ε}), central charge Z_{γ_f} = −2πi, **deterministic Stokes constant
  Ω(±γ_f) = 1**, Borel poles at minus the central charges of 4d BPS states (8–15 digit verified). Iwaki–Koike–
  Takei (arXiv:1805.10945) give the Voros coefficients (Bernoulli-number form, topological recursion). **This is
  the exact deterministic baseline the stochastic Stokes constant is a shift OF.**
- **The program's own deterministic connection is closed-form and machine-verified.** Both factors are Weber
  Γ-ratios (confining × oscillatory) with an e^{3iπ/4} inversion phase; the resonance root is **λ₀ = 0.8896 −
  0.8896i (exactly −45°)**. Solvers: `coupled-atlas/connection_closed_form.py` (Jost/Γ-equation),
  `complex_scaling.py` (resonance = zero of rotated Wronskian D_θ(λ)). See `W_CLOSED_FORM_LEADS.md` Step 1c.
- **The stochastic half (rigorous but resurgence-blind).** The q=1 template — stochastic Airy / TW_β (Ramírez–
  Rider–Virág, arXiv:math/0607331) — has a PROVED sharp tail via the stochastic operator + Girsanov (Dumaz–Virág
  arXiv:1102.4818), but NEVER via resurgence. All-order tail only heuristic (Borot–Nadal loop equations) or
  two-term rigorous. **No random-operator edge law has ever been given a proved resurgent trans-series — 𝒲 would
  be first-of-kind.** The two halves have never been joined; joining them is the whole task.
- **A validated deterministic first-explosion PDE** (this session, `PROGRAM2_RICCATI_PDE_NOTES.md`,
  `coupled-atlas/_riccati_pde_probe.py`): 𝒲 = first passage of the Prüfer phase to π; the escape law solves a
  Scharfetter–Gummel / implicit Fokker–Planck PDE, validated to ~0.005 against the SDE at β=2. A working
  Monte-Carlo-free handle at the physical point (but NOT good for small-η coefficient extraction — see ruled-out).

## The scoped attack (from the deep-research, ranked)
**The bridge — co-equational parametric resurgence (Aniceto–Crew 2024, arXiv:2410.13690).** It rigorously
(systematic algebraic construction) encodes Borel singularities that MOVE on a parameter-dependent algebraic
curve Σ_z = {P_z(w,φ)=0} and can CROSS the integration contour — **the exact algebraic mechanism a noise-shifted
Borel phase (45°→54–63°) requires.** Caveat: it treats a DETERMINISTIC parameter z; **promoting z to the noise
coupling η / a noise realization, i.e. making Voros symbols random variables on a fixed deterministic Σ, is the
genuinely novel unproven step** — no one has done it. That step is the crux of route (b).

**The mechanism hint — Bureković–Schäfer–Grauer 2024 (PRL 133 077202, arXiv:2401.16264).** In weak-noise
stochastic NLS the bare Freidlin–Wentzell instanton action gives the WRONG tail; the **Gaussian-fluctuation
prefactor (one-loop determinant + zero mode) dominates**. Precedent that the fluctuation sector, not the classical
action, sets a noise-induced tail — plausibly why 𝒲's Borel singularity sits off the classical-action value.

## Concrete sub-attacks (in recommended order)
1. **The minimal target lemma (highest value, most likely to be provable).** Show that averaging the RANDOM Voros
   symbol over the noise — via a Malliavin/Wick (Nourdin–Peccati) expansion — leaves the deterministic Borel-
   singularity LOCATIONS fixed while shifting only the Stokes CONSTANT. That would *define* the stochastic Stokes
   constant rigorously as 𝔼[random Voros datum] and reduce the frontier to a computation. Start analytically on
   the O(η²) term (first noise correction to Ω=1); the program's Green's-function/Wick machinery
   (`weaknoise_greens.py`, the chaos engine) already builds the relevant iterated integrals.
2. **The phase-shift mechanism (decisive numerical experiment).** Compute the COMPLEX instanton / turning-point
   action of the Weber weak-noise problem AND its one-loop Gaussian-fluctuation determinant, and test whether
   including the fluctuation shifts the effective Borel-singularity phase from λ₀'s 45° toward the observed ~54°.
   This directly tests the Bureković–Grauer mechanism for 𝒲. Program tools: `instanton_action.py`,
   `mc_instanton_check.py`, `PERSISTENCE_ITEM2_NOTES.md` (the s⁵/10 tail).
3. **Co-equational Σ_z identification.** Work out what the parameter z and curve P_z(w,φ) ARE for the Weber
   weak-noise series (likely z = the spectral/shift parameter λ, or η²), then test numerically whether the
   measured Borel phase falls on the resulting singularity track. Needs the conceptual bridge first (§ "novel
   step" above).
4. **Closure toward a theorem.** Assemble: coefficients as proved renormalized chaos integrals (Phase 1 — Nourdin–
   Peccati, tools exist) + Borel singularity at the instanton action (Phase 2 — lift Spada/Berry–Howls finite-D
   thimble argument) + the stochastic Stokes constant (Phase 3 — sub-attacks 1–3) + moment-determinacy closure
   (Phase 4 — Tier A already closed). See `PROGRAM2_PROVED_CLOSED_FORM_ROUTE.md` for the full route table.

## Ruled-out routes (do NOT re-attempt — burned this session, see `PROGRAM2_STOCHASTIC_WKB_DEEPRESEARCH.md`)
- **Naive noise-shifted resonance via complex scaling: ill-defined.** You cannot complex-scale a white-noise
  realization (the noise lives on the real axis, not the rotated contour). This is exactly why the resonance route
  is the frontier and not a quick calc.
- **FP-generator complex spectrum: wrong object.** The non-self-adjoint first-explosion generator's complex
  eigenvalues sit at ~2–10° (not 54°) and coalesce at a Y-dependent exceptional point — a LOCAL decay phenomenon,
  not the GLOBAL connection-resonance that carries the Borel pair. Confirmed negative.
- **PDE small-η coefficient re-extraction: defeated.** The deterministic PDE hits the same small-η resolution wall
  (narrow escape peak, width ~η) the chaos engine was built to dodge. Good at β=2, useless for weak-noise v_n.
- **PIV / parabolic-cylinder Fredholm reduction: known-distinct.** Xia–Xu–Zhao (arXiv:2301.05807) give a Weber-
  kernel Fredholm det for Clarkson–McLeod PIV, but the program already PROVED 𝒲 is NOT the standard PIV σ-form
  (`STAGE1_GAP_DETERMINANT_NOTES.md`). A check to run, not a solved reduction.

## What "done" looks like (milestones, not all-or-nothing)
- **Milestone 1 (reachable): ✅ DONE 2026-07-10** — stochastic Stokes constant DEFINED as 𝔼[random Voros datum]
  and its O(η²) value COMPUTED + certified: **Ω₂ = −0.451 + 0.352i** (|.|=0.572, arg 142°≈3π/4). Independent
  exact-Newton validation 0.67%. See `PROGRAM2_STOCHASTIC_STOKES_O_ETA2_NOTES.md`,
  `coupled-atlas/stochastic_stokes_o_eta2.py`. Two structural payoffs: (i) the 45°→54° shift is NOT a mean-shift
  effect at O(η²) (mean rotates only −3.2°/η², stays ≈−45°; the shift is nearly radial) — it lives in the
  FLUCTUATION sector (rms|δλ|≈1.34η, anisotropic), sharpening the Bureković–Grauer hint; (ii) the minimal lemma
  is corrected — the Borel *location* (period) mean shift is δ(0)-divergent (the v3 renorm threshold), so the
  renormalization-CLEAN object is 𝔼[connection root], because the log-derivative is a non-local Ito functional.
- **Milestone 2: ATTACKED 2026-07-10 → the premise is SHAKIER than thought.** A 9-agent adversarial workflow +
  identifiability test (`PROGRAM2_MILESTONE2_FLUCTUATION_BOREL_NOTES.md`) found: (i) the 54°-vs-45° gap is **NOT
  robustly established** from the 6 coefficients — the phase is *unidentifiable* over [35°,56°] and a true 45°
  contaminated by the known real `s⁵/10` instanton pole is read as 54–64° by the very estimators used; (ii) the
  "|ζ|≈|λ₀|" radius clue is an **overfit artifact** (constrained |ζ|=1.45–1.81, not 1.258); (iii) three concrete
  mechanisms were tested and **all fail to derive the magnitude** (A one-loop-det = category error, prefactor≠
  location, REFUTED; B resonance-cloud `E[λ⁻ⁿ]` = assumption falsified by the ladder, REFUTED; C coalescence
  discriminant `ζ²=λ₀²−c·pv` = plausible but a surmise with a free knob). The **one survivor**: the *sign* of any
  shift is structurally forced (arg λ₀²=−90° exact, pseudo-variance off-radial at −22.5°), the magnitude is not.
  **Critical path:** a grid-converged, certified v5/v6 (= Handoff A "v6 grind" via the transfer engine + δ(0)
  counterterm) AND real-pole subtraction before extraction — only then can 45-vs-54 be settled.
- **Milestone 3 (the theorem):** a proved partial resurgent representation of 𝒲 — the "theorem-grade milestone"
  of `PROGRAM2_CONSOLIDATION.md` §4. Remaining bridge: map the resonance-root shift Ω₂ ↔ Hao's deterministic
  Stokes constant Ω=1 residue (is the stochastic Stokes constant `1 + η²·(...)` in Hao's normalization?).

## Primary sources (arXiv)
Nikolaev 2410.17224 · Hao 2507.06922 · Iwaki–Koike–Takei 1805.10945 · Aniceto–Crew 2410.13690 · Bureković–
Schäfer–Grauer 2401.16264 · Ramírez–Rider–Virág math/0607331 · Dumaz–Virág 1102.4818 · Borot–Nadal 1111.2761 ·
Cleri–Dunne 2002.06270 (resurgence survives dropping integrability — but on the wrong axis, add-randomness ≠
drop-integrability) · Spada/Serone (SISSA thesis) · Berry–Howls · Sueishi–Kamata–Misumi–Ünsal JHEP 12(2020)114
(DDP dictionary: saddle ambiguity ↔ Stokes-curve topology change) · van Spaendonck–Vonk 2204.09062.

## Program files to build on
`coupled-atlas/`: `connection_closed_form.py`, `complex_scaling.py`, `connection_poles.py` (deterministic λ₀);
`instanton_action.py`, `mc_instanton_check.py` (instanton sector); `weaknoise_greens.py`, `chaos_diagram.py` +
`chaos_transfer.py` (coefficients as Wick integrals); `_riccati_pde_probe.py` (the working FP-PDE). Notes:
`PROGRAM2_PROVED_CLOSED_FORM_ROUTE.md`, `PROGRAM2_STOCHASTIC_WKB_DEEPRESEARCH.md`, `W_CLOSED_FORM_LEADS.md`,
`PROGRAM2_RICCATI_PDE_NOTES.md`, `PROGRAM2_KICKOFF_STOCHASTIC_EXACT_WKB.md`, `T2_CONNECTION_DATA_DERIVATION.md`.
