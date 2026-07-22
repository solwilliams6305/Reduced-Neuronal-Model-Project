# Program 2, Milestone 2 — the fluctuation→Borel-phase link: a rigorous PARTIAL/NEGATIVE result

_2026-07-10 (Fable 5). Attacked "Thread 4 / the λ₀→Borel bridge" — why the stochastic Borel phase is ~54° not
λ₀'s 45° — via a 9-agent adversarially-verified workflow (robust re-measurement + 3 independent mechanisms) plus
a decisive follow-up identifiability/contamination test. Context brief: `MILESTONE2_CONTEXT.md`. Scripts:
`coupled-atlas/milestone2_*.py`. Figure: `coupled-atlas/figures/milestone2_phase_identifiability.png`. Every
number below traces to code that was run and (for the mechanisms) independently re-checked by a skeptic agent._

## Headline (honest)
**The 54°-vs-45° gap that Milestone 2 set out to explain is NOT robustly established from the 6 available
coefficients.** The most important result is negative and it *reframes* the frontier: before asking "what
rotates 45°→54°," one must first establish that there is a rotation. From v0..v5 there is not — the phase is
**unidentifiable** in [35°, 56°], and the known real instanton pole **biases** the very estimators used
(Borel–Padé [2/2], Darboux) to read a true 45° as 54–64°. Two premises the program had been building on are
**corrected**: (i) "|ζ| ≈ |λ₀|" is an overfit artifact, and (ii) "54° is robustly ≠ 45°" holds only if one
trusts the least-certain coefficient. One genuine structural fact survives (the *sign* of any shift is forced).

## What the ~54° actually is (target, from recon)
arg of a complex-conjugate pair of Borel singularities of the variance series `Var = η²Σ vₙ η^{2n}`, extracted
by Borel–Padé + Darboux from `v = (0.134, 0.111, 0.100, −0.02, −0.45, −1.1)` (v0..v5). A competing **real** Borel
pole (the far-tail instanton `S=s⁵/10`, real, BVP-verified) sits at nearly the same modulus, arg 0. The
deterministic baseline is arg(λ₀)=−45° (λ₀=0.8896−0.8896i). The old "45°" was **hardwired** (a `cos(nπ/4)` fit),
not a floated measurement — so it was never independent evidence for 45° either.

## Result 1 — robust re-measurement: θ = 55 ± 9°, NOT robustly distinct from 45°
Five independent estimators (Borel–Padé across all [L/M]; floated-θ Darboux; γ-deflated Cauchy root; direct
conjugate-pair recurrence), all validated to recover θ=54.00°, |ζ|=1.258 exactly on a synthetic pure-pair
series, and cross-checked with a from-scratch mpmath dps=40 implementation. Joint bootstrap (N=8000) over the
stated coefficient uncertainties:
- **Full v5 range honored** (symbolic −1.1 vs grid −2.29, ~2× apart): θ median **49.5°**, 68% CI [40, 56],
  **P(θ>45°) = only 75%** — 45° sits inside the 68% band. v5 alone contributes sd(θ)=8.95° (dominant; v3 4.98,
  v4 3.51, v2 2.94).
- **Trusting symbolic v5≈−1.1**: θ = 56.7°, P(θ>45°)=98% — distinct *only* in this conditional.
- **Drop-v5 (reliable v0..v4 only)**: trusted estimators cluster tightly at **54.5°**; adding v5 *disperses and
  drifts* the estimate (54.6→61.8°) rather than tightening it — a warning the series is not yet asymptotic.

**Caveat that partially rescues distinctness:** the grid v5=−2.29 is from n=10, and n=10 *fails to reproduce the
known* v3,v4 (gives −0.62, −0.63 vs true −0.02, −0.45) — so the grid v5 is unconverged garbage and the symbolic
−1.1 is the better estimate. Under that reading the phase is likely ~54–57°. But "likely" is not "established."

## Result 2 — the radius clue is REFUTED (correction to the program)
`MILESTONE2_CONTEXT.md`/`PROGRAM2_CONSOLIDATION.md` asserted |ζ|≈1.44≈|λ₀|=1.258, motivating "ζ = noise-dressed
λ₀." **Not reproduced.** Constrained estimators give |ζ| = 1.45–1.81 (15–45% *above* |λ₀|); only 1.6% of
bootstrap samples land within 10% of 1.258. The single matching value (1.286) is a Darboux fit with 5 params on
5 points (max-residual 0.000 = exact interpolation), not a measurement. **|ζ| is essentially unmeasurable from 6
coefficients** (per-step estimates swing 1.2, 4.4, 45, 0.7, 10). The "radius match" was an artifact.

## Result 3 (this session's decisive test) — the phase is UNIDENTIFIABLE, and the real pole biases it
`milestone2_contamination.py`. Two findings:
- **Identifiability [panel A].** Fit `(pair@θ_true + real pole)` to the asymptotic tail v1..v5 and scan θ_true:
  the residual is **flat** (0.0037–0.0057, ~0.4% of |v5|) across **θ_true ∈ [35°, 56°]** — 45° (res 0.0042) fits
  as well as 54° (0.0037). The 6 coefficients **cannot discriminate** 45° from 54°.
- **Contamination [panel B].** A *true* 45° pair plus the known real instanton pole (|ζ|≈1.5) is read by the
  [2/2] Borel–Padé and Darboux estimators — the ones used to claim 54–63° — as an *apparent* phase rising to
  54–64° as the real-pole strength grows. And a true-45°+real-pole model **reproduces** the measured v0..v5 tail.
- **Estimators are unbiased on a pure pair** (45→45, 54→54), so the bias is specifically the real-pole
  contamination + Padé hypersensitivity (a ~2% coefficient change swings the [2/2] phase 45°↔54°).

**Conclusion:** the data are consistent with a true 45° (inherited from λ₀) contaminated by the s⁵/10 instanton
pole. The claimed 9–18° shift is within extraction uncertainty + contamination. This does **not** prove 45°; it
shows 54° is **not established** and the target must be pinned with more/better coefficients + real-pole removal.

## Result 4 — the four candidate mechanisms (all fail to DERIVE the magnitude)
| Mechanism | Prediction | Verdict (adversarial) |
|---|---|---|
| **A. complex instanton + one-loop determinant** (Bureković–Grauer) | −56.25° via `λ₀·D'^{−1/2}` | **REFUTED** — category error: a fluctuation determinant is a multiplicative *prefactor*; its phase is provably **invisible to the Borel singularity LOCATION** (synthetic test: varying arg C over {0,45,90,−90}° left the [2/2] pole pinned at −45°). The `−1/2` power was chosen post hoc to hit the target. |
| **B. resonance-cloud saddle** `E[λ^{−n}]` (built on Milestone 1) | −54° | **REFUTED** — `vₙ ~ E[λ^{−n}]` is *falsified* by the real ladder (E[λ^{−n}] is complex, phase-winding, *decaying* since \|λ₀\|>1; the true vₙ are real, sign-alternating, *growing*). The −54° is a free function of the product n·η² tuned to the target; no isolated Borel singularity (Darboux ratio doesn't converge in n). |
| **C. coalescence discriminant** `ζ²=λ₀²−c·pv` | −53.6° (c=1) | **PLAUSIBLE but a surmise** — reproduces from independent Milestone-1 inputs, but c is a free knob (c=1→−53.6°, 1.4→−57.3°, 2→−62.8°, spanning the whole band) and the *minus* sign (which does the work) is unproven (+pv rotates the wrong way, −38°). Radius 1.21 misses the measured 1.45–1.81. |

**The one surviving structural fact.** arg(λ₀²) = −90.0° *exactly* (λ₀ at −45°), and the Milestone-1
pseudo-variance `pv = E[δλ²]/η² = 0.436 − 0.181i` sits **off-radial** at arg −22.5°. Any fluctuation correction
to λ₀² of the collision form `λ₀² − c·pv` therefore acquires a *negative real part* and rotates below −90°,
pushing the half-angle arg(ζ)=½arg(λ₀²−c·pv) below −45°. So **the SIGN of a fluctuation-induced shift (toward
more-negative, i.e. 45°→54°) is structurally motivated; the MAGNITUDE is not pinned** (needs c, unproven). This
is the genuine, if modest, content of the λ₀→Borel bridge to date.

## What is ESTABLISHED vs CONJECTURAL
**Established** (code re-run + independently checked): the estimator machinery is correct and the Borel–Padé
poles are real (not float artifacts); the Milestone-1 inputs (λ₀, pv, Ω₂) reproduce exactly; the mean-shift
route is ruled out (Milestone 1); the single-kernel-eigenvalue route is ruled out (squared-Green chain is
1/k!-suppressed, plain-Green chain is a real-negative/alternating pole); the one-loop-det-phase route is ruled
out (prefactor ≠ location); the phase is unidentifiable from 6 coeffs and consistent with 45°+real-pole; the
radius clue is an artifact; the *sign* of any shift is structurally forced.
**Conjectural / open:** the exact phase, and whether it is even distinct from 45°; the magnitude of any shift;
the resolvent/spectral identity that would actually *link* the resonance-root fluctuation to the variance
coefficients (Mechanism B's version is falsified, C's is an unproven surmise). This remains the genuine frontier.

## Critical path (the single highest-value next step)
**Obtain a grid-CONVERGED, certified v5 (ideally v6), then re-extract with the real instanton pole SUBTRACTED
first.** v5 drives ~9° of the ~9° phase uncertainty; the n=10 small-grid route is a verified dead end (fails the
known v3,v4); the productive path is the transfer-operator engine (`chaos_transfer.py`) pushed to the continuum
with the δ(0) boundary-renormalization counterterm (= Handoff A, the "v6 grind"). Because a competing real pole
contaminates the extraction, the second half is essential: subtract the known `s⁵/10` instanton contribution
(its Borel location is fixed) before reading the complex pair, so the pair phase is not biased. Only with (a)
converged v5/v6 and (b) real-pole removal can the 45°-vs-54° question be settled and the mechanisms tested
against data rather than tuned to a pre-known target.

## Files
`coupled-atlas/`: `milestone2_measure.py`, `milestone2_bootstrap.py`, `milestone2_budget.py`,
`milestone2_verify_independent.py` (measurement + bootstrap, mpmath cross-check); `milestone2_mechA_instanton.py`,
`milestone2_mechA_oneloop.py`, `milestone2_mechB_cloud.py`, `milestone2_mechC_chaos.py`,
`milestone2_refute_check.py` (the mechanisms + refutations); `milestone2_contamination.py`,
`milestone2_figure.py` (this session's identifiability/contamination test + figure). Depends on Milestone 1
(`stochastic_stokes_o_eta2.py`) and the extractor `_borel_analysis.py`.
