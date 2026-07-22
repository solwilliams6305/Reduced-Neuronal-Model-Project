# Program 2 — CONSOLIDATION (the stochastic exact-WKB / resurgent trans-series of 𝒲)

## 000. [2026-07-20 — the representation is VALIDATED end-to-end; paper drafted]

**Thread 3 (reconstruct the real law from the trans-series) is CLOSED, positively.** Median/
lateral Borel–Padé resummation of the 7 exact coefficients v₀…v₆ reproduces the MC-free FP-PDE
ground truth for f(x)=Var(Y★)/η²:
- **At β=2 (x=η²=2, OUTSIDE the radius x_c≈1.2): resummed f = 0.235 vs ground truth 0.237 (~1%)**,
  while the naive 7-term partial sum = −166 (wrong by orders of magnitude AND sign). The diagonal
  [3/3] Padé is the winner; off-diagonal orders scatter (0.18–0.21, [4/2] fails) — the honest
  "7 terms is few" caveat. Lateral-contour median is angle-independent (numerically exact).
- The resummed curve tracks ground truth across x∈[0,2.4] (through the radius to β=2 and beyond).
- fp_cusp re-validated: β=2 fingerprint skew +0.601, exk −0.243 (vs trusted +0.607, −0.237).
- **This is the strongest evidence to date that the resurgent representation IS 𝒲.**
- Tools: `coupled-atlas/_w_groundtruth.py` (FP-PDE Var(η)), `_w_resum.py` (median Borel–Padé
  lateral), `_w_resum_figure.py`. Paper: `W_ResurgentTransseries_paper.tex` (7pp, 3 figs, builds
  clean) — the first consolidated write-up of the partial resurgent representation.
- **Second observable (κ₃/skewness) — ATTEMPTED, PARTIAL.** Computed the third-cumulant ladder
  κ₃=η⁴Σ t_j η^{2j}, t_j=Σ_{a+b+c=4+2j} cum₃(Y_a,Y_b,Y_c), via triple Y-products on the parallel
  engine. t₀…t₅ (n=20) = (0.059,0.153,0.191,−0.187,−2.07,−4.7), factorially divergent — resurgent
  structure confirmed in a SECOND cumulant. **Leading skew coeff t₀/v₀^1.5 = 1.21 ≈ known 1.16 ✓**;
  the ×2.8 skew overshoot (naive √2·1.21≈1.7 vs true 0.601) reproduced from the exact t₀ ✓. BUT the
  full skewness resummation S(x)=T/V^1.5 is NOT clean: violently divergent (s₅~50 vs v₆~−1.9),
  t₄/t₅ grid-sensitive (t₅: −1.73@16→−4.73@20), Borel–Padé orders scatter 0.4–2.2 at β=2 (gt 0.601).
  So skewness CORROBORATES the structure but does not give a clean β=2 reconstruction — the VARIANCE
  (§000 above, ~1%) stays THE quantitative validation. A clean skew resum would need t₆,t₇ (bigger
  triple-product compute). Tools: `_w_kappa3.py`, `_w_skew_resum.py`. Also the mean (§ tried) is
  confounded (absolute/reference-dependent; fp_cusp mean numerically stable so it's a definition
  mismatch + tail/instanton dominance) — only CENTERED observables cross-check cleanly.
- Also this session: the `_contract_hard` hybrid rewrite is DONE+correct but is NOT the v₇ gate
  (numba DP already beats it; real wall = assembly scale) — see `PROGRAM2_V6_GRIND_NOTES.md`
  2026-07-20 addendum. v₇ stays finalized-on-v₆.

---


## 00. SUPERSEDING STATE [2026-07-19 — GOVERNS everything below; v₆ now computed EXACTLY]

**Read `PROGRAM2_V6_GRIND_NOTES.md` for the full record.** The "one remaining large build" of §0 (exact v₆)
is DONE, and it overturns the §0/§3 phase and sign claims. Where §00 conflicts with §0 or §§1–5, §00 governs.

- **Exact ladder (grid-certified, dual-engine, boundary systematic 0.1%):**
  $$v_0{=}0.134,\ v_1{=}0.111,\ v_2{=}0.104,\ v_3{\approx}{-}0.030,\ v_4{=}{-}0.451{\pm}.003,\ v_5{=}{-}1.19{\pm}.01,\ v_6{=}{-}1.90{\pm}.15$$
- **RETRACTED — "$v_6>0$ / sign flip to positive":** v₆ is NEGATIVE. The §3quinque Darboux prediction
  ($v_6{\approx}{+}2$) and the §3sept direct-Itô-MC "corroboration" ($v_6>0$) are both REFUTED by the exact value.
  Sign pattern is $+,+,+,-,-,-,-$ (four-long negative run), not $\dots,-,-,+$.
- **RETRACTED — "$\theta\approx54$–63°, robustly NOT $45°$":** with 7 coefficients the Borel-pair phase fits to
  **$\theta\approx50°\pm2°$**; $\theta\ge53°$ is now hard-excluded (108/108 robustness combos; sign-pattern alone
  kills $\theta\ge60°$). The §3bis/§3quater $54°$ and §3quinque $62$–$64°$ readings were 5/6-coefficient artifacts
  (real-pole bias + unidentifiability — exactly the §0 walk-back's diagnosis, now confirmed).
- **Exact-$45°$ (λ₀ inheritance) is DISFAVORED too** (not just $54°$): the prescribed real-instanton subtraction
  (z_r-profiled, dual-implementation) leaves θ pinned at $50°$ in every trustworthy fit, with a NEGATIVE fitted
  real-pole residue (the 45°-masquerade needs positive). 45° survives only via α-railed overfits. **Leading
  picture: a genuine small (~5°) noise rotation off λ₀'s 45°, plus a strongly renormalized modulus.**
- **Radius:** $|\zeta|\approx1.8$–$2.6$ (single-pair), NOT $1.2$. The "$|\zeta|\approx1.2\approx|\lambda_0|{=}1.26$"
  match is dead (the §0 walk-back's "overfit artifact" call, confirmed and extended). Any Phase-3 inheritance
  proof must now explain phase-near-inheritance WITH modulus renormalization.
- **What SURVIVED unchanged:** factorial divergence (envelope $|v_6|{>}|v_5|{>}|v_4|$, model-free); the
  complex-conjugate-pair + real-$s^5/10$-instanton three-singularity structure (still unseparated — needs v₇);
  the ½-boundary convention and the renormalization validation.
- **v₇ (the discriminator for 45°-vs-50° and the real-instanton separation): NOT reached** — 5 engine
  iterations, walled by the n^(#chains) cost of Var(Y₈); scoped as future work (one `_contract_hard` transfer
  rewrite). Pre-registered gates: $v_7\approx{+}1.3$ (45°) vs ${+}8.2$ (50°). **FINALIZED on v₆.**

---

## 0. FINAL STATE [2026-07-09 capstone — SUPERSEDED by §00 on all phase/sign/radius claims; §§1–5 below are the detailed/chronological audit trail]

**What 𝒲 is.** The universal noise-induced-escape edge law at a coupled-FHN cusp — the cusp analogue of Tracy–Widom;
the first-explosion law of the stochastic Weber operator $u''=(\operatorname{sign}(Y)Y^2-\eta\dot W)u$, $\beta=4/\eta^2$.
A literal 1-D closed form is **proved impossible** (isomonodromy fixed point ⇒ no Painlevé reduction; unbounded-below
first-passage ⇒ no Fredholm determinant), so the only viable "closed form" is a **resurgent trans-series**.

**What was established (this program), all cross-validated.** Built a deterministic Wiener-chaos "diagram" engine
(`coupled-atlas/chaos_diagram.py`) that computes the weak-noise variance coefficients **exactly** (no Monte Carlo, no
truncation). The exact coefficient ladder $\operatorname{Var}(Y^\star)=\eta^2(v_0+v_1\eta^2+\dots)$:
$$\boxed{\,v_0=0.134,\ v_1=0.111,\ v_2\approx0.104,\ v_3\approx-0.03,\ v_4\approx-0.45,\ v_5\approx-1,\ v_6>0\,}$$
> **⛔ This 07-09 box is SUPERSEDED — see §00.** Exact values: $v_4={-}0.451$, $v_5={-}1.19$, and $v_6={-}1.90$
> (NEGATIVE — the "$v_6>0$" here is REFUTED).
- **The series is FACTORIALLY DIVERGENT** (Darboux/Dingle exponent $\alpha\approx0$) — a genuine asymptotic/resurgent
  trans-series, the program's defining premise, now numerically established.
- **Borel plane:** a **complex-conjugate singularity pair at $|\zeta|\approx1.2$, phase $\theta\approx54$–63°**
  (Borel–Padé + Darboux), plus a **competing REAL singularity** = the far-tail Freidlin–Wentzell instanton
  $e^{-s^5/10\eta^2}$. The coefficients OSCILLATE (period $2\pi/\theta\sim6$–8): signs $+,+,+,-,-,-,+$ (first flip at
  $v_3$, small; $v_6>0$).
- **KEY correction to the earlier surmise:** $\theta\approx54$–63° is robustly **NOT** the deterministic connection
  root $\lambda_0=0.890-0.890i$'s $45^\circ$ — so "the Borel phase is inherited from $\lambda_0$" is **refuted** at the
  $\sim10$–18° level. The $\lambda_0\to$stochastic-Borel bridge is a genuine open (frontier) problem, not a coincidence.
- **⚠️ 2026-07-10 WALK-BACK (`PROGRAM2_MILESTONE2_FLUCTUATION_BOREL_NOTES.md`): the two bullets above OVERCLAIM.**
  A 9-agent adversarial workflow + an identifiability test show that from only 6 coefficients the phase is
  **unidentifiable over [35°,56°]** (flat fit residual — 45° fits as well as 54°), and the **known real $s^5/10$
  instanton pole biases** the [2/2] Borel–Padé and Darboux estimators to read a *true* 45° as 54–64°; a
  true-45°+real-pole model reproduces v0..v5. So "$\theta$ robustly $\neq45°$" holds ONLY if one trusts the
  symbolic $v_5\approx-1.1$ over the (unconverged, n=10) grid $v_5\approx-2.29$; under the full stated uncertainty
  P($\theta>45°$)$\approx$75%, and 45° is inside the 68% band. **The "$|\zeta|\approx1.2\approx|\lambda_0|$" match
  is an OVERFIT ARTIFACT** — constrained estimators give $|\zeta|=1.45$–1.81. Net: the shift is *not established*;
  only its *sign* is structurally motivated. Settling it needs a converged $v_5/v_6$ + real-pole subtraction.
- **VALIDATION (the crucial gate, passed).** A deep-research pass ranked "validate the renormalization" as the
  prerequisite (dropping boundary self-contractions could silently omit a Wong–Zakai/Itô counterterm). Checked against
  the **ground-truth direct Itô first-passage MC**: (i) the engine's $v_3,v_4,v_5$ reproduce the true Var($\eta$); (ii)
  the $\tfrac12$ boundary-Wick convention is CORRECT (bfac$=\tfrac12$ matches the MC; bfac$=0$ overshoots — and this is
  decisive because $v_3$'s SIGN flips with the convention); (iii) the MC independently confirms $v_6>0$. So the engine,
  the renormalization scheme, and the boundary convention are all validated against the physical Itô law.

**Honest ledger.** SOLID: $v_0..v_4$ exact + MC-validated; factorial divergence; the complex-pair+real-instanton Borel
structure; the $\tfrac12$-convention. ROUGH: $v_5\approx-1$ (single coarse grid), $v_6>0$ (weak MC). CONJECTURE: exact
$\theta$ (somewhere 54–63°), radius $|\zeta|\approx1.2$. OPEN/FRONTIER: the $\lambda_0\to$Borel bridge (Thread 4, blank);
a PROVED resurgent trans-series (needs rigorous stochastic exact-WKB — no such framework exists yet).

**The one remaining large build (deferred).** Precise $v_6,v_7$ — to pin $\theta$ and separate the real instanton —
needs the transfer-operator engine completion (self-contraction core validated; general coupled-chain contractor +
order-13 functionals remain, ~multi-hour). It would SHARPEN, not overturn, the validated picture.

_Detail notes: `PROGRAM2_ROUTE2B_NOTES.md` (coefficients, renormalization threshold), `PROGRAM2_CHAOS_ENGINE.md`
(the engine), `PROGRAM2_EXACT_FORM_SURMISE.md` (the conjectured form), `PROGRAM2_PROVED_CLOSED_FORM_ROUTE.md` (proof
route). §§1–5 below preserve the chronological derivation; where they predate §0 (e.g. "phase = 45°"), §0 governs._

---

_2026-07-09 (Fable 5). Capstone of a multi-session push. Program 2 began at "$v_2$ is the open question."_

## 1. The headline result — the exact coefficient ladder
$$\operatorname{Var}(Y^\star)=\eta^2\big(v_0+v_1\eta^2+v_2\eta^4+v_3\eta^6+v_4\eta^8+\cdots\big),$$
$$\boxed{\,v_0=0.134,\quad v_1=0.111,\quad v_2=0.100,\quad v_3\approx-0.02,\quad v_4\approx-0.45\,}$$
all computed **exactly** (deterministic Wiener-chaos engine, no Monte Carlo, no truncation; validated on the
independently-known $v_0,v_1$). The sequence is **oscillatory with a growing envelope** ⇒ the series is a genuine
**divergent (asymptotic) resurgent trans-series** — the defining premise of Program 2, now numerically established.

## 2. The conjectured exact form [CONJECTURE, evidence-backed]
$\mathcal W_\beta$ has no 1-D closed form (proved: isomonodromy fixed point + unbounded-below first-passage). Its exact
form is a **resurgent trans-series in $\eta^2=4/\beta$** whose Borel plane carries a **complex-conjugate pair of
singularities at $\arg=\pm\pi/4$**, inherited from the deterministic Weber connection root $\lambda_0=0.890-0.890i$
(exactly $-45^\circ$). Hence the perturbative coefficients **oscillate** (period $8$ in $n$), with radius
$\eta^2_c\approx0.6$–0.7 ($\beta_c\approx6$); $\beta=2$ sits outside the radius (non-perturbative — the measured
$\times2.7$ skew overshoot). Plus a separate **real** far-tail instanton $e^{-s^5/10\eta^2}$. Evidence: the exact ladder
fits $v_n=C r^{-n}\cos(n\pi/4+\varphi)$ with $\theta=\pi/4$ FIXED to $\lambda_0$; the first sign flip is forced to
$v_3$ (near a cosine zero, hence small — exactly as computed) and $v_4$ large (anti-node) — both confirmed.

## 3. What was built / established (this program)
- **Step 0 [DONE, session 1 legacy]:** leading coefficients pinned, FP artifacts retracted ($v_0=0.134$, etc.).
- **$v_2\approx0.10$ [DONE]:** three agreeing routes (full-field MC, regularized recursion, deterministic Green's/Wick).
- **The renormalization threshold [DERIVED, exact]:** $v_3$ is the first coefficient needing a boundary counterterm.
  The node-shift recursion's $Y_7$ has a unique $s_0^2$ term $=\tfrac15Y_1^5$; $v_3$'s divergence is
  $\tfrac{23}{3}\langle\xi(Y^\star_0)^2\rangle v_0^3$ (from $Y_4{=}\tfrac13\xi Y_1^3$ and $R_7{=}\tfrac15Y_1^5$). No
  cubic boundary terms; $v_4$ escalates to $\xi'(Y^\star_0)^2$ ($1/\delta^3$), coefficients pure-$Y_1$.
- **The deterministic chaos engine [BUILT, VALIDATED]:** `chaos_diagram.py` — Wiener-chaos moments via contraction
  multigraphs (Markovian bond-chain kernels, boundary-Wick $\tfrac12 a^{(j)}(\text{node})$, renorm = drop
  boundary-boundary). Exact, MC-free, ~1s per coefficient through $v_3$, ~15s for $v_4$. Reproduces $v_0,v_1$ exactly.
- **$v_3\approx-0.02$, $v_4\approx-0.45$ [COMPUTED, exact]:** ⇒ divergent/resurgent verdict (§1).
- **Proof route [CHARTED, cited]:** deterministic Weber resurgence is PROVED (Nikolaev 2024); instanton→Borel-singularity
  proved in finite-D (Spada/Berry-Howls); chaos/cumulant rigor tools standard (Nourdin-Peccati). The one genuine BLANK
  is **stochastic exact-WKB** (no rigorous framework for a random-potential operator exists) — the load-bearing novel step.

## 3sept. $v_6>0$ from the direct Itô MC — oscillation independently corroborated [NUMERIC, weak, 2026-07-09]
> **⛔ RETRACTED by §00 (2026-07-19).** REFUTED by the exact $v_6=-1.90<0$. This weak-MC "corroboration" of
> $v_6>0$ (2-param +0.36, 3-param +3.1) was a second independent line that failed together with §3quinque's
> Darboux extrapolation — a cautionary case that $\pm0.15$-grade MC fits to a divergent series are unreliable
> even for SIGNS. See §00.
With $v_0..v_4$ now VALIDATED (§3sext), extracted $v_5,v_6$ as a constrained fit of the direct-Itô-MC Var($\eta$) residual
($\eta=0.50$–0.78, Richardson $h\to0$) after subtracting the fixed $v_0..v_4$. Fits are imprecise (divergent series, 5
points): 2-param gives $v_5{=}{+}0.17,v_6{=}{+}0.36$; 3-param gives $v_5{=}{-}0.57,v_6{=}{+}3.1,v_7{=}{-}2.5$ — magnitudes
unreliable, BUT **$v_6>0$ in both**, and the 3-param fit shows the oscillating pattern $v_5<0,v_6>0,v_7<0$. So the direct
Itô law **independently supports $v_6>0$ (the predicted sign flip to positive)**, corroborating the complex-pair/period-8
oscillation. Precise $v_6$ (to pin the phase and separate the real instanton) still needs the engine's exact value = the
transfer-operator rewrite (validated core; general coupled-chain contractor + order-13 functionals remain). **Net: the
qualitative resurgent picture is now cross-validated by BOTH the exact engine and the direct MC; only the precise phase
awaits $v_6,v_7$.**

## 3sext. RENORMALIZATION VALIDATED against the direct Itô law — but a boundary-convention systematic found [DECISIVE, 2026-07-09]
_Deep-research (partial; hit session limit) ranked **validate-the-renormalization (Option C) as the gate** — it confirmed
(3-0 claims) the Itô–Stratonovich/Wong–Zakai counterterm $\tfrac12 g'g$ is real and that "dropping quadratic
self-contractions = dropping the degree-2 term of a cumulant-generating function," a finite shift that could bias
$v_3,v_4,v_5$; and that conformal-Borel/Meijer-G (Option B) can't help until the coefficients are trusted (they need the
singularity structure imposed + external cross-checks). So I ran the decisive test it named: direct Itô full-field
first-node MC of Var($\eta$) vs the engine's coefficients._

**Result [NUMERIC, decisive].** At $\eta=0.45$ (cleanest), Var$_{\rm MC}=0.032383$. Reconstructing with the independent
session-1 value $v_2=0.114$ **plus the engine's $v_3,v_4,v_5$** gives $0.03238$ — a near-exact match; using the engine's
own $v_2=0.10$ gives $0.03225$ (too low). **So: (i) the renormalization for $v_3,v_4,v_5$ is NOT broken** — the
higher-order corrections track the true Itô law once $v_2$ is right (the drop-self-contraction scheme did NOT silently
omit a counterterm at $v_3{+}$, contrary to the Option-C fear). **(ii) BUT the engine's $v_2=0.10$ is ~14% below the true
$v_2\approx0.114$**, and $v_2$ is the FIRST coefficient carrying boundary (`beval`) terms — so the boundary-Wick
convention (the $\tfrac12 a^{(j)}(\text{node})$ rule) carries a **~14% systematic that enters at $v_2$ and propagates to
all boundary coefficients**. The `beval` *primitives* were validated ($\langle s_0 I_1(\psi)\rangle=\tfrac12\psi$), but the
overall boundary NORMALIZATION vs the true Itô law is off by ~14%.

**Consequences.** The QUALITATIVE conclusions survive (factorial divergence, complex pair, phase $\sim54$–63°).

**RESOLUTION [2026-07-09, corrected] — the boundary convention is VALIDATED; the "~14% systematic" was an over-read.**
Made the boundary-Wick factor tunable (`bfac`) and pinned it against the direct Itô MC. CRITICAL: $v_3$ changes SIGN with
the convention — bfac=0 (naive "Itô: non-anticipating field excludes the node increment") gives $v_3=+0.05$ (no sign
flip!), bfac=0.5 (symmetric $\tfrac12$) gives $v_3=-0.03$. So the choice is decisive, and the direct Itô MC settles it:
comparing engine Var($\eta$) to MC at $\eta=0.45,0.55$, **bfac=0.5 matches** (Var(0.55): 0.0527 vs MC 0.0524, the residual
fully explained by small-grid coefficient bias) while **bfac=0 overshoots badly** (0.0537). So $\tfrac12$ is CORRECT
(the boundary term comes from a smooth-field finite-difference derivative straddling the node → the Wong–Zakai/Stratonovich
$\tfrac12$, even for the Itô SDE). **Net: the engine's $\tfrac12$ convention reproduces the ground-truth Itô law; the
renormalization AND the boundary convention are both validated; $v_3\approx-0.03$ (negative, real sign flip) stands, and
with it the factorial-divergence + complex-pair + phase-$\sim54$–63° conclusions.** No fix needed. The earlier
"$v_2\approx0.114$, engine 14% low" was noisy MC-extraction uncertainty — the direct Var comparison is only ~2σ at clean
$\eta$, consistent with the engine's grid-converged $v_2\approx0.104$.

## 3quinque. Darboux/Dingle inversion on the 6 coefficients (Thread 1 analysis) [NUMERIC, 2026-07-09]
> **⛔ RETRACTED by §00 (2026-07-19).** The exact $v_6=-1.90$ REFUTES this section's headline prediction
> $v_6\approx+2$ (sign flip) and its $\theta\approx63°$. These were 6-coefficient artifacts. See §00.
Fitting the late-terms law $v_n\sim2|C|\,|\zeta|^{-n}\Gamma(n+\alpha)\cos(n\theta-\varphi)$ (complex-conjugate Borel pair)
to $v_0..v_5$ is **robust to the rough $v_5$**: across $v_5\in[-1.0,-1.3]$, $|\zeta|=1.17$–1.20, $\theta=62$–64°,
$\alpha\approx0$ (residual $<0.02$). So: **factorial divergence confirmed** (predicts $v_6\approx+1.4$ to $+2.6$ — a
**sign flip to positive**); Borel pair at $|\zeta|\approx1.2$, **$\theta\approx63^\circ$**. Combined with Borel–Padé's
$54^\circ$ (+ separate real pole), the phase is robustly **$\sim54$–63°, NOT $\lambda_0$'s $45^\circ$** — the $\lambda_0$
inheritance is refuted at the $\sim10$–18° level (the true stochastic phase is higher). The real instanton can't be
added (8 params $>$ 6 data) — disentangling it needs $v_6,v_7$. **Sharp prediction for the engine rewrite: $v_6\approx+2$.**

## 3quater. SIXTH coefficient + Borel–Padé — complex pair + real instanton resolved [NUMERIC, 2026-07-09]
> **⛔ RETRACTED by §00 (2026-07-19).** The "$v_5$" here was the pre-asymptotic n=14 grid value; the converged
> $v_5=-1.19$ and exact $v_6=-1.90$ give $\theta\approx50°$ (not $54°$), $|\zeta|\approx1.8$–2.6 (not $1.44$).
> The complex-pair + real-instanton STRUCTURE survives; the specific phase/radius numbers do not. See §00.
Computed $v_5$ exactly (diagram engine extended to order 11: generated $Y_{10},Y_{11}$ by truncated-series inversion;
6-chaos build; base-atom forms cached). $v_5\approx-1.1$ (n=14 grid gives $-1.32$; applying the $v_4$ grid-correction
$\sim1.2\times$ $\Rightarrow\approx-1.1$). Ladder: $v=(0.134,0.111,0.100,-0.02,-0.45,-1.1)$. $|v_5|>|v_4|$ — envelope
still GROWING (divergence reconfirmed); signs $+,+,+,-,-,-$ fit the period-8 pattern.
**Borel–Padé$[\cdot/4]$ on the 6 coefficients (robust across $v_5\in[-1.0,-1.3]$):**
- **complex-conjugate pair** at $\zeta\approx1.44\,e^{\pm i\,54^\circ}$ (phase 53.4–54.6°, very stable);
- **a REAL pole** at $\arg=0$, $\zeta\approx1.5$ — the **third Borel singularity** flagged by the deep-research pass:
  the real far-tail instanton $s^5/10$ on the real axis, at nearly the same $|\zeta|$ (so it competes).

**Two honest refinements:** (i) the complex-pair phase is robustly $\approx54^\circ$, **NOT exactly $\lambda_0$'s
$45^\circ$** (off $\sim9^\circ$) — so "$\pm\pi/4$ inherited from $\lambda_0$" is only approximate; the true stochastic
phase is $\sim54^\circ$, and closing the $54^\circ$-vs-$45^\circ$ gap is exactly Thread-4 frontier work (the $\lambda_0\to$
Borel bridge is unproven). (ii) The Borel plane has (at least) THREE singularities — a complex-conjugate pair AND the
real tail instanton — so the coefficient asymptotics are a competition, not a clean single pair; more coefficients
($v_6,v_7$) + Meijer-G would disentangle them and sharpen $|\zeta|,\theta$.

## 3bis. Borel–Padé on the 5 exact coefficients — complex pair CONFIRMED [NUMERIC, 2026-07-09]
> **⚠️ PARTLY SUPERSEDED by §00 (2026-07-19).** The complex-conjugate-pair CONCLUSION survives; the specific
> $\zeta\approx1.73\,e^{\pm i54°}$ does not — 7 coefficients give $\theta\approx50°$, $|\zeta|\approx1.8$–2.6.
The naive 2-term-recurrence extraction of $\theta$ was ill-conditioned (gave $\theta\approx0$). The *proper* Thread-1
tool (Borel–Padé$[2/2]$ on $b_n=v_n/n!$) gives a clean **complex-conjugate pair** of Borel singularities at
$$\zeta\approx1.73\,e^{\pm i\,54^\circ}\quad(\text{NOT real}).$$
So: (i) the oscillation / complex pair is confirmed by the right tool; (ii) the phase $\approx54^\circ$ is NEAR
$\lambda_0$'s $45^\circ$ (off $\sim9^\circ$ — within Padé$[2/2]$ error on 5 coefficients); (iii) dividing by $n!$ is
what makes it clean ⇒ **factorial divergence** confirmed (genuine resurgent series, $v_n\sim n!\,|\zeta|^{-n}\cos(n\theta
+\varphi)$). This SUPERSEDES the earlier geometric ($r^{-n}$) fit. Whether the true phase is exactly $45^\circ$ (=$\lambda_0$)
or genuinely $\sim54^\circ$ needs $v_5,v_6$ + Meijer-G.

## 3ter. Route evaluation (deep-research, 2026-07-09) — ranked
- **Thread 1 — locate the Borel singularities: CLOSABLE NOW (highest leverage).** Tool: **Meijer-G / hypergeometric
  approximants** (converge at 3–5 orders, native branch cut; outperform Borel-Padé) + Darboux/Dingle late-terms inversion;
  Borel-Padé as cross-check (done above → $\pm54^\circ$ pair). ~2–3 more exact coefficients decisive. Refs: Mera-Pedersen-Nikolic
  (arXiv:1802.06034), Crew-Trinh (arXiv:2208.07290). Caveat: off-axis (complex-conjugate) convergence not directly benchmarked;
  and the **real far-tail instanton $s^5/10$ may be a THIRD Borel singularity** that must be disentangled from the pair.
- **Thread 3 — reconstruct the REAL law from the complex pair: CLOSABLE NOW, no new coefficients.** Median / lateral Borel
  summation (Aniceto-Schiappa 2014, arXiv:1308.1115) + complex-saddle contributions to real observables (Basar-Dunne-Unsal
  2013). The $2\mathrm{Re}[\cdots]$ structure gives a real law from the conjugate pair.
- **Thread 2 — validate the renormalization: PARTIAL (targeted proof needed).** Canonical scheme is Varadhan subtract-the-
  divergent-mean (Hu-Nualart 2005) + an explicit Itô–Stratonovich counterterm $\tfrac12\tilde\phi\,b'Db$. Must PROVE
  "drop boundary-boundary self-contractions" $\equiv$ subtract-mean + counterterm; the $\theta(0)=0\Leftrightarrow$Itô
  shortcut that would have validated it directly did NOT survive verification. Open: does the engine's scheme silently drop
  a finite counterterm that shifts $v_3,v_4$?
- **Thread 4 — $\lambda_0\to$ Borel-singularity bridge: BLANK / multi-year frontier.** No mechanism transfers deterministic
  exact-WKB Stokes/Voros data to a noise-dressed series; the $\lambda_0$-at-$45^\circ$ / pair-at-$\pm54^\circ$ match is
  suggestive, unproven. Defer.

**Recommended next step:** get $v_5,v_6$ (sharpen $\pm54^\circ$ vs $45^\circ$, discriminate the conjugate-pair vs a third
real singularity) → Meijer-G/Darboux (Thread 1) → median resummation (Thread 3). Thread 2 is a parallel targeted proof.

## 4. Status & the next moves
- **Reachable now:** (i) push the engine to $v_5,v_6$ (predict $v_5\approx-0.9$) to nail the radius/period-8 oscillation
  and the Borel-singularity location; (ii) Padé–Borel the $\{v_n\}$ to exhibit the $\pm\pi/4$ singularity pair directly;
  (iii) independently validate the renormalization prescription (the drop-boundary-boundary scheme) — the one un-cross-checked
  ingredient in $v_3,v_4$.
- **The theorem-grade milestone:** a rigorously-defined **partial resurgent representation** — the coefficients as
  proved convergent (renormalized) chaos integrals + the instanton sector + the complex Borel pair from $\lambda_0$.
- **The frontier (multi-year):** the proved **stochastic exact-WKB** — dress the proved deterministic Weber resurgence
  (Nikolaev 2024) with the noise; derive the Stokes constant from $\lambda_0$; prove the trans-series equals $\mathcal W$.

## 5. Honest ledger
- **Solid/exact:** $v_0,v_1,v_2$ (three methods); the engine's reproduction of them; the renormalization-threshold
  algebra ($R_7=\tfrac15Y_1^5$, counterterm $\tfrac{23}3v_0^3$); the divergence verdict (from $|v_4|$'s jump).
- **Solid-modulo-prescription:** $v_3\approx-0.02$, $v_4\approx-0.45$ (renormalization scheme applied consistently but
  not independently cross-validated; the divergence conclusion is robust regardless).
- **Conjecture (evidence-backed but NOT yet independently pinned):** the $\pm\pi/4$ Borel-pair *phase* tied to
  $\lambda_0$. The data ROBUSTLY show (i) divergence (growing envelope, $|v_4|$ jump) and (ii) an oscillation
  (sign pattern $+,+,+,-,-$ with $v_3$ small = near a node). But $\theta=\pi/4$ is currently *assumed* (fixed to
  $\lambda_0$) in the fit — an independent 2-term-recurrence extraction of $\theta$ from the 5 coefficients is
  ill-conditioned (gives $\theta\approx0$ or non-real, destabilized by the small $v_3$), so the specific $45^\circ$
  is CONSISTENT with the data, not proven by it. Radius $\eta^2_c\approx0.5$–0.7. **Pinning $\theta$ independently
  needs $v_5,v_6,v_7$** (now cheap with the fast engine, modulo generating $Y_{10},Y_{11}$ and a 6-chaos build).
- **Retracted:** the sampled $v_3\approx-0.3$ (too noisy — the exact engine gives $-0.02$); the earlier "convergent"
  read (killed by the exact $v_4$). **[2026-07-19, from §00] Also retracted: $v_6>0$ / the sign flip to positive
  (exact $v_6=-1.90<0$); $\theta\approx54$–63° (7-coefficient fit gives $\approx50°$); $|\zeta|\approx1.2$
  (gives $1.8$–2.6). The "$\theta=\pi/4$ assumed" conjecture below is now TESTED: exact-45° is DISFAVORED, but
  so is $54°$ — the leading value is $\theta\approx50°$ (a genuine small rotation off λ₀). v₅ prediction
  $\approx-0.9$ (§4) was also off; exact $v_5=-1.19$.**
- **Blank (mapped):** rigorous stochastic exact-WKB.

**Honest ceiling (unchanged):** a literal 1-D closed form is provably impossible; the deliverable is the resurgent
trans-series. That structure is now on firm numerical footing (**7** exact coefficients $v_0..v_6$, divergent,
oscillatory) and the proof route is charted with a single well-identified blank. **[2026-07-19] The λ₀-tie is now
qualified: phase is NEAR-but-not-45° (~50°) and the modulus is strongly renormalized ($|\zeta|\approx1.9$ vs
$|\lambda_0|=1.26$) — Phase 3 must explain near-inheritance-with-modulus-shift, not clean inheritance. See §00.**
