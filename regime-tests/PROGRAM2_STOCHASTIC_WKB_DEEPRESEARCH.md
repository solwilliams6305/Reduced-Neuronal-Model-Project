# Program 2 — deep-research synthesis: the stochastic exact-WKB / Stokes-constant frontier (2026-07-10)

_Deep-research pass (105 agents, 23 primary sources, 25 claims adversarially 3-vote verified → 23 confirmed,
2 killed). Targets the load-bearing BLANK of `PROGRAM2_PROVED_CLOSED_FORM_ROUTE.md` Phase 3: exact-WKB /
resurgence / Stokes for the noise-dressed Weber operator, and a proved "stochastic Stokes constant" for 𝒲.
Full output: `tasks/wm6sl807z.output`. Tags [confirmed 3-0]/[medium]/[refuted]._

## Headline
The blank is REAL (no stochastic exact-WKB framework exists — confirmed), but the two halves needed to build one
are both mature and independently proved, and the pass identified **the specific recent machinery that supplies the
missing algebraic mechanism** (a noise-moved Borel singularity) plus **two concrete first moves**. This converts
"multi-year frontier" → "scoped attack with a named tool and a decisive first experiment."

## The seven verified findings
1. **[confirmed 3-0] The core gap is genuinely blank.** Random-operator spectral results attack via pure
   stochastic analysis (stochastic Prüfer, Volterra, BDG, Borel–Cantelli, Girsanov) and package noise as adapted
   Brownian integrals / Gaussian variance — never as trans-series / Voros / Stokes. Dong–Jian–Yuan 2026
   (arXiv:2606.22426, 1D white-noise Schrödinger high-energy asymptotics); Dumaz–Virág (arXiv:1102.4818, TW tail
   via Girsanov NOT exact-WKB).
2. **[confirmed 3-0] The deterministic Weber half is fully closed — with a NEW closed-form Stokes datum.**
   Nikolaev 2024 (arXiv:2410.17224, WKB resurgence proved on Riemann surfaces, meromorphic potentials). **Hao 2025
   (arXiv:2507.06922)** gives the closed-form **Weber Stokes datum**: disc $=-\tfrac12\log(1+e^{-2\pi i/\varepsilon})$,
   central charge $Z_{\gamma_f}=-2\pi i$, **deterministic Stokes constant $\Omega(\pm\gamma_f)=1$**, Borel poles at
   minus central charges of 4d BPS states, verified 8–15 digits. **This is the exact baseline the stochastic
   Stokes constant is a shift OF.** (New since the program's last research pass — a concrete deterministic anchor.)
3. **[medium, single source] THE BRIDGE — co-equational parametric resurgence (Aniceto–Crew 2024,
   arXiv:2410.13690).** Rigorously (systematic algebraic construction) encodes Borel singularities that **MOVE on a
   parameter-dependent algebraic curve $\Sigma_z=\{P_z(w,\varphi)=0\}$ and can CROSS the integration contour**, so
   exponentially-small corrections grow to dominate. **This is the exact algebraic mechanism a noise-shifted Borel
   phase (λ₀ 45° → stochastic 54–63°) requires.** Caveat: it treats a DETERMINISTIC external parameter $z$;
   promoting $z$ to a noise realization / making Voros symbols random variables on a fixed $\Sigma$ is the genuinely
   novel unproven step — no source has done it.
4. **[confirmed 3-0] Fluctuation-prefactor lesson (Bureković–Schäfer–Grauer, PRL 133 077202 2024, SNLS).** The
   Freidlin–Wentzell instanton action ALONE gives the WRONG tail; the Gaussian-fluctuation prefactor + zero mode
   dominates the quasi-singular tail (anomalous $a^{-1.88}$ the bare instanton can't produce). **Concrete precedent
   that the one-loop fluctuation sector, not the bare action, controls a noise-induced tail** — plausibly why the
   stochastic Weber Borel singularity sits off the classical-action value. (Refuted 1-2: a UNIVERSAL prefactor-
   dominance factorization; dominance holds specifically in the large-dissipation tail, not everywhere.)
5. **[confirmed 3-0] Route 3 (FP-PDE) is viable — and directly actionable given our new validated PDE.**
   Exact-WKB/Stokes CAN be applied to non-self-adjoint Fokker–Planck first-passage generators: Dao Duc–Schuss–
   Holcman (arXiv:1312.6992, 1405.7821) compute the **complex higher-order eigenvalues** of the non-self-adjoint FP
   operator in the small-noise limit (→ oscillatory survival probabilities) via asymptotic spectral analysis. Novel
   step: give those complex eigenvalues a resurgent/Borel structure. **This is the route that uses our just-built
   first-explosion PDE (`PROGRAM2_RICCATI_PDE_NOTES.md`).**
6. **[confirmed 3-0] Painlevé-IV ↔ Weber-kernel bridge (route 3 alternative, but a KNOWN-distinct object).**
   Xia–Xu–Zhao (arXiv:2301.05807): σ-form of Clarkson–McLeod PIV has a Fredholm-determinant rep via the
   **parabolic-cylinder (Weber) kernel**, governed by $D_{\alpha-1/2}^2$; van Spaendonck–Vonk (arXiv:2204.09062) +
   Sueishi–Kamata–Misumi–Ünsal (JHEP 12(2020)114, DDP dictionary) supply machinery to attach resurgent Stokes data
   to a Weber/PIV kernel by exact-WKB. **Caveat: the program already proved 𝒲 is NOT the standard PIV σ-form**
   (`STAGE1_GAP_DETERMINANT_NOTES`), so this is related-but-distinct — a check to run, not a solved reduction.
7. **[confirmed 3-0] Closest completed template = TW_β / stochastic Airy (q=1), and it exposes the gap exactly.**
   TW_β sharp tail PROVED via the stochastic operator + Girsanov (Dumaz–Virág), all-order tail only heuristic
   (Borot–Nadal loop equations, arXiv:1111.2761 — a claimed explicit all-order β-series was REFUTED 0-3) or two-term
   rigorous — **no proved resurgent trans-series for any random-operator edge law exists (𝒲 would be first-of-kind).**
   Deterministic Painlevé II (Hastings–McLeod) HAS a proved resurgent trans-series and Cleri–Dunne (arXiv:2002.06270)
   proved "integrability is not essential" — it survives dropping integrability. BUT that is the **drop-integrability
   axis, NOT the add-randomness axis** (rated 2-1): it motivates, does not support, the noise-dressing route.

## The single highest-leverage first move (from the synthesis + open questions)
**Promote the deterministic parameter of Aniceto–Crew co-equational resurgence to the noise coupling η (or a noise
realization), test whether the observed ~54–63° Borel phase falls on a computable co-equational singularity track
$\Sigma_z$ off the deterministic 45°.** If the phase shift lands on an algebraic curve, that IS the mechanism and
the "stochastic Stokes constant" becomes a computable object rather than a mystery.

**The minimal provable lemma to aim for:** averaging the random Voros symbol over the noise (Malliavin/Wick
expansion) leaves the deterministic Borel-singularity LOCATIONS fixed while shifting only the Stokes constant —
defining the stochastic Stokes constant rigorously as $\mathbb E[\text{random Voros datum}]$. (Open question 4.)

## Honest caveats (from the verifier)
- The bridge (Finding 3) rests on ONE source; "rigorous" = systematic construction w/ worked examples, not a
  theorem monograph; and it is deterministic-parameter, not random-field. The promotion-to-noise step is unproven.
- No source computed or proved a noise-induced Borel-phase shift; the 45°→54–63° shift is OUR numerical
  observation, uncorroborated in the literature (which is consistent with it being genuinely new).
- Robustness evidence (Cleri–Dunne) is on the wrong axis (drop-integrability, not add-randomness) — a 2-1 heuristic.
- Fast-moving frontier (all enabling pieces are Oct 2024 – Jun 2026); could be partially closed by unindexed 2026
  preprints.

## Attempted first moves (2026-07-10) — two cheap shortcuts RULED OUT [NEGATIVE, useful]
Before committing to the expensive routes, tested the two cheapest independent attacks on the 45°→54° puzzle.
Both fail — which confirms there is no free lunch; the 54° needs either exact coefficients (v6 grind, engine
built) or the frontier resonance theorem. Scripts `coupled-atlas/_fp_spectrum_probe.py`.
- **FP-generator complex spectrum (Finding-5 route): NEGATIVE.** Eigendecomposed the non-self-adjoint
  first-explosion generator (frozen-Y, from the validated PDE). Its complex eigenvalues sit at |arg|~2–10°
  (NOT ~54°) and coalesce to purely real above an exceptional point. The exceptional point IS real operator
  structure but is strongly Y-dependent ($\eta^2_{\rm EP}$: 0.11 at Y=−1.0 → 1.04 at Y=−2.19), overlapping the
  perturbative radius $\eta^2_c\!\approx\!0.6$–0.7 only near Y≈−1.8 — a curve $\eta_{\rm EP}(Y)$, not the single
  radius. Conclusion: the complex-pair Borel structure is NOT a local (frozen-slice) FP property; it is a GLOBAL
  connection-resonance (noise-dressed $\lambda_0$) object — which is the frontier itself (defining a resonance
  for a random oscillatory potential = the stochastic-exact-WKB blank; complex-scaling a white-noise realization
  is not well-defined).
- **PDE re-extraction of the Borel phase: DEFEATED at small η.** The deterministic PDE is superb at β=2 but
  extracting weak-noise $v_n$ needs Var$(\eta)$ at small η, where the escape peak (width ~η) is grid-resolution-
  limited: Var$/\eta^2$ converges only slowly (η=0.2: 0.187→0.152 as N: 1000→2800, target ~0.138). This is the
  SAME small-η wall the chaos engine was built to circumvent. Realistically cross-checks $v_0,v_1$ with heavy
  grid effort, cannot independently pin the phase.

**Net of the probes:** the 54° cannot be reached cheaply. The decisive moves remain (a) exact $v_6,v_7$ from the
transfer engine (built; slow), or (b) the noise-shifted global resonance (the frontier theorem). Neither is a
quick experiment — the deep-research verdict that this is a genuine multi-year frontier is empirically reconfirmed.

## Net
Two independently-proved halves — deterministic Weber exact-WKB (resurgent, closed-form Stokes datum, Hao 2025) and
the stochastic-operator edge law (rigorous tail, but resurgence-blind, Dumaz–Virág) — have **never been joined**.
The join has a candidate mechanism (co-equational moving Borel singularities), a candidate concrete route through
our validated first-explosion FP-PDE (non-self-adjoint complex spectrum → resurgence), and a candidate minimal
lemma (noise-averaged Voros datum). The frontier is now scoped.
