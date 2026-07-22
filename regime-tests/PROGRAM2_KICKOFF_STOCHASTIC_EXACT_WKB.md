# Kickoff — Program 2: the resurgent trans-series / stochastic exact-WKB of 𝒲 (for a fresh chat + /deep-research)

_Paste this into a new chat, then run `/deep-research` with the research questions in §6. This is the second of
the two multi-year cores of the coupled-FHN cusp program (the first — Tier B universality-with-rate — is ~75%
done). Program 2 is ~40% done. **UPDATE 2026-07-09 — see `PROGRAM2_CONSOLIDATION.md` §0 (FINAL STATE):** the first
~six weak-noise coefficients are computed by an exact deterministic Wiener-chaos engine (`chaos_diagram.py`) and
**cross-validated against the ground-truth direct Itô MC**: $v=(0.134,0.111,0.104,-0.03,-0.45,\sim{-}1,\,v_6{>}0)$.
The series is **FACTORIALLY DIVERGENT (genuine resurgent trans-series)**; Borel plane = a **complex-conjugate pair at
$|\zeta|\approx1.2$, phase $\theta\approx54$–63°** PLUS a competing **real tail-instanton** ($e^{-s^5/10\eta^2}$).
**Correction:** $\theta$ is robustly **NOT** $\lambda_0$'s $45^\circ$ — the "$\lambda_0$-inherited phase" surmise is
refuted; the $\lambda_0\to$Borel bridge is a genuine frontier blank. The renormalization + $\tfrac12$ boundary
convention are validated against the Itô law. Remaining: precise $v_6,v_7$ (transfer-engine rewrite) to pin $\theta$;
and the PROVED trans-series (rigorous stochastic exact-WKB — Nikolaev 2024 gives the deterministic backbone). Repo:
`Reduced Neuronal Model Project/`; trackers `coupled-atlas/PROGRESS.md`, notes in `regime-tests/`._

---

## 1. The object 𝒲

Noise-induced escape at a **cusp** in a coupled FitzHugh–Nagumo system has a universal edge law $\mathcal W_\beta$
(the cusp analogue of Tracy–Widom, $\beta=4/\eta^2$ the Dyson index). It is the **first-explosion law of the
stochastic Weber Riccati** $dp=(\operatorname{sign}(Y)Y^2-p^2)\,d(-Y)+\tfrac2{\sqrt\beta}dW$, $p\sim+|Y|$
recessive — equivalently the first-node law of the **stochastic Weber field** $u''=(\operatorname{sign}(Y)Y^2-
\eta\dot W)u$. At $\beta=2$ its fingerprint is skew $+0.607$, exk $-0.237$; tails are $(2q{+}1,\,3q/2)=(5,3)$.
Tier A (a rigorous characterization: well-posed, moment-determined β-family, backward-Kolmogorov PDE) is closed;
`fp_cusp.py` is the trusted MC-free ground truth (β=2: skew +0.601, exk −0.244).

## 2. Why the closed form MUST be a resurgent trans-series (provable obstructions)

A literal 1-D closed form is *provably* unavailable — these are theorems, not gaps, so **do not pursue them**:
- **Not a Painlevé-σ-ODE.** The deterministic skeleton is an **isomonodromy FIXED POINT** (the only relevant
  PC-class deformation, the linear $\Delta Y$ term, flows cusp→fold), so there is no Painlevé τ-flow to reduce
  onto. (`T2_CONNECTION_DATA_DERIVATION`, `CONNECTION_CLOSED_FORM_NOTES`.)
- **Not a Fredholm/soft-edge determinant.** $\mathcal W$ is a first-passage law of an **unbounded-below**
  operator (escape down the $-Y^2$ side, no bottom eigenvalue) — not an eigenvalue-gap of a determinantal
  process. (`STAGE1_GAP_DETERMINANT_NOTES`.)
- The τ-function / RH-with-flow route hits the same wall. A clean "no determinantal form" negative is itself a
  theorem worth stating.

So the **only viable "closed form" is a resurgent trans-series** — it embraces the irreducibility. The
deterministic connection is already closed-form (both Weber Γ-factors, machine-verified: $e^{3i\pi/4}\Gamma(3/4-
i\lambda/4)/\Gamma(1/4-i\lambda/4)=\Gamma(3/4-\lambda/4)/\Gamma(1/4-\lambda/4)$, root $\lambda_0=0.890-0.890i$).

## 3. What's DONE (Step 0 + structure)

- **Leading weak-noise coefficients pinned** by three agreeing converged methods (perturbative expansion off the
  exact Weber backbone; small-η MC; closed-form Green's function): mean $\langle Y^\star\rangle=-2.188+0.212\,
  \eta^2+\cdots$; variance $=\eta^2(0.1337+0.110\,\eta^2+\cdots)$; skew $=1.16\,\eta+\cdots$. **The old FP
  small-η values (C_V=0.20, skew-coeff 0.66) were RETRACTED as boundary-layer artifacts**, along with the
  FP-based "resurgent-growth ratios 1.23→2.33" (the FP had the wrong *sign* for $v_1$). (`PROGRAM2_ROUTE2B_NOTES`,
  `weaknoise_v1.py`, `weaknoise_sectors.py`.)
- **Instanton (non-perturbative) sector:** the left tail is the Freidlin–Wentzell instanton $-\log\mathbb P
  (Y^\star<-s)\to s^5/(10\eta^2)=\beta s^5/40$ — the Borel-plane singularity sits at the instanton action
  $S=s^5/10$. (`PERSISTENCE_ITEM2_NOTES`, `instanton_action.py`; and now *rigorously* upper-bounded at exponent
  5, `PROGRAM1_TIERA_THEOREM.md` §4.2.)
- **Resurgence signature (robust part):** at $\beta=2$ ($\eta=\sqrt2$) the leading weak-noise skew extrapolates
  to $1.16\times1.414=1.64$ vs the true $0.605$ — a **×2.7 overshoot** ⇒ β=2 is non-perturbative, the series
  needs resummation, and the perturbative radius is $\eta^2\approx v_0/v_1\approx1.2$ (β≈3.3).

## 4. The FIRST TACKLE (this session): v₂ reframed from wall to a finite Wick computation

The perturbative recursion $Y^\star=Y^\star_0+\eta Y_1+\eta^2Y_2+\cdots$ (node-shift about the deterministic node
$Y^\star_0=-2.188$) is **noise-safe through $Y_3$** but at 4th order needs $u_1'''(Y^\star_0)=V'^\star u_1+V^\star
u_1'+\xi(Y^\star_0)\,u_0'^\star$ — the **white noise at the fixed node**. Reframed:
- **It is a limitation of the recursion method, not of the law** (the field's first node always exists).
- **$\xi(Y^\star_0)$ is pathwise-singular but finite in the moments:** in $\langle\cdots\rangle$ it Wick-contracts
  with the lower-$Y_k$ noise, $\langle\xi(Y^\star_0)\xi(s)\rangle=\delta(Y^\star_0-s)$ **δ-localizes** to a finite
  **boundary value**. So $v_2$ is a *finite* boundary-Wick computation.
- **Two routes:** (a) **boundary-Wick** — treat $\xi(Y^\star_0)$ as an Itô/Malliavin object, collect the
  δ-localized boundary terms (the clean analytic path = the stochastic-exact-WKB step; generalizes the
  deterministic Weber Voros/exact-WKB machinery); (b) **regularized full-field** — mollify the noise, node-cumulants,
  extrapolate $\delta\to0$ (validation, precision-limited for the $O(\eta^6)$ coefficient).

So the "named hard problem" is now a **well-posed subproblem**, not a distributional wall — this is the on-ramp.

## 5. The target milestone (reachable) vs the far frontier

- **Reachable (publishable): a PARTIAL resurgent representation** — a few converged perturbative orders (via §4's
  boundary-Wick) + the instanton sector ($S=s^5/10$) + the leading resurgence relation / Stokes constant, i.e.
  "the defining resurgent structure of $\mathcal W$."
- **Far frontier (multi-year): full stochastic exact-WKB** — all-orders Voros symbols for the *stochastic* Weber
  operator, the Borel transform, the Stokes automorphism connecting sectors. Deterministic Weber Voros
  coefficients are known (Iwaki–Koike–Takei, Bernoulli numbers); the stochastic generalization is new.

## 6. RESEARCH QUESTIONS for the /deep-research pass

_Applied-math / mathematical-physics — real theorems, authors/titles/years, "known/partial/blank" verdicts; not
a general survey. Frame: does the machinery to compute the perturbative sectors + Borel/Stokes structure of a
noise-induced escape law exist, and what is the single most transferable technique?_

1. **Exact WKB / Voros symbols for the (deterministic) Weber / parabolic-cylinder equation** — Iwaki–Koike–Takei
   "Voros coefficients for the Weber equation" and the Bernoulli-number closed form; the Borel summability and
   Stokes structure of the Weber Voros coefficients. What exactly is known, and in what form?
2. **STOCHASTIC / random exact-WKB or resurgence** — is there ANY exact-WKB, Voros-symbol, or Borel–resurgence
   framework for a Schrödinger operator with a RANDOM (white-noise) potential, or for an SDE first-passage/
   explosion functional? (Likely a blank — that itself is informative.) Adjacent: resurgence in stochastic PDEs,
   Borel summation of weak-noise / small-η loop expansions, instanton–anti-instanton in Langevin escape.
3. **Resurgence of FIRST-PASSAGE / large-deviation functionals** — Borel-resummed weak-noise expansions of exit
   distributions / Kramers rates where the perturbative series is asymptotic and the Borel singularity sits at
   the instanton action; trans-series for escape-time/location laws (Freidlin–Wentzell prefactor expansions,
   Bender–Wu-type, ecalle resurgence applied to path-integral instantons).
4. **The boundary-Wick / Malliavin computation of perturbative escape coefficients** — is the "white-noise-at-the-
   node" term (§4) a known object? Malliavin calculus for first-passage functionals with pointwise-noise boundary
   contributions; the local-time / δ-localization structure; any method that computes the $O(\eta^{2n})$ cumulant
   coefficients of an SDE first-explosion law to high order.
5. **Borel plane + Stokes constant from the instanton action** — given the instanton action $S=s^5/10$ (β s^5/40)
   and leading coefficients $v_0=0.134$, $v_1=0.110$, what does resurgence predict for the large-order growth
   $v_n\sim n!/S^n$ and the leading Stokes constant? Any worked example connecting an $O(1)$ instanton action to
   the coefficient asymptotics of an escape law.
6. **Worked analogues** — any transcendent with NO ODE reduction whose "closed form" was given as a resurgent
   trans-series (a precedent for the deliverable), especially in an isomonodromy-fixed-point / non-integrable
   setting.

**Deliverable from the pass:** for each thread, is the machinery available (cite) or a blank; the single most
transferable technique to (i) compute the perturbative sectors beyond $v_1$ via boundary-Wick, and (ii) assemble
the partial resurgent representation (few orders + instanton + Stokes constant). Flag genuine blanks.

## 7. Practical (for the new chat)

- Trusted tools: `fp_cusp.py` (MC-free β=2 ground truth), `weaknoise_sectors.py`/`weaknoise_v1.py` (converged
  perturbative coefficients — the machinery that hits the $v_2$ wall), `instanton_action.py` (left-tail instanton,
  const 1/10), `connection_closed_form.py`/`resonance_gamma.py` (closed-form connection, machine-verified). scipy
  available (`pbdv`, `pbwa`, `loggamma`, `solve_bvp`).
- Program culture (load-bearing): tag everything proved / cited / numerically-validated / conjectural; never
  dress a numerical check as a proof; cheap decisive checks vs `fp_cusp` every time; prefer non-confounded
  signatures (coefficient growth, Borel-singularity location, β-independence) over bulk cumulant-fitting; retract
  overclaims openly; a precise obstruction + a map is an acceptable outcome. Honest ceiling: a literal 1-D-ODE
  closed form is provably impossible — aim for the **partial resurgent representation**, frame full
  stochastic-exact-WKB as the mapped frontier.
