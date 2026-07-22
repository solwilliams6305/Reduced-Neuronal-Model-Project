# Channel (state-dependent, multiplicative) noise — first probe: robust CLASS, but the law is a FUNCTIONAL of the noise profile (not a frozen β_eff), plus a degenerate fork

_July 2026 (incoming agent, Fable 5). Cheap decisive test of the channel-noise extension: replace the additive
constant inner noise $\eta$ in the trusted `fp_cusp` escape solver with a **state-dependent** $\eta(Y)$ (channel
noise: multiplicative, $\propto\sqrt{\text{gating activity}}$, $\sigma\sim1/\sqrt N$). **Result CORRECTS the
prior prediction** ("easy corollary: localized at the escape node → frozen $\beta_{\rm eff}$"): the cusp law is a
**functional of the noise profile over the whole approach** (turning + confining side both matter, landing point
matters least), the **class is robust** (skew+/exk− cusp fingerprint survives mild state-dependence), but the
**effective β is shifted** (~15–30% cumulant change for a realistic ~10% profile variation) — not a single
frozen value. And a **degenerate fork** (noise vanishing near the turning) leaves the β-family (exk → +3).
Scripts `coupled-atlas/channel_noise_probe.py`, `channel_noise_probe2.py`. Tags **[NUMERIC probe]/[prediction
CORRECTED]**. Preliminary (FP-based, one profile family; the degenerate exk needs finer numerics)._

---

## Setup
Inner Riccati $dp=(\operatorname{sign}(Y)Y^2-p^2)\,d\tau+\eta(Y)\,dW$, $Y=Y_0-\tau$; $\eta$ depends on the sweep
$Y$ (leading channel-noise model: intensity set by the deterministic passage $v\approx v(Y)$; no Itô/Strat drift
since $\eta$ is not $p$-dependent). FP with $D(Y)=\eta(Y)^2/2$. Reference constant-η β-family (FP): β=1/2/4/8 →
skew +0.26/+0.60/+0.85/+0.85, exk −0.78/−0.24/+0.64/+1.20.

## Findings

**(1) The prior "localized at the escape node" claim is WRONG.** Changing $\eta$ *inside* the landing region
[−2.5,−1] (case D) barely moved the law (mean −1.599 vs −1.594); changing it *outside* (cases B,C, which include
the turning) moved it drastically (mean → −0.76, skew even flipping). Corrected localization test: a ×1.5 bump in
the **turning** [−0.5,0.5] shifts skew −0.29 and exk to −1.0; a bump in the **confining tail** [1,3] shifts skew
+0.24; a bump in the **landing** [−2.2,−1.2] shifts skew only −0.14. So sensitivity is **distributed over the
approach (turning + confining side), NOT localized at $Y^\star$** — consistent with the Malliavin sensitivity
kernel $\Phi(s)$ being spread over $[Y^\star,Y_0]$, peaking near the turning, not at the landing. **[NUMERIC ✓]**

**(2) The CLASS is robust to realistic (mild) state-dependence.** A $\pm10\%$ profile variation across the
turning ($\eta=\sqrt2(1+0.1\tanh Y)$, the scale of real channel noise's $O(\sqrt\varepsilon)$ variation at
$\varepsilon\sim0.01$) keeps the cusp fingerprint: skew +0.70 (vs +0.60), exk −0.16 (vs −0.24) — same signs, same
family. So the $\mathcal W$-universality survives channel noise qualitatively. **[NUMERIC ✓]**

**(3) The effective β is a FUNCTIONAL of the profile — and it's DERIVED.** By Itô isometry (Tier A §3.4) the
escape variance with state-dependent noise is $\operatorname{Var}(Y^\star)=\int\eta(s)^2\Phi(s)^2\,ds$, $\Phi$
the sensitivity kernel. Hence the effective noise is the **$\Phi^2$-weighted average** and
$$
\boxed{\ \beta_{\rm eff}=\frac{4\int\Phi(s)^2\,ds}{\int\eta(s)^2\Phi(s)^2\,ds},\qquad \eta(s)=\sigma\,g(v(s)).\ }
$$
This *is* the "functional, not frozen value," and it is now **VERIFIED to 3 digits at weak noise**
(`channel_noise_beta_eff_weak.py`, converged perturbative $Y_1$): for $\eta=\sqrt2(1+a\tanh Y)$ the variance ratio
$\operatorname{Var}(Y_1,{\rm prof})/\operatorname{Var}(Y_1,{\rm const})$ = **0.873** ($a{=}0.10$) and **0.813**
($a{=}0.15$), matching the $u_c^4$-weighted predictions **0.873 / 0.813** exactly. So the weak-noise escape
variance *is* $\int\eta^2\Phi^2$ with $\Phi^2\propto u_c^4$.

**Two corrections this forces (honest):** (i) The variance kernel $u_c^4$ actually peaks on the **oscillatory /
landing side** ($s\approx-1.1$, 82% of the weight on $s<0$) — **NOT the turning**. My earlier "turning-weighted"
statement was from the *skew* shifts, which have a *different* kernel; the *variance* is landing-weighted. So the
"localization" is **cumulant-dependent**. (ii) The $\beta_{\rm eff}$ formula is a **weak-noise** statement. At the
physical $\beta=2$ (strong noise, where $\operatorname{Var}/\eta^2$ is already 74% above $C_V$), it does **not**
apply: the FP profile tests there show the sensitivity migrated to the **turning**, and the law deforms **off the
additive β-family** (std↑ *and* skew↑ are incompatible with any single β). So the physical-regime channel-noise
effect is a genuinely **nonlinear bulk deformation**, beyond the weak-noise $\beta_{\rm eff}$. Biophysically
$\beta_{\rm eff}\sim\varepsilon N$ *scaling* holds; the constant is the $u_c^4$-weighted gating profile at weak
noise, with a nonlinear correction at $\beta=O(1)$. **[DERIVED + VERIFIED weak-noise 3-digit; strong-noise
deformation NUMERIC]**

**(4) The CLASS (tail exponent) is ROBUST — including the degenerate case; the degenerate exk is a BULK effect
only.** The non-confounded class marker is the **left-tail exponent** (=2q+1=5, the FW instanton over
$\operatorname{sign}(Y)Y^2$). Deep-FP measurement gives **5.22 (const), 5.28 (mild ±10%), 5.22 (degenerate
$\eta\!\to\!0$ at the turning)** — all ≈5. So the FW tail is **potential-driven and invariant to the noise
profile**, even when the noise vanishes at the turning. This **corrects the earlier "degenerate fork = new
class"**: the degenerate case's exk→+3.1 is a **bulk/shape deformation, NOT a tail-class change** — the cusp
universality class survives channel noise *in every case tested*. **[NUMERIC ✓, non-confounded]**

## Verdict (final for this probe)
**The cusp universality CLASS is robust to channel noise; the noise only reshapes the BULK fingerprint through a
$\Phi^2$-weighted effective $\beta$.** Concretely:
- **Class marker (left-tail exponent 5) is invariant** to the noise profile — mild *and* degenerate — because it
  is the FW instanton over $\operatorname{sign}(Y)Y^2$, a *potential*-driven quantity. So $\mathcal W_\beta$ is
  genuinely universal against state-dependent/multiplicative channel noise. This is the strongest and cleanest
  outcome, and it's non-confounded (tail exponent, not a cumulant).
- **Bulk fingerprint (cumulants / effective $\beta$) is a derived functional** $\beta_{\rm eff}=4\int\Phi^2/\int
  \eta^2\Phi^2$, dominated by the gating profile *at the turning* (where $\Phi^2$ peaks) — not the frozen
  $g(v_{\rm cusp})$ I first predicted, and not the landing point.
- **Two prior predictions corrected:** (i) "localized at $Y^\star$ → frozen $\beta_{\rm eff}$" — *wrong*
  (turning-weighted functional); (ii) "degenerate fork = new class" — *wrong at the class level* (the tail
  exponent stays 5; only the bulk exk leaves the additive range).
- **Difficulty (revised):** *moderate but tractable* — the class-robustness is now established (tail exponent), and
  $\beta_{\rm eff}$ is a Tier-A Itô-isometry corollary; Tier A/B extend with $\eta(Y)$ bounded+floored. The
  biophysical $\beta\sim\varepsilon N$ scaling holds with a $\Phi^2$-weighted constant.

## Prior art + capstone (deep-research verdict, 110 agents, verified)

**Our result is a NOVEL instance of a KNOWN general principle** — which is the reassuring outcome (we're on
solid ground *and* contributing something new).

- **The principle is classical.** The large-deviation/Eyring–Kramers exit-rate *exponent* is set by the
  deterministic quasipotential, and state-dependent/multiplicative noise enters only through the prefactor +
  an $O(\sigma^2)$ barrier correction (Bouchet–Reygner 2015, arXiv:1507.02104; Yang–Potter–Cameron/Dahiya 2018,
  1806.05321; **Moreno–Barci–González-Arenas, PRE 101 062110, 2020**, 1912.00514). So our tail-exponent
  invariance is the expected FW behavior. Independently, β-edge/Tracy–Widom universality is robust to the
  confining potential and entry distribution, with all β-dependence in one scalar $2/\sqrt\beta$ (RRV 2011; KRV
  2016; Liu–Zou 2025, 2508.17838) — the RMT analogue of "edge law depends on the noise through one parameter."
- **Our β_eff has a direct precedent to bridge to.** Moreno et al. show the multiplicative-noise rate reshaping
  is captured by *local* noise-function values at the escape extrema (their two-point $|g_0|^{2\alpha}
  |g_a|^{2(1-\alpha)}$ structure) — the discrete, rate-level analogue of our continuum edge-law functional
  $\beta_{\rm eff}=4\int\Phi^2/\int\eta^2\Phi^2$. Deriving ours as the continuum limit of theirs would rigorously
  connect the two (a clean next theorem).
- **But the specific combined statement is NOVEL [BLANK]:** no source proves *degenerate* multiplicative-noise
  robustness of a *Weber edge law* $\mathcal W_\beta$ with the *tail exponent pinned at 5* and the *explicit
  $\Phi^2$-weighted $\beta_{\rm eff}$*. Every analogue falls short on ≥1 axis (non-degenerate diffusion; exit-RATE
  not edge-LAW; potential/entry universality not multiplicative-SDE universality). Bouchet–Reygner explicitly
  assume a **nondegenerate** diffusion matrix — so **our degenerate case (noise vanishing at the turning) is
  genuinely open and is where the novelty concentrates.**

**Our findings line up precisely with the literature's caveats** (a good consistency check):
- Moreno et al.'s exponent-invariance holds only for $\sigma^2\ll\Delta U$ (weak noise); above a threshold noise
  destroys the barrier — **this is exactly our weak-noise/strong-noise split** (the β=2 off-family deformation is
  the strong-noise regime where our $\beta_{\rm eff}$ fails).
- Their $O(\sigma^2)$ barrier correction is **Itô/Stratonovich-prescription-dependent** ($\alpha$; zero only for
  anti-Itô) — **this is precisely the noise-induced drift** I flagged for B1d; it's real and convention-dependent.

**Capstone (coupled Morris–Lecar + channel noise) — method exists, target unbuilt:**
- The route is published: WKB / most-probable-escape-path large deviations through the slow-fast conductance
  model, with concrete numerics (gMAM for MPPs, ordered-upwind for the density) — **Newby–Bressloff–Keener 2013**
  (1304.6952), **Newby 2014** (1406.2914).
- **Cautionary + validating:** the frozen-slow-variable Kramers reduction *fails* for channel-noise escape in
  slow-fast conductance models — one must use the full large-deviation/instanton formulation. **This validates
  our instanton approach** and warns against naive reductions.
- **Model choice:** for the escape *tail*, use a stochastic-hybrid/**PDMP** channel-noise model, not the Fox–Lu
  diffusion approximation (flagged inaccurate for metastable/escape dynamics). Open: does tail-exponent-5 survive
  the PDMP formulation?
- **[BLANK] the coupled-ML cusp is unbuilt:** no source exhibits a cusp (two coalescing folds / synchrony-loss
  singularity) in single or coupled stochastic Morris–Lecar. So the capstone is an **open construction**: (a) find
  or engineer a cusp in coupled ML's synchrony mode (does its SNIC/Hopf/fold structure even admit one?), (b) PDMP
  channel noise + instanton, (c) verify the $\mathcal W_\beta$ / tail-5 signature.

## The three follow-up pieces — RESULTS (all decisive)

**Piece 1 — the bridge (channel noise enters only through the noise intensity along the escape path).** Unifies
bulk + tail:
- *Bulk* (weak-noise variance, the body of $\mathcal W_\beta$): $\operatorname{Var}(Y^\star)=\int\eta^2\Phi^2$,
  $\Phi^2\propto u_c^4$ (landing-weighted) — derived + 3-digit verified (above).
- *Tail* (deep escape, the FW instanton): with state-dependent noise the action is $\tfrac12\int\phi^2/\eta(Y)^2
  \,d\tau$, so the noise enters as $1/\eta^2$ **along the instanton path**; the exponent stays 5 (potential-driven,
  path unchanged) and the *rate* is set by $\eta$ at the escape depth. **Confirmed** (`channel_noise_pieces12.py`):
  a profile with $\eta{\times}1.4$ for $|Y|>2.5$ (deeper noise) *lowers* the deep-tail $-\log P$ — ratio at
  $s{=}5,5.5,6$ = 0.61, 0.59, 0.58, trending to the predicted $(\eta_{\rm shallow}/\eta_{\rm deep})^2=1/1.4^2=0.51$
  as the instanton goes deeper. **The bridge to Moreno et al.:** their multiplicative-noise Kramers rate depends on
  $g$ at two points (well + barrier, $|g_0|^{2\alpha}|g_a|^{2(1-\alpha)}$); ours is the *swept/edge-law* version —
  $\eta$ along the most-probable escape path (bulk kernel $\Phi^2=u_c^4$; tail instanton), with the exponent
  invariant. Same principle, generalized from static Kramers to the swept first-passage edge law. **[DERIVED +
  NUMERIC ✓]**

**Piece 2 — the strong-noise deformation is a 2-parameter family; the 2nd parameter is the ODD moment of $\eta(Y)$.**
At β=2, testing odd vs even $O(0.15)$ profiles (`channel_noise_pieces12.py`), against the β-family curve
(std,skew): (.90,.26)(.69,.60)(.48,.85):
- **ODD** $\eta=\sqrt2(1\pm0.15\tanh Y)$: (0.730, **+0.745**) and (0.658, **+0.435**) — skew shifted $\pm0.15$
  *off* the curve (at those std's the family gives skew ≈ 0.53 / 0.65). **OFF-family, both directions.**
- **EVEN** $\eta=\sqrt2(1+0.15\,g_{\rm even})$: (0.696, 0.626) and (0.678, 0.639) — essentially *on* the curve.
  **In-family (just a rescaled $\beta_{\rm eff}$).**
So the **even part of the noise profile rescales $\beta_{\rm eff}$ (stays in-family); the odd part drives the
deformation off the additive family** — it couples to the odd potential $\operatorname{sign}(Y)Y^2$. The
channel-noise law at $\beta=O(1)$ is thus (to leading order) a **2-parameter family: $(\beta_{\rm eff},\,
m_{\rm odd}[\eta])$**, $m_{\rm odd}$ the odd moment of the gating profile. This is the genuinely-new nonlinear
content. **[NUMERIC ✓, decisive]**

**Piece 3 — coupled Morris–Lecar DOES admit an antisymmetric-mode cusp (the capstone target exists).** The
research flagged this as an unbuilt blank; it's not. Structural argument: the FHN cusp needs only *swap symmetry +
a fold + the cubic term*, none FHN-specific; coupled ML has the swap symmetry and (relaxation regime) an N-shaped
fast nullcline with a fold, so the same mechanism applies. Confirmed numerically (`morris_lecar_cusp.py`, standard
ML, $I=90$): the $v$-nullcline is N-shaped (folds at $v=-28.4, 10.6$); the antisymmetric linear coefficient
$\mu(w)=\partial_v f(v_s(w),w)-2g_c$ has $\partial_v f$ with an interior extremum $0.4185$ → **cusp at
$g_{\rm crit}=0.209$**; and the fold separation obeys $\Delta\sim\sqrt{g_{\rm crit}-g_c}$ ($\Delta/\sqrt{\Delta g}$
= 61.7, 61.7, 62.5, 65.9 over a 64× range — constant to ~7%) — the **same √ fold-coalescence signature as FHN's
$\Delta=2\sqrt{-2g/3}$**. So the coupled-cusp/$\mathcal W_\beta$ scenario is realizable in a genuine conductance
model. **[DERIVED-structural + NUMERIC ✓]** *(Remaining for the full capstone: the noise-induced escape at this
cusp with PDMP channel noise → verify the $\mathcal W_\beta$/tail-5 fingerprint; $g_{\rm crit}=0.209>0$ here, so
the regime/orientation needs the noisy sim — the geometric prerequisite is now established.)*

## The noisy ML capstone — ATTEMPTED, INCONCLUSIVE, precisely diagnosed [honest negative]

I ran the noisy coupled-ML sweep (`morris_lecar_capstone_attempt.py`: 2 units, diffusive $v$-coupling, additive
noise, sweep $g_c$ through the cusp, measure per-cycle $\max|v_-|$ desync amplitudes). **It does NOT cleanly
confirm the $\mathcal W_\beta$ fingerprint**, and the failure is diagnostic:
- **The exk trend is not robust — it reverses with noise level.** At $\sigma=3$ the desync-amplitude exk *falls*
  ($+2.0\to-0.2$) as $g_c$ weakens; at $\sigma=1.5$ it *rises* ($+0.9\to+3.1$). Opposite trends ⇒ the observable
  is confounded. The only clean signal is the spread amplification (~×3.6, matching FHN's ×3.4–3.9), but that is
  *generic* (weaker coupling → more desync → more spread), not cusp-specific.
- **Diagnosis (the real lesson):** the per-cycle desync **amplitude** is the wrong observable — the $\mathcal
  W_\beta$ fingerprint lives in the escape **location/timing** distribution (where/when the peel-off occurs), not
  its amplitude. Rare noise-induced peel-offs make the amplitude heavy-tailed (leptokurtic, exk>0), which is not
  the sub-Gaussian escape-location signature. Additionally the linearized $g_{\rm crit}=0.209$ (a fast-subsystem
  cusp) is **not** the dynamical (Floquet) synchrony-loss coupling, so the sweep wasn't centered on the transition.

**What the proper capstone needs (mapped):** (i) the escape-**location** observable matching the FHN
physical-signature protocol (`PHYSICAL_SIGNATURE_NOTES.md`, `crossover_fold_to_cusp.py`), not desync amplitude;
(ii) the *dynamical* synchrony-loss coupling (Floquet analysis of the coupled limit cycle), not the linearized
cusp; (iii) contrast against a genuine *fold* regime (exk cusp < exk fold); (iv) ideally PDMP channel noise.
So the capstone remains **open** — the geometric cusp is solid (§Piece 3), but demonstrating the noise-induced
$\mathcal W_\beta$ signature in full coupled ML is a careful protocol, not a quick sweep. My quick version is
retracted as inconclusive.

## Remaining (lower priority)
- Full ML capstone via the correct escape-location protocol (mapped above) — the real biophysical deliverable.
- Rigorous derivation of the 2-parameter $(\beta_{\rm eff}, m_{\rm odd})$ family (Piece 2 is numeric).
