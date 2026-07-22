# Detailed scoping of the two top-tier open programs

_July 2026 (incoming agent, Fable 5). Research-map-level scoping of the two programs that would lift the 𝒲 work
to top-tier: (1) a rigorous operator-limit / universality theorem; (2) the exact integrable identity (closed
form). For each: target theorem, proof architecture, what's in hand vs open, the precise obstruction, minimal
publishable path, difficulty. Grounded in this session's results (`CONNECTION_CLOSED_FORM_NOTES`,
`FRONTIER_PROGRAMS_NOTES`, `PERSISTENCE_ITEM2_NOTES`, T1). Script `coupled-atlas/transseries_orders.py`._

---

## PROGRAM 1 — Rigorous operator-limit / universality theorem

Two tiers; Tier A is tractable and self-contained, Tier B is the hard universality statement.

### Tier A — characterization (the tractable half; ~1 paper)
**Target.** $\mathcal W_\beta$ is rigorously the first-explosion law of the stochastic Weber Riccati
$dp=(\operatorname{sign}(Y)Y^2-p^2)\,d(-Y)+\tfrac2{\sqrt\beta}dW$, $p\sim+|Y|$ recessive; it is a well-defined,
moment-determined, monotone $\beta$-family whose CDF solves the backward-Kolmogorov PDE.

**Architecture** (status): **A1** SDE well-posed up to explosion [standard]; **A2** explosion a.s. finite via
comparison $\dot p\le M-p^2$ [in hand; numerically fraction $=1.00000$]; **A3** explosion-law = Feynman–Kac
solution of the backward-FP PDE [validated by `fp_cusp`]; **A4** moment-determinacy from rigorous two-sided tail
bounds (left exp 5, right exp 3, both $>2$ ⇒ Carleman) [**the crux** — needs the FW tails made rigorous; the T1
confined-Bernstein bound is the tool for the right/barrier tail, the left/persistence tail needs the instanton
made rigorous]; **A5** $\beta$-monotonicity via coupling [doable].
**Verdict:** provable modulo **A4** (rigorous tails), which *shares the T1 estimate machinery*. Highest-ROI next
theorem.

### Tier B — universality/convergence (the hard half; multi-paper)
**Target.** The rescaled coupled-FHN antisymmetric-mode escape at the cusp converges in law to $\mathcal W_\beta$
($\beta=4/\eta^2$) as $\varepsilon\to0$, **with a rate $\varepsilon^{p}$**.

**Architecture** (status): **B1** stochastic reduction FHN → antisymmetric swept SDE near the cusp, error uniform
through the merge [GSPT/blow-up; deterministic backbone Kristiansen–Pedersen; stochastic error = STOCHASTIC_BLOWUP
work]; **B2** uniform sub-Gaussian tube [**T1 — closed this session** modulo the re-derived Malliavin regularity];
**B3** upgrade the tube (moment bounds) to convergence-in-distribution with a rate [**NEW, HARD** — needs a
quantitative coupling / Stein–Lindeberg argument]; **B4** show the outer (finite-$\varepsilon$) corrections vanish
at rate $\varepsilon^p$ [**the crux** — this is exactly the contamination that makes the physical fingerprint only
*qualitative* at finite $\varepsilon$, per `PHYSICAL_SIGNATURE_NOTES`].

**The precise obstruction (Tier B):** a *quantitative rate* for a singular-perturbation-with-multiplicative-noise
limit **through a cusp**. RRV is an RMT tridiagonal limit (different structure); Berglund–Gentz cover the noisy
*fold* (sample-path tubes) but not the *cusp merge with a rate*. B3–B4 (the swept + merging + multiplicative-noise
convergence rate) are genuinely new. **The measured outer contamination (exk overshoots $-1.1$ vs $\mathcal W$'s
$-0.24$) is direct evidence the rate is slow / $p$ is small — the central quantitative unknown.**

**Minimal path:** prove **Tier A** (publishable: "$\mathcal W_\beta$ rigorously characterized"); state **Tier B**
as the universality conjecture with T1 + the numerical fingerprint as evidence. **Most valuable next step:**
complete **A4** (rigorous two-sided tails) — it closes Tier A and is the shared core of Tier B.
**Difficulty:** A moderate (tools exist); B hard (multi-year; the cusp-rate is the crux).

---

## PROGRAM 2 — Exact integrable identity / closed form

Three routes; the fixed-point theorem kills the first, falsifies the third, and forces the second.

### Route 2a — Riemann–Hilbert / τ-function (the "TW = PII" analog) — **LIKELY OBSTRUCTED**
The analog of TW = Fredholm det of the Airy kernel would be a determinant from the RH problem with the asymmetric
Weber Stokes data we found ($\Gamma(3/4-\lambda/4)$ real-Weber × $e^{3i\pi/4}\Gamma(3/4-i\lambda/4)$ imaginary-Weber).
**Obstruction:** TW's tractability comes from the *isomonodromic deformation* (the PII flow) that makes the
determinant a $\tau$-function. We **proved** the cusp is an isomonodromy **fixed point** — no deformation, no
$\tau$-flow — and that $\mathcal W$ is **not** a soft-edge determinant (unbounded-below, no bottom eigenvalue). So
2a yields the *deterministic* connection (**done, closed-form** this session) but **not** a determinantal structure
for the stochastic law. Same wall as the no-σ-ODE theorem. **Do not pursue for the full law** (a clean *negative* —
"no determinantal form" — would itself be a theorem worth writing).

### Route 2b — resurgent trans-series / stochastic exact-WKB — **THE PROMISING ROUTE**
The modern "closed form" for a transcendent with no ODE reduction: a **resurgent trans-series** — Borel-summable
perturbative sectors (in $\eta^2=4/\beta$) + non-perturbative instanton sectors ($\sim e^{-s^5/10\eta^2}$),
connected by resurgence (Stokes constants). It **embraces** the irreducibility instead of fighting it.
**In hand:** the structure; leading perturbative coefficients; the Borel singularity at the instanton action;
and now (MC-free FP, `transseries_orders.py`) the **asymptotic/resurgent signature** — variance-series coefficient
ratios **grow** ($1.23\to2.33$ per order) ⇒ divergent series ⇒ Borel singularity at an $O(1)$ instanton action.
**Open / obstruction:** requires developing **stochastic exact WKB** — the Voros-symbol/Borel machinery for the
*stochastic* Weber operator. Deterministic Weber Voros coefficients are known in closed form (Iwaki–Koike–Takei,
Bernoulli numbers); the stochastic generalization is new (a genuine subfield contribution).
**Honest prerequisite [⚑]:** the leading trans-series *coefficients* are not yet pinned — FP and MC disagree at
small $\eta$ (skew/$\eta\approx0.7$–$0.8$ FP vs $\approx1.2$ MC; var-coeff $0.20$ vs $0.134$), because *both* are
numerically fragile there (FP boundary layer; MC sparse deep tail). They **agree at $\beta=2$** (skew $0.607$ vs
$0.61$). **Reconciling the small-$\eta$ numerics is step 0** of a quantitative resurgent representation.
**Minimal path:** (0) reconcile small-$\eta$ FP/MC; (i) perturbative sectors to $\sim5$ orders; (ii) the Borel-plane
singularity structure (at the instanton actions we have); (iii) the leading Stokes constant. A *partial* resurgent
representation (few orders + instanton + the resurgence relation) is publishable as "the defining resurgent
structure of $\mathcal W$." **Difficulty:** partial — moderate/feasible; full stochastic exact-WKB — hard.

### Route 2c — resonance/spectral expansion — **FALSIFIED for the law**
The closed-form resonance condition (Γ-equation, $\lambda_0=0.89-0.89i$) characterizes the *operator*, but the
Gamow/resonance-survival expansion gives the **wrong skew sign** (`PERSISTENCE_GAMOW_NOTES`) — resonances are
deterministic operator data, not the noisy law. Not a route to $\mathcal W$'s law.

**Program 2 net:** the fixed-point theorem forces the answer — **the closed form, if it exists, is a resurgent
trans-series (2b), not a Painlevé/Fredholm object.** Achieving it means developing stochastic exact-WKB; a partial
resurgent representation is the feasible publishable milestone.

---

## Recommendation (priority order)

1. **Program 1, Tier A** (rigorous characterization) — most tractable, self-contained, highest ROI. Gated on the
   rigorous two-sided tails (A4), which also feed Tier B. **Do this first.**
2. **Program 2, Route 2b step 0** (reconcile small-$\eta$ FP/MC) then the perturbative sectors — the feasible half
   of the closed-form program; instantiates the resurgent representation.
3. **Program 1, Tier B** (universality with rate) and **full stochastic exact-WKB** — the hard, multi-year cores;
   pursue after 1–2 establish the foundations. The cusp-convergence *rate* (Tier B) is the single deepest unknown.
