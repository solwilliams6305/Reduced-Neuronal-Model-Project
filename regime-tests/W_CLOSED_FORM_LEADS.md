# Cross-disciplinary leads for a closed form / exact representation of 𝒲 (deep-research synthesis)

_July 2026 (incoming agent, Fable 5). Synthesis of a 5-angle deep-research sweep (104 agents, 22 primary
sources, 25 claims adversarially verified — 20 confirmed, 5 killed) for the open problem #17: a closed-form or
exact-representation characterization of 𝒲 (the asymmetric inverted-Weber first-passage law). Tags
**[confirmed]/[refuted]/[template]/[gap]**. Sources are peer-reviewed arXiv/journal primary papers._

## Headline (honest)

**No discipline has solved an object equivalent to 𝒲.** Every solved analogue found is **deterministic,
symmetric, and equilibrium/bounded** — none simultaneously has 𝒲's three obstructions (asymmetric
confining-vs-oscillatory connection, unbounded-below *first-passage* rather than an eigenvalue gap, no 1-D-ODE
reduction). This **corroborates the novelty claim** (§4 of `COUPLED_CUSP_RESULTS.md`) but means there is no
ready-made closed form to lift — the path is to *graft* two mature toolsets onto the stochastic-operator setting.

## Lead 1 — Exact WKB / resurgence on the Weber skeleton (strongest technique match) [confirmed]

This is the technique class that computes connection/Stokes data of **irregular-singularity ODEs via Borel-summed
Voros periods/symbols and Stokes-graph topology, with NO Painlevé reduction** — exactly what our fixed-point
theorem says we need (there is no σ-ODE to reduce onto). Concretely:
- **The Weber equation's Voros coefficients are already known in closed form** (Bernoulli numbers, via
  Iwaki–Koike–Takei topological recursion, [arXiv:1805.10945], AHP 2022). A ready connection-data template for our
  deterministic skeleton. **[confirmed, unanimous]**
- A dedicated **"Weber-type exact WKB"** for *inverted / unbounded-below* potentials exists (Kamata–Misumi–
  Pazarbasi–Taya, [arXiv:2604.05878], 2026): uses $[-h^2\partial_y^2+y^2/4-hE]\phi=0$ as the saddle local model and
  derives Weber connection formulae, monodromy, Voros multipliers, and trans-series across **PT-symmetric /
  resonance / anti-resonance** boundary conditions on equal footing. Closest existing machinery to our skeleton.
  **[confirmed]**

**Honest caveats (verified):** (i) standard exact-WKB Assumption 1.1 needs $Q(x)$ with zeros/poles; pure
$Q=-x^2$ *violates* it ([arXiv:2512.17599]) — must use the Weber-type variant, not the generic theory
**[refuted 0-3 that generic EWKB applies directly]**; (ii) EWKB yields **trans-series requiring Borel resummation,
not literal elementary closed form** **[refuted 0-3 that it "solves" the connection problem in closed form]**;
(iii) all of it is **deterministic** — the genuine gap is the *stochastic* generalization.

## Lead 2 — Inverted oscillator / c=1 matrix model (strongest template) [confirmed]

The $-Y^2$ skeleton is **solved exactly** in 2D string theory:
- single-fermion reflection coefficient $R(\omega)=\tfrac1{\sqrt{2\pi}}(i\alpha\mu)^{i\alpha\omega}\Gamma(\tfrac12+i\alpha\mu-i\alpha\omega)$
  — a **closed-form Gamma function of complex argument** ([arXiv:hep-th/0309148]);
- closed-form scattering phase $e^{i\phi_0(E)}=\sqrt{1/2\pi}\,e^{-(\pi/2)(E-i/2)}\Gamma(iE+\tfrac12)$, density
  $\rho(E)=-\tfrac1{2\pi}\phi'(E)$ ([arXiv:hep-th/0208034], Kostov);
- exact **Fredholm-determinant** grand partition function $Z(\mu)=\mathrm{Det}(1+e^{-\beta(\mu+H_0)})$ (free fermions).

This is the **template for the confining-side connection data** (the $\Gamma(\tfrac12-\lambda/2)$-type factors we
already see in `T2_CONNECTION_DATA_DERIVATION`) and possibly a Fredholm form. **Honest caveats:** symmetric,
deterministic, *equilibrium* (a partition function, not a first-passage law); the scattering phase is **not** a
complete characterization **[refuted 0-3]**. Relevance is by shared skeleton/template only.

## "No free lunch" on the other directions (useful negatives) [confirmed]

- **Persistence / first-passage** (Bray–Majumdar–Schehr; random-acceleration [arXiv:2302.04029]): non-Markovian
  first-passage of *this class* has itself resisted closed form; only **weak-noise / instanton tails** are
  tractable — which **matches our derived left tail** ($-\log P\sim s^5/10\eta^2$). Confirms our tail method is
  the right (and only currently-tractable) piece; won't give the full law.
- **Hyperuniformity** (Ghosh–Lebowitz [arXiv:1608.07496]): rigorous class-I rigidity exists only as *sufficient*
  structure-factor conditions (pair correlation decay $\lesssim|x|^{-2}$ in 1-D), no clean necessity — characterizes
  our node process's *class* but yields no exact structure factor.
- **Non-Hermitian / complex-scaling** (Morikawa–Ogawa [arXiv:2508.09211]): technique matches (unbounded-below,
  S-matrix + resonances via EWKB+complex scaling) but worked examples are **symmetric Rosen–Morse barriers**, not
  Weber/asymmetric.

## Specific gem to check [template]

**PIV σ-form ↔ parabolic-cylinder (Weber) *kernel* Fredholm determinant** ([arXiv:2301.05807]): the σ-form of the
Clarkson–McLeod PIV solutions has a determinantal representation with an *integrable operator whose kernel is the
Weber kernel*. We proved 𝒲 is **not** the standard PIV σ-form, so this is a *related-but-distinct* object — but it
is the closest existing "Weber-kernel Fredholm" and should be checked directly against our FP 𝒲 (tails, β-family)
to see what deformation separates them.

## Recommended concrete path (the synthesis) — the four open questions

1. **Graft the c=1 Γ-function (confining side) onto the oscillatory side** to build 𝒲's *asymmetric* connection
   formula; test whether it reproduces the derived $(5,3)$ tails and the β-family. *(Most direct shot at a closed
   form.)*
2. **Stochastic exact-WKB:** does the Weber-type EWKB Voros/period structure admit a noise generalization mapping
   our validated 2-variable backward-Kolmogorov PDE onto a Borel-summable period — an operator-limit analogue of
   RRV's Riccati diffusion, but for $q=2$?
3. **First-passage Fredholm vs. obstruction:** is there a Fredholm/Pfaffian representation for a *first-passage*
   law of an unbounded-below operator (vs. the c=1 *equilibrium* $\mathrm{Det}(1+\dots)$), or is the missing bottom
   eigenvalue a fundamental obstruction to any determinantal form? *(A clean negative here would itself be a
   theorem — cf. our existing "not determinantal" results.)*
4. **Node process:** upgrade class-I hyperuniform rigidity to an exact structure factor by connecting
   Dyson–Schmidt / Frisch–Lloyd 1-D-random-operator machinery to the exact-WKB connection data.

**Net:** the standard-tools room is confirmed empty (no one solved this there); the light now points at a
concrete new frame — **Weber-type exact WKB + c=1 Γ-function connection data, generalized to the stochastic
operator** — which is a place to *work*, and where either a closed form emerges or a clean determinantal-obstruction
theorem does. Both are publishable. Sources are logged in the deep-research output (`tasks/wfa984k5a.output`).

---

## STEP 1 EXECUTED — the confining-side Γ-factor is CONFIRMED in our skeleton [NUMERIC ✓✓, non-confounded]

_July 2026. First concrete step on the deterministic skeleton $u''=(\operatorname{sign}(Y)Y^2-\lambda)u$.
Scripts `coupled-atlas/connection_gamma_test.py`, `connection_poles.py`; figure `figures/connection_poles.png`._

- **Anchor confirmed:** the canard (attracting, $p=+\sqrt V$) solution's first oscillatory node is
  $Y^\star_{\rm det}(\lambda{=}0)=-2.188$ — matches the program's quoted $-2.19$ exactly.
- **Confining-side connection coefficient = c=1/Weber $\Gamma(3/4-\lambda/4)$ — pinned by the DECISIVE pole test.**
  The first node depth $\Theta_1(\lambda)$ decreases and hits $0.000$ **exactly at $\lambda=3.00$ and $7.00$**
  (mean deviation $0.01$), where $u(0)=0$ and a $Y>0$ node appears. These are precisely the poles of
  $\Gamma(3/4-\lambda/4)$ (at $\lambda=3,7,11$), i.e. the zeros of the textbook recessive Weber value
  $U(a,0)=\sqrt\pi\,2^{-a/2-1/4}/\Gamma(3/4+a/2)$ with $a=-\lambda/2$ (DLMF §12.2.6). So the confining side of
  𝒲's connection **is** the c=1 Γ-function the research pointed to — confirmed in-model. **[NUMERIC ✓✓]**
- **Honest negative (a real methodological finding):** the *smooth* connection phase $\phi(\lambda)$ is
  **confounded** — after de-linearizing, its nonlinear residual correlates $\approx0.90$ with *every* Γ-candidate
  AND with a WKB-only null. So the smooth phase cannot discriminate the Γ-form; **only the non-confounded pole
  structure is decisive** (consistent with the program's recurring "smooth fits pass vacuously" theme). The
  precise Γ-identification therefore comes from the pole/bound-state test, not phase-fitting.

**Step (1b) — oscillatory-side factor: a structural obstruction on the real axis [NUMERIC, negative].** Tested
whether the recessive $U$-normalized oscillatory amplitude $M(\lambda)^2$ follows the c=1 modulus law
$|\Gamma(1/2-i\lambda/2)|^2=\pi/\cosh(\pi\lambda/2)$ (`oscillatory_gamma.py`). **It does not:** large-$\lambda$
slope of $\log M^2$ is $-0.04$ (predicted $-\pi/2$), ratio $M^2/\cosh$ ranges $677\to0.68\to21$ — the cosh law is
**refuted for this observable.** **Diagnosis [DERIVED, clean]:** the first-passage solution is **real**, so on
$Y<0$ it is a standing wave whose two counter-propagating amplitudes are complex conjugates of equal modulus —
**there is no genuine connection modulus on the real axis; the oscillatory connection is pure phase (which §Step-1
already found confounded).** The $\Gamma(1/2-i\lambda/2)$ / cosh factor exists only in the **complex** extension.
So the oscillatory factor is **not** accessible by any real-axis amplitude; it is encoded in the **resonances**
(M1: $\lambda_0=0.86-0.82i,\,2.30-1.22i,\,4.14-1.08i$).

**Step (1c) — the resonance probe CLOSED both factors [NUMERIC ✓✓, machine precision].** See
`CONNECTION_CLOSED_FORM_NOTES.md`. The resonance condition (recessive on $Y>0$, outgoing on $Y<0$; match
log-derivatives at $Y=0$) is, in **closed form**:
$$\underbrace{e^{3i\pi/4}\tfrac{\Gamma(3/4-i\lambda/4)}{\Gamma(1/4-i\lambda/4)}}_{\text{oscillatory Jost }L(\lambda)/2\ \text{(imag Weber)}}=\underbrace{\tfrac{\Gamma(3/4-\lambda/4)}{\Gamma(1/4-\lambda/4)}}_{\text{confining }R(\lambda)/2\ \text{(real Weber)}}.$$
Verified to **machine precision**: $L(\lambda)/[2\Gamma\text{-ratio}]=e^{3i\pi/4}$ exactly (mod 1.0000, arg
$0.75\pi$) across 9 real & complex $\lambda$; the Γ-equation root is $\lambda_0=0.8896-0.8896i$
($|{\rm cond}|=10^{-16}$) = M1's independent FD resonance $0.890-0.890i$. **Both connection factors are now
closed-form Weber Γ-ratios** — confining (real $\lambda$) × oscillatory (imaginary $\lambda$) + the $e^{3i\pi/4}$
inversion Stokes phase. This realizes the deep-research recommendation exactly and **supersedes the falsified M3**
(a real-axis meromorphic $D(\lambda)$; the correct object is the real×imaginary Weber Γ-*equation*).

**Net of Step 1 (complete):** the **deterministic** asymmetric connection / resonance data of 𝒲's skeleton is in
**closed form** (both Γ-factors, machine-verified). Remaining for a closed form of the stochastic law 𝒲:
**(iii) the stochastic generalization** — integrate the noise against this now-exact connection (the q=2 analogue
of RRV's Riccati diffusion). The skeleton it lives on is no longer a mystery; the noise step is the frontier.
