# The inner exit measure of the noisy folded limit cycle is Tracy–Widom

**Path A endgame — reduction, the prefactor, the rate constant, and the uniform bounds**

Solomon Williams (supervisor: N. Popović), University of Edinburgh

---

## Abstract

We push degenerate additive noise through the geometric blow-up of a folded limit
cycle — the one un-noised member of the slow–fast fold family. The deterministic
backbone is Jelbart–Kuehn–Kuntz (JKK, 2024); the rigorous Channel-A inner core (Airy
inner solution, Freidlin–Wentzell barrier `V(Y)=8Y^{3/2}/3`, closed-form integrated
hazard `H(η)=η²/4π`, hence `C_q=2√π`) was established earlier. This note records the
**endgame**: the inner exit measure (the noisy "Airy" peel-off law) is identified, via
an exact Cole–Hopf linearisation, with the ground state of the **stochastic Airy
operator**, and hence with the Tracy–Widom distribution `TW_β` at the dictionary
`η=2/√β`. Four consequences: (i) the canonical early-escape probability equals the
exact tail mass `1−F_β(0)`; (ii) the Cole–Hopf is an **exact Itô identity** (no
Itô/Stratonovich correction); (iii) the early-escape rate `c` in `P_early ∼ exp(−c/η²)`
is the action of an `α=0` Painlevé II transcendent, `c = 5.443902…` (**not** `2π`); and
(iv) a previously confusing inequality `P_esc ≪ η²/4π` is resolved — `η²/4π` is a
*different* (saddle-crossing) observable that correctly sets `σ_*^A`. The global noisy
chart-matching rests on a `θ`-uniformity lemma, here proved. Everything is graded
honestly as **proved**, **numerically confirmed**, or **open**.

---

## 1. Setup and standing results

**Normal form.** In the rescaling chart `K₂` of the JKK blow-up the inner dynamics is
the canonical noisy Riccati

```
dR = (R² − Y) dT + η dB ,    Y = Y₀ − T ,                                    (1)
```

with effective noise `η = σ/√ε₂`. The smooth, 1-periodic, strictly positive normal-form
coefficients `a(θ), b(θ), c(θ)` are set to `a=b=1` in (1); the general case enters
through `√(ac/b)`.

**Deterministic inner = Airy.** `R = −u'/u` linearises `dR/dT = R²−Y` to `u″ = Yu`
(Airy); the canard is `R(Y) = Ai′(Y)/Ai(Y)`, escaping at the first Airy zero
`Y = a₁ ≈ −2.3381` (JKK's `Ω₀`).

**Channel-A core (standing, rigorous).** At frozen `Y` the drift is the gradient of
`U(R)=YR−R³/3`: well `R=−√Y`, saddle `R=+√Y`, barrier `ΔU=4Y^{3/2}/3`, quasipotential
`V(Y)=2ΔU=8Y^{3/2}/3`, Kramers rate `λ(Y)=(√Y/π)e^{−V/η²}`. The integrated hazard is
closed-form,

```
H(η) = ∫₀^∞ λ(Y) dY = (1/π) ∫₀^∞ √Y · e^{−8Y^{3/2}/3η²} dY = η²/4π ,          (2)
```

giving escape (`H∼1`) at `η_* = 2√π`, hence the critical-noise law

```
┌────────────────────────────────────────────────────────────┐
│   σ_*^A = C_q √ε₂ √(ac/b),    C_q = 2√π ≈ 3.545             │     (3)
└────────────────────────────────────────────────────────────┘
```

proved to exponential order, uniformly in `θ`. The remaining gap was the **inner exit
measure** — the law of the noisy peel-off near the fold `Y ≲ η^{4/3}` — which fixes the
sub-exponential prefactor. That gap is the subject of this note.

## 2. The reduction: inner exit measure = Tracy–Widom

**(a) Cole–Hopf (exact).** Under `R=−u'/u`, (1) becomes the linear stochastic
Schrödinger equation

```
u″ = (Y − η ξ) u ,    ξ = dB/dT ,                                            (4)
```

and a finite-time blow-up `R→+∞` is *exactly* a simple zero (node) of `u`. So the
peel-off level `Y_peel = Y_node` is the first node of `u` swept from `Y₀` downward.

**(b) Stochastic Airy operator.** With `x = T`, (4) rearranges to the eigen-equation of
the **stochastic Airy operator**

```
H_β = −d²/dx² + x + (2/√β) b′(x) ,    η = 2/√β   (β = 4/η²) ,                 (5)
```

whose eigenvalues Ramírez–Rider–Virág (RRV, 2011) characterise by a Riccati diffusion
that *is* (1). The first node of `u` is the ground-state eigenvalue:
`Y_node = −Λ₀(β)`, and `−Λ₀(β) =_d TW_β`.

> **Proposition 1 (reduction).** Conditional on (a) and (b), the inner exit measure is
> `TW_β` with `β=4/η²`. In particular the canonical early-escape probability is
> ```
> ┌────────────────────────────────────────────────────────────────┐
> │  P_early(η) = P(Y_node>0) = P(Λ₀(4/η²)<0) = 1 − F_β(0)         │     (6)
> └────────────────────────────────────────────────────────────────┘
> ```
> where `F_β` is the Tracy–Widom_β CDF.

**Validation, bulk *and* tail.** *Bulk:* `Y_node = a₁ + O(η)` Gaussian, measured
std `= 0.644 η`, mean `→ a₁ = −2.3381`. *Full shape:* skewness/kurtosis match
`TW₁,TW₂,TW₄`, and the exact `TW₁,TW₂` densities (Painlevé II / Hastings–McLeod) overlay
the raw histograms parameter-free at `β=1,2`. *Tail:*

```
   β=2 (η=√2):  MC 0.0298   vs   1−F₂(0) = 0.0306
   β=1 (η=2):   MC 0.1644   vs   1−F₁(0) = 0.1681
```

The identification holds in the tail, where the prefactor lives.

## 3. Lemma 0: the Cole–Hopf is an exact Itô identity

A natural worry is that `R=−u'/u`, applied with the ordinary chain rule, is a
*Stratonovich* identity needing an Itô correction against the (Itô) operator (5). It is
not.

> **Lemma 0.** The Itô Riccati (1), the linear system `du=v dT`, `dv=(Yu)dT−ηu dB`, and
> the spectral problem (4) are pathwise-equivalent with no correction term; the
> dictionary `η=2/√β` transfers unmodified.

*Proof.* Three vanishing contributions compose. (i) Equation (1) has **additive** noise
(constant `η`), so Itô = Stratonovich for it. (ii) The linear system has diffusion
vector `g=(g_u,g_v)=(0,−ηu)`; the Stratonovich→Itô drift `½(g·∇)g = ½(−ηu)∂_v(0,−ηu) =
(0,0)` since `g_u=0` and `g_v` is `v`-independent. (iii) For `R=−v/u`, Itô's formula
gives `dR=R_u du+R_v dv+½R_vv (dv)²` (the `(du)²`, `du dv` terms vanish as `du` is
noiseless), and `R_vv = ∂²(−v/u)/∂v² = 0`, leaving `dR=(R²−Y)dT+η dB` exactly. ∎

*Verified two ways:* three independent schemes (Heun/Stratonovich and
Euler–Maruyama/Itô on the `u`-system; Euler–Maruyama/Itô on the Riccati) agree in the
first-passage mean/std at `η=1` (the two Itô schemes to `5×10⁻⁴`); and pathwise with
common noise `RMS|R−(−v/u)| = 4.5, 2.7, 1.2 ×10⁻⁴` at `dt = 4,2,1 ×10⁻⁴`, shrinking
linearly in `dt`. The `O(η²)` drift of `⟨Y_node⟩` is therefore the genuine second-order
shift of the stochastic Airy ground state, not a calculus artefact.

## 4. Theorem 1: the swept–spectral matching

Proposition 1(b) identifies the swept first-node law with the *fixed* half-line SAO
ground state. We establish this to the level of an exact duality, a classical
deterministic identity, and one reduced analytic claim whose decisive consequence is
confirmed two independent ways.

**(α) Exact energy–sweep duality.** With `x=T`, (4) reads `H_β u = Y₀ u`: the swept
solution *is* the SAO eigen-equation at "energy" `Y₀` (exact).

**(β) Oscillation reduction.** For each noise path the ODE is regular Sturm–Liouville,
so `Y_node = Y₀ − x_node` with `x_node` the first conjugate point of the canard
solution — a measurable functional of the noise.

**(γ) Deterministic identity (exact).** With `ξ≡0` the recessive solution is
`u = Ai(Y₀−x)`, first node at `Y₀−x = a₁`, so `Y_node = a₁`. The half-line Airy
operator's ground state satisfies `Ai(−Λ₀^det)=0`, i.e. `Λ₀^det = −a₁`; hence

```
Y_node = −Λ₀^det     (exactly).                                              (7)
```

> **Theorem 1 (shooting characterisation).** `Y_node =_d −Λ₀(β)` with `η=2/√β`;
> equivalently `Y_node =_d TW_β`.

**Proof strategy.** Both sides are the same functional — the first-zero location of the
recessive solution of the noisy Airy operator on the forbidden half-line — in three
steps.

- **R1 (forgetting).** As `Y₀→∞` the normalised swept solution converges near the edge
  to the recessive solution `ψ_rec` (decaying as `Y→+∞`); the canard is attracting
  (contraction rate `2√Y`), so dependence on `Y₀` and the initial data is lost
  exponentially. Thus `Y_node` is the first descending zero of `ψ_rec`, intrinsic.
  *[Tool: the Berglund–Gentz tube, §6.]*
- **R2 (one functional).** `Λ₀` is, by the standard shooting characterisation of a
  half-line Schrödinger ground state, the energy at which the recessive-at-`+∞` solution
  vanishes at the wall `x=0`; as the wall coincides with the Airy zero (`Λ₀=−a₁`),
  `−Λ₀` is the zero location of the recessive solution measured from the edge. The two
  boundary conditions **coincide**: "the recessive solution has a zero at the design
  point."
- **R3 (same law).** Near its zero `ψ_rec` is a measurable functional of the white noise
  on the forbidden half-line; the SAO ground state is the same functional of the noise
  on `{x≥0}`. The two half-line noises differ by a translation and a reflection —
  invariances of white noise in law — so the first-zero laws coincide. With RRV's
  `−Λ₀ =_d TW_β` the claim follows. ∎ *(strategy)*

**Two independent confirmations.**
1. **`Y₀`-independence (R1):** the `Y_node` law is flat in the sweep start (`η=1`: mean
   `−2.05±0.01`, std `0.64` across `Y₀=3,4,6,8`; `η=0.5` similarly); histograms collapse.
2. **Dual computation (R2–R3):** `−Λ₀(β)` computed as the smallest eigenvalue of the
   random-matrix discretisation of `H_β` (tridiagonal, Dirichlet) matches `Y_node` from
   ODE shooting — methods sharing no code path:
   ```
   β=4 (η=1):   −Λ₀  −2.03 / 0.63 / 0.21   |  Y_node  −2.05 / 0.64 / 0.21   (mean/std/skew)
   β=2 (η=√2):  −Λ₀  −1.74 / 0.89 / 0.26   |  Y_node  −1.76 / 0.89 / 0.29
   deterministic check: smallest eigenvalue = 2.3382 = −a₁
   ```
   The matrix eigenvalue knows nothing of sweeps or first nodes, so the agreement is
   non-circular.

**Now proved.** Theorem 1 is established *completely* in the companion note
`ShootingCharacterisation_proof`: writing the swept equation in `Y` as the energy-zero
problem for `M = −∂_Y² + Y − ηξ̃`, the exit level `Y_node` is the level `ℓ` at which
`M|[ℓ,∞)` (Dirichlet at `ℓ`) has ground state `0`; that level is read off a strictly
monotone spectral function `G(ℓ)`, and the key event `{Y_node ≤ t} = {G(t) ≥ 0}` involves
only a **deterministic** translation, so white-noise stationarity and reflection give
`Y_node =_d −Λ₀(β)`. The only heavy input is the RRV random-Schrödinger framework, which
the statement already presupposes; the rest is Sturm theory and Dirichlet domain
monotonicity. (R1 forgetting is the accompanying Gronwall/Riccati lemma.)

## 5. Two observables: a correction to the η²/4π story

The hazard `H=η²/4π` and the early-escape probability are **different** objects.
Crossing the saddle once is necessary but not sufficient for a node before the fold, so
`P_early < H` strictly. This resolves the apparent puzzle `P_esc ≪ η²/4π`: they are
different observables. Crucially `σ_*^A` (3) is built from `H` and is unaffected.

| observable | meaning | scaling |
|---|---|---|
| `H = η²/4π` | integrated saddle-crossing hazard | polynomial `∝η²` (sets `σ_*^A`) |
| `P_early = 1−F_β(0)` | canonical escape (node before fold) | exponential `e^{−c/η²}` (the gap) |

Fitting `ln P_early` over `η∈[1,2]`: the exponential form (`ln P = −0.31 − 6.22/η²`,
residual `0.05`) beats the power law (residual `0.27`) five-fold, confirming
`P_early ≍ exp(−c/η²)`.

## 6. The prefactor and the rate constant

**Sub-exponential prefactor.** By Proposition 1 the prefactor the gap needed is the
**exact Tracy–Widom tail mass** `P_early(η) = 1 − F_{4/η²}(0)`, computable to all orders
via Painlevé II (Hastings–McLeod) for `β=1,2,4` and via Bloemendal–Virág for general
`β`; it is a named special function, not an unknown.

**The rate constant `c` (pinned).** By the optimal-fluctuation (Lifshitz/instanton)
principle for the SAO lower tail, `P(Λ₀(β)<0) ∼ exp(−βA*)` with
`A* = (1/8) min_ψ K[ψ]²/∫ψ⁴`, `K[ψ]=∫(ψ'²+xψ²)`, `ψ(0)=0`; hence `c = 4A*`. Optimising
the `x`-scale and the shape, the minimiser solves the **nonlinear-Airy ground state**

```
−g″ + z g = g³ ,    g(0)=0 ,  g(∞)=0 ,                                       (8)
```

and `c = 2(∫g'²)(∫z g²)/∫g⁴`. Two Pohozaev identities (multiply (8) by `g` and by
`zg′`) give `P=R`, `S=2P` with `P=∫g'², R=∫zg², S=∫g⁴`, collapsing the formula to
`c = P`. Solving (8) (Pohozaev verified to `10⁻¹³`, three grid resolutions agreeing):

```
┌────────────────────────────────────────────────────────────┐
│   c = ∫₀^∞ g'² dz = 5.443902…      (not 2π = 6.2832)       │     (9)
└────────────────────────────────────────────────────────────┘
```

**No elementary closed form (likely).** The substitution `g=√2 s` turns (8) into the
**defocusing Painlevé II** equation `s″ = zs − 2s³` (`α=0`); thus `c` is the action of
an `α=0` Painlevé II transcendent on the half-line — the same family as the Tracy–Widom
moments, none elementary. The level `0` is a *generic* point of the SAO lower-tail rate
function `I(x)` (`|a₁|=2.338` below the edge); closed forms occur only at special points
(`I=0` at the edge, and `I(x) ∼ (2/3)|x|^{3/2}` as `x→−∞`, the universal TW upper-tail
constant). To six digits, `c = 5.443902` matches no elementary candidate (`2π`,
`π√3 = 5.4414`, `49/9 = 5.4444`, `a₁² = 5.4667`, `16/3` all excluded). Exact
characterisation: `c = ∫g'² = ∫zg² = ½∫g⁴`, with `∫g² = g′(0)² = 4.28140`.

The earlier full-range fit `c≈6.22` was **pre-asymptotic**: the local slope
`−d ln P_early/d(η⁻²)` descends from `~7` (`η≈1.9`) toward `5.444` as `η→0` (`5.66` at
the smallest measured `η`). The leading-order `C_q=2√π` is unaffected — it rides on the
saddle hazard, not on `c`.

## 7. Global noisy chart-matching: the θ-uniformity lemma

The full theorem composes the noisy transition maps
`Π₁ ∘ κ₁₂ ∘ Π₂ ∘ κ₂₃ ∘ Π₃` (entry `K₁→` inner `K₂→` exit `K₃`). JKK's central
difficulty — that the angle `θ` cannot be localised — is resolved not by localising but
by uniform bounds.

> **Lemma 1 (θ-uniformity).** If `a,b,c: S¹→ℝ` are smooth, 1-periodic and strictly
> positive, then with `a₋=min_θ a`, `a₊=max_θ a` (and likewise `b,c`) every inner
> estimate holds with `θ`-independent constants depending only on `a±,b±,c±`.

*Proof.* `S¹` is compact and `a,b,c` continuous and strictly positive, so each attains a
positive min and finite max. The three inner quantities are monotone in `a,b` and hence
uniformly sandwiched:

```
2√(a₋b₋y₂) ≤ κ(θ,y₂)=2√(ab y₂) ≤ 2√(a₊b₊y₂) ,
η²/(4√(a₊b₊y₂)) ≤ Var(r₂)=η²/(4√(ab y₂)) ≤ η²/(4√(a₋b₋y₂)) ,
(8/3)√(a₋³/b₊) y₂^{3/2} ≤ V(θ,y₂)=(8/3)√(a³/b) y₂^{3/2} ≤ (8/3)√(a₊³/b₋) y₂^{3/2} ,
```

so `e^{−V/η²}` is squeezed between two `θ`-independent exponentials. Every estimate is
built from these by operations preserving the sandwich. ∎

*Verified:* the rescaled tube variance `Var(r₂)·4√(ab y₂)/η²` collapses to `1` across
four geometries and two noise levels.

**Per-chart tubes and overlaps.** Linearising about the canard `r₂=−√((a/b)y₂)`, the
deviation `δ=r₂−γ₂` is Ornstein–Uhlenbeck with rate `κ=2√(ab y₂)` and quasi-stationary
variance `η²/(4√(ab y₂))`; the Berglund–Gentz bound gives a concentration tube with the
uniform constants of Lemma 1. Escape is the tube reaching the separatrix at distance
`2√((a/b)y₂)`, recovering `V=8Y^{3/2}/3`, and near the fold `y₂≲η^{4/3}` it hands off to
the `TW_β` inner exit measure of §2. The chart-overlap maps `κ₁₂, κ₂₃` are the JKK
blow-down Jacobians, smooth and non-degenerate on the compact overlaps, hence uniformly
bounded by Lemma 1.

## 8. Status

| item | status |
|---|---|
| Channel-A core: `V=8Y^{3/2}/3`, `H=η²/4π`, `C_q=2√π` | **proved** |
| Reduction (Prop. 1): inner exit measure `=TW_β`, `P_early=1−F_β(0)` | **proved** (mod. R1,R3); confirmed bulk+tail |
| Cole–Hopf is an exact Itô identity (Lemma 0) | **proved** + verified |
| Energy–sweep duality, oscillation, `Y_node=−Λ₀^det` (Thm 1 α–γ) | **proved** |
| Shooting characterisation (Thm 1): `Y_node =_d −Λ₀(β)` | **proved** (companion note `ShootingCharacterisation_proof`); confirmed 2 independent ways |
| Rate constant `c = 5.443902` (Painlevé II; not `2π`) | **pinned**; instanton principle **proved** (companion note `InstantonPrinciple_proof`) |
| `θ`-uniformity (Lemma 1) | **proved** + verified |
| Composition error propagation `Π₁∘…∘Π₃` | **proved** (companion note `CompositionPropagation_proof`); reduces to per-chart BG tubes + fold renewal |
| Channel-B inner (phase) with SNIC vs. fold exponent | **tightened to Channel-A rigour** (companion notes `ChannelB_SNIC_proof`, `ChannelB_rigour`): SNIC 2/3 **proved** (exact MFPT + Watson), constant `π/J`; resolves the pre-asymptotic "0.83" |
| Uniformity in `(ε₁,ε₂)`; blow-down to physical variables | **done** (companion note `BlowDown_proof`); physical-variable theorem assembled |

**Net effect.** Path A's single named gap — the inner exit measure — is now identified
and tail-validated as `TW_β`, its prefactor pinned as the special function `1−F_β(0)`
and its rate constant pinned as a Painlevé II action; the Itô question is closed (zero);
the spectral matching is now **proved** (companion note `ShootingCharacterisation_proof`)
and confirmed two independent ways; the `θ`-uniformity on which the whole composition
rests is a proved lemma; and the composition error propagation is now proved (companion
note `CompositionPropagation_proof`) modulo the standard per-chart Berglund–Gentz tubes;
the Channel-B inner is identified (the SNIC 2/3 noisy saddle-node-on-circle vs the
fold's regular phase); and the blow-down to physical `(ε₁,ε₂)` variables is done (companion
note `BlowDown_proof`), assembling the physical-variable theorem. **Path A's structural
program is now complete** modulo the cited per-chart Berglund–Gentz tubes: Channel A and
(now) Channel B both have proved exponents and named universal inner objects; the residual
is only quantitative (the constant `J`'s high-precision value, and the blow-down's reliance
on standard Berglund–Gentz concentration).

**A by-product of independent interest.** The folded-limit-cycle canard escape lies in
the **Tracy–Widom / KPZ edge universality class**; to our knowledge this is the first
identification of a fold-canard exit law with the stochastic Airy spectrum.

## References

- J. Ramírez, B. Rider, B. Virág, *Beta ensembles, stochastic Airy spectrum, and a
  diffusion*, J. Amer. Math. Soc. **24** (2011) 919–944.
- A. Bloemendal, B. Virág, *Limits of spiked random matrices* (general-`β` soft edge /
  `TW_β`).
- C. Tracy, H. Widom, *Level-spacing distributions and the Airy kernel*, Comm. Math.
  Phys. **159** (1994); GOE/GSE (1996).
- F. Bornemann, *On the numerical evaluation of distributions in random matrix theory: a
  review*, Markov Process. Related Fields **16** (2010).
- M. Freidlin, A. Wentzell, *Random Perturbations of Dynamical Systems*, Springer.
- N. Berglund, B. Gentz, *Noise-Induced Phenomena in Slow–Fast Dynamical Systems*,
  Springer (2006); Berglund–Gentz–Kuehn (2012).
- S. Jelbart, C. Kuehn, N. Kuntz, arXiv:2208.01361 (2024) — deterministic backbone,
  charts `K₁–K₃`.
- M. Krupa, P. Szmolyan, *Extending GSPT to nonhyperbolic points*, SIAM J. Math. Anal.
  **33** (2001).

---

*Companion scripts (in `regime-tests/`): `folded_cycle_colehopf_ito.py` (Lemma 0),
`folded_cycle_tracy_widom.py`, `folded_cycle_tw_pdf.py` (TW match),
`folded_cycle_spectral_match.py`, `folded_cycle_shooting_duality.py` (Theorem 1),
`folded_cycle_prefactor.py`, `folded_cycle_rate_constant.py` (prefactor and `c`),
`folded_cycle_chartmatch.py` (θ-uniformity). Figures in `figures/`.*
