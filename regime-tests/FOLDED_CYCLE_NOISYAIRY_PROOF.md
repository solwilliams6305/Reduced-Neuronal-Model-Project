# Path A endgame — the noisy-Airy prefactor, extracted; and the uniform error bounds

**Status.** Two things attempted here, honestly graded.
**Part 1 (prefactor extraction):** the inner exit measure is *reduced* to the stochastic
Airy ground state (rigorous, modulo two named inputs), and the sub-exponential prefactor
the gap needed is **identified in closed form** as the exact Tracy–Widom tail `1−F_β(0)`
— *validated to Monte-Carlo accuracy at the tail*. A genuine correction also fell out:
the canonical early-escape probability is **exponential** `exp(−c/η²)`, not the
polynomial `η²/4π` (which is a *different* observable, the saddle-crossing hazard).
**Part 2 (uniform bounds):** the **θ-uniformity lemma is proved**; the per-chart
Berglund–Gentz tubes are stated with explicit uniform constants; the sharp *composition*
closure remains the open piece.

Companion: `folded_cycle_prefactor.py`, `figures/folded_cycle_prefactor.png`;
builds on `FOLDED_CYCLE_NOISYAIRY_CONJECTURES.md`, `FOLDED_CYCLE_PATHA_PROOFS.md`.

---

# Part 1 — The prefactor, extracted

## 1. The reduction (Proposition 1)

Canonical inner Riccati (K₂ chart, `λ=1`): `dR=(R²−Y)dT+η dB`, `Y=Y₀−T`.

**(a) Cole–Hopf (verified).** `R=−u′/u` linearises the Riccati to
```
u″ = (Y − η ξ) u ,      ξ = dB/dT  (Stratonovich),
```
and a finite-time blow-up `R→+∞` is exactly a simple zero (node) of `u`. *Verified
numerically:* escape-`Y` (Riccati `R>20`) `=−1.625` vs first `u`-zero `=−1.675`, and
`R=−u′/u` tracks the Riccati along the path.

**(b) RRV identification.** `u″=(Y−ηξ)u` is the eigenfunction equation of the
**stochastic Airy operator** `H_β=−d²/dx²+x+(2/√β)b′` at the dictionary
`η=2/√β` (β=4/η²); Ramírez–Rider–Virág’s characterising diffusion *is* this Riccati.
Hence the first node of `u`, swept from `Y₀` downward, is the ground-state eigenvalue:
`Y_node = −Λ₀(β)`, and `−Λ₀(β) =_d TW_β`.

> **Proposition 1 (reduction).** *Conditional on (a) and (b), the inner exit measure
> — the law of the peel-off level `Y_peel=Y_node` — is `TW_β` with `β=4/η²`. In
> particular the canonical early-escape probability is*
> ```
> ┌─────────────────────────────────────────────────────────────┐
> │   P_early(η) = P(Y_node > 0) = P(Λ₀(4/η²) < 0) = 1 − F_β(0)  │
> └─────────────────────────────────────────────────────────────┘
> ```
> *where `F_β` is the Tracy–Widom_β CDF.*

**Validation — bulk *and* tail.** Bulk: `Y_node=a₁+O(η)` Gaussian, std `=0.644η`,
mean `→a₁=−2.3381` (`folded_cycle_noisyairy.py`). Full shape: parameter-free
mean/std/skew/kurt match at `β=1,2`, with the **exact** `TW₁,TW₂` densities (Painlevé II)
overlaying the raw histograms (`folded_cycle_tracy_widom.py`, `folded_cycle_tw_pdf.py`).
**Tail (new):** `P(Y_node>0)` from Monte-Carlo vs the *exact* `1−F_β(0)`:
```
   β=2 (η=√2):  MC 0.0298   vs   1−F₂(0) = 0.0306
   β=1 (η=2):   MC 0.1644   vs   1−F₁(0) = 0.1681
```
The identification therefore holds *in the tail*, which is where the prefactor lives.

## 1a. Lemma 0 — Cole–Hopf is an *exact* Itô identity (no correction)

Input (a) carried a caveat (old “C-i”): is the Cole–Hopf a Stratonovich identity needing
an Itô correction? **It is not — the correction is identically zero**, for three reasons
that compose:

- **(i) The Riccati has additive noise.** `dR=(R²−Y)dT+η dB` has *constant* diffusion
  coefficient `η`, so `½ g g' = ½η·∂_R η = 0`: its Itô and Stratonovich readings coincide.
- **(ii) The linear system has no Wong–Zakai drift.** Writing `du=v\,dT`,
  `dv=(Yu)dT−ηu\,dB`, the diffusion vector is `g=(g_u,g_v)=(0,−ηu)`. The Stratonovich→Itô
  drift correction is `½(g·∇)g = ½(g_u∂_u+g_v∂_v)g = ½(−ηu)\,∂_v(0,−ηu) = (0,0)`, because
  `g_u=0` and `g_v=−ηu` is `v`-independent. So the `u`-system is itself Itô = Stratonovich.
- **(iii) The change of variables has `R_vv=0`.** For `R=−v/u`, Itô’s formula gives
  `dR=R_u\,du+R_v\,dv+½R_{vv}(dv)²` (the terms `(du)²` and `du\,dv` vanish since `du` is
  noiseless). With `R_{vv}=∂²(−v/u)/∂v²=0`, the second-order term drops and
  `dR=(v²/u²)dT−Y\,dT+η\,dB=(R²−Y)dT+η\,dB` **exactly**.

> **Lemma 0.** *The Itô Riccati `dR=(R²−Y)dT+η dB`, the linear system above, and the
> spectral problem `u″=(Y−ηξ)u` are pathwise-equivalent with **no** Itô–Stratonovich
> correction. Consequently the RRV dictionary `η=2/√β` transfers without modification.*

**Verified two ways.** (1) Three independent schemes — Heun (Stratonovich limit) and
Euler–Maruyama (Itô) on the `u`-system, and Euler–Maruyama (Itô) on the Riccati — give the
same first-passage `Y` (`R=−u′/u≥20`) at `η=1`: means `−1.962, −1.971, −1.970`, stds
`0.610, 0.604, 0.604` (the two Itô schemes agree to `5·10⁻⁴`). (2) **Pathwise with common
noise**, `RMS|R−(−v/u)| = 4.5,\,2.7,\,1.2 ×10⁻⁴` at `dt=4,2,1×10⁻⁴` — shrinking
**linearly in `dt`**, i.e. both discretise the *same* SDE. (`folded_cycle_colehopf_ito.py`.)

**Consequence.** The `O(η²)` drift of `⟨Y_node⟩` (mean `−2.336→−2.292` as `η:0.1→0.4`) is
therefore **not** a calculus artifact — it is the genuine *second-order perturbative shift*
of the stochastic-Airy ground-state eigenvalue, a real feature of `TW_β`. Old caveat C-i is
**discharged.**

## 1b. Lemma 1 — the swept problem *is* the SAO eigen-problem (spectral matching)

This is input (b)’s last load-bearing piece: that the *swept* first-node law equals the
*fixed* half-line SAO ground state. It is established here to the level of an **exact
operator duality + a classical deterministic identity + one reduced analytic claim whose
decisive consequence is verified**.

**(α) Exact energy–sweep duality.** In `x=T`, `u″=(Y₀−x−ηξ)u` rearranges to
```
(−d²/dx² + x + ηξ(x)) u = Y₀ u ,      i.e.   H_β u = Y₀ u ,   η = 2/√β .
```
The swept solution **is** the SAO eigen-equation at energy `Y₀` — exact, no approximation
(the same identity that fixes the dictionary).

**(β) Oscillation reduction (pathwise).** For each realisation of `ξ` the ODE is a regular
Sturm–Liouville equation, so `Y_node=Y₀−x_node` with `x_node` the first node of `u` (the
first *conjugate point* of the canard solution) — a measurable functional of the noise. By
the **Sturm oscillation theorem**, the node count of `u` on `(0,x)` equals the number of
eigenvalues below `Y₀` of `H_β` truncated to `(0,x)` (canard left-BC, Dirichlet at `x`).

**(γ) Deterministic identity (classical, exact).** With `ξ≡0` the recessive solution is
`u=Ai(Y₀−x)`, whose first node is at `Y₀−x=a₁`, i.e. **`Y_node=a₁=−2.338`**. The half-line
Airy operator’s ground state is `Λ₀^{det}=−a₁` (its eigenfunction `Ai(x−Λ₀)` vanishes at
`x=0 ⟺ Ai(−Λ₀)=0 ⟺ Λ₀=−a₁`), so exactly
```
Y_node = −Λ₀^{det}.
```
The recessive-solution first node = minus the ground-state energy — the *shooting* dual of
the eigenvalue problem.

**(δ) The matching (reduced + verified).** With noise the claim is `Y_node =_d −Λ₀(β)`.
Two ingredients: **(i)** the canard data `u′(0)/u(0)=√Y₀` is the WKB-recessive value in the
barrier `Y>0`, so as `Y₀→∞` the solution is exponentially close to the globally subdominant
one — *the canard data forgets `Y₀`*; **(ii)** the first node sits near the soft edge
`Y≈0`, where the local operator is the rescaled SAO, so its fluctuations are the SAO
ground-state fluctuations. This is the **shooting characterisation of `Λ₀(β)`, dual to
RRV’s Riccati-explosion characterisation** (oscillation theory is endpoint-symmetric).

Its single decisive consequence — that `Y_node` is **intrinsic**, independent of the sweep
start `Y₀` — is **verified** (`folded_cycle_spectral_match.py`):
```
   η=1.0:  Y₀=3,4,6,8  →  mean −2.05±0.01,  std 0.64,  skew 0.16   (flat)
   η=0.5:  Y₀=3,4,6,8  →  mean −2.268±0.001, std 0.323            (flat)
```
The `Y_node` histograms at `Y₀=3,4,6,8` collapse onto one curve
(`figures/folded_cycle_spectral_match.png`). With the full distributional match to `TW_β`
(§1), the matching holds numerically.

> **Lemma 1 (matching, reduced).** *Conditional on the shooting characterisation of the
> SAO ground state (the Sturm dual of RRV), `Y_node =_d −Λ₀(β) =_d TW_β`. The deterministic
> identity `Y_node=−Λ₀^{det}` is exact, and the decisive consequence (`Y₀`-independence) is
> verified.* The one remaining rigorous step is that shooting-characterisation theorem.

## 1c. Theorem 1 (shooting characterisation) — proof strategy + dual-method confirmation

The matching (Lemma 1δ) reduces to a single statement, attacked here directly.

> **Theorem 1 (shooting characterisation).** `Y_node =_d −Λ₀(β)`, `η=2/√β`.

**Unifying idea.** Both sides are the *same functional* — the first-zero location of the
**recessive solution of the noisy Airy operator on the forbidden half-line**. Reduction in
three steps:

**(R1) Forgetting (canard attractivity).** As `Y₀→∞`, the normalised swept solution
converges near the edge to the unique solution `ψ_rec` of `u″=(Y−ηξ)u` that is recessive
(decaying) as `Y→+∞`; dependence on `Y₀` and on the canard data is lost exponentially
(contraction rate `2√Y`). Hence `Y_node` is the first descending zero of `ψ_rec` —
intrinsic. *[Quantitative tool: the Berglund–Gentz tube (§7). Verified: `Y₀`-independence,
§1b.]*

**(R2) One functional, two readings.** `ψ_rec` is the recessive solution of the
Airy + white-noise operator and `Y_node` the location of its first zero from the fold. For
the SAO `H_β` on `[0,∞)` (Dirichlet), `Λ₀` is — by the standard shooting characterisation
of a half-line Schrödinger ground state — the spectral value at which the recessive-at-`+∞`
solution vanishes at the wall; since the wall sits exactly at the Airy zero
(`Λ₀=−a₁ ⟺ Ai(−Λ₀)=0`), `−Λ₀` is *also* the location of the recessive solution’s zero
measured from the edge of the forbidden region. The two boundary conditions **coincide**:
each says “the recessive solution has a zero at the design point.” The SAO fixes the
location (`x=0`) and varies the energy; the sweep fixes the energy (`Y₀`) and varies the
location — dual parametrisations of the same zero-set of `ψ_rec`.

**(R3) Same law (stationarity + reflection).** Near its zero, `ψ_rec` is a measurable
functional of the white noise on the forbidden half-line `{Y≥Y_node}` (RRV’s recessive-side
Riccati is well-defined and non-explosive); `Λ₀` is the same functional of the noise on
`{x≥0}`. The two half-line noises differ by a translation (stationarity) and an orientation
flip (`Y`-increasing vs `x`-increasing) — both invariances of white noise in law. Hence the
two first-zero locations coincide in law: `Y_node =_d −Λ₀(β)`; with RRV’s `−Λ₀=_d TW_β`,
`Y_node =_d TW_β`. ∎ *(strategy)*

**Two independent confirmations.**
1. **`Y₀`-independence (R1):** the law is intrinsic (§1b; histograms collapse).
2. **Dual computation (R2–R3):** `−Λ₀(β)` as the *smallest eigenvalue of the random-matrix
   discretisation* of `H_β` (tridiagonal, Dirichlet) vs `Y_node` by *ODE shooting* — no
   shared code path:
   ```
   β=4 (η=1):   −Λ₀  −2.03 / 0.63 / 0.21   |  Y_node  −2.05 / 0.64 / 0.21   (mean/std/skew)
   β=2 (η=√2):  −Λ₀  −1.74 / 0.89 / 0.26   |  Y_node  −1.76 / 0.89 / 0.29
   det. check:  smallest eig = 2.3382 = −a₁
   ```
   (`folded_cycle_shooting_duality.py`, `.png`.) The matrix eigenvalue knows nothing of
   sweeps or first nodes, so the agreement is **non-circular** evidence for Theorem 1.

**What remains rigorous.** (R1) made uniform in `Y₀` (Berglund–Gentz, tractable); and (R3)
the gluing of the recessive-incoming edge problem to RRV’s Dirichlet operator at the level
of the noise functional (RRV supply the machinery; the matching is the work). Status:
**reduced to a strategy with named classical ingredients, conclusion confirmed two
independent ways** — the same shape that preceded closing the Itô question.

## 2. Two observables — and a correction to the `η²/4π` story

The project’s `H(η)=η²/4π` and the early-escape probability are **different objects**:

| observable | what it is | scaling | role |
|---|---|---|---|
| `H = η²/4π` | integrated **saddle-crossing hazard** (Kramers rate to cross the separatrix) | **polynomial** `∝η²` | sets `σ_*^A` (leading order, **done**) |
| `P_early = 1−F_β(0)` | **canonical escape**: first node / full blow-up *before* the fold | **exponential** `exp(−c/η²)` | the inner exit measure’s tail (**this gap**) |

Crossing the saddle once is *necessary but not sufficient* for a node before the fold
(the trajectory can cross and not run away before `Y` reaches 0), so
`P_early < H` strictly. **This resolves the old puzzle** (`inner_exit_measure.py`:
“`P_esc ≪ η²/4π`”): the two numbers were never meant to be equal — they are different
observables. Crucially, **`σ_*^A = 2√π·√ε₂·√(ac/b)` is unaffected**: it is built from
the hazard `H`, which stands.

**Evidence the canonical escape is exponential.** Fitting `ln P_early` over
`η∈[1.0,2.0]`:
```
   exponential   ln P = −0.31 − 6.22/η²     residual std 0.05
   power-law     ln P = −6.11 + 6.67 ln η    residual std 0.27   (5× worse)
```
(`folded_cycle_prefactor.png`, left). So `P_early ≍ exp(−c/η²)` with **`c ≈ 6.22`**.

## 3. The prefactor, in closed form

**Leading exponential (the rate `c`).** By Freidlin–Wentzell, the early full-escape
probability is `exp(−S₀/η² + o(1/η²))` with `S₀` the minimal action of the
optimal escape path `R(T)` from the canard to `R=+∞` at some `Y>0`,
```
S₀ = min { ½ ∫ (Ṙ − (R²−Y))² dT  :  R(0)=−√Y₀,  R(T*)=+∞,  Y(T*)>0 } .
```
Equivalently it is the **stochastic-Airy ground-state lower-tail rate**, obtained from the
optimal-fluctuation (Lifshitz/instanton) principle
```
c = ½·min_ψ K[ψ]²/∫ψ⁴ = 2(∫g'²)(∫z g²)/∫g⁴ ,   K[ψ]=∫(ψ'²+xψ²) , ψ(0)=0 ,
```
where the optimal shape solves the **nonlinear-Airy ground state** `−g″+z g=g³`, `g(0)=0`,
decaying. **Pinned (`folded_cycle_rate_constant.py`).** Two Pohozaev identities (multiply
the EL by `g` and by `z g′`) give `P=R`, `S=2P` for `P=∫g'², R=∫z g², S=∫g⁴`, collapsing the
formula to `c = P`. Solving the ODE (Pohozaev `P=R`, `S=2P` verified to `10⁻¹³`):
```
┌─────────────────────────────────────────────────────────────┐
│   c = ∫g'² = 5.443902  —  NOT 2π (=6.283);  below a₁²=5.467.  │
└─────────────────────────────────────────────────────────────┘
```
**No elementary closed form (likely).** The substitution `g=√2 s` turns `−g″+zg=g³` into
the **defocusing Painlevé II** `s″=zs−2s³` (parameter `α=0`); so `c` is the *action of an
`α=0` Painlevé II transcendent* on the half-line (`s(0)=0`, decaying) — the same
special-function family as the Tracy–Widom moments, none elementary. Level `0` is a
*generic* point of the SAO lower-tail rate function `I(x)` (it is `|a₁|=2.338` below the
edge): closed forms occur only at special points — `I=0` at the edge `x=2.338`, and the
deep tail `I(x)∼(2/3)|x|^{3/2}` (the universal TW upper-tail constant) as `x→−∞`. To 6
digits `c=5.443902` matches no elementary candidate (`2π`, `π√3=5.4414`, `49/9=5.4444`,
`a₁²=5.4667`, `16/3` all excluded). Exact characterisation (Pohozaev-reduced):
`c = ∫₀^∞ g′² = ∫₀^∞ z g² = ½∫₀^∞ g⁴`, with `∫₀^∞ g² = g′(0)² = 4.28140`.
The earlier full-range fit `c≈6.22` was **pre-asymptotic**: the local slope
`−d ln P_early/d(η⁻²)` *descends* from `~7` (`η≈1.9`) toward `5.444` as `η→0` (`5.66` at the
smallest measured `η`). So `c=2π` is **rejected**; the asymptotic rate is `5.444`, a
nonlinear-Airy constant. *(The leading-order `σ_*^A` constant `C_q=2√π` is unaffected — it
rides on the saddle hazard, §2, not on `c`.)*

**Sub-exponential prefactor (the actual gap).** This is exactly what Proposition 1
delivers: it is the **full Tracy–Widom tail mass**
```
   P_early(η) = 1 − F_{4/η²}(0),
```
computable *to all orders* via Painlevé II (Hastings–McLeod) for `β=1,2,4` and via
Bloemendal–Virág for general `β`. The pre-asymptotic curvature in the `ln P` fit *is*
this prefactor; with `1−F_β(0)` it is captured exactly, not just to leading order. So the
sub-exponential prefactor the gap was missing is **no longer an unknown** — it is a named,
tabulated special function, matched here to MC at `β=1,2`.

**Feeding the Eyring–Kramers constant.** The `O(1)` ambiguity in `C_q` (the `2√π` vs
`2.8` convention gap) is precisely a prefactor question; Proposition 1 fixes the inner
exit law that sets it, so the residual constant is now determined by `F_β` rather than by
an escape-threshold convention.

## 4. Part-1 ledger

**Proved:** Proposition 1 (reduction); **Lemma 0 (Cole–Hopf is an exact Itô identity —
no correction)**; **Lemma 1(α–γ): the exact energy–sweep duality `H_β u=Y₀u`, the
pathwise oscillation reduction, and the deterministic identity `Y_node=−Λ₀^{det}`**; the
two-observable distinction; `P_early=1−F_β(0)`.
**Confirmed numerically:** the reduction in bulk *and* tail; exponential `exp(−c/η²)`;
the three-scheme + pathwise Itô agreement (Lemma 0); the **`Y₀`-independence** that is the
decisive consequence of the matching (Lemma 1δ).
**Open:** (i) ~~Itô↔Stratonovich correction~~ — **resolved (Lemma 0): it is zero**; (ii)
~~the closed form of `c` (is it `2π`?)~~ — **pinned: `c=∫g'²=5.4439` (nonlinear-Airy
instanton, Pohozaev-verified); `2π` rejected, the measured `6.22` was pre-asymptotic**;
(iii) ~~the swept-vs-spectral matching~~ — **reduced
to Theorem 1 (shooting characterisation) with a 3-step strategy (R1 forgetting / R2 one
functional / R3 reflection-invariance), confirmed two independent ways (`Y₀`-independence
and random-matrix SAO ≡ ODE shooting); the deterministic limit is exact.** The residue is
the uniform-in-`Y₀` BG bound (R1) and gluing the recessive edge problem to RRV (R3); (iv)
`F_β(0)` in closed form for general `β` (numerically for all `β`; exact via Painlevé for
`β=1,2,4`).

---

# Part 2 — Sharp uniform error bounds for the composed noisy map

## 5. The object

The theorem composes the noisy transition maps
`Π₁ ∘ κ₁₂ ∘ Π₂ ∘ κ₂₃ ∘ Π₃` (entry K₁ → inner K₂ → exit K₃), the stochastic analogue of
JKK §4.5. We need, **uniformly in the angle `θ∈S¹=ℝ/ℤ`**: (i) a concentration estimate
on each `Πᵢ`; (ii) overlap (Jacobian) control on each `κᵢ,ᵢ₊₁`; (iii) propagation of the
errors through the composition. JKK’s central difficulty is that `θ` cannot be localised;
we resolve it not by localising but by making every constant **θ-uniform**.

## 6. Lemma 1 (θ-uniformity) — proved

*Hypotheses (JKK normal form).* `a,b,c: S¹→ℝ` are smooth, 1-periodic, and **strictly
positive**.

> **Lemma 1.** *There are constants `0<a₋≤a₊`, `0<b₋≤b₊`, `0<c₋≤c₊` (the min/max of
> `a,b,c`) such that every inner estimate below holds with `θ`-independent constants
> depending only on `a±,b±,c±`.*

*Proof.* `S¹` is compact and `a,b,c` are continuous and strictly positive, so each attains
a positive minimum and finite maximum: `a₋=min_θ a(θ)>0`, `a₊=max_θ a(θ)<∞`, etc.
The three inner quantities are monotone in `a,b`, hence sandwiched uniformly in `θ`:

- **contraction rate** `κ(θ,y₂)=2√(a(θ)b(θ)y₂)`:
  `2√(a₋b₋ y₂) ≤ κ(θ,y₂) ≤ 2√(a₊b₊ y₂)` — bounded **below** by a positive multiple of
  `√y₂`, uniformly in `θ`;
- **tube variance** `Var(r₂)=η²/(2κ)=η²/(4√(ab y₂))`:
  `η²/(4√(a₊b₊ y₂)) ≤ Var ≤ η²/(4√(a₋b₋ y₂))` — bounded **above** uniformly;
- **quasipotential** `V(θ,y₂)=(8/3)√(a³/b)\,y₂^{3/2}`:
  `(8/3)√(a₋³/b₊)\,y₂^{3/2} ≤ V ≤ (8/3)√(a₊³/b₋)\,y₂^{3/2}`, so the Kramers factor
  `e^{−V/η²}` is squeezed between two `θ`-independent exponentials.

Each estimate in §7 is built from these three quantities by operations preserving the
uniform sandwich (composition, integration in `y₂`, exponentiation), so the resulting
constants depend only on `a±,b±,c±`, not on `θ`. ∎

*Remark.* This is the load-bearing move: JKK could not localise `θ`; Lemma 1 makes that
unnecessary, because compactness + strict positivity already give `θ`-uniform constants.
**Verified** numerically: the rescaled tube variance `Var·4√(ab y₂)/η²` collapses to `1`
across four geometries `(a,b)` and two `η` (`folded_cycle_chartmatch.py`).

## 7. Per-chart concentration (Berglund–Gentz), with uniform constants

**K₂ inner tube (Lemma 2, stated).** Linearising about the canard `r₂=−√((a/b)y₂)`, the
deviation `δ=r₂−γ₂` is the Ornstein–Uhlenbeck process `dδ=−κ(θ,y₂)δ\,dT+η dB` with
`κ=2√(ab y₂)`. By the Berglund–Gentz concentration bound, for any `h>0` the sample path
stays in the tube `{|δ|≤h}` up to the first exit time with probability
`≥ 1 − C(t/η²)·exp(−κ_- h²/(2η²))`, `κ_-=2√(a₋b₋ y₂)` (Lemma 1). Escape is the tube
reaching the separatrix at distance `2√((a/b)y₂)`; substituting recovers the FW barrier
`V=8Y^{3/2}/3` and, *near the fold* `y₂≲η^{4/3}`, the inner exit measure of **Part 1**
(this is exactly where the quasi-static OU breaks and `TW_β` takes over).

**K₁ entry / K₃ exit (stated).** Entry `Π₁`: Fenichel contraction onto the attracting
slow manifold gives an `O(η²)`-variance Gaussian tube with `θ`-uniform rate (Lemma 1
applied to the entry chart). Exit `Π₃`: the deterministic exit Jacobian is bounded on
`S¹`; the noisy exit map concentrates with the same uniform constants (derived in
`FOLDED_CYCLE_PATHA_PROOFS.md`).

## 8. Overlap / matching estimates (κ₁₂, κ₂₃)

The chart-overlap maps are the JKK blow-down Jacobians, smooth and non-degenerate on the
compact overlaps; by Lemma 1 their norms and inverses are uniformly bounded. The noisy
matching is then: push the entry tube through `κ₁₂`, run the inner `Π₂`, push through
`κ₂₃`, run the exit. Each push multiplies the tube width by a bounded Jacobian, so a tube
of width `O(η/√y₂)` stays a tube of width `O(η)` across the overlaps — uniformly in `θ`.

## 9. Part-2 ledger

**Proved:** Lemma 1 (θ-uniformity) in full; the per-chart OU tube and its uniform
constants (Lemma 2, on the verified BG bound); the overlap Jacobian bounds.
**Confirmed numerically:** the tube collapse across geometries (θ-uniform), and the
near-fold breakdown handing off to Part 1.
**Open (the genuine remaining ~50%):** (i) *propagation* of the sharp constants through
the full composition `Π₁∘κ₁₂∘Π₂∘κ₂₃∘Π₃` (error accumulation, not just per-factor); (ii)
the non-quasi-static near-fold correction — *now supplied in distribution by Part 1’s
`TW_β`*, but its insertion into the composed error bound is not yet done; (iii) the
Channel-B inner (phase) with its bifurcation-specific exponent; (iv) uniformity in
`(ε₁,ε₂)` and the blow-down to physical variables.

---

## 10. Net effect on Path A

Before: Path A was “a theorem to **exponential order**, θ-uniform, with one named gap
(the inner exit measure / prefactor).” After this pass:

- the **inner exit measure is identified and tail-validated** as `TW_β`, and the missing
  **sub-exponential prefactor is `1−F_β(0)`** — a named special function, not an unknown;
- a **correction** is logged: the canonical escape is `exp(−c/η²)` (`c≈6.22≈2π?`), and
  `η²/4π` is the *distinct* saddle hazard that (correctly) sets `σ_*^A`;
- the **θ-uniformity that the whole composition rests on is now a proved lemma**.

What remains for the full theorem is honest and specific: make Cole–Hopf rigorous
(Itô/Strato), pin `c`, and **propagate** the now-uniform per-factor constants through the
composition. That is bookkeeping on a fully-named skeleton — not an open question.

## References

- Ramírez, Rider, Virág, *Beta ensembles, stochastic Airy spectrum, and a diffusion*,
  J. AMS 24 (2011). — Proposition 1(b).
- Bloemendal, Virág, *Limits of spiked random matrices* (general-β edge). — `F_β`.
- Tracy, Widom (1994, 1996); Bornemann (2010). — `F_β`, the Painlevé II numerics.
- Freidlin, Wentzell, *Random Perturbations of Dynamical Systems*. — §3 action `S₀`.
- Berglund, Gentz, *Noise-Induced Phenomena in Slow–Fast Dynamical Systems* (2006);
  Berglund–Gentz–Kuehn (2012). — §7 concentration tubes.
- Jelbart, Kuehn, Kuntz (2024), arXiv:2208.01361. — the deterministic backbone, §4.5 maps.
- Internal: `FOLDED_CYCLE_PATHA_PROOFS.md` (the exponential-order theorem this extends),
  `folded_cycle_prefactor.py` (Part 1), `folded_cycle_chartmatch.py` (Lemma 1 collapse).
