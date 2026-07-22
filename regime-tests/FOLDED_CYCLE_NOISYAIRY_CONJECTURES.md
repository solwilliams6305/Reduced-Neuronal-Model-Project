# The noisy-Airy gap — creative approaches (conjectures)

**Status:** Path A bucket 1 (`FOLDED_CYCLE_PATHA_PROOFS.md`) proved the Channel-A
critical-noise law `σ_*^A = 2√π·√ε₂·√(ac/b)` **to exponential order, uniformly in
θ**. The single named gap is the **inner exit measure**: the law of where the
noisy canard peels off in the noisy-Airy region `Y ≲ η^{4/3}`, which sets the
**sub-exponential prefactor**. This document speculates how to close it.

**The headline:** it almost certainly does **not** need a brand-new idea — it needs
a **cross-field connection that already exists**. The inner problem is *exactly
linear* (Cole–Hopf), and that linear object is the **stochastic Airy operator**,
whose edge is **Tracy–Widom / KPZ universality**. Conjecture 1 is now **numerically
confirmed**: the inner exit measure matches Tracy–Widom₍β₎ with `β = 4/η²`,
*parameter-free* in mean, variance, skewness and kurtosis at `β = 1, 2` — and the
**raw** histograms sit on the **exact** `TW₁`/`TW₂` densities (Painlevé II) with no
fit. The remaining work is turning the identification into the rigorous prefactor.

Companion: `folded_cycle_noisyairy.py` (bulk `O(η)`), `folded_cycle_tracy_widom.py`
(skewness vs `β`), `folded_cycle_tw_pdf.py` (exact-pdf overlay);
`figures/folded_cycle_{noisyairy,tracy_widom,tw_pdf}.png`.

---

## UPDATE (2026-06-07) — C2 closed, tail reconciled, prefactor now a number

**C2 (Frisch–Lloyd) — DONE.** The inner escape current is computed exactly (MFPT /
log-sum-exp quadrature) and confirmed by Monte-Carlo to a few percent across four decades.
Under the exact inner rescaling `R=η^{2/3}ρ, Y=η^{4/3}y` the frozen-`Y` SDE is
`η`-independent (`dρ=(ρ²−y)dτ+dW`), so there is ONE universal current `𝒥(y)`, with the
physical rate `λ(Y,η)=η^{2/3}𝒥(Y/η^{4/3})`. Findings: `𝒥(y)/𝒥_qs(y)→1` as `y→∞`
(Kramers recovered, ratio 0.98 at `y=3`); Kramers **over**estimates by ~18% near `y≈0.6`;
and `𝒥(0)=0.158>0` is finite where quasi-static Kramers gives exactly **zero** — escape
persists at and past the fold. Scripts `frisch_lloyd_current.py`,
`figures/frisch_lloyd_current.png`.

**Closing H (the outer matching, the old C3 worry) — DONE.** The early-escape integrated
hazard `H=∫_0^∞ λ dY = η²·∫_0^∞ 𝒥(y)dy` is finite: `𝒥` matches the outer quasi-static rate
at large `y`, and the divergent `y→−∞` branch (`𝒥~√|y|`) is the *certain late escape past
the fold*, correctly excluded by the early-escape (`Y>0`) definition — so there is no
divergent outer piece; the inner integral IS the uniformly-valid answer. Numerically
`∫_0^∞ 𝒥_qs = 0.07956 ≈ 1/4π = 0.07958` (check) and

> **H = K·η²/4π  with  K = 1.039** (a +3.9% prefactor; the finite rate at the fold slightly
> outweighs the mid-range deficit). The sub-exponential prefactor on the integrated hazard
> is now a closed number. Script `close_H_matching.py`, `figures/close_H_matching.png`.

**Tail reconciliation (falsifier 2) — DONE.** `P_early(η)=P(Y_node>0)` measured by the
recessive Cole–Hopf sweep (integrate `u''=(Y−ηξ)u` on the *repelling* canard branch, first
node of `u`). `ln P_early` is **perfectly linear in `1/η²`** (`R²=1.00`) — the exponential
SAO/TW left tail — while the polynomial hazard `η²/4π` has the wrong shape and over-counts
`P_early` by 1.9× (`η=2`) up to 59× (`η=1`): confirming `P_early < H` (a separatrix crossing
is necessary but not sufficient for a committed node). The MC lands on the **exact** TW
values: `β=2 (η=√2): 0.0312` vs `0.0306`; `β=1 (η=2): 0.165` vs `0.168`. Effective tail
constant over `η∈[1,2]` is `c≈6.35`, sitting just above and trending toward the asymptotic
Painlevé II value `c=5.4439` as `η→0`. (Methodological note: the canard is the recessive /
repelling `+√Y` branch — seeding the *attracting* `−√Y` branch gives the wrong anchor
`−2.13`; the recessive sweep nails the first Airy zero `−2.338`.) Scripts
`tail_reconciliation.py`, `figures/tail_reconciliation.png`.

---

## 0. The gap, precisely

The inner (K₂) Riccati at `λ=1` is `dR = (R²−Y) dT + η dB`, `Y = Y₀−T`. The
deterministic canard escapes at the first Airy zero `Y = a₁ ≈ −2.338` (`Ω₀`).
The **quasi-static** integrated hazard gives the leading exponential and the
closed form `H(η) = η²/4π`. But near the fold the barrier `8Y^{3/2}/3` **vanishes**,
so the quasi-static OU approximation breaks in the window `Y ≲ η^{4/3}`
(the inner scaling `R_* ∼ η^{2/3}`, `Y_* ∼ η^{4/3}`). There, `H = η²/4π`
**overestimates** the true early-escape fraction (measured in
`inner_exit_measure.py`). The missing object is the **inner exit measure**
`μ_η(Y_peel, R_peel)` — equivalently the **non-quasi-static prefactor** that
corrects `η²/4π`. That is the whole remaining gap.

## 1. The unifying observation — the inner problem is EXACTLY linear

The noisy Riccati Cole–Hopf-linearises. With `R = −u′/u`,

```
dR = (R²−Y)dT + η dB     ⇔     u″ = (Y − η ξ) u ,      ξ = dB/dT.
```

(Check: `R′ = −u″/u + (u′/u)² = −(Y−ηξ) + R² = R²−Y+ηξ`. Stratonovich reading;
see caveat C-i.) **Escape `R→+∞` ⇔ first node of `u`.** So the inner exit measure
is the **first-node law of a linear stochastic ODE** — a Schrödinger equation with
the deterministic Airy potential `Y(T)` plus a **white-noise potential** `−ηξ(T)`.

**Verified (Cole–Hopf, earlier):** escape `Y` (Riccati, `R>20`) `= −1.625` vs first
zero of `u` `= −1.675`, and `R` and `−u′/u` track along the path. **MATCH.**

This is the pivot. Everything below is a different way to read the *same* linear
operator. The reason to be optimistic: a 1D Schrödinger operator with a
deterministic ramp `+` white noise is one of the **most-studied objects in
mathematical physics** (Anderson localization, random matrices, KPZ). We are not
short of an idea; we are short of the right *dictionary entry*.

---

## 2. Conjecture 1 (strongest) — the stochastic Airy operator → Tracy–Widom edge

`u″ = (Y − ηξ)u` **is** the eigenvalue equation of the **stochastic Airy operator**

```
H_β = −d²/dx² + x + (2/√β)·b′(x),        dictionary:  x = −Y,   η = 2/√β   (β = 4/η²).
```

Sweeping `Y` downward = sweeping the spectral level; the **first node of `u`** is the
**ground-state eigenvalue** `λ₀(β)` of `H_β`. Ramírez–Rider–Virág (2011) prove
exactly this: the soft edge of β-ensembles converges to the spectrum of `H_β`,
characterised by a **Riccati diffusion** `p′ = x − λ − p² + (2/√β)ξ` — which, under
`p ↔ −R`, `x−λ ↔ Y`, is **our inner Riccati**. The conjecture:

```
┌────────────────────────────────────────────────────────────────────────────┐
│  Inner exit measure μ_η  =  ground-state law of the stochastic Airy operator │
│  H_{β=4/η²}.  Two regimes:                                                    │
│    BULK (typical peel):   Y_node = a₁ + O(η)·Gaussian   [large-β SAO edge]    │
│    TAIL (early escape):   left large-deviation of λ₀ → the prefactor on η²/4π │
└────────────────────────────────────────────────────────────────────────────┘
```

**Grounding — bulk (NEW, this script).** First-order perturbation theory predicts
`Y_node = a₁ + η·∫ψ₀² dW + O(η²)`, Gaussian with std `η·√(∫ψ₀⁴)`. Measured:
std `~ η^{1.00}` (fit `p = 1.001`), constant `std/η = 0.644 = √(∫ψ₀⁴)` (an explicit
Airy integral), mean `→ −2.336 → a₁ = −2.3381`. The bulk is **exactly the SAO
perturbative edge.** (Figure, left + histogram.)

**Grounding — full shape (NEW, `folded_cycle_tracy_widom.py`).** The standardized
skewness/kurtosis are affine invariants, so they fingerprint the law with **no
centring/scaling freedom**. Measured against `TW_β` at `β = 4/η²`:

```
   β=4/η²   η      meas (mean, std, skew, kurt)        Tracy–Widom_β
   4.00   1.000    -2.06  0.640  0.173  ~0   (shape ✓; 2^{-1/6} GSE-convention scale)
   2.00   1.414    -1.778 0.904  0.244  0.13   TW2: -1.771 0.902 0.224 0.093
   1.00   2.000    -1.209 1.273  0.306  0.16   TW1: -1.207 1.268 0.293 0.165
```

At `β = 1, 2` the match is **parameter-free in all four moments** (well-resolved,
`dt=5e-4`, `N=3·10⁴`). The skewness rises from `0` (Gaussian, `β→∞`) and tracks
`TW_β` along `η = 2/√β` (figure, left); the standardized law is visibly right-skewed
— the heavy tail is the **early-escape** side (figure, right). `β=4` (GSE) matches in
*shape* (skew `0.17` vs `0.166`); its location/scale carry the well-known `2^{-1/6}`
GSE edge-convention factor, not a failure of the identification.

**This confirms Conjecture 1:** the inner exit measure **is** the stochastic-Airy /
Tracy–Widom edge law, dictionary `η = 2/√β`. The early-escape tail (`Y_node > 0` ⇔
`λ₀ < 0`) is the `TW_β` left large-deviation `∼ exp(−βs³/24) = exp(−s³/6η²)`
(Dumaz–Virág); the **true** prefactor on the quasi-static `H = η²/4π` is the SAO
ground-state left-deviation density. The remaining analytic task is to *read off*
that prefactor — now known to be a `TW_β` edge quantity, not an open unknown.

**Why this closes the gap (not just renames it):** RRV give a *rigorous* diffusion
characterisation of `λ₀(β)` — our Riccati — and the edge density is a *known*
special function (Painlevé II / Hastings–McLeod for `β=2`, the `TW_β` family
generally). The prefactor we need is therefore *already computed in the literature*,
once the dictionary `η = 2/√β` and the sweep↔spectrum identification are nailed.

**Caveats (technical, not conceptual).**
- **C-i (Itô/Stratonovich) — RESOLVED.** It looked as if `R=−u′/u` (ordinary chain rule
  ⇒ Stratonovich) would need an Itô correction against the RRV operator. It does **not**:
  the correction is identically zero (Lemma 0, `FOLDED_CYCLE_NOISYAIRY_PROOF.md`) —
  additive Riccati noise, `(g·∇)g=0` for the linear system, and `R_{vv}=0` in Itô’s
  formula all vanish. Verified pathwise (`RMS|R−(−v/u)|∝dt→0`). The `O(η²)` mean drift is
  the genuine second-order shift of the stochastic-Airy ground state, not a calculus
  artifact.
- **C-ii (boundary/sweep).** RRV is the half-line operator with a fixed soft edge;
  ours has the *moving* level `Y(T)=Y₀−T` and a finite passage. Identifying
  "first node during the sweep" with "lowest eigenvalue" needs a matching argument
  (the entry from `Y₀≫1` supplies the decaying Airy initial condition `u∼Ai`).
- **C-iii (multiplicative vs additive).** SAO has the potential *additive* in `x`;
  our noise multiplies `u`. They coincide for the *spectral* (eigenvalue) question
  but differ for pathwise statements — fine here, since we only want the exit law.

## 3. Conjecture 2 — Frisch–Lloyd / stationary Riccati density (1D localization)

Read `R(T)` as a diffusion and ask for its **quasi-stationary density with an
absorbing escape current**. The stationary Fokker–Planck for `dR=(R²−Y)dT+η dB`
at frozen `Y`,

```
(η²/2) ρ″ + ∂_R[(R²−Y)ρ] = −J  (const probability current = escape rate),
```

is the **Frisch–Lloyd equation** of 1D Anderson localization: the escape rate `J`
is the **integrated density of states / rotation number** of `u″=(Y−ηξ)u`, and the
Riccati blow-up *counts eigenvalues*. The classical Frisch–Lloyd (1960) / Halperin
current formula gives `J` in closed quadrature (an Airy-function integral at the
relevant `Y`). The prefactor on the escape rate is then the Frisch–Lloyd current,
**not** the quasi-static `√Y/π`.

**Promise:** the most *elementary and directly computable* route — it is the
standard tool for exactly "count the blow-ups of a Riccati with white-noise
forcing," and it does not require the full RRV machinery. **Falsifier:** integrate
the stationary current for `u″=(Y−ηξ)u` at fixed `Y` and compare to the Kramers
rate `(√Y/π)e^{−8Y^{3/2}/3η²}`; the *deviation* near `Y∼η^{4/3}` is the inner
correction we want.

## 4. Conjecture 3 — exact backward-Kolmogorov committor in special functions

The **committor** `q(R,Y) = P(escape before the fold | R,Y)` solves

```
(η²/2) q_RR + (R²−Y) q_R − q_Y = 0,      q→1 (R→+∞),  q→0 (R→−∞ branch).
```

This is a parabolic PDE: a heat term `(η²/2)q_RR` over the cubic-potential drift.
**Conjecture:** in the inner scaling `R = η^{2/3} ρ`, `Y = η^{4/3} y`, `T = η^{... }τ`
the PDE collapses to an **η-independent canonical PDE**, whose boundary flux **is**
the universal inner exit measure. The deterministic part is Airy; with the heat
term, seek similarity / separated solutions in **parabolic-cylinder or Airy
functions**. If the similarity reduction exists, this gives the *cleanest possible*
closed form — a single canonical profile with `(Y_peel,R_peel)=(η^{4/3},η^{2/3})×`
(universal law), which is precisely the scaling `inner_exit_measure.py` already saw
(`⟨Y_peel⟩ ∼ η^{4/3}`, fitted exponent ≈ 4/3).

**Promise:** best if the similarity reduction works — it would make the prefactor a
*number* (a quadrature of the canonical profile). **Risk:** the reduction may fail
to fully decouple `η` because of the matching to the outer quasi-static region; then
this becomes a boundary-layer (matched-asymptotics) problem, not a clean similarity.

## 5. Conjecture 4 — instanton + Gelfand–Yaglom fluctuation determinant

The early-escape tail is a large deviation: leading exponential = the FW action
(reproducing `η²/4π`), prefactor = the **one-loop fluctuation determinant** about the
**instanton** (optimal escape path), computed by **Gelfand–Yaglom**. Two features
make this both standard and interesting:

- **Near-zero mode.** Time-translation along the sweep gives an approximate zero mode
  of the fluctuation operator; it must be treated collectively (a Faddeev–Popov /
  collective-coordinate factor), and it typically *supplies* the `√` and the `η`-power
  in the prefactor.
- **Vanishing-barrier caustic.** The barrier `8Y^{3/2}/3 → 0` at the fold, so as the
  escape time approaches the fold the instanton **merges with the saddle**
  (instanton–anti-instanton / caustic). A naive Gelfand–Yaglom determinant
  **diverges** there; the cure is a **uniform (Airy-type) treatment** of the
  determinant near the caustic — which *loops straight back to Conjecture 1's Airy
  edge structure.*

**Promise:** the most physically transparent and the standard way to get prefactors.
**Risk/interest:** the vanishing-barrier caustic is the crux, and the uniform
treatment there is *exactly* the stochastic-Airy edge — i.e. C4 and C1 are the same
mountain from two sides.

---

## 6. Ranking, and do we need a new idea?

```
  C1  SAO → Tracy–Widom edge       payoff ●●●●●  grounding ●●●●●  (TW match confirmed, β=1,2)
  C2  Frisch–Lloyd current         payoff ●●●○○  grounding ●●●●●  (DONE: K=1.039 on H; see UPDATE)
  C3  backward-Kolmogorov similar. payoff ●●●●○  grounding ●●○○○  (cleanest if it reduces)
  C4  instanton + Gelfand–Yaglom   payoff ●●●○○  grounding ●●○○○  (caustic = C1)
```

**Recommended attack:** run **C2 and C1 together.** Frisch–Lloyd gives a *computable*
current (the prefactor as a quadrature) fast; the SAO/TW identification tells you
*what universal object* that current is and supplies the rigorous backbone (RRV's
diffusion = our Riccati). C3 is the high-reward side-bet (try the similarity
reduction — it is one afternoon to find out if `η` decouples). C4 is the physicist's
sanity check and, via its caustic, *certifies* that C1 is unavoidable.

**Do we need a NEW idea?** **Probably not.** The gap yields to a **connection**, not
an invention: the inner exit measure is (conjecturally) the **stochastic Airy
operator's ground-state law**, with `η = 2/√β`; the prefactor is its **edge density**;
and three independent classical tools (RRV diffusion, Frisch–Lloyd current,
Gelfand–Yaglom determinant) all target it. The *genuinely* novel content — if any is
needed — is the **uniform treatment at the vanishing-barrier caustic `Y→0`**, where
instanton-merging meets the Airy edge. That precise junction (a fold-canard escape
read as a Tracy–Widom edge) does not appear in the literature in this form, and is
the one place a new idea might be required.

## 7. Numerical falsifiers (next, no scipy needed)

1. **TW shape — DONE, pointwise (`folded_cycle_tracy_widom.py`, `folded_cycle_tw_pdf.py`).**
   Skewness/kurtosis match `TW_β` at `β = 4/η²`; parameter-free in mean+std+skew+kurt at
   `β = 1, 2`. The exact `TW₁` and `TW₂` densities (Painlevé II / Hastings–McLeod,
   bare-numpy RK4, validated to 4 dp against tabulated moments) overlay the **raw**
   simulated histograms with **no centring, scaling, or fit** — see
   `figures/folded_cycle_tw_pdf.png`. ✓✓
2. **Tail reconciliation — DONE (see UPDATE).** `ln P_early` is linear in `1/η²` (`R²=1.00`,
   the exponential SAO tail, `c≈6.35→5.444` as `η→0`), NOT the polynomial hazard `η²/4π`
   (which over-counts up to 59×); MC matches exact TW at `β=1,2`. `tail_reconciliation.py`. ✓✓
3. **Frisch–Lloyd — DONE (see UPDATE).** Universal inner current `𝒥(y)` computed (quadrature
   = MC); matches `(√Y/π)e^{−8Y^{3/2}/3η²}` for `y≫1`, departs at `y=O(1)`, finite at the fold;
   integrated to `H=1.039·η²/4π`. `frisch_lloyd_current.py`, `close_H_matching.py`. ✓✓

## 8. What this buys Path A

Closing this gap upgrades `FOLDED_CYCLE_PATHA_PROOFS.md` from "critical-noise law to
**exponential** order, θ-uniform" to "**with the sharp sub-exponential prefactor**"
— i.e. the full Eyring–Kramers constant `A₀` for the folded-cycle canard escape,
and with it the *complete* leading-order-plus-prefactor inner exit measure that the
global chart composition `Π₁∘κ₁₂∘Π₂∘κ₂₃∘Π₃` needs. And since C1 is **numerically
confirmed**, it exhibits the folded-limit-cycle canard escape as a member of the
**Tracy–Widom / KPZ edge universality class** — a result of independent interest
beyond the JKK program, and (to our knowledge) the first time a fold-canard exit law
has been identified with the stochastic Airy spectrum.

## 9. References

- **Ramírez, Rider, Virág**, *Beta ensembles, stochastic Airy spectrum, and a
  diffusion*, J. AMS 24 (2011) — the stochastic Airy operator and its Riccati
  characterisation (= our inner Riccati; the `2/√β` dictionary).
- **Dumaz, Virág**, *The right tail of the largest eigenvalue / SAO ground-state
  tails* — the `exp(−βs³/24)` left tail used in §2.
- **Frisch, Lloyd** (1960); **Halperin** (1965) — 1D Schrödinger with white-noise
  potential, the IDOS/rotation-number current (C2).
- **Gelfand, Yaglom** — functional determinants for the one-loop prefactor (C4);
  instanton/caustic uniformization.
- **Tracy, Widom** — the `TW_β` edge laws; Painlevé II / Hastings–McLeod.
- JKK 2024 §4.3 (Airy inner, `Ω₀`); Krupa–Szmolyan 2001 (fold blow-up).
- Internal: `FOLDED_CYCLE_PATHA_PROOFS.md` (the proof this completes),
  `inner_exit_measure.py` (the `η^{4/3}` inner scale), `folded_cycle_noisyairy.py`
  (the bulk SAO verification here).
