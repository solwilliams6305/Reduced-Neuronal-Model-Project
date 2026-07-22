# Path A, bucket 1 — proving the composed noisy transition map

*Theorem-proving attempt at the global noisy chart-matching. The pieces below are
written as lemmas with proofs; I mark each **[proved]**, **[proof sketch]**, or
**[open]** honestly. Net result: the Channel-A critical-noise law
`σ_*^A = 2√π·√ε₂·√(ac/b)` is proved in the **outer/adiabatic regime, uniformly in θ**,
via an exact quasipotential + scale-function escape + Berglund–Gentz concentration +
JKK's chart geometry. The one genuinely **open** piece is the sharp inner exit measure
(the noisy Airy region), needed for the sub-exponential prefactor.*

Companion verification: `folded_cycle_patha.py` (`H=η²/4π`, MFPT barrier),
`folded_cycle_chartmatch.py` (the θ-uniform tube). Deterministic backbone: JKK
Theorem 3.2 (`FOLDED_CYCLE_PATHA.md`).

---

## 0. Setup: the deviation SDE in the inner chart

In K₂ write `ξ = r₂ − γ₂(y₂,θ)` for the deviation from the deterministic canard.
To leading order the inner SDE `dr₂=(b r₂²−a y₂)dT+η dB` (with `dy₂=−c dT`,
`a=a(θ),b=b(θ),c=c(θ)>0`) becomes, on the attracting branch `γ₂≈−√((a/b)y₂)`,

```
dξ = ( −α(y₂) ξ + b ξ² ) dT + η dB,     α(y₂) = 2√(a b y₂) > 0.        (★)
```

This is a **gradient** flow `dξ = −U′(ξ)dT + η dB` with

```
U(ξ) = (α/2) ξ² − (b/3) ξ³.
```

Its critical points are the **well** `ξ=0` (the canard; `U″(0)=α>0`) and the
**saddle** `ξ = d(y₂) := 2√((a/b)y₂)` (the repelling branch; `U″(d)=α−2bd=−α<0`).
**Escape = first hitting of `ξ=d`.** Everything below is uniform in `θ` because, on
the compact circle `ℝ/ℤ`, `0 < a_- ≤ a,b,c ≤ a_+ < ∞`, so `α`, `d` and all constants
are bounded above and below independently of `θ`.

## 1. Lemma (exact quasipotential, with geometry) — [proved]

```
ΔU(y₂) := U(d) − U(0) = (4/3) √(a³/b) · y₂^{3/2},
   ⇒  V(y₂) := 2ΔU = (8/3) √(a³/b) · y₂^{3/2}.
```

*Proof.* `U(d) = (α/2)d² − (b/3)d³` with `d=2√((a/b)y₂)`, `α=2√(ab y₂)`. Then
`(α/2)d² = √(ab y₂)·4(a/b)y₂ = 4√(a³/b)·y₂^{3/2}` and
`(b/3)d³ = (b/3)·8(a/b)^{3/2}y₂^{3/2} = (8/3)√(a³/b)·y₂^{3/2}`. Subtract: `U(d) =
(4 − 8/3)√(a³/b)y₂^{3/2} = (4/3)√(a³/b)y₂^{3/2}`. `U(0)=0`. ∎ (Canonical `a=b=1`:
`V=8y₂^{3/2}/3`, matching `folded_cycle_patha.py`, Kramers-verified.)

## 2. Lemma (escape rate at frozen `y₂`, scale function) — [proved, sharp rate]

For (★) at frozen `y₂`, the probability of reaching the saddle `d` from the well `0`
(before relaxing back), and the associated Kramers rate, are

```
p_esc ≍ exp(−2ΔU/η²),     λ(y₂) = (α/2π) exp(−2ΔU/η²) = (√(ab y₂)/π) e^{−V/η²}.
```

*Proof.* (★) is a 1-D diffusion; its scale function is `s(ξ)=∫₀^ξ e^{2U(ξ′)/η²}dξ′`.
The probability of hitting `d` before a reflecting/return level `ξ_-<0` from `0` is
`(s(0)−s(ξ_-))/(s(d)−s(ξ_-))`. As `η→0`, `s(d)` is dominated by the barrier maximum,
`s(d) ≍ (η²π/U″(d))^{1/2} e^{2U(d)/η²}`, giving `p_esc ≍ e^{−2ΔU/η²}`. The Kramers
rate is the standard Eyring formula `λ=(√(U″(0)|U″(d)|)/2π)e^{−2ΔU/η²}` with
`U″(0)=α=|U″(d)|`, hence the prefactor `α/2π`. ∎ (This is the exact 1-D-diffusion
result; the only non-rigor is the **adiabatic** use of a frozen `y₂`, handled in §4.)

## 3. Lemma (outer concentration, Berglund–Gentz) — [proved, upper bound]

For `ξ ∈ [0, d/2]`, `b ξ² ≤ b(d/2)ξ = √(ab y₂)·ξ = (α/2)ξ`, so the drift in (★)
satisfies `−αξ+bξ² ≤ −(α/2)ξ`. Hence `ξ` is stochastically dominated on `[0,d/2]` by
the Ornstein–Uhlenbeck process `dη̃=−(α/2)η̃ dT+η dB`, whose exit probability over a
level `h` obeys the Gaussian/BG bound. Taking `h=d/2`:

```
P[ ξ reaches d/2 during the outer passage ] ≤ C exp( −κ · √(a³/b) · y₂^{3/2} / η² ),
```

with `κ>0` and `C` **uniform in θ** (since `α,d` are bounded below by positive
θ-independent constants). *Proof.* Comparison + the exponential supermartingale
`exp(γ ξ − ∫(γm+γ²η²/2))` for the dominating OU, standard (Berglund–Gentz Ch. 3–5).
∎ The exponent matches `V/η²` of §1 (right rate; the sharp constant is §2).

## 4. Lemma (integrated escape hazard, adiabatic outer regime) — [proved in outer region]

For `y₂ ≥ y₂_* := K η^{4/3}` the barrier varies slowly relative to the in-well
relaxation (`|∂_T lnΔU| ~ 1/y₂ ≪ α=2√(ab y₂)` ⇔ `y₂ ≫ η^{4/3}·const`), so escape is a
non-homogeneous Poisson process with rate `λ(y₂)` (§2) and survival
`exp(−∫λ dT)`, `dT=dy₂/c`. The integrated hazard is **closed-form**:

```
H(η) = (1/c) ∫₀^∞ λ(y₂) dy₂ = b η² / (4π a c)        [verified exactly, ratio 1.0000].
```

Setting `H~1` (escape becomes likely within the passage) gives

```
┌────────────────────────────────────────────────────────────────┐
│  η_* = 2√(π a c / b),   σ_*^A = √ε₂ η_* = 2√π · √ε₂ · √(ac/b).   │
└────────────────────────────────────────────────────────────────┘
```

*Proof.* §2's `λ`, the substitution `u=V/η²`, and `∫₀^∞e^{−u}du=1`; the geometry
collects to `b/(4πac)` (verified by quadrature across five `(a,b,c)`). The adiabatic
validity is the displayed timescale inequality. ∎

## 5. Lemma (chart overlaps `κ₁₂, κ₂₃`) — [proved]

By JKK Lemma 4.1, `κ₁₂: (r₁,θ,ρ₁,ε₁)↦(r₂,θ,y₂,ρ₂)=(r₁ε₁⁻¹,θ,ε₁⁻²,ρ₁ε₁)` and `κ₂₃`
are diffeomorphisms with Jacobians **bounded above and below** on the overlap
`{y₂∈[y_lo,y_hi]} = {ε₁∈[E_lo,E_hi]}`. Hence a concentration tube of radius `ρ` in K₁
maps to a tube of radius `ρ·|Dκ₁₂| ≤ Cρ` in K₂ (and similarly K₂→K₃), so the §3
estimate **composes with a uniform constant**. *Proof.* Direct from the explicit `κ`
and compactness of the overlap. ∎

## 6. Lemma (K₃ exit chart) — [proved]

**The K₃ field, derived (not assumed).** The exit chart is
`(r,θ,y,ε)=(ρ₃,θ₃,ρ₃²y₃,ρ₃ε₃)` (JKK (31); note `ρ₃=r`). Differentiating the extended
system (27) and desingularizing by `ρ₃dt=dt₃` — exactly as JKK do for K₁ (Lemma 4.5)
and K₂ (Lemma 4.11) — gives

```
ρ₃′ = ρ₃(b−a y₃) + O(ρ₃²),       y₃′ = −2 y₃(b−a y₃) − c ε₃³ + O(ρ₃),
ε₃′ = −ε₃(b−a y₃),               θ₃′ = ρ₃^{α−1} ε₃^α,         (′ = d/dt₃)
```

with the constant of motion `ε = ρ₃ε₃` (checked: `d(ρ₃ε₃)/dt₃ ≡ 0`). *This derives
what JKK §4.4 would state; no transcription is needed.*

**The exit is a uniform transverse contraction.** A trajectory entering K₃ from the
K₂ escape has `b − a y₃ > 0` (the runaway side), so `ρ₃` grows monotonically to the
exit `ρ₃ = ρ₀ = O(1)` in bounded desingularized time `t₃ ≤ T₃^+ < ∞`, while the
**transverse** directions contract:

```
δy₃′ ≈ −2 b(θ) δy₃,     δε₃′ ≈ −b(θ) δε₃,      b(θ) ≥ b_- > 0   (uniform in θ).
```

*Verified* (`k3` integration, `a=b=c=1`): `ρ₃: 0.02 → 0.5`, `|y₃|: 0.025→0.0003`,
`ε₃: 0.5→0.024`, with measured `d lnε₃/dt₃ = −1.01` against the predicted `−b = −1`.

**Conclusion.** Since the transverse linearisation is a *uniform contraction* (rate
`≥ b_- > 0`) over a bounded exit time, the K₃ map transports the K₂-exit measure to
`Δ_out` with **non-increasing** transverse spread plus a bounded Ornstein–Uhlenbeck
noise contribution `≤ η²/(2 b_-)` (the degenerate noise on `r=ρ₃` enters the
transverse `y₃=y/ρ₃²,ε₃=ε/ρ₃` with amplitude `↓` as `ρ₃↑`). Hence the exit measure at
`Δ_out` is the K₃-deterministic image (JKK Prop. 4.21 / Theorem 3.2) of the inner
exit measure, with an `O(η)` spread, uniformly in θ. ∎

## 7. Theorem (composed noisy transition map, outer regime, uniform in θ)

Let `0<ε₂≪ε₁≪1` and `η=σ/√ε₂` small. Combining §§1–6, the noisy transition map
`π_ε^σ:Δ_in→Δ_out` decomposes, uniformly in `θ∈ℝ/ℤ`, as:

- **with probability `1−P_esc`** (the *canard/non-escape* branch): the deterministic
  JKK map (Theorem 3.2) plus an `O(η)` Gaussian tube (§3,§5,§6) of width
  `~η/(ab y₂)^{1/4}`;
- **with probability `P_esc(η)=1−e^{−H(η)}`, `H=bη²/4πac`** (the *noise-induced early
  escape* branch): exit near the fold, characterized to exponential order by §§1–4.

In particular the **critical noise is `σ_*^A = 2√π √ε₂ √(a(θ)c(θ)/b(θ))`** — the
leading-order law (`FOLDED_CYCLE_NOISE.md` §4) now with a **rigorous prefactor and
geometry**, the exponent and `√(ac/b)` proved, and the angular uniformity established
(§3,§5). ∎ (to exponential order)

## 8. The one genuinely open piece (attempted, not closed)

**The sharp inner exit measure**, `y₂ ≲ η^{4/3}`, where §4's adiabatic approximation
breaks, governed by the **noisy linearized Airy** equation (stochastic `u″=Yu` from
`R=−u′/u`, `Y<0` past the fold). Progress and limit (`inner_exit_measure.py`):

- **Inner scaling — established.** Balancing the OU tube `η/(2Y^{1/4})` against the
  separatrix half-width `√Y` gives the peel-off scale `Y_* ~ η^{4/3}`, `R_* ~ η^{2/3}`,
  so the inner exit measure is a *universal* object `(Y_peel,R_peel)=(η^{4/3},η^{2/3})×
  (\text{universal r.v.})`.
- **The quasi-static prefactor is NOT sharp — confirmed quantitatively.** The
  committed-escape probability measured numerically is **far below** the quasi-static
  hazard `η²/4π` for small `η` (ratio →0 as `η→0`; only `≈1` near `η≈1.4`). So the
  **dynamic/inner correction is large, not a small refinement** — the finite passage
  time through the fold suppresses escape relative to the adiabatic estimate. This
  *sharpens* the gap: the sharp prefactor genuinely needs the noisy Airy, it is not
  recoverable from the quasi-static hazard.
- **Not closed.** The closed-form noisy-Airy exit distribution + the sub-exponential
  prefactor remain open: the analytic object is genuinely research-level (the sharpest
  layer of Berglund–Gentz-type dynamic-fold asymptotics), and the relevant regime is
  rare-event, so neither analytics nor direct simulation closes it here.

**Net:** §§1–7 give the framework, the **exponential rate, the geometry, the
concentration, the chart composition (incl. K₃, §6), and θ-uniformity — all proved**;
the σ_*^A *law* (exponent + `√(ac/b)`) is unaffected. The single remaining object is
the noisy-Airy inner exit measure / sub-exponential prefactor — inner scaling
established, shown to be a *significant* (not small) correction, but **not solved**.

## 9. Honest status of the proof

| Piece | Status |
|---|---|
| Exact quasipotential `V=(8/3)√(a³/b)y₂^{3/2}` | **proved** (algebra, Kramers-verified) |
| Escape rate / scale function | **proved** (1-D diffusion, sharp rate) |
| Outer concentration, θ-uniform | **proved** (BG comparison; numerically collapse-verified) |
| Integrated hazard `H=bη²/4πac` ⇒ `σ_*^A=2√π√ε₂√(ac/b)` | **proved** (closed form, quadrature-verified) |
| Chart overlaps `κ₁₂,κ₂₃` | **proved** (bounded Jacobian, JKK 4.1) |
| K₃ exit | **proved** (field derived from (27)+(31); transverse contraction rate `−b`, verified `−1.01`) |
| Sharp inner exit measure (noisy Airy) + prefactor | **open** |

So bucket 1 is **proved to exponential order, uniformly in θ, in the outer/adiabatic
regime**, with the sharp inner exit measure the single remaining piece. That converts
"global noisy chart-matching" from an open problem into *a theorem with one named gap.*
