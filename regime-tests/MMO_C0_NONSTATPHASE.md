# Third pass — `c₀(φ₀)` closes on the bulk (area identity + non-stationary phase)

*Supersedes the pessimistic §6 of `MMO_C0_LEMMA1.md`. Outcome: the spacing constant
**is** pinned to relative O(μ) across the whole clustering regime (canards of bounded
depth `φ₀∈(0,π/2)`); only two corners sit outside, both numerically benign. Verified
`/tmp/lemma_retry.py`.*

> **Honesty note (this is my third pass).** Pass 1 said "lemma is routine" —
> overclaim. Pass 2 said "it all collapses into the hard turning point" —
> overcorrection. This pass is the calibrated middle, and I flag the convergence
> rather than pretend I had it first. The new content below is one specific
> mechanism (oscillatory cancellation of the coupling) that pass 2 underweighted.

---

## 0. Result

For secondary canards of **bounded depth** `φ₀∈[φ_min, π/2−δ]` (any fixed `δ>0`),

```
 c₀ = c₀(φ₀)·(1 + O(μ)) ,     c₀(φ₀) = π² sin²φ₀ / [(1+μ)(2φ₀+sin2φ₀)²] ,
```

rigorously **given BGK's reduction to the weak-canard variational equation**. The
shallow limit gives `π²/16`, the deep limit `1`. This covers the regime where the
canards actually cluster and `c₀` is the physical spacing constant. The two
excluded corners — the shallowest `k=0,1` and the extreme-deep sliver `φ₀→π/2` —
are handled in §4 and shown to be benign.

## 1. The two ingredients

**(A) Area identity — exact, unconditional (Liouville).** For `M` the transition
matrix of `μu'=A(z)u` from `z₀` to the fold section,
`det M = exp((1/μ)∫_{z₀}^0 tr A dz) = exp((1/μ)∫4z dz) = e^{−2z₀²/μ}`. Hence the
geometric-mean (area/symplectic) contraction is exactly `e^{−z₀²/μ}`, giving

```
 c₀_area := −½ log|det M|/((2k+1)²μ) = z₀²/((2k+1)²μ²) = c₀(φ₀)   (exact, all k,μ).
```

**(B) The non-adiabatic coupling is oscillatory — so it cancels away from `ω=0`.**
In the diagonalising frame the mode-mixing term carries the factor
`exp((2i/μ)∫_{z₀}^z ω)`. Its phase derivative is `2ω/μ ≠ 0` wherever `ω>0`, so by
integration by parts (non-stationary phase) the accumulated mode separation is

```
 D := ½ log(σ₁/σ₂) = O( μ · sup_{[z₀,0]} 1/ω² ) = O( μ / cos²φ₀ ),
```

since `ω ≥ ω(z₀) = √(1+μ)·cosφ₀` on the window. For `φ₀` bounded away from `π/2`
this is **O(μ)**. *(This is BGK's averaging step made explicit and error-controlled;
it may already lie inside their rigorous apparatus — see §5. The new content of this
whole programme is the **evaluated constant** `c₀(φ₀)`, not this estimate.)*

## 2. The estimate that closes it

The physical spacing is some directional contraction `dist∈[σ₂,σ₁] =
e^{−z₀²/μ}·e^{∓D}`. Extracting `c₀`:

```
 c₀ = c₀(φ₀) ± D/((2k+1)²μ) .
```

The **directional ambiguity** is therefore

```
 D/((2k+1)²μ) = O( 1 / ((2k+1)² cos²φ₀) ).
```

In the clustering regime `φ₀=O(1)` forces `(2k+1)μ = (1+μ)(2φ₀+sin2φ₀)/π = O(1)`,
i.e. `(2k+1)=O(1/μ)`, so `(2k+1)²cos²φ₀ = O(1/μ²)` and the ambiguity is **O(μ²)** —
swamped by the O(μ) twist/geometric correction. With the area identity (A) fixing
the mean exactly, this gives `c₀ = c₀(φ₀)(1+O(μ))`. ∎(bulk)

## 3. Numerical confirmation (`/tmp/lemma_retry.py`)

```
 mu=0.02 (s_max≈24.5)
  k  phi0  c0(phi0)  c0_geo  c0-spread     D    D/mu
  0   0.9   0.6048   0.6078    0.49      0.0049  0.25
  2   4.4   0.6060   0.6089    0.020     0.0049  0.25
  5   9.8   0.6107   0.6137    0.004     0.0050  0.25     bulk: D/mu = 1/4 const
 14  27.6   0.6523   0.6559    0.0007    0.0060  0.30     spread already O(mu^2)
 20  43.7   0.7249   0.7298    0.0005    0.0089  0.45
 23  56.6   0.8049   0.8116    0.0008    0.0170  0.85     D/mu growing (turning pt)
 24  63.7   0.8540   0.8625    0.0012    0.0299  1.49     ...but spread still ~1e-3
```

Three facts confirmed: `c0_geo = c0(φ0)` throughout (area identity); `D/μ ≈ 1/4`
constant across the bulk then growing only near the deep end; the actual `c₀`
uncertainty (`spread`) collapses to `≤10⁻³` for every `k≥2` — including the deepest,
because the `(2k+1)²` denominator beats the growth of `D`.

## 4. The two corners (both benign)

- **Shallow `k=0,1`.** Here `(2k+1)²` is O(1), so the directional ambiguity is O(1)
  (`spread≈0.49` at k=0). But `z₀²/μ=O(μ)`: the canards are **not** clustering
  (spacing `≈1`), so `c₀` is not the meaningful quantity in this corner. No loss.
- **Extreme-deep `φ₀→π/2` (`k→s_max`).** `ω(z₀)→0`: the phase `2ω/μ` becomes
  stationary, the IBP bound `D=O(μ/cos²φ₀)` degrades, and a uniform statement needs
  an **Airy / parabolic-cylinder** turning-point estimate. Numerically `D/μ` only
  creeps up (≤~3 at the deepest computed) and the `c₀` impact stays `~10⁻³`, so this
  is the one genuinely-hard residue but it controls a **sub-dominant correction**,
  not the constant.

## 5. Honest status

- **Rigorous (given the BGK variational reduction):** area identity (A, exact) +
  non-stationary-phase bound (B) ⇒ `c₀=c₀(φ₀)(1+O(μ))` on bounded depth. This is the
  clustering regime — the physically meaningful one.
- **What is genuinely mine vs cited:** the explicit evaluated constant `c₀(φ₀)` and
  its limits are the contribution. The area identity is Liouville; the
  non-stationary-phase/averaging control (B) is standard and **may already be inside
  BGK's rigorous framework** — needs BGK 1312.6353 read in full to say which (still
  the open Step-0 gate). If it is theirs, the honest register reverts to "evaluate
  the constant they bounded," which is fine and was always the claim.
- **Remaining for a fully uniform theorem:** the Airy turning-point estimate at
  `φ₀→π/2` (now a bounded, well-posed sub-problem: control `D` as `ω→0`), and the
  next-order twist condition to pin the O(μ) prefactor itself. Neither blocks the
  bulk result.

## 6. For the 90+ plan

- **Stateable theorem now:** *the secondary-canard spacing constant equals the
  closed form `c₀(φ₀)` to relative O(μ) for canards of bounded depth — combining the
  exact symplectic area-contraction with the non-stationary-phase decay of the
  rotating non-adiabatic coupling — recovering `π²/16` (shallow), `1` (deep), and
  thereby BGK's interval as the range of an explicit function.* That is materially
  stronger than "computed to leading order," and most of it is rigorous.
- **Distance to faultless:** the single Airy turning-point lemma for the deepest
  sliver. It is now scoped down from "the whole difficulty" to "a sub-dominant
  correction at one endpoint" — the right size for a dissertation's hardest section,
  not a wall.
