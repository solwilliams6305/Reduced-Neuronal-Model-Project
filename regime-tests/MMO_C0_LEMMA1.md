# Lemma 1 attempt — the O(μ) bound, and what it actually shows

*Targets gap (1) of `MMO_C0_PROOF.md` §6: upgrade "conditional on BGK's averaging
+ leading-order twist" to "unconditional, fixed k." Outcome: **gap (1) does not
close on its own.** The non-adiabatic remainder is O(μ) as hoped — but at fixed k
that is not enough to pin `c₀`, and the honest consequence is that gaps (1) and (2)
are the same problem. What survives unconditionally is an exact area-contraction
identity. Verified numerically (`/tmp/lemma1.py`, `/tmp/lemma1b.py`).*

---

## 0. Bottom line (so the rest is read in the right light)

- The remainder **is** O(μ): in the canonical frame the singular-value spread is
  `σ₁/σ₂ = 1 + ½μ + O(μ²)` (measured). So the easy estimate I promised does hold.
- **But it does not pin `c₀` at fixed k.** For fixed k the exponent `z₀²/μ` is
  *itself* O(μ), so an O(μ) directional ambiguity becomes an **O(1)** ambiguity in
  the extracted `c₀` (measured spread ≈ 0.48 at k=0, around a value ≈ 0.6).
- The one **unconditional, exact** statement is the **area (symplectic) contraction
  identity**: `−½ log|det M| / ((2k+1)²μ) = z₀²/((2k+1)²μ²) = c₀(φ₀)`, for all k, μ,
  by Liouville — no averaging, no WKB.
- The directional ambiguity vanishes only as **k grows** (spread 0.48→0.02→0.002 for
  k=0,2,8) — i.e. `c₀` is robustly defined exactly in the large-k regime, which is
  the O(1)-window, `ω→0` **turning-point** regime of gap (2). **Gaps (1) and (2)
  merge.** My earlier "(1) routine, (2) hard" split was wrong.

---

## 1. The structural fact that helps — and the one that bites

Variational equation `μu' = A(z)u`, eigenvalues `λ± = 2z ± 2iω(z)`, `ω=√(1+μ−z²)`.

- **Helps:** `λ±` are complex conjugates ⇒ *common* real part `2z` ⇒ a common
  contraction `e^{(1/μ)∫2z}=e^{−z₀²/μ}` and a rotation; there is **no growing mode**
  to amplify the non-adiabatic coupling.
- **Bites:** for fixed k the twist condition forces `|z₀| = O(μ)` (the window is
  O(μ) wide), hence `z₀²/μ = c₀(2k+1)²μ = O(μ)`: the contraction is only `1−O(μ)`.

## 2. Canonical frame — the leading flow is an isotropic contraction × rotation

With `S̃(z)=[[1,0],[−z,ω]]` (real parts / imag parts of the eigenvector):

```
 S̃⁻¹ A S̃ = [[ 2z,  2ω],
             [−2ω,  2z]] ,        S̃⁻¹S̃' = [[0,0],[−1/ω, −z/ω²]] = O(1).
```

So in `w=S̃⁻¹u` the leading generator `(1/μ)[[2z,2ω],[−2ω,2z]]` produces exactly
`e^{−z₀²/μ}·R(Φ)` — an **isotropic** contraction times a rotation. *All* singular-
value spread comes from the O(1) non-adiabatic term `−S̃⁻¹S̃'`.

## 3. The remainder is O(μ) (the part that worked)

Interaction picture `w=Φ₀ w̃`, `w̃' = −B w̃`, `B = Φ₀⁻¹(S̃⁻¹S̃')Φ₀`. The off-diagonal
entries of `B` carry the factor `exp((1/μ)∫(λ_j−λ_i)) = exp((1/μ)∫(∓4iω))`, of
**modulus 1** (conjugate eigenvalues) — no exponential amplification. Hence
`‖B‖=O(1)` uniformly, and over a window `|z₀|=O(μ)`, Gronwall gives transition
matrix `= I + O(μ)`. Measured (canonical frame, k=0):

```
   mu     sig1/sig2 - 1    (ratio-1)/mu
  0.100      4.66e-2          0.466
  0.020      9.85e-3          0.493
  0.005      2.49e-3          0.498     ->  sig1/sig2 = 1 + (1/2) mu + O(mu^2)
```

and `det M = e^{−2z₀²/μ}` reproduced to integrator accuracy (Liouville).

## 4. Why that is not enough at fixed k (the bite, quantified)

`σ₁σ₂ = |det M| = e^{−2z₀²/μ}` is exact, but `σ₁,σ₂ = e^{−z₀²/μ}(1 ± ¼μ)`.
Extracting `c₀ = −log(dist)/((2k+1)²μ)` from a *directional* distance gives

```
 c₀ = c₀(φ₀)  ∓  (¼μ)/((2k+1)²μ)  =  c₀(φ₀) ∓ 1/(4(2k+1)²) ,
```

an **O(1)** shift at k=0. Measured (exact variational eq, three extractions):

```
  k   mu    c0(phi0)  c0_geo   c0_lo   c0_hi   hi-lo
  0  0.020   0.6048   0.6078   0.363   0.853   0.490     <- O(1) ambiguous
  2  0.020   0.6060   0.6089   0.599   0.619   0.020
  8  0.020   0.6193   0.6225   0.622   0.623   0.002     <- robust
```

`c0_geo` (area/geometric mean) tracks `c₀(φ₀)`; the `lo/hi` spread is the honest
uncertainty. At k=0 it is ≈0.49 — `c₀` is simply not pinned by the fixed-k window.

## 5. The clean invariant that survives, and the merge of the two gaps

**Exact, unconditional (Liouville):**

```
 c₀_area := −½ log|det M(0,z₀)| / ((2k+1)²μ) = z₀²/((2k+1)²μ²) = c₀(φ₀),   ∀ k, μ.
```

No averaging, no WKB, no fixed-k restriction. So **if** the secondary-canard
spacing is the area/symplectic contraction (which is what BGK's averaged radial
coordinate is), then `c₀ = c₀(φ₀)` *exactly*. The whole question reduces to that
identification — and §4 shows it holds to relative O(μ) only when the directional
spread `1/(4(2k+1)²)` is negligible against `c₀(φ₀)`, i.e. for **large k**. Large k
is precisely the O(1)-window, `ω(z₀)→0` turning-point regime of gap (2).

**Therefore gaps (1) and (2) are one problem,** not two. The fixed-k regime is
analytically easy but informationally too weak to carry `c₀`; the large-k regime
carries `c₀` robustly but needs the turning-point analysis.

## 6. Honest status and revised plan

- **Proved unconditionally:** (i) the area-contraction identity
  `c₀_area = c₀(φ₀)` (all k, exact); (ii) the non-adiabatic remainder is O(μ) over
  the fixed-k window (`σ₁/σ₂ = 1+½μ`).
- **Not closed:** pinning the *secondary-canard spacing* to `c₀(φ₀)` for fixed k —
  the directional ambiguity is O(1) there.
- **Correction to `MMO_C0_PROOF.md` §6/§7:** the "two lemmas, one routine" framing
  was too optimistic. There is effectively **one** substantive theorem left:
  show the canard spacing equals the area-contraction to leading order, which is
  automatic at large k but requires the `ω→0` turning-point (Airy/PCF) analysis.
- **Best unconditional dissertation statement available now:** *the area
  (symplectic) contraction of the weak-canard variational flow between entry and
  fold section is exactly `e^{−z₀²/μ}`, so the area-measured secondary-canard
  spacing constant is exactly the closed form `c₀(φ₀)` (→ π²/16 fixed-k, → 1 deep)*
  — with the identification of the geometric spacing with the area contraction, in
  the clustering (large-k) regime, as the one remaining theorem.

This is a real sharpening (the area identity is exact and unconditional) and an
honest demotion (fixed-k alone cannot pin the constant). It also says where the
single remaining difficulty lives: the turning point, not a Gronwall estimate.
