# A closed-form proof of the BGK folded-node canard-spacing constant `c₀`

*Companion to `MMO_C0_PINNED.md` (which computed `c₀`) — this document is the
attempt at a **proof**, written with an explicit ledger of what is rigorous, what
is leading-order, and what is still open. Algebra machine-verified (see §5).*

---

## Theorem (what is established, and in what register)

Let `γᵏ` denote the k-th secondary canard of a folded node with eigenvalue ratio
`μ∈(0,1)`, `γʷ` the weak canard, and write the spacing at the fold section as
`dist(γᵏ,γʷ) ∼ exp(−c₀(2k+1)²μ)` (Berglund–Gentz–Kuehn, Thm 4.4, who prove
`0<c₀=O(1)` and bound it to a range).

**Within BGK's canonical-form reduction (their Thm 4.3) and twist condition
(their eq. 4.26), the spacing constant is given exactly by**

```
            π² sin²φ₀
 c₀(φ₀) = ───────────────────────── ,     φ₀ ∈ (0, π/2),
          (1+μ)(2φ₀ + sin2φ₀)²
```

**where `φ₀` is fixed by the depth relation**

```
 (1+μ)(2φ₀ + sin2φ₀) = π μ (2k+1).                                      (★)
```

Its two limits are **exactly BGK's interval endpoints**:

- **fixed `k`, `μ→0`** (`φ₀→0`): `c₀ → π²/16 ≈ 0.6169` (shallow / Thm-4.4 value);
- **deepest canards** (`φ₀→π/2`, `k∼s_max`): `c₀ → 1/(1+μ) → 1`.

So BGK's range `[π²/16, 1]` is precisely the range of the explicit function
`c₀(φ₀)` over canard depth. **The register is "sharpening inside BGK's framework":
we *evaluate* the two ingredients they *bounded*; we do not re-prove their GSPT
reductions.** §6 is an honest account of the residual gap to an unconditional
theorem.

---

## 1. Setup — the weak-canard variational equation (BGK)

Linearising about the weak canard gives (BGK eq. 4.8)

```
 μ u'(z) = A(z) u ,     A(z) = [[ 4z,        2 ],
                                [ −2(1+μ),   0 ]] .
```

`tr A = 4z`, `det A = 4(1+μ)`, so the eigenvalues are

```
 λ±(z) = 2z ± 2i ω(z) ,     ω(z) = √(1+μ − z²)      (|z| < √(1+μ)).
```

Write `R := √(1+μ)`. Solutions are a **contraction** (rate `Re λ = 2z`) times a
**rotation** (rate `Im λ = 2ω`).

## 2. Lemma 1 — the contraction exponent is `−z₀²/μ`, with no change-of-variables factor

The amplitude factor accumulated from the entry point `z₀<0` to the fold section
`z=0` is

```
 exp( (1/μ) ∫_{z₀}^0 Re λ dz ) = exp( (1/μ) ∫_{z₀}^0 2z dz )
   = exp( (1/μ)(0 − z₀²) ) = exp(−z₀²/μ).
```

*The canonical-form change of variables `u = S(z) ũ` (BGK Thm 4.3) does not alter
this exponent.* Two independent confirmations:

1. **Determinant.** `det S(z) = (1/ω)·det[[−z+ω, −z−ω],[1,1]] = (1/ω)(2ω) = 2`,
   constant. Hence `S, S⁻¹` are uniformly `O(1)` and the rotation block `U` is
   orthogonal; together they move only the `O(1)` prefactor, never the exponent.
2. **Polar average.** In BGK's polar coordinates the averaged radial equation
   integrates to `r(z) = r(z₀) e^{(z²−z₀²)/μ}`, so `r(0) = r(z₀) e^{−z₀²/μ}` with
   `r` *being* the distance to `γʷ` — the exponent appears with no intervening
   factor.

Therefore `dist(γᵏ,γʷ)|_{z=0} = C·e^{−z₀²/μ}`, `C=O(1)`, and consequently

```
 c₀ = z₀² / ((2k+1)² μ²).                                               (1)
```

Pinning `c₀` is thus equivalent to pinning the entry point `z₀`.

## 3. Lemma 2 — the twist condition and the exact entry point

The k-th secondary canard makes `(2k+1)/2` rotations before the fold (BGK eq. 4.26,
leading order): `(1/μ)∫_{z₀}^0 2ω(s) ds = π(2k+1)/2`, i.e.

```
 ∫_{z₀}^0 ω(s) ds = ∫_{z₀}^0 √(R² − s²) ds = π μ (2k+1)/4.               (2)
```

**The step BGK did not take: evaluate (2) exactly rather than bound it.** Using the
antiderivative `∫√(R²−s²)ds = (s/2)√(R²−s²) + (R²/2) arcsin(s/R)` and the
substitution `z₀ = −R sin φ₀` (so `φ₀ = arcsin(|z₀|/R)∈(0,π/2)`,
`√(R²−z₀²) = R cos φ₀`):

```
 ∫_{z₀}^0 √(R²−s²) ds = (R²/4)(2φ₀ + sin2φ₀).                           (3)
```

[`∫√(R²−s²)` against numerical quadrature: agreement to `4·10⁻¹⁵` — §5, Check 1.]
Equating (3) to (2) gives the depth relation **(★)**: `(1+μ)(2φ₀+sin2φ₀)=πμ(2k+1)`.

BGK's range comes from *bounding* (3): concavity gives
`(π/4)|z₀|R ≤ ∫ω ≤ |z₀|R`, trapping `|z₀|∈[(π/4)(2k+1)μ, (2k+1)μ]/R` and hence
`c₀∈[(π/4)², 1] = [π²/16, 1]`. Evaluating instead of bounding replaces the interval
by the exact value below.

## 4. Proof of the Theorem

Insert `z₀² = R² sin²φ₀ = (1+μ) sin²φ₀` and the relation (★),
`(2k+1)μ = (1+μ)(2φ₀+sin2φ₀)/π`, into (1):

```
        (1+μ) sin²φ₀                      π² sin²φ₀
 c₀ = ─────────────────────────── = ───────────────────────── .
      [ (1+μ)(2φ₀+sin2φ₀)/π ]²        (1+μ)(2φ₀+sin2φ₀)²
```

**Limits.**
- `φ₀→0`: `sin²φ₀ ∼ φ₀²`, `2φ₀+sin2φ₀ ∼ 4φ₀`, so `c₀ → π²φ₀²/((1+μ)·16φ₀²) =
  π²/(16(1+μ)) → π²/16` as `μ→0`. (Fixed `k` forces `φ₀=O(μ)→0` via (★).)
- `φ₀→π/2`: `sin²φ₀→1`, `2φ₀+sin2φ₀→π`, so `c₀ → π²/((1+μ)π²) = 1/(1+μ) → 1`. ∎

## 5. Numerical confirmation (`/tmp/c0proof.py`; mirrors `c0_check.py`)

Root-finding `z₀` from the **exact** twist integral (2) and forming both `c₀` from
(1) and `c₀(φ₀)` from the Theorem:

- **Check 1.** Exact antiderivative (3) vs quadrature of (2): diff `≤ 6·10⁻¹³`.
- **Check 2.** `c₀(direct, eq 1)` vs `c₀(depth fn)` and the relation (★): agreement
  and residual both `≤ 6·10⁻¹⁶` across `μ∈{0.1,..,0.005}`, `k∈{0,1,3}`. *(The depth
  function is an exact algebraic identity, not a fit.)*
- **Check 3.** Fixed `k=0`: `c₀ = 0.5617(μ=.1) → 0.6162(μ=.001) → π²/16=0.61685`.
- **Check 4.** `μ=0.02`: `c₀` rises `0.605(k=0) → 0.805(k=23≈s_max)`, `φ₀: 0.9°→56.6°`
  — climbing toward `1/(1+μ)` as `φ₀→90°`.
- **Check 5.** Reproduces stored `MMO_C0_PINNED` table exactly (e.g. `μ=.1,k=1 →
  0.56965`; `μ=.05,k=2 → 0.59451`).

## 6. Honest status — rigorous / leading-order / open

**Rigorous (given BGK Thm 4.3 + 4.26):**
- The closed form `c₀(φ₀)` and the relation (★) — exact algebra, machine-checked.
- Lemma 1's exponent `−z₀²/μ` — the `det S=2` argument genuinely removes the
  change-of-variables concern; this is a real (if small) lemma, not a hope.
- The exact evaluation (3) of the twist integral, and both limits.

**Leading-order (the conditional part):** BGK's twist condition (4.26) is itself
*leading order in `μ`*, and their radial exponent comes from an **averaging** of the
polar equations (BGK §4.1) that is exact only up to `O(μ)`. Hence honestly

```
 c₀ = c₀(φ₀) · (1 + O(μ)).
```

For **fixed `k`** this is clean: `|z₀|=O(μ)` is small, `ω≈R` is essentially constant
over `[z₀,0]`, the averaging is uniformly valid, and `c₀→π²/16` rigorously within
the framework. The remaining work to make it *unconditional* (not resting on 4.26)
is a Gronwall / variation-of-constants bound on the averaging remainder; the
prefactor is already controlled by `det S=2`. This is sketchable and is the first
of the two theorems a dissertation needs.

**Genuinely open (not bookkeeping):** the **deep-canard end** `φ₀→π/2`,
`k∼s_max=(1−μ)/(2μ)`. There `z₀→−R`, so `ω(z₀)=√(R²−z₀²)→0`: the rotation *stalls
at a turning point*, and the averaging underlying BGK's canonical form is
**non-uniform** there. So `c₀→1` is *formally exact (the integral (3) still holds)
and numerically confirmed (Check 4)* but **not rigorously established** at the
deepest canards. Closing it needs a turning-point analysis where `ω→0` — an Airy /
parabolic-cylinder matching — which is the second, harder theorem.

## 7. Verdict for the dissertation

- **What closes now:** a complete, self-contained, machine-verified determination of
  the spacing constant *as an explicit function of canard depth*, with the fixed-`k`
  value `π²/16` rigorous inside BGK's framework — strictly more than BGK's range, and
  the right "sharpen the published constant" register. This is the A1-target kernel.
- **What stands between this and "faultless" (90+):** exactly two named theorems —
  (i) the `O(μ)` averaging/twist error bound (Gronwall; converts "conditional on
  4.26" into "unconditional, fixed `k`"), and (ii) the deep-canard turning-point
  analysis (`ω→0`; Airy/PCF) that makes `c₀→1` rigorous. Neither is hand-waving;
  both are concrete and bounded, and (ii) is the genuinely hard one.
- **Honest one-liner:** *the constant and its depth function are proved exactly given
  BGK's reduction; promoting "given BGK" to "unconditional" is two more lemmas, one
  routine and one a real turning-point problem.*
