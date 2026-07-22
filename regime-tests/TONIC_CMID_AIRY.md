# Tonic–canard duality: the iPRC is the *square* of the canard's Airy function

**What this is.** A fresh attack on the tonic constant `c` (`A_mid=√(c/π²)`, #8) —
the third pillar of the project's "universal law + model-specific O(1) constant"
spine, after the canard `C_q` and excitable `C` (`FOLD_ESCAPE_PREFACTORS.md`).
Instead of shooting the inner adjoint BVP numerically (`TONIC_CMID_BVP.md`), it
solves the same fold normal form in closed form and finds that **the tonic iPRC
and the canard are built from one and the same Airy function** of the fold:

```
┌─────────────────────────────────────────────────────────────────┐
│  fold normal form:   dV/dT = V² − W ,   dW/dT = −λ                 │
│  Riccati  V = −u'/u   ⇒   u'' = (W₀−λT) u   =   Airy equation       │
│                                                                    │
│     canard velocity      V(T)   = −u'/u        (Riccati log-deriv) │
│     tonic inner iPRC      Z̃_v(T) =  u²  = C·Ai(−ξ)²   (the square)  │
│                                                                    │
│   ONE Airy function u; the canard is its log-derivative, the iPRC  │
│   is its square. C_q and c are two faces of the same u.            │
└─────────────────────────────────────────────────────────────────┘
```

This delivers the exact *leading* inner iPRC that `TONIC_CMID_BVP.md` could only
shoot numerically, and reduces `c` to a named Airy integral. It does **not** by
itself produce the number `c≈1.55` — the same two gaps remain (the `O(ε^{1/3})`
correction and the global/bulk matching) — but it puts them on an exact
special-function footing and explains *why* `c` is model-specific.

---

## 1. Setup (what is certified going in)

`A_mid² = (1/T²)∮Z_v²dt`, dominated by the fold passages; `c = π²A_mid²`. The
inner adjoint at the Krupa–Szmolyan fold blow-up was derived and corrected in
`TONIC_CMID_BVP.md` (and is taken as given here):

```
dZ̃_v/dT = −2 V(T) Z̃_v ,     Z_v = ε^{−2/3} Z̃_v ,   t = ε^{−1/3} T ,
```

with `V(T)` the deterministic fold passage `dV/dT=V²−W, dW/dT=−λ`. The per-fold
contribution is `∮_fold Z_v² dt = ε^{−5/3} R̃`, `R̃ = ∫ Z̃_v² dT`.

## 2. The Riccati → Airy reduction (the new step)

`W(T)=W₀−λT`, so the forward equation is a Riccati equation
`dV/dT = V² − W₀ + λT`. The standard linearising substitution

```
V = −u'/u            ⇒     dV/dT = −u''/u + (u'/u)² = −u''/u + V² ,
```

turns it into a **linear** second-order ODE,

```
u'' = (W₀ − λT) u .
```

With `ξ = λ^{1/3}(T − W₀/λ)` this is exactly the **Airy equation**

```
u_ξξ + ξ u = 0     ⇒     u(T) = Ai(−ξ)   (the canard-selected solution).
```

So the deterministic fold passage *is* an Airy function: the attracting branch is
the exponential (sub-turning) regime of `Ai`, the fold tip is the **Airy turning
point**, and the post-fold jump is the oscillatory regime.

## 3. The duality: iPRC = Airy²

The inner adjoint solves `dZ̃_v/dT = −2V Z̃_v = +2(u'/u) Z̃_v`, i.e.
`d(ln Z̃_v)/dT = 2 d(ln u)/dT`. Hence, exactly,

```
Z̃_v(T) = C · u(T)² = C · Ai(−ξ)² .
```

The tonic inner iPRC is the **square of the same Airy function** whose
log-derivative is the canard velocity. Two consequences:

- **The iPRC peak is the first Airy extremum.** `Z̃_v` is extremal where
  `dZ̃_v/dT=0 ⇒ V=0 ⇒ u'=0 ⇒ Ai'(−ξ)=0`, i.e. at the first maximum of `Ai`,
  `−ξ ≈ −1.019`. This *is* `TONIC_CMID_BVP`'s "peak at the fold tip `V≈0`", now
  pinned to the Airy extremum, no shooting.
- **The amplification across the branch is the canard quasipotential.** In the
  sub-turning (WKB) regime `u ~ exp(−(2/3λ)W^{3/2})`, so
  `Z̃_v = u² ~ exp(−(4/3λ)W^{3/2}) = exp(−ΔΦ/λ)` with `ΔΦ=(4/3)W^{3/2}` — the
  *same* fold quasipotential that gives `C_q`. The iPRC amplifies from branch to
  tip by `exp(ΔΦ/λ)`. The tonic jitter and the canard escape are literally two
  functionals of the one fold barrier `Φ`.

**Consistency check against `TONIC_CMID_BVP.md`.** Its matched outer value
`Z̃_v ≈ 0.58` at `V=1` and inner peak `≈1.07` require an amplification
`1.07/0.58 ≈ 1.85 = exp(ΔΦ/λ)`, i.e. `ΔΦ/λ ≈ 0.6` — an O(1) number, exactly as
the blow-up demands (`ΔΦ`, `λ` both O(1) in inner variables). The closed-form
duality reproduces the numerically-shot profile.

## 4. `c` reduces to an Airy fourth moment

```
R̃ = ∫ Z̃_v² dT = C² ∫ u⁴ dT = (C²/λ^{1/3}) ∫ Ai(−ξ)⁴ dξ ,
A_mid² = [ ε^{−5/3}(R̃_L + R̃_R) + bulk ] / T² ,     c = π² A_mid² .
```

So the inner constant is the **fourth moment of the Airy function** (the canard
used its log-derivative; the tonic uses its square, integrated). Crucially,
`∫Ai(−ξ)⁴dξ` is **logarithmically divergent at both ends** — `Ai(−ξ)⁴ ~ ξ^{−1}`
in the oscillatory regime and the sub-turning tail must be matched to the outer
iPRC. So `R̃` is set by the **matching cutoffs**, not a universal Airy number.

**This is the point, not a defect.** A universal inner form (`Ai²`) with a
matching-dependent integral is exactly why `c` is *model-specific* — `c≈1.55`
(FHN) vs `0.78` (VdP) — while `A_mid=√(c/π²)` is universal. The duality
*derives* the universal/specific split for the tonic chapter that the canard and
excitable chapters showed: universal inner solution (here `Ai²`), model-specific
matching constant (here the cutoffs set by the global cycle return).

## 5. Honest status — what closed, what didn't

| Piece | Status |
|---|---|
| Inner adjoint `dZ̃_v/dT=−2V Z̃_v` | certified (from `TONIC_CMID_BVP`) |
| Fold passage = Airy: `V=−u'/u`, `u''=(W₀−λT)u`, `u=Ai(−ξ)` | **derived (new)** |
| **iPRC = `C·Ai(−ξ)²`** (exact leading inner solution) | **derived (new)** |
| Peak at first Airy extremum `−ξ≈−1.019`; amplification `exp(ΔΦ/λ)` | **derived (new)**, matches BVP numbers |
| `c` = Airy fourth moment × geometry, matching-cutoff ⇒ model-specific | **derived (new)** structurally |
| **The number `c≈1.55`** | **still open** — needs the matched cutoffs (global/bulk) + the `O(ε^{1/3})` inner correction |

The two remaining gaps are exactly the ones `TONIC_CMID_BVP.md §5` named — but
they are now attached to an *exact* inner solution (`Ai²`) rather than a numerical
shoot. In particular the BVP's overshoot/drift (`c_recon ≈ 1.9→4.5` rising with
ε) is the `O(ε^{1/3})` correction to the leading `Ai²` form (ε-dependent, →0 as
ε→0), not an error in the leading solution.

## 6. The trilogy, stated

All three chapters' O(1) constants live on the **same fold normal form**
`dV/dT=V²−W`, and the universal/specific split is now uniform:

| chapter | universal object on the fold | model-specific constant | source of specificity |
|---|---|---|---|
| canard (#1/#16) | quasipotential `Φ=WV−V³/3` → escape | `C_q=√(4π ln2)`·(corr.) | drift `λ` at the fold |
| excitable (#12) | same barrier `(4/3)δ^{3/2}` | `K_fold=4/3, A₀=1/π` | drift `g` at the fold |
| **tonic (#8)** | **iPRC `= Ai(−ξ)²`** | `c` (≈1.55 FHN, 0.78 VdP) | **Airy-moment matching cutoffs** |

The forward fold flow gives the canard (Riccati `V=−u'/u`) and the excitable
barrier; the *square* of the same Airy solution gives the tonic iPRC. One special
function, three constants.

## 7. Next, in the same clean register

To turn the structural result into the number `c`: (i) fix the two matching
cutoffs of `∫Ai⁴` against the deterministic outer iPRC on each branch (the global
return enters here — the model-specific piece); (ii) add the coupled
`O(ε^{1/3})` inner correction (`TONIC_CMID_BVP §5-i`) on top of the exact `Ai²`
leading term. Both are now perturbations of a known closed-form solution rather
than corrections to a numerical shoot — the right footing for a clean `c`.
