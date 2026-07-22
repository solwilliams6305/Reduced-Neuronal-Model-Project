# Fold-escape noise prefactors — `C_q` certified, excitable `C` derived, the two unified

**What this is.** A certification pass on the canard `C_q` derivation
(`DERIVATION_STRATEGIES.md` §1) — the same scrutiny that caught the Route-B sign
bug — followed by harvesting the *identical* Freidlin–Wentzell fold-escape method
to close the excitable constant `C = K_fold/A₀` (#12), which the audit flagged
"stated computable, never computed." The payoff is a single formula for both
chapters' noise thresholds:

```
┌──────────────────────────────────────────────────────────────┐
│   σ_crit = √( 4π · D · ln2 ) · √ε ,    D = slow drift AT the fold │
│   D = λ = b(I−I_fold)        (canard, #1/#16)                    │
│   D = g  = |a−1−b(I−2/3)|    (excitable, #12)                    │
│   the √(4π ln2) is universal; D is the regime-specific constant. │
└──────────────────────────────────────────────────────────────┘
```

This is the universal-law + model-specific-constant split — now *derived* for the
noise threshold itself, not just asserted.

---

## Part 1 — Certifying `C_q = √(4π ln2)`

Independent re-derivation of every step. The blow-up canard SDE (Krupa–Szmolyan):

```
dV/dT = V² − W + η dB_T ,   dW/dT = −λ ,   η = σ/√ε .
```

**Step 1 — the frozen-W barrier (re-derived).** At fixed `W>0`, `V²−W = −Φ'(V)`
with `Φ(V) = WV − V³/3`. Critical points `Φ'=W−V²=0 ⇒ V=±√W`; `Φ''=−2V`, so
`V₋=−√W` is the well (`Φ''=+2√W`) and `V₊=+√W` the saddle (`Φ''=−2√W`).

```
Φ(±√W) = ±W^{3/2} ∓ W^{3/2}/3 = ±(2/3)W^{3/2}
⇒  ΔΦ = Φ(V₊) − Φ(V₋) = (4/3) W^{3/2} .            ✓ CERTIFIED
```

**Step 2 — the Kramers rate (re-derived).** For `dV=−Φ'dT+η dB` the stationary
density is `∝ exp(−2Φ/η²)` (diffusion `D=η²/2`). Overdamped Kramers,
`k = √(Φ''_min |Φ''_max|)/(2π) · exp(−2ΔΦ/η²)`, with
`√(2√W·2√W)=2√W`:

```
k(W) = (√W / π) · exp( −(8/3) W^{3/2} / η² ) .       ✓ CERTIFIED
```

**Step 3 — the accumulated-hazard collapse (re-derived).** `dW/dT=−λ ⇒
dT=−dW/λ`; integrating the hazard over the passage `W: ∞→0`:

```
H = (1/λ) ∫₀^∞ (√W/π) e^{−(8/3)W^{3/2}/η²} dW .
```

Substitute `u=(8/3)W^{3/2}/η²` ⇒ `√W dW = (η²/4) du`; the `√W` Kramers prefactor
cancels the Jacobian *exactly*:

```
H = (1/λ)(1/π)(η²/4) ∫₀^∞ e^{−u} du = η² / (4π λ) .   ✓ CERTIFIED
```

**Step 4 — the constant.** Escape is a non-homogeneous Poisson process, survival
`e^{−H}`. The measured criterion (median `R_hit` crosses the window) is 50 %
escape, `e^{−H}=½ ⇒ H=ln2`:

```
η²/(4πλ)=ln2 ⇒ η_* = √(4π ln2)·√λ ⇒ C_q = √(4π ln2) = 2.9514… ✓
general quantile p:  C_q(p) = √( −4π ln(1−p) ).
```

### Certification verdict

No algebra bug (unlike Route B). The structure is **rigorous**: the barrier
`(4/3)W^{3/2}` is exact, and the hazard collapses to `η²/(4πλ)` exactly given the
Kramers-form integrand. But the *number* `2.95` is a **leading-order** value, not
3-figure-rigorous, for three named reasons — stated honestly:

1. **Kramers prefactor is large-barrier asymptotic.** The hazard integral is
   dominated by `W` where `(8/3)W^{3/2}/η² = O(1)` — i.e. the barrier is *not*
   ≫ η². There the Kramers prefactor `√W/π` carries an `O(1)` correction (the
   exact escape rate over the cubic barrier differs from Kramers at finite
   barrier). This is the main uncertainty in the prefactor.
2. **Adiabatic (frozen-W) approximation.** Valid because the barrier changes
   slowly vs the escape-attempt rate — the ratio is `O(λ)` near the window
   (`barrier-change/attempt-freq ~ λ^{7/6}/λ^{1/6} = λ → 0`). So the
   non-adiabatic correction is `O(λ)`, vanishing as `λ→0`, consistent with the
   ~5 % residual at the swept `λ ≲ 0.04`.
3. **Definitional.** The measured `2.8` (normal form) uses `V_cross=1`, *not* the
   threshold-free separatrix `V₊=√W` this derivation uses; part of `2.95 vs 2.8`
   is this criterion difference, not error.

**Status:** `C_q = √(4π ln2) ≈ 2.95` certified as a clean, parameter-free,
*leading-order* result; the structure is exact; the prefactor inherits an `O(1)`
correction from the finite-barrier regime (item 1) — the well-posed next step to
make it 3-figure-rigorous is the exact (non-Kramers) cubic-barrier escape rate.
The universal half of #16 (VdP) rides along unchanged (identical blow-up).

---

## Part 2 — Deriving the excitable `C = K_fold/A₀` (#12), same method

The audit: excitable `σ_crit = √(ε/C)`, `C = K_fold/A₀`, "computable from
geometry, no free parameters" — but `K_fold, A₀` never computed. They are exactly
what the Part-1 method delivers.

**Fold geometry (no free parameters).** FHN fast nullcline `f(v)=v−v³/3`,
`f'(−1)=0` (left fold), `f''(−1)=2`. With `x=v+1`, `Δw=w−w_f`, `w_f=I−2/3`:

```
ẋ = f(v) − w + I = x² − Δw + O(x³)      (½f''(−1)=1 ⇒ clean x² fold).
```

**`K_fold` — the barrier coefficient.** Identical to Part 1 with `W→δ=|Δw|`
(distance to fold): `U(x)=δx−x³/3`, barrier between `x=∓√δ`:

```
ΔU(δ) = (4/3) δ^{3/2}     ⇒     K_fold = 4/3   (exact; = (4/3) for any cubic
                                  fold with ½f''=1, so FHN and VdP alike).
```

**`A₀` — the Kramers prefactor.** `√(U''_min|U''_max|)/(2π) = 2√δ/(2π)=√δ/π`:

```
λ_esc(δ) = (√δ / π) · exp( −(8/3) δ^{3/2}/σ² )   ⇒   A₀ = 1/π .
```

So the two constants the audit left open are **`K_fold = 4/3`, `A₀ = 1/π`** — pure
geometry, exactly as promised, now actually computed.

**Assembly (same collapse).** Slow drift carries `δ` at rate `|dδ/dt|=ε g`,
`g=|v+a−bw|_fold=|a−1−b(I−2/3)|`. Then `dt=dδ/(εg)` and the *same* substitution
cancels the `√δ` prefactor:

```
H = (1/(εg)) ∫₀^∞ (√δ/π) e^{−(8/3)δ^{3/2}/σ²} dδ = σ²/(4π ε g) ,
H=ln2 ⇒  σ_crit = √(4π g ln2) · √ε .
```

For `(a,b,I)=(0.7,0.8,−0.1)`: `g = |0.7−1−0.8(−0.767)| = 0.314`, giving
`σ_crit ≈ √(4π·0.314·0.693)·√ε = 1.65 √ε`. In the audit's `√(ε/C)` form,
`C = 1/(4π g ln2) ≈ 0.37`.

### Scope and honesty for #12

- **Rock-solid:** `K_fold=4/3` and `A₀=1/π` — the near-fold barrier coefficient
  and Kramers prefactor, derived from the cubic geometry, no fitting. These close
  the literal "computable but uncomputed" gap.
- **The assembled `σ_crit=√(4π g ln2)√ε` holds in the drift-toward-fold regime**
  (a net `O(ε)` slow passage through the fold — near-SNIC / ramped). For the
  *canonical* stable-rest excitable (`I=−0.1`, FP at `v=−1.26`, finite distance
  from the fold) the escape is Kramers from the FP over the *global* barrier
  (`ΔU(W_FP)≈0.0254`); the near-fold `δ^{3/2}` law is then the `I→I_SNIC`
  asymptotic, where the FP approaches the fold and `δ→0`. State which regime when
  quoting the number.
- **Same prefactor caveat as Part 1 item 1** (Kramers vs exact rate near the
  dominant `δ`). Validation against the measured excitable failure boundary is a
  (numerical) follow-up.

---

## The unification (the reason this was the high-return move)

The canard (#1), VdP (#16), and excitable (#12) noise thresholds are **one
formula**:

```
σ_crit = √(4π · D · ln2) · √ε ,
```

with `D` the deterministic slow-drift speed *at the fold* — `λ` for the canard,
`g` for the excitable — and `√(4π ln2)` a universal constant from the cubic
fold's quasipotential. The exponent `√ε` and the `√(4π ln2)` prefactor are
regime-independent; only the geometric drift `D` changes. One Freidlin–Wentzell
calculation, three audit rows.

| Item | `D` (drift at fold) | `σ_crit` | status |
|---|---|---|---|
| #1 canard | `λ = b(I−I_fold_L)` | `√(4π ln2)·√(ε λ)` | certified (leading order) |
| #16 VdP | `λ = 1+a` | same form | universal half rides #1 |
| #12 excitable | `g = a−1−b(I−2/3)` | `√(4π g ln2)·√ε` | `K_fold=4/3`, `A₀=1/π` derived |

**Remaining to make it 3-figure-rigorous (one computation, shared by all three):**
the exact cubic-barrier escape rate replacing the Kramers prefactor near
barrier`/η² = O(1)` — it renormalises the universal `√(4π ln2)` by a computable
`O(1)` factor and is the single thing standing between "clean leading-order" and
"asymptotic theory" for the whole fold-escape cluster.
