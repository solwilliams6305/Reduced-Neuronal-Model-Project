# PROGRAM2 — TW_β halo: the exact Hastings–McLeod Stokes constant (T2c closed form)

**Date:** 2026-07-22. **Executes:** `regime-tests/PROGRAM2_HANDOFF_C_TW_BETA_HALO.md`, open
tasks **#1** (verify 5 classical citations) and **#2** (exact `C`; the backbone half). **Memory:**
`program2-tw-beta-halo.md`.

---

## Result (one line)

The large-order / Stokes constant of the Hastings–McLeod Painlevé II backbone — which by Theorem T1
governs the resurgent `TW_β` tail — is, **certified**,

$$\boxed{\;a_k \sim -C\,\Gamma\!\big(2k-\tfrac12\big)\,(9/8)^k,\qquad \gamma=-\tfrac12,\qquad
C=\sqrt{\tfrac{2}{3\pi^{3}}}=\tfrac1\pi\sqrt{\tfrac{2}{3\pi}}=0.146632271193848478\ldots\;}$$

**This refutes the paper's conjectured closed form `C ≟ 1/(4√π) = 0.14105…` (3.81 % too small).**
The earlier "`C ≈ 0.139`" (`_hm_stokes.py`, `N=10`) was simply an under-converged estimate; the true
value is `0.14663…`, *further* from `1/(4√π)`, not closer.

The value **equals the published HM Stokes constant** (Dunne, *Introductory Lectures on Resurgence*,
arXiv:2511.15528, eq. (2.16): `S = (1/π)√(2/(3π)) = √(2/(3π³))`). Independent computation and
independent literature agree to **28 digits**.

---

## What was done

### 1. Exact coefficients — a closed-form rational recursion (replaces the N=10 sympy solve)

Substituting `q = 2^{-1/2} Σ_k a_k u^{1/2-3k}` (`u=-s`) into the HM equation `q_uu = 2q³ - u q`
(α=0) and collecting the coefficient of `u^{3/2-3m}` gives, for `m ≥ 1`,

```
a_m = ½ [ a_{m-1} · p · (p-1)  −  T(m) ],     p = ½ − 3(m−1),
T(m) = Σ_{i+j+l=m, 0≤i,j,l≤m−1} a_i a_j a_l    (cubic self-convolution, m-th index excluded),
a_0 = 1.
```

Hand-checks: `a_1 = −1/8`, `a_2 = −73/128` (matches the paper and `_hm_coeffs.json` exactly). This is
`O(N³)` in exact `fractions.Fraction`, fast to `N = 120` (and extendable). Coefficients are
**non-sign-alternating** (`a_{k≥1} < 0`) — the Cleri–Dunne signature ⇒ positive-real-axis Borel
singularity ⇒ imaginary lateral ambiguity ⇒ median = real tail. Certified against the old
`_hm_coeffs.json` (`k<10`, rel error `0`).

### 2. Method A — Richardson/Neville extraction (the certified numbers)

- **γ:** consecutive-ratio estimator `γ_k = [3+√(1+4ρ_k)]/2 − 2k`, `ρ_k = (a_k/a_{k−1})/(9/8)`,
  Neville-extrapolated in `1/k`. → `γ* = −0.5` with `|γ*+½| ≈ 6·10⁻³¹`, window spread `~10⁻²⁵`.
  **γ = −1/2 exactly** (25+ digits).
- **C:** `C_k = −a_k / [Γ(2k+γ)(9/8)^k]`, Neville-extrapolated. → `C* = 0.146632271193848477888…`,
  window spread `1.9·10⁻²⁸`. Identical for `γ=γ*` and `γ=−1/2` (as it must).

### 3. Closed form + three independent validations

- **(LIT)** `√(2/(3π³)) = 0.146632271193848477888509285742`; `C* − √(2/(3π³)) = 5·10⁻³⁰` (the
  extrapolation floor) ⇒ **match to ~28 digits**. Also `C·π = √(2/(3π))` to 20 digits. Published
  constant = Dunne 2511.15528 eq. (2.16); classical origin Its–Kapaev / Kapaev nlin/0411009; textbook
  FIKN ch. 11; connection review Bothner arXiv:2003.14374.
- **(BOR)** Borel-space cross-check. `c_k = a_k/(2k)! ~ −C·2^{-3/2}k^{-3/2}(9/8)^k`; the Borel
  transform `Σ c_k x^k` (`x=ξ²`) has a **bounded** `(x₀−x)^{1/2}` branch at `x₀ = 8/9 = 1/A²`.
  Extracted branch exponent `−1.5000000001` (⇒ γ=−1/2) and amplitude `C = 0.1466322711939…` (matches
  to 13 digits) — a different reprocessing (Borel space, via Stirling) than method A.
- **(ODE)** Direct HM PII integration (scipy DOP853, boundary `q ~ Ai(s)` at `s→+∞`). At
  `s = −2,…,−6` the optimally-truncated series reproduces the true transcendent with remainder
  tracking the least term and **decaying at the action rate `e^{−A(−s)^{3/2}}`, A = 2√2/3** — the
  `a_k` are the genuine HM asymptotics (median resummation reproduces HM up to the exp-small ambiguity).

### 4. Citation audit (handoff task #1) — DONE

Web-verified against primary sources (Springer/Euclid/ADS/AMS/publisher). Two need fixing (both now
patched in `TWbeta_Resurgent_skeleton.tex`):

| ref | verdict |
|---|---|
| Tracy–Widom 1994, CMP **159** 151–174 | ✅ correct |
| **Flaschka–Newell 1980**, CMP **76** 65–116 | ⚠️ title punctuation: "…deformations**. I**" (was "…deformations I") |
| **Jimbo–Miwa–Ueno 1981**, Physica D **2** 306–352 | ⚠️ dropped subtitle: "…rational coefficients**. I. General theory and τ-function**" |
| FIKN 2006, Math. Surveys Monogr. **128** | ✅ correct |
| Clarkson 2003, JCAM **153** 127–140 | ✅ correct (em-dash cosmetic) |

Bonus: **HM Stokes multipliers `(s₁,s₂,s₃) = (−i, 0, i)` confirmed** (Bothner §8, γ=1 boundary of the
Ablowitz–Segur family; PII constraint `s₁−s₂+s₃+s₁s₂s₃ = −2sin πα = 0` at α=0 ✓). Exponentially-small
correction: `δq(s) ~ (2√(2π))^{-1}(−s)^{1/4} exp(−(2√2/3)(−s)^{3/2})`, amplitude `1/(2√(2π))`.

---

## Impact on the paper

`TWbeta_Resurgent_skeleton.tex` updated (rebuilds clean, 6 pp):
- Theorem T2(c) amplitude `≈0.14` → **`C = √(2/(3π³))` (exact)**.
- §T2c-proof: the `C≈0.139≈1/(4√π)` line → the certified recursion + `γ=−1/2` + `C=√(2/(3π³))` +
  the refutation of `1/(4√π)` + the two cross-checks. New cites `\cite{Dunne2025,KapaevHM,BothnerRHP}`.
- Remark (i): "exact `C`" removed from the open list (now fixed); only the **a-dependence of `S(a)`**
  remains as the refinement.
- Status table: T2c row → `C=√(2/(3π³))` (28 digits).
- Bibliography: FN1980 + JMU1981 corrected; added Dunne2025, KapaevHM, BothnerRHP.

**Net:** handoff task #1 fully closed; task #2 **backbone half closed exactly** (the `a`-independent
`C`); only the tail-specific `a`-dependence `S(a)` (transport hierarchy) remains of task #2.

---

## Files

- `coupled-atlas/_hm_stokes_exact.py` — exact rational recursion, method-A Richardson, PSLQ probe,
  saves `_hm_coeffs_exact.json` (121 exact coeffs + `C_str`).
- `coupled-atlas/_hm_stokes_validate.py` — (LIT) closed-form cert, (BOR) Borel-space, (ODE) direct HM.
- `coupled-atlas/_hm_coeffs_exact.json` — 121 exact `a_k` (as `num/den`) + `C` to 30 digits.
- Supersedes the guess in `coupled-atlas/_hm_stokes.py` / `O5d_stokes_computed.tex` (`C≈0.139`).

## Gotchas (this session)

- **`float(Fraction)` overflows** for `k ≳ 50` (`a_k` has hundreds of digits) — keep exact as
  `num/den` strings; convert to `mpf` via `mpf(num)/mpf(den)`.
- **Borel–Padé via the derivative is unreliable here**: `γ=−½<0` makes the branch point *bounded*
  (`B` finite, only `B'` diverges), so `(x₀−x)^{1/2}B'` from a `[60/60]` Padé is noisy near the cut.
  The robust Borel check is the **coefficient** asymptotic `c_k ~ k^{-3/2}` (leg BOR), not Padé.
- **mpmath `odefun` choked** on the long `+9 → −7` integration at `tol=1e-60`; scipy `DOP853`
  (`rtol=1e-13`) is plenty for `s∈[−2,−6]` (remainder `10⁻³…10⁻⁸`).
- **PSLQ needs full-precision input**: feeding `C*` (good to ~28 digits) at `dps=80` failed to find
  `[2,−1,1,3]` (the relation `2 lnC − ln2 + ln3 + 3 lnπ = 0`) because the trailing digits are
  extrapolation noise. The literature match is what nailed the closed form; PSLQ would need `C*`
  trusted to ≥ its own precision.

## Remaining (unchanged priority from handoff)

2′. **a-dependence of `S(a)`** — the tail-specific Stokes constant beyond the backbone `C`; needs the
    transport hierarchy `c_n(a)` (`O1_reduction_note.tex` §3). The `a`-independent piece is now done.
3.  α=1/2 vs α=0 assignment (saddle vs distribution). 4. rough-potential rigor. 5. direct exact-WKB.
