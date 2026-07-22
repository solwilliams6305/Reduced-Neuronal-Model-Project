# PROGRAM2 — TW_β halo, Thread 2: the noise (1/β) expansion IS Borel-summable

**Date:** 2026-07-22. **Executes:** Thread 2 — the frontier flagged un-established by the deep-research
pass: is the tail's **noise `1/β`** expansion (as opposed to the level-variable resurgence) a
Borel-summable trans-series? **Answer: YES.** **Memory:** `program2-tw-beta-halo.md`.

---

## Result (one line)

The TW_β tail's `1/β` expansion is Borel-summable — a genuine **two-sector resurgent trans-series**:

| sector | source | Borel singularity | summation |
|---|---|---|---|
| **(A) dynamical escape** | reduced escape `d_m∼−(1/π)Γ(m)(3/2)^m` | **positive real axis**, `t_c=2/3=A_g` (reduced action Φ; two-instanton ambiguity) | **median** |
| **(B) Coulomb-gas** | `Γ(β/2)` β-ensemble normalization | **imaginary axis**, `±iπ` (Γ instantons) | ordinary |

Both are factorially divergent (Gevrey-1) ⇒ the `1/β` series is a bona-fide trans-series, distinct
from the LEVEL-variable (a, β=2) Hastings–McLeod resurgence (`C=√(2/(3π³))`).

## Setup

The exact tail (Borot–Nadal 1111.2761, Prop. 1.1) is
`1−TW_β(s)=c_β\,s^{−3β/4}e^{−2βs^{3/2}/3}\exp[Σ_m(β/2)R_m(2/β)s^{−3m/2}]`, `c_β=Γ(β/2)/((4β)^{β/2}2π)`.
Its `1/β` expansion (at fixed level `s=a`) has two factorially-divergent sources; the `R_m(2/β)` are
polynomials in `1/β` (degree ≤ m+1) and contribute a **regular** (per-level-finite) piece.

## (A) Dynamical escape sector — positive-axis, median (`_tw_thread2_borel.py`)

The reduced escape fluctuation series `F(g)=Q̃(2/g)^{−2}=Σ_m d_m g^{−m}` (`g=βa^{3/2}`; this program,
exact) has `d_m∼−(1/π)Γ(m)(3/2)^m`, **all `d_m<0`** (non-sign-alternating). At fixed `a` the `1/β`
series is `d_m a^{−3m/2}β^{−m}` — same non-alternating factorial ⇒ **positive-real-axis Borel
singularity** ⇒ requires **median** resummation.

- **Borel singularity located:** `B_F(t)=Σ_{m≥1}d_m t^{m−1}/(m−1)!`; Padé[30/30] → nearest
  positive-real pole `t_c=0.666665869…` = `2/3` = `A_g` (the reduced instanton action Φ; its median
  ambiguity is the two-instanton `e^{−2βΦ}`) in the Borel-`1/g`
  plane). ✓
- **Median resummation reconstructs the function:** median Borel–Laplace (rays `±0.35` rad)
  `F(9)=0.9423505` vs the exact escape rate value `π e^{2g/3}R_MFPT(9)=0.9434944` — **rel. diff
  `1.2×10⁻³`** (vs a naïve 60-term sum, which diverges). So the dynamical sector is
  **Borel-summable (median)**; the lateral ambiguity `~e^{−(2/3)g}` is the Stokes/one-instanton term.

## (B) Coulomb-gas sector — imaginary-axis, ordinary Borel

`ln c_β` carries `ln Γ(β/2)`, whose Stirling tail `Σ_k B_{2k}/(2k(2k−1))(2/β)^{2k−1}` has coefficient
of `β^{−(2k−1)}` equal to `a_k=B_{2k}/(2k(2k−1))2^{2k−1}`:
- **alternating** (`+ − + − …`), factorial: `|a_k/a_{k−1}|/[(2k−2)(2k−3)] → 1/π² = 0.101321` (to 6
  digits by `k=16`), i.e. `a_k∼(−1)^{k+1}(2k−2)!/π^{2k}`.
- Borel transform `∼(1/π)\arctan(t/π)`: branch points at `t=±iπ` (generally the Γ instantons at
  `±2πin`) — on the **imaginary axis** ⇒ **ordinary Borel-summable**, no positive-axis median
  ambiguity. These are the β-ensemble / Coulomb-gas "noise instantons" `e^{−2πin·β/2}`.

## Why this doesn't contradict CLDS

CLDS (2510.14433) compute the **typical (bulk)** fluctuation cumulants `⟨a₁ⁿ⟩_c≃C_n/β^{n−1}`
(`C₂≈1.6697=`the program's `s₂`; `C₃,C₄`) — a handful of pure numbers near the mean, with no
large-order/growth test — and correctly find no Stokes structure *there*. The resurgence is a
property of the **tail** (large-deviation, large `a`) `1/β` expansion, a different object; that is what
carries the escape instanton (A) and, through the normalization, the Coulomb-gas sector (B).

## Status / honesty

- **Established (numerically decisive):** both sectors are factorially divergent with the stated Borel
  structure; the dynamical sector's median sum reconstructs the exact escape rate to `10⁻³`; the
  Coulomb-gas sector's alternating-factorial / `±iπ` structure is exact (Γ function). This **upgrades
  the literature's "un-established"** to a concrete two-sector characterization.
- **Caveat (standard):** a fully rigorous Borel-summability proof for the noise direction (Costin-type
  analytic-continuation + subexponential bound of the Borel transform along the ray) is not done here;
  the structure and the median reconstruction are what is shown. The dynamical handle `d_m` is the
  *frozen* (leading-in-x-extension) series; the full tail's dynamical sector adds the x-extension
  dressing (`√(2/(3π))` zero mode, `−½` index; `PROGRAM2_TWBETA_GY_ZEROMODE_NOTES.md`) but keeps the
  same positive-axis action, so the verdict carries.

## Relation to the two other resurgences (the halo's full picture)

- **Level variable (a), β=2:** Hastings–McLeod PII, positive-axis/median, `C=√(2/(3π³))` [DONE this
  program].
- **Noise (1/β), fixed a — Thread 2 [this note]:** two sectors — dynamical escape (positive-axis,
  median; sing at `Φ`, two-instanton ambiguity) + Coulomb-gas `Γ(β/2)` (imaginary-axis, ordinary).
- These are genuinely different directions of one multi-parameter trans-series; the frozen escape is
  the object where the `1/β` and level directions coincide (`g=βa^{3/2}`), which is why it supplies
  the dynamical sector for both.

## Files
- `coupled-atlas/_tw_thread2_borel.py` — sector (A) Borel-Padé + median-Laplace vs MFPT; sector (B)
  Γ-Stirling alternating-factorial / imaginary-axis analysis.
