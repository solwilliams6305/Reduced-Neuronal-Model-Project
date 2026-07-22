# PROGRAM2 — TW_β halo: the √(2/(3π)) "x-extension factor" is a zero-mode Jacobian

**Date:** 2026-07-22. **Executes:** the residual of `PROGRAM2_HANDOFF_C_TW_BETA_HALO.md` task #2 —
derive the elementary factor relating the frozen-escape Stokes constant `S₀=1/π`
(`PROGRAM2_TWBETA_S_OF_A_NOTES.md`) to the full Hastings–McLeod backbone
`C=√(2/(3π³))` (`PROGRAM2_TWBETA_STOKES_CONSTANT_NOTES.md`). **Memory:** `program2-tw-beta-halo.md`.

---

## Result (one line, and a correction)

The frozen→full change is **not** a lone multiplicative factor (as the previous note loosely put it).
It is a **joint shift of the Stokes data** — constant *and* resurgent index together:

$$\Big(S_0=\tfrac1\pi,\ \gamma'=0\Big)\ \xrightarrow{\ \text{x-extension}\ }\ \Big(C=\sqrt{\tfrac{2}{3\pi^3}},\ \gamma=-\tfrac12\Big),\qquad
\frac{C}{S_0}=\sqrt{\tfrac{2}{3\pi}}=\sqrt{\tfrac{A_g}{\pi}},\ \ A_g=\tfrac23 .$$

The constant factor is exactly `√(A_g/π)` with `A_g=2/3` the **one-instanton action** (`e^{−(2/3)g}`),
i.e. the **x-translation zero-mode (collective-coordinate) Jacobian**; and it comes *bundled* with the
`−½` index shift. Constant + index shift together are the signature of **one continuous zero mode** —
precisely the instanton-translation mode the frozen (single-x) escape lacks and the x-extension adds.
**Confirmed numerically to 9 digits.**

## The concrete diagnostic (`coupled-atlas/_tw_frozen_vs_hm.py`)

Dictionary (matched two-instanton exponentials, `PROGRAM2_TWBETA_S_OF_A_NOTES.md`):
`g = βa^{3/2} = √2\,(−s)^{3/2}`, so `(−s)^{−3}=2/g²`. Re-expanding the HM series
`q∼√(−s/2)Σ_k a_k(−s)^{−3k}` in `g` lands it at **even** orders `g^{−2k}` with coefficient `a_k 2^k`.
Compare to the frozen loop series `F(g)=Σ_m d_m g^{−m}` at `m=2k`. Both established growths:

```
frozen :  d_{2k}   ~ -(1/π) Γ(2k)     (9/4)^k       (γ'=0 ; base 3/2 ⇒ (3/2)^{2k}=(9/4)^k)
HM     :  a_k 2^k  ~ -C     Γ(2k-1/2) (9/4)^k       (γ=-1/2 ; base 9/8 ⇒ (9/8)^k 2^k=(9/4)^k)
```

**Same action `(9/4)^k`** — so the x-extension does *not* move the Borel singularity. It changes only
the index (`Γ(2k)→Γ(2k−½)`) and the constant. Hence the ratio

$$r_k:=\frac{d_{2k}}{a_k 2^k}\sim \frac1{\pi C}\,\frac{\Gamma(2k)}{\Gamma(2k-\tfrac12)}\sim \sqrt{\tfrac{3\pi}{2}}\,\sqrt{2k},
\qquad \frac{r_k}{\sqrt{2k}}\to \frac1{\pi C}=\sqrt{\tfrac{3\pi}{2}}=\frac1{\sqrt{2/(3\pi)}}.$$

Numerics (`k` up to 80): `r_k/√(2k) → 2.170803762`, versus `1/(πC)=√(3π/2)=2.170803764` — **agree to
`1.4×10⁻⁹`**. So `C/S₀=√(2/(3π))` and the index shift `−½` are both verified.

## Why `√(A_g/π)` — the zero-mode reading

`√(2/(3π)) = √((2/3)/π) = √(A_g/π)`, `A_g=2/3` the one-instanton action, is the standard
**collective-coordinate Jacobian** for a single translation zero mode: trading the zero mode of the
x-extended instanton for its collective coordinate multiplies the Stokes constant by `√(A_g/π)` and
lowers the resurgent index by `½` (one continuous mode ⇒ one half-integer in the Γ). This is exactly
the `(√(2/(3π)), −½)` pair we measure.

**This is the SAME zero mode the O1 note already flagged.** `O1_reduction_note.tex` §3: "`A₀` … its zero
mode (instanton translation) is removed within the Riccati formulation as in Schorlepp–Grafke–Grauer,
leaving the finite `a^{−3β/4}` prefactor." So the single x-translation zero mode supplies **both** the
`a^{−3β/4}` eikonal prefactor (O1 note) **and** the `√(2/(3π))` Stokes-constant factor (here). The frozen
escape has no x-direction, hence no such mode, hence its cleaner `(1/π, γ'=0)`.

## Literature confirmation (deep-research pass this session)

The decomposition gives a *physical reading of the published constant*. **Dunne, Introductory Lectures
on Resurgence (arXiv:2511.15528), eq. (2.16)** writes the HM Stokes constant in exactly the factored form
```
c_n^{(0,-)} ~ -(1/π)√(2/(3π)) · Γ(2n-1/2)/((2√2/3)^{2n}) · (1 - (17/72)/(2n-3/2) + (1513/10368)/((2n-3/2)(2n-5/2)) - …)
```
i.e. `C=(1/π)√(2/(3π))` is the literature's own form — and our two pieces have physical content:
`1/π` = the cubic-barrier **escape-rate** Stokes constant (computed from scratch, 25 digits), `√(2/(3π))`
= the **x-translation zero-mode Jacobian** `√(A_g/π)`. Our `a_n` reproduce (2.16) not only at leading
order but at the **first subleading** order: forming `R_n = a_n/[-C\,Γ(2n-1/2)(9/8)^n]`, the extracted
`(R_n-1)(2n-3/2)` Richardson-extrapolates to **−0.23611111111111 = −17/72 to 16 digits** (Dunne's
coefficient). So `a_n` ARE Dunne's `c_n^{(0,-)}` — the HM level-variable coefficients.

**Level vs. noise — a caveat the literature makes sharp.** This resurgence lives in the **level
variable at the fixed classical point β=2** (Hastings–McLeod PII; Dunne, Cleri–Dunne 2002.06270). It is
**not** a large-β/noise object. Comtet–Le Doussal–Smith (2510.14433) treat the `1/β` corrections as
ordinary **cumulants** (pure numbers `C₂≈1.6697` — matching the program's weak-noise variance `s₂=1.67`
from `_o1_falsifier.py`; `C₃≈0.7438`, `C₄≈0.5576`), with **no** Borel/Stokes/factorial structure; and
Borot–Nadal (1111.2761, Prop. 1.1) give the exact all-β soft tail
`1−TW_β(s)=[Γ(β/2)/((4β)^{β/2}2π)]\,s^{-3β/4}e^{-2βs^{3/2}/3}\exp[Σ_m(β/2)R_m(2/β)s^{-3m/2}]` — a
**single-exponential level series, no second instanton**. So the program's "`1/β` trans-series"
(Thread 2) is genuinely un-established in the literature; the solid, literature-confirmed resurgence is
the level-variable HM one, whose constant we have now both computed (C) and physically decomposed
(1/π × zero-mode Jacobian). This scopes the frontier honestly.

## Status: what is derived, what remains

- **Derived / confirmed (9 digits):** the frozen→full map is `(1/π,0)→(C,−½)`, with constant ratio
  exactly `√(A_g/π)=√(2/(3π))` (the one-instanton-action zero-mode Jacobian) and index shift `−½`
  (one zero mode). This *explains the value* of the HM Stokes constant `C=√(2/(3π³))=(1/π)·√(A_g/π)`
  as (cubic-barrier escape constant) × (x-translation zero-mode Jacobian).
- **Already rigorous, independently:** `C=√(2/(3π³))` itself is fixed by the T1 Painlevé-II reduction
  (`PROGRAM2_TWBETA_STOKES_CONSTANT_NOTES.md`), so the tail Stokes constant needs no separate
  derivation — this note *interprets* its value.
- **Residual — NOW CLOSED (`PROGRAM2_TWBETA_GY_ZEROMODE_NOTES.md`).** The Gelfand–Yaglom zero-mode
  Jacobian is `√(‖ψ₀‖²/2π)` with `ψ₀=−sech²τ` the tanh-kink escape instanton (`P⋆=−tanh τ`) and
  `‖ψ₀‖²=∫sech⁴=4/3=ΔV=2A_g` (virial). So `√(A_g/π)=√(2/(3π))` — the factor is `√(A_g/π)` **because
  `‖ψ₀‖²=ΔV=2A_g`**, i.e. the `2` is the virial factor (not a boundary artifact). Assembles to
  `C=(1/π)√(2/(3π))=√(2/(3π³))` to 40 digits.

## Correction to `PROGRAM2_TWBETA_S_OF_A_NOTES.md`

That note said the x-extension "dresses `1/π → C`, contributing the elementary `√(2/(3π))`." Accurate
refinement: it is a **joint** shift of the Stokes data `(1/π, γ'=0) → (C, γ=−½)`; the constant part is
`√(2/(3π))=√(A_g/π)` (zero-mode Jacobian) and it is inseparable from the `−½` index shift — one
continuous (x-translation) zero mode supplies both. The action is unchanged (`(9/4)^k` at matched
orders), so this is a determinant/index effect, not a Borel-singularity shift.

## Files
- `coupled-atlas/_tw_frozen_vs_hm.py` — the matched-order diagnostic (`r_k/√(2k)→1/(πC)=√(3π/2)`, 9 digits).
