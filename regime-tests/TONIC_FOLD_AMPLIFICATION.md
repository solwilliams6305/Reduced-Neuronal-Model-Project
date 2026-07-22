# Mid-tonic fold-passage amplification: a first-principles audit of `A_fold ≈ 2.7`

**Status:** derived + numerically settled. This closes the first-principles
half of TONIC_PHASE.md §20 / §21.1-item-2, but with a **correction**: the
clean cycle-geometric constant `A_fold ≈ 2.7` asserted in §20 does **not**
survive a direct first-principles computation. The σ = 0.02 mid-tonic CV
enhancement is the **onset of non-perturbative fold escape**, not a
deterministic multiplicative factor — vindicating §18.2 over §20.

Companion script: `tonic_fold_amplification.py`
Companion figure: `figures/tonic_fold_amplification.png`

---

## 1. The claim under audit

TONIC_PHASE.md §20 closed "Open Q3" (the σ = 0.02 mid-tonic CV ratio of 2.7)
by asserting a **deterministic Lyapunov amplification**

```
A_fold  =  exp( ∮ max(λ⊥(t), 0) dt ),   λ⊥(t) = (1 − γ_v(t)²) − εb,      (20.3)
CV^eff / CV_leading  =  A_fold  =  2.7   (σ-independent).                  (20.5)
```

with `A_fold ≈ 2.7` obtained by **inverting** the measured CV ratio
(`∮max(λ⊥,0)dt = ln 2.7 ≈ 0.99`). The §21.1 roadmap then flagged the
*first-principles* evaluation of (20.3) from the cycle's `V(T)` profile as a
"one-page calculation." This document does that calculation. It does not give
2.7. Here is why, and what the correct statement is.

## 2. The cycle is overwhelmingly net-contracting

`λ⊥(t) = (1 − γ_v²) − εb` is the transverse Floquet rate (the non-phase
eigen-rate; the tangential direction contributes zero to `tr J` by Floquet's
theorem). Integrating it around the deterministic mid-tonic cycle gives the
non-trivial Floquet exponent `λ2 = (1/T) ∮ tr J dt`:

| ε | T_cycle | λ2 = ∮trJ dt / T | ∮ trJ dt (= λ2·T) | Floquet multiplier e^{λ2 T} |
|---|---|---|---|---|
| 0.04 | 62.84 | −1.186 | −74.5 | 4.3×10⁻³³ |
| 0.08 | 36.42 | −1.025 | −37.3 | 6.2×10⁻¹⁷ |
| 0.16 | 22.24 | −0.837 | −18.6 | 8.2×10⁻⁹ |

The cycle is **extraordinarily stable transversally**: a transverse
perturbation is contracted by 16–33 orders of magnitude per period. This is
the first problem with (20.3): it reads a *permanent* amplification off the
expanding part of `λ⊥` while discarding the (far larger) contracting part.

## 3. The exponentiated expanding integral is O(10–10²), not 2.7

Computing the positive part of `λ⊥` directly from the cycle:

| ε | ∮ max(λ⊥,0) dt | exp(∮ max(λ⊥,0) dt) |
|---|---|---|
| 0.04 | 5.00 | 149 |
| 0.08 | 4.11 | 61 |
| 0.16 | 3.21 | 25 |

So the literal evaluation of (20.3) gives `A_fold` between 25 and 149 — not
2.7, and *ε-dependent*. A transverse fluctuation is not permanently amplified
by this factor; it is amplified during the fast fold transit and immediately
re-contracted by the strongly stable slow branches (§2). The Floquet
decomposition already accounts for the net effect, and the net is contraction.

## 4. The "post-fold-tip" inner integral is cutoff-dependent

A charitable reading of §20.1–20.2 restricts the integral to the brief
*expanding* sliver just past the fold tip, where in Krupa–Szmolyan blow-up
coordinates (`v = −1 + ε^{1/3}V`) one has `λ⊥ = 2ε^{1/3}V − εb + …` so that
`λ⊥ > 0` for `V > 0`. The amplification over that sliver,

```
∫ max(λ⊥,0) dt  =  2 ∫_0^{V_out} V dT     (ε-independent; dt = ε^{−1/3}dT),
```

is computed from the deterministic inner passage `dV/dT = V² − W, dW/dT = −λ_rel`
(here `λ_rel = b(I_mid − I_fold_L) = 0.467`):

| V_out | 2∫₀^{V_out} V dT | exp |
|---|---|---|
| 0.5 | 0.272 | 1.31 |
| 1.0 | 0.750 | 2.12 |
| 1.5 | 1.218 | 3.38 |
| 2.0 | 1.632 | 5.12 |
| 3.0 | 2.304 | 10.0 |

The integral **grows without bound** as the upper cutoff `V_out` increases
(the runaway branch has `V → ∞`). It equals `ln 2.719 = 0.99` only at the
tuned value `V_out ≈ 1.26`. There is no first-principles reason to stop at
`V_out = 1.26` rather than 2 or 3, so **2.7 is not a convergent
cycle-geometric number** — it is the inversion target dressed as a cutoff.

## 5. The measured CV ratio is σ-dependent → not a deterministic factor

A genuine deterministic `A_fold` (eq. 20.5) must be **σ-independent**: the
ratio `CV_meas/(σ·A)` would be flat in σ. It is not. The §10 data already
show it climbing (1.0 → 2.7 → 14 at σ = 0.005, 0.02, 0.05). A fresh stochastic
ISI sweep at the §10 operating point (I = 0.83, ε = 0.08, A = 0.3989) sharpens
this:

| σ | CV_meas | CV/(σ·A) | regime |
|---|---|---|---|
| 0.005 | 0.0020 | 0.99 | leading order CV = σ·A ✓ |
| 0.015 | 0.0058 | 0.98 | leading order CV = σ·A ✓ |
| 0.020 | 0.0332 | 4.16 | **broken** |
| 0.025 | 0.0406 | 4.08 | broken |
| 0.030 | 0.1448 | 12.1 | broken |

The leading-order phase reduction holds *cleanly* (ratio ≈ 1) up to
σ ≈ 0.015, then the ratio **jumps sharply** to 4–12. This is the fingerprint
of rare, large, non-perturbative excursions (escape-like fold-passage events)
contaminating the ISI sample — not a smooth multiplicative correction. The
exact ratio just above onset is **sampling-dependent**: §10 reported 2.72 at
(I = 0.83, σ = 0.02); this run gives 4.16 at the same point and an earlier
private run gave 0.91 at I = 0.875. A deterministic constant would not move
with the seed or with a 0.045 shift in I.

## 6. The correct statement

Putting §2–§5 together:

> **There is no deterministic, σ-independent `A_fold ≈ 2.7`.** The transverse
> dynamics of the mid-tonic cycle are governed by the (large, negative)
> Floquet exponent `λ2` (§2); the principled finite-noise correction is the
> second-order phase-reduction coefficient `K₂` of §16, which TONIC_PHASE.md
> §17.2 computed and found *convergent but small* (`σ²K₂/⟨Z_v²⟩ ≈ 0.15` at
> σ = 0.02). The factor-of-a-few CV enhancement that appears at σ ≈ 0.02 is
> the **onset of non-perturbative escape from the cycle during fold passages**,
> at a threshold `σ_onset ≈ 0.015–0.02`. Below `σ_onset`, CV = σ·A holds with
> no fold correction; above it, the ISI distribution acquires a heavy tail
> from rare excursions and CV = σ·A·(deterministic factor) is the **wrong
> functional form**. The right framework is Freidlin–Wentzell escape from the
> limit cycle (§18.2): a contribution scaling as `exp(−A_action/σ²)` times a
> rate prefactor, which is exactly what produces a *sharp threshold* in σ
> rather than a smooth power-law correction.

This is the same escape mechanism as the canard chapter's `σ_*`, now applied
to the tonic cycle's transverse direction. It **unifies** the tonic and canard
chapters under the structural finding already in the README §1: *noise destroys
the deterministic structure where transverse stability is weakest* — here, at
the fold passages, where `λ⊥` momentarily turns positive even though the cycle
is globally hyper-stable.

## 7. What is genuinely first-principles here

| Quantity | First-principles value | Status |
|---|---|---|
| Transverse Floquet exponent `λ2` | −1.19 / −1.03 / −0.84 (ε = 0.04/0.08/0.16) | computed from `∮ trJ dt` |
| Expanding integral `∮ max(λ⊥,0) dt` | 5.0 / 4.1 / 3.2 | computed; exponentiation overcounts |
| Post-tip inner integral `2∫V dT` | cutoff-dependent (diverges) | computed; not convergent |
| `σ_onset` (escape onset, mid-tonic) | ≈ 0.015–0.02 | measured (fresh sweep) |
| `A_fold = 2.7` as a deterministic constant | **not reproduced** | refuted (§3–§5) |

The deliverable of "first-principles `A_fold`" is therefore a **negative
result with a mechanism**: the clean Lyapunov factor does not exist; the
correct object is the FW escape action that sets `σ_onset`. Computing that
action (the cheapest path leaving the cycle near a fold and returning on the
next pass) is the genuine PhD-scope analytical step, and it is the *tonic*
instance of the canard chapter's blow-up + large-deviation calculation.

## 8. Consequence for the σ-validity table

TONIC_PHASE.md §20.5 / README §5.6 list a three-regime structure with the
linearised+`A_fold` band running to σ ≈ 0.1. That band should be **deleted**.
The corrected structure at mid-tonic is two regimes:

```
σ ≤ σ_onset ≈ 0.015–0.02 :  leading-order phase reduction, CV = σ · A  (no fold factor)
σ > σ_onset              :  non-perturbative fold escape (FW action; heavy-tailed ISI)
```

The off-centre cells (I = 0.475, 1.275) have `σ_onset` larger (their cycles do
not make sharp fold passages), which is why §10 saw them stay perturbative at
σ = 0.02 while mid-tonic did not. `σ_onset(I)` — smallest at mid-tonic, larger
toward the Hopf edges — is the clean, defensible observable, and it is set by
fold-passage transverse geometry exactly as the "ghost of the canard" picture
(§5) predicts.

## 9. Reproduce

```
python3 regime-tests/tonic_fold_amplification.py
```

Computes the Floquet exponent `λ2`, the expanding integral and its
exponentiation, the cutoff-dependence of the post-tip inner integral, the §10
σ-ratios, and a fresh stochastic ISI onset sweep; writes
`results/tonic_fold_amplification/summary.txt` and
`figures/tonic_fold_amplification.png`.
