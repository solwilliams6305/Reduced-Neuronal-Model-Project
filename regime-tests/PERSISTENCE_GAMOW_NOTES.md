# 𝒲 as a first-passage law: Gamow expansion (A) and persistence exponent (B) — what holds, what doesn't

_June 2026. Executes the reframe of `W_IDENTITY_NOVEL_ANGLES.md` (𝒲 = Gamow survival with a Kramers right tail
and a persistence left tail). Outcome: **the first-passage decomposition is validated** (right = Kramers, left =
persistence — both named), but **the Gamow-survival *bulk* frame is falsified** and the exponents are
**β-dependent** (no single universal number). Figure `coupled-atlas/figures/persistence_decomposition.png`;
scripts `persistence_check.py`, `standardized_tails.py`, `persistence_figure.py`. Tags
**[DERIVED]/[NUMERIC]/[FALSIFIED]/[CITED]**; ⚑ load-bearing._

---

## 0. Outcome

| piece | result | status |
|---|---|---|
| **Reframe** (𝒲 = first-passage, not spectral) | confirmed: two clean mechanisms, right vs left | **VALIDATED** |
| **C — right tail = Kramers** | heavier, low exponent (~1.2 at β=2, exponential-like) — the early-escape barrier LDP | **VALIDATED (recognition)** |
| **B — left tail = persistence** | *is* a persistence probability $P(\text{no zero of }u)$ — **named**; exponent ~2.0 at β=2 (lighter) | **IDENTIFIED ✓ (but β-dependent)** |
| **A — Gamow survival produces 𝒲** | dominant-resonance bulk = Rayleigh, **skew $-0.63$** vs 𝒲 $+0.60$; left rate η-dep, $\neq1.64$ | **FALSIFIED ⚑** |
| net | 𝒲 is a **β-family** of first-passage laws (Kramers right ⊕ persistence left); **no single universal closed form** | **DERIVED** |

**Verdict.** The note's *direction* is right — 𝒲 lives in first-passage/persistence/Kramers language, not RMT.
But neither A nor B delivers a *universal number*: the Gamow expansion fails to produce the bulk (wrong skew
sign), and the persistence exponent is **β-dependent** (part of the 𝒲_β family), not a single Bray–Majumdar–Schehr
constant. Honest: a productive reframe with named mechanisms, not a closed form.

---

## 1. A — Gamow / resonance survival expansion: falsified as the producer of 𝒲 ⚑

The corrected map (𝒲 = survival $\sum_n c_n e^{-i\lambda_n\Theta}$) was the natural candidate. Two decisive,
cheap checks (`persistence_check.py`) **falsify** it:

- **Bulk skew, wrong sign.** The dominant resonance $\lambda_0=0.86-0.82i$ gives survival
  $S(\Theta)\sim e^{-2|\mathrm{Im}\lambda_0|\Theta}=e^{-1.64\Theta}$; with phase-depth $\Theta=Y^2/2$ this is
  $f(Y^\star)\propto|Y^\star|e^{-0.82Y^{\star2}}$ — a **Rayleigh** in $|Y^\star|$, **skew $-0.631$**. 𝒲 has skew
  **$+0.60$** — *opposite sign*. The Gamow survival cannot produce the bulk. **[FALSIFIED, analytic ⚑]**
- **Left-tail rate, η-dependent and $\neq1.64$.** The measured survival rate in $\Theta=Y^2/2$ is
  **$4.72$ (η=1.0), $6.84$ (η=0.7)** — i.e. $\approx 4.75/\eta$, **not** the η-independent Gamow $1.64$. **[FALSIFIED, NUMERIC ⚑]**

**Diagnosis.** The resonances $\lambda_n$ are a **deterministic** (η-independent) spectral property of the
operator (M1). 𝒲 is the **noise-induced** (η/β-dependent) escape-location fluctuation. The Gamow expansion
describes the *operator's metastable decay*, not the *noisy fluctuation law* — they are different objects. So the
resonances **characterize the operator but do not produce 𝒲**. **[DERIVED.]**

## 2. B — the left tail *is* a persistence probability (named), but β-dependent

The escape is the first zero of the stochastic Weber field $u$ ($p=-u'/u$ explodes at a node), so the
**late-escape left tail is literally $P(\text{no sign change of }u\text{ up to depth }\Theta)$ — a persistence
probability** in the precise Bray–Majumdar–Schehr sense. This **identification holds** and is the genuine
reframe: the "anomalous, no-closed-form" left tail is the **persistence of the stochastic Weber field**, a named
object. **[IDENTIFIED ✓]**

**But the exponent is β-dependent.** Standardized left-tail exponent (`standardized_tails.py`,
$-\log P\sim|s|^{a}$): $a_L=1.99,\,1.85,\,1.68$ at $\beta=2,4,8$. It **drifts with β** (≈Gaussian, $a_L\!\approx\!2$,
at β=2; lighter than that for larger β). So it is **not a single universal persistence exponent** — it is a
**β-family** of persistence exponents. **Reason [DERIVED]:** the stochastic Weber field is **non-stationary**
(instantaneous frequency $|Y|$ grows) *and* noise-scaled by β; Bray–Majumdar–Schehr universal exponents are for
**stationary** processes, so ours has no single universal value. The persistence *framework* is correct; the
*number* is β-dependent. **[IDENTIFIED + honest caveat.]**

## 3. C — right tail = Kramers (validated as a mechanism)

The early-escape right tail is the heavier, **low-exponent** side ($a_R=1.20,\,0.87,\,0.80$ at $\beta=2,4,8$ —
exponential-like, the barrier-crossing LDP). It is the **+skew driver** (right heavier ⇒ skew $>0$), and a
**different mechanism** from the left (so a single closed form was never going to cover both tails — retroactively
explaining the retracted "clean collapse"). At β=2: right exp $\approx1.2$, left exp $\approx2.0$, skew $+0.61$
(matches FP-𝒲). **[VALIDATED.]**

## 4. The honest synthesis — 𝒲 is a β-family of first-passage laws

> **𝒲_β is a first-passage law: a Kramers (barrier-crossing) right tail ⊕ a persistence (no-zero) left tail,
> with a β-dependent crossover bulk. No single universal closed form — it is a genuine one-parameter (β) family
> of first-passage laws.**

This is the corrected, validated picture. It **confirms** the note's reframe (first-passage, not spectral; right
Kramers + left persistence) and **names** both tails, while being honest that (i) the Gamow-survival *bulk* frame
is falsified, and (ii) the persistence exponent is β-dependent, not universal. The mechanisms are named; the law
is a β-family, not a closed form.

## 5. Net + tracker

**Validated:** the first-passage decomposition (right Kramers, left persistence) and the **persistence-exponent
identification** of the left tail (a named object: persistence of the stochastic Weber field). **Falsified:** the
Gamow-survival expansion as the producer of 𝒲 (wrong bulk skew; η-dependent left rate). **Diagnosis:** 𝒲 is a
β-family; the resonances are operator data, not the noisy law; the persistence exponent is β-dependent
(non-stationary process).

**Tracker: B stays 96%.** All three pieces of 𝒲 are now **named** — right tail Kramers (proved earlier), left
tail persistence (identified here), bulk β-family crossover — which fully articulates B's *honest ceiling*: 𝒲_β
is a β-family of first-passage laws with no single closed-form, a settled characterization rather than a missing
calculation. Remaining frontier (unchanged): a β-dependent persistence-exponent computation (non-stationary
Gaussian-process theory) would pin the left tail quantitatively; the right tail's Kramers constant is the proved
no-early-escape rate.
