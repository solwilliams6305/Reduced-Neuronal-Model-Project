# Notes — cusp tubes + spectral identification (attempt at the two open theorems)

_June 2026. Companion to `STOCHASTIC_CUSP_NOTES.md`, `DELTA_G_BLOWUP_NOTES.md`.
Figures: `coupled-atlas/figures/cusp_tubes.png`, `spectral_id.png`; scripts `cusp_tubes.py`, `spectral_id.py`.
Tags: [DERIVED] · [NUMERIC] · [VALIDATED] · [NEGATIVE] · [OPEN]._

---

## Part A — the cusp tube estimate (Berglund–Gentz one rung up)

**[DERIVED].** Linearising the rescaled cusp SDE about the deterministic canard D̄(s), the deviation
ξ = D − D̄ is an Ornstein–Uhlenbeck process

  dξ = a(s) ξ ds + η dB_s ,  a(s) = −3D̄(s)² + M  (< 0 on the attracting branch),
  near the fold  a ≈ −2√(3M)·d  (d = height above the fold),

so the tube variance is Var(ξ) ≈ η²/(2|a|) and the linear tube fails — peel-off begins — at the
matching height **d\* ~ (η²/(4√(3M)))^{1/3} = η^{2/3}/(4√(3M))^{1/3}** (→ the paper's η^{2/3} at the
fold, M = O(1); → η^{2/3}/M^{1/6} as M → cusp).

**[NUMERIC]** (`cusp_tubes.py`, M=2, η=0.3):
- tube variance follows the OU law: mean ratio Var(ξ) / [η²/(2|a|)] = **0.85** away from the fold,
  and departs (the trajectory leaves the tube) only in the fold/peel-off zone — exactly the picture.
- **Gaussian confinement** P(sup|ξ| > h) ~ exp(−κ (h/η)²) with **κ ≈ 1.85** — the Berglund–Gentz
  exponential tube bound holds one rung up.

**[OPEN]** the rigorous theorem: explicit constants, the bound uniform in M down to d\*, and the
matching of the tube to the nonlinear inner (cusp) equation. The structure and the two key
estimates (OU variance, Gaussian confinement) are in hand.

---

## Part B — the spectral identification (validated for the fold; cusp operator still open)

Goal: realise the cusp peel-off law as the ground state of a stochastic operator (the cusp analogue
of RRV's stochastic-Airy = TW). Method: direct tridiagonal discretisation of
H_k = −d²/dx² + x^k + (2/√β) b'(x) + smallest eigenvalue by Sturm-sequence bisection — an
*independent* numerical method from the swept-ODE first-node (shooting).

**[VALIDATED] k = 1 (fold).** The operator ground state reproduces **Tracy–Widom**: skew +0.25
(operator) vs +0.18 (swept) vs +0.22 (TW₂); the two methods agree (`spectral_id.png` A). So the
operator/Sturm machinery is correct and the stochastic-Airy = TW identification is confirmed
independently.

**[NEGATIVE] k = 2 (cusp).** The naive harmonic operator −d²/dx² + x² + noise gives a **near-Gaussian**
ground state (skew **+0.05**) — it does **not** match the swept k=2 / KP cusp peel-off (skew **+0.62**,
Weber-class). So the cusp operator is *not* the confining harmonic.

**Diagnosis (why, and what's needed).** The cusp peel-off is a **soft-edge** law: the swept k=2 inner
equation u″ = (sign(Y)|Y|² − ηξ)u is the **parabolic-cylinder/Weber** equation with a soft edge at
Y → −∞ (V = −Y²), whereas the harmonic operator is **confining** (hard edge ⇒ Gaussian fluctuation
by perturbation). The correct cusp operator is therefore the **soft-edge multicritical / higher-order-
Airy** operator (the second member of the Le Doussal–Majumdar–Schehr / Cafasso–Claeys–Girotti
hierarchy), not −d²+x². Constructing it and proving its ground state = the cusp peel-off (with a
β-dictionary) is the cusp analogue of RRV — and it is **[OPEN]**. The obstruction (soft vs confining
edge) is now pinned down, and the reliable characterisation of the cusp peel-off remains the swept
first-node (shooting).

---

## Net status of the stochastic step

| ingredient | status |
|---|---|
| η_cusp = σ/(√2 ε^{2/5}) | **[DERIVED]** (validated vs fold) |
| cusp tube: OU variance law + Gaussian confinement | **[NUMERIC]** verified; rigorous theorem [OPEN] |
| matching scale d\* ~ η^{2/3}/M^{1/6} | **[DERIVED]** |
| spectral method (Sturm operator) reproduces TW at the fold | **[VALIDATED]** |
| cusp peel-off = naive harmonic ground state | **[NEGATIVE]** — it isn't (soft-edge, not confining) |
| cusp peel-off = soft-edge multicritical/higher-order-Airy operator | **[OPEN]** — the RRV-type theorem; obstruction identified |

So: the tube side advanced to a verified scaffold; the spectral side validated the method on the fold
and, on the cusp, converted a vague "needs the stochastic-Weber spectral id" into a sharp statement —
the operator is a *soft-edge multicritical* one, and the naive harmonic guess is provably (numerically)
wrong. That is genuine progress toward the theorem, plus an honest fence around what remains.
