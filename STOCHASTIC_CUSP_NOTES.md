# Notes — the stochastic cusp step (toward the TW→Weber-edge theorem)

_June 2026. Companion to `DELTA_G_BLOWUP_NOTES.md`, `NOISY_CUSP_CROSSOVER_NOTES.md`.
Figure: `coupled-atlas/figures/stochastic_cusp.png`; script `coupled-atlas/stochastic_cusp.py`.
Tags: [DERIVED] · [NUMERIC] · [HEURISTIC] · [OPEN]._

This is the noise layer on top of the deterministic Δ(g)/blow-up scaffold. Three pieces.

---

## 1. Noise covariance in the cusp chart — [DERIVED]

Additive noise enters the antisymmetric mode as σ_δ = σ/√2 (from independent σ dW on each v_i).
Push it through the A₃ cusp blow-up (δ,μ,δw) = (r D, r² M, r³ W), inner time s = r² t, with
r ~ ε^{1/5} fixed by balancing the slow load drift d(δw)/dt = O(ε) against the inner cusp flow.
The Itô-on-the-diffeomorphism calculation (no Stratonovich correction, additive noise) gives

  dD = (−D³ + M D − W) ds + **η_cusp** dB_s ,   **η_cusp = σ_δ / r² = σ / (√2 · ε^{2/5}).**

**Validation:** the identical method applied to the *fold* chart (r=ε₂^{1/3}, s=ε₂^{1/3}t) reproduces
the paper's η = σ/√ε₂ exactly — so the η_cusp exponent 2/5 (vs the fold's 1/2) is trustworthy.

---

## 2. Berglund–Gentz tube, one rung up — [NUMERIC]

Simulating the rescaled cusp SDE, the noisy trajectory stays in an **O(η)-tube around the
deterministic canard** (the upper branch) until the peel-off: at M=2, η=0.6 the pre-escape tube
half-width is ≈ 0.21 (≈ 0.36 η), and it scales with η. So the Berglund–Gentz confinement picture
transfers from the fold/folded-node to the cusp — the prerequisite for sample-path control.
(Rigorous tube *estimates* in the cusp chart remain [OPEN]; this is numerical confirmation.)

---

## 3. The stochastic uniform connection (TW ↔ Weber) — [NUMERIC, with caveats]

Sweeping the chart parameter M = μ/ε^{2/5} = −2g/ε^{2/5} and measuring the inner exit measure
(peel-off load at escape): the law's shape **shifts toward the Weber/cusp class as the two folds
coalesce (M → cusp)**. Skew rises from ≈ −0.4 at resolved M (M ≳ 2.6) to ≈ +0.13 at M = 0.5
(cusp) — the inner stochastic Airy↔Weber connection is visible (panel B/C).

**Honest caveats:**
- The skew(M) curve is **non-monotone** on the resolved side: the local fold curvature √(3M)
  sets an M-dependent β, so the "TW value" itself drifts with M. The robust, clean signal is the
  rise toward the cusp, not a single TW plateau.
- **Observable reconciliation (important).** This inner SDE is the *antisymmetric* mode's OWN
  escape, which is Weber near *its* cusp (small M, i.e. g near 0). The full-model peel-off I measured
  earlier (`kp_cusp_noise.py`: w₁ at the *symmetric* fold release) is a DIFFERENT observable — the
  symmetric escape **modulated by** the antisymmetric funnel — whose Weber signature is strongest
  when the funnel is *active and strong*, i.e. at large μ (g very negative). The two are consistent
  once read as different observables; the inner antisym cusp escape and the modulated symmetric
  escape have opposite μ-dependence by construction. Making this reconciliation rigorous is part of
  the theorem.

---

## What the theorem still needs — [OPEN]

1. **Rigorous tube estimates** in the cusp chart (the Berglund–Gentz tubes one rung up): prove the
   O(η_cusp) confinement and the exit-time/exit-measure control that §2 shows numerically.
2. **The uniform connection, proved:** the inner exit measure is the ground state of the stochastic
   **Weber/parabolic-cylinder** operator for M = O(1) and of the stochastic **Airy** operator (→ TW)
   for M ≫ 1, *uniformly in M* — the noise analogue of Olver's coalescing-turning-points asymptotics.
   §3 is the numerical target; the proof needs the RRV-type spectral identification for the Weber
   operator (the cusp analogue of the paper's stochastic-Airy = TW step).
3. **Observable map:** relate the inner antisymmetric exit measure to the physical symmetric-fold
   release (the full-model peel-off), i.e. how the funnel imprints on the escape — closing the loop
   to the measured KP-cusp law (`weber_vs_pearcey.py`, Part 1).

---

## Status ledger (stochastic step)

| statement | status |
|---|---|
| η_cusp = σ/(√2 ε^{2/5}); method validated on the fold (→ σ/√ε₂) | **[DERIVED]** |
| O(η) Berglund–Gentz tube around the canard in the cusp chart | **[NUMERIC]** (rigorous estimates [OPEN]) |
| inner exit measure shifts TW → Weber as folds coalesce (M→cusp) | **[NUMERIC]** (non-monotone; clean signal = rise toward cusp) |
| inner exit measure = stochastic-Weber / stochastic-Airy ground state, uniform in M | **[OPEN]** — the theorem |
| reconciliation of inner antisym escape vs full-model symmetric release | **[HEURISTIC]** now; [OPEN] rigorously |

Net: the stochastic step's *scaffolding* is in place — the noise scale is derived (and validated),
the tube confinement holds one rung up, and the inner exit measure numerically shows the TW→Weber
shift. What remains is genuinely the hard analysis: the cusp tube estimates and the stochastic-Weber
spectral identification (the cusp analogue of RRV's stochastic-Airy = TW), plus the observable map.
