# Two-time folded medium: what physical drift gives Airy₂ (capstone)

*The final step: a two-time disordered folded medium where the slow drift is a
**physical slow drive**, not a hand-imposed OU. Question: does any slow drive give the
Airy₂ process, or specifically a stochastic one? Script: `twotime_folded_pde.py` ·
figure: `figures/twotime_folded_pde.png`.*

## The discriminator is the increment-variance EXPONENT, not the marginal

The leading escape edge `ξ(τ)=N^{1/6}(λ_top(C(τ))-2√N)` under two physical slow drives of
the disordered coupling `C(τ)`:

| slow drive | increment-variance small-lag exponent | plateau | verdict |
|---|---|---|---|
| **stochastic coupling drift** (noisy plasticity / parameter noise) | **0.77 ≈ 1** (Brownian) | **1.67 ≈ 2·Var(TW₂)=1.63** | **Airy₂** ✓ |
| **deterministic slow drive** (smooth rotation of the disorder) | **1.96 ≈ 2** (ballistic) | turns over | **NOT Airy₂** ✗ |
| detuning-only / uncoupled | — | — | Gumbel marginal (no TW) ✗ |

**Both** coupling drivers share the *same instantaneous marginal* `=TW₂` (var 0.81,
skew +0.17 from independent edge draws) — because `cos(ωτ)A+sin(ωτ)B` of independent GUE
is GUE at every `τ`. The **marginal cannot tell them apart.** The increment-variance
exponent can: Brownian (slope 1) vs ballistic (slope 2).

## What this pins down

Airy₂ is the edge of Dyson **Brownian** motion — the *Brownian* is load-bearing. A
deterministic slow drive, even slow and even over a fully disordered medium with the
exact TW₂ one-point law, produces **smooth (ballistic) spectral flow** — a *different
process* with the *same marginal*. So the hand-imposed OU of the previous note was
secretly supplying exactly one essential ingredient: **stochasticity of the slow drift.**

Physical reading: a two-time folded medium realises the Airy₂ process **iff** the slow
drift of its disordered couplings is *stochastic* (noisy plasticity / slow parameter
noise). Deterministic drives → ballistic, non-KPZ. Detuning-only → no edge. This is the
sharp, falsifiable condition the whole arc was after.

## The full ladder, now complete

1. **TW marginal** — one escape (paper). ✓
2. **Airy point process** — sorted spectrum of one swept disordered network (rung 2). ✓
3. **Airy₂ process** — leading escape under a *stochastic* slow drift of a *disordered
   long-range* coupling (rungs: disorder ⇒ RMT edge; stochastic drift ⇒ Brownian ⇒ KPZ
   class). ✓ Realised here; exponent 1 + plateau 1.63.

Failure modes mapped along the way: nearest-neighbour → Anderson (no edge); uniform
mean-field → BBP outlier; deterministic drive → ballistic (wrong process); detuning-only
→ Gumbel. The needle: **disordered + long-range + stochastically drifting.**

## Caveats / honest scope
- Leading escape `= -2.338 - λ_min(C)` via the validated nonlinear-field correspondence
  (corr 0.95, `FOLDED_PDE_NETWORK_SWEEP.md`); the "medium" is the disordered transverse
  operator (the faithful canard linearisation), not yet a real-space reaction–diffusion
  PDE. Building that PDE is engineering; this note identifies the **mechanism it must
  contain** (stochastically-drifting long-range disorder).
- Single-trajectory marginal skews are unreliable (few independent samples); the clean
  marginal uses independent edge draws. Finite-N throughout (N=96).
- Real (GOE) couplings → the Airy₁/TW₁ sibling.
