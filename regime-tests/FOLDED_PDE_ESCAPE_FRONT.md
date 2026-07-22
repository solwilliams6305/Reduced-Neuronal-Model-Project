# Escape front in space-time — dynamical SDE field (follow-up to FOLDED_PDE_AIRY_PROCESS)

*The operator note settled the mechanism statically. This is the faithful dynamical
test: integrate the actual field of noisy folded cycles through the fold and look at
the escape front in real space-time. Script: `folded_pde_escape_front.py`
(engine reused from `folded_pde_airy_process` ideas) · figure:
`figures/folded_pde_escape_front.png`.*

## Engine + a regime lesson learned the hard way

**[R] Stable canard engine.** Real-space front needs the *canard* (repelling-branch)
escape, not the generic one. Integrating the **linear Cole–Hopf** field downward on the
**recessive branch** (`U_YY=(Y+C-ηΞ)U`, `R=-U'/U`, peel-off = first node of `U`) is the
self-correcting route (the handover's branch lesson). Validated: deterministic first
node `Y=-2.340` (first Airy zero −2.338); uncoupled noisy marginal reproduces
**TW₂** (−1.76/0.81/0.19) and **TW₁** (−1.22/1.56/0.26).

**[R] The trap, confirmed.** A naive *nonlinear* field with a slow *physical* ramp
escapes at the **saddle-node** (`Y≈0`, σ²≈0.006), **not** the canard — the generic
trajectory rides the attracting branch and pops at the fold. The TW edge is an
**inner/blow-up** object living on the repelling branch; it is invisible to a face-value
physical-regime simulation. (This is exactly why the first space-time run gave a marginal
at −0.2, not −2.34.)

## What the space-time front shows

**[N] Local coupling → local front.** Diffusively coupled canard field: the peel-off
front `Y_node(s)` has a **short correlation length** (~1 site in Cole–Hopf coupling; ~6
sites for physical-`R` diffusion) and a marginal that stays TW-shaped but narrows
(mean −2.08, var 0.62, skew +0.24). Uncoupled → iid front (zero correlation). So a
**homogeneous diffusive folded medium produces a smooth *local* escape front** — no
long-range, KPZ-class structure.

**[N] Leading-edge law — the discriminator.** The first site/mode to peel (leading edge
of the front) over many realizations:

| coupling            | leading-edge skew | law |
|---------------------|-------------------|-----|
| uncoupled           | **+1.18**         | Gumbel (extreme-value of iid; no repulsion) |
| diffusive (local)   | +0.77             | still extreme-value-like — local coupling adds no genuine repulsion |
| **disordered long-range** | **−0.35** ( \|0.35\| ) | **Tracy–Widom-class** (TW₁ \|skew\| 0.29; reflected — the lower spectral edge of `C` sets the first peel-off) |

This **dynamically reproduces the static operator result**: only disordered long-range
coupling turns the extreme-value (Gumbel) leading edge into a Tracy–Widom edge — the
prerequisite for the Airy process.

## Verdict (consistent across both notes)

The **literal real-space escape front is local/saddle-node** physics; the **Tracy–Widom /
Airy edge lives in the eigenmode (disordered-coupling) structure**, not in the
real-space front of a homogeneous medium. Two regimes, cleanly separated:

- homogeneous + local (diffusive PDE) → smooth saddle-node front, finite correlation
  length, **no Airy**;
- disordered + long-range → Gumbel→**TW** leading edge, genuine repulsion → the Airy
  route (operator note).

So "escape front in space-time" and "Airy₂ process" are **not the same picture**: the
front you can watch sweep through space is the wrong place to look for TW; the right
place is the spectrum of the disordered transverse operator. The genuine Airy₂ target
remains a **disordered, long-range** folded medium (quenched-random network), per
`FOLDED_PDE_AIRY_PROCESS.md`.

## Caveats / next
- Cole–Hopf coupling ≠ physical-`R` coupling; a fully nonlinear disordered-network field
  (escape front in eigenmode time, swept) is the next, heavier confirmation and where
  the Airy₂ **process** covariance (not just the TW marginal) would be measured.
- Skew estimates at N=40, R=250; the contrast (1.18 vs 0.35) is robust, the exact values
  are finite-size.
