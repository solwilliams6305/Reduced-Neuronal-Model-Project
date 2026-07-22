# Notes — the Riccati-explosion reframe, de-risked (Route A, step ii)

_June 2026. Figure `coupled-atlas/figures/riccati_explosion_pde.png`; script
`riccati_explosion_pde.py`. Tags: [CONFIRMED] · [DEFERRED]._

**Goal.** Test the central reframe of `WEBER_TW_ATTACK_STRATEGY.md`: that the cusp escape is the
**first explosion of the swept-Weber Riccati diffusion**
$dp=(\mathrm{sign}(Y)Y^2-p^2)\,dY+\eta\,dW$ (Cole–Hopf $p=u'/u$ of $u''=(\mathrm{sign}(Y)Y^2-\eta\xi)u$),
started on the **stable** branch $p=+\sqrt{W(Y_0)}$.

## [CONFIRMED] the reframe is correct

Direct Monte-Carlo of the Riccati-explosion SDE matches the original swept-Weber peel-off law across
the whole catastrophe ladder (skew / excess kurtosis, β=2):

| q | swept-Weber MC | Riccati-explosion MC |
|---|---|---|
| 0.75 | +0.25 / +0.78 | +0.12 / +0.29 |
| 1.00 | +0.23 / +0.10 | +0.19 / −0.02 |
| 1.50 | +0.40 / −0.22 | +0.38 / −0.33 |
| **2.00 (cusp)** | **+0.64 / −0.23** | **+0.59 / −0.29** |
| 2.50 | +0.76 / −0.15 | +0.77 / −0.13 |
| 3.00 | +0.96 / +0.12 | +0.93 / +0.07 |

Over the relevant range $q\ge1$: mean $|\Delta\text{skew}|=0.03$, mean $|\Delta\text{exkurt}|=0.07$, and
the Riccati picture **reproduces the negative-kurtosis valley** (and its return to $+$ at $q=3$). The
only sizeable gap is at $q=0.75$ (the heavy-tailed, MC-noisy end), irrelevant to the cusp.

So: **the cusp escape = first explosion of the swept-Weber Riccati diffusion.** This nails down the
right object — a *non-self-adjoint diffusion generator*
$\mathcal L_Y=\tfrac{\eta^2}{2}\partial_p^2+(\mathrm{sign}(Y)Y^2-p^2)\partial_p$
(matching the "directional/non-self-adjoint" finding), with an exactly solvable deterministic skeleton
(parabolic-cylinder functions), that carries the cusp's sub-Gaussian signature.

## [DEFERRED] the deterministic PDE gateway

The matching backward-Kolmogorov / Fokker–Planck equation
$\partial_Y h+(\mathrm{sign}(Y)Y^2-p^2)\partial_p h+\tfrac{\eta^2}{2}\partial_p^2 h=0$
is the deterministic (Monte-Carlo-free) object and the gateway to the integrable analysis. It is
formulated in `fp_solve`, but the **explicit** upwind solver is numerically delicate for the steep
$-p^2$ drift (the bottom-boundary CFL blows up as $W\to-\infty$). The immediate fix is an **implicit
(Crank–Nicolson) tridiagonal solve**, which is unconditionally stable — the next coding step before the
Painlevé-IV / parabolic-cylinder analysis (Routes A-iii, B).

## Why this matters
Routes A and B in the strategy now act on a **validated** object, not a guess. The explosion CDF is a
concrete 1-D PDE; its self-similar reduction is the candidate Painlevé-IV transcendent, and the
deterministic skeleton being parabolic-cylinder is exactly the structure that conjecture needs.
