# Canard + noise: conceptual grounding before the instanton

Written after the dynamic-passage experiment showed the measured "exponent" is
protocol-bound (ramp 0.5->0.25 moves it 0.94->1.60; (a,b) moves it 0.78->0.94) and the
ramp->0 extrapolation has no clean limit (sigma_crit -> 0 because the autonomous canard
strip is exponentially thin). That outcome forced the question: do we have the right
theory frame? Partly. The engine is right; the observable and the quoted exponent were not.

## What is solid (keep)

The rigorous backbone for noise in slow-fast FHN is Freidlin-Wentzell large deviations
organized through a **quasipotential**, with the slow variable frozen. Doss-Thieullen make
this explicit for FHN: freeze the slow variable y, the fast subsystem is a 1-D gradient
flow with two wells, and the well depths

    V_pm(y) = -2 * integral_{x*_pm(y)}^{x*_0(y)} (-y + f(u)) du

are the FW transition costs; escape times obey exp(V/eps~) with eps~ = noise^2/timescale.
This is the **same machinery** as our Arrhenius / B(delta) work — the resonator barrier B
is a quasipotential well-depth and exp(-B/sigma^2) is the FW rate. So the core grounding
is correct and already validated. We simply aimed it at the wrong probe (ramp-induced
delay) for the canard question.

## What was shaky (fix)

"Canard" was covering three distinct objects:

1. **Autonomous canard window** — exponentially thin in I, O(exp(-c/eps)). Not power-law.
2. **Slow passage / delayed Hopf** — what we simulated. Delay observable is intrinsically
   protocol-dependent; the ramp experiment proved it.
3. **Noise threshold** separating "noise tracks the canard" from "noise jumps early."

The `sigma_crit ~ eps^{3/4}` was a **project conjecture, not a Berglund-Gentz theorem**.
BG do not quote a bare power law; they organize the canard-noise problem through a scaling
ratio (their c proportional to noise^2 * |log| / timescale) and through **blown-up
(desingularized) coordinates** near the folded singularity, where the natural scales are
eps^{1/3} in space/time and the noise enters as sigma / eps^{1/4}. The "3/4" was the
heuristic sigma/eps^{1/4} reaching an eps^{1/2} tube. It lives inside a scaling region that
ordinary FW theory cannot resolve, because the deterministic drift vanishes at the fold and
the naive quasipotential degenerates there.

## The correct observable for the instanton

The clean quantity that survives sigma, ramp -> 0 is the large-deviation cost to **peel off
the repelling slow manifold** onto the opposing attracting branch: the quasipotential gap
between the attracting and repelling Fenichel manifolds. Compute it in **blow-up coordinates
near the v = -1 fold**, keeping noise **degenerate** (only in v — Berglund-Gentz's natural
fast-variable-noise case, and the case our |q_v|^2 projection already respects). The
eps-dependence of that scaled gap is the exponent, derived rather than fit.

Concretely, the FHN critical manifold is w = v - v^3/3 + I (fast nullcline), folding at
v = +-1; |v|>1 attracting, |v|<1 repelling. The lower Hopf is a **singular Hopf** (Hopf
within O(eps) of the v=-1 fold), i.e. the canard-generating configuration. Near the fold
set v = -1 + eps^{1/3} u, rescale time by eps^{1/3}, and the leading problem is a
constant-drift + cubic fold normal form; the FW action in these scaled variables yields the
sigma-vs-eps threshold directly.

## Recommendation

Do **not** rebuild the instanton solver against the dynamic-passage framing. Build it to
compute the quasipotential gap between attracting/repelling slow manifolds in the fold
blow-up, with degenerate noise. The exponent is then a property of the scaled barrier, not
of any transport protocol — and it can be checked against BG's scaling-region results
rather than against a power-law fit.

## Sources

- Berglund, N. & Gentz, B., *Noise-Induced Phenomena in Slow-Fast Dynamical Systems: A
  Sample-Paths Approach*, Springer (2006) — quasipotential / scaling-region framework.
- Doss, C. & Thieullen, M., *Oscillations and Random Perturbations of a FitzHugh-Nagumo
  System*, arXiv:0906.2671 — frozen-slow quasipotential V_pm(y), LDP exit times.
- Berglund, Gentz, Kuehn, "Hunting French ducks in a noisy environment" (canards + noise).
- Muratov, Vanden-Eijnden, E., *Self-induced stochastic resonance in excitable systems*,
  Physica D (2005) — noise/timescale matching in fast-variable-noise FHN.
