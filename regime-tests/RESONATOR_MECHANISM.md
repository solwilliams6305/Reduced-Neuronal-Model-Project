# The Resonator Regime: a mechanistic picture of phase-gated stochastic commitment

**Model.** Stochastic FitzHugh–Nagumo with noise in the voltage channel only:

```
dv = ( v - v^3/3 - w + I ) dt + sigma * sqrt(dt) * xi
dw = eps ( v + a - b w ) dt
```

Parameters held fixed at `a = 0.7`, `b = 0.8`, `eps = 0.08`. The deterministic Hopf
bifurcation (lower branch) sits at `I_H = 0.331281`. The resonator regime is the
band just below it, where the resting fixed point is a **stable spiral (focus)**:
small perturbations ring and decay rather than relaxing monotonically. The
canonical operating point used throughout is `I = 0.30` (a deficit of `delta =
0.031` below the Hopf), with fixed point `(v*, w*) = (-0.993, -0.367)`, decay rate
`alpha = -0.0253`, spiral frequency `omega = 0.280`, and damping-to-frequency
ratio `kappa = |alpha|/omega = 0.090`.

In this regime the cell does not fire on its own. Noise jiggles the trajectory
around the focus; it rings, decays back, gets re-kicked, and occasionally — on one
lucky orbit — is thrown across the separatrix into a full spike. The question this
work answers is: **when the cell finally commits to a spike, is there a preferred
direction it leaves in, and what sets that direction?**

## The phenomenon

Working in the local spiral coordinates (the eigenbasis of the focus, where the
trajectory's angular position is a well-defined "phase"), the escapes are sharply
phase-concentrated. But concentration measured *at the spiking threshold itself* is
a tautology — escape is *defined* as crossing `v ≈ 1`, so of course all escapes
share that surface (the measured resultant `R ≈ 0.997` there is geometric, not
informative).

The real signal is upstream. Measured on a **mid shell** deep inside the spiral
(at `0.25` of the threshold radius), the phase of the *final outbound commitment* —
the last time a trajectory crosses that shell before spiking — concentrates at an
angle that is offset from the deterministic "steepest-outward" direction (the
geometric-danger angle) by roughly **half a radian**. This offset is signed (the
commitment sits to one angular side of the geometric direction) and, crucially, it
depends on the noise amplitude. That noise-dependence is the first sign that the
offset is a *selected* feature, not a property of the deterministic flow geometry.

## The chain of evidence

Four tests, each ruling something in or out.

**1. The gate is real, not geometric** (`test_commitment_gate.py`). Comparing the
escape-commitment phase to the deterministic steepest-outward direction at a range
of radii shows they coincide at the threshold (tautological) but separate by ~0.5
rad on the mid shell, and the separation grows as noise shrinks. A purely geometric
gate would show a fixed, noise-independent offset. This one moves with sigma, so it
is noise-selected.

**2. It is not linear resonance** (`resonance_scaling.py`). If the offset were set
by the focus ringing at its natural frequency, it should scale with the linear
damping-to-frequency ratio `kappa = |alpha|/omega`. Sweeping operating points that
independently move `alpha` (distance below the Hopf) and `omega` (via `eps`), the
offset stays essentially invariant at ~ -0.55 rad; the correlation with `kappa` is
weak (`-0.40`) and there is no clean scaling. So the gate is **not** the cell's
linear resonance.

**3. It is set by the nonlinear escape path** (`instanton_exit.py`). Freidlin–
Wentzell large-deviation theory says that in the small-noise limit, escapes
overwhelmingly follow a single most-probable path — the instanton — the minimizer
of the action functional. Because the noise is degenerate (voltage only), the
w-channel must follow its noiseless ODE exactly, which we impose with a large
anisotropic metric penalty. Computing this path and reading off the phase at which
it crosses the mid shell gives a *predicted* commitment offset. The leading-order
instanton predicts an offset of the **same sign** as the measurement.

**4. The prediction is robust to winding** (`instanton_exit.py --method gmam`). The
worry — a good one — was that typical resonator escapes orbit the focus several
times before the knockout, so maybe the relevant instanton should wind, and a
single fixed-duration path was an artifact. The geometric (arc-length) minimum-
action method removes the path-duration freedom entirely. Started from a straight
line, a 2-turn spiral, and a 4-turn spiral, it converges to the **same** committed
crossing phase to within `0.008` rad. The committed offset is `-1.06` rad,
*independent of how many times the path winds*. Winding near the focus is cost-free
(it just follows the drift); the final outbound leg is geometrically pinned.

**5. Finite-noise measurements converge to it** (`sigma_lead_scaling.py`). The
gMAM instanton is the `sigma → 0` ceiling; the ensemble measurements are at finite,
near-threshold noise. Lowering sigma from `0.06` to `0.02`, the measured offset
climbs monotonically and with steepening slope: `-0.39 → -0.43 → -0.49 → -0.55 →
-0.56 → -0.60 → -0.65`, heading straight for the `-1.06` ceiling. Escapes vanish
below `sigma ≈ 0.02` (only 25% of trajectories spike even at `T = 600`), which is
precisely why the instanton calculation is needed in that limit. The earlier
apparent "discrepancy" between `-1.06` (theory) and `-0.55` (measurement) was never
a discrepancy — it is the finite-noise correction, and it behaves exactly as that
interpretation demands.

## The mechanistic picture

Putting it together, the resonator regime commits to a spike like this.

The resting state is a stable focus, so the cell's response to a kick is to ring:
spiral inward at frequency `omega`, decaying at rate `alpha`. Under continuous
noise the trajectory therefore lives as a rotating, breathing cloud around the
focus — winding many times, each orbit a fresh attempt at escape, almost always
failing and decaying back. The *waiting* is governed by how rare a large enough
fluctuation is; the *winding* is just the relaxation dynamics playing out between
attempts and carries no directional information about the eventual exit.

When an escape does happen, it is not a kick in the geometrically easiest
direction, nor a resonant build-up at the natural frequency. It is a single
large-deviation event that follows the **minimum-action escape path** — the
cheapest route, in the sense of the noise, from the focus across the separatrix to
the spiking branch. That path is a property of the *nonlinear* vector field and the
*degenerate* (voltage-only) noise geometry, not of the linearized spiral. Its final
approach to the separatrix crosses the mid shell at a definite angle, and that angle
— offset ~1 rad to one side of the deterministic steepest-outward direction in the
`sigma → 0` limit — is the commitment gate.

At realistic, finite noise the cell does not have to pay the full action, so the
realized escapes are pulled partway back toward the geometric direction; the gate
sits at a smaller offset (~0.4–0.65 rad over the measured range) and slides toward
the instanton angle as noise is reduced. The gate is sharp (concentration `R ≈
0.99`) because the action landscape is steep around its minimizer: paths even
slightly off the optimal exit angle are exponentially less likely.

So the one-line mechanism: **in the resonator regime the spike-commitment direction
is set by the nonlinear, degenerate-noise large-deviation escape path (the
instanton), robustly independent of the orbital winding that precedes it, with the
finite-noise gate sitting at a reduced offset that converges to the instanton's
prediction as noise shrinks.** The focus geometry (steepest descent) and the linear
resonance (`kappa`) are both *not* the controlling quantities — they are ruled out
directly by tests 1 and 2.

## Numbers, at a glance

| Quantity | Value |
|---|---|
| Canonical point | `I = 0.30`, `a = 0.7`, `b = 0.8`, `eps = 0.08` |
| Hopf current | `I_H = 0.331281` (deficit `delta = 0.031`) |
| Focus | `alpha = -0.0253`, `omega = 0.280`, `kappa = 0.090` |
| Threshold-surface concentration | `R ≈ 0.997` (tautological) |
| Mid-shell gate offset (finite sigma) | ~ `-0.4` to `-0.65` rad, sharp (`R ≈ 0.99`) |
| Offset vs `kappa` | invariant; corr `-0.40` → not linear resonance |
| gMAM instanton offset (`sigma → 0`) | `-1.06` rad, init-independent (spread `0.008`) |
| Measured offset vs sigma | `-0.39 (0.06) → -0.65 (0.02)`, monotone toward `-1.06` |
| Arrhenius fit at canonical point | `B = 2.5e-3`, `A = 0.063`, `R^2 = 0.984` |
| `B(delta)` low-delta slope | `0.32 ± 0.09` (bootstrap over seeds) |
| `B(delta)` high-delta slope | `0.63 ± 0.04`; crossover positive in `99.8%` of resamples |
| Single fitted exponent | protocol-dependent (`~0.52` vs `~0.71`) → not a constant |

## The barrier vs. Hopf deficit: a crossover, not a universal exponent

The escape rate follows a clean Arrhenius law `k(sigma) ~ A exp(-B / sigma^2)`
(fit `R^2 = 0.984` at the canonical point), so the barrier `B = ΔV` is the single
load-bearing quantity controlling how rare escapes are. The natural next question
is how `B` scales with the Hopf deficit `delta = I_H - I` as the cell is pushed
deeper below the bifurcation.

A multi-seed measurement (`multiseed_delta.py`: deltas `0.01–0.08`, seeds `1,2,3`,
sigma-grid `0.04–0.08`) settles what is and is not real here:

- **Real (locked in):** `B` rises with `delta`, and the *local* slope steepens as
  you move away from the Hopf. Bootstrapping over seeds, the high-`delta` slope
  exceeds the low-`delta` slope (`0.63 ± 0.04` vs `0.32 ± 0.09`; difference
  positive in `99.8%` of resamples). So there is a genuine **crossover**, not a
  single clean power law, and it survives Monte-Carlo error bars.
- **Not real / not a constant:** the *value* of any single fitted exponent. It is
  protocol-dependent — `~0.52` on this sigma-grid versus `~0.71` on a lower one,
  each with tight error bars. The shift of `~0.2` with the fit window is the
  signature of a pre-asymptotic regime (shallow near-Hopf barrier; the frozen-w
  well and saddle nearly merged at the fold, so the 1D barrier `≈ 0`), not of a
  universal scaling law waiting to be measured more precisely.

The defensible claim is therefore the qualitative one: **B(delta) has a positive,
delta-dependent effective slope with a robust low-to-high crossover, and there is
no universal `B ~ delta^α`.** Pinning a single `α` numerically is not meaningful
in this regime; giving the exponent physical meaning would require an analytic
near-Hopf prediction (normal-form / degenerate-noise WKB), not more seeds.

## Chapter synthesis: one regime, two questions, one answer

The resonator regime poses two questions, and they turn out to have the same answer
read from two angles of the same object — the nonlinear, degenerate-noise escape
geometry.

*Where does the spike leave?* Not in the geometrically easiest direction and not at
the focus's resonant frequency, but along the large-deviation instanton — the
minimum-action escape path — whose committed crossing angle is pinned independent of
the orbital winding that precedes it (`-1.06` rad as `sigma -> 0`, finite-noise gate
sliding toward it).

*How rare is the spike?* Arrhenius, `k ~ A exp(-B/sigma^2)`, with the barrier
`B = ΔV` the whole `sigma`-dependence. And `ΔV` is `1/2` the geometric action of that
same instanton — the gate's *direction* and the escape's *rate* are the two readouts
of one quasipotential.

The barrier's parameter-dependence then closes the loop with the linear theory that
tests 1–2 ruled out. The Hopf is subcritical (`a_3 = +0.268`), so the unstable limit
cycle is the threshold and the normal form gives `B -> alpha^2/(a_3 |q_v|^2) ~ delta^2`
in the `(delta, sigma) -> 0` corner. The measured `~delta^0.5` is the pre-asymptotic,
floor-dominated barrier *below* that corner, and the empirical crossover is the data
bending up toward `delta^2` as `delta` grows. So the resonator regime is, end to end:
**a focus whose linear ringing carries no directional information, in which both the
direction and the rarity of commitment are set by a single nonlinear large-deviation
object, whose near-Hopf scaling is `delta^2` but is only observable pre-asymptotically.**
The linear focus geometry (`kappa`) and the 1D frozen-w barrier are both explicitly
not the controlling quantities. That is the complete picture for this regime.

## Caveats and the one open thread

The instanton is leading-order (`sigma → 0`); the w-constraint is enforced by a
finite metric penalty, not exactly; and the action values across the duration-based
solver were not a clean single minimization (which is why the geometric method
was the one to trust). The `sigma → 0` extrapolation of the measured curve cannot
be taken all the way numerically because escapes become too rare — the trend and
its acceleration are unambiguous, but the precise limiting value rests on the
instanton, not on extrapolated data. A naive quadratic extrapolation lands near
`-0.87`; the true curve steepens faster than quadratic at the small-sigma end, so
`-1.06` is consistent but should be read as "the trend points there," not "the data
reach there." The clean way to tighten the last step would be a sub-exponential
prefactor (Eyring–Kramers-type) correction to the instanton, predicting the *rate*
at which the finite-sigma offset relaxes toward the ceiling.

## Scripts

The argument is reproducible from four scripts in `regime-tests/`:

- `test_commitment_gate.py` — gate is real and noise-selected, not geometric.
- `resonance_scaling.py` — gate does not scale with `kappa`; not linear resonance.
- `instanton_exit.py` (`--method gmam`) — gate is the large-deviation instanton,
  robust to winding (`-1.06` rad).
- `sigma_lead_scaling.py` — finite-noise gate converges toward the instanton.
- `arrhenius_escape.py` — escape rate is Arrhenius; measures the barrier `B = ΔV`.
- `multiseed_delta.py` — `B(delta)` crossover survives error bars; no universal exponent.
- `normal_form_barrier.py` — subcritical-Hopf normal form gives `B ~ delta^2`; the
  measured shallow exponent is the pre-asymptotic floor (see `BARRIER_NORMAL_FORM.md`).
- `multiseed_delta.py` — `B(delta)` crossover survives error bars; exponent is not
  a universal constant.
