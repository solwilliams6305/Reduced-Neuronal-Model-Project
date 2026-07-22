# The two edges at the collective level — neural-mass theory track (opening probe)

> **CORRECTION (see `MPR_FINITE_SIZE_LANGEVIN.md` §5).** The "macroscopic *phase* edge" claimed below
> was **wrong**. On closer diagnosis the macroscopic bifurcation at this operating point is a **fold of
> limit cycles** (the *amplitude* edge): the period stays finite ($\sim$17) and the amplitude steady
> ($\sim$2.1), then the cycle vanishes abruptly, with a coexisting rest state. The collective
> $\sigma^{2/3}$ test fails (slope $-0.18$, CV $\approx0.78$): the noise-induced bursts are
> **rate-governed bistable hopping**, not the phase-edge first-passage law; the classifier's "SNIC" call
> was a misfire on skewed escape intervals. What *survives* is the system-size result — the macroscopic
> noise is $O(N^{-1/2})$, validated directly ($\mathrm{Var}(r)\propto1/N$, slope $-0.98$). Read this file
> as the (partly mistaken) exploration; the corrected account is `MPR_FINITE_SIZE_LANGEVIN.md`.

*The first probe of the parallel theory track: do the single-neuron two-edge laws govern the COLLECTIVE
dynamics of a next-generation neural-mass model? Code: `neural_mass_edge.py` (+ `neural_mass_figure.py`).
Figure: `figures/neural_mass_edge.png`. Tags: **[R]** proved / **[N]** numerical / **[H]** heuristic.*

## The bridge

The single-neuron edge theory lives in the noisy QIF/SNIC. The **Montbrió–Pazó–Roxin** (2015) exact
mean-field of an all-to-all QIF network is a 2D macroscopic system in the firing rate `r` and mean
potential `v`:
```
  r' = Δ/π + 2 r v ,   v' = v² + η̄ + J r − π² r² .
```
Pure 2D MPR has only up/down bistability. Adding spike-frequency **adaptation** — the standard
next-generation neural-mass route to collective rhythms — makes it a 3D slow–fast system
(`τ_a a' = −a + α r`, with `−a` in the `v` equation) whose fast `(r,v)` subsystem is swept by the slow
adaptation: a **network-level folded cycle**.

## The result

- **[N] A macroscopic phase edge.** Drifting the mean drive `η̄`, the collective-burst period
  **diverges** as `η̄ → η̄* ≈ −2.32` (period 12.4 → 17.2, then the cycle is destroyed) — a *network-level*
  SNIC/homoclinic of the rate dynamics. The whole population has a phase edge, just like the single
  neuron's SNIC.
- **[N] The collective bursts sit at that same edge.** Just below `η̄*`, with noise on the macroscopic
  drive, the **noise-induced inter-burst intervals** are the universal first-passage shape:
  positive-skew, classified **Type-I (SNIC)** by the *same atlas* that classifies single-neuron ISIs,
  read by the *same inversion*. CV falls toward the critical value as `η̄ → η̄*` (1.5 → 0.95) — the
  phase-edge crossover, now of the collective variable.
- **The bridge, stated:** the single neuron and the whole network sit at the **same universality edge**.
  The two-edge dichotomy is not only a single-cell statement; it organises the fluctuation statistics of
  the macroscopic neural-mass dynamics — which the next-gen neural-mass literature models
  deterministically. That is where these results are novel to that community.

## Why this is the right wedge into the field

Next-generation neural-mass models (MPR; Coombes–Byrne) are a live, growing modelling program, and they
are largely **deterministic** or add noise ad hoc. A principled account of **finite-size / dynamic-noise
fluctuations near their bifurcations** is a genuine gap — and the two-edge theory is exactly a theory of
fluctuation universality at those bifurcations. This slots into models people are already building,
rather than competing with the (crowded) early-warning or seizure literatures.

## Honest status — this is an opening probe, not a theory

- **The noise is phenomenological**, added to the macroscopic `v` (a fluctuating global drive), not
  derived. The rigorous object is the **finite-size Langevin correction** to MPR (a `1/√N` system-size
  expansion; the Lorentzian/Ott–Antonsen ansatz that makes MPR exact breaks under dynamic noise). Deriving
  that — and showing its near-bifurcation reduction is the stochastic-Airy / quartic-FPT inner equation —
  is the actual theorem, and it is not done here.
- **Burst extraction is imperfect.** Inter-burst intervals are read from threshold-crossings of `r` with a
  refractory; near the edge the bursts re-cluster and a small short-interval residual survives (a cleaner
  order-parameter / phase definition would remove it). The class call and the bulk shape are robust; the
  precise constants are not claimed.
- **Only the phase edge is shown.** The macroscopic **amplitude edge** (a network Hopf / fold-of-cycles,
  Tracy–Widom collective amplitude statistics) is the natural companion and is untouched.
- **One operating point.** `(J, Δ, τ_a, α) = (15, 1, 15, 5)`; the macroscopic bifurcation structure across
  parameter space is unmapped.

## The program this opens

(1) Derive the finite-size Langevin MPR and its near-bifurcation inner reduction → a *theorem* that the
macroscopic phase fluctuations are the quartic-FPT universal law. (2) The collective amplitude edge at a
network Hopf → collective Tracy–Widom. (3) `σ^{2/3}` vs `σ^{1/2}` scaling of the collective fluctuations
vs system size `N`. Each is a concrete, citable contribution to the next-gen neural-mass community.
