# α=1 (rotating phase) — attempting the four approaches, and what the numerics decided

*Goal: try the creative routes to the coupled rotating-phase folded-cycle escape
(A two-zone decoherence, B Riccati→Schrödinger Lyapunov, C Fourier-in-phase
averaging, + Kapitza wildcard) and let a discriminating simulation pick the winner.
Scripts: `outputs/alpha1_test{,2,3}.py`.*

---

## The discriminating prediction

All four routes agree on the *machinery* but split on one sharp, testable question:
**when the phase rotates fast through the passage, is the escape threshold set by the
weakest phase (min G) or by an average of G(θ) = √(a c / b)(θ)?**

- **A/B/C as I first framed them ⇒ min-barrier.** Reasoning: a fast-rotating cycle
  revisits its weakest phase every turn, and a rare escape is a Laplace integral over
  the visited phases ⇒ dominated by `min_θ G`. (B: escape = band edge of the swept
  potential; A: decohered core sees the easiest phase.)
- **Kapitza ⇒ average-barrier.** High-frequency averaging replaces the modulated
  barrier by its mean ⇒ `⟨G⟩`.

So I built the α=1 model SDE
`dR=(R²−a(θ)Y)dT+η dB_R,  dθ=ω dT+(η/2πR)dB_θ,  dY=−dT`,
measured the noise-induced early-escape threshold `η_*` vs rotation rate `ω`, and used
a **narrow weak-window** `a(θ)` so min, mean, median are well separated.

## What the numerics said — my prediction was wrong

```
narrow weak-window:  G_min=0.632  G_median=1.265  G_mean=1.193
 predicted  η_*(fast)/η_*(frozen):  min=0.500   mean=0.943   no-effect=1.000

 ω:      0     1     3    10    40   150        (phase noise OFF, pure rotation)
 η_*:  2.45  2.47  2.48  2.48  2.48  2.48
 ratio: 1.00  1.01  1.01  1.01  1.02  1.01
```

`η_*` is **flat in ω**: ratio ≈ 1.01 out to ω=150. The **min-barrier prediction (0.50)
is decisively refuted**; the result sits at the typical/averaged barrier, not the
minimum. (The earlier symmetric test couldn't tell mean from median — `mean a = median
a = 1` for a cosine — but it already showed the same flatness; the narrow window
separates them and kills the min hypothesis cleanly.)

## Why min is wrong — the mechanism the test revealed

The Laplace/min argument ignores that **canard escape is not instantaneous**: the noise
must push `R` across the separatrix, which takes a finite local Riccati time
`T_loc ~ 1/√(aY)`. A fast-rotating phase dwells in the weak window only for
`δ/ω ≪ T_loc`, so the trajectory **never has time to complete an escape while passing
through the weak phase.** The weak phase is therefore *invisible* to the rare escape;
the threshold is set by the barrier sustained over `T_loc`, i.e. the typical/averaged
value. Finite escape time beats brief weak-phase exposure.

## Corrected verdict on the four approaches

| Approach | Status after the test |
|---|---|
| **C — Fourier-in-phase / averaging** | **Vindicated.** Escape sees the phase-*averaged* barrier (n=0 mode dominates); α=1 ⇒ an effective autonomous canard with averaged coefficients. The right reduction. |
| **Kapitza (wildcard)** | **Vindicated — I had it backwards.** High-frequency averaging → effective *average* barrier is correct, because escape is average-dominated, not min. |
| **A — two-zone decoherence** | **Half right.** The amplitude-gated decoherence (crossover `R_dec ~ η`) is real, but the inner zone uses the *averaged* barrier, not `min G`. |
| **B — Riccati→Schrödinger→Lyapunov** | **Framework correct, my "band-edge=min" reading wrong.** The Riccati→Hill linearization (`u″ = a(θ)Y·u`) is exact and is the route to the prefactor + the log; but the Lyapunov exponent of the swept potential *averages*, it does not select the band edge. The `|ln ε₂|` comes from the ramp/turning-point (Airy) structure, not a localization anomaly. |

## The actual upshot (better than feared)

The α=1 coupled problem is **gentler than the brief made it sound.** Fast rotation
*averages* the phase modulation rather than exploiting its weak point, so at leading
order

```
σ_*(α=1)  ≈  σ_* of an effective autonomous canard with phase-AVERAGED √(ac/b),
```

robust to the rotation rate. So:

- **Methodology that works:** stochastic averaging (Khasminskii) / Kapitza
  effective-barrier (**C + wildcard**), made rigorous by the Fourier-hierarchy
  truncation, with the **Riccati→Schrödinger (B)** linearization as the route to the
  O(1) prefactor and the log. The min/localization idea (**A/B-as-min**) is out.
- **Honest status:** this is a normal-form model test, leading order. It refutes the
  min hypothesis robustly and points to averaging; pinning *which* average (and the
  prefactor/log) is the Path-A inner-solution job, now with a clear target.

## The meta-point

The creative min-barrier idea was natural and wrong, and the simulation caught it in
one clean design. The reason it's wrong — finite escape time vs weak-phase dwell — is
itself the real insight, and it tells you the α=1 case reduces to averaging, which is
tractable. That is the value of testing a conjecture instead of writing it up.

Reproduce: `python3 outputs/alpha1_test3.py` (and `_test{,2}.py`).
