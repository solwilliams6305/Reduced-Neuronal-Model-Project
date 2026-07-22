# A normal-form derivation of the escape barrier B(delta)

**Question.** The escape rate obeys a clean Arrhenius law `k(sigma) ~ A exp(-B/sigma^2)`,
with the barrier `B = ΔV` the single quantity that controls how rare spikes are. How
should `B` depend on the Hopf deficit `delta = I_H - I` as the cell is pushed below the
bifurcation — and does that explain why the *measured* exponent refuses to settle on a
universal value?

The answer, derived below and tested against the multi-seed data, is:

> **There is a clean theoretical scaling, `B ~ alpha^2 / (a_3 |q_v|^2) ~ delta^2`, but it
> is a `(sigma -> 0, delta -> 0)` corner law. The measured window sits *below* that
> corner, in a pre-asymptotic regime where the barrier is set by the global distance to
> the spike separatrix, not by the shrinking local unstable cycle. The theory therefore
> explains the *absence* of a universal measured exponent rather than supplying one — and
> it identifies the empirical "crossover" as the data bending up toward `delta^2`.**

## 1. The reduction

The deterministic field has a single nonlinearity — the cubic `-v^3/3` in the voltage
equation; the recovery equation is linear. Near the lower-branch Hopf the resting focus
loses stability at `I_H` where `trace J = 0`, i.e. `v_H = -sqrt(1 - eps*b)`. With
`a=0.7, b=0.8, eps=0.08` this gives

```
I_H = 0.331281 ,  v_H = -0.96747 ,  omega_0 = 0.27551 .
```

Writing the deficit `delta = I_H - I`, the linearization has eigenvalues
`lambda = alpha +- i omega` with `alpha = 1/2 trace J`. Because only one nonlinear term
is present, the planar Hopf normal form in the complex amplitude `z` is exact at cubic
order:

```
z' = (alpha + i omega) z + c_1 |z|^2 z ,   c_1 = a_3 + i b_3 .
```

The radial amplitude `r = |z|` is a gradient flow `r' = -U'(r)` in the amplitude
potential

```
U(r) = -1/2 alpha r^2 - 1/4 a_3 r^4 .
```

## 2. The cubic coefficient and criticality

Computing `a_3` by the standard planar Hopf (Kuznetsov) formula at the Hopf point — using
`f_vv = -2 v_H`, `f_vvv = -2`, all other second/third derivatives zero — gives

```
a_3 = +0.2678   (first Lyapunov coefficient l_1 = +0.972 > 0)  ->  SUBCRITICAL Hopf.
```

The positive sign is the load-bearing fact: a subcritical Hopf with a stable focus
(`alpha < 0`) has a single well at `r = 0` plus a barrier maximum at the **unstable limit
cycle**

```
r_u^2 = |alpha| / a_3 ,
```

and that unstable cycle is exactly the spiking threshold — cross it and the trajectory is
thrown out to the large excursion (the spike). The barrier height is

```
ΔU = U(r_u) - U(0) = alpha^2 / (4 a_3) .
```

## 3. Degenerate (voltage-only) noise

Noise enters only through `v`, so it projects onto the amplitude coordinate with the
`v`-weight of the Hopf eigenvector, `|q_v|^2`, and the `cos^2` phase average `1/2`. The
phase-averaged radial diffusion is `D_r = 1/2 sigma^2 |q_v|^2 * 1/2`. Kramers' law
`k ~ exp(-ΔU / D_r)` then maps onto the measured `exp(-B/sigma^2)` form with

```
B_nf = ΔU / (1/4 |q_v|^2) = alpha^2 / (a_3 |q_v|^2) ,   |q_v|^2 = 0.926 .
```

## 4. The scaling, and why the data don't show it

Numerically `alpha(delta)` is linear over the whole measured range —
`d ln|alpha| / d ln delta = 0.992` — i.e. `alpha ~ -0.81 delta`. The leading normal-form
prediction is therefore unambiguous:

```
B_nf ~ alpha^2 ~ delta^2 .   (slope 2 in log-log)
```

The measurement is nowhere near this:

```
d ln B / d ln delta  = 0.517      (measured)
d ln B / d ln|alpha| = 0.521      (measured)
                     = 2          (normal form)
```

The diagnostic is the ratio `alpha^2 / B`, which the normal form says should be the
constant `a_3 |q_v|^2 = 0.248`. Instead it sweeps across the data:

| delta | alpha | B_meas | alpha^2 / B | B_nf = alpha^2/(a_3 |q_v|^2) |
|---|---|---|---|---|
| 0.010 | -0.0081 | 1.95e-3 | 0.034 | 2.67e-4 |
| 0.015 | -0.0122 | 2.13e-3 | 0.070 | 5.99e-4 |
| 0.020 | -0.0162 | 2.44e-3 | 0.108 | 1.06e-3 |
| 0.030 | -0.0243 | 2.77e-3 | 0.213 | 2.38e-3 |
| 0.050 | -0.0403 | 3.77e-3 | 0.431 | 6.55e-3 |
| 0.080 | -0.0640 | 5.85e-3 | 0.701 | 1.65e-2 |

Read this as follows. At **small** `delta` the normal form *under*-predicts the barrier by
~7x: it lets the barrier collapse as `delta^2` because the unstable cycle radius
`r_u ~ sqrt(|alpha|) ~ sqrt(delta)` shrinks to zero, but the *real* threshold — the
distance the trajectory must travel to the spike separatrix — does not vanish that fast.
The true barrier hits a **floor** set by the global geometry, so the measured curve is far
flatter (`~ delta^0.5`). At **large** `delta` the quartic truncation *over*-predicts
(higher-order amplitude terms, sextic and up, bend the unstable cycle back). The two
curves cross near `delta ~ 0.03` — which is exactly the canonical operating point, and
exactly why the canonical measured barrier `B = 2.5e-3` agrees with the leading normal
form there.

## 5. What this resolves

The earlier empirical findings now have a mechanism rather than a fit:

- **Why no universal exponent.** The clean exponent (`2`) lives in the
  `delta -> 0, sigma -> 0` corner. The measurable window is pre-asymptotic — the barrier
  there is floor-dominated, so the *effective* slope is small and slides with the fit
  window. The protocol-dependence of the measured exponent (0.52 vs 0.71) is not noise;
  it is the signature of measuring a curve in its bend, not on its asymptote.
- **What the empirical "crossover" is.** The local slope rising with `delta`
  (`0.32 -> 0.63`, robust at 99.8%) is the data **bending up toward the `delta^2` normal
  form** as `delta` grows away from the floor. The crossover is the approach to the
  asymptotic regime, seen from below.
- **Why the canonical point is special.** `delta ~ 0.03` is where the floor regime and the
  normal-form `delta^2` regime cross, so the leading normal form is quantitatively right
  there (`2.4e-3` predicted vs `2.5-2.8e-3` measured) and nowhere else.

## 6. The honest statement

We *do* have a normal-form derivation of `B(delta)`: subcritical Hopf (`a_3 = +0.268`),
unstable limit cycle as threshold, degenerate-noise projection `|q_v|^2 = 0.926`, giving

```
B(delta)  ->  alpha^2 / (a_3 |q_v|^2)  ~  delta^2     as (delta, sigma) -> 0.
```

But this is a corner asymptotic, and the data live below it. The measured `~delta^0.5` is
the pre-asymptotic, floor-dominated barrier, and the only way to *measure* the `delta^2`
law is to go to `delta` so small (and `sigma` so small) that escapes become
unobservable — which is precisely the regime where the instanton/normal-form calculation
is needed and the simulation cannot follow. The theory's contribution is not a number to
replace the fit; it is the reason the fit was never going to give a universal number.

## 7. Reproduce

```
python3 normal_form_barrier.py
```

Reads `results/multiseed_delta/barriers.csv`, computes the fixed point, `alpha(delta)`,
the Lyapunov cubic coefficient `a_3`, and the eigenvector v-weight; overlays `B_nf` on the
measured `B(delta)` and reports the log-log slopes. Outputs
`results/normal_form_barrier/{normal_form_barrier.png, .csv}`.
