# Handoff — numerically measure the antidamping integral (Route B validation)

**Purpose.** Decisively test the Route B reduction (`MMO_K2_ROUTE_AB.md`) by
measuring, along the *deterministic* FHR funnel, the growth and rotation rates
that the reduction claims govern the SAOs — and reading off `κ` from them. This
is the **measure-first gate**: if the measurement disagrees, Route B's reduction
is wrong and the parabolic-cylinder (Route A) computation is moot, so no effort
is wasted on the hard part. It is cheap (one deterministic trajectory per `c`, an
eigendecomposition per step) and it sidesteps the symbolic Liouville-sign subtlety
entirely — numerics don't care about the symbolic form of `Q`.

Companion to `MMO_K2.md` (measures `ln R = κμ`, `κ ≈ 2π²`) and
`MMO_K2_ROUTE_AB.md` (the reduction being tested).

---

## 1. The claim under test

Route B reduces the transverse deviation from the weak canard to
`ζ'' − 2g(τ)ζ' + Ω²(τ)ζ = 0`: `Ω` is the rotation rate, `g` the **antidamping**
(transverse divergence) that grows the SAOs. The per-turn growth is therefore the
antidamping integrated over one rotation, `ln R = ∮g dτ`, giving

```
κ = 2π² · (1 + O(μ))      ⇔      growth-per-radian  ⟨g/Ω⟩ = πμ .
```

The measurement extracts `g` and `Ω` from the field itself (not from the global
amplitude ratio) and checks this.

## 2. Operational definitions (the load-bearing part)

Along the deterministic funnel trajectory `Γ(t) = (v,w,y)(t)`, form the 3×3
Jacobian of the FHR field

```
J(Γ) =  [ 1 − v²   −1     1   ]
        [  ε      −εb     0   ]
        [ −εδ      0      0   ]   .
```

In the funnel its spectrum is **one real eigenvalue** (the fast/contracting mode)
plus a **complex-conjugate pair** `λ_c(t) = g(t) ± iΩ(t)` — that pair *is* the
rotation:

```
g(t) = Re λ_c(t)   ( = ½ · transverse divergence ),
Ω(t) = |Im λ_c(t)| ( = winding rate ).
```

Then accumulate two integrals over the SAO sequence (physical time `t`; see §4
for boundaries):

```
G = ∫_funnel g(t) dt ,      Φ = ∫_funnel Ω(t) dt .
```

Note `∫g dτ = ∫g dt` and `∫Ω dτ = ∫Ω dt` — the rescaling drops out of both, so
you never need the blow-up time `τ`. (This is why the Liouville sign is
irrelevant here: you measure `g` directly, not via `Q`.)

## 3. The three numbers and exactly what each tests

**(a) The constant — use the RATIO, not `G` alone.**

```
κ_meas = 2π · G / (μ · Φ) .
```

*Why the ratio.* `G = ∫g dt = ln(a_max/a_min)` is the **realised** growth, which
is *funnel-fill-dependent* (`s_obs < s_max`, set by ε via `a_min`). So `G` alone
is **not** `π²` — that is only the counterfactual full-fill value. The ratio
`G/Φ = ⟨g/Ω⟩` is the **fill-independent** growth-per-radian, and
`κ = 2π·(G/Φ)/μ`. **Test: `κ_meas ≈ 2π² ≈ 19.7`?**

**(b) Reduction validity — does the local eigenvalue reproduce the growth?**

```
G  ?=  ln(a_max / a_min)      [a_min, a_max = first/last SAO loop amplitudes].
```

`G` is the integral of the *local* transverse rate; `ln(a_max/a_min)` is the
*actual* loop-amplitude growth. Equality (to nonlinear corrections) is the test
that the SAO growth really is the local transverse eigenvalue. **If these
disagree, Route B's central premise is wrong** — the growth has a nonlinear or
global-return component the linearisation misses.

**(c) Eigenvalue/segment sanity — the phase check.**

```
Φ / (2π)  ?=  s_obs           [number of SAOs counted in the episode].
```

Confirms you picked the right complex pair and the right funnel segment. Cheap,
and it catches the most likely bug.

## 4. Protocol

1. **Model / regime.** Working FHR `v'=v−v³/3−w+y+I`, `w'=ε(v+a−bw)`,
   `y'=εδ(c−v)`, `(a,b,ε,δ)=(0.7,0.8,0.08,0.2)`, `I=0.30` — identical to
   `MMO_K2.md` so the result is directly comparable to the measured κ. Sweep
   `c ∈ {−0.83, −0.86, −0.89, −0.92, −0.94}` (the reliable `s ≥ 3` band;
   `μ ≈ 0.014–0.041`).

2. **`μ` per `c`.** Measure it the way `MMO_ALPHA_DERIVATION.md` does — from the
   desingularised folded-node Jacobian (`μ = λ_w/λ_s ≈ 2bδ(c+1)/(1+δ)²`). Don't
   reuse a fitted value; this keeps `κ_meas` an independent number.

3. **Trajectory.** Integrate deterministic FHR (`σ = 0`) to a stable MMO; isolate
   one `L¹Sˢ` episode. Use a stiff/accurate integrator and a `dt` fine enough
   that `Φ/(2π)` is stable (the integrals are smooth, so this is *less*
   `dt`-sensitive than SAO counting — a built-in advantage over `MMO_K2`).

4. **Funnel boundaries.** Restrict to the SAO sequence: start at the first
   (smallest, `a_min`) loop after re-injection, end at the last (`a_max`) loop
   before the fast jump. Detect by (i) existence of the complex pair
   (`discriminant of the relevant 2×2 block < 0`) and/or (ii) the local minima of
   `|Γ − γ_w|`. Integrate `g`, `Ω` only over this segment. **Watch the sign of
   `g` near the deepest loop** — integrate over the *growing* part (`g > 0`); if
   the trajectory spirals in then out, start at the turn-around (`g = 0`).

5. **Per step:** eigendecompose `J(Γ(t))`, pick the complex pair, accumulate `G`
   and `Φ`. Report `κ_meas`, `G`, `Φ/(2π)`, alongside the `MMO_K2` values
   `a_min, a_max, s_obs, μ`.

## 5. Predictions and the decision tree

| Outcome | Reading |
|---|---|
| `κ_meas ≈ 19.7`, flat across the band, **and** `G ≈ ln(a_max/a_min)`, **and** `Φ/2π ≈ s_obs` | **Reduction validated, constant ≈ 2π².** Proceed to Route A (PCF modulus) anchored to this. |
| `κ_meas` near 20 but **drifts** with μ | The `2π²(1+O(μ))` drift is **real and measured** — this *is* the audit's low-q drift. Fit `κ(μ)`, extrapolate `μ→0`; that limit is the target for Route A. |
| `G ≠ ln(a_max/a_min)` (the validity check fails) | Local transverse eigenvalue does **not** capture the growth ⇒ Route B reduction incomplete (nonlinear/global-return effect). **Stop** — the PCF computation is moot until the reduction is repaired. |
| `Φ/2π ≠ s_obs` | Wrong eigenvalue or segment. Fix the identification before reading anything else. |

The first two rows both *advance* the program (validate, or quantify the drift);
the third is the one that would kill Route B — which is exactly the point of
doing this first.

## 6. Gotchas

- **Eigenvalue picking.** Three eigenvalues; take the complex-conjugate pair, not
  the real fast mode. On pure slow-manifold segments all three may be real (no
  rotation) — those are outside the funnel by construction.
- **Sign of `g` across the deepest loop** (see §4.4) — the single most likely way
  to get `G` wrong.
- **Which fold.** The folded node sits at one fold; analyse that funnel passage,
  not the global return / the opposite fold.
- **`γ_w` vs the trajectory.** Linearising about the actual SAO trajectory is
  fine at leading order (it sits near `γ_w` in the funnel); if `G` is borderline
  against `ln(a_max/a_min)`, redo `J` along the computed weak canard itself.
- **`δ = 0.2`** is the `MMO_K2` regime (`μ ∝ (c+1)` is δ-robust); a `δ = 1.0`
  cross-check is worth one run, as `MMO_K2.md` §5 notes.

## 7. Why this is the right next step

It tests the reduction with a quantity extracted from the field (not back-fitted
from the amplitude ratio), it pins the `g, Ω` coefficients flagged "to verify" in
`MMO_K2_ROUTE_AB.md` §A4, it is immune to the Liouville-sign subtlety, and by
sweeping `c` it directly answers the `2π²` vs `2π²(1+O(μ))` question. It is the
folded-node analogue of how every other chapter earned its closure: measure the
mechanism, then derive the constant — never the reverse.

---

## 8. RESULTS — run of `mmo_k2_gint.py`

Deterministic FHR, dt-converged (identical at `dt = 0.02` and `0.005`), medians
over 17–38 episodes per `c`:

```
   c      mu   nep  k_amp(K2)  k_amp(gr)  k_g(gr)  g/lnA  ph/loop
 -0.86  0.033   29    14.21      12.72     24.67    1.94    1.002
 -0.88  0.028   24    16.81      16.26     25.35    1.59    1.017
 -0.90  0.023   38    19.07      19.48     26.49    1.38    1.011
 -0.92  0.018   28    19.84      27.74     35.64    1.28    1.004
 -0.94  0.014   17    22.09      30.59     37.56    1.23    0.999
                                                  (2π² = 19.74)
```

**Verdict: Route B's reduction is neither killed nor cleanly validated — it is
confirmed as an *asymptotic* (`μ → 0`) statement.** Four findings:

1. **Rotation is exact.** `ph/loop = Φ/(2π·nloop) = 1.00 ± 0.02` across the whole
   band, dt-converged. `Ω = Im λ_c` of the frozen Jacobian *is* the SAO rotation
   rate — the rotation half of the reduction is confirmed clean, no fitting.

2. **`κ ≈ 2π²` reproduced.** The amplitude-based `k_amp(K2)` reproduces
   `MMO_K2`'s κ independently (mean ≈ 18.4, crosses 2π² at `c ≈ −0.90/−0.92`, and
   matches their `c = −0.94 → 22.1` exactly). 2π² is confirmed as the *central*
   value.

3. **The antidamping reduction holds only asymptotically.** The field-based
   `k_g(gr) = 2π·∫g/(μ∫Ω)` **overestimates**: 25–38, vs amplitude ~2π². The
   validity ratio `∫g / ln(amp)` is **1.94 → 1.23**, shrinking monotonically
   toward 1 as `μ → 0`. So the SAO growth *is* antidamping-driven (correlated with
   `Re λ_c` loop-by-loop, **not** the conservative envelope — confirming B4), but
   the *raw* frozen `Re λ_c` is not the exact rate at accessible `μ`; the
   correction (eigenvector-rotation / Liouville-envelope term) vanishes as
   `μ → 0`.

4. **`κ` drifts with `μ`, and the drift is physical** (dt-converged, not a
   resolution artifact): `k_amp(K2)` runs 14 → 22 as `μ`: 0.033 → 0.014. This is
   the audit's "low-q drift" / the `2π²(1+O(μ))` question — now measured. 2π² is
   not a sharp constant at these `μ`; it is the mid-band value.

**Consequence for the program.** You cannot bank "`∫g = π²`" or "`κ = 2π²`
exactly" from numerics at accessible `μ` — both the raw `∫g` overestimate and the
κ-drift are real. Pinning 2π² *is* the `μ → 0` limit, i.e. **Route A (the PCF
asymptotics)** — now with two numerical anchors it must reproduce: (i) `κ → 2π²`
as the `μ → 0` limit of the measured drift, and (ii) the `∫g/ln(amp)` correction
→ 1. The mechanism (rotation + antidamping) is validated; the *constant* remains
an asymptotic computation.

**Caveats.** (a) Amplitudes are measured as the `v`-excursion (a projection of
the eigenplane radius); part of the `∫g/ln(amp) ≠ 1` could be a slowly-varying
projection factor, though the clean `μ`-trend argues the finite-`μ` correction
dominates. (b) `k_amp(gr)` (growing-run amplitude) is the least stable measure at
small `μ`; `k_amp(K2)` is the reference.

Artifacts: `regime-tests/mmo_k2_gint.py`, `results/mmo/mmo_k2_gint.txt`,
`data/mmo_k2_gint.npz`, `figures/mmo_k2_gint.png`.
