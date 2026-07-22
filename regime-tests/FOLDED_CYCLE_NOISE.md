# Folded limit cycle — stochastic blow-up / Channel-A derivation

**Status:** **Both channels derived and validated on the normal form.**
Channel A (amplitude escape): effective noise scale `η = σ/√ε₂` and
`σ_*^A(θ) = C_q √ε₂ √(a c / b)(θ)`, with the `√(ac/b)` geometry factor
confirmed by a geometry collapse and the `C_q` prefactor shown to be the
fold/canard constant (cross-chapter inheritance). Channel B (phase
diffusion): the iPRC `~ 1/r` makes phase diffusion blow up *in the same K₂
chart*, giving a log-divergent phase-variance integral and
`σ_*^B ≈ 2π√3 √(ac/b) √(ε₂/|ln ε₂|)` — the **same leading order** as
Channel A, separated by a logarithm. The A-vs-B race is **geometry-independent**
and set by `ε₂` alone. **Verdict (§11): a complete *leading-order* picture and an
explicit A-vs-B regime diagram — but read honestly: at α = 2 the channels are
*separable*, so this is two thresholds *raced*, not a coupling *solved* (the
genuine fusion is α = 1, open); the striking "B wins for full models" line rests
on a *borrowed* FHN `C_q`, so the worked Liénard model is load-bearing, not
polish; and novelty stays gated on the Popović email. Not closed.**

This is the one singularity in the fold family nobody had noised. The
deterministic scaffold is Jelbart–Kuehn–Kuntz 2024 (JKK, arXiv:2208.01361,
*J. Nonlinear Sci.* 34:17). Channel A transfers the §3–§4 move of
`CANARD_BLOWUP.md` (push `σ dW` through the rescaling chart, then
Freidlin–Wentzell accumulated variance) onto the JKK folded-cycle charts.

Companion scripts:
- `folded_cycle_normal_form_map.py` — Channel-A regime map on the JKK K₂
  inner SDE: geometry collapse + falsification control + canard-inheritance.
- `folded_cycle_phase_diffusion.py` — Channel-B phase-diffusion integral, the
  log law, and the A-vs-B race.

Companion figures:
- `figures/folded_cycle_normal_form_regime_map.png` (Channel A).
- `figures/folded_cycle_phase_diffusion.png` (Channel B + A-vs-B race).

---

## 0. Step 0 — novelty gate (passed)

- **Scaffold confirmed.** JKK 2024 is the deterministic blow-up of a *regular
  folded limit-cycle manifold* in three-time-scale "semi-oscillatory" systems
  with two small parameters. It is the object we add noise to.
- **Pipeline check.** No published *stochastic* version exists. The nearest
  work — Ahsan–Dankowicz–Kuehn 2025 (SIADS, online May 2025, arXiv:2404.13429)
  — treats noisy limit cycles by **adjoint/covariance**, not blow-up. A
  near-miss, not a collision.
- **Outstanding human check (load-bearing, not automatable):** the Popović email
  to Kuehn — "is the stochastic/blow-up folded limit cycle free, or in Christian's
  pipeline?" Everything in this chapter is *standard ingredients recombined* (JKK
  blow-up + canard accumulated-variance + Floquet/iPRC phase reduction), so a
  referee or Kuehn could know adjacent or unpublished work that a literature search
  cannot see. The novelty claim is **gated** on this email — send it before
  investing in the hard parts. Draft: `POPOVIC_NOVELTY_EMAIL.md`.

## 1. Setup

The JKK prototypical normal form (their eq. (14), tildes dropped) is a
three-time-scale system in cylindrical coordinates `(r, θ, y) ∈ ℝ_{≥0} ×
ℝ/ℤ × ℝ`, with `r` the fast *radius* (oscillation amplitude), `θ` the fast
limit-cycle *phase*, and `y` a slowly drifting *parameter*:

```
r' = −a(θ) y + b(θ) r² + 𝓡_r,            (fast,   amplitude)
θ' = ε₁,                                  (phase,  rotation)
y' = ε₂ ( −c(θ) + 𝓡_y ),                  (drift,  parameter)
```

with `a(θ), b(θ), c(θ) > 0` smooth and **1-periodic in θ**, and two small
parameters ordered

```
0 < ε₂ ≪ ε₁ ≪ 1            ("semi-oscillatory" regime).
```

`ε₁` sets the phase rotation rate; `ε₂` sets the parameter drift. The
remainders are `𝓡_r = 𝓞(r³, y², ry, ε₁r², ε₁y, ε₂)`, `𝓡_y = 𝓞(r, y, ε₁, ε₂)`.

**The folded limit cycle.** At `ε₂ → 0` the system has a 2-D limit-cycle
manifold `S₀ = {F(r,θ,y)=0}`; its fold is the non-hyperbolic *circle*
`S₀ᶜ = {r = 0, θ ∈ ℝ/ℤ, y = 0}` (regular-fold conditions JKK (11)–(12)). The
critical manifold near the fold is `y = (b/a) r² + 𝓞(r³)`, i.e. **`y ~ r²`** —
the fold geometry that fixes the blow-up weights below. Crucially the
singularity is *global in the angle θ*: the folding object is a whole fast
oscillation carrying a phase, not a point.

**Noise (project convention).** Degenerate noise on the fast variable only.
In the reduced coordinates the fast variable is the radius `r`, so

```
dr = ( −a(θ) y + b(θ) r² ) dt + σ dW_t,
dθ = ε₁ dt,
dy = −ε₂ c(θ) dt.
```

(In a Cartesian realisation — e.g. a forced Liénard oscillator — the physical
degenerate fast-variable noise projects onto **both** `r` and `θ`: the radial
projection drives Channel A below, the tangential projection is the Channel-B
source. The O(1) radial-projection factor is absorbed into `C_q`.)

## 2. JKK blow-up

Set `ε₁ = ε^α`, `ε₂ = ε³` (JKK bookkeeping; `ε_KSW = ε₂ = ε³` matches
Krupa–Szmolyan so no fractional powers appear). The weighted blow-up of the
fold circle (JKK (30)) is

```
(r, θ, y, ε) = ( ρ r̄, θ, ρ² ȳ, ρ ε̄ ),     (r̄, ȳ, ε̄) ∈ S²,
```

so the **weights are `r ~ ρ¹`, `y ~ ρ²`, `ε ~ ρ¹`** (the `y ~ r²` fold forces
weight(y) = 2·weight(r)); `θ` is *not* blown up. The whole fold circle blows
up to a **torus of spheres `S² × S¹`**.

The central **rescaling chart K₂** (`ε̄ > 0`, JKK (31)) is
`(r, θ, y, ε) = (ρ₂ r₂, θ₂, ρ₂² y₂, ρ₂)` with `ρ₂ = ε`. After the time
rescaling `dt₂ = ε dt` the desingularised field (JKK K₂, their (39)–(40)) is

```
r₂' = −a(θ₂) y₂ + b(θ₂) r₂² + 𝓞(ε),
θ₂' = ε^{α−1},
y₂' = −c(θ₂) + 𝓞(ε),
```

a **regular** perturbation in `(r₂, θ₂, y₂)`. Its `ε → 0` limit, with the
canonical scaling `T₂ = (abc)^{1/3} t₂`, `r₂ = (ac/b²)^{1/3} R₂`,
`y₂ = (c²/ab)^{1/3} Y₂`, becomes the **Krupa–Szmolyan / canard Riccati fold**

```
dR₂/dT₂ = R₂² − Y₂,
dY₂/dT₂ = −1,                              (Airy-solvable; the canard normal form)
```

— *identical* to `CANARD_BLOWUP.md` §2, but with `θ`-dependent coefficients
`a(θ), b(θ), c(θ)`. Two sub-cases of the phase:

- **α = 2** — `θ₂' = ε`: phase *frozen* in the chart ⇒ a `θ`-family of planar
  canard problems. **This chapter does α = 2 (Channel A).**
- **α = 1** — `θ₂' = 1`: phase *rotates* through the passage (`𝓞(ε⁻²)`
  rotations, JKK Cor. 3.3) ⇒ a non-autonomous, periodically-forced Riccati.
  This is where Channel B lives.

## 3. Degenerate noise pushed through the blow-up

Carry `σ dW_t` into the K₂ chart. With `r = ε r₂` (so `dr = ε dr₂`) and
`dt₂ = ε dt` (so `dW_t = ε^{−1/2} dB_{t₂}`):

```
ε dr₂ = ε²(−a y₂ + b r₂²) dt + σ dW_t
   dr₂ = (−a y₂ + b r₂²) dt₂ + (σ / ε^{3/2}) dB_{t₂}.
```

So the **effective noise amplitude in the rescaling chart is**

```
η = σ / ε^{3/2} = σ / √ε₂.
```

This is the folded-cycle analogue of the fold's `η = σ/√ε`. The content is
*which* small parameter appears: **the drift parameter `ε₂`, not the angular
`ε₁`.** The fast-oscillation timescale drops out of amplitude escape at
leading order (the phase enters Channel A only through the rotating-phase
case α = 1, §7). The result is parametrisation-independent — the bare KSW
balance `δ_r = ε₂^{1/3}`, `δ_y = ε₂^{2/3}`, `δ_t = ε₂^{−1/3}` gives
`η = σ·δ_r^{−1}·δ_t^{1/2} = σ ε₂^{−1/2}` directly.

The drift equation `dy₂ = −c(θ₂) dt₂` is deterministic at leading order
(degenerate noise; `y` only feels `σ` through the higher-order leak, dropped
here).

## 4. Freidlin–Wentzell accumulated variance ⇒ the Channel-A law

Frozen phase (α = 2), constants `a = a(θ)`, `b = b(θ)`, `c = c(θ)`. The
inner SDE is the canard problem with these coefficients:

```
dr₂ = ( b r₂² − a y₂ ) dt₂ + η dB_{t₂},     dy₂ = −c dt₂.
```

**Natural scales** (the `CANARD_BLOWUP.md` §4 balance, carrying a, b, c). The
attracting branch is `r₂ ≈ −√((a/b) y₂)`, so the amplitude in the window is
`r₂_* = √((a/b) y₂_*)`. The Riccati time at `r₂_*` is `1/(2√(ab y₂_*))`;
equating to the drift transit time `y₂_*/c` pins the window:

```
y₂_* = (c² / (a b))^{1/3},      (drift window;  analogue of W_* = λ^{2/3})
r₂_* = (a c / b²)^{1/3},        (escape amplitude;  analogue of V_* = λ^{1/3})
T_win = (a b c)^{−1/3}.         (passage time in t₂;  analogue of λ^{−1/3})
```

These reproduce the JKK canonical-rescaling constants exactly — an internal
check that the noise lives on the same scales as the deterministic passage.

**Accumulated Brownian variance** (not a single jump). The fluctuating `r₂`
accrues `Var(r₂) = η²·T_win` across the window. Escape when this reaches the
escape amplitude `r₂_*²`:

```
η² · T_win  ~  r₂_*²
η² · (abc)^{−1/3}  ~  (ac/b²)^{2/3}
η²  ~  a c / b          ⇒    η_*  =  C_q · √( a c / b ).
```

Hence, in physical units (`η = σ/√ε₂`):

```
┌─────────────────────────────────────────────────────────────┐
│  σ_*^A(θ)  =  C_q · √ε₂ · √( a(θ) c(θ) / b(θ) ).             │
└─────────────────────────────────────────────────────────────┘
```

`C_q` is the same O(1) accumulated-variance constant as the fold/folded node
(§6 shows the inner SDE *is* the canard normal form, so the constant is
inherited, not new).

**Two things are genuinely new relative to the point fold (BG):**

1. **Exponent in the *drift* parameter.** `σ_* ~ √ε₂` — the half-power of the
   smallest (drift) parameter; the angular `ε₁` is absent at leading order
   (frozen phase). Setting `(a,b,c) = (1,1,λ)` recovers the fold law
   `σ_* = C_q √ε₂ · λ^{1/2}` exactly (`√(ac/b) = √λ`), i.e. BG with `ε → ε₂`.

2. **Phase modulation.** Because `a, b, c` are 1-periodic in `θ`, the escape
   threshold is *modulated around the cycle*:

   ```
   σ_*^A(θ) = C_q √ε₂ · G(θ),     G(θ) := √( a(θ) c(θ) / b(θ) ).
   ```

   The cycle escapes most easily at the **weakest phase** `θ_* = argmin_θ G(θ)`.
   This phase-dependence of the critical noise is the fingerprint of folding a
   *cycle* rather than a point; it vanishes (→ a single BG threshold) iff
   `a, b, c` are `θ`-independent. It is also the first place the phase couples
   into Channel A, and the natural hand-off to Channel B.

## 5. Scaled variables for the regime map

Exactly as `CANARD_BLOWUP.md` §5:

```
Θ      =  η / η_*  =  η / √(a c / b)            ( = σ / (√ε₂ · √(ac/b)) ),
R_hit  =  y₂_hit / y₂_*  =  y₂_hit / (c²/ab)^{1/3}.
```

If `√(ac/b)` captures all the geometry dependence, the distribution of
`R_hit` at fixed `Θ` is **independent of `(a, b, c)`** — the collapse
criterion. Operationally: `R ≪ 0` = late / fold-edge escape; `R = 𝓞(1)` =
canonical-window escape; `R ≫ 0` = early noise-driven escape.

## 6. Numerical regime map (normal form)

`folded_cycle_normal_form_map.py` integrates the inner SDE

```
dr₂ = ( b r₂² − a y₂ ) dt₂ + η dB,     dy₂ = −c dt₂,
```

from `r₂ = −√((a/b) y₂)` at `y₂ = 5 y₂_*`, recording `y₂_hit` at the first
up-crossing of the geometry-scaled spike threshold `r₂ = r₂_*` (canonical
amplitude `R = 1`). Resolution is fixed in canonical time so all geometries
are integrated identically (`dT_canon = 10⁻³`, `T_canon_max = 40`, 500
trajectories).

**Sweep.** Six geometries `(a,b,c) ∈ {(1,1,1), (2,1,1), (1,2,1), (1,1,2),
(0.5,1.5,1), (1.7,0.8,1.3)}` × `Θ ∈ logspace(−0.5, 1.0, 10)`.

**Collapse across geometry (median `R_hit` at each `(a,b,c)`, fixed `Θ`):**

```
   Theta  (1,1,1)  (2,1,1)  (1,2,1)  (1,1,2)  (0.5,1.5,1)  (1.7,0.8,1.3)
   0.316   -1.618   -1.634   -1.596   -1.621   -1.634      -1.627
   1.000   -1.334   -1.293   -1.286   -1.350   -1.275      -1.274
   1.468   -0.879   -0.723   -0.903   -0.852   -0.815      -0.771
   2.154   +0.419   +0.367   +0.387   +0.357   +0.322      +0.432
   3.162   +2.493   +2.338   +2.548   +2.374   +2.502      +2.294
  10.000   +4.764   +4.741   +4.790   +4.763   +4.793      +4.769
```

The curves lie on top of one another: **mean cross-geometry std of median R =
0.035** across an interquartile-scale range of geometry. That is the
universality predicted by `η_* = C_q √(ac/b)`.

**Falsification control.** Re-running with the *wrong* normalisation
`η_* = √(a/b)` (the `c`-factor dropped) breaks the collapse — **mean
cross-geometry std = 0.224, ~6× larger**, splaying maximally on the geometries
where `c ≠ 1`. So the `c` genuinely belongs in `√(ac/b)`.

**Critical noise / prefactor.** Pooled median `R_hit` crosses zero at

```
Θ_crit ≈ 1.9       (geometry-scaled spike threshold, canonical R = 1).
```

**`C_q` inheritance (cross-chapter consistency).** The canard/fold normal form
is the special case `(a,b,c) = (1,1,λ)` of this inner SDE. Run with the canard
chapter's *absolute* threshold convention (`V_cross = 1`, i.e. canonical
`R = λ^{−1/3}`), the same script returns

```
λ = 0.01 → Θ_crit = 2.82
λ = 0.02 → Θ_crit = 2.83        ⇒  reproduces CANARD_BLOWUP.md C_q ≈ 2.8.
λ = 0.04 → Θ_crit = 2.69
```

So the folded-cycle Channel-A escape constant is the *same* O(1) constant as
the fold and folded node (the brief's `MMO_CROSSOVER` C_q-inheritance check);
the `Θ_crit ≈ 1.9` of the main sweep is that same constant read off a
different (geometry-scaled) spike threshold — the prefactor's
threshold-dependence is exactly the one documented in `CANARD_BLOWUP.md` §6.
For a full worked model expect `C_q ≈ 8–10`, as for full FHN.

See `figures/folded_cycle_normal_form_regime_map.png`.

## 7. What is new, what is next

**Positioning.**

- vs **Berglund–Gentz (fold)**: their fold is a *point* with no phase; here the
  fold is a *circle*, and `σ_*^A(θ) = C_q √ε₂ G(θ)` is **phase-modulated**.
  BG is the `θ`-independent, single-parameter (`ε₂`) reduction.
- vs **Berglund–Gentz–Kuehn (folded node)**: there the rotation is *slow* and
  the fast subsystem still jumps; here the folding object is a *sustained fast
  limit cycle*. Different singularity.
- vs **JKK (deterministic folded cycle)**: this is its stochastic extension —
  `η`, `σ_*^A`, and the geometry collapse are the new objects.

**The honest scope of this chapter.** This is **Path B** (accumulated variance
+ numerics), Channel A only, frozen phase (α = 2). It is the bounded first
result the brief's §6 asks for.

**Next (not done here):**

1. **Channel B — phase diffusion.** The fast cycle carries a phase; degenerate
   noise gives `D_φ = σ² ⟨Z_v²⟩` (the `TONIC_PHASE.md` engine). Through the
   fold the cycle amplitude → 0, so `Z_v ~ 1/amplitude` becomes singular *in
   the same K₂ inner chart* where Channel A lives. Tracking `D_φ`'s blow-up
   scaling and racing it against the Channel-A escape time is the genuinely
   new coupling — and the place the analysis must *not* collapse back to pure
   amplitude escape (or it re-collides with BG/BGK).
2. **Rotating phase (α = 1).** Then `θ` sweeps `𝓞(ε⁻²)` times through the
   passage, so the escape threshold sees a **phase average** of `G(θ)` rather
   than a frozen value — the A/B interplay entering Channel A directly. The
   inner problem is the non-autonomous Riccati JKK (45)–(47); this is the
   Path-A analog of the `MMO_K2_ROUTE_AB.md` inner solution and is **open**.
3. **Worked example.** Reproduce the law on a periodically-forced Liénard
   system (JKK §5.1) — the cross-model map à la `CANARD_BLOWUP.md` §10 /
   `VDP_CROSSMODEL.md` — to fix the model `C_q` and confirm survival of
   finite-`ε` corrections. (JKK §5.1's explicit Liénard equations were beyond
   reach in this pass; pull them from the v3 PDF, pp. ~22–26.)

## 8. Result summary (for the σ_crit table)

| Singularity | σ_crit law | Mechanism | New vs prior |
|---|---|---|---|
| Fold (BG) | `C_q √ε λ^{1/2}` | accumulated variance, point fold | done (`CANARD_BLOWUP.md`) |
| Folded node (BGK) | `C_q √ε F(μ)`, `F ~ μ^{3/2}` | K₂ staircase dissolution | done (`MMO_NOISE.md`) |
| **Folded limit cycle (this)** | **`C_q √ε₂ √(a(θ)c(θ)/b(θ))`** | **accumulated variance in JKK K₂ chart; `η = σ/√ε₂`; phase-modulated** | **Channel A new; Channel B + interplay open** |

```
η = σ/√ε₂        (effective noise, JKK rescaling chart)
σ_*^A(θ) = C_q √ε₂ √(a(θ)c(θ)/b(θ)),   C_q ≈ 2.8 (normal form, canard threshold)
```

## 9. Channel B — phase diffusion through the fold (derivation)

The folding object is a *sustained fast oscillation*, so it carries a phase `θ`.
Degenerate noise on the fast Cartesian variable (the "noise in v" convention)
projects onto polar `(r, θ)`:

```
x = r cos 2πθ;   noise σ dW on x  ⇒
   radial:  dr += σ cos(2πθ) dW              ⇒ Channel A  (rms σ/√2)
   phase :  dθ += −σ sin(2πθ)/(2π r) dW      ⇒ Channel B  (1/r blow-up)
```

The phase noise carries a **1/r** factor — the standard iPRC `~ 1/(cycle
amplitude)` (`TONIC_PHASE.md` §6). As the cycle folds, `r → 0`, so phase
diffusion becomes singular *exactly where Channel A's escape lives*. Pushing
through the K₂ chart (`r = ε r₂`, `η = σ/√ε₂`, `dt₂ = ε dt`):

```
dθ₂ = ε^{α−1} dt₂ + [ η / (2π r₂) ] · (proj) dB_{t₂},
```

the *same* `η` as Channel A, divided by `2π r₂`. Accumulated phase variance
through the passage (`⟨sin²⟩ = 1/2` over the fast phase):

```
Var(θ) = (η²/8π²) ∫ dt₂ / r₂².
```

On the attracting branch `r₂² = (a/b) y₂`, `dt₂ = −dy₂/c`, so

```
∫ dt₂/r₂² = (b/ac) ∫ dy₂/y₂ = (b/ac) · ln( y₂_in / y₂_* ),
```

a **log-divergent** integral (range of `r₂` from O(1) to the inner scale), cut
off at the fold layer `y₂_*`. With `y₂_in ~ ε₂^{−2/3}` (physical entry O(1)),
`ln ~ (2/3)|ln ε₂|`. Phase randomization (`Var ~ 1`) gives

```
┌─────────────────────────────────────────────────────────────┐
│  σ_*^B  ≈  2π√3 · √( a c / b ) · √( ε₂ / |ln ε₂| ).          │
└─────────────────────────────────────────────────────────────┘
```

**The structural punchline.** Both channels carry the *same* leading order
`√ε₂ · √(ac/b)`:

```
σ_*^A = C_q · √ε₂ · √(ac/b),            C_q ≈ 2.8  (normal form),
σ_*^B = C_B · √(ε₂/|ln ε₂|) · √(ac/b),  C_B = 2π√3 ≈ 10.9,
σ_*^B / σ_*^A = (C_B/C_q) / √|ln ε₂|.
```

The geometry `√(ac/b)` **cancels in the ratio** — so the A-vs-B race is
*geometry-independent*, controlled by `ε₂` alone (through the log) against the
prefactor ratio. B drops below A only when `|ln ε₂| > (C_B/C_q)²`.

## 10. Channel B — numerical validation

`folded_cycle_phase_diffusion.py`; figure `figures/folded_cycle_phase_diffusion.png`.

**(i) The log law.** The (regularized) integral `I_φ = ∫dt₂/r₂²` vs
`ln(y₂_in/y₂_*)`:

```
(1,1,1),  b/ac=1.00:     ε₂=1e-2 → I_φ=3.20 vs (b/ac)lnL=3.07  (ratio 1.04)
                         ε₂=1e-6 → I_φ=9.34 vs 9.21            (ratio 1.01)
(1.7,0.8,1.3), b/ac=0.36: same ratios; slope tracks (b/ac).
```

`I_φ → (b/ac) ln(y₂_in/y₂_*)`, ratio → 1 as `ε₂ → 0`. The `1/r²` blow-up and the
`(b/ac)` prefactor are confirmed (figure, left panel).

**(ii) The A-vs-B race (geometry-independent).** `η_*^B/η_*^A` at `C_q = 2.8`:

```
   ε₂     η_*^B/η_*^A   (identical for both geometries)
  1e-2       1.77
  1e-3       1.46
  1e-4       1.27
  1e-6       1.04
```

Same for both geometries (the `√(ac/b)` cancels), → 1 as `ε₂ → 0`, crossing 1
near `ε₂ ≈ 10⁻⁷` for the normal-form `C_q` (figure, right panel):

```
ε₂ ≳ 10⁻⁷ :  σ_*^A < σ_*^B   — amplitude escape sets σ_*  (Channel A wins)
ε₂ ≲ 10⁻⁷ :  σ_*^B < σ_*^A   — phase randomizes first     (Channel B wins)
```

For a *full* model the crossover moves up — but this was the soft spot, and the
worked model (Part 2, `FOLDED_CYCLE_WORKED_MODEL.md`) has now **deflated it**. The
earlier line *assumed* `C_q ≈ 8–10` (FHN, borrowed) ⇒ `C_B/C_q ≈ 1.2` ⇒ crossover
at `ε₂ ≈ 0.2` ⇒ "B wins across the accessible range." But the **measured**
fold-of-cycles value is `C_q ≈ 5.7`, giving `C_B/C_q ≈ 1.9` and a crossover at

```
|ln ε₂| ≈ (C_B/C_q)² ≈ 3.6   ⇒   ε₂ ≈ 0.03.
```

So with the real number, **Channel B wins only for `ε₂ ≲ 0.03`**, not across the
whole range; amplitude escape (A) sets `σ_*` for `ε₂ ≳ 0.03`. The "B wins for full
models" headline is therefore a *small-ε₂* statement, not a generic one. (A fully
stiffened folded-cycle model could lift `C_q` toward 8 as FHN's finite-ε
corrections did, nudging the crossover back up; `5.7` is the clean fold-of-cycles
value.) This is exactly what the worked model was for — turning an extrapolation
into a measured, and more conservative, statement.

**(iii) Diffusion law + interplay.** The stochastic `(r₂, θ₂)` run gives
`Var(θ) = η² I_φ/8π²` at small `η` (`η=0.25`: 0.0037 measured vs 0.0037
predicted), saturating below the linear extrapolation at large `η` (inner-scale
cutoff + radial-noise-shortened approach). The channels are entangled: whichever
threshold is crossed first depends on `ε₂`, and via the phase-dependent
`σ_*^A(θ)` a randomized phase makes amplitude escape select the weakest phase
`θ_*`.

## 11. Verdict — does Channel B close the question?

**Partly.** It delivers the full *leading-order* picture and the regime diagram
the brief asked for — but it does **not** close the problem.

Settled (Path B, leading order, normal form):

- `η = σ/√ε₂`, `σ_*^A(θ) = C_q √ε₂ √(ac/b)(θ)` — derived + collapse-validated.
- `σ_*^B ≈ 2π√3 √(ac/b) √(ε₂/|ln ε₂|)` — derived + log-law-validated.
- **Two channels raced — not the coupling solved.** At α = 2 the phase is frozen
  in the inner chart, so A and B are derived *independently* and compared; this is
  not a solved coupled problem. Both share `√ε₂ √(ac/b)`, so the A-vs-B boundary is
  geometry-independent, set by `ε₂` vs the prefactor ratio (`|ln ε₂| ≷ (C_B/C_q)²`).
  Normal form: A wins for `ε₂ ≳ 10⁻⁷`. The *shape* of this regime diagram is solid;
  its absolute crossover for full models is the extrapolation flagged in §10(ii).
  The genuine *fusion* of the channels is α = 1 (open, below).
- It does **not** reduce to BG/BGK — both the phase channel and the phase
  modulation of A are absent there.

Open (so it is *not* closed), in the order worth doing them — cheap/grounding/
gating first, hardest last, so the most effort isn't sunk into the riskiest part
before the result is shown to survive on a real system and to be free:

1. **Worked Liénard model — grounds the headline (do first).** Compute `C_q` and
   `C_B` on a periodically-forced Liénard system (JKK §5.1) so the "B wins for full
   models" line stops borrowing FHN's constant. Cheap relative to its payoff;
   without it the headline (§10(ii)) is extrapolation, not result.
2. **Popović email — gates the novelty (do first).** Everything here is standard
   ingredients recombined, so adjacent/unpublished work is entirely plausible.
   Confirm the stochastic/blow-up folded cycle is free before investing further.
   Load-bearing, not a formality.
3. **Rotating phase α = 1 — the genuine fusion (after 1–2).** When the phase
   rotates `𝓞(ε⁻²)` turns, A sees a phase *average/min* of `√(ac/b)` and the two
   noises (same `dW`) no longer separate; the honest object is the non-autonomous
   Riccati JKK (45)–(47) with both noises — outlined, not solved. **Highest ceiling
   — this is the real new mathematics, what makes it "a coupled theory" rather than
   "two thresholds compared" — but the hardest and riskiest; don't sink the most
   effort here until 1–2 land.**
4. **Rigour (Path A) + prefactor sharpening.** A rigorous inner solution through
   the JKK charts (analog of the Wechselberger K₂ solution left open for the folded
   node), and pinning the O(1)/sub-log constants in `C_B = 2π√3`. The deepest item;
   the project's fold/node chapters also stop short of it.

**Bottom line.** Channel B is the genuinely new mechanism, and with it the
folded-limit-cycle stochastic problem now has a complete *leading-order* theory
and an explicit regime diagram — a real, self-contained result. Closing it in the
project's full sense (rigorous inner solution + the rotating-phase coupling + a
worked model) is the next phase, not this one.

## References

- Jelbart, Kuehn, Kuntz (2024). *Geometric Blow-Up for Folded Limit Cycle
  Manifolds in Three Time-Scale Systems.* J. Nonlinear Sci. 34:17.
  arXiv:2208.01361. *(deterministic scaffold)*
- Berglund & Gentz (2006). *Noise-Induced Phenomena in Slow–Fast Dynamical
  Systems*, Ch. 5. *(fold = done)*
- Berglund, Gentz, Kuehn (2015). arXiv:1312.6353. *(folded node = done)*
- Ahsan, Dankowicz, Kuehn (2025). *Adjoint-Based Projections … Stochastically
  Perturbed Limit Cycles and Tori.* SIADS. arXiv:2404.13429. *(noisy cycles by
  covariance, NOT blow-up — the near-miss)*
- Krupa & Szmolyan (2001). *SIAM J. Math. Anal.* 33. *(fold blow-up)*
- Project: `CANARD_BLOWUP.md`, `TONIC_PHASE.md`, `MMO_NOISE.md`,
  `MMO_K2_ROUTE_AB.md`, `MMO_CROSSOVER.md`, `PRIOR_BLOWUP_WORK_AUDIT.md`,
  `FOLDED_CYCLE_NOISE_PROMPT.md`.
```
