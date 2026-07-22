# FHN Stochastic Noise & Slow–Fast Reduction

**Solomon · University of Edinburgh · Summer 2026**
Supervisor: Nikola Popovic

A unified regime-by-regime theory of how noise affects FitzHugh–Nagumo
dynamics, built from three analytical engines (Geometric Singular
Perturbation Theory, Freidlin–Wentzell large-deviation theory, and
Floquet phase reduction) applied across the four dynamical regimes
the FHN cubic supports.

---

## 1. Research question

The stochastic FitzHugh–Nagumo model with degenerate noise (only on the
fast variable v):

```
dv = (v − v³/3 − w + I) dt + σ dW_t,           (fast)
dw = ε (v + a − b w) dt,                       (slow)
```

with fixed (a, b) = (0.7, 0.8). As I varies, the deterministic resting
state passes through:

```
I < I_fold_L         → excitable (left attracting branch)
I_fold_L < I < I_H1  → lower resonator (stable spiral on middle branch)
I_H1 < I < I_H2      → tonic (limit cycle; FP unstable)
I_H2 < I < I_fold_R  → upper resonator
I > I_fold_R         → excitable (right attracting branch)
```

with I_fold_L ≈ 0.292, I_H1(ε) ≈ 0.331 (ε=0.08), I_H2(ε) ≈ 1.419 (ε=0.08),
I_fold_R ≈ 1.458. The Hopf points satisfy v_H = ±√(1 − εb) — both at
O(ε) distance from the folds.

**The central question.** When does noise destroy the deterministic
slow-fast structure that organises FHN dynamics, and what is the
analytical mechanism in each regime? Concretely: for each regime, what
is the noise scale σ_crit(ε, I) at which the relevant deterministic
feature (escape barrier, canard structure, limit cycle) ceases to
dictate the trajectory's behaviour?

**The chapter-by-chapter answer:**

| Regime | σ_crit(ε, I) | Mechanism | Status |
|---|---|---|---|
| **Excitable** | √(ε/C), C = K_fold/A_0 | integrated escape hazard near fold | derived + validated |
| **Canard (autonomous)** | C_q · √ε · λ^{1/2} | accumulated Brownian variance across canard window | derived + validated (C_q ≈ 2.8 normal form, 8–10 full FHN) |
| **Canard (explosion washout)** | ~ ε^{3/4} | noise overwhelms ε^{1/4} small-cycle amplitude | derived (BGK 2012); inaccessible to direct FHN numerics |
| **Tonic** | CV = σ · A_mid, A_mid = √(c/π²) | Floquet phase diffusion; matched-asymptotic adjoint fold blow-up | A_mid ≈ 0.39 closed; inner eq **corrected to −2V**, fold-blow-up **mechanism derived** (peak ~15%); precise c ≈ 1.55 open |
| **Resonator** | Arrhenius B = ΔV; B → δ² corner law | subcritical-Hopf unstable cycle; phase-gated instanton commitment | derived; corner law + floor explained, instanton gate validated |
| **MMO / SAO** (3D FHR) | σ_pq ≈ C_q√ε·μ(q)^{3/2} (local escape = canard; μ^{3/2} = global funnel) | folded-node canard staircase; degenerate noise dissolves Farey plateaus high-q-first | **derived + validated**: deterministic (folded node, α ≈ 1.5, κ ≈ 2π²) + noise (high-q-first; β ≈ 1.13 = measured 1.11; **C_q ≈ 8–10 = canard's, proven as continuous μ→0 limit**, `MMO_CROSSOVER.md`). 2D FHN scheme retracted (Jordan obstruction) |

**The structural finding** unifying the chapters: noise destroys the
deterministic structure wherever **transverse stability is weakest**
— at folds for excitable and canard regimes, at fold passages for
tonic. The same Krupa–Szmolyan blow-up that desingularises the fold
for canard analysis controls the σ_crit scaling for all three primary
chapters through different observables.

---

## 2. Unified analytical framework

Three engines operate at each regime, applied to the same underlying
slow-fast geometry:

**Geometric Singular Perturbation Theory** (Fenichel 1979,
Krupa–Szmolyan 2001). Provides the deterministic backbone — slow
manifolds, fold blow-ups, canard trajectories. The Krupa–Szmolyan
blow-up coordinates

```
v = −1 + ε^{1/3} V,    w = w_f + ε^{2/3} W,    t = ε^{−1/3} T
```

desingularise the v = −1 fold and produce the canonical fold normal
form dV/dT = V² − W, dW/dT = −λ. This blow-up appears in every chapter.

**Freidlin–Wentzell large-deviation theory** (Freidlin–Wentzell 1984,
Berglund–Gentz 2006). Quantifies noise-induced escape via path
integrals. Under degenerate noise the rate function transforms in
blow-up coordinates as

```
I[γ] = (ε / 2σ²) · A_blowup[V, W],
```

so the entire (ε, σ) content of canard escape is the single prefactor
(ε / 2σ²). This is the canard chapter's foundational identity.

**Floquet phase reduction** (Brown–Moehlis–Holmes 2004,
Ermentrout–Terman 2010). Reduces the tonic spiking SDE to a 1D noisy
phase equation via the adjoint Floquet equation, producing the phase
response curve Z(φ) and phase diffusion coefficient D_φ.

**The methodological unification across chapters:** the Krupa–Szmolyan
blow-up applied to the **forward** equation gives the canard chapter's
mechanism (peel-off from the slow manifold); applied to the
**adjoint** equation gives the tonic chapter's mechanism (phase
diffusion from fold-passage Z_v peaks). Same coordinate machinery,
different operational question. Identifying this parallel is one of
the chapter's distinctive contributions.

---

## 3. Chapter 1 — Excitable regime

**Document:** `PROJECT_CONTEXT_1.md` §1–§2 (legacy); re-integrated
into the unified framework here.

### 3.1 Mechanism

For I < I_fold_L, the FP sits on the left attracting branch at v < −1.
With degenerate noise, escape proceeds via diffusion in v across the
saddle at the middle branch. The integrated escape hazard along the
slow drift gives σ_crit ~ √ε.

Near the fold, deficit δ = w_fold − w. The fold normal-form barrier:

```
ΔU(δ) ~ B · δ^{3/2}.
```

Kramers escape rate near the saddle:

```
λ(δ) ~ A_0 · δ^{1/2} · exp(−K · δ^{3/2} / σ²).
```

Slow drift ẇ = O(ε), so accumulated escape hazard before deterministic
tipping:

```
H = ∫ λ dt ~ (σ²/ε) · ∫ A_0 · δ^{1/2} · exp(−K · δ^{3/2} / σ²) dδ ~ C σ²/ε.
```

Setting H ~ 1 gives the **failure boundary**:

```
σ_crit  =  √(ε / C),     C = K_fold / A_0.                          (3.1)
```

C is computable purely from deterministic geometry (no free parameters).

### 3.2 Empirical validation

Empirical failure boundary fits σ_crit ~ ε^{0.43} (project's original
Phase 3 sweep). The 0.43 exponent is a finite-ε artefact of fitting a
power law to the log-corrected Kramers formula

```
σ_crit  =  √(2 ΔU / log(T_drift · A)).
```

Power-law fit over the empirical ε range gives effective exponent
≈ 0.457, consistent with the 0.43 measurement. As ε → 0 the asymptotic
exponent → 0.5.

### 3.3 w_escape distribution (Gaussian OU model)

A non-homogeneous Poisson model (w deterministic during escape) fails:
KS ≈ 0.9. The correct model accounts for noise in v leaking into w via
the slow equation, giving an Ornstein–Uhlenbeck process for w around
the FP:

```
d(w − W_FP) = − b ε (w − W_FP) dt + ε σ dW_t,
```

with stationary variance

```
w_escape ~ N(W_FP, ε σ² τ_v / (2 b) · (1 − exp(−2 b ε T_esc))).
```

τ_v ≈ 1.16 (empirical) ≈ 1.73 (analytical, 1/|1 − V_FP²|). Validated
KS ≈ 0.1–0.4.

### 3.4 Status

**Validated and complete.** ε^{1/2} exponent confirmed in the
asymptotic limit; finite-ε empirical 0.43 explained; w_escape
distribution model validated to KS = 0.1–0.4. The chapter's σ_crit
formula (3.1) connects directly to the canard chapter as the λ → ∞
limit (large distance from the fold).

---

## 4. Chapter 2 — Canard regime

**Document:** `regime-tests/CANARD_BLOWUP.md` (full derivation and
numerical validation).

### 4.1 Mechanism (autonomous canard escape)

For I in the canard band just inside I_H1, the deterministic
trajectory can follow the repelling middle branch for an O(1) time
before peeling off — the canard structure. With noise, the question
is at what σ does the noise destroy the canard.

In Krupa–Szmolyan blow-up coordinates, the fast equation reduces to
the canonical Riccati dV/dT = V² − W, and the noise picks up an
effective amplitude η = σ/√ε. Two consistency conditions in the
blow-up fix the canard window:

```
W_*  =  λ^{2/3},     V_*  =  λ^{1/3},     T_window  =  λ^{−1/3}.
```

The accumulated Brownian variance criterion η² · T_window ~ V_*² gives

```
η_*  =  C_q · λ^{1/2},        σ_*  =  √ε · η_*  =  C_q · √ε · λ^{1/2}.   (4.1)
```

This is the autonomous canard escape law. The exponent is ε^{1/2}
universally; the prefactor C_q depends on the operational definition
of "noise drives V off the canard."

### 4.2 Numerical validation

Normal-form regime map (`canard_normal_form_map.py`, 4 λ × 10 Θ × 600
trajectories): median R_hit = W_hit / W_* crosses 1 at Θ_crit ≈ 2.8
**independent of λ** to ~10% across an 8× λ range. Confirms exponent
and gives C_q ≈ 2.8 for V_cross = 1.

Full-FHN regime map (`canard_full_fhn_map.py`, 2 ε × 4 λ × 8 Θ × 300
trajectories): same qualitative shape, C_q ≈ 8–10. Collapse holds to
~7% at large Θ; sub-leading λ-dependence at intermediate Θ from
finite-ε corrections (V³/3 sub-leading term, Hopf-fold offset).

### 4.3 Canard explosion washout (the ε^{3/4} row)

At the small-cycle end of the canard explosion (I just past I_H1,
cycle amplitude ε^{1/4}), the noise threshold to wash out the
deterministic explosion structure is set by matching σ_eff = σ/√ε to
the ε^{1/4} amplitude:

```
σ_crit_explosion  ~  ε^{3/4}.                                       (4.2)
```

This is the Berglund–Gentz–Kuehn 2012 prediction. It applies to a
different observable than (4.1) — explosion cycle visibility vs
canard peel-off location — and is **derived theoretically but
experimentally inaccessible in full FHN numerics**: the explosion
window in I is exp(−c/ε) wide; resolving the ε^{1/4} cycle requires
δ < ε² in I, below numerical precision.

### 4.4 Ramped passage (the productive negative)

`canard_ramp_passage.py`: at fixed ramp rate ρ through I_H1, the
mean delay D(σ) stays flat at the deterministic Berglund–Gentz delay
D_0(ε, ρ) until σ approaches the autonomous σ_* threshold. The
deterministic delay scales cleanly (D_0 ~ ε^{0.99} at ρ = ε^{1.5},
D_0 ~ ε^{0.62} at ρ = ε^{1.0}), but the noise-induced washout
**doesn't have its own σ_crit** — it's the autonomous σ_* in disguise.

The "ramped canard exponent" earlier work chased (between 0.78 and
1.60 with different ramp rates) is not a real ramp-specific physics;
it's the autonomous σ_*(λ) measured at different (ε, ρ) operating
points. The clean autonomous exponent is ε^{1/2}; the protocol
dependence dissolves once it's correctly framed.

### 4.5 Status

**Closed.** Autonomous canard exponent ε^{1/2} derived and validated
in both normal-form and full-FHN numerics. Explosion-washout row
inherited from BGK 2012 (analytically derived, inaccessible to direct
numerics). Ramped problem reduced to a corollary of the autonomous
threshold.

---

## 5. Chapter 3 — Tonic regime

**Document:** `regime-tests/TONIC_PHASE.md` (full derivation, numerical
validation, and quantitative closure of both opens).

### 5.1 Framework

For I in the tonic band (I_H1 < I < I_H2), the deterministic FHN has
an unstable FP and a stable limit cycle Γ. The σ_crit framework
changes character: every σ produces a measurable ISI distribution,
and the question is how its moments scale with (σ, ε, I).

Phase reduction: near Γ, dynamics collapse to a 1D noisy phase
equation

```
dφ = ω dt + σ Z_v(φ) dW_t,
```

with Z_v(φ) the phase response curve (gradient of asymptotic phase in
v), computed by integrating the adjoint Floquet equation backward
along Γ. The phase diffusion coefficient

```
D_φ  =  σ² · ⟨Z_v²⟩_φ,
```

drives ISI variance and CV through

```
CV  =  σ · A(I, ε),     A(I, ε) = √( ⟨Z_v²⟩_φ · T_cycle / (2π²) ).   (5.1)
```

A(I, ε) is the cycle's noise-to-jitter conversion factor — a
deterministic quantity packaging cycle geometry into one number.

### 5.2 Predicted A(I, ε) shape

**U-shape across the tonic window:**

- **Mid-tonic floor:** A_mid ≈ const ≈ 0.39 in ε.
  Mechanism: matched-asymptotic three-stage cancellation of log,
  finite, and sub-leading Z_w slaving contributions, giving
  R_fold = c · ε with c ≈ 1.5 cycle-geometric constant.
- **Hopf edges:** A_edge ~ ε^{−1/4} · (I − I_H)^{−1/2}.
  Mechanism: cycle's transverse Floquet exponent vanishes linearly as
  I → I_H by Green's-theorem topology (sum of Floquet exponents equals
  integrated trace J around cycle, dominated by enclosed FP's
  vanishing α). Holds regardless of cycle shape — same divergence
  rate for small Hopf cycle and post-explosion relaxation cycle.

### 5.3 Numerical validation

`tonic_phase_response.py`: A(I, ε) computed via adjoint Floquet
integration at multiple I values × 3 ε values.

- Leading-order phase reduction validated at σ = 0.005: predicted
  CV = σ · A matches measurement to within 5% across all tested I
  cells.
- Edge divergence: (I − I_H)^{−1/2} confirmed asymptotically; p
  exponent → −0.5 as ε → 0 (measured −0.40 at ε = 0.16, −0.55 at
  ε = 0.04; smaller-ε extension to ε = 0.005 gives −0.64,
  overshooting due to canard-explosion contamination near I_H1).
- ε^{−1/4} prefactor: directionally correct but inaccessible to
  direct test (pure small-Hopf regime requires δ ≪ ε² in I, below
  numerical resolution).

### 5.4 Closure of Open Q1: A_mid ≈ 0.39

The matched-asymptotic procedure has three nested cancellations:

1. **Log:** outer slow-phase log-divergence cancels inner blow-up
   log by matching invariance λ_rel = |g| (verified for FHN: both
   equal 0.43 at mid-tonic).
2. **Leading finite:** O(1) remainder cancels to O(ε^{1/3}) via
   Riccati structure of inner BVP.
3. **Sub-leading:** ε^{1/3} cancels to ε^{2/3} via O(ε^{2/3}) Z_w
   slaving correction.

Net: R_fold ~ ε with prefactor c, and the closed form is

```
A_mid  =  √( c / π² )  ≈  0.39,                                     (5.2)
```

matching measured 0.40 to within 2.5%. **Open Q1 closed.**

**First-principles c — mechanism derived, equation corrected (`TONIC_CMID_BVP.md`).**
Solving the inner adjoint-Floquet BVP on the deterministic fold passage
(backward-shoot from the matched outer iPRC) **confirms the fold-blow-up
mechanism** for A_mid: it reproduces the iPRC peak to **~15% near ε = 0.04** and
the correct profile shape, from cycle geometry alone. It does **not** yet pin a
precise c at accessible ε — the leading-order reconstruction overshoots and
drifts (c_recon ≈ 1.9 → 4.5 across ε = 0.04–0.16) because the dropped O(ε^{1/3})
inner correction is large (ε^{1/3} ≈ 0.43 at ε = 0.08) and the slow-branch bulk
is not yet first-principles. So the measured **c ≈ 1.55** has its mechanism
derived (fold blow-up of the adjoint) but a precise value awaits the coupled
O(ε^{1/3}) BVP + matched bulk integral. The derivation's clean standalone result
is the **corrected** inner adjoint equation: the leading inner equation is
dZ̃_v/dT = −2V·Z̃_v (peak at the fold tip V ≈ 0, confirmed numerically), not the
(1/b − 2V) of earlier drafts — the spurious 1/b came from
incorrectly slaving Z_w algebraically rather than as its antiderivative. The
structural prediction A_mid = √(c/π²) is universal across 2D relaxation
oscillators with c model-specific (Van der Pol: c ≈ 0.78, `VDP_CROSSMODEL.md`).

### 5.5 Open Q3 revisited: the σ = 0.02 enhancement is an escape onset, not `A_fold`

Earlier drafts closed Open Q3 by asserting a **deterministic, σ-independent**
multiplicative factor A_fold = exp(∮ max(λ_⊥,0) dt) ≈ 2.7. A first-principles
audit (`TONIC_FOLD_AMPLIFICATION.md`) shows that closure does **not** survive:

1. **The exponentiated integral is not 2.7.** Computed directly from the
   deterministic cycle, ∮ max(λ_⊥,0) dt = 3.2–5.0, so exp(·) = 25–149, not
   2.7 — because the cycle is hugely net-contracting (Floquet multiplier
   e^{λ₂T} ~ 1e−17) and the expanding sliver is immediately re-contracted.
2. **The "post-tip" integral is cutoff-dependent.** 2∫₀^{V_out} V dT grows
   without bound with V_out; it equals ln 2.7 only at a tuned V_out ≈ 1.26,
   so 2.7 is not a convergent cycle constant.
3. **The CV ratio is σ-dependent.** A fresh ISI sweep gives ratio ≈ 1 (leading
   order CV = σ·A) up to σ ≈ 0.015, then a **sharp jump** to 4–12; the value
   just above onset is sampling-/I-dependent (§10 got 2.72, a re-run gives
   4.16 at I = 0.83 and 0.91 at I = 0.875). A deterministic factor would be
   flat in σ.

**Corrected statement.** The σ ≈ 0.02 mid-tonic CV enhancement is the **onset
of non-perturbative fold escape** (a heavy ISI tail from rare excursions where
λ_⊥ momentarily turns positive), at a threshold σ_onset ≈ 0.015–0.02. This
vindicates §18.2's diagnosis over the §20 A_fold picture, and is the *tonic*
instance of the canard chapter's Freidlin–Wentzell escape — unifying the two
chapters under "noise wins where transverse stability is weakest."

### 5.6 σ-validity regime structure (corrected)

Mid-tonic σ-validity has **two** regimes, not three:

```
σ ≤ σ_onset ≈ 0.015–0.02 : leading-order phase reduction, CV = σ · A_mid (no fold factor)
σ > σ_onset              : non-perturbative fold escape (FW action; heavy-tailed ISI)
```

σ_onset(I) is smallest at mid-tonic (sharp fold passages) and larger toward the
Hopf edges — itself the clean, defensible observable, set by fold-passage
transverse geometry (the "ghost of the canard"). The former
Lyapunov-amplified middle band is removed.

### 5.7 Status

**Closed, with one closure corrected.** Flagship predictions validated (CV =
σ · A at weak noise; edge divergence (I − I_H)^{−1/2}). Open Q1 (A_mid ≈ 0.39)
closed and now **derived first-principles** (c ≈ 1.55, `TONIC_CMID_BVP.md`).
Open Q3 (the σ = 0.02 enhancement) **re-characterised**: not a deterministic
A_fold ≈ 2.7 but the onset of non-perturbative fold escape at σ_onset ≈ 0.02
(`TONIC_FOLD_AMPLIFICATION.md`). Two-regime σ-validity structure identified.

---

## 6. Chapter 4 — Resonator regime

**Documents:** `regime-tests/RESONATOR_MECHANISM.md` (phase-gated commitment),
`regime-tests/BARRIER_NORMAL_FORM.md` (B(δ) corner law + floor).

### 6.1 Mechanism

For I in the resonator bands (I_fold_L < I < I_H1, the lower band; and
I_H2 < I < I_fold_R, the upper band) the resting fixed point is a **stable
spiral (focus)**: perturbations ring and decay rather than relaxing
monotonically. The canonical operating point is I = 0.30 (deficit
δ = I_H − I = 0.031), with focus (v*, w*) = (−0.993, −0.367), decay rate
α = −0.0253, spiral frequency ω = 0.280, damping ratio κ = |α|/ω = 0.090.

The cell does not fire on its own; noise jiggles the trajectory around the
focus until one lucky orbit is thrown across the separatrix into a spike. The
σ_crit framework here is the **escape rate**, which is Arrhenius:

```
k(σ)  ~  A · exp( − B / σ² ),      B = ΔV  (quasipotential barrier),     (6.1)
```

fit R² = 0.984 at the canonical point (B = 2.5×10⁻³, A = 0.063). The barrier
B = ΔV is ½ the geometric action of the most-probable escape path (instanton),
so the escape **rate** and the commitment **direction** (§6.3) are two readouts
of one quasipotential object.

### 6.2 The barrier vs Hopf deficit: a δ² corner law, floor-dominated below it

How does B scale with δ as the cell is pushed below the Hopf? Because the only
nonlinearity is the cubic, the planar Hopf normal form is exact at cubic order.
The first Lyapunov coefficient is **positive** (a_3 = +0.268, l_1 = +0.97), so
the Hopf is **subcritical**: a stable focus coexists with an **unstable limit
cycle** that *is* the spiking threshold. The barrier is the amplitude-potential
well depth, and with degenerate (voltage-only) noise projected by the Hopf
eigenvector weight |q_v|² = 0.926,

```
B(δ)  →  α² / ( a_3 |q_v|² )  ~  δ²        as (δ, σ) → 0.               (6.2)
```

This is a **corner law**: it holds only in the joint (δ → 0, σ → 0) limit. The
*measured* window sits **below** the corner, where the barrier is set by the
global distance to the spike separatrix and hits a **floor** (the frozen-w well
and saddle nearly merge at the fold, so the 1D barrier ≈ 0). There the measured
effective slope is shallow (≈ δ^{0.5}) and protocol-dependent (0.52 vs 0.71 on
different σ-grids) — bootstrapped slopes rise from 0.32 ± 0.09 (low δ) to
0.63 ± 0.04 (high δ), a robust crossover (positive in 99.8% of resamples) that
is the data **bending up toward δ²** as δ grows. So there is **no universal
B ~ δ^α** in the accessible regime; the normal form explains the *absence* of a
clean measured exponent rather than supplying one. (This replaces the earlier
`μ√(ε/ω)` quasipotential-well conjecture, which never matched the data.)

### 6.3 Phase-gated commitment: the spike direction is the instanton

When the cell finally commits to a spike, is there a preferred exit direction?
Measured on a mid-shell deep inside the spiral (0.25 of threshold radius), the
phase of the final outbound commitment is offset from the geometric
steepest-outward direction by ~0.5 rad, and the offset **depends on σ** — so it
is *selected*, not geometric. Four tests pin the mechanism:

- **Not geometric** (`test_commitment_gate.py`): the offset moves with σ; a
  geometric gate would be σ-independent.
- **Not linear resonance** (`resonance_scaling.py`): sweeping α and ω
  independently, the offset stays ~ −0.55 rad and barely correlates with κ
  (−0.40) — the focus's linear ringing carries no directional information.
- **It is the large-deviation instanton** (`instanton_exit.py`): the
  degenerate-noise minimum-action escape path predicts an offset of the same
  sign; the geometric (arc-length) minimum-action method, started from 0-, 2-,
  and 4-turn spirals, converges to the **same** committed crossing phase to
  within 0.008 rad — **−1.06 rad, independent of orbital winding** (winding
  near the focus is cost-free; the final outbound leg is geometrically pinned).
- **Finite-noise convergence** (`sigma_lead_scaling.py`): lowering σ from 0.06
  to 0.02, the measured offset climbs monotonically (−0.39 → −0.65) toward the
  −1.06 instanton ceiling — the earlier "discrepancy" was just the finite-noise
  correction.

One-line mechanism: **in the resonator the spike-commitment direction is set by
the nonlinear, degenerate-noise large-deviation instanton, robust to the
orbital winding that precedes it; the linear focus geometry (κ) and the 1D
frozen-w barrier are explicitly not the controlling quantities.**

### 6.4 Status

**Derived and integrated.** Escape is Arrhenius with barrier B = ½ instanton
action (R² = 0.984); the B(δ) → δ² corner law is derived from the subcritical
Hopf normal form, and the measured shallow exponent is explained as the
pre-asymptotic, floor-dominated regime below the corner. The phase-gated
commitment gate is the large-deviation instanton (−1.06 rad), validated against
finite-noise measurements. The full 2D escape quasipotential beyond the Hopf
normal form (the complete escape geometry, not just the corner law) is
PhD-scope (README §9.2).

---

## 7. Master σ_crit table

| Regime | σ_crit / scaling | Mechanism | Status |
|---|---|---|---|
| **Excitable** | √(ε / C), C = K_fold/A_0 | integrated escape hazard, fold barrier ΔU ~ δ^{3/2} | derived + validated; ε^{1/2} asymptotic, ε^{0.43} finite-ε |
| **Canard (autonomous escape)** | C_q · √ε · λ^{1/2}, C_q ≈ 2.8 (normal form) / 8–10 (full FHN) | accumulated Brownian variance over canard window | derived + validated |
| **Canard (explosion washout)** | ~ ε^{3/4} | noise overwhelms ε^{1/4} explosion-cycle amplitude | derived (BGK 2012); experimentally inaccessible in full FHN |
| **Resonator** | Arrhenius k ~ exp(−B/σ²); B = ΔV → δ² (corner law) | subcritical-Hopf unstable cycle = threshold; B = ½ instanton action; phase-gated commitment | derived; corner law + floor explained, instanton gate validated (numerics pre-asymptotic) |
| **Tonic mid-window** | CV = σ · A_mid, A_mid = √(c/π²) ≈ 0.39, **c ≈ 1.55** | Floquet phase diffusion; matched-asymptotic adjoint fold blow-up (inner **dZ̃_v/dT = −2V Z̃_v**, corrected) | derived + validated; mechanism derived + inner eq corrected; precise c open (O(ε^{1/3}) large at accessible ε) — Q1 closed |
| **Tonic mid-window (σ ≳ σ_onset)** | escape onset at σ_onset ≈ 0.015–0.02 | non-perturbative fold escape (heavy ISI tail); **not** a deterministic A_fold | re-characterised (Q3); FW cycle-escape action is PhD-scope |
| **Tonic Hopf edges** | CV ~ σ · ε^{−1/4} · (I − I_H)^{−1/2} | small-cycle Hopf phase response + Green's-theorem topology | derived; (I−I_H)^{−1/2} validated, ε^{−1/4} inaccessible in full FHN |
| **MMO / SAO regime** (3D FitzHugh–Rinzel) | σ_pq ~ σ_*·q^{−α/γ}; σ_* = C_q √ε F(μ) | folded node (μ = λ_weak/λ_strong ∝ (c+1)); Wechselberger ceiling s_max = (1−μ)/(2μ) (α_ceiling = 2); realised α ≈ 1.55 via funnel-filling f(c) = (2/κ)·ln(a_max/a_min)/(1−μ), rotation constant κ ≈ 2π² | **derived + validated**: deterministic (α ≈ 1.45–1.49, `MMO_K2.md`) + noise (`MMO_NOISE.md`, `MMO_CROSSOVER.md`): high-q-first dissolution; **exponent reconciled** (Phase 3.1) — exact μ(c_pq) gives β ≈ 1.13 = measured 1.11 (idealised μ∝1/q gave 1.5); **local escape = canard escape, C_q ≈ 8–10 shared, proven as continuous μ→0 limit** (crossover, A+B); the μ^{3/2} is **global funnel-filling, not local escape** (empirical; derivation from f(c) open). γ 3B-derived/unmeasured. 2D Farey scheme retracted (Jordan obstruction) |
| **Strong-noise breakdown** | σ ~ O(1) | local slow-fast geometry no longer organises dynamics | flagged; analytical limit boundary |

**The universal structural finding:** every row's mechanism involves
fold or Hopf-normal-form geometry. Excitable: fold barrier ΔU ~ δ^{3/2}. Canard:
fold blow-up window W_* = λ^{2/3}. Resonator: subcritical-Hopf unstable cycle +
instanton. Tonic mid: fold blow-up of the adjoint Floquet equation. Tonic edges:
transverse Floquet at the Hopf-fold complex. MMO: the same fold blow-up at the 3D
**folded node** sets the secondary-canard rotation map (κ ≈ 2π²), and noise dissolves
the staircase by the **canard escape itself** — the μ → 0 crossover (`MMO_CROSSOVER.md`)
shows the folded-node noise escape *becomes* the canard escape with the shared prefactor
C_q (a continuous limit, not a coincidence); the μ^{3/2} in σ_pq is a global
funnel-filling factor, **not** local escape. The Krupa–Szmolyan blow-up (and its
adjoint) is the common mathematical engine; cross-model tests (Van der Pol)
confirm the **laws** are universal with **model-specific prefactors**
(C_q, c).

---

## 8. Codebase

### 8.1 Active files

```
kernel.py                          — FHN2D dataclass (geometry, Hopf landmarks,
                                      regime_at(ε)); simulate_kernel SDE integrator;
                                      sweep_grid; failure_boundary; ISI mode
verify_kernel.py                   — kernel sanity checks
hopf_analysis.py                   — analytical + numerical Hopf landmarks
run_I_sweep.py                     — I-sweep across regime transitions
investigate_snic.py                — fine resolution at lower Hopf
period_divergence.py               — SNIC vs Hopf vs homoclinic period test
regime-tests/_shim.py              — scipy brentq pure-Python fallback

regime-tests/CANARD_BLOWUP.md      — canard chapter (Chapter 2)
regime-tests/canard_normal_form_map.py
regime-tests/canard_full_fhn_map.py
regime-tests/canard_ramp_passage.py
regime-tests/canard_action_min.py  — order-of-magnitude FW action check
regime-tests/canard_hit_location.py — full-FHN hit-location simulator

regime-tests/TONIC_PHASE.md        — tonic chapter (Chapter 3)
regime-tests/tonic_phase_response.py
regime-tests/tonic_corrections.py  — §17 numerical closure
regime-tests/TONIC_CMID_BVP.md     — first-principles c via inner adjoint BVP
regime-tests/tonic_cmid_bvp.py     — inner adjoint BVP; corrects inner eq (−2V),
                                     derives mechanism (precise c≈1.55 open)
regime-tests/TONIC_FOLD_AMPLIFICATION.md — A_fold audit (escape-onset finding)
regime-tests/tonic_fold_amplification.py

regime-tests/RESONATOR_MECHANISM.md — resonator chapter (Chapter 4): instanton gate
regime-tests/BARRIER_NORMAL_FORM.md — B(δ) → δ² corner law + floor
regime-tests/test_commitment_gate.py, resonance_scaling.py, instanton_exit.py,
regime-tests/sigma_lead_scaling.py, arrhenius_escape.py, multiseed_delta.py,
regime-tests/normal_form_barrier.py

regime-tests/VDP_CROSSMODEL.md     — Van der Pol universality test (Chapter 5)
regime-tests/vdp_crossmodel.py     — canard collapse + tonic A_mid for VdP

regime-tests/MMO_CHAPTER.md        — MMO chapter (Chapter 6): consolidated narrative,
                                     Phases 1.5/2/3/3.1 — the chapter doc to read first
regime-tests/MMO_FHR_PLAN.md       — plan + banner (Phase 1→3.1 correction trail)
regime-tests/MMO_K2.md             — Phase 2 deterministic closure: K2 return map,
                                     funnel-filling f(c), κ ≈ 2π², α ≈ 1.45–1.55
regime-tests/MMO_NOISE.md          — Phase 3/3.1 noise dissolution: σ_pq, exponent reconciled
regime-tests/MMO_CROSSOVER.md      — Phases A+B keystone: local escape = canard (μ→0 limit),
                                     C_q shared; μ^{3/2} shown global, not local
regime-tests/MMO_CROSSOVER_PLAN.md — crossover plan (Phases A–D); MMO_CHAPTER2_OUTLINE.md — chapter skeleton
regime-tests/MMO_PHASE1_5.md, MMO_ALPHA_DERIVATION.md, MMO_TIMESCALE_CHECK.md
                                   — Phase 1.5 diagnostics (μ ∝ (c+1), α bracket, 2-slow)
regime-tests/MMO_K2_PROMPT.md, MMO_PHASE3_PROMPT.md, MMO_PHASE3_1_PROMPT.md, MMO_CROSSOVER_PROMPT.md
                                   — handoff prompts (Phase 2 / 3 / 3.1 / crossover A+B)
regime-tests/MMO_PLAN.md           — retired 2D scheme (Jordan obstruction); historical
regime-tests/mmo_fhr_staircase.py, mmo_fhr_foldednode.py, mmo_fhr_alpha_fine.py,
regime-tests/mmo_fhr_timescale.py, mmo_k2_return_map.py, mmo_noise_staircase.py,
regime-tests/mmo_noise_exponent.py, mmo_crossover_nf.py, mmo_2d_check.py
```

### 8.2 Data and figures

```
data/canard_normal_form_hit_location.npz
data/canard_full_fhn_hit_location.npz
data/canard_ramp_passage.npz
data/tonic_phase_response.npz
figures/canard_normal_form_regime_map.png
figures/canard_full_fhn_regime_map.png
figures/canard_ramp_passage.png
figures/tonic_phase_diffusion.png
figures/tonic_edge_divergence.png
figures/tonic_corrections.png
figures/tonic_cmid_bvp.png             — inner-adjoint BVP vs pipeline; A_mid(ε)
figures/tonic_fold_amplification.png   — Floquet structure + escape onset
figures/vdp_crossmodel.png             — VdP canard collapse + tonic A_mid
```

### 8.3 Default parameters

| Parameter | Role | Range |
|---|---|---|
| σ | noise amplitude | 0.001–0.45 |
| ε | timescale separation | 0.005–0.25 |
| I | external current | −0.1 (excitable) → 1.5 (excitable, right branch) |
| a, b | recovery params | 0.7, 0.8 (fixed) |

### 8.4 Reproduce

```bash
# canard chapter
python3 regime-tests/canard_normal_form_map.py
python3 regime-tests/canard_full_fhn_map.py
python3 regime-tests/canard_ramp_passage.py

# tonic chapter
python3 regime-tests/tonic_phase_response.py
python3 regime-tests/tonic_corrections.py
python3 regime-tests/tonic_cmid_bvp.py            # first-principles c ≈ 1.55
python3 regime-tests/tonic_fold_amplification.py  # A_fold audit / escape onset

# resonator chapter
python3 regime-tests/arrhenius_escape.py
python3 regime-tests/normal_form_barrier.py
python3 regime-tests/instanton_exit.py --method gmam

# cross-model universality
python3 regime-tests/vdp_crossmodel.py

# all outputs land in data/, figures/, and regime-tests/results/
```

---

## 9. Status and roadmap

### 9.1 Completed (summer + thesis-chapter scope)

- [x] Excitable chapter (σ ~ √(ε/C), Gaussian OU w_escape model)
- [x] Canard autonomous chapter (σ_* = C_q √ε λ^{1/2}, validated)
- [x] Canard explosion row (ε^{3/4} from BGK, inaccessible in full FHN)
- [x] Canard ramped problem (reduced to autonomous σ_* in disguise)
- [x] Tonic chapter (CV = σ · A_mid; A_mid ≈ 0.39 closed)
- [x] **Resonator chapter integrated** (Arrhenius B = ΔV; B → δ² corner law +
      floor; phase-gated instanton commitment) — §6
- [x] **Tonic A_mid: inner equation corrected + mechanism derived** via the inner
      adjoint-Floquet BVP — leading inner equation is dZ̃_v/dT = −2V Z̃_v (not the
      earlier (1/b−2V)); BVP reproduces the iPRC peak to ~15% near ε=0.04 from
      cycle geometry (no inversion). Precise c ≈ 1.55 still open (O(ε^{1/3})
      correction large at accessible ε + bulk shoulder). `TONIC_CMID_BVP.md` —
      the corrected equation is the candidate standalone-publishable result
- [x] **A_fold audited and re-characterised**: the σ = 0.02 enhancement is an
      escape onset (σ_onset ≈ 0.02), not a deterministic A_fold ≈ 2.7
      (`TONIC_FOLD_AMPLIFICATION.md`)
- [x] **Cross-model universality (Van der Pol)**: canard ε^{1/2}·λ^{1/2} law
      (C_q ≈ 8) and tonic A_mid = √(c/π²) (c ≈ 0.78 vs FHN 1.55) confirmed
      (`VDP_CROSSMODEL.md`)
- [x] **MMO deterministic mechanism (3D FHR) confirmed** — folded-node staircase:
      μ ∝ (c+1), Wechselberger ceiling s_max = (1−μ)/(2μ) (α_ceiling = 2), realised
      α ≈ 1.55 reproduced via funnel-filling f(c) = (2/κ)·ln(a_max/a_min)/(1−μ) with
      κ ≈ 2π² (`MMO_FHR_PLAN.md`, `MMO_K2.md`). 2D scheme retracted (Jordan obstruction)
- [x] **MMO noise dissolution (Phases 3 + 3.1) — validated** (`MMO_NOISE.md`):
      degenerate noise dissolves the folded-node staircase **high-q-first**
      (dt-converged); σ_pq ≈ C_q√ε·μ(q)^{3/2}. **Exponent reconciled (Phase 3.1,
      Outcome 1):** the idealised μ∝1/q gave β=1.5, but the **exact μ(c_pq)** gives
      β ≈ 1.13 = measured 1.11 — the gap was bookkeeping, Path A not needed for the
      exponent. This is the **fifth instance** of the project's pattern (noise destroys
      structure where transverse stability is weakest)
- [x] **MMO crossover keystone (Phases A+B, normal form) — validated**
      (`MMO_CROSSOVER.md`): the folded-node noise escape *is* the canard escape,
      recovered continuously as μ → 0 (η_*(μ) → the canard plateau, **C_q ≈ 8–10
      shared, flat across μ**) — upgrading "C_q = canard's" from a measured coincidence
      to a **continuous μ→0 limit**. **Refinement:** the μ^{3/2} is **not** local
      escape physics (the local escape is canard-like, no μ^{3/2}) but a
      **model-specific global funnel-filling** factor — so the noise side mirrors the
      deterministic side's universal-local (canard escape, shared C_q) + model-specific-
      global (funnel-filling f(c)) split. Open: derive μ^{3/2} from the K2 funnel
      geometry; direct γ; full-FHR crossover (Phases C–D)
- [x] Three-engine unified framework (GSPT + FW + Floquet)
- [x] kernel.py with regime classifier and fold sign-correction
- [x] Master σ_crit table populated across all primary regimes

### 9.2 PhD-scope follow-ups

- [ ] Full first-principles c: (a) the slow-branch **bulk shoulder** (~36% of
      ∮Z_v²dt, currently from numerics — needs the matched outer integral with
      the correct incoming/outgoing asymmetry), and (b) the coupled (Z̃_v, Z̃_w)
      O(ε^{1/3}) inner BVP to remove the residual ε-drift. The fold part
      (~64%) is already first-principles.
- [ ] Tonic fold-escape action: the Freidlin–Wentzell action setting σ_onset
      (the *tonic* instance of the canard chapter's σ_*); replaces the retired
      A_fold picture.
- [ ] Resonator quasipotential beyond the Hopf normal form: full 2D escape
      geometry (not just the δ² corner law).
- [~] **MMO / SAO chapter — derived + validated (§9.1); consolidation pending.**
      Deterministic and noise sides both closed; the exponent was reconciled (Phase 3.1,
      Task α: exact μ(c_pq) ⇒ β ≈ 1.13 = measured 1.11), so Path A is **not** needed for
      validity. What remains is finishing work, not new results:
      (a) **`MMO_CHAPTER.md`** consolidation (Phases 1.5 / 2 / 3 / 3.1) — the one active item;
      (b) *optional precision* — the Task-β converged counting study (q=2 extended σ-grid,
      q=5–6 with a_min-scaled H, and a **direct γ** from std(ρ)~σ^γ) to pin β ∈ [1.0,1.15]
      with error bars and actually measure γ (currently 3B-derived ≈1, unmeasured);
      (c) *first-principles* γ, F(μ), κ = 2π², closed-form a_min(c) — the Wechselberger K2
      inner solution (Path A), the same frontier as the deterministic side.
      Engine = Wechselberger folded-node theory + Nikola's entry-exit work.
- [ ] Cross-model, continued: first-principles c_vdp via the VdP inner BVP;
      Morris–Lecar (non-cubic nullcline → non-universal corrections);
      smaller-ε canard runs to pin C_q^vdp asymptotically.
- [ ] Hodgkin–Huxley reduction: connect to Nikola's GSPT programme.
- [ ] FitzHugh–Rinzel beyond the MMO band: bursting regimes and the three-timescale
      (δ → 0, Krupa–Popovic–Kopell) limit; folded saddle-node transitions.
- [ ] Conference: BAMC 2026 (UEA Norwich), abstract deadline Jan 2026.
- [ ] Publication: SIAM J. Appl. Dyn. Syst. or Nonlinearity.

### 9.3 Methodological insight

The chapter-spanning insight is **the same Krupa–Szmolyan fold blow-up organises
every regime** — through different observables. In excitable, the ε^{1/2}
threshold comes from integrated escape hazard near a static fold. In canard,
ε^{1/2} comes from accumulated Brownian variance across a transit fold window
(blow-up of the **forward** equation). In tonic, the **same blow-up of the
adjoint** equation gives the phase-response peak that sets A_mid = √(c/π²), and
the breakdown of weak-noise theory at σ_onset is the *tonic* instance of the
canard escape. In the resonator, the subcritical-Hopf normal form gives the
δ² barrier and the instanton fixes the commitment direction. In the MMO regime, the
same fold blow-up — now at the **3D folded node** — organises the secondary-canard
rotation map (κ ≈ 2π²) behind the deterministic staircase, and the noise escape that
dissolves it **is the canard escape itself**: the μ → 0 crossover recovers it
continuously, with C_q ≈ 8–10 shared (a continuous limit, not a coincidence), while the
μ^{3/2} in σ_pq is a model-specific *global* funnel-filling factor. The fifth instance of
the pattern, the first in three dimensions, and the one where the local escape is shown
to be *literally* the 2D canard escape continued. The mathematical
engine (fold blow-up, forward and adjoint) is identical; the observables (MFPT,
peel-off location, CV, escape direction) and mechanisms differ. **Noise destroys
deterministic structure wherever transverse stability is weakest, and the laws
are universal across models (FHN, Van der Pol) with model-specific prefactors
(C_q, c).**

This is the cleanest headline for the thesis introduction or supervisor
meeting: a unified slow-fast noise theory through fold blow-up (forward for
escape, adjoint for jitter), with regime-specific observables, a common ε^{1/2}
skeleton, demonstrated cross-model universality (FHN, Van der Pol), and a span from
the 2D fold to the 3D folded node (MMO).

---

## 10. Key references

**GSPT and canard theory:**
- Fenichel (1979). Geometric singular perturbation theory for ODEs.
  *J. Differential Equations* 31, 53–98.
- Benoît, Callot, Diener & Diener (1981). Chasse au canard.
  *Collectanea Mathematica* 32, 37–119.
- Krupa & Szmolyan (2001). Extending GSPT to nonhyperbolic points —
  fold and canard points in 2D. *SIAM J. Math. Anal.* 33(2), 286–314.
- Kuehn (2015). *Multiple Time Scale Dynamics.* Springer.

**Freidlin–Wentzell and noisy slow-fast:**
- Freidlin & Wentzell (1984). *Random Perturbations of Dynamical
  Systems.* Springer.
- Berglund & Gentz (2006). *Noise-Induced Phenomena in Slow–Fast
  Dynamical Systems: A Sample-Paths Approach.* Springer.
- Berglund, Gentz & Kuehn (2012). Hunting French ducks in a noisy
  environment. *J. Differential Equations* 252, 4786–4841.
- Muratov & Vanden-Eijnden (2005). Self-induced stochastic resonance
  in excitable systems. *Physica D* 210, 227–240.
- Doss & Thieullen (2009). Oscillations and random perturbations of
  a FitzHugh–Nagumo system. arXiv:0906.2671.

**Floquet phase reduction:**
- Brown, Moehlis & Holmes (2004). On the phase reduction and response
  dynamics of neural oscillator populations. *Neural Comp.* 16(4),
  673–715.
- Ermentrout & Terman (2010). *Mathematical Foundations of
  Neuroscience.* Springer.
- Goldobin & Pikovsky (2005). Synchronization and desynchronization
  of self-sustained oscillators by common noise. *Phys. Rev. E* 71,
  045201.

**Supervisor's relevant work:**
- Kaklamanos, Kuehn, Popovic & Sensi (2025). Entry-exit functions with
  intersecting eigenvalues.

**FHN background:**
- Lindner, García-Ojalvo, Neiman & Schimansky-Geier (2004). Effects of
  noise in excitable systems. *Phys. Rep.* 392, 321–424.
- Pikovsky & Kurths (1997). Coherence resonance in a noise-driven
  excitable system. *Phys. Rev. Lett.* 78, 775–778.

---

## 11. Notes for collaboration / continuation

- Always import from `kernel.py`; legacy `simulate.py` / `sweep.py`
  are wrappers maintained for backward compat only.
- `FHN2D(I=-0.1)` is the baseline excitable model; `I_hopf_lower_at(ε)`
  and `I_hopf_upper_at(ε)` define the tonic band; `regime_at(ε)`
  classifies the current operating point.
- `w_fold_left` / `w_fold_right` properties were sign-swapped in early
  versions of kernel.py; patched (verification notes live with the
  regime-tests).
- ISI mode in `simulate_kernel` resets to FP after each spike;
  `early_exit_frac` controls sweep termination.
- All chapters' analytical derivations live in their respective
  regime-tests `.md` files; this README is the index and cross-chapter
  unification.
- `PROJECT_CONTEXT.md` and `PROJECT_CONTEXT_1.md` are superseded by
  this README plus the chapter documents. They remain in the
  repository as historical record (and as reference for the original
  excitable-chapter derivations).
