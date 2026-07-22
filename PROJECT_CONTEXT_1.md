# FHN Stochastic Reduction Project — Context Document

**Solomon Williams · University of Edinburgh · Summer 2026**  
**Supervisor: Nikola Popovic**

---

## Research Question

When does the slow-manifold reduction of the stochastic FitzHugh–Nagumo (FHN) model give qualitatively wrong predictions, and what is the precise analytical mechanism? The broader goal is a unified classification of reduction failure across all dynamical regimes, using GSPT normal forms combined with Freidlin–Wentzell large deviation theory.

---

## Model

**Full 2D SDE (Euler–Maruyama):**
```
dv = (v - v³/3 - w + I) dt + σ dW
dw = ε(v + a - bw) dt
```

**Default parameters:** I = −0.1, a = 0.7, b = 0.8  
**Fixed point:** (V_FP, W_FP) ≈ (−1.2563, −0.6954) on left stable branch  
**Potential barrier:** ΔU(W_FP) ≈ 0.0254  
**Left fold:** w_fold = I + 2/3 ≈ −0.767

**Reduced 1D slow flow (σ → 0):**
```
dw/dt = ε(v_s(w) + a - bw)
```
where v_s(w) is the stable-branch root of the cubic. Fires when w reaches w_fold.

**Reduced MFPT:**
```
T_drift(ε) = ∫_{W_FP}^{w_fold} dw / F_w(w),   F_w = ε(v_s + a - bw)
```
Scales as 1/ε.

---

## Correct Bifurcation Structure (SETTLED — analytic, see hopf_analysis.py)

The earlier "excitable → resonator → tonic" sketch was muddled. The exact
structure follows from the Jacobian of the deterministic flow:

```
J = [[1 - v², -1], [ε, -ε b]]
trace = 1 - v² - ε b ,    det = ε(1 - b + b v²) > 0   (always positive)
```

A Hopf bifurcation requires trace = 0, i.e. **v* = ±√(1 − ε b)**. As I increases
the operating fixed point climbs the middle branch and crosses BOTH roots, so
there are **two Hopf points** bracketing a single wide tonic window:

```
I < I_fold_L            v* < -1                 Excitable (left branch)
I_fold_L < I < I_H1     -1 < v* < -√(1-εb)      Stable spiral  (lower resonator)
I_H1 < I < I_H2         -√(1-εb) < v* < √(1-εb) TONIC limit cycle (FP unstable)
I_H2 < I < I_fold_R     √(1-εb) < v* < 1        Stable spiral  (upper resonator)
I > I_fold_R            v* > 1                  Excitable (right branch)
```

**Bifurcation values (ε = 0.08, a = 0.7, b = 0.8):**
```
I_fold_L  ≈ 0.292   (left-branch FP annihilates — excitable → lower resonator)
I_H1      ≈ 0.331   (lower Hopf — TONIC SPIKING BORN)   v* = -√(1-εb) = -0.9675
I_H2      ≈ 1.419   (upper Hopf — tonic dies, spiral re-stabilises) v* = +0.9675
I_fold_R  ≈ 1.458   (right-branch FP appears — upper resonator → excitable)
```

**Hopf is confirmed (not SNIC):** at trace = 0 the eigenvalues are ±iω with
ω = √det = √(ε(1 − ε b²)) ≈ 0.276, so the limit cycle is born with a **finite**
nascent period 2π/ω ≈ 22.8. A SNIC would give a diverging period. The earlier
period_divergence.py gave a poor fit because it anchored I_c to the left-fold
value (0.292) rather than to the true Hopf at I_H1 ≈ 0.331.

**Correction to the old reading:** the "resonator" band is NOT a single window
on the right branch. There are TWO narrow stable-spiral windows (just inside
each fold, 0.292–0.331 and 1.419–1.458), and the tonic limit cycle fills the
whole wide gap 0.331–1.419 between the two Hopf points. The code's existing
`I_hopf_at` (negative root) was in fact correct for I_H1; the defect was the
`regime` classifier, now fixed via `regime_at(eps)` in kernel.py.

**Canard strip**: exponentially thin O(exp(−C/ε)) region near each Hopf (esp.
I_H1) where trajectories can follow the unstable slow manifold. Stochastic
canards dominate there.

> NOTE: landmark currents above are analytic (exact for the Hopf condition).
> The deterministic limit-cycle amplitude/period confirmation in
> hopf_analysis.py still needs one execution run (sandbox was unavailable when
> this was written) to lock in the I_H2 cycle-death numerics.

---

## Key Results

### 1. Excitable Regime — Integrated Escape Hazard (COMPLETE)

The failure boundary exponent σ* ~ ε^0.43 empirically, explained analytically:

Near the fold, δ = w_fold − w, and:
```
ΔU(δ) ~ B·δ^(3/2)          (fold normal form — universal)
λ(δ)  ~ A₀·δ^(1/2) · exp(−K·δ^(3/2)/σ²)
H     ~ (σ²/ε) · ∫ A₀·δ^(1/2) · exp(−K·δ^(3/2)/σ²) dδ  ~  C·σ²/ε
```

Setting H ~ 1: **σ_crit ~ √(ε/C), C = K_fold/A₀**

The 0.43 exponent is a finite-ε artefact of forcing a power law onto
σ* ~ √(2ΔU/log(T_drift·A)). True asymptotic exponent is 0.5. Power-law fit
to the Kramers formula over the finite ε range gives effective exponent ~0.457,
consistent with empirical 0.43. As ε → 0, exponent → 0.5.

### 2. w_escape Distribution — Gaussian OU Model (VALIDATED)

Non-homogeneous Poisson model (w deterministic) fails: KS ~ 0.9.

Correct model: noise in v leaks into w via dw = ε(v + a − bw)dt, giving OU process:
```
d(w − W_FP) = −bε(w − W_FP) dt + εσ dW_t
```

w_escape distribution:
```
w_escape ~ N(W_FP,  εσ²τ_v/(2b) · (1 − exp(−2bε·T_esc)))
```
τ_v ≈ 1.16 (empirical), ≈ 1.73 (analytical: 1/|1−V_FP²|)

Gaussian OU achieves KS ~ 0.1–0.4, beating Poisson (KS ~ 0.6–0.9) consistently.
Residual gap ~0.2 KS: covariance between OU fluctuations and escape path integral
(second-order correction, noted but not pursued).

### 3. I-Sweep Results (IN PROGRESS)

- Excitable regime: α_LR increases from 0.60 → 1.15 as I → I_left_fold.
  ΔU → 0 destabilises Kramers picture, σ_crit → 0 for all ε.
- Resonator regime (0.292–0.324): trajectories settle to right-branch FP,
  no tonic spiking observed. ISI mode inappropriate here.
- Tonic regime (I > 0.324): CV maps show clean structure.
  At I=0.352: CV boundary at σ ~ 0.06, α_CV ≈ 0.31.
  At I=0.440: CV boundary flat at σ ~ 0.02, weak ε dependence.
  At I=0.690: smooth CV gradient, D_φ ~ σ²/T_cycle (phase diffusion).

### 4. Bifurcation Type — Period Divergence Test (IN PROGRESS)

Period divergence test reveals:
- Period is unmeasurable (>>3000) for I ∈ (0.292, 0.324) — no limit cycle
- Limit cycle born at I ≈ 0.324 with finite period ~50 (Hopf, not SNIC)
- Power law fit: T ~ ΔI^(−0.133), R² = 0.98 — consistent with Hopf (finite period birth)
- SNIC (T ~ ΔI^(−0.5)) and homoclinic (T ~ log(1/ΔI)) both fit poorly (R² negative)
- **The lower transition is not a SNIC — it is the annihilation of the left-branch FP**
  followed by a Hopf at a separate I value that creates the limit cycle.

---

## Unified Analytical Framework

The programme is: derive σ_crit for each bifurcation normal form using
GSPT + Freidlin–Wentzell, then verify across FHN parameter space and other models.

| Normal form | σ_crit | Status | Exponent |
|-------------|--------|--------|---------|
| Fold (excitable) | √(ε/C), C=K_fold/A₀ | DERIVED + VALIDATED | ε^0.5 |
| Hopf (tonic birth) | finite period, phase diffusion | IN PROGRESS | — |
| Resonator bulk | μ√(ε/ω) quasipotential | PLANNED | — |
| Canard strip | ε^(3/4) | CONJECTURED (BG) | ε^0.75 |
| Node transition | — | NOT STARTED | — |

**Key insight:** The fold normal form result (ε^0.5) is universal — it applies to any
model with a fold bifurcation regardless of the global vector field. The prefactor C
encodes model-specific geometry. Testing Van der Pol and Morris-Lecar should confirm
the same exponent with different C.

---

## Codebase

### Active Files (use these)
```
kernel.py                    — core machinery (Layer 1/2/3, USE THIS)
verify_kernel.py             — kernel sanity checks
run_I_sweep.py               — I-sweep across excitable→tonic transition
investigate_snic.py          — fine investigation around SNIC/Hopf transition
period_divergence.py         — SNIC vs homoclinic vs Hopf period test
validate_commitment.py       — w_escape distribution validation
plot_predicted_boundaries.py — analytical boundary overlays on regime map
fhn_animation.html           — interactive phase portrait (open in browser)
```

### Legacy Files (backward compat only)
```
simulate.py    — legacy EM simulator
sweep.py       — legacy sweep
plot.py        — visualisation functions
run_sweep.py   — legacy main driver
run_diagnostics.py — Kramers/fold-layer diagnostics
verify.py      — legacy sanity checks
```

### Kernel Architecture (kernel.py)

**Layer 1 — Model:**
```python
@dataclass
class FHN2D:
    I: float = -0.1
    a: float =  0.7
    b: float =  0.8
    # properties: V_FP, W_FP, w_fold_left, V_SADDLE, delta_U,
    #             I_snic, I_hopf_at(eps), regime
    # methods: dv(), dw(), barrier_at(w)
```

**Layer 2 — Simulator:**
```python
simulate_kernel(model, sigma, eps,
    mode='first_passage'|'isi_sequence',
    n_trajectories, n_isi, dt, T, ...)
# first_passage → mfpt, fpts, w_escapes, w_spikes
# isi_sequence  → isis, isi_mean, isi_std, cv
```

**Layer 3 — Analysis:**
```python
reduced_prediction(model, eps, mode)   # T_drift or T_cycle
sweep_grid(model, sigma_vals, eps_vals, mode, ...)
failure_boundary(sweep_result, delta_ratio=2.0)
compute_cv(isis)
log_ratio(mfpt_full, mfpt_reduced)
```

**NOTE:** FHN2D.regime property still uses old classification. Needs updating
to reflect correct excitable/resonator/tonic structure.

### Key Parameters
| Parameter | Role | Typical range |
|-----------|------|---------------|
| σ | Noise amplitude | 0.02 – 0.45 |
| ε | Timescale separation | 0.01 – 0.25 |
| I | External current | −0.1 (excitable), 0.30 (resonator), 0.40+ (tonic) |
| a, b | Recovery params | 0.7, 0.8 (fixed) |

---

## Roadmap

### Completed
- [x] Phase 1: FHN implementation, Euler–Maruyama, deterministic behaviour
- [x] Phase 2: Full simulation pipeline, MFPT sweep, regime map
- [x] Phase 3: Discrepancy metric, failure boundary, power-law fit
- [x] Analytical mechanism: integrated escape hazard → σ ~ √ε (chapter 1)
- [x] w_escape distribution: Gaussian OU model validated
- [x] Kernel refactor: FHN2D dataclass, ISI mode, modular sweep
- [x] I-sweep infrastructure: run_I_sweep.py with CV-based boundary
- [x] Bifurcation structure corrected: excitable → resonator → tonic

### In Progress
- [ ] Period divergence test: confirm Hopf birth of limit cycle at I≈0.324
- [ ] Correct FHN2D.regime classification in kernel.py
- [ ] Z(θ) phase sensitivity function: compute numerically for tonic regime
- [ ] Phase diffusion coefficient D_φ: derive CV scaling analytically
- [ ] Resonator bulk: quasipotential computation, σ_crit ~ μ√(ε/ω)

### Next (Summer scope)
- [ ] Confirm Hopf vs SNIC vs homoclinic at I≈0.324 definitively
- [ ] Z(θ) peaks at ghost — show fold passage dominates phase diffusion
- [ ] D_φ(I) scaling as I → I_Hopf_right — extract exponent
- [ ] Resonator quasipotential across spiral window [0.324, 1.888]
- [ ] Canard strip: numerical exponent measurement near I_Hopf_right
- [ ] Van der Pol: swap vector field, verify same fold exponent
- [ ] Morris-Lecar: non-cubic nullcline, test universality

### Future (PhD scope)
- [ ] Canard strip rigorous treatment (Berglund–Gentz–Kuehn)
- [ ] Fitzhugh–Rinzel: 3-variable, 2D slow manifold, n^T GG^T n projection
- [ ] Full (a, b) parameter space map: cusp catastrophe, all regime types
- [ ] Hodgkin–Huxley reduced: connect to Nikola's existing GSPT work
- [ ] Publication: SIAM Journal on Applied Dynamical Systems or Nonlinearity
- [ ] Conference: BAMC 2026 (UEA Norwich) — abstract deadline ~Jan 2026

---

## Supervisor Context

**Nikola Popovic** (Edinburgh) works on:
- GSPT and fast-slow systems with fold singularities
- Entry-exit functions with intersecting eigenvalues (2025 paper with Kaklamanos, Kuehn, Sensi)
- Hodgkin-Huxley GSPT reduction (3D reduction, mixed-mode oscillations)
- Collaborator with Christian Kuehn (TU Munich) — the K in Berglund-Gentz-Kuehn

This project is conceived as the **stochastic complement** to Nikola's deterministic
GSPT programme. The entry-exit function (his recent paper) is the deterministic backbone
of the commitment distribution we've been computing stochastically.

Meeting outcomes: extend to resonator, bursting (second slow variable), canards,
coherence measure (CV) as unified diagnostic. Summer = chapter 1 of a PhD programme.

---

## Key References

- Berglund & Gentz (2006). *Noise-Induced Phenomena in Slow-Fast Dynamical Systems*
- Kuehn (2015). *Multiple Time Scale Dynamics*
- Fenichel (1979). GSPT foundational theorem
- Pikovsky & Kurths (1997). Coherence resonance in FHN
- Lindner et al. (2004). Effects of noise in excitable systems
- Freidlin & Wentzell (1984). *Random Perturbations of Dynamical Systems*
- Kaklamanos, Kuehn, Popovic, Sensi (2025). Entry-exit functions, intersecting eigenvalues

---

## Notes for New Conversations

- Always import from `kernel.py`, not `simulate.py` or `sweep.py`
- **Bifurcation structure is excitable → resonator → tonic, NOT excitable → tonic**
- Left-branch FP annihilates at I ≈ 0.292 (not a SNIC in the classical sense)
- Tonic limit cycle born at I ≈ 0.324 via Hopf of right-branch FP
- Resonator window: I ∈ [0.292, 0.324] for ε=0.08 (wide spiral window up to I≈1.888)
- The empirical failure boundary uses `delta_ratio=2.0` (2× MFPT error)
- ISI mode resets to the fixed point (V_FP, W_FP) after each spike
- τ_v ≈ 1.16 (empirical), ≈ 1.73 (analytical: 1/|1−V_FP²|)
- The 0.43 exponent is explained — true asymptotic exponent is 0.5
- w_escape is Gaussian around W_FP, not shifted by deterministic drift
- Spiral window width = I_node − I_Hopf ≈ 0.47 for ε=0.08, grows with ε
- Canard strip is exponentially thin O(exp(−C/ε)) near I_Hopf_right ≈ 0.324
- FHN2D.regime in kernel.py needs updating to reflect correct classification
