# FHN Stochastic Reduction Project — Context Document

**Solomon Williams · University of Edinburgh · Summer 2026**  
**Supervisor: Nikola Popovic**

---

## Research Question

When does the slow-manifold reduction of the stochastic FitzHugh–Nagumo (FHN) model give qualitatively wrong predictions, and what is the precise analytical mechanism?

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

## Key Results So Far

### 1. Regime Map
- Empirical failure boundary (where MFPT_full = T_drift/2) fits σ* ~ C·ε^0.43
- Berglund–Gentz prediction σ ~ ε^(3/2) is wrong for this regime
- The failure boundary exponent 0.43 is explained analytically (see below)

### 2. Analytical Mechanism — Integrated Escape Hazard

The correct framework is **not** Kramers at fixed W. It is the integrated escape hazard over the slow passage near the fold.

Near the fold, define δ = w_fold − w (distance to fold). The barrier scales as:
```
ΔU(δ) ~ B·δ^(3/2)
```

The Kramers/Freidlin–Wentzell escape rate:
```
λ(δ) ~ A₀·δ^(1/2) · exp(−K·δ^(3/2) / σ²)
```

Since ẇ = O(ε), the accumulated escape hazard before deterministic tipping:
```
H = ∫ λ dt ~ (σ²/ε) · ∫ A₀·δ^(1/2) · exp(−K·δ^(3/2)/σ²) dδ
```

Substituting z = K·δ^(3/2)/σ², the integral collapses to ~Cσ², giving:
```
H ~ C·σ²/ε
```

Setting H ~ 1 gives the **failure boundary**:
```
σ_crit ~ √(ε/C),   C = K_fold / A₀
```

where K_fold is the near-fold barrier coefficient and A₀ is the Kramers prefactor — both computable purely from the deterministic geometry with no free parameters.

**The 0.43 exponent** is a finite-ε artefact: forcing a power law onto the log-corrected Kramers formula σ* ~ √(2ΔU/log(T_drift·A)) over a finite ε range gives effective exponent ~0.457, consistent with empirical 0.43. As ε → 0 the true exponent approaches 0.5.

### 3. w_escape Distribution — Gaussian OU Model

The non-homogeneous Poisson model (treating w as deterministic during escape) fails to describe the w_escape distribution — KS ~ 0.9.

The correct model: noise in v leaks into w through dw = ε(v + a − bw)dt, giving an **Ornstein–Uhlenbeck process** for w near W_FP:
```
d(w − W_FP) = −bε(w − W_FP) dt + εσ dW_t
```

The **w_escape distribution** is approximately:
```
w_escape ~ N(W_FP,  εσ²τ_v / (2b) · (1 − exp(−2bε·T_esc)))
```

where τ_v ≈ 1.16 is the v autocorrelation time near the fixed point (empirically estimated).

**Validated:** Gaussian OU model achieves KS ~ 0.1–0.4 across all test pairs, consistently beating the Poisson model (KS ~ 0.6–0.9). The mean is W_FP (not shifted by the deterministic drift) because escape is fast relative to the slow drift timescale in the Kramers regime.

**Residual gap** (~KS 0.2): attributable to the covariance between OU fluctuations and the escape path integral — the noise that drives escape also perturbs w during the fast crossing from V_FP to V_SADDLE. This is a second-order correction noted but not pursued.

### 4. Generalised Failure Boundary Formula

For the excitable regime, the closed-form result is:
```
σ_crit ~ √(ε/C),   C = K_fold / A₀
```

where:
- **K_fold** = coefficient in ΔU(δ) ~ K_fold·δ^(3/2) near the fold
- **A₀** = O(1) prefactor in λ(δ) ~ A₀·δ^(1/2)·exp(...)

Both extractable from `potential_barrier(w)` by fitting near w_fold.

---

## Bifurcation Structure

As I increases from −0.1:

```
Excitable  →(SNIC)→  Tonic spiking  →(Hopf)→  Stable spiral (resonator)
I < I_SNIC           I_SNIC < I < I_Hopf        I > I_Hopf
```

**I_SNIC:** Fixed point annihilates at left fold (v* = −1):
```
I_SNIC = (v_fold + a)/b − v_fold + v_fold³/3  ≈ −0.3333  (for a=0.7, b=0.8)
```

**I_Hopf (finite ε):** tr(J) = 0 → v*_Hopf = −√(1 − bε):
```
I_Hopf(ε) such that fixed point sits at v*_Hopf on right branch
```

**At SNIC:** ΔU → 0, Kramers breaks down entirely, any σ > 0 causes failure.

**At Hopf:** stable spiral with purely imaginary eigenvalues, ΔU = 0, canard territory.

---

## σ_crit Bounds by Regime

| Regime | Mechanism | σ_crit scaling |
|--------|-----------|----------------|
| Excitable | Integrated hazard near fold | ~√(ε/C) |
| SNIC | Barrier vanishes, ISI diverges | → 0 as I → I_SNIC |
| Tonic spiking | Phase diffusion on limit cycle | ~√(ε·T_cycle) |
| Resonator | 2D escape near spiral | ~μ√(ε/ω), μ = Re(eigenvalue) |

---

## Codebase

### Files
```
kernel.py            — core machinery (USE THIS for new code)
simulate.py          — legacy EM simulator (backward compat wrapper)
sweep.py             — legacy sweep (backward compat wrapper)
plot.py              — visualisation functions
run_sweep.py         — main driver (legacy, uses simulate/sweep)
run_diagnostics.py   — Kramers/fold-layer crossover diagnostics
run_I_sweep.py       — NEW: sweeps I across excitable→tonic transition
validate_commitment.py — w_escape distribution validation
plot_predicted_boundaries.py — analytical boundary overlays
verify.py            — legacy sanity checks
verify_kernel.py     — NEW: kernel sanity checks
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

### Key Parameters
| Parameter | Role | Typical range |
|-----------|------|---------------|
| σ | Noise amplitude | 0.05 – 0.45 |
| ε | Timescale separation | 0.01 – 0.25 |
| I | External current | −0.1 (excitable) → 0.5+ (tonic) |
| a, b | Recovery params | 0.7, 0.8 (fixed) |

---

## Roadmap

### Completed
- [x] Phase 1: FHN implementation, Euler–Maruyama, deterministic behaviour
- [x] Phase 2: Full simulation pipeline, MFPT sweep, regime map
- [x] Phase 3: Discrepancy metric, failure boundary, power-law fit
- [x] Analytical mechanism: integrated escape hazard → σ ~ √ε
- [x] w_escape distribution: Gaussian OU model validated
- [x] Kernel refactor: FHN2D dataclass, ISI mode, modular sweep

### In Progress
- [ ] I-sweep: regime maps across excitable → tonic transition (run_I_sweep.py)
- [ ] Tonic spiking: CV-based failure boundary, phase diffusion framework
- [ ] SNIC: failure boundary behaviour as I → I_SNIC

### Future (PhD scope)
- [ ] Resonator regime: 3-variable FHN, canards, Berglund–Gentz–Kuehn machinery
- [ ] Bursting: second slow variable, 2D slow manifold
- [ ] Coherence measure: CV as unified diagnostic across all regimes
- [ ] Higher-dimensional generalisation: n^T GG^T n noise projection

---

## Key References

- Berglund & Gentz (2006). *Noise-Induced Phenomena in Slow-Fast Dynamical Systems*
- Kuehn (2015). *Multiple Time Scale Dynamics*
- Fenichel (1979). GSPT foundational theorem
- Pikovsky & Kurths (1997). Coherence resonance in FHN
- Lindner et al. (2004). Effects of noise in excitable systems
- Freidlin & Wentzell (1984). *Random Perturbations of Dynamical Systems*

---

## Notes for New Conversations

- Always import from `kernel.py`, not `simulate.py` or `sweep.py`
- `FHN2D(I=-0.1)` is the baseline excitable model
- The empirical failure boundary uses `delta_ratio=2.0` (2× MFPT error)
- ISI mode resets to the fixed point (V_FP, W_FP) after each spike
- τ_v ≈ 1.16 (empirical), ≈ 1.73 (analytical: 1/|1−V_FP²|)
- The 0.43 exponent is explained — true asymptotic exponent is 0.5
- w_escape is Gaussian around W_FP, not shifted by deterministic drift
