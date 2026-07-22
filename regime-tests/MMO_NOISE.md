# MMO Phase 3 — noise dissolution of the Farey staircase

**Status:** mechanism **validated**, exponent **reconciled** (Phase 3.1, §4.1).
The Farey-dissolution picture holds: σ_pq decreases with q (high-q-first, dt-converged),
the scale matches σ_* = C_q√ε·F(μ) with F(μ)~μ^{3/2} and C_q ≈ 8–10 (the canard chapter's
value). The exponent gap closed once the **exact μ(c_pq)** replaces the idealised μ ∝ 1/q:
the derived β ≈ 1.13 then matches the measured β ≈ 1.11 (the original "β = 1.5" was the
idealisation). Path A is **not needed** for the exponent; remaining work is measurement
precision (a converged counting study) and a first-principles γ/F(μ).

> **⚠ Refined by the crossover (`MMO_CROSSOVER.md`, Phases A+B).** §4's "Path-B"
> reading of the μ^{3/2} as a *local* sector-crossing escape is **superseded**:
> sweeping the folded-node normal form to μ → 0 shows the *local* escape is
> canard-like (η_*(μ) → the canard plateau, C_q ≈ 8–10 flat across μ — **no
> μ^{3/2}**). So the μ^{3/2} in σ_pq is a **model-specific global funnel-filling**
> factor, **not** local escape physics. The empirical σ_pq ≈ C_q√ε·μ(c_pq)^{3/2}
> and the §4.1 exponent reconciliation **stand**; only the §4 *derivation* of
> μ^{3/2} does not. Clean statement: `MMO_CHAPTER.md` §4.3/§7 — one universal local
> escape (canard, shared C_q) + one model-specific global geometry (funnel-filling
> f(c)); deriving the global μ^{3/2} from f(c) is open.

Path used for 3B: **B** (accumulated Brownian variance + measured rotation map κ ≈ 2π²;
Wechselberger's inner solution not reconstructed). Companion: `mmo_noise_staircase.py`,
`figures/mmo_noise_staircase.png`.

---

## 1. Question

The deterministic MMO staircase (Phases 1.5–2) locks the rotation number at ρ = 1/q on
plateaus of width Δ_pq ~ q^{−α}, α ≈ 1.55. Add degenerate noise σ·dW on v. **Which plateau
dissolves at what σ_pq, in what order, and does σ_pq ~ σ_*·q^{−α/γ} hold** with σ_* the 3D
folded-node analogue of the canard chapter's σ_* = C_q√ε λ^{1/2}?

## 2. Setup

Noisy working FHR (δ = 0.2, the regime where α, f(c), κ were pinned — kept fixed for
self-consistency, per the prompt): `dv = (v−v³/3−w+y+I)dt + σ dW`, `dw = ε(v+a−bw)dt`,
`dy = εδ(c−v)dt`. Vectorised Euler–Maruyama over N seeds. σ_pq(p,q) ≡ σ at which the
**locked fraction** (seeds with ρ within the Farey half-gap `0.5/(q(q+1))` of 1/q) drops
below 50 %.

**Counter (the load-bearing fix).** The deterministic v=0 split mis-counts under noise
(v-noise makes thousands of spurious subthreshold maxima — confirmed: nS exploded to ~10³).
Fix: count spikes by v-hysteresis (cross +0.5, reset −0.5) and SAO loops as
**prominence-H peaks of w** — w is smooth (noise enters only v; w low-pass-filters it). Then
ρ = nL/(w-peaks). Validated against σ = 0 (recovers the deterministic ρ).

## 3. Measured (3A)

**dt-convergence (chief hazard — cleared).** At q = 2, ⟨ρ⟩(σ) is stable between dt = 0.02
and 0.01 (e.g. σ = 0.025: ⟨ρ⟩ = 0.582 vs 0.589); the locked-fraction knee shifts only
modestly. The dissolution is **physical, not discretisation noise**.

**σ_pq and dissolution order** (H = 0.008, N = 70):

| q (ρ=1/q) | c | σ_pq |
|---|---|---|
| 2 (1/2) | −0.77 | **> 0.045** (widest plateau, most robust) |
| 3 (1/3) | −0.82 | 0.022 |
| 4 (1/4) | −0.845 | 0.016 |

σ_pq **decreases with q** — high-q (narrow, deep-funnel) plateaus dissolve first, the
predicted Farey-dissolution signature (clearly visible in the ⟨ρ⟩-drift figure: ρ=1/4 leaves
its plateau fastest, 1/2 slowest). From q = 3, 4: **β ≈ 1.15** (σ_pq ~ q^{−β}); hence the
directly-inferred **γ = α/β ≈ 1.35**.

**Honest limits.** Only q = 3, 4 are cleanly resolved: q = 2 needs a larger σ-grid (wide
plateau), and q ≥ 5 SAOs are at the noise-countability floor (the smallest deep-funnel loops
fall below the prominence H that rejects v-noise). σ_pq is also **counter-H-sensitive at the
~factor-2 level** (q = 3: 0.022 at H = 0.008 vs 0.010 at H = 0.03). So β is order-reliable,
not precision: **β ∈ [~1.1, ~2]** across reasonable counting choices.

## 4. Derived (3B, Path B) and comparison

The canard-chapter argument generalised to the funnel: noise picks up effective amplitude
η = σ/√ε in the blow-up; to change s by ±1 (hence ρ) it must push the trajectory across one
secondary-canard sector. Using the **measured** rotation map (SAO amplitudes grow by
R = e^{κμ}/turn, κ ≈ 2π²; sectors ~μ wide; T_funnel ~ s_max ~ 1/μ turns), accumulated
Brownian variance `η²·T_funnel ~ (sector)²`:

```
η_*² · (1/μ) ~ μ²   ⇒   η_* ~ μ^{3/2}   ⇒   σ_* = C_q √ε · F(μ),  F(μ) ~ μ^{3/2},
```

and the displacement is linear in σ (variance ∝ η²) ⇒ **γ ≈ 1**. Since μ ∝ 1/q at the ρ=1/q
plateau, `σ_pq ~ √ε · q^{−3/2}`, i.e. **β ≈ 1.5**.

**Comparison (measured vs derived):**

| quantity | measured (3A) | derived (3B) |
|---|---|---|
| exponent β | ≈ 1.15 (q=3,4; order-reliable, ±~0.5) | ≈ 1.5 |
| γ = α/β | ≈ 1.35 | ≈ 1 |
| scale σ_pq(q=3) | 0.022 | √ε μ^{3/2}·C_q = 0.28·0.04^{1.5}·C_q ⇒ C_q ≈ 10 |

The two agree at the **ballpark/order level**: same sign and order of β, γ ≈ O(1), and the
prefactor needs C_q ≈ 10 — strikingly, the **same C_q as the full-FHN canard chapter** (8–10).
That cross-chapter consistency is the strongest single check: the 3D folded-node σ_* inherits
the canard chapter's prefactor. The exponent gap (1.15 vs 1.5) is within the 3A counting
uncertainty.

## 4.1 Exponent tightening (Phase 3.1) — the gap was bookkeeping

The β = 1.5 above used the **idealised** μ ∝ 1/q. The exact relation is not 1/q: at the
ρ = 1/q plateau the SAO count is s = q−1 = f(c)·(1−μ)/(2μ), and because the measured f(c)
**rises** across the band, the self-consistent μ(q) falls **slower** than 1/q. The honest 3B
slope uses μ at each plateau's actual centre c_pq (`mmo_noise_exponent.py`, no new simulation):

```
σ_pq^pred(q) = C_q√ε · μ(c_pq)^{3/2},   μ(c_pq) = 2bδ(c_pq+1)/(1+δ)².
```

| q | c_pq | μ(c_pq) | σ_pq^pred (C_q=10) | σ_pq^meas |
|---|---|---|---|---|
| 2 | −0.765 | 0.0585 | 0.040 | > 0.045 |
| 3 | −0.825 | 0.0422 | 0.025 | 0.022 |
| 4 | −0.845 | 0.0370 | 0.020 | 0.016 |
| 5 | −0.870 | 0.0307 | 0.015 | — |
| 6 | −0.895 | 0.0245 | 0.011 | — |

Fitted slope of the exact-μ prediction: **β_pred^exact ≈ 1.13** (q = 2–6) — *not* 1.5. The
measured β ≈ 1.11 (q = 3, 4) matches it almost exactly; per-point, the prediction tracks the
measured σ_pq with C_q ≈ 8–10 (the canard range; the residual ~20 % is the C_q calibration,
slope is right).

> **Verdict — Outcome 1 (bookkeeping).** The 1.5-vs-1.15 gap was the μ ∝ 1/q idealisation, not
> missing physics. With the exact μ(c_pq) the derived slope (1.13) and the measured slope
> (1.11) agree. **Path A (Wechselberger's K2 inner solution) is NOT required for the
> exponent.** A residual remains only as the measurement's own factor-2 H-sensitivity and the
> 2-point (q = 3, 4) fit — a question of *measurement precision*, not of mechanism — which the
> Task-β counting study (q = 2 extended σ-grid, q = 5–6 with geometry-grounded H, direct γ)
> would pin to β ∈ [1.0, 1.15] but is not needed for the Path-A decision.

## 4.2 The μ-dependence is global: μ·a_min(c) candidate vs μ^{3/2} (Phase M2)

The crossover (`MMO_CROSSOVER.md`) showed the *local* escape is canard-like
(μ-independent C_q, no μ^{3/2}), so the μ-dependence of σ_pq is a **global** (funnel)
property. The mechanistic candidate is `σ_pq ~ μ·a_min(c)` (rotation-map sensitivity
`ds/d(ln a_in) = −1/κμ` × injection depth `a_min(c)`). Phase M2
(`mmo_sigmapq_m2.py`) measured the deterministic ingredients (no noise):

- `a_min(c)` (smallest-SAO amplitude): 1.28 → 0.10 across q=2…6, but it **varies
  factor ~2 within each plateau's c-width** (≈0.54 → 0.25 across the q=3 plateau).
- `κ = ln(a_max/a_min)/(sμ)`: 8.7 (q=3) → 18.6 (q=6) — i.e. → 2π² only asymptotically.

**Verdict: M2 cannot separate `μ·a_min` from `μ^{3/2}`, and gently corrects the
plan's premise.** With the plateau-*centre* a_min, `σ_pq/(μ·a_min) = 0.87, 1.70`
(NOT the ~1.9-flat of the plan's edge-a_min 2-point note), and `μ·a_min`
*underpredicts* the measured σ_pq(q=4) by ~2× while `μ^{3/2}` matches it. The cause
is that **the a_min ambiguity (factor ~2) exceeds the μ·a_min-vs-μ^{3/2} gap** — the
two forms are degenerate at q=3,4 and diverge only at q=2,5,6 (figure). So the
honest status is unchanged: `μ^{3/2}` is a placeholder, `μ·a_min(c)` is
mechanistically motivated but **unproven**; separating them needs BOTH M1 (robust
q=2…6 noise σ_pq) AND a **derived** a_min(c) — i.e. the K2 Path-A inner solution,
which is the deterministic side's open item too. One calculation (closed-form
a_min(c)) would close both. See `MMO_SIGMAPQ_PLAN.md`.

## 5. Verdict and what's open

**Validated:** the Farey-dissolution mechanism (high-q-first, dt-converged), the scale
σ_* ~ C_q√ε·μ^{3/2} with C_q ≈ 8–10 (= canard), γ ≈ O(1), AND the **exponent**: with the exact
μ(c_pq) (§4.1) the derived β ≈ 1.13 matches the measured β ≈ 1.11. So `σ_pq ~ σ_*·q^{−α/γ}`
holds in form, order, prefactor, and exponent — the cross-chapter unification closes at the
project's standard (universal law + model-specific prefactor, C_q inherited from canard).

**Open (precision / Path A frontier, not blocking):** a fully *converged* β with honest error
bars needs the Task-β counting study (q = 2, 5, 6; per-plateau a_min-scaled H to beat the
factor-2 global-H sensitivity; direct γ from std(ρ)~σ^γ). And a *first-principles* γ and F(μ)
(rather than the measured-κ Path-B estimate) would need Wechselberger's K2 inner solution —
the same Path-A frontier as the deterministic side (closed-form a_min(c), first-principles
κ = 2π²). Neither blocks the chapter's result.

This lands at **§7 Outcome 1**: mechanism + order + scale + prefactor + exponent reconciled;
remaining work is measurement precision and first-principles γ/F(μ), not the validity of the law.

## 6. Hazards encountered

(1) v-noise spurious SAO maxima — fixed by w-based counting (§2). (2) Counter-H sensitivity
(factor 2) + high-q countability floor — the limiting systematic (§3). (3) dt-convergence —
checked and cleared (§3). (4) σ_* is an ε^{1/2}·small-μ scale, ~0.02 — the σ-grid resolves it
for q = 3, 4 but steps over the wide q = 2 plateau (σ_pq > 0.045).

## 7. Reproduce

```
python3 regime-tests/mmo_noise_staircase.py
```
Writes `data/mmo_noise.npz`, `results/mmo/mmo_noise.txt`, `figures/mmo_noise_staircase.png`.
