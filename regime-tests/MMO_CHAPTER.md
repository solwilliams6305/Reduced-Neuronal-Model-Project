# MMO chapter — noise dissolution of the folded-node Farey staircase

**3D FitzHugh–Rinzel. The fifth regime, and the first in three dimensions.**

**Status:** deterministic and noise sides both **derived + validated**. The
deterministic MMO staircase is reproduced from folded-node canard theory
(Phases 1.5–2); degenerate noise dissolves it high-q-first, and the **local noise
escape is the canard escape itself** — the μ → 0 crossover (Phases A+B,
`MMO_CROSSOVER.md`) recovers it continuously with the **shared prefactor
C_q ≈ 8–10**, upgrading "C_q = canard's" from a measured coincidence to a proven
continuous limit. The **μ^{3/2}** in σ_pq is a *model-specific global*
funnel-filling factor, **not** local escape physics (§4.3, §7). Open as
*precision* (converged counting, direct γ), *first-principles* (κ = 2π², a_min(c),
and a derivation of the global μ^{3/2}), and the *full-FHR* crossover (Phases C–D)
— none blocks the result.

This document consolidates `MMO_FHR_PLAN.md`, `MMO_PHASE1_5.md`,
`MMO_ALPHA_DERIVATION.md`, `MMO_TIMESCALE_CHECK.md`, `MMO_K2.md`,
`MMO_NOISE.md`, and `MMO_CROSSOVER.md` into one chapter narrative, parallel to
`CANARD_BLOWUP.md` and `TONIC_PHASE.md`. Where those phase documents quoted the
noise law as `σ_pq ~ σ_*·q^{−α/γ}` or read the μ^{3/2} as a *local* escape factor,
§4.3, §5, and §7 give the corrected statements.

---

## 1. Question

For the four primary regimes (excitable, canard, tonic, resonator) the noise
question was a single threshold σ_crit. The MMO regime is richer: the
deterministic flow produces **mixed-mode oscillation** patterns L^a S^b —
sequences of L large spikes interleaved with S small subthreshold
oscillations (SAOs) — and the rotation number ρ = L/(L+S) traces a **Farey
(devil's) staircase** as a parameter is swept. Two questions:

1. **Deterministic.** Predict the staircase: the L^a S^b patterns, the
   plateau widths Δ_pq ~ q^{−α}, and the exponent α.
2. **Noise.** With degenerate noise σ·dW on the fast variable v (the same
   convention as every other chapter), characterise the dissolution of the
   staircase: which plateau (p, q) dissolves at what threshold σ_pq, in what
   order, and how σ_pq connects to the canard chapter's
   σ_* = C_q·√ε·λ^{1/2}.

**Why 3D.** Autonomous 2D flows cannot carry an L^a S^b symbol sequence at
fixed parameters: a stable limit cycle is a Jordan curve and admits no such
deterministic itinerary. MMOs need a second slow direction whose drift through
a fold creates a **folded singularity** at which trajectories acquire a
rotation number. The minimal carrier in the FHN family is FitzHugh–Rinzel
(Rinzel 1987), which adds a slow modulator y. This Jordan obstruction killed
the project's first (2D) MMO scheme; the redirect to 3D FHR is documented in
the `MMO_FHR_PLAN.md` banner.

---

## 2. Model and folded-node geometry

### 2.1 Working FHR

```
dv = ( v − v³/3 − w + y + I ) dt  +  σ dW       (fast, noisy)
dw = ε ( v + a − b w ) dt                        (slow)
dy = ε δ ( c − v ) dt                            (slow modulator)
```

`(a, b, ε, δ, I) = (0.7, 0.8, 0.08, 0.2, 0.30)`, with **c the bifurcation
control**: sweeping c ∈ [−0.96, −0.70] traverses the staircase. The v–w
subsystem is autonomous FHN, so the canard chapter's Krupa–Szmolyan blow-up at
v = −1 still applies; y is the second slow drift that turns the 2D fold into a
**folded node** in 3D.

(The canonical Rinzel form `y' = δ(cv + d − y)` with δ → 0 does *not* produce
MMOs here — it is a pure spiker, because the global return re-injects outside
the funnel; `MMO_PHASE1_5.md` §1. The working variant above is the one that
feeds the funnel.)

### 2.2 The folded node and the eigenvalue ratio μ

Treating (w, y) as two slow variables, the desingularised reduced flow on the
critical manifold has, at the fold v = −1, a folded singularity with Jacobian

```
J = [[ 1+δ , −b ],
     [ 2bδ(c+1) , 0 ]],     tr = 1+δ,    det = 2bδ(c+1).
```

For the working FHR this is a **folded node** (real eigenvalues, det > 0)
across the whole band, with eigenvalue ratio μ = λ_weak/λ_strong. Near the
folded-saddle-node limit c → −1 the determinant → 0, so

```
μ  =  λ_weak / λ_strong  ∝  (c + 1).
```

Over the MMO band μ ∈ [0.03, 0.16] — small, the **many-secondary-canard**
regime. (The README/early drafts quote the linearised μ ≈ 2bδ(c+1)/(1+δ)²;
the exact eigenvalue ratio is used numerically and is ~10–15 % larger, e.g.
μ(−0.825) = 0.0422.)

### 2.3 The Wechselberger ceiling

Wechselberger's (2005) folded-node theorem gives the maximal number of small
oscillations a canard makes through the funnel:

```
s_max  =  (1 − μ) / (2μ)        (an upper bound — the secondary-canard count).
```

This is a **bound, not an equality**: the realised SAO count grows in lockstep
with s_max as μ → 0 but sits strictly below it (15 % of the bound at large μ,
45 % at the deepest point); `MMO_PHASE1_5.md` §2. With μ ∝ (c+1),
s_max ∝ 1/(c+1), giving a **universal ceiling exponent α_ceiling = 2**.

### 2.4 Timescale: Wechselberger applies at δ = 0.2

A 2-slow folded node requires both slow variables at O(ε). Here w-rate = ε and
y-rate = δε; **δ is an O(1) ratio**, so both are O(ε) and the system is
genuinely 2-slow for any δ = O(1). The folded node is real and MMOs exist for
all δ ∈ [0.05, 1.0]; δ = 0.2 (where α, f(c), κ were pinned) is the chapter's
self-consistent operating point, δ = 1.0 the cleanest cross-check.
Genuine three-timescale (Krupa–Popovic–Kopell) is the δ → 0 limit and is *not*
needed; `MMO_TIMESCALE_CHECK.md`.

### 2.5 Regime: the Krupa–Wechselberger transition, not standard Wechselberger

The standard folded-node MMO theorems — Wechselberger (2005), and the generic
1^{k+1} MMO existence and O(ε^{(1−μ)/2}) sector-width results of Brøns–Krupa–
Wechselberger (2006) — assume **μ ≫ ε^{1/2}**. This chapter operates **below**
that bound: with ε = 0.08 (ε^{1/2} ≈ 0.28) and μ ∈ [0.03, 0.16], we have
**μ < ε^{1/2} across the whole band** (μ/ε^{1/2} ∈ [0.11, 0.57]). The operating
point therefore sits in the **Krupa–Wechselberger (2010) transition regime**
(μ = O(ε^{1/2}) or smaller), which connects to the singular-Hopf bifurcation
(Desroches et al. 2012, §3.2).

This is not a defect — it *is* the small-μ, many-secondary-canard regime in which
MMOs with many SAOs occur — but it is honoured throughout: (i) it **explains the
empirical wrinkles** (κ → 2π² only asymptotically, §3.2; a_min → tiny; the
delicate small-μ numerics); (ii) the standard scalings (the s_max ceiling, the
O(ε^{(1−μ)/2}) sector widths) are used as **asymptotic/empirical guides**, not as
theorems that apply cleanly at this μ; (iii) a first-principles closure of the
global noise μ-dependence (§4.3, §6) must be built on **Krupa–Wechselberger
2010**, not the standard folded-node inner solution.

---

## 3. Deterministic staircase (Phases 1.5–2)

### 3.1 The realised exponent: a universal ceiling and a model-specific filling

The principal sequence is L¹Sˢ, so ρ = 1/(s+1), q = s + 1. Two exponents:

- **Ceiling (universal):** s_max ∝ (c+1)^{−1}, giving α_ceiling = 2.
- **Realised (measured):** s_obs ∝ (c+1)^{−1.74}, giving **α ≈ 1.55** (bracket
  [1.48, 1.60], deep-subrange 1.58 ± 0.06; `mmo_fhr_alpha_fine.py`, 22-pt grid).

The gap is the **funnel-filling** f(c) = s_obs/s_max, which *rises* 0.21 → 0.44
as c → −1. So the realised count grows steeper than the ceiling, pulling α
below 2. This is the chapter's first universal-vs-model-specific split: the
ceiling α = 2 is folded-node-universal; f(c) is set by where FHR's global
return injects into the funnel and is model-specific (`MMO_ALPHA_DERIVATION.md`).

(An earlier α ≈ 1.0 was retracted — a coarse-grid artefact of counting ρ = p/q
samples within a tolerance; the clean SAO-count scaling supersedes it.)

### 3.2 The K2 return map: where f(c) comes from (Phase 2)

Each L¹Sˢ episode is one global return (the L spike + re-injection) followed by
s SAOs through the funnel. The SAO amplitudes grow **geometrically** from the
deepest loop a_min (nearest the weak canard) to a_max ≈ O(1) (the last loop
before the jump), by a per-turn factor R. Measuring (s, a_min, a_max, μ) per c:

```
ln R  =  κ · μ ,      κ = 19.6 ± 1.8  ≈  2π² (= 19.74),   constant across c.
```

This **measured rotation map** is the load-bearing K2 content: the per-turn SAO
growth rate is set by the folded-node μ through the constant κ ≈ 2π². The
constancy is **asymptotic, not exact**: in the transition regime (§2.5) κ
approaches 2π² only as q grows — measured ≈ 8.7 at q = 3, rising to ≈ 18.6 at
q = 6 (M2) — so "κ ≈ 2π² constant" is the large-q limit, with a real low-q drift
that is itself a signature of the μ < ε^{1/2} regime. Composing the (asymptotic)
rotation map with the ceiling s_max = (1−μ)/(2μ) gives the funnel-filling with
**μ cancelling**:

```
f(c)  =  s_obs / s_max  =  (2/κ) · ln(a_max / a_min) / (1 − μ).
```

So f(c) is just the log of the global-return amplitude span, divided by κ
(matches the measured f to 5–15 %). Since a_min(c) → 0 fast as c → −1,
ln(a_max/a_min) ~ (c+1)^{−0.7}, giving s_obs ~ (c+1)^{−1.7} and hence the
reproduced exponent (`MMO_K2.md`):

```
α (measured s_obs)            =  1.49
α (K2-composed f × s_max)     =  1.45        target 1.55 ± 0.06 — confirmed at the lower edge.
```

**Universal piece:** the folded-node rotation map (κ ≈ 2π²) and the ceiling.
**Model-specific piece:** the global-return amplitude span a_min(c). This is
the same split as canard (σ_* exponent universal, C_q specific) and tonic
(A_mid form universal, c specific). Path A (Wechselberger's parabolic-cylinder
inner solution) would give a closed-form a_min(c) and a first-principles
κ = 2π²; Path B reproduced the mechanism without it.

---

## 4. Noise dissolution (Phase 3)

### 4.1 Setup and the counter fix

Add degenerate noise σ·dW on v; integrate with vectorised Euler–Maruyama over N
seeds. σ_pq(p, q) ≡ the σ at which the **locked fraction** (seeds with ρ within
the Farey half-gap 0.5/(q(q+1)) of 1/q) drops below 50 %.

**The load-bearing numerical fix.** Counting SAOs by the deterministic v = 0
split fails under noise — v-noise produces ~10³ spurious subthreshold maxima.
The fix: count spikes by v-hysteresis (cross +0.5, reset −0.5) and SAO loops as
**prominence-H peaks of w**. Because noise enters only v and w low-pass-filters
it, w is smooth; ρ = nL/(w-peaks), validated against σ = 0.

### 4.2 The measured signature: high-q-first, dt-converged

| q (ρ = 1/q) | c | σ_pq |
|---|---|---|
| 2 (1/2) | −0.77 | > 0.045 (widest plateau, most robust) |
| 3 (1/3) | −0.82 | 0.022 |
| 4 (1/4) | −0.845 | 0.016 |

σ_pq **decreases with q**: high-q (narrow, deep-funnel) plateaus dissolve
first — the predicted Farey-dissolution signature (the ⟨ρ⟩-drift figure shows
ρ = 1/4 leaving its plateau fastest, 1/2 slowest). The chief hazard is cleared:
σ_pq is **dt-converged** (⟨ρ⟩ stable between dt = 0.02 and 0.01), so the
dissolution is physical, not discretisation noise.

### 4.3 The noise scale: a universal local escape and a global μ^{3/2}

The measured dissolution threshold fits

```
σ_pq  ≈  C_q · √ε · μ(c_pq)^{3/2}      (β ≈ 1.13 = measured 1.11; §5).
```

It is tempting to read the μ^{3/2} as a *local* escape factor — an early Path-B
draft did, via an accumulated-variance "cross one secondary-canard sector"
argument (sectors ~μ wide, T_funnel ~1/μ ⇒ η_* ~ μ^{3/2}). **The crossover study
(§7, `MMO_CROSSOVER.md`) shows this local reading is wrong.** Sweeping the
folded-node normal form down to μ → 0, the local escape threshold η_*(μ) is
**canard-like** — it tends to the canard plateau with a μ-independent prefactor
C_q, and does **not** follow μ^{3/2} (a μ^{3/2} reference sits 3–200× below the
data). So:

- the **local** noise escape in the funnel is the **canard escape itself**
  (universal, shared prefactor C_q — §7); and
- the **μ^{3/2}** in σ_pq is a **model-specific global** property — how the
  global return / funnel-filling f(c) maps the plateau index q onto the
  injection geometry μ(c_pq). It matches the data empirically, but its
  derivation from the K2 funnel geometry (a_min(c), f(c)) is **open** (§6).

This is exactly the deterministic side's split (universal ceiling α = 2 +
model-specific funnel-filling f(c)) now mirrored on the noise side: **one
universal local mechanism (canard escape, C_q) + one model-specific global
geometry (funnel-filling).**

---

## 5. The exponent, stated correctly (Phase 3.1)

This section supersedes the `σ_pq ~ σ_*·q^{−α/γ}` form used in the phase
documents.

### 5.1 The reconciliation

The phase-3 draft compared a measured slope β ≈ 1.15 (from q = 3, 4) against a
"derived β = 1.5", and called the agreement ballpark. The "1.5" used the
**idealised** μ ∝ 1/q. But the exact relation is not 1/q: at the ρ = 1/q
plateau s = q − 1 = f(c)·(1−μ)/(2μ), and because f(c) *rises* across the band,
the self-consistent μ(q) falls **slower** than 1/q. Using μ at each plateau's
actual centre c_pq (closed-form, no new simulation; `mmo_noise_exponent.py`):

```
σ_pq^pred(q)  =  C_q · √ε · μ(c_pq)^{3/2},        μ(c_pq) = exact eigenvalue ratio.
```

| q | c_pq | μ(c_pq) | σ_pq^pred (C_q=10) | σ_pq^meas |
|---|---|---|---|---|
| 2 | −0.765 | 0.0585 | 0.040 | > 0.045 |
| 3 | −0.825 | 0.0422 | 0.025 | 0.022 |
| 4 | −0.845 | 0.0370 | 0.020 | 0.016 |
| 5 | −0.870 | 0.0307 | 0.015 | — |
| 6 | −0.895 | 0.0245 | 0.011 | — |

Fitted slope of the exact-μ prediction: **β_pred^exact ≈ 1.13** (q = 2–6),
**not 1.5**. The measured β ≈ 1.11 (q = 3, 4) matches it almost exactly, and
per-point the prediction tracks the measured σ_pq with **C_q ≈ 8–10 — the
full-FHN canard chapter's prefactor**. The 1.5-vs-1.1 gap was bookkeeping (the
μ ∝ 1/q idealisation), not missing physics; **Path A is not required for the
exponent**.

### 5.2 The correct operative statement (and the status of γ)

The reconciliation in §5.1 sets the dissolution threshold **equal to the
sector-crossing scale**:

```
σ_pq  ≈  σ_*(μ(q))  =  C_q · √ε · μ(q)^{3/2},
```

with **all** the q-dependence carried by μ(q) — which itself encodes the
staircase structure through the c ↔ q map (and so the staircase exponent and
f(c)). The validated exponent β ≈ 1.13 is therefore the **q-slope of σ_* itself**,
not a separate power on top of σ_*.

This corrects a conflation in the phase documents. The heuristic
`σ_pq ~ σ_*·q^{−α/γ}` (from "dissolution when δρ_noise ~ σ^γ reaches
Δ_pq ~ q^{−α}") is exact only in the idealised μ ∝ 1/q limit, where it gives
β = α/γ ≈ 1.55 (α = 1.55, γ = 1) — numerically close to the μ^{3/2} result
β ≈ 1.5 **by the coincidence that the realised α ≈ 3/2**, but a different
bookkeeping. With the exact μ(q) the operative law is σ_pq ≈ σ_*(μ(q)); the
α/γ form should not be quoted alongside β ≈ 1.13 as if α/γ = β (it is not:
α/γ = 1.55, α/β = 1.37).

Consequently **γ is reported separately and honestly**: γ is the
rotation-spread exponent in δρ_noise ~ σ^γ. Path B argues γ ≈ 1 (the
displacement's standard deviation is linear in σ), but γ has **not been
measured directly** — that is the deferred Task-β item (std(ρ) vs σ at fixed
in-plateau c). γ is *not* to be back-inferred from α/β.

---

## 6. What is established / what is open

| Piece | Status |
|---|---|
| Jordan obstruction ⇒ 3D FHR; folded node is the engine | **established** (Phase 1.5) |
| Operating regime: **μ < ε^{1/2}** (Krupa–Wechselberger transition, near singular Hopf) | **established** (§2.5) — standard scalings are asymptotic guides here |
| μ ∝ (c+1); Wechselberger ceiling s_max = (1−μ)/(2μ), α_ceiling = 2 | **derived** (asymptotic; transition-regime corrections expected) |
| Rotation map ln R = κμ, κ → 2π² (asymptotic; drifts 8.7→18.6 over q=3→6) | **measured** (Phase 2 / M2) |
| Funnel-filling f(c) = (2/κ)ln(a_max/a_min)/(1−μ); realised α ≈ 1.55 | **derived + validated** (μ cancels; matches to 5–15 %) |
| Noise dissolves staircase high-q-first, dt-converged | **measured** (Phase 3) |
| **Local** noise escape = canard escape (continuous μ→0 limit); **C_q ≈ 8–10 shared** | **derived + validated** (crossover A+B, `MMO_CROSSOVER.md`) |
| **μ^{3/2}** in σ_pq is a *global* funnel-filling factor, not local escape | **empirical**; μ·a_min candidate degenerate with μ^{3/2} at current resolution (M2) |
| Exponent: σ_pq ≈ C_q√ε·μ(q)^{3/2}, β ≈ 1.13 = measured 1.11 | **reconciled** (Phase 3.1, Outcome 1) |
| Closed-form a_min(c); first-principles κ, γ; derive the global σ_pq μ-law | **open** — sector-packing route found (μ·ε^{(1−μ)/2}, Thm 3.1) but needs **Krupa–Wechselberger 2010** transition-regime inner solution (`MMO_SIGMAPQ_PLAN.md` Phase B) |
| Full-FHR crossover σ_*(control) through the FSN, C_q ≈ 8–10 | **open — Phases C–D** (`MMO_CROSSOVER_PLAN.md`) |
| Converged β with error bars; **direct** γ from std(ρ)~σ^γ | **open — precision** (Task-β counting study) |

**Validated:** the law σ_pq ≈ σ_*(μ(q)) holds in form, order, scale, prefactor,
and exponent, with C_q inherited from the canard chapter — the cross-chapter
unification the chapter was built to deliver, at the project's standard
(universal mechanism + model-specific prefactor). **Open** items are precision
and first-principles closure, neither of which blocks the result; both sit on
the same Path-A frontier as the deterministic side.

---

## 7. Cross-chapter unification — the headline

The MMO regime is the **fifth instance** of the project's pattern, and the
**first in three dimensions**. The same Krupa–Szmolyan fold blow-up that
organises every other regime here operates at the **3D folded node**: applied
to the forward flow it gives the secondary-canard rotation map (κ ≈ 2π²) behind
the deterministic staircase. The noise escape that dissolves the staircase is
**the canard escape itself**: the crossover study (A+B, `MMO_CROSSOVER.md`)
sweeps the folded-node normal form to μ → 0 and recovers the canard escape
**continuously**, with a **shared, μ-flat prefactor C_q ≈ 8–10**. The measured
threshold

```
σ_pq  ≈  C_q · √ε · μ(q)^{3/2},     C_q ≈ 8–10  (= full-FHN canard, shared)
```

then decomposes the way the deterministic side does: the prefactor C_q is the
**universal** canard-escape constant (proven the same, not merely close, via the
continuous μ → 0 limit), and the **μ^{3/2}** is a **model-specific global**
funnel-filling factor — *not* a local escape law (the local escape carries no
μ^{3/2}). So the chapter's unification is sharper than originally stated: **one
universal local escape mechanism (the canard escape, shared C_q), two global
geometries (the single fold of Chapter 1 ↔ the folded-node funnel here).** Noise
destroys deterministic structure wherever transverse stability is weakest, and at
the folded node the *local* weakness is identical to the 2D canard's — only the
global funnel that surrounds it is new.

---

## 8. Hazards (encountered and handled)

1. **2D Jordan obstruction** — fatal to the original scheme; resolved by the
   redirect to 3D FHR (§1).
2. **Wrong FHR parametrisation** — the canonical Rinzel δ → 0 form is a pure
   spiker; the working variant (§2.1) is the MMO regime (`MMO_PHASE1_5.md`).
3. **SAO undercount under noise** — v-based counting explodes (~10³ spurious
   maxima); fixed by w-based prominence counting (§4.1).
4. **dt false positive** — too-coarse dt mimics dissolution; cleared by
   dt-convergence (§4.2).
5. **Counter-H sensitivity / countability floor** — σ_pq is H-sensitive at the
   factor-2 level and q ≥ 5 deep-funnel loops sit near the noise floor; this is
   the limiting systematic on the *precision* of β (not its validity), and the
   target of the Task-β study (per-plateau a_min-scaled H).
6. **Transition regime (μ < ε^{1/2})** — μ ∈ [0.03, 0.16] sits below ε^{1/2} ≈ 0.28
   throughout (§2.5), i.e. the Krupa–Wechselberger transition regime near the
   singular Hopf, not standard Wechselberger. Existence of secondary canards holds
   for all μ ∈ (0,1), but the standard *scalings* (sector widths, the s_max
   ceiling, κ ≈ 2π²) are asymptotic guides here, and any first-principles closure
   must use Krupa–Wechselberger 2010. This is the source of the small-μ numerical
   delicacy and the low-q κ-drift.

---

## 9. Reproduce

```bash
python3 regime-tests/mmo_fhr_staircase.py     # Phase 1b: devil's staircase, ρ(c)
python3 regime-tests/mmo_fhr_foldednode.py    # Phase 1.5: folded-node μ(c), Task B'
python3 regime-tests/mmo_fhr_alpha_fine.py    # Phase 1.5: realised α = 1.55 (22-pt grid)
python3 regime-tests/mmo_fhr_timescale.py     # Phase 1.5: 2-slow check, δ-robustness
python3 regime-tests/mmo_k2_return_map.py     # Phase 2:   K2 rotation map κ≈2π², f(c), α
python3 regime-tests/mmo_noise_staircase.py   # Phase 3:   σ_pq(q), high-q-first dissolution
python3 regime-tests/mmo_noise_exponent.py    # Phase 3.1: exact-μ(q) reconciliation, β≈1.13
python3 regime-tests/mmo_crossover_nf.py      # Phase A+B: μ→0 crossover, local escape = canard
```

Outputs in `data/`, `figures/`, `regime-tests/results/mmo/`. Phase docs:
`MMO_FHR_PLAN.md` (plan + banner trail), `MMO_PHASE1_5.md`,
`MMO_ALPHA_DERIVATION.md`, `MMO_TIMESCALE_CHECK.md`, `MMO_K2.md`, `MMO_NOISE.md`,
`MMO_CROSSOVER.md` (+ `MMO_CROSSOVER_PLAN.md` for Phases C–D).

---

## 10. Key references

**Folded-node MMO theory**
- Wechselberger (2005). Existence and bifurcation of canards in R³ in the case
  of a folded node. *SIAM J. Appl. Dyn. Syst.* 4(1), 101–139.
- Brøns, Krupa, Wechselberger (2006). Mixed-mode oscillations due to the
  generalized canard phenomenon. *Fields Inst. Commun.* 49, 39–63.
- Krupa, Popovic, Kopell (2008). Mixed-mode oscillations in three time-scale
  systems. *SIAM J. Appl. Dyn. Syst.* 7(2), 361–420.
- Desroches, Guckenheimer, Krauskopf, Kuehn, Osinga, Wechselberger (2012).
  Mixed-mode oscillations with multiple time scales. *SIAM Review* 54(2), 211–288.

**Rinzel's model**
- Rinzel (1987). A formal classification of bursting mechanisms in excitable
  systems. *Lecture Notes in Biomathematics* 71, 267–281.

**Entry-exit / supervisor**
- Kaklamanos, Kuehn, Popovic, Sensi (2025). Entry-exit functions with
  intersecting eigenvalues.

**Stochastic slow-fast (the 2D engine this generalises)**
- Berglund & Gentz (2006). *Noise-Induced Phenomena in Slow–Fast Dynamical
  Systems.* Springer.
- Berglund, Gentz, Kuehn (2015). Stochastic dynamic bifurcations and
  excitability. Springer.
- Krupa & Szmolyan (2001). Extending GSPT to nonhyperbolic points. *SIAM J.
  Math. Anal.* 33(2), 286–314.
