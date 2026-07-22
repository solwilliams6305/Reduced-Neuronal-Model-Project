# MMO chapter — 3D FitzHugh–Rinzel plan

**Status:** active plan, supersedes the 2D version in `MMO_PLAN.md` (which
failed at the Jordan-curve obstruction caught in Phase 1; see banner there
and `MMO_PHASE1.md` for the diagnostic numerics). This document scopes the
chapter using the corrected model, with the same noise question but the
right phase-space dimension.

> **⚠ Phase 1.5 executed — see `MMO_PHASE1_5.md`.** Three updates: (1) **the §0
> canonical FHR (Rinzel form) does NOT produce MMOs** — it is a pure spiker at
> every δ∈[0.005,0.04] tested (folded node exists but the global return bypasses
> its funnel). The working MMO model is the Phase-1 variant `v'=v−v³/3−w+y+I,
> y'=εδ(c−v), δ=0.2, c≈−0.7…−0.96`; §0's equations should be replaced with it.
> (2) **Folded-node μ is verified as the engine** (μ(c) computed from the
> desingularised reduced flow; observed SAO counts grow with μ→0), **but
> Wechselberger's `(1−μ)/(2μ)` is the SAO _upper bound_, not an equality** — the
> §2.1 "≈" should read "≤". (3) **α ≈ 1.0** (finer sweep, q=2–5), **below** the
> circle-map [2,3] — consistent with folded-node (secondary-canard) organisation
> rather than a critical circle map.
>
> **(3-update — α derivation attempted, `MMO_ALPHA_DERIVATION.md`): α ≠ 1.** The
> α≈1.0 above was a coarse-grid artifact (retracted). From the derived `μ ∝ (c+1)`
> + Wechselberger's `s_max=(1−μ)/(2μ)`, the folded-node *ceiling* is **α = 2**; the
> *realised* staircase is **α ≈ 1.6** (`s_obs~(c+1)^{−1.74}`), reduced below 2 by
> the global-return funnel-filling f(c). Phase 2's K2/entry-exit return map must
> supply f(c) to close the realised α; the noise corollary becomes
> `σ_pq ~ σ_*·q^{−α/γ}` (the plan's `σ_*/√q` is superseded).
>
> **(4) Phase-2 prerequisite cleared (`MMO_TIMESCALE_CHECK.md`): Wechselberger
> applies.** δ=O(1) (0.2 and 1.0) is a genuine 1-fast/2-slow system (both w,y are
> O(ε)); the folded node is real ∀δ∈[0.05,1.0]; MMOs exist at the cleanest 2-slow
> limit δ=1.0. **Launch Phase 2 K2-chart at δ=1.0**, NOT at the plan's §0 δ→0
> (that is the three-timescale / Krupa–Popovic–Kopell regime and is unnecessary).
> Caveat: μ small (≈0.03–0.16) ⇒ many-secondary-canard regime, K2 matching needs
> care (hazard §6.4).
>
> **(5) α tightened (`mmo_fhr_alpha_fine.py`, 22-pt grid): realised α = 1.55,
> bracket [1.48, 1.60]** (deep-subrange headline 1.58±0.06; ±0.02 within-fit,
> ±0.06 systematic across fit choices). Both Phase-1.5 prerequisites are now
> cleared. **K2 target set:** does the K2/entry-exit f(c) reproduce α≈1.55? A
> result of 1.55 is a hit, 1.3 or 1.8 a miss — precise enough to separate
> "matching right, numerics off" from "matching wrong" in the small-μ regime.
> Trajectory to close the chapter: K2 return map → derive f(c)/α → FW action for
> γ → σ_pq ~ σ_*·q^{−α/γ} → write-up.
>
> **(6) Phase 2 K2 calculation done (Path B, `MMO_K2.md`): mechanism CONFIRMED,
> α reproduced ≈1.45–1.49** (target 1.55±0.06 — confirmed at the lower edge; the
> ~0.05–0.10 gap is within the flagged SAO-count dt-sensitivity). Key result: the
> folded-node rotation map is `ln(R)=κμ` with **κ≈2π²** (=19.7; measured 19.6±1.8,
> constant across c) — the per-turn SAO growth rate set by μ. This gives the
> funnel-filling with μ **cancelling**: `f(c)=(2/κ)·ln(a_max/a_min)/(1−μ)`, i.e. the
> global-return amplitude span / κ (matches measured f to ~5–15%). Universal piece
> = folded-node rotation (κ≈2π²) + ceiling; model-specific piece = numerical global
> return a_min(c). Open (Path A): closed-form a_min(c) and a first-principles κ=2π².
> Deterministic side closes; next = σ_pq (γ via FW action) → noise staircase → write-up.
>
> **(7) Phase 3 noise dissolution done (Path B, `MMO_NOISE.md`): mechanism
> VALIDATED, exponent ballpark.** Degenerate noise on v dissolves the staircase
> high-q-first (σ_pq: q=2 >0.045, q=3 ≈0.022, q=4 ≈0.016; dt-converged — physical,
> not discretisation). `σ_pq ~ σ_*·q^{−α/γ}` holds in form/order/prefactor: 3B
> (accumulated Brownian variance + measured κ≈2π²) gives `σ_*=C_q√ε·F(μ)`,
> **F(μ)~μ^{3/2}**, γ≈1, β≈1.5; prefactor needs **C_q≈10 — same as the full-FHN
> canard chapter** (key cross-chapter check). Measured β≈1.15, γ≈1.35 — ballpark
> consistent. NOT precision-tested: SAO-counting-under-noise (the flagged limiting
> hazard) → factor-2 counter-H sensitivity, only q=3,4 cleanly resolved. Open
> (Path A): first-principles γ/F(μ); finer-counting β. Remaining: `MMO_CHAPTER.md`
> + README MMO row → "derived + validated (mechanism)".
>
> **(8) Phase 3.1 exponent tightening done (Task α, `mmo_noise_exponent.py`):
> OUTCOME 1 — the gap was bookkeeping.** The β=1.5 used the idealised μ∝1/q; the
> exact μ(c_pq) (closed form, no new sim) gives **β_pred^exact ≈ 1.13** (q=2–6),
> matching the **measured β ≈ 1.11** — slope AND prefactor (C_q≈8–10) reconciled.
> So `σ_pq ~ σ_*·q^{−α/γ}` is validated in form/order/prefactor/exponent. **Path A
> NOT needed for the exponent.** Open only as precision (Task β: converged counting
> at q=2,5,6 + direct γ → β∈[1.0,1.15]) and first-principles γ/F(μ) (the K2 Path-A
> frontier, shared with the deterministic side). Exponent reconciled.
>
> **(9) Crossover A+B done (normal form, `MMO_CROSSOVER.md`): keystone inheritance
> PROVEN as a continuous limit + a refinement.** The folded-node K2 SDE
> `dW/dT=μV−(1+μ)/2` → canard `dW/dT=−1/2` (λ=1/2) as μ→0; measuring the escape
> with the canard criterion verbatim, **η_*(μ) → the canard plateau** as μ→0 and
> **C_q=Θ_crit is flat ≈2.0 for μ≲0.1** (rising mildly to 2.5 by μ=0.5, onset of
> funnel near μ_c≈0.6). So MMO and canard share ONE local escape mechanism + C_q,
> as a continuous limit (not a coincidence). **Refinement:** η_* is canard-like
> (constant-drift scale) throughout — **no μ^{3/2} in the local normal-form escape**;
> the μ^{3/2} is the GLOBAL-return funnel-filling (MMO_K2 f(c)). Same
> universal-local + model-specific-global split as the deterministic side. Absolute
> C_q (2.0 here vs 2.8 canard map) is dt/settings-dependent; constancy + μ→0 limit
> are robust. Open: full-FHR Phases C–D (FSN scaling, C_q≈8–10) — separate handoff.

---

## 0. Why FHR, not FHN

2D autonomous flows have Jordan-curve attractors: a stable limit cycle is
a simple closed curve, so it cannot carry a deterministic L^a S^b symbol
sequence at fixed parameters. MMOs in autonomous slow-fast systems need
**at least one extra slow direction** — a second slow variable whose drift
through a fold creates a **folded singularity** (folded node, folded
saddle, folded saddle-node) at which trajectories acquire a rotation
number.

The minimal carrier model in the FHN family is FitzHugh–Rinzel (Rinzel
1987), which adds a single slow modulator y to the FHN equations:

```
dv/dt  =  v − v³/3 − w − y + I,             (fast)
dw/dt  =  ε ( v + a − b w ),                 (slow)
dy/dt  =  δ ( c v + d − y ).                 (slow modulator)
```

with two small parameters ε, δ and a third slow timescale set by δ.
Typical parameter values for clean MMOs (from Rinzel 1987 / Desroches
et al. 2012): a ≈ 0.7, b ≈ 0.8, ε ≈ 0.08, δ ≈ 0.001–0.01, c ≈ 1.0,
d ≈ 0.8, with I as the bifurcation control. The slow modulator y dragges
the effective excitability slowly through the canard region, producing
sequences of L spikes interleaved with S subthreshold loops.

**This keeps the chapter intellectually continuous with the canard
chapter:** the v–w subsystem is identical to autonomous FHN, the
Krupa–Szmolyan blow-up at v = −1 still applies, and y appears as a
second slow drift that turns the 2D fold into a **folded node** in 3D.

---

## 1. Question

For FHR in the MMO regime (parameter window where the deterministic
flow produces L^a S^b patterns):

1. **Deterministic structure.** Predict the L^a S^b pattern at each (I, ε, δ)
   in the MMO band. Equivalently, predict the rotation number ρ = L/(L+S)
   and the mode-locked plateau widths Δ_pq.

2. **Noise dissolution.** With degenerate noise σ on v (same noise
   convention as the FHN chapters), characterise the dissolution of the
   Farey staircase: P(L^a S^b | I, ε, δ, σ), order of dissolution as σ
   grows, q-dependent threshold σ_pq.

3. **Cross-chapter unification.** Connect σ_pq to the canard chapter's
   σ_* = C_q · √ε · λ^{1/2}. The conjecture is σ_pq ~ σ_*(λ_folded_node) · q^{−β}
   where λ_folded_node is the effective slow-drift parameter at the
   folded singularity and β is set by the noise's propagation through the
   3D return map.

---

## 2. Analytical framework — folded-node geometry

### 2.1 The folded node

In 3D autonomous slow-fast systems, the slow manifold is 2-dimensional
and generically has fold curves rather than fold points. At a generic
point on a fold curve the dynamics is 2D fold-like (canard chapter
machinery applies), but at **isolated points** on the fold curve the
slow flow develops a critical structure — a **folded singularity** —
classified by Wechselberger (2005) as folded saddle, folded node, or
folded focus depending on the eigenvalues of the reduced (desingularised)
slow flow.

The **folded node** is the MMO-producing case. It has two real negative
eigenvalues with ratio μ = λ_weak / λ_strong ∈ (0, 1). The eigenvalue
ratio μ is the **single most important quantity** in the chapter — it
controls the rotation number of canard trajectories passing through the
folded-node neighbourhood:

```
Number of small oscillations per spike  ≈  (1 − μ) / (2μ)    (Wechselberger 2005, Thm 3.1)
```

So the deterministic L^a S^b structure is set by μ(I, δ), which is in
turn a known function of the FHR parameters via the linearisation at the
folded singularity.

### 2.2 The chapter's analytical engine

Wechselberger's folded-node theorem provides:

- **A 3D extension of the Krupa–Szmolyan blow-up** (Wechselberger 2005,
  §4). The blow-up has multiple charts (the "K1, K2, K3" charts) covering
  different scaling regimes around the folded singularity. The canard
  chapter's 2D K-S blow-up is the special case of the K2 chart for the
  folded saddle-node limit.
- **2k+1 secondary canards** for μ ∈ (1/(2k+3), 1/(2k+1)). These
  organise the L^a S^b patterns: trajectories that pass through different
  numbers of small loops are separated by canard surfaces in 3D phase
  space.
- **An explicit return map** at the folded node, derived from matching
  the K2 chart's local solution to the global flow. The map is a 1D map
  on a Poincaré section transverse to the cycle, with rotation number
  determined by μ.

This is exactly the supervisor-aligned machinery: **Wechselberger's
folded-node theory + Krupa–Popovic–Kopell's 3D-MMO framework + Nikola's
entry-exit functions** combine into a single analytical engine. The
chapter doesn't need to invent new mathematics — it needs to apply the
existing canon to FHR specifically and then layer the noise question on
top.

### 2.3 Connection to the canard chapter's σ_*

The canard chapter's σ_* = C_q · √ε · λ^{1/2} was derived for the K-S 2D
blow-up. In 3D, the analogous calculation lives in the K2 chart of
Wechselberger's blow-up, with the slow-drift parameter λ replaced by the
folded-node's μ-dependent local drift. The noise scaling should still be

```
σ_* (3D)  =  C_q · √ε · F(μ)
```

with F(μ) a folded-node-geometry function that reduces to λ^{1/2} in the
2D limit. **Computing F(μ) explicitly is the natural cross-chapter
analytical deliverable** — it extends σ_*(λ) from 2D fold to 3D folded
node and would be a clean Section in the MMO chapter that directly
parallels the canard chapter's σ_* derivation.

---

## 3. Noise framework

### 3.1 The dissolution mechanism is unchanged

The Farey-staircase dissolution argument from the original MMO_PLAN.md
§3 carries over without modification: plateaus of width Δ_pq dissolve
when noise-induced rotation-number fluctuation δρ_noise(σ) ~ Δ_pq.
High-q plateaus dissolve first, the dissolution sequence is the same as
before, the σ_pq ~ σ_* · q^{−β} conjecture stands.

What changes is the σ_* scale that sets the dissolution threshold:

- 2D canard chapter: σ_*(λ) = C_q √ε λ^{1/2}
- 3D MMO chapter: σ_*(μ, ε, δ) = C_q √ε F(μ, δ)

with the folded-node function F(μ, δ) replacing λ^{1/2} as the
geometric factor.

### 3.2 What the MMO chapter measures

- P(L^a S^b | I, ε, δ, σ): joint distribution of MMO patterns at given
  parameters under noise.
- σ_pq: threshold at which plateau (p, q) dissolves (each plateau a
  separate observable, indexed by p/q).
- Dissolution order: empirical sequence of plateau dissolution as σ
  grows. Compared to the prediction (high-q first).
- σ_pq / σ_* collapse: testing the cross-chapter unification.

---

## 4. Task sequence (3D version)

### Phase 1.5 — Deterministic FHR foundations (1–2 weeks)

Sonnet's Phase 1 has already started this with `mmo_fhr_staircase.py`,
showing a clean ρ = L/(L+S) staircase locking at 1, 2/3, 1/2, 1/3, …
The remaining Phase 1.5 work:

**Task A':** Finer sweep across the MMO band at fixed ε, δ. Resolve the
Farey mediants (L^2 S^3, L^3 S^2 between L^1 S^1 and L^1 S^2 / L^2 S^1)
by using ~200 c-grid points instead of the current ~50. Confirm Δ_pq
scaling on at least 5 denominators (currently 4).

**Task B':** Compute the folded-node μ as a function of (I, δ) at the
FHR canonical operating point. Verify Wechselberger's prediction
"#small oscillations ≈ (1-μ)/(2μ)" against the measured ρ in the
staircase.

**Task C':** Honest α-fit for Δ_pq ~ q^{−α}. Sonnet's preliminary
α ≈ 1.1 is unreliable on the coarse grid; the finer sweep should
either confirm α ∈ [2, 3] (standard circle-map) or surface the
canard-modified value. Either is a substantive result.

### Phase 2 — Folded-node blow-up + return map (3–4 weeks)

This is the rate-limiting analytical phase. Pre-meeting prep with
Nikola on the right approach is essential here.

**Task D':** Apply Wechselberger's 3D blow-up at the FHR folded node.
Compute the inner solution in the K2 chart. Match to the global slow
flow on both sides.

**Task E':** Derive the explicit 1D return map R from the matched
asymptotic. R is parameterised by μ and the global slow flow. Predict
Δ_pq from R's mode-locked windows.

**Task F':** Connect to entry-exit. Nikola's framework gives the right
way to describe the trajectory's passage through the folded node;
should reduce to or extend the Wechselberger return map.

### Phase 3 — Noise dissolution (2–3 weeks)

**Task G':** σ_*(μ, ε, δ) derivation. Apply the canard chapter's
accumulated-Brownian-variance argument inside the K2 chart of the
folded-node blow-up. Get F(μ, δ) explicitly. This is the chapter's
cross-chapter unification deliverable.

**Task H':** Numerical σ_pq sweep. For three plateaus (L^1 S^1, L^1 S^2,
L^2 S^1) measure the σ at which each dissolves. Compare with the
F(μ, δ) prediction.

**Task I':** Devil's-staircase noise broadening. Same I-sweep as Task A'
but with σ > 0; show plateau dissolution as σ increases.

### Phase 4 — Write-up (1–2 weeks)

**Task J':** `MMO_FHR_CHAPTER.md` with derivation + numerical
validation. Structure parallel to `CANARD_BLOWUP.md` and
`TONIC_PHASE.md`.

**Task K':** Integrate MMO row into master README §1 table:
`σ_pq ~ σ_*(μ) · q^{−β}` with status updated from "planned (PhD scope)"
to "derived + validated."

**Total: 7–11 weeks** (about a week longer than the 2D plan since the
3D blow-up is more involved).

---

## 5. Pre-meeting questions for Nikola

The 2D-to-3D redirect is a substantive plan change and worth flagging at
the supervisor meeting. Specific questions:

1. **Wechselberger vs entry-exit as the analytical engine.** Both
   describe the folded-node passage; do they agree, or do they give
   complementary information? Is there a place where entry-exit gives
   sharper predictions than Wechselberger's 2005 framework?

2. **Cross-chapter unification.** The conjecture σ_* (3D) =
   C_q √ε F(μ, δ) reducing to the 2D σ_* in the appropriate limit —
   is this the cleanest framing, or should the MMO chapter break with
   the canard chapter's prefactor convention and use a μ-based
   parameterisation throughout?

3. **FHR vs more general 3D systems.** The MMO chapter as planned is
   specifically for FHR. Should it instead be framed as a general
   folded-node MMO chapter with FHR as a worked example, parallel to
   how the canard chapter has FHN + Van der Pol? That would set up
   Morris–Lecar with a second slow variable (or similar) as the
   universality test for this chapter.

4. **Publication target.** Wechselberger's group publishes folded-node
   MMO work in *SIAM J. Appl. Dyn. Syst.* The noise extension is the
   chapter's distinctive contribution. Is SIADS the right target or
   would *Journal of Nonlinear Science* be a better fit?

---

## 6. Hazards specific to 3D

1. **The folded-node blow-up has multiple charts.** Wechselberger's 2005
   paper carries 3 coordinate patches (K1, K2, K3) covering different
   scaling regimes, with matching between them. This is more involved
   than the 2D K-S blow-up's single chart. Expect Phase 2 to take
   longer than the analogous canard-chapter derivation.

2. **Numerical FHR is more expensive.** Three variables and the y-drift
   timescale δ ≪ ε mean long simulation times to resolve the slow
   manifold structure. The σ-sweeps in Phase 3 will need careful
   timestep / integrator choice to avoid spurious dissolution from
   numerical noise.

3. **Multiple slow timescales** (ε from w, δ from y) create new
   scaling regimes that don't appear in the canard chapter. The right
   ε/δ ratio for clean MMOs is ε ≫ δ; in the singular limit ε → 0
   first, then δ → 0, but at finite values the two-scale geometry has
   to be tracked.

4. **The Wechselberger framework assumes generic folded nodes.** FHR's
   specific algebraic structure might produce non-generic folded
   singularities (folded saddle-nodes, in particular) where the
   standard 2k+1-canards theorem doesn't apply directly. Check for
   non-generic structure before applying the theorem blindly.

---

## 7. Why this is the right way to go

Three reasons:

1. **Mathematical correctness.** 2D autonomous flows can't carry MMO
   structure (Jordan obstruction). The chapter has to live in 3D or
   higher. FHR is the minimal carrier in the FHN family.

2. **Continuity with existing chapters.** FHR contains FHN as a
   subsystem (set y = 0 and ignore the y dynamics). The canard
   chapter's σ_* derivation extends to FHR's folded-node case in
   the K2 chart of Wechselberger's blow-up. The methodological
   parallel that the existing tonic chapter set up
   (Krupa–Szmolyan applied to adjoint Floquet) generalises here as
   well.

3. **Supervisor alignment.** Wechselberger's folded-node theory is
   directly continuous with Nikola's entry-exit work and his PhD
   advisor lineage (Wechselberger). The MMO chapter in 3D FHR sits
   inside Nikola's research line in a way the 2D version did not.
   Co-authorship and publication strategy align cleanly.

---

## 8. References

**Folded-node MMO theory:**
- Wechselberger (2005). Existence and bifurcation of canards in R³ in
  the case of a folded node. *SIAM J. Appl. Dyn. Syst.* 4(1), 101–139.
- Brøns, Krupa, Wechselberger (2006). Mixed mode oscillations due to
  the generalized canard phenomenon. *Fields Inst. Commun.* 49, 39–63.
- Krupa, Popovic, Kopell (2008). Mixed-mode oscillations in three
  time-scale systems: a prototypical example. *SIAM J. Appl. Dyn. Syst.*
  7(2), 361–420.
- Desroches, Guckenheimer, Krauskopf, Kuehn, Osinga, Wechselberger
  (2012). Mixed-mode oscillations with multiple time scales.
  *SIAM Review* 54(2), 211–288.

**Rinzel's model:**
- Rinzel (1987). A formal classification of bursting mechanisms in
  excitable systems. *Lecture Notes in Biomathematics* 71, 267–281.

**Entry-exit / supervisor:**
- Kaklamanos, Kuehn, Popovic, Sensi (2025). Entry-exit functions with
  intersecting eigenvalues.

**Stochastic 3D slow-fast (less developed than the 2D theory):**
- Berglund, Gentz, Kuehn (2015). Stochastic dynamic bifurcations and
  excitability. (Chapter in *Stochastic Processes, Multiscale Modeling,
  and Numerical Methods for Computational Cellular Biology*, Springer.)

---

## 9. Status

**Planning document.** Phase 1.5 work is the immediate next step
(finer FHR sweep, Wechselberger-μ verification, honest α fit).
The chapter is otherwise queued behind the supervisor meeting and the
canard + Van der Pol writeup. Reread before starting Phase 1.5 to
verify the FHR parameter choices and the operational definitions in
§4.

The 2D MMO_PLAN.md document is retained as historical record + Phase 1
diagnostic context. The 2D banner there points here.
