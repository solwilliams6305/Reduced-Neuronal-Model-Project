# MMO crossover — the folded-node noise law as the canard law in disguise

**Goal.** Establish, analytically and numerically, that the 3D folded-node
noise threshold `σ_*(μ) = C_q·√ε·μ^{3/2}` and the 2D canard escape threshold
`σ_* = C_q·√ε·λ^{1/2}` are two regimes of a **single** crossover function, with
the **same** prefactor C_q, connected as the folded node degenerates to a
folded saddle-node (μ → 0). This converts the C_q "inheritance" from a measured
coincidence into a proven continuous limit, and — because it is literally the
bridge to the canard chapter — supplies the unification pillar that turns the
MMO write-up into a self-justifying PhD chapter.

**Why this is the right pillar.** It pre-empts the two sharpest examiner
objections at once ("is C_q *really* the same constant?" and "is the μ^{3/2}
law connected to anything, or just fitted?"), and it does double duty as the
Chapter-1 → Chapter-2 connective tissue. It is also tractable: the analytical
core is a clean normal-form statement, and the numerics reuse both chapters'
existing simulators.

---

## 1. The claim

The canard chapter (Ch.1) and the MMO chapter (Ch.2) study the **same local
object** — a fold/canard threaded by degenerate noise — at two ends of one
family of folded singularities:

```
folded saddle  ——  folded saddle-node (FSN)  ——  folded node
(single canard,        (μ → 0, the                (funnel of s_max
 Ch.1's 2D passage)     crossover)                 secondary canards, Ch.2)
```

Both noise laws have the form `σ_* = C_q·√ε·[geometry]`, with C_q the *local*
escape constant (how many canard-amplitudes of accumulated Brownian noise it
takes to leave the canard) and `[geometry]` the regime-specific factor
(`λ^{1/2}` for the single passage, `μ^{3/2}` for the funnel). **The claim:** C_q
is one constant across the whole family; only the geometric factor changes
form, and it does so continuously through the FSN.

---

## 2. The structural spine — the normal forms meet at μ → 0

This is the fact that makes the crossover rigorous rather than rhetorical.

**Canard K-S normal form (Ch.1):**
```
dV/dT = V² − W,        dW/dT = −λ.
```
Constant slow drift −λ; no rotation; a single fold passage. Accumulated
variance over the window (V_* = λ^{1/3}, T_window = λ^{−1/3}) gives
η_* = C_q·λ^{1/2}, i.e. σ_* = C_q·√ε·λ^{1/2} (CANARD_BLOWUP.md §4).

**Folded-node K2 normal form (Ch.2):**
```
dV/dT = V² − W,        dW/dT = μV − (1+μ)/2.
```
(Wechselberger's K2 chart; convention to be verified against Wechselberger 2005
§4 — see Hazard H1.) The **μV term is what generates the funnel rotation** and
the s_max = (1−μ)/(2μ) secondary canards.

**The bridge.** Taking μ → 0 in the folded-node normal form sends
`dW/dT = μV − (1+μ)/2 → −1/2`: the V-dependent (rotation) term vanishes and the
drift becomes **constant**, recovering the canard normal form with λ = 1/2.
So the canard passage is the **μ = 0 face of the same normal form**, the noise
enters identically as η = σ/√ε in both, and the local escape constant C_q is
therefore the *same object computed in the same coordinates*. The two chapters
are not analogous — they are one calculation evaluated at two values of μ.

**The crossover question (the analytical deliverable).** As μ grows from 0, the
rotation term μV competes with the constant drift. There is a crossover scale
μ_c at which rotation takes over: below it the dynamics is canard-like and
η_* tracks the constant-drift `λ^{1/2}` regime; above it the funnel forms and
η_* tracks the `μ^{3/2}` regime. The deliverable is the single crossover
function η_*(μ; λ) interpolating the two power laws, and the scale μ_c where the
geometry switches.

---

## 3. Task sequence

### Phase A — Normal-form crossover (analytical core)

**A1.** Confirm the folded-node K2 normal form and its μ → 0 reduction to the
canard K-S form (against Wechselberger 2005; do not reconstruct the
parabolic-cylinder inner solution — only the normal form and its limit are
needed). State the constant-drift limit cleanly (λ_eff = 1/2 in NF units).

**A2.** Identify the crossover scale μ_c by balancing the rotation term μV
against the constant drift over the canard window: μ·V_* ~ O(drift). With
V_* ~ λ^{1/3} this gives μ_c as a function of the (normal-form) drift. Below
μ_c: single-passage, η_* ~ λ^{1/2}; above: funnel, η_* ~ μ^{3/2}.

**A3.** Write the unified geometric factor F(μ, λ) with
F → λ^{1/2} (μ ≪ μ_c) and F → μ^{3/2} (μ ≫ μ_c). This is the matched-asymptotic
crossover function; C_q is the (μ-independent) prefactor multiplying it.

*Deliverable:* a one- to two-page derivation. The clean, defensible result even
if the full F(μ,λ) is hard is **A1** (the normal forms coincide at μ → 0) plus
**A2** (the crossover scale) — that alone proves continuity and shared C_q.

### Phase B — Normal-form SDE numerics (validate the crossover, cheaply)

**B1.** Integrate the folded-node normal-form SDE
`dV = (V² − W)dT + η dB, dW = (μV − (1+μ)/2)dT`, sweeping μ from O(1) down
toward 0 (e.g. μ ∈ {0.5, 0.2, 0.1, 0.05, 0.02, 0.01}). For each μ measure the
escape threshold η_* by the **same operational criterion** the canard chapter
used (noise drives V off the canard by V_*; CANARD_BLOWUP.md gives Θ_crit ≈ 2.8
for V_cross = 1).

**B2.** Plot η_*(μ). Show it traces the predicted crossover: the `μ^{3/2}` slope
at larger μ bending to a μ-independent plateau (the canard `λ^{1/2}` value with
λ = 1/2) as μ → 0, with **C_q ≈ 2.8 constant throughout** — the same
normal-form C_q as Ch.1.

*This is the lowest-risk confirmation and the figure that sells the result:* one
curve, one constant, two regimes. Reuses the canard normal-form simulator
(`canard_normal_form_map.py`).

### Phase C — Full FHR: locate the transition and the μ ↔ λ relation

**C1.** In FHR, map the folded singularity across the control parameter:
compute the desingularised-flow eigenvalues (Jacobian
`J = [[1+δ, −b],[2bδ(c+1), 0]]`) and locate where μ → 0 (the folded
saddle-node), the folded-node band (μ ∈ (0,1)), and any folded-saddle/
folded-focus boundaries. Identify the control that drives μ → 0 (c, and/or the
y-coupling).

**C2.** Establish the **μ ↔ λ relation** in the FSN scaling regime — how the
folded-node eigenvalue ratio μ and the effective 2D slow-drift λ co-vary as the
singularity is approached. This is the crux that lets the two *physical* laws
(`μ^{3/2}`, `λ^{1/2}`) be read as two regimes of one σ_*(control). Expect to
need the folded-saddle-node scaling theory (Krupa–Wechselberger 2010); if a
closed form is out of reach, extract μ(control) and λ(control) numerically and
form the relation empirically.

### Phase D — Full FHR SDE across the transition (the headline figure)

**D1.** Define ONE noise observable that specialises correctly on both sides:
the threshold at which noise displaces the trajectory off the (secondary-)canard
by its local amplitude. In the folded-node band this is the staircase σ_pq
(`mmo_noise_staircase.py`); approaching the FSN it should continuously become
the canard-escape σ_* (`canard_escape_autonomous.py`). Use the common
"off-canard by V_*" definition so the two are the same measurement.

**D2.** Sweep the control through the transition with degenerate noise on v,
measure σ_*(control), and show it is **continuous through the FSN** with
**C_q ≈ 8–10 throughout** (the full-FHN/FHR value, the finite-ε counterpart of
the normal-form 2.8). Overlay the Phase-A crossover function.

*Deliverable:* the chapter's keystone figure — σ_* one continuous curve from the
canard regime (λ^{1/2}) to the MMO regime (μ^{3/2}), one prefactor C_q.

### Phase E — Integrate as the chapter's unification section

**E1.** Write `MMO_CROSSOVER.md` (derivation + numerics), then fold it into the
Chapter-2 draft as the §5 unification pillar (§9 of the chapter skeleton).

**E2.** Update README: the MMO §7 structural-finding clause and §9.3
methodological insight from "C_q inherited from canard" to "C_q proven
continuous across the folded saddle-node," and add the crossover figure to the
figure index. Banner line on `MMO_FHR_PLAN.md`.

---

## 4. Success criteria / falsifiability

- **Clean success:** Phase A gives the normal forms coinciding at μ → 0 and a
  crossover scale μ_c; Phase B shows η_*(μ) crossing from μ^{3/2} to the canard
  plateau with C_q ≈ 2.8 constant; Phase D shows the full-FHR σ_* continuous
  through the FSN with C_q ≈ 8–10. The inheritance is then a theorem-backed
  continuous limit, not a coincidence.
- **Partial success (still chapter-worthy):** Phases A + B close (normal-form
  crossover proven and numerically confirmed) but the full-FHR μ ↔ λ relation
  (C2) stays empirical. The unification is then established at the universal
  (normal-form) level with the model-specific FHR mapping numerical — the same
  universal-vs-model-specific split as the rest of the chapter.
- **Informative failure:** if η_*(μ) does **not** approach the canard plateau as
  μ → 0 — e.g. C_q drifts with μ — then C_q is *not* one constant and the
  "inheritance" was a near-coincidence. That is a real (negative) result and
  must be reported as such; it would mean the funnel escape carries physics the
  single-passage canard does not.

---

## 5. Hazards

**H1. Normal-form convention.** The K2 normal form `dW/dT = μV − (1+μ)/2` is
quoted from the project's own notes; verify the exact coefficients against
Wechselberger (2005) §4 / Brøns–Krupa–Wechselberger (2006) before building A2–A3
on it. The structural point (μ → 0 ⇒ constant drift ⇒ canard NF) is convention-
robust, but μ_c's prefactor is not.

**H2. The μ → 0 divergence.** s_max = (1−μ)/(2μ) → ∞ as μ → 0: the *local*
theory predicts infinitely many secondary canards exactly at the FSN, while the
physical trajectory does finitely many (funnel-filling). The crossover must be
taken in the FSN *scaling* regime, not at the singular point; expect the
small-μ numerics to be delicate (the same many-canard hazard as Phase 3).

**H3. Which folded saddle-node.** FSN type I (folded saddle ↔ folded node,
transcritical) vs type II (folded node ↔ regular fold / singular Hopf) have
different local scalings. Determine which transition FHR undergoes (Phase C1)
before importing Krupa–Wechselberger scalings; the canard chapter's 2D fold
must be matched to the correct degenerate case.

**H4. A common observable.** σ_pq (staircase dissolution) and σ_* (canard
escape) are *different* measurements as stated; the crossover only makes sense
if both are reduced to the one "noise drives the trajectory off the canard by
its local amplitude V_*" threshold (D1). Define this first; do not compare a
dissolution σ_pq on one side to an escape σ_* on the other.

**H5. Geometric factor changes form.** Because μ^{3/2} and λ^{1/2} are different
powers, do **not** expect a single straight line on a log-log plot across the
whole transition — the crossover function bends. A failed single-power fit is
expected, not a problem; the bend *is* the result.

---

## 6. References

- Wechselberger (2005). Existence and bifurcation of canards in R³ in the case
  of a folded node. *SIAM J. Appl. Dyn. Syst.* 4(1), 101–139.
- Krupa & Wechselberger (2010). Local analysis near a folded saddle-node
  singularity. *J. Differential Equations* 248(12), 2841–2888.
- Brøns, Krupa, Wechselberger (2006). Mixed-mode oscillations due to the
  generalized canard phenomenon. *Fields Inst. Commun.* 49, 39–63.
- Krupa & Szmolyan (2001). Extending GSPT to nonhyperbolic points. *SIAM J.
  Math. Anal.* 33(2), 286–314.
- Berglund & Gentz (2006); Berglund, Gentz, Kuehn (2012, 2015) — the noise side.
- Internal: `CANARD_BLOWUP.md` §4 (the λ^{1/2} law + Θ_crit ≈ 2.8),
  `MMO_K2.md` / `MMO_NOISE.md` (the μ^{3/2} law + C_q ≈ 8–10), `MMO_CHAPTER.md`
  §5, §7.

---

## 7. How it slots into Chapter 2

This crossover is the chapter's **§5 (unification)** in the skeleton from the
chapter-scoping discussion: §1 intro + Jordan obstruction, §2 background, §3
deterministic staircase, §4 noise dissolution, **§5 the crossover — σ_pq is the
canard σ_* continued through the folded saddle-node**, §6 discussion. It is the
movement that elevates the chapter from "one good noise result" to "the noise
law is universal across the fold→folded-node family," and it is the natural
place to cite and connect to Chapter 1 explicitly.

**Effort estimate.** Phases A + B (the normal-form crossover and its numerics)
are ~1–2 focused sessions and deliver the defensible core. Phases C + D (full
FHR) are the heavier lift (~2–3 weeks) and carry the H2–H4 hazards. A or A+B
alone already establishes the keystone claim; C + D upgrade it from
universal-normal-form to demonstrated-in-the-model.
