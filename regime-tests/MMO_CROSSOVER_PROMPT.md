# MMO crossover — Phase A+B (normal-form core): handoff to Opus

This is a research handoff. You are establishing the **keystone unification
result** of the MMO chapter: that the 3D folded-node noise law
`σ_*(μ) = C_q√ε·μ^{3/2}` and the 2D canard escape law `σ_* = C_q√ε·λ^{1/2}` are
two regimes of **one** crossover function with the **same** prefactor C_q,
connected as the folded node degenerates to a folded saddle-node (μ → 0).

You are doing **Phases A and B only** — the normal-form analysis and its SDE
numerics. This is the cheap, high-value core (~1–2 sessions) that proves the
claim at the universal (normal-form) level. The full-FHR confirmation (Phases
C–D) is a separate, heavier follow-on; do not attempt it here.

## 1. Brief context

The project is a regime-by-regime theory of how degenerate noise (σ on the fast
variable v only) affects slow–fast neuronal dynamics. Two chapters are central
here:

- **Canard (Ch.1):** a trajectory follows a repelling slow manifold past a 2D
  fold; noise knocks it off. σ_* = C_q√ε·λ^{1/2}, with λ the slow drift through
  the fold. Derived by accumulated Brownian variance in the Krupa–Szmolyan
  blow-up; C_q ≈ 2.8 in the normal form, 8–10 in full FHN. (`CANARD_BLOWUP.md`.)
- **MMO (Ch.2):** in 3D FitzHugh–Rinzel a folded node creates a funnel of
  secondary canards; degenerate noise dissolves the resulting Farey staircase.
  σ_*(μ) = C_q√ε·μ^{3/2}, μ = λ_weak/λ_strong the folded-node eigenvalue ratio;
  the fitted prefactor is C_q ≈ 8–10 — **the same as the full-FHN canard value**.
  (`MMO_NOISE.md`, `MMO_CHAPTER.md`.)

That shared C_q is currently a *measured consistency*. This task makes it a
*continuous limit*: the canard law is the μ → 0 face of the folded-node law.

## 2. Read first (in order)

1. `regime-tests/MMO_CROSSOVER_PLAN.md` — the full plan; you are executing its
   §3 Phases A and B, with §2 (the normal-form spine) as the analytical target
   and §5 (hazards) as the guardrails.
2. `regime-tests/CANARD_BLOWUP.md` §4 — the canard λ^{1/2} law: η = σ/√ε,
   W_* = λ^{2/3}, V_* = λ^{1/3}, T_window = λ^{−1/3}, η_* = C_q·λ^{1/2}, and the
   normal-form criterion Θ_crit ≈ 2.8 for V_cross = 1. **You will reuse this
   exact escape criterion** so the two sides are the same measurement.
3. `regime-tests/MMO_NOISE.md` §4 — the μ^{3/2} law and the C_q ≈ 8–10 fit.
4. `regime-tests/canard_normal_form_map.py` — the canard normal-form SDE
   simulator you will adapt (the folded-node normal form differs only in the
   dW/dT drift term).

## 3. What's established (rely on, don't re-derive)

- **Canard K-S normal form:** `dV/dT = V² − W,  dW/dT = −λ`. Constant drift,
  single fold passage, no rotation. η_* = C_q·λ^{1/2}, C_q ≈ 2.8 (V_cross = 1).
- **Folded-node K2 normal form:** `dV/dT = V² − W,  dW/dT = μV − (1+μ)/2`
  (Wechselberger's K2 chart — **verify the coefficients** against Wechselberger
  2005 §4 before relying on the constant; Hazard §7.1). The μV term generates
  the funnel rotation and the s_max = (1−μ)/(2μ) secondary canards.
- **The bridge:** as μ → 0, `dW/dT → −1/2` — the rotation term vanishes and the
  drift becomes constant, recovering the canard normal form with λ = 1/2. The
  noise enters identically as η = σ/√ε in both. So C_q is the *same object in
  the same coordinates*; only the geometric factor (λ^{1/2} ↔ μ^{3/2}) changes.

## 4. Task A — the normal-form crossover (analytical)

**A1.** Confirm the folded-node K2 normal form and its μ → 0 reduction to the
canard K-S form. Use **only** the normal form and its limit — do **not**
reconstruct the parabolic-cylinder inner solution (that is Path A of the
deterministic side and is not needed here). State the constant-drift limit
cleanly (λ_eff = 1/2 in normal-form units).

**A2.** Find the crossover scale μ_c by balancing the rotation term μV against
the constant drift over the canard window: μ·V_* ~ O(drift), with V_* ~ λ^{1/3}.
Below μ_c the dynamics is a single canard passage (η_* set by the λ^{1/2},
constant-drift regime); above μ_c the funnel forms (η_* ~ μ^{3/2}).

**A3.** Write the unified geometric factor F(μ, λ) interpolating
F → λ^{1/2} (μ ≪ μ_c) and F → μ^{3/2} (μ ≫ μ_c), with C_q the μ-independent
prefactor. If a closed F is hard, **A1 + A2 alone** (normal forms coincide at
μ → 0, plus the crossover scale) already prove continuity and shared C_q — that
is the defensible core; do not force A3.

## 5. Task B — normal-form SDE numerics (validate the crossover)

**B1.** Integrate the folded-node normal-form SDE
`dV = (V² − W) dT + η dB,   dW = (μV − (1+μ)/2) dT`,
sweeping μ ∈ {0.5, 0.2, 0.1, 0.05, 0.02, 0.01} (down toward the canard limit).
For each μ, measure the escape threshold η_* using the **same operational
criterion as the canard chapter** (noise drives V off the canard by V_*;
Θ = η/η_* with Θ_crit ≈ 2.8 for V_cross = 1). Reuse the integrator and
threshold logic from `canard_normal_form_map.py` — the *only* change is the
dW/dT drift term.

**B2.** Plot η_*(μ). The prediction: at larger μ it follows the `μ^{3/2}` slope;
as μ → 0 it bends to a **μ-independent plateau** at the canard value
(η_* = C_q·λ^{1/2} with λ = 1/2), with **C_q ≈ 2.8 constant throughout**. This
one curve — two regimes, one constant — is the result.

## 6. Falsifiability

- **Success:** η_*(μ) crosses continuously from μ^{3/2} to the canard plateau
  with C_q ≈ 2.8 constant across the whole μ-sweep. The inheritance is then a
  continuous limit, not a coincidence — the chapter's keystone.
- **Partial:** the bridge (A1) and the qualitative crossover (B) hold but C_q
  drifts mildly with μ at the small-μ (delicate) end — report the band and
  attribute to the many-canard numerics (Hazard §7.2), still supporting
  continuity.
- **Informative failure:** η_*(μ) does **not** approach the canard plateau as
  μ → 0, or C_q drifts systematically. Then C_q is *not* one constant and the
  "inheritance" was a near-coincidence — a real negative result; document it,
  because it would mean the funnel escape carries physics the single-passage
  canard does not.

## 7. Hazards

1. **Normal-form convention.** The `−(1+μ)/2` constant is quoted from project
   notes; verify against Wechselberger (2005) §4 / Brøns–Krupa–Wechselberger
   (2006). The structural μ → 0 ⇒ constant-drift fact is convention-robust; the
   μ_c prefactor is not.
2. **Small-μ delicacy.** s_max = (1−μ)/(2μ) → ∞ as μ → 0 — many accumulating
   secondary canards make the small-μ runs delicate (longer T, finer dt). Take
   the crossover in the scaling regime, not at the singular point.
3. **One observable, both sides.** Use the single "noise drives V off the canard
   by V_*" threshold throughout, so the small-μ measurement is literally the
   canard-chapter measurement. Do not compare different criteria across regimes.
4. **Don't expect a single power law.** μ^{3/2} and λ^{1/2} are different powers,
   so η_*(μ) **bends** on log-log — a failed single-power fit is expected; the
   bend *is* the crossover.

## 8. Deliverables

1. `MMO_CROSSOVER.md` (Phases A+B portion) — the normal-form reduction (A1),
   crossover scale (A2), optional F(μ,λ) (A3), and the η_*(μ) numerics (B) with
   the C_q-constant verdict. Structure parallel to `CANARD_BLOWUP.md`. Be
   explicit about the normal-form convention check (§7.1).
2. `mmo_crossover_nf.py` — self-contained; adapts `canard_normal_form_map.py`'s
   integrator + threshold; sweeps μ; writes `data/` + `figures/`.
3. `figures/mmo_crossover_nf.png` — η_*(μ) on log-log with the μ^{3/2} line, the
   canard plateau, and the fitted C_q annotated. The keystone figure.
4. Banner line on `MMO_FHR_PLAN.md`: "(9) Crossover A+B: [result] — η_*(μ)
   crosses μ^{3/2} → canard plateau, C_q = … constant; inheritance [proven /
   what's missing]."

## 9. If A+B close

The keystone is established at the normal-form level. The follow-on (a separate
handoff) is Phases C–D from `MMO_CROSSOVER_PLAN.md`: the full-FHR folded-
singularity map, the μ ↔ λ relation in the FSN scaling regime
(Krupa–Wechselberger 2010), and the continuous full-FHR σ_*(control) through the
folded saddle-node with C_q ≈ 8–10. That upgrades the result from
universal-normal-form to demonstrated-in-the-model and supplies §5.4 of the
chapter. A+B alone already supply §5.1–5.3 and the keystone figure.

## 10. Reading priority if time-constrained

1. `MMO_CROSSOVER_PLAN.md` §2 (the normal-form spine — the whole idea).
2. `CANARD_BLOWUP.md` §4 (the λ^{1/2} law + the η_* criterion you reuse).
3. `canard_normal_form_map.py` (the simulator you adapt).

Start Task B as soon as A1 is confirmed — the numerics are the persuasive
deliverable and reuse existing code; A2/A3 can be finished alongside.
