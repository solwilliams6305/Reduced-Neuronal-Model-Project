# MMO / SAO chapter — scope and task plan

**Status:** planning document. The MMO chapter is the natural next regime
after the canard, tonic, resonator, and excitable chapters close. This
document scopes the framework, tasks, and analytical machinery so that the
chapter can be picked up as a focused project (estimated 6–10 weeks of
focused work, depending on supervisor input).

> **⚠ Phase 1 executed — premise corrected (see `MMO_PHASE1.md`).** The
> deterministic `L^a S^b` devil's staircase this plan attributes to **autonomous
> 2D FHN** (§1, §2.1) **cannot exist there**: a planar limit cycle is a simple
> closed (Jordan) curve and cannot realise a>1 large + b>1 small loops in one
> period. Numerics confirm it (`mmo_2d_check.py`: 0/18 I-values show L+S
> coexistence; ρ(I) is a 0→1 step, no plateaus). The staircase **is** real once a
> third slow variable is added: **3D FitzHugh–Rinzel** (`mmo_fhr_staircase.py`)
> shows a clean MMO devil's staircase (ρ locks at 1, 2/3, 1/2, 1/3, …; Δ_pq falls
> with q). **Recommended redirect: take the chapter into FitzHugh–Rinzel**
> (autonomous, folded node, on the roadmap, keeps the canard-blow-up/entry–exit
> machinery). The §5 Phases 2–4 carry over with the model swapped; §2.1's "2D
> limit cycle" framing is the one part that must be rewritten. Hazard §8.1
> ("return-map smoothness fails at the maximal canard") is subsumed by this
> deeper, topological obstruction.

---

## 1. Question

Inside the canard explosion window (I just above I_H1, width O(ε)), the
deterministic limit cycle exhibits **mixed-mode oscillations**: sequences
of L^a S^b patterns, where L is a large relaxation spike and S is a small
canard-amplitude oscillation around the unstable Hopf FP. Across the
window, the ratio (a per b) traces a **devil's staircase** organised by
Farey arithmetic: plateaus at every rational p/q rotation number, with
mediant patterns L^(p₁+p₂) S^(q₁+q₂) appearing between L^p₁ S^q₁ and
L^p₂ S^q₂.

The chapter's questions:

1. **Deterministic structure.** Predict the L^a S^b patterns as a
   function of (I, ε) — equivalently, the rotation number ρ = L/(L+S)
   across the canard explosion window. Pin down plateau widths Δ_pq
   as a function of denominator q.

2. **Noise dissolution.** Add degenerate noise σ. The discrete L^a S^b
   classification becomes a probability distribution P(k | I, ε, σ)
   over SAO counts. The σ at which plateau p/q dissolves into its
   neighbours gives σ_pq. The chapter predicts the **dissolution
   sequence** — high-q plateaus first, then mediants, then low-q —
   and the q-dependent σ_pq scaling.

3. **Cross-chapter unification.** Connect σ_pq to the canard chapter's
   σ_* = C_q · √ε · λ^{1/2}. The natural conjecture is σ_pq ~ σ_* / q^α
   for some α set by mode-lock width geometry, making the MMO chapter
   a refinement of the canard chapter at the level of integer-valued
   observables.

---

## 2. Deterministic framework

### 2.1 Setup

For I ∈ (I_H1, I_H1 + O(ε)), the FHN limit cycle inherits structure from
the canard explosion:
- L = full relaxation excursion (v sweeps from −2 to +2, both attracting
  branches).
- S = small canard oscillation, amplitude ε^{1/4} → O(1) across the
  explosion.
- L^a S^b = a large spikes plus b small loops per period of the symbolic
  return map.

The rotation number is ρ = L / (L + S), where L counts large excursions
and S small loops in one full pattern repetition.

### 2.2 Return map and Farey staircase

The deterministic MMO structure is captured by a **return map** R: state
just after a large spike → state just before the next large spike,
crossing the canard region in between. R is a 1D map on a Poincaré
section transverse to the cycle. The mode-locked windows L^a S^b
correspond to **rotation numbers ρ = a/(a+b) ∈ Q** of R.

Standard circle-map theory (cf. Boyland 1986, Sturmian sequences) gives:

- Each rational ρ = p/q has a mode-locked window of width Δ_pq in I.
- Plateau widths satisfy Δ_pq ~ 1/q^α for some α ∈ [2, 3].
- Mediant rule: between p₁/q₁ and p₂/q₂ (adjacent in Stern–Brocot),
  the Farey mediant (p₁+p₂)/(q₁+q₂) has its own plateau.
- Total plateau measure approaches 1 in the singular limit ε → 0; the
  rising portion of ρ(I) is supported on a Cantor set of measure zero.

### 2.3 The MMO-specific question

For FHN specifically, the return map R is determined by the canard
explosion's blow-up. **The chapter's deterministic deliverable** is:

- Construct R analytically near the canard fold via the Krupa–Szmolyan
  blow-up (inherited from Chapter 2).
- Predict Δ_pq scaling: α = ? (likely α ≈ 2 from naive circle-map
  arguments; the canard structure may modify this).
- Validate numerically across the canard explosion window.

---

## 3. Noise framework

### 3.1 Mode-lock dissolution

With degenerate noise σ, each fold passage inside an MMO pattern is
subject to the canard-chapter peel-off mechanism. Noise blurs the
mode-locked windows by an amount δρ_noise(σ) set by the noise scale's
contribution to the return map. A plateau dissolves when

```
δρ_noise(σ_pq)  ~  Δ_pq.
```

### 3.2 Order of dissolution

Smaller plateaus (high q) dissolve first because Δ_pq ~ 1/q^α decreases
faster than δρ_noise depends on the cycle's complexity. The predicted
order as σ grows:

```
σ small : depth-3+ mediants dissolve  →  fine Farey structure invisible
σ moderate : depth-2 mediants (L^2 S^3, L^3 S^2) dissolve
σ larger : L^1 S^2 and L^2 S^1 dissolve
σ largest : even L^1 S^1 (the widest plateau) dissolves
σ >> dissolution : staircase becomes smooth ρ(I)
```

This sequence is the chapter's headline noise observable.

### 3.3 σ_pq scaling

Combining the canard chapter's σ_* = C_q · √ε · λ^{1/2} with the
plateau width Δ_pq ~ 1/q^α and a return-map-noise propagation factor,
the natural prediction is

```
σ_pq  ~  σ_*  ·  q^{−β}   for some β > 0,
```

with β determined by how the noise propagates through one pass of the
return map. This is a single-parameter family of thresholds indexed by
the plateau denominator.

---

## 4. Analytical engines

Three machineries combine for the MMO chapter:

1. **Krupa–Szmolyan blow-up** (Chapter 2). The fold passage during each
   SAO is the same fold blow-up that gave σ_* in the canard chapter.
   The MMO return map's geometry near the fold is determined by this
   blow-up, with the parameter being the SAO count k.

2. **Return-map / circle-map theory.** Standard mode-locking analysis
   (Boyland 1986, Glass & Mackey 1988): rotation numbers, Farey
   sequences, devil's staircase. Apply to the MMO return map.

3. **Freidlin–Wentzell action on the return map.** Noise-induced
   transitions between adjacent mode-locked windows are large-
   deviation events governed by an action computed on the return map's
   1D dynamics. The action gives σ_pq via Arrhenius scaling.

4. **Entry-exit functions** (Nikola's 2025 paper with Kaklamanos,
   Kuehn, Sensi). The return map's structure at the canard fold is
   exactly the kind of slow-manifold passage that entry-exit
   relations describe. This is the most direct technical bridge to
   the supervisor's current research line.

---

## 5. Task sequence

### Phase 1 — Deterministic foundations (1–2 weeks)

**Task A:** Build an L^a S^b pattern detector. Numerical FHN
trajectories at fixed (ε, I) inside the canard explosion window;
detect large vs small excursions by amplitude threshold (e.g., v > 0
for L, v ∈ [−1.2, −0.8] only for S); output the symbolic sequence.

**Task B:** Sweep I across the canard explosion window at three ε
values; compute rotation number ρ(I, ε) = L count / total. Plot the
devil's staircase. Identify visible plateaus and measure their widths
Δ_pq numerically.

**Task C:** Fit Δ_pq ~ 1/q^α; report α. Compare to standard circle-map
α = 2 and to canard-specific predictions.

### Phase 2 — Analytical return map (2–3 weeks)

**Task D:** Construct the deterministic return map R via Krupa–Szmolyan
blow-up at the canard fold. The map's parameter is k (SAO count); the
map's structure determines which L^a S^b sequences are admissible at
which I. Cross-reference with Nikola's entry-exit framework.

**Task E:** From R, derive Δ_pq analytically; compare to Task C. The
deterministic chapter closes here.

### Phase 3 — Noise-broadened staircase (2–3 weeks)

**Task F:** Add degenerate noise σ to the FHN integrator; re-sweep I
at fixed (ε, σ) for several σ values; measure P(k | I, ε, σ). Plot
the noise-broadened staircase.

**Task G:** Identify σ_pq for the visible plateaus (the σ at which
each plateau's bin probability drops below 0.5). Fit σ_pq vs q;
predict and verify the q^{−β} scaling.

**Task H:** Connect to canard chapter: show σ_pq / σ_*(λ_explosion)
collapses across I to a function of q alone, confirming the
cross-chapter unification.

### Phase 4 — Write-up (1–2 weeks)

**Task I:** Draft `MMO_CHAPTER.md` with derivation + numerical
validation. Same structure as `CANARD_BLOWUP.md` and `TONIC_PHASE.md`:
question, framework, mechanism, numerical evidence, status, open
items.

**Task J:** Integrate MMO row into the master `README.md` §1 table
with derived σ_pq scaling, replacing the current "planned (PhD scope)"
status with "derived + validated."

**Task K:** Update the master synthesis to reflect MMO as the fifth
regime chapter; discuss the integer-valued-observable framing as
distinct from σ_crit threshold framings of the other chapters.

---

## 6. Estimated effort

- Phase 1: 1–2 weeks (mostly numerical; deterministic FHN integration,
  pattern detector, devil's staircase plot).
- Phase 2: 2–3 weeks (analytical work on the return map; the most
  technically demanding phase; intersects directly with Nikola's
  expertise).
- Phase 3: 2–3 weeks (noise sweeps; potentially expensive numerics if
  long simulation times needed at low σ).
- Phase 4: 1–2 weeks (write-up).

**Total: 6–10 weeks.** Phase 2 is the rate-limiting analytical work
and is the right place to have a focused session with Nikola.

---

## 7. Supervisor-alignment notes

The MMO chapter sits at the technical intersection of:

- **Canard theory** (Wechselberger, Brøns, Krupa, Szmolyan, Popovic) —
  the chapter inherits the canard chapter's blow-up.
- **Entry-exit functions** (Kaklamanos, Kuehn, Popovic, Sensi 2025) —
  Phase 2's return map is exactly an entry-exit relation across the
  canard fold.
- **Circle-map / mode-lock theory** (Boyland, Glass & Mackey) — the
  standard machinery for the devil's staircase.
- **Stochastic slow-fast** (Berglund, Gentz, Kuehn 2012) — the noise
  framework on the return map is the Hunting French Ducks framework
  applied to MMO transitions.

**Key supervisor questions for the meeting before starting the
chapter:**

1. Has anyone done noise-broadening of MMO devil's staircases in
   FHN-class systems? If so, what's the right reference to build on
   / cite / distinguish from?
2. Is the entry-exit framework the right analytical engine for the
   return map, or is there a cleaner approach (e.g., direct
   matched-asymptotic in the singular limit)?
3. What's the right journal target for the MMO chapter once it's done?
   *Nonlinearity*, *SIADS*, or a more neuroscience-oriented venue?

---

## 8. Open analytical hazards

Specific places where the chapter could get stuck and require longer
analytical sessions:

1. **The return map's smoothness.** Standard mode-lock theory assumes
   the return map is a smooth diffeomorphism of the circle. For FHN's
   canard return map, smoothness might fail at the maximal canard
   (where the trajectory's behaviour switches between SAO and L). If
   so, the devil's staircase is replaced by a different fractal
   structure.

2. **Non-perturbative regimes.** At larger σ, the noise drives
   trajectories outside the canard window entirely (canard chapter's
   σ_* threshold). Beyond this point, MMO patterns are dominated by
   escape events not by mode-lock dissolution. The chapter needs to
   identify where σ_pq < σ_*(λ_explosion) — i.e. where the MMO
   framework supersedes vs is superseded by the canard framework.

3. **Connection to the tonic chapter.** At the upper edge of the
   canard explosion window, MMO patterns continuously deform into the
   tonic regime's pure L cycle. The chapter needs to capture this
   transition cleanly — probably as a vanishing-S limit of the L^a S^b
   structure.

---

## 9. Key references

- Boyland (1986). Bifurcations of circle maps: Arnol'd tongues, bistability
  and rotation intervals. *Commun. Math. Phys.* 106, 353–381.
- Brøns, Krupa, Wechselberger (2006). Mixed mode oscillations due to the
  generalized canard phenomenon. *Fields Inst. Commun.* 49, 39–63.
- Krupa, Popovic, Kopell (2008). Mixed-mode oscillations in three time-scale
  systems: a prototypical example. *SIAM J. Appl. Dyn. Syst.* 7, 361–420.
- Desroches et al. (2012). Mixed-mode oscillations with multiple time
  scales. *SIAM Review* 54(2), 211–288.
- Kaklamanos, Kuehn, Popovic, Sensi (2025). Entry-exit functions with
  intersecting eigenvalues.
- Berglund, Gentz, Kuehn (2012). Hunting French ducks in a noisy
  environment.
- Glass & Mackey (1988). *From Clocks to Chaos.* Princeton University
  Press.

---

## 10. Status

**Planning document only.** No code written, no numerics run, no
derivations executed. The chapter is queued behind the supervisor
meeting and the writeup of the current canard + tonic + Van der Pol
material. This document exists to capture the framework so that the
chapter can be picked up as a focused project when the time comes.
