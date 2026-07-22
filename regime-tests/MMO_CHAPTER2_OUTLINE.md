# PhD Chapter 2 — Noise and mixed-mode oscillations in 3D FitzHugh–Rinzel

**Outline + slotting map.** Six sections, ~35–45 pp. Each block below lists
what it contains, which existing material feeds it, the **status** of that
material (new-writing / exposition / measured / derived / open / **new
research**), and the figures it needs. The single new research piece is §5 (the
crossover); everything else is either an existing result to write up properly or
exposition.

**Three-act spine:** set up the folded-node staircase (§3) → derive the noise
law that dissolves it (§4) → prove that law *is* the canard escape law continued
through the folded saddle-node (§5). Acts 1–2 exist; Act 3 is the pillar to add.

---

## §1 Introduction (~3–4 pp) — *all new writing*

- **1.1 Mixed-mode oscillations.** L^a S^b patterns, rotation number
  ρ = L/(L+S), the Farey/devil's staircase; where MMOs occur (neurons, chemical
  oscillators, lasers). Set the intuition before any equations.
- **1.2 The question.** How does degenerate noise (σ on the fast variable v)
  dissolve the staircase — which plateau at what σ_pq, in what order? Why it
  matters: real neurons are noisy, so which rhythms are robust signals vs
  artefacts of a noise-free idealisation.
- **1.3 Why three dimensions.** The Jordan obstruction (2D autonomous flows
  cannot carry an L^a S^b itinerary); FitzHugh–Rinzel as the minimal carrier
  that keeps the v–w subsystem identical to Chapter 1's FHN.
- **1.4 Contributions.** (i) the folded-node staircase reproduced from cycle
  geometry (κ ≈ 2π², funnel-filling f(c)); (ii) the noise law
  σ_pq ≈ σ_*(μ(q)) = C_q√ε·μ^{3/2}; (iii) the crossover proving C_q is the
  **same constant** as Chapter 1's canard escape — noise universality across
  dimension.
- **1.5 Relation to Chapter 1.** State that the σ_* accumulated-variance
  machinery is inherited and recapped in §2.4.

*Figure:* a teaser — sample FHR MMO time series + the deterministic staircase.

---

## §2 Mathematical background (~6–8 pp) — *exposition; cite heavily*

- **2.1 GSPT essentials.** Slow–fast systems, critical manifold, Fenichel,
  the fold. (Krupa–Szmolyan 2001, Kuehn 2015.)
- **2.2 Folded singularities.** Desingularised reduced flow; folded saddle /
  node / focus; the eigenvalue ratio μ = λ_weak/λ_strong. (Wechselberger 2005;
  Desroches et al. 2012.)
- **2.3 Wechselberger's theorem.** Secondary canards, s_max = (1−μ)/(2μ), the
  K1/K2/K3 blow-up charts. State the K2 normal form (used in §3 and §5).
- **2.4 Recap: the canard escape law (Chapter 1).** The accumulated-Brownian-
  variance argument, the K-S blow-up, η = σ/√ε, σ_* = C_q√ε·λ^{1/2}, Θ_crit ≈ 2.8.
  *This subsection is load-bearing for §4 and §5* — write it as a clean,
  self-contained recap, not just a cross-reference.
- **2.5 Noisy slow–fast essentials.** Freidlin–Wentzell action, sample-paths
  (Berglund–Gentz 2006); the degenerate-noise convention.

*Figures:* folded-node geometry schematic (the funnel + secondary canards);
the blow-up coordinate sketch.

---

## §3 The deterministic FHR staircase (~8–10 pp) — *measured/derived; write derivations + figures*

- **3.1 The model.** Working FHR `dv=(v−v³/3−w+y+I)dt, dw=ε(v+a−bw)dt,
  dy=εδ(c−v)dt`; why the canonical Rinzel δ→0 form is a pure spiker (global
  return bypasses the funnel). *[MMO_PHASE1_5 §1 — measured]*
- **3.2 The folded node.** Jacobian `J=[[1+δ,−b],[2bδ(c+1),0]]`,
  μ(c) ∝ (c+1); the 2-slow justification (δ = O(1) ratio, not δ→0).
  *[MMO_PHASE1_5 §2, MMO_TIMESCALE_CHECK — derived]*
- **3.3 The staircase and its exponent.** ρ(c); plateau widths; the universal
  ceiling α = 2 (from μ∝(c+1) + s_max) vs the realised α ≈ 1.55; the
  funnel-filling f(c) as the model-specific gap. *[MMO_ALPHA_DERIVATION —
  derived ceiling + measured realised]*
- **3.4 The K2 return map.** The rotation map `ln R = κμ`, κ ≈ 2π² constant in
  c; f(c) = (2/κ)·ln(a_max/a_min)/(1−μ) with μ cancelling; α reproduced
  (1.45–1.49). *[MMO_K2 — measured + derived]*
- **3.5 What's first-principles vs measured.** Ceiling derived; κ and a_min(c)
  measured (Path-A frontier). *[flag open]*

*Figures:* (a) critical manifold + nullclines + folded node; (b) sample MMOs at
several c; (c) the ρ(c) staircase; (d) the κ ≈ 2π² collapse (three quantities,
one line); (e) f(c) measured vs (★). *Mostly exist in `figures/` — reuse.*

---

## §4 Noise dissolution of the staircase (~8–10 pp) — *the core; results exist*

- **4.1 Setup + the counter.** Degenerate noise on v; Euler–Maruyama; the
  load-bearing fix — count SAO loops in **w** (smooth; v-noise gives ~10³
  spurious maxima). *[MMO_NOISE §2 — measured]*
- **4.2 The signature.** σ_pq decreases with q (high-q-first), dt-converged —
  the predicted Farey-dissolution order. *[MMO_NOISE §3]*
- **4.3 The noise scale.** The measured σ_pq ≈ C_q√ε·μ(q)^{3/2}. The μ^{3/2} is
  **not** a local escape factor (the early Path-B sector-crossing reading is
  superseded by §5): the *local* escape is the canard escape (universal, shared
  C_q), and the μ^{3/2} is a *global* funnel-filling factor (model-specific,
  derivation from f(c) open). *[MMO_NOISE §4 + MMO_CROSSOVER — empirical/global]*
- **4.4 The exponent, stated correctly.** Idealised μ∝1/q gives β=1.5; exact
  μ(c_pq) gives β≈1.13 = measured 1.11; the operative law is σ_pq ≈ σ_*(μ(q))
  (β is the q-slope of σ_* itself, **not** an α/γ on top); γ reported as
  3B-derived/unmeasured, not back-inferred. *[MMO_NOISE §4.1, §5.2 — reconciled]*
- **4.5 The inheritance, observed.** C_q ≈ 8–10 = Chapter 1's full-FHN value;
  state it here as an observation, prove it in §5. *[MMO_NOISE — measured
  consistency]*

*To strengthen before submission:* the converged counting study (error-barred
β, **direct** γ from std(ρ)~σ^γ) — removes the "two-point slope" objection.
*[open — precision]*

*Figures:* (a) ⟨ρ⟩-drift dissolution (high-q-first); (b) σ_pq(q) with the
exact-μ reconciliation + idealised line. *Exist — reuse.*

---

## §5 Unification: the local escape *is* the canard escape (~5–7 pp) — ***NEW RESEARCH (A+B done)***

- **5.1 The claim.** The folded-node noise escape and the canard escape are the
  same local mechanism with the same C_q, connected continuously as μ → 0
  (folded saddle-node). *[MMO_CROSSOVER §4]*
- **5.2 The normal-form spine.** Folded-node K2 form `dW/dT = μV − (1+μ)/2`
  → canard K-S form `dW/dT = −1/2` as μ → 0; same noise η = σ/√ε; so C_q is one
  constant. Crossover scale μ_c ≈ 0.6. *[MMO_CROSSOVER §2, A1–A2 — derived]*
- **5.3 Normal-form numerics (the result, not the originally-guessed one).**
  η_*(μ) → the canard plateau as μ → 0 with **C_q flat ≈ 2.0** (settings-
  dependent absolute value); η_* is canard-like and carries **no μ^{3/2}**.
  *[MMO_CROSSOVER §3, B — validated]*
- **5.4 The refinement = the cleaner story.** Because the local escape carries no
  μ^{3/2}, that factor is **global funnel-filling**, not local. So the noise side
  mirrors the deterministic side exactly: universal local escape (canard, C_q) +
  model-specific global geometry (funnel-filling f(c)). *[MMO_CROSSOVER §4]*
- **5.5 Full-FHR crossover.** σ_*(control) continuous through the FSN, C_q ≈ 8–10
  — the follow-on. *[open — Phases C–D, MMO_CROSSOVER_PLAN.md]*

*Status: A+B done (`MMO_CROSSOVER.md`); C–D open.*
*Figure (keystone, exists): η_*(μ) → canard plateau, C_q flat, μ^{3/2} reference
sitting far below the data (`figures/mmo_crossover_nf.png`).*

---

## §6 Discussion (~3–4 pp) — *new writing*

- **6.1 Established / open.** The honest table from MMO_CHAPTER.md §6.
- **6.2 Universal vs model-specific.** κ ≈ 2π² and C_q universal; f(c), a_min(c)
  model-specific — the same split as canard/tonic; situate in the thesis arc.
- **6.3 Outlook.** First-principles κ and Path-A a_min(c); cross-model (Koper);
  noise-induced MMOs (the pure-spiker FHR kicked into the funnel); the
  entry-exit reformulation (Kaklamanos–Kuehn–Popovic–Sensi) as the
  supervisor-aligned next step.
- **6.4 Limits.** Degenerate-noise assumption; δ = O(1) (not three-timescale);
  small-μ many-canard numerical limits.

---

## Front/back matter

- **Notation table** (v, w, y, ε, δ, μ, λ, σ, η, ρ, q, s, κ, C_q, σ_*, σ_pq).
- **"What's new vs Wechselberger" box** in §1 or §3 — make explicit that the
  deterministic theory is *applied* Wechselberger and the noise law + crossover
  are the original contributions. Examiners will look for this.
- **Reproducibility appendix** — the `mmo_*` scripts per section.

---

## Workload summary (so the page count isn't intimidating)

| Section | New research | New writing | Figures | Status |
|---|---|---|---|---|
| §1 intro | — | full | 1 (reuse) | write |
| §2 background | — | full (exposition) | 2 (new schematics) | write |
| §3 deterministic | — | derivations | 5 (mostly reuse) | results exist |
| §4 noise | precision only | core write-up | 2 (reuse) | results exist |
| §5 crossover | **yes (A–D)** | full | 1 (the keystone, new) | **to do** |
| §6 discussion | — | full | — | write |

**Bottom line:** one new research movement (§5), one cheap precision tidy (§4's
β/γ), and the rest is writing/derivation/figures around results you already
have. The chapter is closer to done than its length implies.
