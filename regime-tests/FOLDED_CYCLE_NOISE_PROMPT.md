# Prompt — Stochastic blow-up of the folded limit cycle

*Self-contained brief for a fresh Opus chat. This extends the project's existing
stochastic-blow-up engine (built for the fold and the folded node) to the one
singularity in the family nobody has noised. **Do Step 0 before any derivation.***

---

## 0. Orientation, and what to read first

**Goal.** Push degenerate noise through the *folded-limit-cycle* blow-up and derive a
critical-noise law `σ_*`, exactly as the project already did for the fold
(`CANARD_BLOWUP.md`) and the folded node (`MMO_NOISE.md`) — but for the one case BG and
BGK never treated.

**Read these project files first — they are the engine you are extending, not background:**
- `CANARD_BLOWUP.md` — the fold template. §3 "degenerate noise pushed through the
  blow-up" (`η = σ/√ε`), §4 FW accumulated-Brownian-variance ⇒ `σ_* = C_q√ε·λ^{1/2}`,
  §6/§10 numerical regime-map validation. **Copy this structure.**
- `MMO_NOISE.md` + `MMO_FHR_PLAN.md` + `MMO_K2_ROUTE_AB.md` — the folded-node version
  (Wechselberger K2-chart blow-up return map, `κ ≈ 2π²`, staircase dissolution
  `σ_pq ~ σ_*·q^{−α/γ}`). Same engine, harder chart; `MMO_K2_ROUTE_AB.md` is your
  inner-chart return-map precedent (the Path-A analog for §5).
- `MMO_CROSSOVER.md` — proves the folded-node *local* escape → the canard local escape as
  `μ→0` (one universal local-escape mechanism, shared `C_q ≈ 8–10`). Use it: Channel A's
  `σ_*` should **inherit that same `C_q`** — a built-in cross-chapter consistency check.
- `TONIC_PHASE.md` + `TONIC_CMID_AIRY.md` — the **Floquet / iPRC phase-diffusion** engine.
  You will need it (see §3, Channel B) — this is what makes the limit-cycle case different.
- `VDP_CROSSMODEL.md` — Van der Pol (the canonical relaxation / limit-cycle oscillator) as
  the cross-model universality check; the natural second testbed for §4.
- `PRIOR_BLOWUP_WORK_AUDIT.md` — the map of what's already done and why it collided with BG/BGK.

**Step 0 — gate the novelty (do before deriving; the project has been burned 3× by skipping this):**
1. **Read Jelbart–Kuehn–Kuntz 2024**, *Geometric Blow-Up for Folded Limit Cycle Manifolds
   in Three Time-Scale Systems* (arXiv:2208.01361, J. Nonlinear Sci. 34:17). This is the
   **deterministic scaffold** you add noise to. Reproduce their normal form, weights, and
   charts. Do **not** trust the summaries below for constants — get them from the paper.
2. **Pipeline check.** As of the last audit: no published stochastic version exists; the
   natural extender (Kuntz) moved to ML; Ahsan–Dankowicz–Kuehn 2025 do noisy limit cycles
   but by *adjoint/covariance*, not blow-up. Re-check arXiv for anything newer.
3. **The Popović email** (he collaborates with Kuehn): "is the stochastic/blow-up folded
   limit cycle free, or in Christian's pipeline?" Gate the headline novelty claim on this.

## 1. Why this is genuinely new (the hook — keep it central or you re-collide with BG/BGK)

The fold (Berglund–Gentz) and the folded node (BGK) are done, and this project already
re-derived both. **The folded limit cycle is different in one essential way: the folding
object is a *sustained fast oscillation*, so it carries a phase.** Therefore noise must be
tracked in **two channels at once**:

- **Channel A — amplitude escape** off the folded-cycle manifold (blow-up + FW accumulated
  variance — the `CANARD_BLOWUP.md` move).
- **Channel B — phase diffusion** of the fast oscillation through the passage (Floquet /
  iPRC — the `TONIC_PHASE.md` move).

Neither BG (fold: no phase) nor BGK (folded node: the rotation is *slow*, the fast
subsystem still jumps — not a sustained fast limit cycle) handles a folding fast
oscillation. **This is the first place the project's two engines (canard blow-up +
Floquet phase reduction) must combine.** That combination, and the A-vs-B interplay, is
the contribution. The moment your analysis reduces to pure amplitude escape you are back
in BG/BGK — so keep the fast-oscillation phase central.

## 2. Setup

- **Scaffold (from JKK 2024):** three-timescale "semi-oscillatory" systems with two small
  parameters; a regular folded limit-cycle manifold; blow-up charts covering the relative
  timescale of the *angular* (phase) dynamics vs the parameter drift. Extract their exact
  weights and the rescaling-chart normal form.
- **Testbed:** their worked class — Liénard equations with periodic forcing in the slow
  equation. State one clean minimal example to simulate.
- **Noise:** add `σ dW` to the fast variable(s) only (project convention; matches FHN
  noise-in-v and your `|q_v|²` degenerate-noise treatment).

## 3. The derivation — two channels

**Channel A (amplitude escape — the canard template).**
1. Carry `σ dW` through the JKK rescaling chart. **Derive the effective noise scale**
   `η(σ, ε₁, ε₂)` in blown-up coordinates from their weights (the analog of `η = σ/√ε`;
   do NOT assume the fold's exponent — the folded-cycle weights differ).
2. Identify the inner-chart passage-window scales (the analogs of `W_*, V_*, T_window`
   from `CANARD_BLOWUP.md` §4) from the rescaling-chart dynamics.
3. FW accumulated-Brownian-variance: `η²·T_window ~ (escape amplitude)²` ⇒ the Channel-A
   `σ_*` law in terms of `ε₁, ε₂` and the folded-cycle geometry.

**Channel B (phase diffusion — the tonic template).**
4. The fast limit cycle has a phase; degenerate noise gives phase-diffusion rate
   `D_φ = σ²·⟨Z_v(φ)²⟩` (iPRC `Z_v`, exactly the tonic-chapter CV machinery).
5. **Through the fold the cycle amplitude → 0**, so `Z_v` and `D_φ` become singular /
   acquire a blow-up scaling — track `D_φ` in the inner chart (this is the new coupling).
6. Compare timescales: does phase diffusion randomise the passage before amplitude escape
   fires, or after? Which channel sets `σ_*`, and how does the answer move as `ε₂ → 0`?

**Deliverable:** `σ_*(folded cycle)` plus the **regime diagram** of A-vs-B dominance.

## 4. Numerical validation (the project's standard — mirror `CANARD_BLOWUP.md`)

- **Normal-form regime map:** integrate the blown-up SDE, sweep noise × geometry, build
  `R_hit` vs rescaled noise `Θ`, check the collapse (distribution independent of the
  geometric parameter at fixed `Θ`) — exactly `CANARD_BLOWUP.md` §6.
- **Worked-example map:** full Liénard-with-slow-forcing model, à la §10, to get the
  prefactor `C` and confirm the form survives finite-`ε` corrections.

## 5. Routes and honesty

- **Path B first (your proven route):** accumulated variance + numerics. This is the
  bounded first result; it is what landed the fold and folded-node chapters.
- **Path A (hard, flag as open):** the rigorous inner solution + error control through the
  JKK charts (the analog of the Wechselberger K2 inner solution you left open for the
  folded node). Attempt only after Path B lands.
- **Positioning:** state explicitly against BG (fold), BGK (folded node), JKK (deterministic
  folded cycle). The claim is the *stochastic* folded cycle + the amplitude/phase interplay.

## 6. First deliverable (bounded — "even if small")

The effective noise scale `η` in the JKK rescaling chart **plus** the Channel-A `σ_*` law
for the simplest folded-limit-cycle normal form, validated on one normal-form regime map.
That alone is a real, self-contained first result and the natural first lemma — structured
exactly like `CANARD_BLOWUP.md` §3–§6. Channel B and the interplay come next.

## 7. References

- **Jelbart, Kuehn, Kuntz (2024)** — folded limit cycle blow-up, arXiv:2208.01361. *(scaffold)*
- Berglund & Gentz (2006), *Noise-Induced Phenomena in Slow–Fast Dynamical Systems*, Ch. 5. *(fold = done)*
- Berglund, Gentz, Kuehn (2015), arXiv:1312.6353. *(folded node = done)*
- Krupa & Szmolyan (2001), SIAM J. Math. Anal. 33. *(fold blow-up)*
- Wechselberger (2005), SIADS 4. *(folded-node blow-up)*
- Project: `CANARD_BLOWUP.md`, `MMO_NOISE.md`, `MMO_FHR_PLAN.md`, `MMO_K2_ROUTE_AB.md`,
  `MMO_CROSSOVER.md`, `TONIC_PHASE.md`, `TONIC_CMID_AIRY.md`, `VDP_CROSSMODEL.md`,
  `PRIOR_BLOWUP_WORK_AUDIT.md`.

---

*First action: Step 0 (read JKK 2024; pipeline check; Popović email). Do not derive before the scaffold is in hand.*
