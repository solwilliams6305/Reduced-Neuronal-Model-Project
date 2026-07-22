# Kickoff — coupled-FHN cusp program: the two top-tier theorems

You're picking up a mature research program on noise-induced escape at a **cusp** in a coupled FitzHugh–Nagumo
system. The abstract law is 𝒲 (the cusp analogue of Tracy–Widom). The program is ~92% done; the spine is complete
(fold/TW, cusp/𝒲, the uniform tube T1, the ladder unification a(q)=3q/(q+2), the closed-form deterministic
connection, the physical signature). **Your job is the two remaining top-tier programs** — the ones that would take
this from "strong applied-math paper" to "top-tier."

## Read first (in order)
1. `coupled-atlas/PROGRESS.md` — the tracker (read the top ~8 log entries; they are the recent work).
2. `regime-tests/FRONTIER_SCOPING_NOTES.md` — **the detailed scoping of your two targets. Read this most carefully.**
3. `regime-tests/COUPLED_CUSP_RESULTS.md` — consolidated results (𝒲, tails, connection, ladder, T1).
4. `regime-tests/CONNECTION_CLOSED_FORM_NOTES.md` — the closed-form connection (your Program-2 foundation).
5. `regime-tests/PERSISTENCE_ITEM2_NOTES.md`, `FRONTIER_PROGRAMS_NOTES.md` — the derived tails + first-pass scoping.

## The goal (prioritized — do 1, then 2)

**Program 1, Tier A** (most tractable, highest ROI): prove 𝒲_β is rigorously the first-explosion law of the
stochastic Weber Riccati `dp=(sign(Y)Y²−p²)d(−Y)+(2/√β)dW`, p~+|Y| recessive — well-defined, moment-determined,
monotone β-family, characterized by the backward-Kolmogorov PDE. Four of five steps are standard or in hand
(well-posedness up to explosion; a.s.-finite explosion by comparison ṗ≤M−p²; Feynman–Kac ⇒ the FP PDE;
β-monotonicity by coupling). **The crux is A4: rigorous two-sided tail bounds** ⇒ moment-determinacy (Carleman).
Right/barrier tail = the T1 confined-Bernstein bound (in hand); left/persistence tail = the Freidlin–Wentzell
instanton (−logP→s⁵/(10η²), derived) made into a rigorous bound. Closing A4 also feeds Tier B. This is a
self-contained, publishable theorem — build it end to end and tag every step proved/cited.

**Program 2, Route 2b** (the closed form, as a resurgent trans-series — the ONLY viable route):
- **Step 0 first:** reconcile the small-η FP vs MC discrepancy in the weak-noise coefficients (skew/η ≈0.7–0.8 FP
  vs ≈1.2 MC; variance-coefficient 0.20 FP vs 0.134 MC/Green's-function). Both methods are numerically fragile at
  small η (FP: thin diffusive boundary layer; MC: sparse deep tail). Pin the leading coefficients with a converged
  method — boundary-layer-resolved FP (adaptive/finer p-grid near the absorbing boundary), or tail-importance-
  sampled MC. They already AGREE at β=2 (skew 0.607 vs 0.61); the disagreement is purely small-η numerics.
- Then: perturbative sectors to ~5 orders (MC-free FP is cleanest); the Borel-plane singularity structure (located
  at the instanton actions you already have); the leading Stokes constant. The FP data already shows the resurgent
  signature (variance-series coefficient ratios GROW, 1.23→2.33 per order ⇒ divergent/asymptotic series). A partial
  resurgent representation ("the defining resurgent structure of 𝒲") is publishable.

**Do NOT pursue (provably obstructed — these are theorems, not gaps):** a Painlevé-σ-ODE or a Fredholm/
determinantal closed form for 𝒲 (isomonodromy fixed point; unbounded-below operator has no bottom eigenvalue). The
τ-function / RH-with-flow route hits the same wall. A clean "no determinantal form" negative would itself be a
theorem worth writing — but don't chase a positive there.

## How this program works (honor these — they are load-bearing, and they set the bar)
- **Tag everything:** proved / cited / numerically-validated / conjectural. Mark load-bearing steps (⚑).
- **Never dress a numerical check as a proof. Never declare a theorem closed unless it is.**
- **Attempt the hard reasoning yourself.** These are proof/derivation tasks — do the analysis, don't punt them.
  Work a load-bearing step (a tail bound, a resurgence derivation) with real effort before concluding anything.
- **But don't grind a wall.** If a route genuinely resists after honest effort — the estimate won't bound, the
  Stokes constant won't pin — STOP, diagnose the obstruction precisely, and map the minimal missing ingredient /
  cleaner sub-problem. **A precise obstruction + a map is a fully acceptable outcome, often better than a forced
  result.** This is the single most important cultural rule; the best results here were clean diagnoses.
- **Cheap decisive checks** against the FP-𝒲 (trusted ground truth) and the MC sampler, every time. Prefer
  **non-confounded** signatures (exponents, pole structure, β/θ-independence, tail mechanism, resurgent coefficient
  growth) over bulk cumulant-fitting — the latter is **confounded** (a smooth Gaussian "passes").
- **Retract overclaims openly.** This program retracts and gains trust by it (recent: the "β^0.6 persistence rate,"
  the √ε-onset "opposite trend," each corrected honestly). Do the same.
- Save notes + figures alongside `regime-tests/` and `coupled-atlas/`; update `PROGRESS.md` with **honest**
  percentages (don't inflate). A prior-entries-below chain keeps the history.

## Practical / sandbox
- **scipy IS available** — parabolic-cylinder (`pbdv`, `pbwa`), `loggamma`, `solve_bvp`, `solve_ivp`. Use it.
- **Trusted tools:** `fp_cusp.py` = MC-free ground truth for 𝒲 (β=2: skew +0.607, exk −0.237); the additive
  Riccati `escape()` pattern = trusted sampler; `connection_closed_form.py` / `resonance_gamma.py` = the closed-form
  connection (both Weber Γ-factors, machine-verified); `instanton_action.py` = the left-tail instanton (const 1/10).
- Vectorize heavy sims over realizations (loop only over time steps); split long runs; background them.
- matplotlib mathtext has **no `\mathcal`** (use plain "W"); avoid `tight_layout` with a mathtext suptitle
  (use `fig.subplots_adjust`); reset/renormalize to avoid overflow in deep field integrations.

## Honest ceiling (state it; don't pretend past it)
A literal closed form for 𝒲 at β=2 is **provably not a 1-D ODE** — a theorem about 𝒲, not a missing calculation.
The realistic top-tier wins are: (1) the rigorous characterization (Program 1 Tier A) — reachable; (2) a partial
resurgent representation (Program 2 Route 2b) — reachable. The full universality-with-rate (Tier B) and full
stochastic exact-WKB are genuine multi-year programs, each now reduced to ONE named hard problem (the cusp
tail/convergence estimate; stochastic exact-WKB). Aim for (1) and (2); frame (3) as the mapped frontier.

The physics is real, the standards are high, and the honest negatives are as valuable as the positives.
