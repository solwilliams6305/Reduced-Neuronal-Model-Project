# Notes — Δ(g) and the coupled cusp blow-up (deterministic step toward the crossover theorem)

_June 2026. Companion to `NOISY_CUSP_CROSSOVER_NOTES.md`, `PROBLEM_STATEMENT_NOISY_CUSP.md`.
Figure: `coupled-atlas/figures/delta_g_blowup.png`; script `coupled-atlas/delta_g_blowup.py`.
Tags: [PROVED] exact / classical · [DERIVED] asymptotic balance · [NUMERIC] simulated here · [HEURISTIC] argued · [OPEN] needed for a theorem._

Model (Kristiansen–Pedersen): v_i' = −v_i³ + 3v_i − w_i + g(v_j − v_i), w_i' = ε(v_i − c), i=1,2.

---

## 1. Δ(g) — derived from the slow manifold (replaces the geometry-argued map)

**[PROVED — exact algebra].** Symmetric/antisymmetric coordinates v_{1,2} = v_s ± δ,
w_{1,2} = w_s ± δw. Subtracting the two fast equations:

  δ̇ = −δ³ + μ(v_s,g)·δ − δw ,   **μ(v_s,g) = 3(1 − v_s²) − 2g.**

So the antisymmetric mode is **exactly the cusp (A₃) normal form**, with unfolding μ. (The
symmetric mode keeps the plain fold w_s = −v_s³ + 3v_s, folds at v_s = ±1; the coupling
cancels in the sum.) The fast Jacobian determinant confirms it:
det J = a₁a₂ − g(a₁+a₂), a_i = 3(1−v_i²); on synchrony the (1,−1)-eigenvalue is a − 2g = μ.

The antisymmetric cubic δw = −δ³ + μδ has folds at δ = ±√(μ/3) when μ > 0, so the
**fold separation** is

  **Δ(g, v_s) = 2√(μ/3) = 2√( (3(1−v_s²) − 2g)/3 ).**

Evaluated at the symmetric fold v_s = 1 (where the peel-off happens), μ = −2g and

  **Δ(g) = 2√(−2g/3)   (repulsive g < 0);   μ ≤ 0 for g ≥ 0 ⇒ synchrony stable, no antisym folds.**

The **cusp** is μ = 0, i.e. **g = 0** at the symmetric fold. For g < 0 the antisymmetric mode
is unstable (a folded singularity with rotation = the SAO funnel); for g ≥ 0 it is the bare fold.
This is the exact, model-derived version of the fold↔cusp picture.

---

## 2. The coupled cusp blow-up and the Airy↔Weber uniform connection in the charts

Treat δw and μ as slowly varying (driven by the symmetric slow flow, ẇ_s = ε(v_s − c)).
Near the cusp point (δ, μ, δw) = 0, blow up with the **A₃ quasi-homogeneous weights**

  δ = r·D,  μ = r²·M,  δw = r³·W,   (one rung above the fold's (1,2,3) weights).

**[DERIVED].** In the rescaling (central) chart the inner equation is the *rescaled cusp*
D' = −D³ + M·D − W, and the **variational equation about the connecting (canard) solution is
the parabolic-cylinder / Weber equation** — this is the folded-node inner equation (Wechselberger).
The two antisymmetric turning points sit at D = ±√(M/3):

- **M ≫ 1** (μ ≫ r²): the turning points are *resolved/isolated* → each is an **Airy** turning
  point → (with noise) Tracy–Widom. This is the "two separated folds" / fold-class limit.
- **M = O(1)** (μ ~ r²): the turning points *coalesce* → **parabolic cylinder / Weber** — this is
  exactly **Olver's uniform connection** (isolated turning → Airy; two coalescing → Weber),
  now realized inside the blow-up chart. This is the cusp/Weber-class.

So the deterministic Airy↔Weber crossover is the uniform connection of the variational equation
as M passes through O(1). **[PROVED for the connection itself: Olver, *Asymptotics & Special
Functions*, coalescing-turning-points chapter; DERIVED that it appears at M~O(1) in this chart.]**

### Which scale sets the crossover

The *operative* transition for the measured SAO/Weber signature is the **folded-saddle →
folded-node onset**: the folded singularity's eigenvalues scale as √ε (standard folded-node
result), so the eigenvalue ratio ~ μ/√ε, and SAOs become resolved when

  **μ ~ √ε   ⇔   Δ(g) ~ ε^{1/4}   ⇔   g_crit ~ −k√ε.**

(Equivalently r ~ ε^{1/4}, M = μ/√ε.) The √ε here is the folded-node scale, not the pure
A₃ blow-up radius (r~ε^{1/5}, μ~ε^{2/5}); the data below selects √ε.

---

## 3. Tie-back to the numerics (predicted onset & scaling vs measured)

**[NUMERIC]** (`delta_g_blowup.py`). Measuring the peel-off spread crossover at two ε:

| ε | predicted onset (−¾·… see note) | measured g_crit | g_crit/√ε |
|---|---|---|---|
| 0.015 | −(3/8)√ε = −0.046 … −½√ε = −0.061 | **−0.071** | −0.58 |
| 0.030 | −0.065 … −0.087 | **−0.099** | −0.57 |

- **g_crit/√ε is constant (−0.58, −0.57)** and the **ratio g_crit(0.030)/g_crit(0.015) = 1.40 ≈ √2**
  — confirming **g_crit ∝ √ε** (ε^{2/5} would give 1.32; √ε gives 1.41 — the data picks √ε).
- The measured constant k ≈ 0.58 is O(1), bracketed by the two heuristic estimates (3/8, 1/2);
  pinning it exactly is the folded-node onset constant (needs the chart computation).
- g_crit(ε=0.015) = −0.071 sits inside the independently-measured crossover band **−0.05 … −0.11**
  (`crossover_fold_to_cusp.py`). ✓

So the derived Δ(g) + the folded-node scale predict the onset *and* its ε-dependence, both verified.

---

## 4. Status ledger

| statement | status |
|---|---|
| Antisym mode = cusp normal form; μ = 3(1−v²)−2g; Δ(g)=2√(−2g/3) | **[PROVED]** exact |
| Cusp at g=0; fold-class g≥0, folded-node (Weber funnel) g<0 | **[PROVED]** |
| Variational eq. in the rescaling chart = Weber; turning points coalesce as M→0 (Airy↔Weber) | **[DERIVED]** + [PROVED] (Olver connection) |
| Crossover at μ~√ε ⇔ Δ~ε^{1/4} ⇒ g_crit∝√ε | **[DERIVED]** scaling; **[NUMERIC]** confirmed (ratio 1.40≈√2) |
| O(1) onset constant k≈0.58 | **[NUMERIC]**; exact value [OPEN] |

---

## 5. What the stochastic step still needs (Berglund–Gentz tubes → TW→Weber-edge)

The deterministic blow-up above is the scaffold. To turn the crossover into a *stochastic*
theorem (TW for g≳g_crit, Weber-edge for g≲g_crit) one still needs, on top:

1. **[OPEN] Noise covariance in the cusp chart.** Carry additive σ dW through the A₃ rescaling;
   get the effective inner noise η_cusp = σ / ε^{p} (the cusp analogue of the fold's η = σ/√ε₂).
   The exponent p follows from the (δ,δw)=(r,r³) weights with r~ε^{1/4}.
2. **[OPEN] Berglund–Gentz tubes for the cusp.** Sample-path tube estimates around the canard
   one rung above the fold/folded-node they treated — confinement of the noisy trajectory to an
   O(η_cusp)-tube through the rescaling chart, with the exit measure controlled.
3. **[OPEN] The stochastic uniform connection** — the genuine new theorem: the inner exit measure
   is the ground state of the **stochastic Weber (parabolic-cylinder) operator** for M=O(1)
   (g≲g_crit) and of the **stochastic Airy operator** (→ TW) for M≫1 (g≳g_crit), uniformly in M.
   I.e. the noise analogue of Olver's coalescing-turning-points asymptotics. Part 1
   (`weber_vs_pearcey.py`) is the numerical target this would explain.

The deterministic uniform connection (§2) is the classical, tractable first rigorous step; the
stochastic version (1–3) rests on the coupled blow-up and is the project's standing analytic gap.
