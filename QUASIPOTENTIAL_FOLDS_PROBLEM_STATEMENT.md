# Problem statement — Quasipotential landscapes at slow-fast folded singularities

*A kickoff brief for a fresh investigation. Self-contained: it assumes no prior
conversation. Read §0 first — your first action is a literature check that can
kill or pivot the whole thing.*

---

## 0. How to use this document (and the one thing to do first)

You are being asked to investigate whether a specific gap is genuinely open and,
if so, to start attacking it. The gap, in one line:

> **No one has computed the Freidlin–Wentzell quasipotential *landscape* at a
> slow-fast folded singularity, nor tested whether the standard quasipotential
> solvers survive the singular limit ε→0 there.**

**STEP 0 — confirm the gap before deriving anything.** This brief was written
from a targeted-but-shallow literature pass. Before investing, verify the
intersection is empty:

- Search recent (2023–2026) output of **Maria Cameron** (OUM/OLIM quasipotential
  solvers), **Tobias Grafke / Eric Vanden-Eijnden** (instantons, gMAM for
  multiscale / rare events), and **Freddy Bouchet** (large deviations for
  slow-fast / averaging).
- Search **Berglund–Gentz–Kuehn**, **Wechselberger**, **Desroches** for any
  "quasipotential near a fold / folded node" or "large-deviation landscape for
  canards."
- Search terms: `quasipotential singularly perturbed`, `large deviations slow-fast
  folded node`, `instanton canard fold`, `mean exit time folded singularity`,
  `ordered upwind method stiff multiscale`.
- **Human step (highest value):** ask the supervisor (Popović) directly whether
  Kuehn's group has touched quasipotential-near-fold. Popović co-authors with
  Kuehn (Kaklamanos–Kuehn–Popović–Sensi, *JDDE* 2025), so this is a one-email
  check that could save a month.

If a paper already bridges quasipotential computation and folded singularities,
**stop and report it** — the contribution then shifts to whatever sub-question
they left open. If the two literatures are still sitting apart (the working
assumption), proceed to §1.

**Project files with the closed-form ground truth you will need as benchmarks**
(read these if the workspace is mounted):
`NOISE_ERROR_BOUNDS.md`, `regime-tests/FOLD_ESCAPE_PREFACTORS.md`,
`regime-tests/MMO_C0_PINNED.md`, `README.md` (§3, Excitable).

---

## 1. Thesis

Two mature literatures do not talk to each other:

- **Quasipotential numerics / rare-event ML.** Cameron's Ordered Upwind (OUM)
  and Ordered Line Integral (OLIM) methods; neural-net and SINDy methods that
  learn the orthogonal (gradient + rotational) decomposition of the drift
  (arXiv:2012.09111, 2306.11418, 2409.06886). These compute the global
  quasipotential landscape and mean-exit-times around **fixed-point attractors**,
  in 2–3D, on **smooth O(1) systems**, benchmarked on **Maier–Stein**.

- **Geometric singular perturbation theory (GSPT) of folds and canards.**
  Krupa–Szmolyan blow-up, folded nodes, canards; on the stochastic side,
  Berglund–Gentz–Kuehn **sample-path** concentration bounds near canards/folds.
  This community works with the **fold/folded-node normal forms** and their
  ε-asymptotics, but treats noise via path-concentration, **not** via a global
  quasipotential landscape.

The quasipotential solvers have **never been pointed at a folded singularity**,
where (i) the drift loses normal hyperbolicity (its Jacobian has an eigenvalue
→ 0), and (ii) in the physically relevant FitzHugh–Nagumo case the noise is
**degenerate** (rank-deficient diffusion). Both break the assumptions the solvers
rest on. The opportunity is to bring the quasipotential-landscape object into the
folded-singularity setting, where closed-form ground truth already exists to
benchmark against — and to characterize precisely where and why the standard
tools fail as ε→0.

---

## 2. Setting and notation (self-contained)

**Running model — stochastic FitzHugh–Nagumo, degenerate noise.**

```
dv = ( v − v³/3 − w + I ) dt + σ dW_t      (fast, noisy)
dw = ε ( v + a − b w ) dt                   (slow, deterministic)
```

with `(a,b) = (0.7, 0.8)`, timescale separation `0 < ε ≪ 1`, noise `0 < σ ≪ 1`.
Noise enters **only** the fast variable ⇒ diffusion matrix `D = diag(σ², 0)` is
**rank-deficient (degenerate / hypoelliptic).**

**The quasipotential.** For fixed ε, in the small-noise limit σ→0,
Freidlin–Wentzell theory assigns a quasipotential `V_ε(x)` relative to the stable
rest state `x*`, solving the stationary Hamilton–Jacobi equation

```
H(x, ∇V) = 0 ,     H(x,p) = ½ ⟨p, D̂ p⟩ + ⟨b(x), p⟩ ,
```

where `b` is the drift `(v−v³/3−w+I, ε(v+a−bw))` and `D̂ = diag(1,0)` (noise
direction). `V_ε(x)` measures the σ²-scaled cost of being driven from `x*` to `x`;
escape rate `~ exp(−2 ΔU_ε / σ²)` with barrier `ΔU_ε = V_ε(escape point) − V_ε(x*)`.

**The standard decomposition the solvers use.** Away from degeneracies one writes

```
b = −D̂ ∇V + ℓ ,        with     ⟨∇V, ℓ⟩ = 0 ,
```

(`ℓ` the transverse/rotational drift). OUM/OLIM solve the HJ equation on a mesh;
the ML/SINDy methods learn this decomposition from trajectories. **Two structural
problems at a fold:** `D̂` is rank-deficient (so the decomposition's standard
elliptic form is invalid), and at the fold `Db` has a zero eigenvalue (loss of
normal hyperbolicity), so `ℓ` and `∇V` become ill-separated.

**Fold normal form (Krupa–Szmolyan blow-up).** Near the lower fold of the
v-nullcline, the canonical desingularized system is, to leading order,

```
dV/dT = V² − W ,     dW/dT = −λ ,
```

obtained by `v = v_f + ε^{1/3} V`, `w = w_f + ε^{2/3} W`, `t = ε^{−1/3} T`
(constants `v_f, w_f, λ` model-specific; **re-derive and verify against the
project files** — do not take them on faith). The key takeaway is the
**anisotropic inner scaling**: the fold region is a box of size
`ε^{1/3}` (fast direction) × `ε^{2/3}` (slow direction).

---

## 3. The precise open problem (three components)

**(A) Closed-form quasipotential at the folded geometries, as ε→0.**
Construct `V_ε` near the fold via matched asymptotics / blow-up, in closed form
where possible, for:
  1. the planar **fold** (FHN excitable escape),
  2. the **folded node** (3D FitzHugh–Rinzel; eigenvalue ratio `μ = λ_w/λ_s`),
  3. (adjacent, different mechanism) the **sub-Hopf resonator** escape, as a
     contrast case.
Output: the barrier `ΔU_ε`, its leading ε-scaling, and the O(1) constant. The
fold case has known closed form (see §8) — start there.

**(B) The structural failure prediction (the heart — see §5).**
Show analytically why the orthogonal-decomposition solvers must lose accuracy in
an ε-dependent neighborhood of the fold, and predict the error scaling.

**(C) Benchmark + numerical demonstration (the bankable core — see §6).**
Run the public solvers on the singularly-perturbed systems; measure error against
the closed-form `ΔU_ε` from (A); exhibit the predicted breakdown.

Components (A) and (B) are derivation-heavy and are the original content;
(C) is publishable on its own (see §9).

---

## 4. The two axes of difficulty (why this is genuinely hard, not just untested)

1. **Loss of normal hyperbolicity.** At the fold the critical manifold is
   tangent to the fast fibers; `Db` has a zero eigenvalue. The gradient/rotational
   split `b = −D̂∇V + ℓ` is ill-conditioned there.
2. **Degenerate noise.** `D̂ = diag(1,0)` is rank-1. The HJ Hamiltonian is
   degenerate in `p_w`; the minimum-action path must "borrow" fluctuations through
   the drift coupling. Most solvers assume non-degenerate (or at least full-rank
   anisotropic, cf. OLIM arXiv:1806.05321) diffusion. **Rank-deficient D at a fold
   appears doubly untested** — verify this in Step 0; if true it is a second,
   independent open axis.

A clean intermediate target that isolates axis (1): use **non-degenerate** noise
(`D̂ = I`) first, nail the fold-geometry breakdown, then switch on degeneracy.

---

## 5. Central falsifiable hypothesis

> **H1 (resolution breakdown).** A uniform-mesh quasipotential solver (OUM/OLIM)
> with spacing `h` computes the fold barrier `ΔU_ε` with controlled error only
> while `h ≲ ε^{2/3}` (the slow-direction inner scale). For `h ≳ ε^{2/3}` the
> relative error in `ΔU_ε` grows like `(h / ε^{2/3})^q` for some `q > 0`.
> Determining `q` (analytically and numerically) is part of the contribution.

> **H2 (decomposition breakdown).** The learned-decomposition methods
> (NN / SINDy-on-instanton) degrade near the fold because (a) `Db` is singular
> there, making `∇V`/`ℓ` non-unique to leading order, and (b) a fixed polynomial
> SINDy basis cannot represent the cubic-fold structure in a shrinking
> `O(ε^{1/3})` neighborhood without growing degree. Predicted signature: error
> concentrated in the fold boundary layer, worsening as ε→0 at fixed model order.

Both are concrete, testable, and **falsifiable** — if the solvers sail through
the fold with ε-independent error, the thesis is wrong and you report that
(still a useful negative result).

---

## 6. Benchmark protocol

**Systems (in order):**
1. Planar fold normal form `(V²−W, −λ)` with non-degenerate noise — cleanest test
   of H1.
2. 2D FHN above, non-degenerate noise — realistic geometry, known `ΔU_ε`.
3. 2D FHN, **degenerate** noise (`D̂=diag(1,0)`) — the real problem.
4. (stretch) 3D FitzHugh–Rinzel folded node.

**Solvers:**
- Cameron's OLIM (software is public from her UMD page) — primary grid solver.
- One learned method: the orthogonal-decomposition NN (arXiv:2012.09111) and/or
  SINDy-on-instanton (arXiv:2409.06886). Reimplement minimally if needed.
- Cross-check `ΔU_ε` independently with a direct **geometric minimum-action**
  (gMAM) computation of the escape instanton.

**Measurements:**
- Relative error in `ΔU_ε` vs mesh `h`, swept across `ε ∈ {10⁻¹,…,10⁻³}`. Fit `q`
  in H1; locate the `h ~ ε^{2/3}` knee.
- Spatial error map: confirm error concentrates in the `ε^{1/3}×ε^{2/3}` fold box.
- For learned methods: error vs SINDy degree / NN capacity at fixed ε.

**Ground truth:** the closed-form `ΔU_ε` from §3(A)/§8, plus direct Monte-Carlo
mean-exit-time `τ ~ exp(2ΔU_ε/σ²)` as an end-to-end sanity check.

---

## 7. What is NOT open — cite, do not claim

- **The integrability principle** ("a closed-form quasipotential exists ⟺ the
  Wentzell–Freidlin Hamiltonian is completely integrable"). This is **folklore**
  (classical Hamilton–Jacobi / Liouville; Maier–Stein is the canonical
  non-integrable example). Use it as framing; claim nothing.
- **Generic "benchmark ML quasipotential against a closed-form case."** Done —
  Maier–Stein is *the* established analytic-ground-truth benchmark. Novelty must
  be the **slow-fast / folded** benchmark class and the **breakdown mechanism**,
  not benchmarking per se.
- **Sample-path bounds near canards/folds** belong to **Berglund–Gentz–Kuehn**.
  This project's object (the global quasipotential landscape + solver behavior) is
  different, but position explicitly against BGK and do not re-claim their bounds.

---

## 8. Closed-form ground truth available (benchmark values)

From the existing project (verify by re-derivation; constants sourced from
`FOLD_ESCAPE_PREFACTORS.md`, `NOISE_ERROR_BOUNDS.md`, `MMO_C0_PINNED.md`):

- **Fold barrier:** `ΔU(δ) ~ B · δ^{3/2}`; fold geometric coefficient
  `K_fold = 4/3`; attempt-frequency prefactor `A₀ = 1/π`.
- **Escape rate:** `λ(δ) ~ A₀ δ^{1/2} exp(−K_fold δ^{3/2}/σ²)`.
- **Critical noise:** `σ_crit = √(ε/C)`, `C = K_fold/A₀`; unified form
  `σ_crit = √(4π D ln2 · ε)` with `D = λ` (canard) or `D = g ≈ 0.31` (excitable).
- **Canard quantile law:** `C_q(p) = √(−4π ln(1−p))`, median `C_q = √(4π ln2) ≈ 2.95`.
- **BGK canard-spacing constant (folded node):** `c₀ = π²/16 ≈ 0.6169` (leading
  order), depth function `c₀(φ₀) = π² sin²φ₀ / [(1+μ)(2φ₀+sin2φ₀)²]`.

These are exactly the analytic targets the solver community lacks for this class —
they make the benchmark in §6 quantitative.

---

## 9. Milestones and the bankable core

1. **Step 0 lit check** (§0). Gate everything on this. (~days)
2. **Fold quasipotential in closed form**, non-degenerate noise; reproduce
   `ΔU(δ) ~ (4/3)δ^{3/2}` from the HJ/instanton side. (component A, fold)
3. **OLIM benchmark on the planar fold**; test H1, fit `q`, find the `h~ε^{2/3}`
   knee. **← bankable core: a clean, publishable result by itself.**
4. **Structural explanation** of the breakdown via blow-up inner scales
   (component B / H1). Turn the empirical `q` into a derived one.
5. **Degenerate-noise fold** (axis 2) and **learned-method breakdown** (H2).
6. **Folded node** (stretch); position against BGK.

The floor is steps 1–3: *"the standard quasipotential solvers break near folds in
the singular limit; here is the closed-form truth and the resolution law."* Even
if 4–6 stall, that is a complete, useful paper bridging two communities.

---

## 10. Risk register

- **Gap already filled.** Mitigation: Step 0. If Cameron/Grafke/BGK did it, pivot
  to their open sub-question.
- **Breakdown is "just caustics are hard."** The novelty must be the *slow-fast-
  specific* mechanism (ε-scaled boundary layer from blow-up; loss of normal
  hyperbolicity) and the *derived* resolution law `q`, not a generic remark.
  Mitigation: tie `q` explicitly to the `ε^{1/3}/ε^{2/3}` inner scales.
- **Folded-node closed form may not exist.** Fine — that becomes a "numerical-
  only" entry, which *is* the point of the closed-form-vs-numerical contrast.
- **Degenerate noise too hard.** Mitigation: the non-degenerate fold (steps 2–4)
  is already a complete result; degeneracy is a second paper.

---

## 11. Key references (fetch these)

- Cameron, *Finding the quasipotential for nongradient SDEs* (Physica D, 2012);
  OUM/OLIM software — UMD page `math.umd.edu/~mariakc/Quasipotential.html`.
- *Ordered Line Integral Methods for the quasi-potential* — arXiv:1706.07509;
  anisotropic-diffusion variant arXiv:1806.05321; 3D nongradient arXiv:1808.00562.
- *A Data Driven Method for Computing Quasipotentials* — arXiv:2012.09111.
- *Computing large deviation prefactors via machine learning* — arXiv:2306.11418.
- *Quasi-potential and drift decomposition by sparse identification (SINDy)* —
  arXiv:2409.06886.
- Krupa–Szmolyan, fold blow-up; Wechselberger / Desroches et al., folded nodes &
  canards (review: *Mixed-Mode Oscillations…*, SIAM Review 2012).
- Berglund–Gentz, *Noise-Induced Phenomena in Slow-Fast Dynamical Systems* (2006);
  Berglund–Gentz–Kuehn, noisy canards/MMO (arXiv:1312.6353 and 2012 work).
- Freidlin–Wentzell, *Random Perturbations of Dynamical Systems* (quasipotential,
  HJ, integrability) — the folklore source for §7.

---

*End of brief. First action: §0 Step 0. Do not derive before confirming the gap.*
