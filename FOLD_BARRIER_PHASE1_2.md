# Phase 1 & 2 — fold quasipotential: closed form + the resolution law (H1)

*Run 2026-06-02. Executes steps 1–3 of `QUASIPOTENTIAL_FOLDS_PROBLEM_STATEMENT.md`
§9 (the bankable core), non-degenerate (isotropic) noise. Code and figures in the
folder: `phase1_fold_barrier.py`, `phase2_solver.py`, `phase2_analyze.py`,
`fig_*.png`. Every constant re-derived, not inherited; the solver validated against
an exact non-gradient quasipotential before being trusted on the fold.*

---

## Headline

1. **Phase 1 (closed form), verified two independent ways and against the stored
   project constants.** The planar-fold escape barrier is `ΔU(δ) = (4/3) δ^{3/2}`
   (`K_fold = 4/3`); the Freidlin–Wentzell instanton action — i.e. the
   quasipotential exponent — is `ΔV = 2ΔU = (8/3) δ^{3/2}`; the Kramers prefactor is
   `A₀ = 1/π`. Symbolic check: all PASS.

2. **Phase 2 (numerics).** A uniform-grid OLIM-class solver, validated to a ~0.5%
   floor on an exact case, then run on a slow–fast fold across ε ∈ [3·10⁻⁴, 3·10⁻²]:

   - **H1 is confirmed.** The mesh must satisfy `h ≲ ε^{2/3}` to compute the fold
     barrier; the measured resolution exponent is **q ≈ 0.69–0.73** (fit 0.728) vs
     the **predicted 2/3 = 0.667**, holding over ~1.5 decades of ε
     (`fig_hcrit_vs_eps.png`).
   - **A finding that reshapes the framing:** the 2D quasipotential barrier is *not*
     the frozen `(4/3)δ^{3/2}`. The minimum-action path hugs the critical manifold
     through the fold (paying only the `O(ε)` slow-drift cost), so the barrier
     scales **`B* ∝ ε^{0.95} ≈ ε`**, ~200× below the frozen value. This is the
     **saddle-avoidance** mechanism — an independent reproduction of Börner–Grafke–
     Feudel (PRR 2024), the paper Step 0 flagged as the must-read.

---

## Phase 1 — closed-form barrier (`phase1_fold_barrier.py`, all checks PASS)

Convention (matches Cameron / the project files): `dx = b dt + σ dW`, isotropic
noise, rate `~ exp(−V/σ²)`, Hamiltonian `H(x,p)=½|p|²+b·p`, `V` solves
`H(x,∇V)=0`. For a gradient fast field `b=−U'`, `V = 2U`.

**Route A — frozen fast-subsystem potential (HJ/gradient side).** At frozen slow
value δ>0, the fast drift is `b(x)=x²−δ=−U'(x)`, `U(x)=δx−x³/3`. Critical points
`x=∓√δ` (well at `−√δ`, `U''=+2√δ`; saddle at `+√δ`, `U''=−2√δ`). Barrier
`ΔU = U(√δ)−U(−√δ) = (4/3)δ^{3/2}` ⇒ `K_fold = 4/3`. Kramers prefactor
`√(U''_min|U''_max|)/(2π)/√δ = 1/π` ⇒ `A₀ = 1/π`.

**Route B — Freidlin–Wentzell instanton (minimum-action side).** On the zero-energy
manifold `H=0`, the non-trivial branch is `p = −2b = 2(δ−x²)`. The escape action
`S = ∫_{−√δ}^{+√δ} p\,dx = (8/3)δ^{3/2} = 2ΔU`. The two routes agree, and the
escape rate is `λ(δ) = (√δ/π) exp(−(8/3)δ^{3/2}/σ²)`.

**Cross-check (the discipline the brief demands).** Identical to the stored
derivation in `FOLD_ESCAPE_PREFACTORS.md` (Part 2): `K_fold=4/3`, `A₀=1/π`, exponent
`8/3`. The **instanton route (B) is new** here — the project file derived the
barrier only via potential+Kramers; route B is the independent second derivation
step 2 of §9 asked for. The exact `(4/3)` is robust; the prefactor `1/π` inherits
the same finite-barrier `O(1)` caveat the file already flags (Kramers is large-
barrier asymptotic).

---

## Phase 2 — the solver, its validation, and H1

**Solver (`phase2_solver.py`).** Uniform-grid Freidlin–Wentzell *geometric action*
(Heymann–Vanden-Eijnden): the path cost is `S=∫(|b||x'|−b·x')ds`, discretised on a
square mesh as a shortest-path problem with nonnegative edge weights
`w(X→Y)=|b_mid||Y−X|−b_mid·(Y−X)`, solved by heap Dijkstra on a 16-neighbour
stencil. This is a minimal but genuine member of the OUM/OLIM ordered-upwind family
(it is *not* Cameron's exact OLIM — see limitations).

**Validation (must pass before trusting the fold).** On the linear non-gradient SDE
`dx=−Ax dt+σ dW`, `A=[[1,β],[0,1]]`, the exact quasipotential is
`V(x)=½xᵀΣ₁⁻¹x` with `AΣ₁+Σ₁Aᵀ=I`. The solver reproduces it to **~0.5% median**
error for both β=0 (gradient) and β=1 (non-gradient). The error is an **h-independent
floor** set by the fixed stencil's angular discretisation — the known signature of
ordered-upwind solvers, and the baseline against which fold errors are judged.

**Test system (isotropic noise, explicit ε, a stable rest state):**
`dv=(v−v³/3−w)dt+σdW₁`, `dw=ε(v−a)dt+σdW₂`, with `a=−1.3` (excitable rest on the
left branch; lower fold at `v=−1`, `w_f=−2/3`; `δ=w*−w_f≈0.099`).

**Finding 1 — saddle avoidance (`fig_Bstar_vs_eps.png`).** The converged barrier to
the threshold scales `B*(ε) ≈ 0.030·ε^{0.95}` — i.e. linear in ε, ~200× below the
frozen `(4/3)δ^{3/2}=0.083`. The minimum-action path does not jump straight across
the fast direction; it crawls along the critical manifold (`b_v=0`) and rounds the
fold, paying only the `O(ε)` slow drift. **The frozen fold barrier is a
fast-subsystem object, not the 2D quasipotential.** This independently reproduces
the Börner–Grafke–Feudel saddle-avoidance / "flat quasipotential" effect, and it
ties Step 0's two adjacent literatures together (Cameron-type solvers meet the
saddle-avoidance phenomenon at the fold).

**Finding 2 — H1 confirmed (`fig_barrier_vs_h.png`, `fig_hcrit_vs_eps.png`).** For
each ε, the computed barrier `B(h,ε)` plateaus once `h≲ε^{2/3}` and rises above it.
Defining `h_crit(ε)` as the mesh at which the relative error reaches 25%:

| ε | ε^{2/3} | B*(=finest h) | B*/ε | h_crit(25%) | h_crit/ε^{2/3} |
|---|---|---|---|---|---|
| 3·10⁻² | 0.0965 | 1.21·10⁻³ | 0.040 | 0.038* | 0.40* |
| 1·10⁻² | 0.0464 | 3.68·10⁻⁴ | 0.037 | 0.036 | 0.77 |
| 3·10⁻³ | 0.0208 | 1.08·10⁻⁴ | 0.036 | 0.015 | 0.70 |
| 1·10⁻³ | 0.0100 | 3.95·10⁻⁵ | 0.040 | 0.0055 | 0.55 |
| 3·10⁻⁴ | 0.0045 | 1.59·10⁻⁵ | 0.053 | 0.0029 | 0.65 |

A log–log fit (excluding the coarse-grid-limited ε=0.03*) gives
**`h_crit ≈ 0.98·ε^{0.73}`**, i.e. exponent **q≈0.73**, against the H1 prediction
`q=2/3=0.667`. With a stricter grid filter the fit gives q≈0.69. Either way the
mesh requirement scales as the **slow-direction inner scale ε^{2/3}**, exactly as
H1 predicts — over ~1.5 decades of ε.

---

## Honest limitations

- **Solver.** Self-implemented OLIM-*class* Dijkstra, not Cameron's OLIM; ~0.5%
  angular floor; modest grids. Swapping in the real OLIM and finer/anisotropic
  meshes is the obvious hardening step (brief §6).
- **q is q≈0.70 ± ~0.06, not 3-figure.** The smallest-ε plateaus are not fully
  reached at the finest affordable h, biasing those B* slightly high and h_crit
  slightly low (hence the 0.55 row). Consistent with 2/3; not yet a clean 2/3.
- **Reframing, not refutation, of the brief.** The barrier that breaks down is the
  along-manifold (saddle-avoidance) quasipotential `∝ε`, **not** the frozen
  `(4/3)δ^{3/2}`. H1's *resolution law* `h~ε^{2/3}` holds; but the object it applies
  to must be stated as the slow-manifold-hugging 2D quasipotential. This is the
  single most important correction to the problem statement.
- **Non-degenerate noise only** (axis 1). Degenerate/hypoelliptic noise — the real
  FHN — is untouched here (brief Phase 5).

---

## What this does to the program

- **Phase 1 banks.** Closed-form fold barrier, two independent derivations, matching
  the stored constants. Step 2 of §9 is done.
- **Phase 2 / H1 banks, with a twist.** The central H1 claim — uniform-mesh solvers
  need `h≲ε^{2/3}` near a fold — is empirically confirmed (q≈0.70). This is the
  bankable result of §9 step 3.
- **The framing tightens around Börner et al. 2024.** Since the breakdown is of the
  *saddle-avoidance* quasipotential, the contribution is best stated as: *the
  mesh-resolution law for the slow-manifold-hugging quasipotential at a fold* —
  uniting Cameron-type solvers with the Börner saddle-avoidance mechanism. Read that
  paper in full before writing (still the top to-do, with the Popović email).
- **Next (Phase 3):** derive q analytically from the blow-up — does the action
  accumulated in the `ε^{1/3}×ε^{2/3}` fold box give exactly 2/3? — and map when the
  along-manifold route vs the frozen fast escape dominates (the σ–ε regime diagram).
