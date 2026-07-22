# Phase 3 — the resolution law derived: why h ≲ ε^{2/3} (component B / H1)

*Run 2026-06-02. Executes step 4 of `QUASIPOTENTIAL_FOLDS_PROBLEM_STATEMENT.md`
§9: turn the empirically measured knee (Phase 2, `FOLD_BARRIER_PHASE1_2.md`) into a
*derived* resolution law, tied to the fold blow-up scales and the loss of normal
hyperbolicity. Verification code: `phase3_verify.py`; figure: `fig_fold_box_error.png`.*

---

## Headline

- **Threshold exponent derived: r = 2/3.** A uniform mesh resolves the fold escape
  only while `h ≲ ε^{2/3}`. Measured in Phase 2: **0.70–0.73**. Derived here from the
  Krupa–Szmolyan blow-up box.
- **Above-knee growth derived: q = 2** (relative error `~ (h/ε^{2/3})²` for
  `h ≳ ε^{2/3}`). Measured **≈1.9–2.2** in the well-resolved cases (pooled 1.5 with a
  named bias).
- **Mechanism, confirmed spatially:** the mesh under-resolves the fold *turn* (the
  `ε^{2/3}` slow scale); that error then propagates downstream along the
  saddle-avoidance valley — so the error map lights up the whole manifold-hugging
  route, not an isolated box.

This is the slow-fast-*specific* breakdown the risk register demanded — it comes
from the non-hyperbolic fold's anisotropic blow-up, not from generic "caustics are
hard."

---

## 1. Derivation of the threshold r = 2/3 (rigorous, from the blow-up)

Near the lower fold `(v_f,w_f)=(−1,−2/3)` write `v=−1+ξ`, `w=−2/3+η`. The drift
reduces to the fold normal form (re-derived, consistent with Phase 1):

```
b_v ≈ ξ² − η ,     b_w ≈ ε g₀ ,    g₀ = −1 − a   (constant in the fold region).
```

Balancing the fast drift against the slow drift fixes the **Krupa–Szmolyan inner
scales**: with `ξ = ε^{1/3}X`, `η = ε^{2/3}Y`, `T = ε^{1/3}t`,

```
dX/dT = X² − Y ,   dY/dT = g₀ ,
```

i.e. the fold region is an **anisotropic box**

```
Δv ~ ε^{1/3}   (fast direction)  ×   Δw ~ ε^{2/3}   (slow direction).
```

The minimum-action escape route (Phase 2: the manifold-hugging / saddle-avoidance
path) rounds the fold *inside this box*. A **uniform** mesh of spacing `h` (same in
`v` and `w`) resolves the box geometry only if it resolves **both** extents:

```
h ≲ ε^{1/3}   AND   h ≲ ε^{2/3}.
```

Since `ε^{2/3} < ε^{1/3}` for `ε<1`, the **slow direction binds**:

```
        h_crit ~ ε^{2/3}        ⇒        r = 2/3 .
```

This is exactly H1's statement ("`h ≲ ε^{2/3}`, the slow-direction inner scale").
**Why this is fold-specific:** the anisotropy `ε^{1/3}:ε^{2/3} = 1:ε^{1/3}` *diverges*
as `ε→0`, and exists only because the fold is non-hyperbolic (the fast Jacobian
eigenvalue → 0 there). At a normally hyperbolic point the manifold carries O(1)
scales — no shrinking box, no resolution wall. The breakdown is caused by the loss
of normal hyperbolicity, as required.

## 2. Derivation of the growth exponent q = 2 (variational)

The barrier is a **minimum** of the Freidlin–Wentzell action, `B = min_path S[path]`.
For `h > ε^{2/3}` the mesh cannot place the rounding path within `ε^{2/3}` of the
optimal turn; the best grid-representable path is displaced by `δ ~ h`. Because the
first variation vanishes at a minimum, the leading action error is **quadratic**:

```
B_h − B ≈ ½ ⟨δπ, S''[path*] δπ⟩ ,   S'' ~ B / (ε^{2/3})²,   δπ ~ h
⇒   (B_h − B)/B  ~  (h / ε^{2/3})²        ⇒        q = 2 .
```

So below the knee the error sits at the solver floor; above it, it climbs
quadratically in `h/ε^{2/3}`.

## 3. Numerical confirmation (`phase3_verify.py`)

| Quantity | Derived | Measured | Notes |
|---|---|---|---|
| Threshold `h_crit ∝ ε^{r}` | r = 2/3 ≈ 0.667 | **0.73** (0.69 stricter filter) | over ~1.5 decades of ε |
| Above-knee growth `~(h/ε^{2/3})^q` | q = 2 | **1.88, 2.16** (ε=10⁻², 3·10⁻³) | clean, well-resolved cases |
| ″ (small ε) | q = 2 | 1.14, 1.62 (ε=10⁻³, 3·10⁻⁴) | biased low — see caveat |

**Spatial error map (`fig_fold_box_error.png`).** Error vs a fine reference is
**not** an isolated box at the geometric fold; it lies along the manifold-hugging
(saddle-avoidance) valley emanating from the fold turn (only ~20% of the error sits
in the immediate fold strip). This is the expected behaviour of a marching solver:
the mis-resolved fold turn injects the `(h/ε^{2/3})²` error, which then **propagates
forward along the action-minimising characteristics** down the valley.

## 4. Honest caveats

- **Small-ε bias.** For ε ≤ 10⁻³ the finest affordable `h` does not fully reach the
  plateau, so the reference `B*` is slightly high → measured errors (and the growth
  exponent) are biased *low*. This is exactly why the pooled `q` is 1.5 while the
  resolved cases give ~2; it also slightly *steepens* the threshold fit (pushing the
  measured `r` to 0.73 above the true 2/3). Finer meshes / Richardson extrapolation
  would tighten both toward the derived values.
- **q = 2 is a scaling (variational) argument,** not a rigorous estimate of `S''`;
  the data support it but do not pin it to 3 figures.
- **The clean `ε^{1/3}×ε^{2/3}` box aspect ratio is not directly recoverable** with
  this Dijkstra-class solver, because the fold-turn error propagates along the
  valley. A higher-order / anisotropic OLIM and an error-localisation diagnostic
  would be needed to image the box itself.
- **Non-degenerate (isotropic) noise only.** Degenerate noise (real FHN) is Phase 5.

---

## 5. What Phase 3 settles, and what's next

**Settles (the §9 step-4 goal):** the empirical knee is now a *derived* law —
`h_crit ~ ε^{2/3}` from the blow-up box (binding on the slow scale), with quadratic
error growth above it from the variational structure. Combined with Phase 1 (closed
form) and Phase 2 (the OLIM-class breakdown + saddle-avoidance), the **bankable core
(steps 1–4) is complete**: *closed-form fold truth → solver breaks at the fold in
the singular limit → here is the derived resolution law, tied to the loss of normal
hyperbolicity.*

**Sharpened thesis (one sentence):** *uniform-mesh quasipotential solvers lose the
slow-manifold-hugging (saddle-avoidance) escape route at a fold once `h ≳ ε^{2/3}`,
the Krupa–Szmolyan slow inner scale, with relative barrier error `~(h/ε^{2/3})²` —
a breakdown caused specifically by the fold's loss of normal hyperbolicity.*

**Next:**
1. **Tighten the constants** — finer meshes / Richardson to pull `r→2/3`, `q→2`
   cleanly; or swap in Cameron's actual OLIM (brief §6) to confirm the law is
   solver-class-independent.
2. **The σ–ε regime boundary** — derive where the along-manifold (∝ε) route gives
   way to the frozen fast escape `(4/3)δ^{3/2}`; this is the bridge between this
   work and the project's existing fold-Kramers results.
3. **Degenerate noise (Phase 5)** — repeat with `D̂=diag(1,0)`; the slow direction
   is then noise-free, which should *sharpen* the box and the law.
4. **Read Börner et al. 2024 in full + send the Popović email** (still the two
   open Step-0 gates) before writing for publication.
