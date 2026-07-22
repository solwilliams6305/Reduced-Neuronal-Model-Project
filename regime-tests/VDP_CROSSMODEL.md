# Van der Pol cross-model test: universality vs model-specificity

**Status:** derived + numerically confirmed. Establishes that the two flagship
structural forms from the FHN chapters are **universal for 2D relaxation
oscillators with a generic cubic fold**, with **model-specific O(1) prefactors**:

| Form (universal) | FHN(0.7,0.8) | Van der Pol | verdict |
|---|---|---|---|
| Canard escape `σ_* = C_q·√ε·λ^{1/2}` | C_q ≈ 8–10 (full) / 2.8 (normal form) | **C_q ≈ 8** | exponent universal; prefactor O(1) |
| Tonic `A_mid = √(c/π²)` | c ≈ 1.55 | **c ≈ 0.78** | form universal; c clearly model-specific |

Companion script: `vdp_crossmodel.py`. Figure: `figures/vdp_crossmodel.png`.

---

## 1. Why Van der Pol

Van der Pol is the minimal control on what is FHN-specific vs generic. It shares
FHN's cubic fast nullcline (`f(x) = x − x³/3`, folds at `x = ±1`) but has the
**simplest possible slow dynamics** — a vertical slow nullcline, i.e. `b = 0`:

```
dx = ( x − x³/3 − y ) dt + σ dW,        (fast, degenerate noise as in FHN)
dy = ε ( x − a ) dt.                    (slow)
```

The fixed point is `x* = a`. Because `tr J = 1 − x²`, the Hopf bifurcations sit
**exactly at the folds** `a = ±1` (unlike FHN, whose Hopf is offset O(ε) from the
fold by the `−εb` term). For `|a| < 1` the FP is unstable and the system
relaxation-oscillates (the tonic analogue). This isolates two questions:

- Does the **canard exponent** survive when the Hopf–fold offset is removed?
- Are the **tonic constants** `c`, `A_mid` genuinely model-specific?

## 2. Part A — canard escape exponent (universal)

Blowing up the left fold (`x = −1 + ε^{1/3}X`, `y = y_f + ε^{2/3}Y`, `t = ε^{−1/3}T`,
`y_f = −2/3`) gives the **identical canonical fold normal form** as FHN,

```
dX/dT = X² − Y,    dY/dT = −λ,    λ_vdp = 1 + a,
```

with the same effective noise `η = σ/√ε`. So the autonomous canard law must be
`σ_* = C_q·√ε·λ^{1/2}` with the **same ε^{1/2}, λ^{1/2} exponents**; only `C_q`
can differ. `vdp_crossmodel.py` Part A tests this on the **full** VdP SDE (not the
normal form): start on the attracting branch at `W_init = 5λ^{2/3}`, integrate,
record `R_hit = W_hit/λ^{2/3}` at the first crossing of `x_cross = −1 + ε^{1/3}`,
and sweep the scaled noise `Θ = σ/(√ε λ^{1/2})` across a grid of `(ε, λ)`.

**Result.** The median `R_hit(Θ)` curves for all six `(ε, λ)` cells collapse onto
a single curve (the σ_* law absorbs all the ε- and λ-dependence): `R_hit` rises
monotonically from the fold-edge (negative, canard survives) through the
canonical window to early noise-driven escape, with the pooled median crossing
the canonical window `R = 1` at

```
Θ_crit = C_q^vdp ≈ 8.
```

This is comparable to the **full-FHN** value (`C_q ≈ 8–10`) and well above the
bare normal-form value (`2.8`): in both models the finite-ε corrections to the
parabolic normal form (the cubic's `O(ε^{1/3})` sub-leading term, here at
ε = 0.02–0.04) stiffen the canard and lift the prefactor by the same ~3×. The
**exponent is identical**; the prefactor is an O(1) model+protocol constant. The
universal ε^{1/2}·λ^{1/2} canard law is confirmed on a second model.

## 3. Part B — tonic `A_mid` (model-specific constant)

For the tonic VdP cycle (`a = 0`, symmetric relaxation oscillation) the same
adjoint-Floquet pipeline as `tonic_phase_response.py` (with the VdP Jacobian
`J = [[1−x², −1],[ε, 0]]`) gives the phase-diffusion amplitude
`A_mid² = (1/T²)∮Z_v²dt`:

| ε | T_cycle | A_mid (VdP) | c = π²A² |
|---|---|---|---|
| 0.04 | 58.1 | 0.302 | 0.90 |
| 0.08 | 33.8 | 0.285 | 0.80 |
| 0.16 | 20.6 | 0.255 | 0.64 |

So `A_mid^vdp ≈ 0.28` and `c_vdp ≈ 0.78` — the **same** `A_mid = √(c/π²)` form as
FHN, but a **clearly different constant** (`c_vdp ≈ 0.78` vs `c_FHN ≈ 1.55`, a
factor ≈ 2). The ε-trend mirrors FHN's (mild decrease, the same sub-leading
`O(ε^{1/3})` correction identified in TONIC_CMID_BVP.md §5). The difference in `c`
is exactly what the first-principles formula predicts: `c = π²·∑_folds ε^{−5/3}R̃_i/T²`
with `R̃_i` depending on the model's fold drift `λ_i` and slow-drift `g_i`, which
differ between FHN (`b = 0.8`, sloped nullcline) and VdP (`b = 0`).

## 4. The universality/specificity split — stated cleanly

> Across FHN and Van der Pol, the **scaling laws are universal** — the canard
> escape threshold scales as `√ε·λ^{1/2}` and the tonic CV amplitude as
> `A_mid = √(c/π²)` in both models, because both are governed by the same
> Krupa–Szmolyan fold normal form `dX/dT = X²−Y` (forward equation → canard;
> adjoint equation → tonic `Z̃_v`). The **prefactors are model-specific O(1)
> constants** set by each model's global cycle geometry: `C_q` (≈ 8 for both full
> models at these ε; 2.8 in the bare normal form) and `c` (1.55 for FHN, 0.78 for
> VdP). This is the same split the canard and tonic chapters asserted for FHN
> alone, now demonstrated to be a genuine universality statement rather than an
> FHN coincidence.

Removing FHN's Hopf–fold offset (VdP's `b = 0`) does **not** change the canard
exponent — confirming the offset is a prefactor effect, not an exponent effect,
as CANARD_BLOWUP.md §10 argued.

## 5. Reproduce

```
python3 regime-tests/vdp_crossmodel.py
```

Part A runs the full-VdP canard escape collapse (≈30 s); Part B runs the VdP
adjoint-Floquet tonic `A_mid`. Writes `results/vdp_crossmodel/summary.txt` and
`figures/vdp_crossmodel.png`.

## 6. Caveats

The full-VdP `C_q ≈ 8` is read from the pooled `R_hit` crossing; at ε = 0.02–0.04
there is finite-ε scatter across cells (collapse spread up to ~1 R-unit near the
crossing), so `C_q^vdp` should be read as "O(8), comparable to full FHN," not a
3-figure number. The clean, unambiguous model-specificity result is the tonic
`c_vdp ≈ 0.78` vs `c_FHN ≈ 1.55`. A first-principles inner-adjoint BVP for `c_vdp`
(the VdP analogue of TONIC_CMID_BVP.md) and smaller-ε canard runs to pin
`C_q^vdp` asymptotically are the natural follow-ups.
