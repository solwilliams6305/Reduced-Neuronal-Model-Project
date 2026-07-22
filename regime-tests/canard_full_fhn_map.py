"""
canard_full_fhn_map.py
----------------------
Full-FHN canard-explosion regime map (Roadmap §9 item 1 / CANARD_BLOWUP.md §9).

Mirrors canard_normal_form_map.py but integrates the full FHN SDE in original
(v, w) coordinates, detecting canard peel-off via the left-fold / middle-saddle
crossing.

Sweep:
  eps   ∈ {0.04, 0.08}
  I     = I_H1(eps) + delta,  delta ∈ {0, 0.25, 0.5, 1.0} * eps
          (delta ~ O(eps), not O(sqrt(eps)), to stay inside canard window)
  Theta = sigma / sigma_*,  sigma_* = sqrt(eps) * sqrt(lambda)
          Theta in logspace(-0.8, 1.0, 8)
  N = 300 trajectories per grid point, dt = 5e-3, T_max = 400.

Detector (left-fold side only):
  (a) v_new > v_M(w_new) AND w_new > w_f   — noise-driven saddle crossing
  (b) w_new < w_f                            — deterministic descent past fold
  Right-fold crossing (w > I+2/3) is NOT a trigger.

Outputs:
  data/canard_full_fhn_hit_location.npz
  figures/canard_full_fhn_regime_map.png
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

# --- shim first, then kernel ---
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _ROOT)
sys.path.insert(0, _HERE)

import _shim  # noqa: F401  — installs brentq shim if scipy missing
from kernel import FHN2D

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
A_FHN = 0.7
B_FHN = 0.8
I_FOLD_L = (A_FHN - 1.0 + 2.0 * B_FHN / 3.0) / B_FHN   # ≈ 0.2917


# ---------------------------------------------------------------------------
# Vectorised Cardano: middle root v_M(w) of  v - v³/3 = w - I
# Only meaningful when |w - I| < 2/3 (between the folds).
# ---------------------------------------------------------------------------

def v_middle(w: np.ndarray, I: float) -> np.ndarray:
    """
    Middle (repelling) root of  v - v^3/3 - (w - I) = 0.
    Clips argument to [-1,1] so it is safe to call for any w,
    but the result is only the true middle root for |w - I| < 2/3.
    """
    arg = np.clip(-1.5 * (w - I), -1.0, 1.0)
    theta = np.arccos(arg)
    roots = np.stack([
        2.0 * np.cos(theta / 3.0),
        2.0 * np.cos(theta / 3.0 - 2.0 * np.pi / 3.0),
        2.0 * np.cos(theta / 3.0 - 4.0 * np.pi / 3.0),
    ], axis=-1)
    return np.sort(roots, axis=-1)[..., 1]   # middle root


def v_left_root(w0: float, I: float) -> float:
    """
    Attracting (left) root of  v - v³/3 = w0 - I  by bisection on [-3, -1].
    """
    def f(v):
        return v - v**3 / 3.0 - (w0 - I)
    lo, hi = -3.0, -1.0
    flo = f(lo)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if abs(fm) < 1e-14 or (hi - lo) < 1e-13:
            return mid
        if flo * fm < 0:
            hi = mid
        else:
            lo, flo = mid, fm
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------------------
# Deterministic sanity check
# ---------------------------------------------------------------------------

def deterministic_check(eps: float = 0.04, df: float = 0.25):
    """
    Plot one noiseless trajectory at (eps, I = I_H1 + df*eps).
    Expect: v ≈ -1.3 descending along attracting branch, first hit is
    condition (b) w_new < w_f, R_hit ≈ 0.
    Prints outcome; no plot produced (headless).
    """
    fhn = FHN2D()
    I_H1 = fhn.I_hopf_lower_at(eps)
    I = I_H1 + df * eps
    lam = B_FHN * (I - I_FOLD_L)
    W_star = lam ** (2.0 / 3.0)
    w_f = I - 2.0 / 3.0
    W_init = 5.0 * W_star
    w0 = w_f + eps ** (2.0 / 3.0) * W_init
    v0 = v_left_root(w0, I)

    print(f"  Det. check: eps={eps}, I={I:.4f}, w_f={w_f:.4f}")
    print(f"    v0={v0:.4f}, w0={w0:.4f} (w0-w_f={w0-w_f:.4f})")

    dt = 5e-3
    v, w = v0, w0
    for step in range(int(400 / dt)):
        v_new = v + (v - v**3 / 3.0 - w + I) * dt
        w_new = w + eps * (v + A_FHN - B_FHN * w) * dt
        # condition (b): left fold crossing
        if w_new < w_f:
            W_hit_blowup = (w_new - w_f) / eps ** (2.0 / 3.0)
            R = W_hit_blowup / W_star
            print(f"    First hit (b) at t={step*dt:.2f}: w={w_new:.4f}, R={R:.4f}")
            return
        # condition (a): middle-saddle crossing while inside fold
        if w_new > w_f:
            vM = v_middle(np.array([w_new]), I)[0]
            if v < vM and v_new >= vM:
                W_hit_blowup = (w_new - w_f) / eps ** (2.0 / 3.0)
                R = W_hit_blowup / W_star
                print(f"    First hit (a) at t={step*dt:.2f}: w={w_new:.4f}, R={R:.4f}")
                return
        v, w = v_new, w_new
    print("    WARNING: no hit in T_max=400")


# ---------------------------------------------------------------------------
# Full-FHN integrator (vectorised Euler-Maruyama)
# ---------------------------------------------------------------------------

def simulate_fhn(
    eps: float,
    I: float,
    sigma: float,
    n_traj: int = 300,
    dt: float = 5e-3,
    T_max: float = 400.0,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """
    Simulate full FHN SDE, record w_hit at first canard peel-off event.

    Detector (left-fold side only; right-fold NOT a trigger):
      (a) v_new > v_M(w_new)  AND  w_new > w_f   — noise-driven crossing
      (b) w_new < w_f                              — deterministic fold descent

    For (b): w_hit is interpolated to w_f (the fold value) since the
    trajectory is at the canard peel-off point.
    For (a): w_hit is linearly interpolated to the v_M crossing.

    Returns
    -------
    w_hit : shape (n_traj,)  — NaN only if no event in T_max.
    """
    if rng is None:
        rng = np.random.default_rng()

    lam   = B_FHN * (I - I_FOLD_L)
    w_f   = I - 2.0 / 3.0
    W_star = lam ** (2.0 / 3.0)

    W_init = 5.0 * W_star
    w0 = w_f + eps ** (2.0 / 3.0) * W_init
    v0 = v_left_root(w0, I)

    v = np.full(n_traj, v0, dtype=float)
    w = np.full(n_traj, w0, dtype=float)
    hit   = np.zeros(n_traj, dtype=bool)
    w_hit = np.full(n_traj, np.nan)

    sqrt_dt = np.sqrt(dt)
    n_steps = int(T_max / dt)

    for _ in range(n_steps):
        if hit.all():
            break

        noise = rng.standard_normal(n_traj)
        active = ~hit

        v_new = np.copy(v)
        w_new = np.copy(w)
        v_new[active] = (v[active]
                         + (v[active] - v[active]**3 / 3.0 - w[active] + I) * dt
                         + sigma * sqrt_dt * noise[active])
        w_new[active] = (w[active]
                         + eps * (v[active] + A_FHN - B_FHN * w[active]) * dt)

        if not active.any():
            break

        act_idx = np.where(active)[0]

        # --- condition (b): left-fold descent ---
        # w drops below w_f; interpolate w_hit = w_f (fold value)
        b_mask = w_new[active] < w_f
        if b_mask.any():
            b_idx = act_idx[b_mask]
            # linear interp: w_hit = w_f at the step fraction where w crosses w_f
            w_old_b = w[b_idx]
            w_new_b = w_new[b_idx]
            frac_b = (w_f - w_old_b) / (w_new_b - w_old_b)
            frac_b = np.clip(frac_b, 0.0, 1.0)
            w_hit[b_idx] = w_old_b + frac_b * (w_new_b - w_old_b)  # ≈ w_f
            hit[b_idx] = True

        # --- condition (a): middle-saddle crossing above w_f ---
        # Only check trajectories still active AND w_new > w_f
        still_active = active & ~hit
        a_candidates = still_active & (w_new > w_f)
        if a_candidates.any():
            ac_idx = np.where(a_candidates)[0]
            vM_new = v_middle(w_new[a_candidates], I)
            a_mask = (v[a_candidates] < vM_new) & (v_new[a_candidates] >= vM_new)
            if a_mask.any():
                a_idx = ac_idx[a_mask]
                dv = v_new[a_idx] - v[a_idx]
                frac_a = np.where(
                    np.abs(dv) > 1e-15,
                    (vM_new[a_mask] - v[a_idx]) / dv,
                    0.5
                )
                frac_a = np.clip(frac_a, 0.0, 1.0)
                w_hit[a_idx] = w[a_idx] + frac_a * (w_new[a_idx] - w[a_idx])
                hit[a_idx] = True

        v, w = v_new, w_new

    return w_hit


# ---------------------------------------------------------------------------
# Sweep
# ---------------------------------------------------------------------------

def run_sweep(
    eps_vals=(0.04, 0.08),
    delta_factors=(0.0, 0.25, 0.5, 1.0),   # delta = df * eps  (not sqrt(eps))
    theta_vals=np.logspace(0.0, 1.5, 10),  # 1 to 31.6; transition at ~8
    n_traj=300,
    seed=42,
    verbose=True,
):
    fhn = FHN2D()
    rng_master = np.random.default_rng(seed)

    n_eps = len(eps_vals)
    n_I   = len(delta_factors)
    n_th  = len(theta_vals)

    W_hit_grid  = np.full((n_eps, n_I, n_th, n_traj), np.nan)
    lambda_grid = np.full((n_eps, n_I), np.nan)
    I_grid      = np.full((n_eps, n_I), np.nan)

    t0    = time.time()
    total = n_eps * n_I * n_th
    done  = 0

    for i, eps in enumerate(eps_vals):
        I_H1 = fhn.I_hopf_lower_at(eps)
        for j, df in enumerate(delta_factors):
            delta = df * eps          # FIX 1: delta ~ O(eps), not O(sqrt(eps))
            I     = I_H1 + delta
            lam   = B_FHN * (I - I_FOLD_L)
            sigma_star = np.sqrt(eps) * np.sqrt(lam)
            I_grid[i, j]      = I
            lambda_grid[i, j] = lam

            for k, theta in enumerate(theta_vals):
                sigma = theta * sigma_star
                rng   = np.random.default_rng(rng_master.integers(0, 2**63 - 1))
                w_hit = simulate_fhn(eps=eps, I=I, sigma=sigma,
                                     n_traj=n_traj, rng=rng)
                W_hit_grid[i, j, k, :] = w_hit
                done += 1
                if verbose:
                    w_f   = I - 2.0 / 3.0
                    W_star = lam ** (2.0 / 3.0)
                    R = (w_hit - w_f) / eps ** (2.0 / 3.0) / W_star
                    print(
                        f"  [{done:3d}/{total}] eps={eps:.2f} "
                        f"delta={df:.2f}*eps={delta:.4f} "
                        f"I={I:.4f} theta={theta:.3f}  "
                        f"fired={np.mean(~np.isnan(R)):.2f}  "
                        f"medR={np.nanmedian(R):+.3f}"
                    )

    if verbose:
        print(f"\n  total wall time: {time.time() - t0:.1f}s")

    return dict(
        eps_vals=np.asarray(eps_vals),
        delta_factors=np.asarray(delta_factors),
        theta_vals=np.asarray(theta_vals),
        I_vals=I_grid,
        lambda_vals=lambda_grid,
        W_hit=W_hit_grid,
    )


# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------

def summarise(data: dict):
    eps_vals    = data["eps_vals"]
    theta_vals  = data["theta_vals"]
    W_hit       = data["W_hit"]
    I_vals      = data["I_vals"]
    lambda_vals = data["lambda_vals"]

    bins_edges = np.array([-np.inf, 0.1, 0.5, 2.0, 4.0, np.inf])
    bin_labels = ["R<0.1", "0.1–0.5", "0.5–2", "2–4", "R>4"]

    per_theta = {}
    for k, theta in enumerate(theta_vals):
        R_all = []
        for i in range(len(eps_vals)):
            for j in range(W_hit.shape[1]):
                eps   = eps_vals[i]
                lam   = lambda_vals[i, j]
                W_star = lam ** (2.0 / 3.0)
                w_f   = I_vals[i, j] - 2.0 / 3.0
                R = (W_hit[i, j, k, :] - w_f) / eps ** (2.0 / 3.0) / W_star
                R_all.append(R[~np.isnan(R)])
        R_all = np.concatenate(R_all) if R_all else np.array([])

        if R_all.size > 0:
            med  = float(np.median(R_all))
            q25  = float(np.quantile(R_all, 0.25))
            q75  = float(np.quantile(R_all, 0.75))
            counts, _ = np.histogram(R_all, bins=bins_edges)
            bin_probs = counts / counts.sum()
        else:
            med = q25 = q75 = np.nan
            bin_probs = np.zeros(len(bin_labels))

        per_theta[float(theta)] = dict(
            n_fired=int(R_all.size),
            median=med, q25=q25, q75=q75,
            bin_probs=bin_probs,
        )

    return dict(per_theta=per_theta, bin_labels=bin_labels)


def collapse_check(data: dict):
    eps_vals    = data["eps_vals"]
    theta_vals  = data["theta_vals"]
    W_hit       = data["W_hit"]
    I_vals      = data["I_vals"]
    lambda_vals = data["lambda_vals"]

    rows = []
    for k, theta in enumerate(theta_vals):
        row = {"theta": theta}
        for i in range(len(eps_vals)):
            eps = eps_vals[i]
            for j in range(W_hit.shape[1]):
                lam   = lambda_vals[i, j]
                W_star = lam ** (2.0 / 3.0)
                w_f   = I_vals[i, j] - 2.0 / 3.0
                R = (W_hit[i, j, k, :] - w_f) / eps ** (2.0 / 3.0) / W_star
                R = R[~np.isnan(R)]
                key = f"e{eps:.2f}_l{lam:.4f}"
                row[key] = float(np.median(R)) if R.size else np.nan
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def make_figure(data: dict, summary: dict, out_path: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    theta_vals = data["theta_vals"]
    per        = summary["per_theta"]
    bin_labels = summary["bin_labels"]

    medians = np.array([per[float(t)]["median"] for t in theta_vals])
    q25s    = np.array([per[float(t)]["q25"]    for t in theta_vals])
    q75s    = np.array([per[float(t)]["q75"]    for t in theta_vals])

    fig, axes = plt.subplots(2, 1, figsize=(7.5, 8))

    ax = axes[0]
    ax.fill_between(theta_vals, q25s, q75s, alpha=0.25, color="C1",
                    label="IQR (25–75%)")
    ax.plot(theta_vals, medians, "o-", color="C1", lw=2, label="median $R_{\\rm hit}$")
    ax.axhline(0.0, color="k", lw=1.0, ls="-",  alpha=0.35)
    ax.axhline(1.0, color="k", lw=0.7, ls="--", alpha=0.45)
    ax.axvline(1.0, color="k", lw=0.7, ls="--", alpha=0.45)
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=0.5)
    ax.set_xlabel(r"$\Theta = \sigma / \sigma_*$,  "
                  r"$\sigma_* = \sqrt{\varepsilon}\,\lambda^{1/2}$")
    ax.set_ylabel(r"$R_{\rm hit} = W_{\rm hit}/W_*$")
    ax.set_title("Full-FHN canard-explosion regime map")
    ax.legend(loc="best", frameon=False, fontsize=9)
    ax.grid(True, which="both", alpha=0.25)

    ax = axes[1]
    bottom = np.zeros_like(theta_vals, dtype=float)
    colors = ["#2c7fb8", "#7fcdbb", "#cccccc", "#fdae6b", "#d94701"]
    for b, label in enumerate(bin_labels):
        probs = np.array([per[float(t)]["bin_probs"][b] for t in theta_vals])
        ax.bar(theta_vals, probs, bottom=bottom, width=theta_vals * 0.18,
               color=colors[b], label=label, edgecolor="white", lw=0.5)
        bottom += probs
    ax.set_xscale("log")
    ax.set_xlabel(r"$\Theta$")
    ax.set_ylabel(r"$P(R_{\rm hit}$ in bin$)$")
    ax.set_title("Stacked bin probabilities vs $\\Theta$")
    ax.set_ylim(0, 1)
    ax.legend(ncol=5, loc="lower center", bbox_to_anchor=(0.5, -0.32),
              frameon=False, fontsize=8)
    ax.grid(True, axis="y", alpha=0.25)

    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    out_data = os.path.join(_ROOT, "data")
    out_fig  = os.path.join(_ROOT, "figures")
    os.makedirs(out_data, exist_ok=True)
    os.makedirs(out_fig, exist_ok=True)

    print("\n=== Full-FHN canard-explosion regime map ===\n")
    print(f"  I_FOLD_L = {I_FOLD_L:.6f}\n")

    # --- FIX 4: deterministic sanity check ---
    print("Deterministic sanity check (sigma=0):")
    deterministic_check(eps=0.04, df=0.25)
    print()

    # --- full sweep ---
    data = run_sweep(
        eps_vals=(0.04, 0.08),
        delta_factors=(0.0, 0.25, 0.5, 1.0),
        theta_vals=np.logspace(0.0, 1.5, 10),
        n_traj=300,
        seed=42,
    )

    np.savez(
        os.path.join(out_data, "canard_full_fhn_hit_location.npz"),
        eps_vals=data["eps_vals"],
        I_vals=data["I_vals"],
        lambda_vals=data["lambda_vals"],
        theta_vals=data["theta_vals"],
        W_hit=data["W_hit"],
    )
    print("  data saved -> data/canard_full_fhn_hit_location.npz")

    summary  = summarise(data)
    fig_path = os.path.join(out_fig, "canard_full_fhn_regime_map.png")
    make_figure(data, summary, fig_path)
    print(f"  figure saved -> figures/canard_full_fhn_regime_map.png\n")

    per        = summary["per_theta"]
    theta_vals = data["theta_vals"]

    print("§10 Pooled summary (across all eps and lambda):")
    print(f"  {'Theta':>8s}  {'n_fired':>8s}  {'median R':>10s}  {'IQR':>24s}")
    for theta in theta_vals:
        d = per[float(theta)]
        print(
            f"  {theta:8.3f}  {d['n_fired']:8d}  {d['median']:+10.3f}  "
            f"[{d['q25']:+8.3f}, {d['q75']:+8.3f}]"
        )

    print("\n§10 Collapse check — median R at each (eps, lambda):")
    rows     = collapse_check(data)
    col_keys = [k for k in rows[0] if k != "theta"]
    header   = "  ".join(f"{k:>16s}" for k in col_keys)
    print(f"  {'Theta':>8s}  {header}")
    for row in rows:
        vals = "  ".join(
            f"{row[k]:+10.3f}      " if not np.isnan(row[k]) else "     NaN      "
            for k in col_keys
        )
        print(f"  {row['theta']:8.3f}  {vals}")


if __name__ == "__main__":
    main()
