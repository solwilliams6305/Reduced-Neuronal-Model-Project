"""
canard_hit_location.py
----------------------
Autonomous canard hit-location regime map (Roadmap §4).

Tests the corrected blow-up prediction
    sigma_* = C_q * sqrt(eps) * lambda^{1/2},
    W_*     = lambda^{2/3},
    eta_*   = lambda^{1/2}.

Simulates full FHN with degenerate noise at fixed I in the canard window,
starting at the deterministic fixed point. At the first time the trajectory
crosses v = 0 ("committed to spike"), records w_hit, then maps to blow-up
coordinate W_hit_blowup = (w_hit - w_f) / eps^{2/3} and finally to the
scaled coordinate

    R_hit = W_hit_blowup / W_*,
    Theta = sigma / (sqrt(eps) * lambda^{1/2}).

Produces:
    figures/canard_regime_map.png   — median(R) vs Theta with IQR
                                     and stacked bin probabilities.
    data/canard_hit_location.npz    — raw (Theta_grid, R_hit array).

Pure numpy (no scipy needed for the simulator itself). Imports kernel.py
through _shim.py so the same script runs in your local environment and
in sandboxes without scipy.
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

# Make this script runnable from any cwd, and shim scipy if needed
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _HERE)
sys.path.insert(0, _ROOT)

import _shim  # noqa: F401, installs brentq shim if scipy is missing
from kernel import FHN2D

# ---------------------------------------------------------------------------
# Model parameters
# ---------------------------------------------------------------------------

A = 0.7
B = 0.8


def lambda_phys(I: float) -> float:
    """Blow-up slow-drift parameter lambda = b (I - I_fold_L)."""
    I_fold_L = (A - 1.0 + 2.0 * B / 3.0) / B   # 0.29166...
    return B * (I - I_fold_L)


def W_star(lam: float) -> float:
    return float(lam) ** (2.0 / 3.0)


def Theta(sigma: float, eps: float, lam: float) -> float:
    return float(sigma) / (np.sqrt(eps) * np.sqrt(lam))


# ---------------------------------------------------------------------------
# Vectorised Euler-Maruyama, first-hit at v = v_cross
# ---------------------------------------------------------------------------

def first_hit_W(
    eps: float,
    sigma: float,
    I: float,
    n_traj: int = 200,
    T_max: float = 600.0,
    dt: float = 5e-3,
    v_cross: float = -0.5,
    W_init_over_star: float = 5.0,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """
    Simulate n_traj trajectories of FHN starting on the attracting slow
    manifold ABOVE the canard window (W_init = W_init_over_star * W_*),
    so the deterministic trajectory descends THROUGH the canard window
    toward the fold. Return w at first time v crosses v_cross upward
    (NaN if never).

    v_cross = -0.5 catches the peel-off from the canard structure
    (just past the canard exit near v = -1) before w has drifted during
    the relaxation phase.

    The initial condition matters: starting at the lower-resonator FP
    puts the trajectory at the blow-up origin and there is no canard
    descent to measure. Starting at W_init = 5 W_* gives the trajectory
    a clean attracting-manifold descent through the canard window.
    """
    if rng is None:
        rng = np.random.default_rng()

    m = FHN2D(I=I, a=A, b=B)
    w_f = m.w_fold_left
    lam = lambda_phys(I)
    W_init_blowup = W_init_over_star * (lam ** (2.0 / 3.0))
    # Original-coord initial condition on attracting branch
    w0 = w_f + (eps ** (2.0 / 3.0)) * W_init_blowup
    # On the attracting branch v = -1 - sqrt(w - w_f) (leading order)
    v0 = -1.0 - np.sqrt(max(w0 - w_f, 0.0))

    v = np.full(n_traj, v0, dtype=float)
    w = np.full(n_traj, w0, dtype=float)
    hit = np.zeros(n_traj, dtype=bool)
    w_hit = np.full(n_traj, np.nan)

    sqrt_dt = np.sqrt(dt)
    n_steps = int(T_max / dt)

    for _ in range(n_steps):
        # Vectorised Euler-Maruyama step
        noise = rng.standard_normal(n_traj)
        dv = (v - v ** 3 / 3.0 - w + I) * dt + sigma * sqrt_dt * noise
        dw_inc = eps * (v + A - B * w) * dt
        v_new = v + dv
        w_new = w + dw_inc

        new_hit = (~hit) & (v < v_cross) & (v_new >= v_cross)
        if new_hit.any():
            # Linear interpolation across the crossing for sub-step accuracy
            frac = (v_cross - v[new_hit]) / (v_new[new_hit] - v[new_hit])
            w_hit[new_hit] = w[new_hit] + frac * (w_new[new_hit] - w[new_hit])
            hit |= new_hit
            if hit.all():
                break

        v, w = v_new, w_new

    # Blow-up W coordinate
    W_hit = (w_hit - w_f) / (eps ** (2.0 / 3.0))
    return W_hit


# ---------------------------------------------------------------------------
# Sweep over (eps, lambda, sigma)
# ---------------------------------------------------------------------------

def run_sweep(
    eps_vals=(0.04, 0.08),
    lam_vals=(0.012, 0.024),       # via I = I_fold_L + lam/B
    theta_vals=np.logspace(-0.7, 0.9, 8),   # Theta in [~0.2, ~8]
    n_traj=300,
    seed=42,
    verbose=True,
):
    """
    For each (eps, lam, Theta) compute n_traj trajectories and return
    a 4-D array of W_hit shape (n_eps, n_lam, n_theta, n_traj).
    """
    I_fold_L = (A - 1.0 + 2.0 * B / 3.0) / B
    rng_master = np.random.default_rng(seed)

    grid = np.full(
        (len(eps_vals), len(lam_vals), len(theta_vals), n_traj),
        np.nan,
    )

    t0 = time.time()
    total = len(eps_vals) * len(lam_vals) * len(theta_vals)
    done = 0

    for i, eps in enumerate(eps_vals):
        for j, lam in enumerate(lam_vals):
            I = I_fold_L + lam / B
            for k, theta in enumerate(theta_vals):
                sigma = theta * np.sqrt(eps) * np.sqrt(lam)
                rng = np.random.default_rng(rng_master.integers(0, 2**63 - 1))
                W_hit = first_hit_W(
                    eps=eps, sigma=sigma, I=I,
                    n_traj=n_traj, rng=rng,
                )
                grid[i, j, k, :] = W_hit
                done += 1
                if verbose:
                    frac_fired = float(np.mean(~np.isnan(W_hit)))
                    print(
                        f"  [{done:3d}/{total}] eps={eps:.3f} lam={lam:.3f} "
                        f"Theta={theta:.3f} sigma={sigma:.4f} "
                        f"fired={frac_fired:.2f}"
                    )

    if verbose:
        print(f"  total wall time: {time.time() - t0:.1f}s")

    return {
        "eps_vals": np.asarray(eps_vals),
        "lam_vals": np.asarray(lam_vals),
        "theta_vals": np.asarray(theta_vals),
        "W_hit": grid,
    }


# ---------------------------------------------------------------------------
# Analysis: R_hit summaries and bin probabilities
# ---------------------------------------------------------------------------

def summarise(data: dict):
    """
    Build (Theta, R_hit) pairs across all (eps, lam) and return:
      Theta_flat, R_flat — concatenated arrays for scatter / boxplot
      per_theta: dict of theta -> dict(median, q25, q75, bins, n_fired, n)
    """
    eps_vals = data["eps_vals"]
    lam_vals = data["lam_vals"]
    theta_vals = data["theta_vals"]
    W = data["W_hit"]   # (n_eps, n_lam, n_theta, n_traj)

    bins_edges = np.array([0.0, 0.1, 0.5, 2.0, 4.0, np.inf])
    bin_labels = ["R<0.1", "0.1-0.5", "0.5-2", "2-4", "R>4"]

    per_theta = {}
    Theta_flat = []
    R_flat = []

    for k, theta in enumerate(theta_vals):
        R_all = []
        for i in range(len(eps_vals)):
            for j in range(len(lam_vals)):
                lam = lam_vals[j]
                W_star = lam ** (2.0 / 3.0)
                W_k = W[i, j, k, :]
                R = W_k / W_star
                # Keep only fired trajectories
                R = R[~np.isnan(R)]
                R_all.append(R)
        R_all = np.concatenate(R_all) if R_all else np.array([])

        if R_all.size > 0:
            med = float(np.median(R_all))
            q25 = float(np.quantile(R_all, 0.25))
            q75 = float(np.quantile(R_all, 0.75))
        else:
            med = q25 = q75 = np.nan

        # Bin probabilities
        bin_probs = np.zeros(len(bin_labels))
        if R_all.size > 0:
            counts, _ = np.histogram(R_all, bins=bins_edges)
            bin_probs = counts / counts.sum()

        per_theta[float(theta)] = dict(
            n_total=int(W.shape[3] * len(eps_vals) * len(lam_vals)),
            n_fired=int(R_all.size),
            median=med, q25=q25, q75=q75,
            bin_probs=bin_probs, bin_labels=bin_labels,
            R_values=R_all,
        )
        Theta_flat.append(np.full(R_all.size, theta))
        R_flat.append(R_all)

    return {
        "Theta_flat": np.concatenate(Theta_flat) if Theta_flat else np.array([]),
        "R_flat":     np.concatenate(R_flat) if R_flat else np.array([]),
        "per_theta":  per_theta,
        "bin_labels": bin_labels,
    }


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def make_figure(data: dict, summary: dict, out_path: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    theta_vals = data["theta_vals"]
    per = summary["per_theta"]
    bin_labels = summary["bin_labels"]

    medians = np.array([per[float(t)]["median"] for t in theta_vals])
    q25s    = np.array([per[float(t)]["q25"] for t in theta_vals])
    q75s    = np.array([per[float(t)]["q75"] for t in theta_vals])
    n_fired = np.array([per[float(t)]["n_fired"] for t in theta_vals])

    fig, axes = plt.subplots(2, 1, figsize=(7, 8))

    # Panel 1: median R_hit vs Theta with IQR shading
    ax = axes[0]
    ax.fill_between(theta_vals, q25s, q75s, alpha=0.25, color="C0",
                    label="IQR (25–75%)")
    ax.plot(theta_vals, medians, "o-", color="C0", lw=2, label="median R_hit")
    ax.axhline(1.0, color="k", lw=0.7, ls="--", alpha=0.6)
    ax.axvline(1.0, color="k", lw=0.7, ls="--", alpha=0.6)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\Theta = \sigma / (\sqrt{\varepsilon}\,\lambda^{1/2})$")
    ax.set_ylabel(r"$R_{\rm hit} = W_{\rm hit} / W_*$, $W_* = \lambda^{2/3}$")
    ax.set_title("Autonomous canard hit-location regime map")
    ax.legend(loc="best", frameon=False, fontsize=9)
    ax.grid(True, which="both", alpha=0.25)

    # Panel 2: stacked bin probabilities vs Theta
    ax = axes[1]
    bottom = np.zeros_like(theta_vals, dtype=float)
    colors = ["#1f77b4", "#9ecae1", "#cccccc", "#fdae6b", "#e6550d"]
    for b, label in enumerate(bin_labels):
        probs = np.array([per[float(t)]["bin_probs"][b] for t in theta_vals])
        ax.bar(theta_vals, probs, bottom=bottom, width=theta_vals * 0.18,
               color=colors[b], label=label, edgecolor="white", lw=0.5)
        bottom += probs
    ax.set_xscale("log")
    ax.set_xlabel(r"$\Theta$")
    ax.set_ylabel("P(bin)")
    ax.set_title("Stacked probability of $R_{\\rm hit}$ bin vs $\\Theta$")
    ax.set_ylim(0, 1)
    ax.legend(ncol=5, loc="lower center", bbox_to_anchor=(0.5, -0.32),
              frameon=False, fontsize=8)
    ax.grid(True, axis="y", alpha=0.25)

    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

def main():
    print("\n=== Autonomous canard hit-location regime map ===\n")
    print("Roadmap §4: build R_hit vs Theta with median, IQR, and binned probabilities.\n")

    out_dir_data = os.path.join(_ROOT, "data")
    out_dir_fig = os.path.join(_ROOT, "figures")
    os.makedirs(out_dir_data, exist_ok=True)
    os.makedirs(out_dir_fig, exist_ok=True)

    # Compact but informative sweep (~ a few minutes)
    data = run_sweep(
        eps_vals=(0.04, 0.08),
        lam_vals=(0.012, 0.024),
        theta_vals=np.logspace(-0.7, 0.9, 8),
        n_traj=300,
        seed=42,
    )

    np.savez(
        os.path.join(out_dir_data, "canard_hit_location.npz"),
        eps_vals=data["eps_vals"],
        lam_vals=data["lam_vals"],
        theta_vals=data["theta_vals"],
        W_hit=data["W_hit"],
    )
    print(f"  data saved -> data/canard_hit_location.npz")

    summary = summarise(data)
    fig_path = os.path.join(out_dir_fig, "canard_regime_map.png")
    make_figure(data, summary, fig_path)
    print(f"  figure saved -> figures/canard_regime_map.png")

    # Quick table for the conversation
    per = summary["per_theta"]
    print("\nTheta sweep summary:")
    print(f"  {'Theta':>8s}  {'n_fired':>8s}  {'median R':>10s}  "
          f"{'IQR':>16s}  {'P(R<0.1)':>9s}  {'P(0.5-2)':>9s}  {'P(R>4)':>9s}")
    for theta in data["theta_vals"]:
        d = per[float(theta)]
        bp = d["bin_probs"]
        print(f"  {theta:8.3f}  {d['n_fired']:8d}  {d['median']:10.3f}  "
              f"[{d['q25']:6.3f}, {d['q75']:6.3f}]  "
              f"{bp[0]:9.3f}  {bp[2]:9.3f}  {bp[4]:9.3f}")


if __name__ == "__main__":
    main()
