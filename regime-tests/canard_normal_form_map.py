"""
canard_normal_form_map.py
-------------------------
Autonomous canard hit-location regime map on the blow-up normal form
(Roadmap §2, §4).

The normal-form system is

    dV = (V^2 - W) dT + eta dB_T,
    dW = -lambda dT,

with initial condition on the attracting slow manifold V = -sqrt(W) at
W_init = (W_init_over_star) * W_*, W_* = lambda^{2/3}.

The trajectory descends along V = -sqrt(W) toward W = 0, then on the
maximal canard would continue onto the repelling branch V = +sqrt(W).
Noise pushes the trajectory off the canard. We record W at the first
time V crosses V_cross from below.

Scaled variables:
    Theta = eta / lambda^{1/2}   (eta_* = lambda^{1/2})
    R_hit = W_hit / W_*

This is the direct test of the analytical prediction
    eta_* = C_q * lambda^{1/2}.

The companion full-FHN test (in the canard-explosion window I near I_H1)
is a separate experiment with extra modelling complexity.
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np


# ---------------------------------------------------------------------------
# Normal-form integrator (vectorised Euler-Maruyama)
# ---------------------------------------------------------------------------

def simulate_normal_form(
    lam: float,
    eta: float,
    n_traj: int = 600,
    W_init_over_star: float = 5.0,
    T_max: float = 80.0,
    dT: float = 5e-4,
    V_cross: float = +1.0,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """
    Simulate the canonical fold-blow-up SDE forward in T.

    Returns
    -------
    W_hit : shape (n_traj,), W-coordinate at first time V crosses V_cross
            from below. NaN if no crossing occurred.
    """
    if rng is None:
        rng = np.random.default_rng()

    W_star = lam ** (2.0 / 3.0)
    W_init = W_init_over_star * W_star
    V_init = -np.sqrt(W_init)

    V = np.full(n_traj, V_init, dtype=float)
    W = np.full(n_traj, W_init, dtype=float)
    hit = np.zeros(n_traj, dtype=bool)
    W_hit = np.full(n_traj, np.nan)

    sqrt_dT = np.sqrt(dT)
    n_steps = int(T_max / dT)

    for _ in range(n_steps):
        noise = rng.standard_normal(n_traj)
        V_new = V + (V * V - W) * dT + eta * sqrt_dT * noise
        W_new = W - lam * dT

        new_hit = (~hit) & (V < V_cross) & (V_new >= V_cross)
        if new_hit.any():
            frac = (V_cross - V[new_hit]) / (V_new[new_hit] - V[new_hit])
            W_hit[new_hit] = W[new_hit] + frac * (W_new[new_hit] - W[new_hit])
            hit |= new_hit
            if hit.all():
                break

        V, W = V_new, W_new

    return W_hit


# ---------------------------------------------------------------------------
# Sweep over (lambda, Theta)
# ---------------------------------------------------------------------------

def run_sweep(
    lam_vals=(0.005, 0.01, 0.02, 0.04),
    theta_vals=np.logspace(-0.8, 1.0, 10),
    n_traj=600,
    seed=42,
    verbose=True,
):
    rng_master = np.random.default_rng(seed)
    grid = np.full((len(lam_vals), len(theta_vals), n_traj), np.nan)

    t0 = time.time()
    total = len(lam_vals) * len(theta_vals)
    done = 0

    for j, lam in enumerate(lam_vals):
        for k, theta in enumerate(theta_vals):
            eta = theta * np.sqrt(lam)
            rng = np.random.default_rng(rng_master.integers(0, 2**63 - 1))
            W_hit = simulate_normal_form(
                lam=lam, eta=eta, n_traj=n_traj, rng=rng,
            )
            grid[j, k, :] = W_hit
            done += 1
            if verbose:
                fired = np.mean(~np.isnan(W_hit))
                print(
                    f"  [{done:3d}/{total}] lam={lam:.4f} theta={theta:.3f} "
                    f"eta={eta:.4f}  fired={fired:.2f}"
                )

    if verbose:
        print(f"  total wall time: {time.time() - t0:.1f}s")

    return dict(
        lam_vals=np.asarray(lam_vals),
        theta_vals=np.asarray(theta_vals),
        W_hit=grid,
    )


# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------

def summarise(data: dict):
    lam_vals = data["lam_vals"]
    theta_vals = data["theta_vals"]
    W = data["W_hit"]   # (n_lam, n_theta, n_traj)

    bins_edges = np.array([-np.inf, 0.1, 0.5, 2.0, 4.0, np.inf])
    bin_labels = ["R<0.1", "0.1–0.5", "0.5–2", "2–4", "R>4"]

    per_theta = {}
    for k, theta in enumerate(theta_vals):
        R_all = []
        for j in range(len(lam_vals)):
            W_star = lam_vals[j] ** (2.0 / 3.0)
            R = W[j, k, :] / W_star
            R = R[~np.isnan(R)]
            R_all.append(R)
        R_all = np.concatenate(R_all) if R_all else np.array([])

        if R_all.size > 0:
            med = float(np.median(R_all))
            q25 = float(np.quantile(R_all, 0.25))
            q75 = float(np.quantile(R_all, 0.75))
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


# ---------------------------------------------------------------------------
# Figure
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

    fig, axes = plt.subplots(2, 1, figsize=(7.5, 8))

    ax = axes[0]
    ax.fill_between(theta_vals, q25s, q75s, alpha=0.25, color="C0",
                    label="IQR (25–75%)")
    ax.plot(theta_vals, medians, "o-", color="C0", lw=2, label="median R_hit")
    # Reference lines
    ax.axhline(1.0, color="k", lw=0.7, ls="--", alpha=0.5)
    ax.axvline(1.0, color="k", lw=0.7, ls="--", alpha=0.5)
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=0.5)
    ax.set_xlabel(r"$\Theta = \eta / \lambda^{1/2}$")
    ax.set_ylabel(r"$R_{\rm hit} = W_{\rm hit}/W_*$,  $W_* = \lambda^{2/3}$")
    ax.set_title("Canard hit-location regime map (blow-up normal form)")
    ax.legend(loc="best", frameon=False, fontsize=9)
    ax.grid(True, which="both", alpha=0.25)

    ax = axes[1]
    bottom = np.zeros_like(theta_vals, dtype=float)
    colors = ["#2c7fb8", "#7fcdbb", "#cccccc", "#fdae6b", "#d94701"]
    for b, label in enumerate(bin_labels):
        probs = np.array([per[float(t)]["bin_probs"][b] for t in theta_vals])
        ax.bar(theta_vals, probs, bottom=bottom, width=theta_vals * 0.16,
               color=colors[b], label=label, edgecolor="white", lw=0.5)
        bottom += probs
    ax.set_xscale("log")
    ax.set_xlabel(r"$\Theta$")
    ax.set_ylabel(r"$P(R_{\rm hit}$ in bin$)$")
    ax.set_title("Stacked probability of $R_{\\rm hit}$ bins vs $\\Theta$")
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
    _HERE = os.path.dirname(os.path.abspath(__file__))
    _ROOT = os.path.dirname(_HERE)
    out_data = os.path.join(_ROOT, "data")
    out_fig = os.path.join(_ROOT, "figures")
    os.makedirs(out_data, exist_ok=True)
    os.makedirs(out_fig, exist_ok=True)

    print("\n=== Canard normal-form hit-location regime map ===\n")
    print("Tests eta_* = C_q * lambda^{1/2} directly on the blow-up SDE.\n")

    data = run_sweep(
        lam_vals=(0.005, 0.01, 0.02, 0.04),
        theta_vals=np.logspace(-0.8, 1.0, 10),
        n_traj=600,
        seed=42,
    )

    np.savez(
        os.path.join(out_data, "canard_normal_form_hit_location.npz"),
        lam_vals=data["lam_vals"],
        theta_vals=data["theta_vals"],
        W_hit=data["W_hit"],
    )
    print("  data saved -> data/canard_normal_form_hit_location.npz")

    summary = summarise(data)
    fig_path = os.path.join(out_fig, "canard_normal_form_regime_map.png")
    make_figure(data, summary, fig_path)
    print(f"  figure saved -> figures/canard_normal_form_regime_map.png")

    per = summary["per_theta"]
    print("\nTheta sweep summary (pooled across lambda):")
    print(f"  {'Theta':>8s}  {'n_fired':>8s}  {'median R':>10s}  "
          f"{'IQR':>20s}")
    for theta in data["theta_vals"]:
        d = per[float(theta)]
        print(
            f"  {theta:8.3f}  {d['n_fired']:8d}  {d['median']:+10.3f}  "
            f"[{d['q25']:+6.3f}, {d['q75']:+6.3f}]"
        )

    # Collapse check: lambda-independence of R_hit vs Theta at fixed Theta
    print("\nCollapse check — median R at each (lambda, Theta):")
    print(f"  {'Theta':>8s}  " + "  ".join(
        [f"lam={l:.4f}" for l in data['lam_vals']]
    ))
    for k, theta in enumerate(data["theta_vals"]):
        row = []
        for j, lam in enumerate(data["lam_vals"]):
            W_star = lam ** (2.0 / 3.0)
            R = data["W_hit"][j, k, :] / W_star
            R = R[~np.isnan(R)]
            if R.size:
                row.append(f"{np.median(R):+.3f}")
            else:
                row.append("  NaN  ")
        print(f"  {theta:8.3f}  " + "  ".join(row))


if __name__ == "__main__":
    main()
