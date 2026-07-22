"""
run_sweep.py  —  Main driver for the FHN regime map computation.

Usage:
    python run_sweep.py [--quick] [--n-traj N] [--seed S] [--delta-ratio R]

--quick        : 8×8 coarse grid for fast testing (~3 min)
--n-traj N     : trajectories per grid point (default 500)
--delta-ratio R: log-ratio failure threshold (default 2.0)

Speed notes
-----------
T is adaptive per grid point: min(8 × T_drift(eps), T_max).
At eps=0.25, T_drift≈9.6 → simulates to T=77.  At eps=0.01, T_drift≈239 → T=400 cap.
dt=5e-3 (vs 1e-3 previously) gives 5× speedup with negligible MFPT error.
Early exit at 90% fired per grid point.
"""

import argparse, os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))
from sweep import sweep_grid, failure_boundary, fit_power_law
from plot  import plot_regime_map, plot_mfpt_comparison, plot_boundary_scaling

os.makedirs("figures", exist_ok=True)
os.makedirs("data",    exist_ok=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick",       action="store_true")
    parser.add_argument("--n-traj",      type=int,   default=500)
    parser.add_argument("--seed",        type=int,   default=42)
    parser.add_argument("--delta-ratio", type=float, default=2.0)
    args = parser.parse_args()

    if args.quick:
        sigma_vals = np.linspace(0.05, 0.45, 8)
        eps_vals   = np.linspace(0.01, 0.25, 8)
        T_max      = 200.0
        print("[Quick mode] 8×8 grid  (adaptive T, max 200)")
    else:
        sigma_vals = np.linspace(0.05, 0.45, 20)
        eps_vals   = np.linspace(0.01, 0.25, 15)
        T_max      = 400.0
        print(f"[Full mode] {len(sigma_vals)}×{len(eps_vals)} grid  (adaptive T, max {T_max})")

    print(f"sigma: {sigma_vals[0]:.3f} – {sigma_vals[-1]:.3f}  ({len(sigma_vals)} pts)")
    print(f"eps:   {eps_vals[0]:.3f} – {eps_vals[-1]:.3f}   ({len(eps_vals)} pts)")
    print(f"N_traj={args.n_traj},  delta_ratio={args.delta_ratio},  seed={args.seed}")
    print()

    result = sweep_grid(
        sigma_vals=sigma_vals, eps_vals=eps_vals,
        n_trajectories=args.n_traj, T=T_max,
        seed=args.seed, verbose=True,
    )

    np.savez("data/sweep_result.npz",
             sigma_vals=sigma_vals, eps_vals=eps_vals,
             mfpt_full=result["mfpt_full"], mfpt_reduced=result["mfpt_reduced"],
             log_ratio=result["log_ratio"], rel_error=result["rel_error"],
             drift_times=result["drift_times"])
    print("Saved: data/sweep_result.npz")

    eps_bd, sigma_bd = failure_boundary(result, delta_ratio=args.delta_ratio)
    alpha_fit, C_fit = fit_power_law(eps_bd, sigma_bd)

    if not np.isnan(alpha_fit):
        print(f"\nFitted boundary:  sigma = {C_fit:.3f} * eps^{alpha_fit:.3f}")
    else:
        print("\nNot enough boundary points to fit — try a wider sigma range or more trajectories.")
    print(f"BG prediction:    sigma ~ eps^1.5")

    fig1 = plot_regime_map(result, delta_ratio=args.delta_ratio,
                           boundary_eps=eps_bd, boundary_sigma=sigma_bd,
                           bg_C=C_fit if not np.isnan(C_fit) else 1.0,
                           save_path="figures/regime_map.png")
    print("Saved: figures/regime_map.png")

    fig2 = plot_mfpt_comparison(result, save_path="figures/mfpt_comparison.png")
    print("Saved: figures/mfpt_comparison.png")

    if len(eps_bd) >= 3:
        fig3 = plot_boundary_scaling(eps_bd, sigma_bd, alpha_fit, C_fit,
                                     save_path="figures/boundary_scaling.png")
        print("Saved: figures/boundary_scaling.png")

    print("\nDone.")


if __name__ == "__main__":
    main()
