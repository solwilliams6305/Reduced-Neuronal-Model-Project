#!/usr/bin/env python3
"""
run_finer_exposure_sweep.py

Runs a finer stochastic FHN sweep concentrated around the square-root exposure
boundary

    sigma ≈ C sqrt(eps)

then re-extracts the log_ratio=0 contour and refits:

    sigma = C sqrt(eps)
    sigma = C eps^p
    sigma^2 = alpha + beta eps

This is designed to test whether the exponent near 1/2 survives higher
resolution around the failure boundary.

Usage
-----
Quick-ish:
    python run_finer_exposure_sweep.py --n-traj 1200

Better:
    python run_finer_exposure_sweep.py --n-traj 2500 --eps-n 28 --sigma-band 0.09 --sigma-n 34

Then plot/analyse using:
    python plot_exposure_scaling_boundary.py --sweep-path data/sweep_finer_exposure.npz --tau 1.1615
"""

from __future__ import annotations

import argparse
import os
import sys
import multiprocessing

import numpy as np

try:
    multiprocessing.set_start_method("fork", force=True)
except RuntimeError:
    pass

sys.path.insert(0, os.path.dirname(__file__))

from sweep import sweep_grid
from simulate import reduced_drift_time


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="Run finer FHN sweep near sqrt exposure boundary.")

parser.add_argument("--out", type=str, default="data/sweep_finer_exposure.npz",
                    help="Output .npz path.")
parser.add_argument("--n-traj", type=int, default=1500,
                    help="Trajectories per grid point.")
parser.add_argument("--seed", type=int, default=123)
parser.add_argument("--eps-min", type=float, default=0.015)
parser.add_argument("--eps-max", type=float, default=0.17)
parser.add_argument("--eps-n", type=int, default=24)
parser.add_argument("--sigma-n", type=int, default=30)
parser.add_argument("--sigma-band", type=float, default=0.075,
                    help="Half-width of sigma band around C sqrt(eps).")
parser.add_argument("--C", type=float, default=1.1002,
                    help="Initial square-root exposure prefactor.")
parser.add_argument("--sigma-min", type=float, default=0.08)
parser.add_argument("--sigma-max", type=float, default=0.48)
parser.add_argument("--T", type=float, default=400.0,
                    help="Simulation horizon passed to sweep_grid.")
parser.add_argument("--global-grid", action="store_true",
                    help=(
                        "Use one rectangular sigma grid covering the full band. "
                        "This is compatible with simple heatmap plotting. "
                        "Default is also rectangular, but chosen from union of local bands."
                    ))

args = parser.parse_args()


# ---------------------------------------------------------------------------
# Construct refined grid
# ---------------------------------------------------------------------------

I, a, b = -0.1, 0.7, 0.8

os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

# Put slightly more resolution where previous contour existed.
eps_vals = np.linspace(args.eps_min, args.eps_max, args.eps_n)

# Build a rectangular sigma grid covering the local band around C sqrt(eps).
sigma_center = args.C * np.sqrt(eps_vals)
sigma_low = np.clip(sigma_center - args.sigma_band, args.sigma_min, args.sigma_max)
sigma_high = np.clip(sigma_center + args.sigma_band, args.sigma_min, args.sigma_max)

smin = float(np.min(sigma_low))
smax = float(np.max(sigma_high))

# Add a small safety margin.
smin = max(args.sigma_min, smin - 0.015)
smax = min(args.sigma_max, smax + 0.015)

sigma_vals = np.linspace(smin, smax, args.sigma_n)

print("\nFiner exposure sweep")
print(f"  eps range:   {eps_vals.min():.4f} to {eps_vals.max():.4f}  n={len(eps_vals)}")
print(f"  sigma range: {sigma_vals.min():.4f} to {sigma_vals.max():.4f}  n={len(sigma_vals)}")
print(f"  grid points: {len(eps_vals) * len(sigma_vals)}")
print(f"  n_traj:      {args.n_traj}")
print(f"  output:      {args.out}")

print("\nPredicted sqrt centres:")
for e, s in zip(eps_vals[::max(1, len(eps_vals)//8)], sigma_center[::max(1, len(eps_vals)//8)]):
    print(f"  eps={e:.4f}, C sqrt(eps)={s:.4f}, T_drift≈{reduced_drift_time(eps=float(e), I=I, a=a, b=b):.2f}")

print("\nRunning sweep...")

result = sweep_grid(
    sigma_vals,
    eps_vals,
    n_trajectories=args.n_traj,
    T=args.T,
    seed=args.seed,
    verbose=True,
)

np.savez(
    args.out,
    sigma_vals=result["sigma_vals"],
    eps_vals=result["eps_vals"],
    mfpt_full=result["mfpt_full"],
    mfpt_reduced=result["mfpt_reduced"],
    log_ratio=result["log_ratio"],
    rel_error=result["rel_error"],
    drift_times=result["drift_times"],
)

print(f"\nSaved: {args.out}")
print("\nNext run:")
print(f"  python plot_exposure_scaling_boundary.py --sweep-path {args.out} --tau 1.1615")
