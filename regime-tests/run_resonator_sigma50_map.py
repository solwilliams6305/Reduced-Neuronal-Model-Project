#!/usr/bin/env python3
"""
run_resonator_sigma50_map.py

Batch driver + analysis for the stochastic FHN resonator regime.

What it does
------------
For each I in a resonator-window grid, it calls
`verify_resonator_phase_difference.py` with a sigma sweep, then fits

    P(excursion before T | sigma) ≈ 1 / (1 + exp(-(sigma - sigma50)/width))

and saves:
  - all_sigma_sweeps.csv
  - sigma50_map_summary.csv
  - excursion_probability_curves_by_I.png
  - sigma50_vs_I.png
  - transition_width_vs_I.png
  - escape_phase_R_by_I.png
  - median_rotations_by_I.png

Recommended first run
---------------------
python run_resonator_sigma50_map.py \
  --verify-script verify_resonator_phase_difference.py \
  --I-values 0.295 0.300 0.305 0.310 0.315 0.320 \
  --sigmas 0.005 0.010 0.015 0.020 0.025 0.030 0.035 0.040 0.050 0.060 \
  --N 500 --T 250 --dt 0.005

Notes
-----
This script is intentionally a wrapper around the existing resonator script,
so you do not have to duplicate the simulation logic. It turns a single-I
sigma sweep into the actual object we need: sigma50(I), the finite-time
resonator failure boundary.
"""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    from scipy.optimize import curve_fit
    SCIPY_AVAILABLE = True
except Exception:
    SCIPY_AVAILABLE = False


def logistic(sigma: np.ndarray, sigma50: float, width: float) -> np.ndarray:
    width = max(float(width), 1e-9)
    z = np.clip((sigma - sigma50) / width, -80, 80)
    return 1.0 / (1.0 + np.exp(-z))


def safe_I_label(I: float) -> str:
    s = f"I_{I:.6f}".replace("-", "m").replace(".", "p")
    return s


def ensure_dir(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def run_one_I(args: argparse.Namespace, I: float, outdir: Path) -> Path:
    """Run the existing verifier for one I value and return the summary CSV path."""
    csv_path = outdir / "sigma_sweep_summary.csv"
    if args.skip_existing and csv_path.exists():
        print(f"[skip] I={I:.6g}: using existing {csv_path}")
        return csv_path

    cmd = [
        sys.executable,
        str(args.verify_script),
        "--I", str(I),
        "--a", str(args.a),
        "--b", str(args.b),
        "--eps", str(args.eps),
        "--N", str(args.N),
        "--T", str(args.T),
        "--dt", str(args.dt),
        "--r0", str(args.r0),
        "--phi0", str(args.phi0),
        "--min-radius", str(args.min_radius),
        "--spike-threshold", str(args.spike_threshold),
        "--seed", str(args.seed),
        "--outdir", str(outdir),
        "--max-paths-plot", str(args.max_paths_plot),
        "--sweep-sigma",
        *[str(s) for s in args.sigmas],
    ]
    print("\nRunning:")
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)
    if not csv_path.exists():
        raise FileNotFoundError(f"Expected {csv_path}, but it was not created.")
    return csv_path


def interpolate_sigma50(sig: np.ndarray, p: np.ndarray) -> float:
    """Linear interpolation sigma where p crosses 0.5. Returns nan if not bracketed."""
    order = np.argsort(sig)
    sig = sig[order]
    p = p[order]
    if np.all(p < 0.5) or np.all(p > 0.5):
        return float("nan")
    for k in range(len(sig) - 1):
        p0, p1 = p[k], p[k + 1]
        if (p0 - 0.5) == 0:
            return float(sig[k])
        if (p0 - 0.5) * (p1 - 0.5) <= 0 and p0 != p1:
            return float(sig[k] + (0.5 - p0) * (sig[k + 1] - sig[k]) / (p1 - p0))
    return float("nan")


def fit_sigma_curve(sig: np.ndarray, p: np.ndarray, N: int) -> Dict[str, float]:
    """Fit logistic curve; robust fallback to interpolation if SciPy/fitting fails."""
    sig = np.asarray(sig, dtype=float)
    p = np.asarray(p, dtype=float)
    ok = np.isfinite(sig) & np.isfinite(p)
    sig = sig[ok]
    p = p[ok]

    sigma50_interp = interpolate_sigma50(sig, p)

    result = {
        "sigma50_interp": sigma50_interp,
        "sigma50_fit": float("nan"),
        "width_fit": float("nan"),
        "rmse_fit": float("nan"),
        "method": "interp_only",
    }

    # Need at least one intermediate point for a meaningful sigmoid fit.
    if len(sig) < 4 or len(np.unique(p)) < 3:
        return result

    if not SCIPY_AVAILABLE:
        return result

    try:
        # Initial guesses.
        s50_guess = sigma50_interp if np.isfinite(sigma50_interp) else float(np.median(sig))
        width_guess = max((sig.max() - sig.min()) / 10.0, 1e-3)

        # Binomial standard error; clip to avoid infinite weights at p=0 or 1.
        p_clip = np.clip(p, 1.0 / (2 * N), 1.0 - 1.0 / (2 * N))
        se = np.sqrt(p_clip * (1.0 - p_clip) / max(N, 1))

        bounds = ([sig.min() - 0.05, 1e-5], [sig.max() + 0.05, 0.2])
        popt, _ = curve_fit(
            logistic,
            sig,
            p,
            p0=[s50_guess, width_guess],
            bounds=bounds,
            sigma=se,
            absolute_sigma=False,
            maxfev=20000,
        )
        pred = logistic(sig, *popt)
        rmse = float(np.sqrt(np.mean((pred - p) ** 2)))
        result.update({
            "sigma50_fit": float(popt[0]),
            "width_fit": float(popt[1]),
            "rmse_fit": rmse,
            "method": "logistic_fit",
        })
    except Exception as e:
        print(f"[warn] logistic fit failed: {e}")

    return result


def load_and_annotate(csv_path: Path, I: float, args: argparse.Namespace) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    df.insert(0, "I", I)
    df.insert(1, "a", args.a)
    df.insert(2, "b", args.b)
    df.insert(3, "eps", args.eps)
    df.insert(4, "N", args.N)
    df.insert(5, "T", args.T)
    return df


def analyse_all(all_df: pd.DataFrame, args: argparse.Namespace, outdir: Path) -> pd.DataFrame:
    rows = []
    for I, g in all_df.groupby("I", sort=True):
        sig = g["sigma"].to_numpy(float)
        p = g["excursion_fraction"].to_numpy(float)
        fit = fit_sigma_curve(sig, p, args.N)

        # Pull useful diagnostics near the transition / across escaping points.
        escaping = g[g["excursion_fraction"] > 0]
        high_escape = g[g["excursion_fraction"] >= 0.5]

        def median_col(col: str, frame: pd.DataFrame) -> float:
            if col not in frame or frame.empty:
                return float("nan")
            vals = pd.to_numeric(frame[col], errors="coerce").dropna()
            return float(vals.median()) if len(vals) else float("nan")

        rows.append({
            "I": float(I),
            "a": args.a,
            "b": args.b,
            "eps": args.eps,
            "N": args.N,
            "T": args.T,
            **fit,
            "max_excursion_fraction": float(np.nanmax(p)),
            "min_excursion_fraction": float(np.nanmin(p)),
            "median_escape_phase_R_when_escaping": median_col("escape_phase_R", escaping),
            "median_escape_phase_diff_R_when_escaping": median_col("escape_phase_diff_R", escaping),
            "median_escape_time_when_p_ge_0p5": median_col("median_spike_time", high_escape),
            "median_rotations_when_p_ge_0p5": median_col("median_rotations_before_escape", high_escape),
            "period_linear": median_col("period_linear", g),
            "omega": median_col("omega", g),
        })
    out = pd.DataFrame(rows).sort_values("I")
    out.to_csv(outdir / "sigma50_map_summary.csv", index=False)
    return out


def plot_results(all_df: pd.DataFrame, map_df: pd.DataFrame, outdir: Path) -> None:
    # Excursion probability curves by I.
    plt.figure(figsize=(8, 5))
    for I, g in all_df.groupby("I", sort=True):
        g = g.sort_values("sigma")
        plt.plot(g["sigma"], g["excursion_fraction"], marker="o", label=f"I={I:.3f}")
        if SCIPY_AVAILABLE:
            row = map_df[map_df["I"] == I].iloc[0]
            if np.isfinite(row["sigma50_fit"]) and np.isfinite(row["width_fit"]):
                xs = np.linspace(g["sigma"].min(), g["sigma"].max(), 200)
                plt.plot(xs, logistic(xs, row["sigma50_fit"], row["width_fit"]), alpha=0.35)
    plt.axhline(0.5, linestyle="--", linewidth=1, alpha=0.6)
    plt.xlabel("sigma")
    plt.ylabel("P(large excursion before T)")
    plt.title("Resonator excursion probability curves")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outdir / "excursion_probability_curves_by_I.png", dpi=180)
    plt.close()

    # sigma50 vs I.
    plt.figure(figsize=(7, 4.5))
    y = map_df["sigma50_fit"].where(np.isfinite(map_df["sigma50_fit"]), map_df["sigma50_interp"])
    plt.plot(map_df["I"], y, marker="o")
    plt.xlabel("I")
    plt.ylabel("sigma50: P(excursion before T)=0.5")
    plt.title("Finite-time resonator failure boundary")
    plt.tight_layout()
    plt.savefig(outdir / "sigma50_vs_I.png", dpi=180)
    plt.close()

    # width vs I.
    plt.figure(figsize=(7, 4.5))
    plt.plot(map_df["I"], map_df["width_fit"], marker="o")
    plt.xlabel("I")
    plt.ylabel("logistic transition width")
    plt.title("Noise-transition width across resonator window")
    plt.tight_layout()
    plt.savefig(outdir / "transition_width_vs_I.png", dpi=180)
    plt.close()

    # escape phase concentration by sigma, per I.
    plt.figure(figsize=(8, 5))
    for I, g in all_df.groupby("I", sort=True):
        g = g.sort_values("sigma")
        if "escape_phase_R" in g:
            plt.plot(g["sigma"], g["escape_phase_R"], marker="o", label=f"I={I:.3f}")
    plt.xlabel("sigma")
    plt.ylabel("R of absolute escape phase")
    plt.ylim(-0.05, 1.05)
    plt.title("Escape phase concentration across noise")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outdir / "escape_phase_R_by_I.png", dpi=180)
    plt.close()

    # median rotations by sigma, per I.
    plt.figure(figsize=(8, 5))
    for I, g in all_df.groupby("I", sort=True):
        g = g.sort_values("sigma")
        if "median_rotations_before_escape" in g:
            plt.plot(g["sigma"], g["median_rotations_before_escape"], marker="o", label=f"I={I:.3f}")
    plt.xlabel("sigma")
    plt.ylabel("median rotations before escape")
    plt.title("How many resonator rotations before commitment?")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outdir / "median_rotations_by_I.png", dpi=180)
    plt.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Map the finite-time resonator excursion threshold sigma50(I)."
    )
    parser.add_argument("--verify-script", type=Path, default=Path("verify_resonator_phase_difference.py"))
    parser.add_argument("--I-values", type=float, nargs="+", default=[0.295, 0.300, 0.305, 0.310, 0.315, 0.320])
    parser.add_argument("--sigmas", type=float, nargs="+", default=[0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.050, 0.060])
    parser.add_argument("--a", type=float, default=0.7)
    parser.add_argument("--b", type=float, default=0.8)
    parser.add_argument("--eps", type=float, default=0.08)
    parser.add_argument("--N", type=int, default=500)
    parser.add_argument("--T", type=float, default=250.0)
    parser.add_argument("--dt", type=float, default=0.005)
    parser.add_argument("--r0", type=float, default=0.06)
    parser.add_argument("--phi0", type=float, default=0.0)
    parser.add_argument("--min-radius", type=float, default=0.004)
    parser.add_argument("--spike-threshold", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--outdir", type=Path, default=Path("results/resonator_sigma50_map"))
    parser.add_argument("--max-paths-plot", type=int, default=20)
    parser.add_argument("--skip-existing", action="store_true", help="Reuse existing per-I sigma_sweep_summary.csv files.")
    parser.add_argument("--analysis-only", action="store_true", help="Do not run simulations; only analyse existing per-I summaries.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = ensure_dir(args.outdir)

    if not args.analysis_only and not args.verify_script.exists():
        raise FileNotFoundError(
            f"Could not find verifier script: {args.verify_script}. "
            "Pass --verify-script path/to/verify_resonator_phase_difference.py"
        )

    frames = []
    for I in args.I_values:
        I_dir = ensure_dir(root / safe_I_label(I))
        csv_path = I_dir / "sigma_sweep_summary.csv"
        if args.analysis_only:
            if not csv_path.exists():
                raise FileNotFoundError(f"Missing {csv_path}; cannot use --analysis-only.")
        else:
            csv_path = run_one_I(args, I=I, outdir=I_dir)
        frames.append(load_and_annotate(csv_path, I=I, args=args))

    all_df = pd.concat(frames, ignore_index=True)
    all_df.to_csv(root / "all_sigma_sweeps.csv", index=False)

    map_df = analyse_all(all_df, args=args, outdir=root)
    plot_results(all_df, map_df, root)

    print("\nDone.")
    print(f"Wrote: {root / 'all_sigma_sweeps.csv'}")
    print(f"Wrote: {root / 'sigma50_map_summary.csv'}")
    print("Main plots:")
    for name in [
        "excursion_probability_curves_by_I.png",
        "sigma50_vs_I.png",
        "transition_width_vs_I.png",
        "escape_phase_R_by_I.png",
        "median_rotations_by_I.png",
    ]:
        print(f"  {root / name}")

    print("\nSigma50 summary:")
    cols = ["I", "sigma50_fit", "sigma50_interp", "width_fit", "max_excursion_fraction", "median_escape_phase_R_when_escaping"]
    print(map_df[cols].to_string(index=False))


if __name__ == "__main__":
    main()
