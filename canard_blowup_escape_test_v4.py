#!/usr/bin/env python3
"""
canard_blowup_escape_test_v4.py

V4 fixes the main numerical issue from V3:

    small lambda values were not using the intended scale-invariant window
        W in [window_lo * W_star, window_hi * W_star],
        W_star = lambda^(2/3),

    because a hard W_cut floor, e.g. W_cut=0.02, clipped the lower window.

This script lets you either:
  1. use the true scale-invariant window with no hard lower clipping;
  2. skip lambdas whose requested window would be clipped;
  3. explicitly set a tiny numerical W_min only for avoiding sqrt(W) at W=0.

It also fits eta50 ~ C lambda^beta over:
  - all valid lambdas;
  - lambdas whose window was not clipped;
  - optional lambda_min / lambda_max restricted ranges.

Recommended first run:

python3 canard_blowup_escape_test_v4.py \
  --overwrite \
  --event-mode window \
  --window-lo 0.5 \
  --window-hi 2.0 \
  --target-mode frac \
  --target-frac 0.75 \
  --lambdas 0.004 0.006 0.008 0.012 0.016 0.024 0.032 0.048 0.064 \
  --etas 0 0.02 0.035 0.05 0.07 0.10 0.14 0.20 0.28 0.40 0.60 0.85 1.10 \
  --n 2000 \
  --outdir results/canard_blowup_escape_v4_q075

If the smallest lambda values are still unstable, fit the clean range:

python3 canard_blowup_escape_test_v4.py \
  --overwrite \
  --event-mode window \
  --window-lo 0.5 \
  --window-hi 2.0 \
  --target-mode frac \
  --target-frac 0.75 \
  --fit-lambda-min 0.008 \
  --lambdas 0.004 0.006 0.008 0.012 0.016 0.024 0.032 0.048 0.064 \
  --etas 0 0.02 0.035 0.05 0.07 0.10 0.14 0.20 0.28 0.40 0.60 0.85 1.10 \
  --n 2000 \
  --outdir results/canard_blowup_escape_v4_q075_fit008
"""

from __future__ import annotations

import argparse
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


@dataclass
class SimConfig:
    lam: float
    eta: float
    n: int
    dt: float
    seed: int
    W_start_factor: float
    W_end: float
    event_mode: str
    window_lo: float
    window_hi: float
    target_mode: str
    target_frac: float
    peel_dist: float
    W_min: float
    hard_W_cut: float | None
    skip_clipped_windows: bool


def sigmoid_interp_x_at_y(x: np.ndarray, y: np.ndarray, target: float = 0.5) -> float:
    """
    Interpolate the first x where y crosses target.
    Assumes x increasing. y need not be perfectly monotone, but should mostly be increasing.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    order = np.argsort(x)
    x = x[order]
    y = y[order]

    if np.any(~np.isfinite(x)) or np.any(~np.isfinite(y)):
        mask = np.isfinite(x) & np.isfinite(y)
        x, y = x[mask], y[mask]

    if len(x) < 2:
        return np.nan

    if np.nanmax(y) < target:
        return np.nan
    if np.nanmin(y) > target:
        return x[0]

    for i in range(len(x) - 1):
        y0, y1 = y[i], y[i + 1]
        if (y0 - target) == 0:
            return x[i]
        if (y0 - target) * (y1 - target) <= 0 and y0 != y1:
            a = (target - y0) / (y1 - y0)
            return x[i] + a * (x[i + 1] - x[i])

    return np.nan


def fit_power_law(lam: np.ndarray, eta50: np.ndarray) -> dict:
    """
    Fit eta50 = C * lambda^beta by OLS in log-log space.
    """
    lam = np.asarray(lam, dtype=float)
    eta50 = np.asarray(eta50, dtype=float)

    mask = np.isfinite(lam) & np.isfinite(eta50) & (lam > 0) & (eta50 > 0)
    lam = lam[mask]
    eta50 = eta50[mask]

    if len(lam) < 2:
        return {
            "C": np.nan,
            "beta": np.nan,
            "r2": np.nan,
            "n_fit": len(lam),
        }

    x = np.log(lam)
    y = np.log(eta50)

    beta, logC = np.polyfit(x, y, 1)
    yhat = logC + beta * x

    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = np.nan if ss_tot == 0 else 1.0 - ss_res / ss_tot

    return {
        "C": float(np.exp(logC)),
        "beta": float(beta),
        "r2": float(r2),
        "n_fit": int(len(lam)),
    }


def simulate_one_probability(cfg: SimConfig) -> dict:
    """
    Simulate dV = (V^2 - W)dT + eta dB_T, dW = -lambda dT.

    Event:
      all_positive:
        count hit anywhere W > W_min / hard W cut.
      window:
        count hit only while W in [window_lo W*, window_hi W*],
        unless hard_W_cut clips the lower window.

    target-mode frac:
        hit if V >= target_frac * sqrt(W).

    target-mode dist:
        hit if V >= sqrt(W) - peel_dist.
    """
    lam = cfg.lam
    eta = cfg.eta
    rng = np.random.default_rng(cfg.seed)

    W_star = lam ** (2.0 / 3.0)
    W_start = cfg.W_start_factor * W_star
    W_end = cfg.W_end

    requested_W_hi = cfg.window_hi * W_star
    requested_W_lo = cfg.window_lo * W_star

    # Tiny numerical floor only; this should not define the asymptotic event.
    numerical_floor = max(cfg.W_min, 0.0)

    if cfg.hard_W_cut is None:
        effective_W_lo = max(requested_W_lo, numerical_floor)
        clipped = False
    else:
        effective_W_lo = max(requested_W_lo, cfg.hard_W_cut, numerical_floor)
        clipped = effective_W_lo > max(requested_W_lo, numerical_floor) + 1e-15

    effective_W_hi = requested_W_hi

    if cfg.event_mode == "window":
        if cfg.skip_clipped_windows and clipped:
            return {
                "lambda": lam,
                "eta": eta,
                "hit_prob": np.nan,
                "n": cfg.n,
                "W_star": W_star,
                "requested_W_lo": requested_W_lo,
                "requested_W_hi": requested_W_hi,
                "effective_W_lo": effective_W_lo,
                "effective_W_hi": effective_W_hi,
                "window_clipped": clipped,
                "skipped": True,
            }
        if effective_W_hi <= effective_W_lo:
            return {
                "lambda": lam,
                "eta": eta,
                "hit_prob": np.nan,
                "n": cfg.n,
                "W_star": W_star,
                "requested_W_lo": requested_W_lo,
                "requested_W_hi": requested_W_hi,
                "effective_W_lo": effective_W_lo,
                "effective_W_hi": effective_W_hi,
                "window_clipped": clipped,
                "skipped": True,
            }

    T_total = (W_start - W_end) / lam
    n_steps = int(np.ceil(T_total / cfg.dt))
    dt = T_total / n_steps

    # Start near attracting branch: V = -sqrt(W).
    V = -np.sqrt(np.full(cfg.n, W_start))
    hit = np.zeros(cfg.n, dtype=bool)

    sqrt_dt_eta = eta * math.sqrt(dt)

    W = W_start
    for k in range(n_steps):
        # Euler-Maruyama update.
        dB = rng.normal(0.0, 1.0, size=cfg.n)
        V += (V * V - W) * dt + sqrt_dt_eta * dB
        W -= lam * dt

        if W <= numerical_floor:
            # Avoid sqrt of nonpositive W.
            continue

        if cfg.event_mode == "window":
            active_W = (effective_W_lo <= W <= effective_W_hi)
        elif cfg.event_mode == "all_positive":
            if cfg.hard_W_cut is None:
                active_W = W > numerical_floor
            else:
                active_W = W > max(cfg.hard_W_cut, numerical_floor)
        else:
            raise ValueError(f"Unknown event_mode: {cfg.event_mode}")

        if not active_W:
            continue

        branch = math.sqrt(W)
        if cfg.target_mode == "frac":
            threshold = cfg.target_frac * branch
        elif cfg.target_mode == "dist":
            threshold = branch - cfg.peel_dist
        else:
            raise ValueError(f"Unknown target_mode: {cfg.target_mode}")

        hit_now = V >= threshold
        hit |= hit_now

    return {
        "lambda": lam,
        "eta": eta,
        "hit_prob": float(np.mean(hit)),
        "n": cfg.n,
        "W_star": W_star,
        "requested_W_lo": requested_W_lo,
        "requested_W_hi": requested_W_hi,
        "effective_W_lo": effective_W_lo,
        "effective_W_hi": effective_W_hi,
        "window_clipped": bool(clipped),
        "skipped": False,
    }


def run_grid(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    counter = 0

    for lam in args.lambdas:
        for eta in args.etas:
            counter += 1
            cfg = SimConfig(
                lam=float(lam),
                eta=float(eta),
                n=int(args.n),
                dt=float(args.dt),
                seed=int(args.seed + 100000 * counter),
                W_start_factor=float(args.W_start_factor),
                W_end=float(args.W_end),
                event_mode=args.event_mode,
                window_lo=float(args.window_lo),
                window_hi=float(args.window_hi),
                target_mode=args.target_mode,
                target_frac=float(args.target_frac),
                peel_dist=float(args.peel_dist),
                W_min=float(args.W_min),
                hard_W_cut=args.hard_W_cut,
                skip_clipped_windows=bool(args.skip_clipped_windows),
            )
            row = simulate_one_probability(cfg)
            rows.append(row)

            print(
                f"lambda={lam:g}, eta={eta:g}, "
                f"P={row['hit_prob'] if np.isfinite(row['hit_prob']) else np.nan:.3f}, "
                f"W*={row['W_star']:.5g}, "
                f"Wwin=[{row['effective_W_lo']:.5g},{row['effective_W_hi']:.5g}], "
                f"clipped={row['window_clipped']}, skipped={row['skipped']}"
            )

    results = pd.DataFrame(rows)

    eta50_rows = []
    for lam, group in results.groupby("lambda"):
        clean = group.dropna(subset=["hit_prob"]).sort_values("eta")
        eta50 = sigmoid_interp_x_at_y(clean["eta"].to_numpy(), clean["hit_prob"].to_numpy(), 0.5)

        first = group.iloc[0]
        eta50_rows.append(
            {
                "lambda": float(lam),
                "eta50_interp": float(eta50) if np.isfinite(eta50) else np.nan,
                "W_star": float(first["W_star"]),
                "requested_W_lo": float(first["requested_W_lo"]),
                "requested_W_hi": float(first["requested_W_hi"]),
                "effective_W_lo": float(first["effective_W_lo"]),
                "effective_W_hi": float(first["effective_W_hi"]),
                "window_clipped": bool(first["window_clipped"]),
                "used_points": int(clean["hit_prob"].notna().sum()),
            }
        )

    sigma50 = pd.DataFrame(eta50_rows).sort_values("lambda")

    # Add reference predictions.
    sigma50["lambda_pow_half"] = sigma50["lambda"] ** 0.5
    sigma50["lambda_pow_2over3"] = sigma50["lambda"] ** (2.0 / 3.0)

    # Fit all.
    all_fit = fit_power_law(sigma50["lambda"], sigma50["eta50_interp"])

    # Fit unclipped.
    unclipped = sigma50[~sigma50["window_clipped"]]
    unclipped_fit = fit_power_law(unclipped["lambda"], unclipped["eta50_interp"])

    # Fit user range.
    ranged = sigma50.copy()
    if args.fit_lambda_min is not None:
        ranged = ranged[ranged["lambda"] >= args.fit_lambda_min]
    if args.fit_lambda_max is not None:
        ranged = ranged[ranged["lambda"] <= args.fit_lambda_max]
    ranged_fit = fit_power_law(ranged["lambda"], ranged["eta50_interp"])

    fit_summary = {
        "all": all_fit,
        "unclipped": unclipped_fit,
        "range": ranged_fit,
        "fit_lambda_min": args.fit_lambda_min,
        "fit_lambda_max": args.fit_lambda_max,
    }

    return results, sigma50, fit_summary


def plot_outputs(results: pd.DataFrame, sigma50: pd.DataFrame, fit_summary: dict, args: argparse.Namespace) -> None:
    outdir = Path(args.outdir)

    # Escape curves.
    plt.figure(figsize=(11, 7))
    for lam, group in results.groupby("lambda"):
        group = group.dropna(subset=["hit_prob"]).sort_values("eta")
        plt.plot(group["eta"], group["hit_prob"], marker="o", label=f"lambda={lam:g}")

    plt.axhline(0.5, linestyle=":", color="black")
    plt.xlabel(r"$\eta = \sigma/\sqrt{\varepsilon}$")
    plt.ylabel("P(hit repelling side before fold)")
    plt.title("Blow-up normal form escape curves")
    plt.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(outdir / "blowup_escape_curves_v4.png", dpi=180)
    plt.close()

    # eta50 vs lambda with fitted half and 2/3 references.
    clean = sigma50.dropna(subset=["eta50_interp"]).copy()
    plt.figure(figsize=(10, 7))
    plt.loglog(clean["lambda"], clean["eta50_interp"], marker="o", label="measured eta50")

    lam_grid = np.linspace(clean["lambda"].min(), clean["lambda"].max(), 200)

    # Scale references through median to compare slopes visually.
    med_eta = np.nanmedian(clean["eta50_interp"])
    med_lam = np.nanmedian(clean["lambda"])
    ref_half = med_eta * (lam_grid / med_lam) ** 0.5
    ref_23 = med_eta * (lam_grid / med_lam) ** (2.0 / 3.0)

    plt.loglog(lam_grid, ref_half, "--", label=r"$\lambda^{1/2}$ reference")
    plt.loglog(lam_grid, ref_23, "--", label=r"$\lambda^{2/3}$ reference")

    beta_all = fit_summary["all"]["beta"]
    beta_range = fit_summary["range"]["beta"]

    title = rf"Blow-up normal form: $\eta_{{50}}$ vs $\lambda$"
    if np.isfinite(beta_all):
        title += rf" | all $\beta={beta_all:.3f}$"
    if np.isfinite(beta_range):
        title += rf" | range $\beta={beta_range:.3f}$"

    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$\eta_{50}$")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "blowup_eta50_vs_lambda_v4.png", dpi=180)
    plt.close()

    # Collapse ratio for lambda^1/2.
    C_half = np.nanmedian(clean["eta50_interp"] / np.sqrt(clean["lambda"]))
    clean["ratio_half"] = clean["eta50_interp"] / (C_half * np.sqrt(clean["lambda"]))

    plt.figure(figsize=(10, 6))
    plt.semilogx(clean["lambda"], clean["ratio_half"], marker="o")
    plt.axhline(1.0, linestyle="--", color="black", label=rf"median-normalised $C={C_half:.3g}$")
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"measured / $C\lambda^{1/2}$")
    plt.title(r"Collapse ratio using $\lambda^{1/2}$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "blowup_collapse_ratio_half_v4.png", dpi=180)
    plt.close()

    # Collapse ratio for lambda^2/3.
    C_23 = np.nanmedian(clean["eta50_interp"] / (clean["lambda"] ** (2.0 / 3.0)))
    clean["ratio_23"] = clean["eta50_interp"] / (C_23 * (clean["lambda"] ** (2.0 / 3.0)))

    plt.figure(figsize=(10, 6))
    plt.semilogx(clean["lambda"], clean["ratio_23"], marker="o")
    plt.axhline(1.0, linestyle="--", color="black", label=rf"median-normalised $C={C_23:.3g}$")
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"measured / $C\lambda^{2/3}$")
    plt.title(r"Collapse ratio using $\lambda^{2/3}$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "blowup_collapse_ratio_2over3_v4.png", dpi=180)
    plt.close()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()

    p.add_argument("--outdir", type=str, default="results/canard_blowup_escape_v4")
    p.add_argument("--overwrite", action="store_true")

    p.add_argument("--lambdas", type=float, nargs="+", required=True)
    p.add_argument("--etas", type=float, nargs="+", required=True)
    p.add_argument("--n", type=int, default=2000)
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--dt", type=float, default=0.002)

    p.add_argument("--event-mode", choices=["window", "all_positive"], default="window")
    p.add_argument("--window-lo", type=float, default=0.5)
    p.add_argument("--window-hi", type=float, default=2.0)

    p.add_argument("--target-mode", choices=["frac", "dist"], default="frac")
    p.add_argument("--target-frac", type=float, default=0.75)
    p.add_argument("--peel-dist", type=float, default=0.05)

    # This starts high enough above W* to settle onto the attracting branch.
    p.add_argument("--W-start-factor", type=float, default=8.0)

    # Integrate until W_end. Negative values allow paths to pass through the fold.
    p.add_argument("--W-end", type=float, default=-0.05)

    # Tiny numerical floor only. Should usually be 0 or tiny, not 0.02.
    p.add_argument("--W-min", type=float, default=1e-8)

    # Optional old-style hard cutoff. Leave omitted for clean scale-invariant windows.
    p.add_argument("--hard-W-cut", type=float, default=None)

    # If using --hard-W-cut, skip lambdas where it clips the requested window.
    p.add_argument("--skip-clipped-windows", action="store_true")

    # Fit range controls.
    p.add_argument("--fit-lambda-min", type=float, default=None)
    p.add_argument("--fit-lambda-max", type=float, default=None)

    return p.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)

    if outdir.exists() and any(outdir.iterdir()) and not args.overwrite:
        raise SystemExit(f"Output directory {outdir} exists and is nonempty. Use --overwrite.")

    outdir.mkdir(parents=True, exist_ok=True)

    print("\n=== V4 blow-up normal-form canard escape test ===")
    print(f"event_mode={args.event_mode}")
    print(f"target_mode={args.target_mode}")
    print(f"target_frac={args.target_frac}")
    print(f"window=[{args.window_lo} W*, {args.window_hi} W*]")
    print(f"W_min={args.W_min}")
    print(f"hard_W_cut={args.hard_W_cut}")
    print(f"skip_clipped_windows={args.skip_clipped_windows}")
    print(f"n={args.n}, dt={args.dt}")
    print("")

    results, sigma50, fit_summary = run_grid(args)

    results_path = outdir / "blowup_results_v4.csv"
    sigma50_path = outdir / "blowup_sigma50_v4.csv"
    summary_path = outdir / "blowup_fit_summary_v4.txt"

    results.to_csv(results_path, index=False)
    sigma50.to_csv(sigma50_path, index=False)

    with open(summary_path, "w") as f:
        f.write("V4 blow-up normal-form fit summary\n")
        f.write("===================================\n\n")
        f.write(f"event_mode={args.event_mode}\n")
        f.write(f"target_mode={args.target_mode}\n")
        f.write(f"target_frac={args.target_frac}\n")
        f.write(f"window=[{args.window_lo} W*, {args.window_hi} W*]\n")
        f.write(f"W_min={args.W_min}\n")
        f.write(f"hard_W_cut={args.hard_W_cut}\n")
        f.write(f"fit_lambda_min={args.fit_lambda_min}\n")
        f.write(f"fit_lambda_max={args.fit_lambda_max}\n\n")

        for name, fit in fit_summary.items():
            if isinstance(fit, dict):
                f.write(
                    f"{name}: eta50 = C lambda^beta, "
                    f"C={fit['C']:.6g}, beta={fit['beta']:.6g}, "
                    f"R2={fit['r2']:.6g}, n_fit={fit['n_fit']}\n"
                )

    plot_outputs(results, sigma50, fit_summary, args)

    print("\n=== Fit summary ===")
    for name, fit in fit_summary.items():
        if isinstance(fit, dict):
            print(
                f"{name:10s}: C={fit['C']:.6g}, "
                f"beta={fit['beta']:.6g}, "
                f"R2={fit['r2']:.6g}, n_fit={fit['n_fit']}"
            )

    print("\nSaved:")
    print(f"  {results_path}")
    print(f"  {sigma50_path}")
    print(f"  {summary_path}")
    print(f"  {outdir / 'blowup_escape_curves_v4.png'}")
    print(f"  {outdir / 'blowup_eta50_vs_lambda_v4.png'}")
    print(f"  {outdir / 'blowup_collapse_ratio_half_v4.png'}")
    print(f"  {outdir / 'blowup_collapse_ratio_2over3_v4.png'}")


if __name__ == "__main__":
    main()
