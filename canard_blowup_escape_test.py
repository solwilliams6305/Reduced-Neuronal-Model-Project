#!/usr/bin/env python3
"""
canard_blowup_escape_test_v2.py

Blow-up normal-form diagnostic for canard escape.

Normal form:
    dV = (V^2 - W) dT + eta dB_T
    dW = -lambda dT

This script deliberately does NOT use a deterministic peel-off baseline, because
in the singular blow-up normal form the deterministic path starting on the
attracting branch does not peel to the repelling branch at positive W.

Instead it measures:
    P(hit)
where "hit" means the noisy path reaches a chosen target near the repelling side
while W is still positive.

Default event:
    V >= sqrt(W) - peel_dist
while W >= W_cut.

For a more local test around the predicted optimal jump scale
    W_star = lambda^(2/3),
use:
    --event-mode window
which only counts hits while W is in
    [window-lo * W_star, window-hi * W_star].

The diagnostic fit is:
    eta_50 ~ C * lambda^beta
with the candidate prediction beta = 2/3.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path

import numpy as np


def simulate_hit_probability(
    lam: float,
    eta: float,
    *,
    n: int,
    dt: float,
    W0: float,
    T: float,
    W_cut: float,
    peel_dist: float,
    target_frac: float,
    target_mode: str,
    event_mode: str,
    window_lo: float,
    window_hi: float,
    seed: int,
) -> dict:
    """
    Simulate the blow-up SDE and return hit probability and hit stats.

    Event:
        V >= sqrt(W) - peel_dist

    Count only while:
        all_positive mode: W >= W_cut
        window mode:       window_lo*W_star <= W <= window_hi*W_star

    Notes:
    - Once W falls below the relevant lower cutoff, remaining unhit paths are
      considered "not hit"; we stop early.
    - If V blows up to +infinity after a counted hit, that is fine. If it blows
      up without a counted hit, it is not counted.
    """
    rng = np.random.default_rng(seed)

    V = np.full(n, -math.sqrt(W0), dtype=float)
    W = np.full(n, W0, dtype=float)

    hit = np.zeros(n, dtype=bool)
    t_hit = np.full(n, np.nan)
    W_hit = np.full(n, np.nan)
    V_hit = np.full(n, np.nan)

    sqrt_dt = math.sqrt(dt)
    W_star = lam ** (2.0 / 3.0)

    if event_mode == "window":
        W_low = max(W_cut, window_lo * W_star)
        W_high = window_hi * W_star
    else:
        W_low = W_cut
        W_high = float("inf")

    n_steps = int(math.ceil(T / dt))
    t = 0.0

    for _ in range(n_steps):
        # Only evaluate events while W is in the desired positive region.
        active_window = (W >= W_low) & (W <= W_high)

        if np.any(~hit & active_window):
            sqrtW = np.sqrt(np.maximum(W, 0.0))
            if target_mode == "frac":
                target = target_frac * sqrtW
            else:
                target = sqrtW - peel_dist
            new_hit = (~hit) & active_window & (V >= target)
            if np.any(new_hit):
                hit[new_hit] = True
                t_hit[new_hit] = t
                W_hit[new_hit] = W[new_hit]
                V_hit[new_hit] = V[new_hit]

        # Stop once all unhit paths have passed below the counting region.
        if np.all(hit | (W < W_low)):
            break

        noise = rng.standard_normal(n) if eta > 0 else 0.0

        # Euler-Maruyama update.
        # Clip huge positive V to avoid overflow poisoning the whole vector.
        drift = V * V - W
        V = V + drift * dt + eta * sqrt_dt * noise
        W = W - lam * dt
        t += dt

        V = np.where(V > 1e6, 1e6, V)
        V = np.where(V < -1e6, -1e6, V)

    valid = hit
    return {
        "hit_prob": float(np.mean(hit)),
        "mean_t_hit": float(np.nanmean(t_hit[valid])) if np.any(valid) else np.nan,
        "mean_W_hit": float(np.nanmean(W_hit[valid])) if np.any(valid) else np.nan,
        "mean_V_hit": float(np.nanmean(V_hit[valid])) if np.any(valid) else np.nan,
        "W_star": float(W_star),
        "W_low": float(W_low),
        "W_high": float(W_high),
    }


def interp_x_at_y(x: np.ndarray, y: np.ndarray, target: float = 0.5) -> float:
    """Linear interpolation for first crossing y >= target."""
    order = np.argsort(x)
    x = x[order]
    y = y[order]

    above = np.where(y >= target)[0]
    if len(above) == 0:
        return float("nan")

    k = int(above[0])
    if k == 0:
        return float(x[0])

    x0, x1 = x[k - 1], x[k]
    y0, y1 = y[k - 1], y[k]

    if not np.isfinite(y0) or not np.isfinite(y1) or y1 == y0:
        return float(x1)

    return float(x0 + (target - y0) * (x1 - x0) / (y1 - y0))


def fit_power_law(xs: np.ndarray, ys: np.ndarray) -> tuple[float, float, float]:
    """Fit y = C x^p. Returns p, C, R2."""
    mask = np.isfinite(xs) & np.isfinite(ys) & (xs > 0) & (ys > 0)
    if int(mask.sum()) < 3:
        return np.nan, np.nan, np.nan

    lx = np.log(xs[mask])
    ly = np.log(ys[mask])
    p, logC = np.polyfit(lx, ly, 1)
    pred = p * lx + logC
    ss_res = float(np.sum((ly - pred) ** 2))
    ss_tot = float(np.sum((ly - np.mean(ly)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return float(p), float(np.exp(logC)), float(r2)


def write_csv(path: Path, rows: list[dict], fieldnames: list[str], overwrite: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    mode = "w" if overwrite or not path.exists() else "a"
    with path.open(mode, newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        if mode == "w":
            w.writeheader()
        for r in rows:
            w.writerow(r)


def make_plots(outdir: Path, result_rows: list[dict], sigma_rows: list[dict]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Escape curves
    plt.figure(figsize=(8, 5.5))
    lambdas = sorted({float(r["lambda"]) for r in result_rows})
    for lam in lambdas:
        pts = sorted([r for r in result_rows if float(r["lambda"]) == lam], key=lambda r: float(r["eta"]))
        eta = np.array([float(p["eta"]) for p in pts])
        hp = np.array([float(p["hit_prob"]) for p in pts])
        plt.plot(eta, hp, "o-", label=f"lambda={lam:g}")
    plt.axhline(0.5, color="k", ls=":", lw=1.5)
    plt.xlabel(r"$\eta = \sigma/\sqrt{\varepsilon}$")
    plt.ylabel("P(hit repelling side before fold)")
    plt.title("Blow-up normal form escape curves")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "blowup_escape_curves_v3.png", dpi=160)
    plt.close()

    # eta50 vs lambda
    if sigma_rows:
        lam = np.array([float(r["lambda"]) for r in sigma_rows])
        eta50 = np.array([float(r["eta50"]) for r in sigma_rows])
        pred = np.array([float(r["pred_eta"]) for r in sigma_rows])

        plt.figure(figsize=(7, 5.3))
        plt.loglog(lam, eta50, "o-", label="measured eta50")
        plt.loglog(lam, pred, "--", label=r"$\lambda^{2/3}$ reference")
        plt.xlabel(r"$\lambda$")
        plt.ylabel(r"$\eta_{50}$")
        plt.title(r"Blow-up normal form: $\eta_{50}$ vs $\lambda$")
        plt.legend()
        plt.tight_layout()
        plt.savefig(outdir / "blowup_eta50_vs_lambda_v3.png", dpi=160)
        plt.close()

        ratio = eta50 / pred
        plt.figure(figsize=(7, 5.3))
        plt.semilogx(lam, ratio, "o")
        plt.axhline(np.nanmedian(ratio), color="k", ls="--", label=f"median={np.nanmedian(ratio):.3g}")
        plt.xlabel(r"$\lambda$")
        plt.ylabel("measured / predicted")
        plt.title("Blow-up collapse ratio")
        plt.legend()
        plt.tight_layout()
        plt.savefig(outdir / "blowup_collapse_ratio_v3.png", dpi=160)
        plt.close()


def main() -> None:
    ap = argparse.ArgumentParser()

    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--overwrite", action="store_true")

    ap.add_argument("--lambdas", type=float, nargs="+", default=None)
    ap.add_argument("--etas", type=float, nargs="+", default=None)

    ap.add_argument("--n", type=int, default=None)
    ap.add_argument("--dt", type=float, default=2e-3)
    ap.add_argument("--W0", type=float, default=1.0)
    ap.add_argument("--T", type=float, default=None)
    ap.add_argument("--T-after", type=float, default=20.0)

    ap.add_argument("--W-cut", type=float, default=0.02)
    ap.add_argument("--peel-dist", type=float, default=0.05)
    ap.add_argument("--target-mode", choices=["frac", "dist"], default="frac")
    ap.add_argument("--target-frac", type=float, default=0.75)
    ap.add_argument("--event-mode", choices=["all_positive", "window"], default="all_positive")
    ap.add_argument("--window-lo", type=float, default=0.5)
    ap.add_argument("--window-hi", type=float, default=2.0)

    ap.add_argument("--seed", type=int, default=123)
    ap.add_argument("--outdir", type=str, default="results/canard_blowup_escape")

    args = ap.parse_args()

    if args.quick:
        lambdas = args.lambdas if args.lambdas is not None else [0.008, 0.016, 0.032]
        etas = args.etas if args.etas is not None else [0.0, 0.02, 0.035, 0.05, 0.07, 0.10, 0.14, 0.20, 0.28, 0.40, 0.60]
        n = args.n or 800
    else:
        lambdas = args.lambdas if args.lambdas is not None else [0.004, 0.006, 0.008, 0.012, 0.016, 0.024, 0.032, 0.048, 0.064]
        etas = args.etas if args.etas is not None else [0.0, 0.015, 0.025, 0.035, 0.05, 0.07, 0.10, 0.14, 0.20, 0.28, 0.40, 0.60, 0.85]
        n = args.n or 2000

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("\n=== V3 blow-up normal-form canard escape test ===")
    print("Event: target-mode frac means V >= target_frac*sqrt(W); dist means V >= sqrt(W)-peel_dist")
    print(f"event_mode={args.event_mode}, W_cut={args.W_cut}, peel_dist={args.peel_dist}")
    if args.event_mode == "window":
        print(f"window=[{args.window_lo}, {args.window_hi}] * lambda^(2/3)")
    print(f"lambdas: {lambdas}")
    print(f"etas   : {etas}")
    print(f"n={n}, dt={args.dt}, W0={args.W0}\n")

    result_rows: list[dict] = []

    for i_lam, lam in enumerate(lambdas):
        pred_eta = lam ** (2.0 / 3.0)
        T_needed = args.W0 / lam + args.T_after
        T_run = args.T if args.T is not None else T_needed
        T_run = max(T_run, T_needed)

        print(f"--- lambda={lam:g}; pred eta~lambda^(2/3)={pred_eta:.6g}; T_run={T_run:.2f} ---")

        for i_eta, eta in enumerate(etas):
            res = simulate_hit_probability(
                lam,
                eta,
                n=n,
                dt=args.dt,
                W0=args.W0,
                T=T_run,
                W_cut=args.W_cut,
                peel_dist=args.peel_dist,
                target_frac=args.target_frac,
                target_mode=args.target_mode,
                event_mode=args.event_mode,
                window_lo=args.window_lo,
                window_hi=args.window_hi,
                seed=args.seed + 1000 * i_lam + i_eta,
            )

            row = {
                "lambda": lam,
                "eta": eta,
                "pred_eta": pred_eta,
                "hit_prob": res["hit_prob"],
                "mean_t_hit": res["mean_t_hit"],
                "mean_W_hit": res["mean_W_hit"],
                "mean_V_hit": res["mean_V_hit"],
                "W_star": res["W_star"],
                "W_low": res["W_low"],
                "W_high": res["W_high"],
                "event_mode": args.event_mode,
                "peel_dist": args.peel_dist,
                "W_cut": args.W_cut,
                "n": n,
                "dt": args.dt,
                "T_run": T_run,
            }
            result_rows.append(row)

            print(
                f" eta={eta:.5f} | hit={res['hit_prob']:.3f} "
                f"meanW={res['mean_W_hit']:.5g} meant={res['mean_t_hit']:.3g}"
            )
        print()

    result_fields = [
        "lambda", "eta", "pred_eta", "hit_prob",
        "mean_t_hit", "mean_W_hit", "mean_V_hit",
        "W_star", "W_low", "W_high",
        "event_mode", "peel_dist", "W_cut",
        "n", "dt", "T_run",
    ]
    write_csv(outdir / "blowup_results_v3.csv", result_rows, result_fields, args.overwrite)

    # Extract eta50
    sigma_rows: list[dict] = []
    for lam in sorted(set(float(r["lambda"]) for r in result_rows)):
        pts = sorted([r for r in result_rows if float(r["lambda"]) == lam], key=lambda r: float(r["eta"]))
        eta_arr = np.array([float(p["eta"]) for p in pts])
        hp_arr = np.array([float(p["hit_prob"]) for p in pts])
        eta50 = interp_x_at_y(eta_arr, hp_arr, 0.5)
        pred_eta = lam ** (2.0 / 3.0)
        sigma_rows.append({
            "lambda": lam,
            "eta50": eta50,
            "pred_eta": pred_eta,
            "ratio": eta50 / pred_eta if np.isfinite(eta50) else np.nan,
        })
        print(f"lambda={lam:g}: eta50={eta50:.6g}, pred={pred_eta:.6g}, ratio={eta50/pred_eta if np.isfinite(eta50) else np.nan:.3g}")

    sigma_fields = ["lambda", "eta50", "pred_eta", "ratio"]
    write_csv(outdir / "blowup_sigma50_v3.csv", sigma_rows, sigma_fields, args.overwrite)

    lam_fit = np.array([float(r["lambda"]) for r in sigma_rows])
    eta_fit = np.array([float(r["eta50"]) for r in sigma_rows])
    beta, C, r2 = fit_power_law(lam_fit, eta_fit)

    print("\n=== Fit: eta50 ~ C lambda^beta ===")
    print(f"beta={beta:.3f}, C={C:.4g}, R2={r2:.3f}   target beta=2/3")
    if np.any(np.isfinite(eta_fit / (lam_fit ** (2.0 / 3.0)))):
        ratios = eta_fit / (lam_fit ** (2.0 / 3.0))
        print(f"median eta50/pred={np.nanmedian(ratios):.3g}, range=[{np.nanmin(ratios):.3g}, {np.nanmax(ratios):.3g}]")

    make_plots(outdir, result_rows, sigma_rows)

    print(f"\n-> wrote {outdir/'blowup_results_v3.csv'}")
    print(f"-> wrote {outdir/'blowup_sigma50_v3.csv'}")
    print(f"-> plots: *_v3.png in {outdir}")


if __name__ == "__main__":
    main()
