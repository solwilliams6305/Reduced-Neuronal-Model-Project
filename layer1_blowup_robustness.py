#!/usr/bin/env python3
"""
layer1_blowup_robustness.py

Layer 1 script:
  Robustness tests for the blow-up normal form

      dV = (V^2 - W)dT + eta dB_T
      dW = -lambda dT

Goal:
  Check whether the window-hit law

      eta50 ~ C(q, window) lambda^beta

  keeps beta ≈ 1/2 across:
      - target fractions q
      - canard-window choices [lo W*, hi W*]
      - optional lambda fitting ranges

Outputs:
  results/layer1_blowup_robustness/
      layer1_all_results.csv
      layer1_eta50_summary.csv
      layer1_fit_summary.csv
      layer1_beta_by_setting.png
      layer1_eta50_grid.png
      layer1_collapse_examples.png

Recommended run:

python3 layer1_blowup_robustness.py \
  --overwrite \
  --target-fracs 0.5 0.75 0.9 \
  --windows 0.5:2.0 0.25:2.0 0.5:3.0 \
  --lambdas 0.004 0.006 0.008 0.012 0.016 0.024 0.032 0.048 0.064 \
  --etas 0 0.02 0.035 0.05 0.07 0.10 0.14 0.20 0.28 0.40 0.60 0.85 1.10 \
  --n 2000 \
  --fit-lambda-min 0.008 \
  --outdir results/layer1_blowup_robustness

Notes:
  This is intentionally standalone. It does not import your V4 script.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


@dataclass
class Setting:
    q: float
    window_lo: float
    window_hi: float


def parse_window(s: str) -> tuple[float, float]:
    if ":" not in s:
        raise argparse.ArgumentTypeError("Window must be formatted like lo:hi, e.g. 0.5:2.0")
    lo, hi = s.split(":", 1)
    lo, hi = float(lo), float(hi)
    if lo <= 0 or hi <= lo:
        raise argparse.ArgumentTypeError("Need 0 < lo < hi")
    return lo, hi


def interp_x_at_y(x: np.ndarray, y: np.ndarray, target: float = 0.5) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    order = np.argsort(x)
    x, y = x[order], y[order]
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]

    if len(x) < 2 or np.max(y) < target:
        return np.nan
    if np.min(y) >= target:
        return float(x[0])

    for i in range(len(x) - 1):
        if y[i] == target:
            return float(x[i])
        if (y[i] - target) * (y[i + 1] - target) <= 0 and y[i] != y[i + 1]:
            a = (target - y[i]) / (y[i + 1] - y[i])
            return float(x[i] + a * (x[i + 1] - x[i]))
    return np.nan


def fit_power(lam: np.ndarray, eta50: np.ndarray) -> dict:
    lam = np.asarray(lam, dtype=float)
    eta50 = np.asarray(eta50, dtype=float)
    mask = np.isfinite(lam) & np.isfinite(eta50) & (lam > 0) & (eta50 > 0)
    lam, eta50 = lam[mask], eta50[mask]
    if len(lam) < 2:
        return {"C": np.nan, "beta": np.nan, "r2": np.nan, "n_fit": len(lam)}

    x, y = np.log(lam), np.log(eta50)
    beta, logC = np.polyfit(x, y, 1)
    yhat = logC + beta * x
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = np.nan if ss_tot == 0 else 1 - ss_res / ss_tot
    return {"C": float(np.exp(logC)), "beta": float(beta), "r2": float(r2), "n_fit": int(len(lam))}


def simulate_prob(
    lam: float,
    eta: float,
    q: float,
    window_lo: float,
    window_hi: float,
    n: int,
    dt: float,
    seed: int,
    W_start_factor: float,
    W_end: float,
    W_min: float,
    V_abs_clip: float,
) -> dict:
    rng = np.random.default_rng(seed)

    W_star = lam ** (2.0 / 3.0)
    W_start = W_start_factor * W_star
    W_lo = max(window_lo * W_star, W_min)
    W_hi = window_hi * W_star

    T_total = (W_start - W_end) / lam
    n_steps = int(np.ceil(T_total / dt))
    dt_eff = T_total / n_steps

    V = -np.sqrt(np.full(n, W_start))
    active = np.ones(n, dtype=bool)
    hit = np.zeros(n, dtype=bool)

    sqrt_dt_eta = eta * math.sqrt(dt_eff)
    W = W_start

    for _ in range(n_steps):
        if not np.any(active):
            break

        idx = active
        dB = rng.normal(size=np.count_nonzero(idx))
        V[idx] += (V[idx] * V[idx] - W) * dt_eff + sqrt_dt_eta * dB

        # Numerical hygiene: if a path has gone huge, it has escaped the local chart.
        # This prevents overflow without changing hit decisions inside the intended window.
        V[idx] = np.clip(V[idx], -V_abs_clip, V_abs_clip)

        W -= lam * dt_eff
        if W <= W_min:
            continue

        if W_lo <= W <= W_hi:
            threshold = q * math.sqrt(W)
            hit_now = active & (V >= threshold)
            hit |= hit_now
            active[hit_now] = False

    return {
        "q": q,
        "window_lo": window_lo,
        "window_hi": window_hi,
        "lambda": lam,
        "eta": eta,
        "hit_prob": float(np.mean(hit)),
        "W_star": W_star,
        "W_lo": W_lo,
        "W_hi": W_hi,
        "n": n,
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", default="results/layer1_blowup_robustness")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--target-fracs", type=float, nargs="+", default=[0.5, 0.75, 0.9])
    p.add_argument("--windows", type=parse_window, nargs="+", default=[(0.5, 2.0), (0.25, 2.0), (0.5, 3.0)])
    p.add_argument("--lambdas", type=float, nargs="+", required=True)
    p.add_argument("--etas", type=float, nargs="+", required=True)
    p.add_argument("--n", type=int, default=2000)
    p.add_argument("--dt", type=float, default=0.002)
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--W-start-factor", type=float, default=8.0)
    p.add_argument("--W-end", type=float, default=-0.05)
    p.add_argument("--W-min", type=float, default=1e-8)
    p.add_argument("--V-abs-clip", type=float, default=1e6)
    p.add_argument("--fit-lambda-min", type=float, default=None)
    p.add_argument("--fit-lambda-max", type=float, default=None)
    args = p.parse_args()

    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()) and not args.overwrite:
        raise SystemExit(f"{outdir} exists and is nonempty. Use --overwrite.")
    outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    counter = 0
    settings = [Setting(q=q, window_lo=w[0], window_hi=w[1]) for q in args.target_fracs for w in args.windows]

    print("\n=== Layer 1: blow-up robustness ===")
    print(f"settings={len(settings)}, lambdas={len(args.lambdas)}, etas={len(args.etas)}, n={args.n}\n")

    for s in settings:
        print(f"\n--- q={s.q:g}, window=[{s.window_lo:g}W*, {s.window_hi:g}W*] ---")
        for lam in args.lambdas:
            for eta in args.etas:
                counter += 1
                row = simulate_prob(
                    lam=lam,
                    eta=eta,
                    q=s.q,
                    window_lo=s.window_lo,
                    window_hi=s.window_hi,
                    n=args.n,
                    dt=args.dt,
                    seed=args.seed + 100000 * counter,
                    W_start_factor=args.W_start_factor,
                    W_end=args.W_end,
                    W_min=args.W_min,
                    V_abs_clip=args.V_abs_clip,
                )
                rows.append(row)
                print(f"q={s.q:g}, win={s.window_lo:g}:{s.window_hi:g}, lambda={lam:g}, eta={eta:g}, P={row['hit_prob']:.3f}")

    all_results = pd.DataFrame(rows)

    eta50_rows = []
    for keys, g in all_results.groupby(["q", "window_lo", "window_hi", "lambda"]):
        q, wlo, whi, lam = keys
        g = g.sort_values("eta")
        eta50 = interp_x_at_y(g["eta"].to_numpy(), g["hit_prob"].to_numpy(), 0.5)
        first = g.iloc[0]
        eta50_rows.append({
            "q": q,
            "window_lo": wlo,
            "window_hi": whi,
            "lambda": lam,
            "eta50": eta50,
            "W_star": first["W_star"],
            "W_lo": first["W_lo"],
            "W_hi": first["W_hi"],
        })
    eta50_df = pd.DataFrame(eta50_rows).sort_values(["q", "window_lo", "window_hi", "lambda"])

    fit_rows = []
    for keys, g in eta50_df.groupby(["q", "window_lo", "window_hi"]):
        q, wlo, whi = keys

        fit_all = fit_power(g["lambda"], g["eta50"])
        fit_rows.append({"q": q, "window_lo": wlo, "window_hi": whi, "fit": "all", **fit_all})

        gr = g.copy()
        if args.fit_lambda_min is not None:
            gr = gr[gr["lambda"] >= args.fit_lambda_min]
        if args.fit_lambda_max is not None:
            gr = gr[gr["lambda"] <= args.fit_lambda_max]
        fit_range = fit_power(gr["lambda"], gr["eta50"])
        fit_rows.append({"q": q, "window_lo": wlo, "window_hi": whi, "fit": "range", **fit_range})

    fit_df = pd.DataFrame(fit_rows)

    all_results.to_csv(outdir / "layer1_all_results.csv", index=False)
    eta50_df.to_csv(outdir / "layer1_eta50_summary.csv", index=False)
    fit_df.to_csv(outdir / "layer1_fit_summary.csv", index=False)

    # Plot beta by setting.
    range_fits = fit_df[fit_df["fit"] == "range"].copy()
    if range_fits["beta"].isna().all():
        range_fits = fit_df[fit_df["fit"] == "all"].copy()

    labels = [f"q={r.q:g}\n[{r.window_lo:g},{r.window_hi:g}]" for r in range_fits.itertuples()]
    x = np.arange(len(range_fits))

    plt.figure(figsize=(max(10, 1.2 * len(x)), 6))
    plt.bar(x, range_fits["beta"])
    plt.axhline(0.5, linestyle="--", color="black", label=r"$1/2$")
    plt.axhline(2/3, linestyle=":", color="black", label=r"$2/3$")
    plt.xticks(x, labels, rotation=0)
    plt.ylabel(r"fitted $\beta$")
    plt.title(r"Layer 1: fitted exponent in $\eta_{50}\sim C\lambda^\beta$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "layer1_beta_by_setting.png", dpi=180)
    plt.close()

    # Grid eta50 plots.
    plt.figure(figsize=(11, 7))
    for keys, g in eta50_df.groupby(["q", "window_lo", "window_hi"]):
        q, wlo, whi = keys
        plt.loglog(g["lambda"], g["eta50"], marker="o", label=f"q={q:g}, win={wlo:g}:{whi:g}")
    lam_grid = np.linspace(min(args.lambdas), max(args.lambdas), 200)
    # purely visual slopes, median-scaled
    med_lam = np.median(args.lambdas)
    med_eta = np.nanmedian(eta50_df["eta50"])
    plt.loglog(lam_grid, med_eta * (lam_grid / med_lam) ** 0.5, "--", label=r"$\lambda^{1/2}$ slope")
    plt.loglog(lam_grid, med_eta * (lam_grid / med_lam) ** (2/3), ":", label=r"$\lambda^{2/3}$ slope")
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$\eta_{50}$")
    plt.title("Layer 1: eta50 curves across target/window choices")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outdir / "layer1_eta50_grid.png", dpi=180)
    plt.close()

    # Collapse examples for q=0.75, first window if present.
    ex = eta50_df[(eta50_df["q"] == 0.75)]
    if len(ex) == 0:
        ex = eta50_df.copy()
    first_keys = ex.groupby(["q", "window_lo", "window_hi"]).size().index[0]
    ex = eta50_df[
        (eta50_df["q"] == first_keys[0]) &
        (eta50_df["window_lo"] == first_keys[1]) &
        (eta50_df["window_hi"] == first_keys[2])
    ].copy()
    C_half = np.nanmedian(ex["eta50"] / np.sqrt(ex["lambda"]))
    C_23 = np.nanmedian(ex["eta50"] / (ex["lambda"] ** (2/3)))

    plt.figure(figsize=(10, 6))
    plt.semilogx(ex["lambda"], ex["eta50"] / (C_half * np.sqrt(ex["lambda"])), marker="o", label=r"$1/2$ collapse")
    plt.semilogx(ex["lambda"], ex["eta50"] / (C_23 * ex["lambda"] ** (2/3)), marker="s", label=r"$2/3$ collapse")
    plt.axhline(1.0, linestyle="--", color="black")
    plt.xlabel(r"$\lambda$")
    plt.ylabel("measured / prediction")
    plt.title(f"Layer 1 collapse example: q={first_keys[0]:g}, window={first_keys[1]:g}:{first_keys[2]:g}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "layer1_collapse_examples.png", dpi=180)
    plt.close()

    print("\n=== Layer 1 fit summary ===")
    print(fit_df.to_string(index=False))
    print(f"\nSaved outputs to {outdir}")


if __name__ == "__main__":
    main()
