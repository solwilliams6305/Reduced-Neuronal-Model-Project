#!/usr/bin/env python3
"""
scaling_law_diagnostics.py

Immediate Test 1: is the Arrhenius barrier scaling exponent stable, or is the
observed ~2/3 just a finite-range crossover?

Input
-----
Use the output from arrhenius_ab_grid_checked.py, preferably:
    arrhenius_ab_summary_aggregated.csv
or, if you ran multiple seeds:
    arrhenius_ab_summary_by_seed.csv

What it tests
-------------
1. Fits B(delta) = C delta^alpha for each (a,b) pair and globally.
2. Compares fixed exponent models:
       alpha = 1/2, 2/3, 1, 3/2, 2
   against a free-alpha fit.
3. Computes leave-one-delta-out prediction error.
4. Computes sliding-window/free-window exponents when enough deltas exist.
5. Makes plots:
       barrier_vs_delta_fits.png
       exponent_by_pair.png
       model_comparison.png
       residuals_by_model.png
       window_exponents.png

Interpretation
--------------
If alpha is stable near 2/3 across pairs and delta windows, the 2/3 conjecture is
strengthened. If alpha drifts strongly with delta window, the 2/3 is probably a
finite-range effective exponent/crossover rather than a true asymptotic law.
"""

from __future__ import annotations

import argparse
import itertools
import math
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


FIXED_ALPHAS = [0.5, 2/3, 1.0, 1.5, 2.0]


def _pick_col(df: pd.DataFrame, candidates: list[str]) -> str:
    for c in candidates:
        if c in df.columns:
            return c
    raise ValueError(f"None of these columns found: {candidates}")


def load_barriers(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)

    # Support both aggregated and by-seed formats.
    b_col = _pick_col(df, ["B_meas_mean", "B_meas", "B_mean", "barrier_B", "B"])
    r2_col = None
    for c in ["arrhenius_r2_mean", "arrhenius_r2", "r2_mean", "R2"]:
        if c in df.columns:
            r2_col = c
            break

    required = ["a", "b", "delta"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Input missing required columns: {missing}")

    out = df.copy()
    out["B"] = pd.to_numeric(out[b_col], errors="coerce")
    out["delta"] = pd.to_numeric(out["delta"], errors="coerce")
    if r2_col is not None:
        out["arrhenius_r2"] = pd.to_numeric(out[r2_col], errors="coerce")
    else:
        out["arrhenius_r2"] = np.nan

    out = out[np.isfinite(out["B"]) & np.isfinite(out["delta"]) & (out["B"] > 0) & (out["delta"] > 0)]
    out = out.sort_values(["a", "b", "delta"]).reset_index(drop=True)
    if out.empty:
        raise ValueError("No usable positive B/delta rows found.")
    return out


def fit_fixed_alpha(delta: np.ndarray, B: np.ndarray, alpha: float) -> dict:
    """Least squares in original B-space for B = C delta^alpha."""
    x = delta ** alpha
    C = float(np.sum(x * B) / (np.sum(x * x) + 1e-300))
    pred = C * x
    resid = B - pred
    sse = float(np.sum(resid ** 2))
    rmse = math.sqrt(sse / max(len(B), 1))
    # log-space error too, useful for power laws.
    lerr = np.log(B) - np.log(np.maximum(pred, 1e-300))
    log_rmse = math.sqrt(float(np.mean(lerr ** 2)))
    return {"alpha": alpha, "C": C, "sse": sse, "rmse": rmse, "log_rmse": log_rmse}


def fit_free_alpha(delta: np.ndarray, B: np.ndarray) -> dict:
    """OLS in log-space: log B = log C + alpha log delta."""
    x = np.log(delta)
    y = np.log(B)
    X = np.vstack([x, np.ones_like(x)]).T
    alpha, logC = np.linalg.lstsq(X, y, rcond=None)[0]
    C = float(math.exp(logC))
    pred = C * delta ** alpha
    resid = B - pred
    sse = float(np.sum(resid ** 2))
    log_resid = y - (alpha * x + logC)
    log_sse = float(np.sum(log_resid ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2)) + 1e-300
    log_r2 = 1.0 - log_sse / ss_tot
    return {
        "alpha": float(alpha), "C": C, "sse": sse,
        "rmse": math.sqrt(sse / max(len(B), 1)),
        "log_rmse": math.sqrt(log_sse / max(len(B), 1)),
        "log_r2": float(log_r2),
    }


def loo_error(delta: np.ndarray, B: np.ndarray, alpha: float | None) -> dict:
    """Leave-one-delta-out log prediction error. alpha=None means free-alpha."""
    errs = []
    n = len(B)
    if n < 3:
        return {"loo_log_rmse": np.nan, "loo_rel_rmse": np.nan}
    for i in range(n):
        tr = np.ones(n, dtype=bool); tr[i] = False
        if alpha is None:
            if tr.sum() < 2:
                continue
            fit = fit_free_alpha(delta[tr], B[tr])
            pred = fit["C"] * delta[i] ** fit["alpha"]
        else:
            fit = fit_fixed_alpha(delta[tr], B[tr], alpha)
            pred = fit["C"] * delta[i] ** alpha
        errs.append((math.log(B[i]) - math.log(max(pred, 1e-300)), (B[i] - pred) / B[i]))
    if not errs:
        return {"loo_log_rmse": np.nan, "loo_rel_rmse": np.nan}
    loge = np.array([e[0] for e in errs])
    rele = np.array([e[1] for e in errs])
    return {
        "loo_log_rmse": float(np.sqrt(np.mean(loge ** 2))),
        "loo_rel_rmse": float(np.sqrt(np.mean(rele ** 2))),
    }


def fit_models_for_group(g: pd.DataFrame, group_name: str) -> pd.DataFrame:
    delta = g["delta"].to_numpy(float)
    B = g["B"].to_numpy(float)
    rows = []

    free = fit_free_alpha(delta, B)
    loo = loo_error(delta, B, alpha=None)
    rows.append({
        "group": group_name,
        "model": "free_alpha",
        "alpha": free["alpha"],
        "C": free["C"],
        "rmse": free["rmse"],
        "log_rmse": free["log_rmse"],
        "log_r2": free["log_r2"],
        "loo_log_rmse": loo["loo_log_rmse"],
        "loo_rel_rmse": loo["loo_rel_rmse"],
        "n": len(g),
    })

    for a in FIXED_ALPHAS:
        fit = fit_fixed_alpha(delta, B, a)
        loo = loo_error(delta, B, alpha=a)
        # log-space pseudo R^2 for fixed alpha with fitted C.
        pred = fit["C"] * delta ** a
        y = np.log(B)
        yp = np.log(np.maximum(pred, 1e-300))
        ss_res = float(np.sum((y - yp) ** 2))
        ss_tot = float(np.sum((y - y.mean()) ** 2)) + 1e-300
        rows.append({
            "group": group_name,
            "model": f"alpha_{a:.6g}",
            "alpha": a,
            "C": fit["C"],
            "rmse": fit["rmse"],
            "log_rmse": fit["log_rmse"],
            "log_r2": 1.0 - ss_res / ss_tot,
            "loo_log_rmse": loo["loo_log_rmse"],
            "loo_rel_rmse": loo["loo_rel_rmse"],
            "n": len(g),
        })
    return pd.DataFrame(rows)


def window_exponents(df: pd.DataFrame, min_points: int = 3) -> pd.DataFrame:
    rows = []
    for (a, b), g in df.groupby(["a", "b"]):
        g = g.sort_values("delta")
        vals = g[["delta", "B"]].to_numpy(float)
        n = len(vals)
        for L in range(min_points, n + 1):
            for i in range(0, n - L + 1):
                sub = vals[i:i+L]
                fit = fit_free_alpha(sub[:, 0], sub[:, 1])
                rows.append({
                    "a": a, "b": b,
                    "delta_min": sub[:, 0].min(),
                    "delta_max": sub[:, 0].max(),
                    "n": L,
                    "alpha": fit["alpha"],
                    "C": fit["C"],
                    "log_r2": fit["log_r2"],
                })
    return pd.DataFrame(rows)


def make_plots(df: pd.DataFrame, model_df: pd.DataFrame, win_df: pd.DataFrame, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    # Global fit curves.
    global_models = model_df[model_df["group"] == "global"].copy()
    dgrid = np.geomspace(df["delta"].min() * 0.9, df["delta"].max() * 1.1, 200)

    fig, ax = plt.subplots(figsize=(8, 5.5))
    for (a, b), g in df.groupby(["a", "b"]):
        ax.plot(g["delta"], g["B"], marker="o", linewidth=1, alpha=0.7, label=f"a={a:g}, b={b:g}")
    for _, r in global_models.iterrows():
        if r["model"] in ["free_alpha", "alpha_0.5", "alpha_0.666667", "alpha_1"]:
            label = f"{r['model']} (α={r['alpha']:.3g})"
            ax.plot(dgrid, r["C"] * dgrid ** r["alpha"], linestyle="--", linewidth=2, label=label)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Hopf deficit δ = I_H - I")
    ax.set_ylabel("measured Arrhenius barrier B")
    ax.set_title("Barrier scaling: data and global power-law fits")
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(outdir / "barrier_vs_delta_fits.png", dpi=180)
    plt.close(fig)

    # Free exponent by pair.
    pair_free = model_df[(model_df["model"] == "free_alpha") & (model_df["group"] != "global")].copy()
    if not pair_free.empty:
        fig, ax = plt.subplots(figsize=(8, 4.8))
        labels = pair_free["group"].to_list()
        ax.bar(np.arange(len(pair_free)), pair_free["alpha"])
        ax.axhline(2/3, linestyle="--", linewidth=1, label="2/3")
        ax.axhline(0.5, linestyle=":", linewidth=1, label="1/2")
        ax.set_xticks(np.arange(len(pair_free)))
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_ylabel("free fitted exponent α")
        ax.set_title("Free exponent by (a,b) pair")
        ax.legend()
        fig.tight_layout()
        fig.savefig(outdir / "exponent_by_pair.png", dpi=180)
        plt.close(fig)

    # Model comparison by global LOO.
    glob = model_df[model_df["group"] == "global"].sort_values("loo_log_rmse")
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.bar(np.arange(len(glob)), glob["loo_log_rmse"])
    ax.set_xticks(np.arange(len(glob)))
    ax.set_xticklabels(glob["model"], rotation=45, ha="right")
    ax.set_ylabel("leave-one-delta-out log RMSE")
    ax.set_title("Global model comparison")
    fig.tight_layout()
    fig.savefig(outdir / "model_comparison.png", dpi=180)
    plt.close(fig)

    # Residual plot for fixed candidates.
    fig, ax = plt.subplots(figsize=(8, 5.2))
    for m in ["alpha_0.5", "alpha_0.666667", "alpha_1", "free_alpha"]:
        r = global_models[global_models["model"] == m]
        if r.empty:
            continue
        r = r.iloc[0]
        pred = r["C"] * df["delta"].to_numpy(float) ** r["alpha"]
        resid = np.log(df["B"].to_numpy(float)) - np.log(np.maximum(pred, 1e-300))
        ax.scatter(df["delta"], resid, label=f"{m} α={r['alpha']:.3g}", alpha=0.8)
    ax.axhline(0, linewidth=1)
    ax.set_xscale("log")
    ax.set_xlabel("δ")
    ax.set_ylabel("log residual")
    ax.set_title("Residuals against global scaling models")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "residuals_by_model.png", dpi=180)
    plt.close(fig)

    # Window exponents if enough delta points exist.
    if not win_df.empty:
        fig, ax = plt.subplots(figsize=(8, 5.2))
        for (a, b), g in win_df.groupby(["a", "b"]):
            centers = np.sqrt(g["delta_min"] * g["delta_max"])
            ax.scatter(centers, g["alpha"], label=f"a={a:g}, b={b:g}", alpha=0.75)
        ax.axhline(2/3, linestyle="--", linewidth=1, label="2/3")
        ax.axhline(0.5, linestyle=":", linewidth=1, label="1/2")
        ax.set_xscale("log")
        ax.set_xlabel("geometric centre of δ window")
        ax.set_ylabel("window-fitted exponent α")
        ax.set_title("Does the exponent drift with δ-window?")
        ax.legend(fontsize=7, ncol=2)
        fig.tight_layout()
        fig.savefig(outdir / "window_exponents.png", dpi=180)
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=str, default="results/arrhenius_ab_grid_checked/arrhenius_ab_summary_aggregated.csv")
    ap.add_argument("--outdir", type=str, default="results/scaling_law_diagnostics")
    ap.add_argument("--min-r2", type=float, default=None,
                    help="Optional: drop Arrhenius barrier estimates with R^2 below this.")
    args = ap.parse_args()

    input_path = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = load_barriers(input_path)
    if args.min_r2 is not None and "arrhenius_r2" in df.columns:
        before = len(df)
        df = df[(df["arrhenius_r2"].isna()) | (df["arrhenius_r2"] >= args.min_r2)].copy()
        print(f"Filtered by min R^2={args.min_r2}: {before} -> {len(df)} rows")

    # If multiple seeds are still present, average them before fitting by pair/delta.
    df = (
        df.groupby(["a", "b", "delta"], as_index=False)
          .agg(B=("B", "mean"), B_sd=("B", "std"), arrhenius_r2=("arrhenius_r2", "mean"))
          .sort_values(["a", "b", "delta"])
    )

    model_parts = [fit_models_for_group(df, "global")]
    for (a, b), g in df.groupby(["a", "b"]):
        model_parts.append(fit_models_for_group(g, f"a={a:g},b={b:g}"))
    model_df = pd.concat(model_parts, ignore_index=True)

    win_df = window_exponents(df, min_points=3)

    df.to_csv(outdir / "barrier_data_used.csv", index=False)
    model_df.to_csv(outdir / "power_model_comparison.csv", index=False)
    win_df.to_csv(outdir / "window_exponents.csv", index=False)

    make_plots(df, model_df, win_df, outdir)

    print("\n================ SCALING DIAGNOSTIC ================")
    print(f"Input: {input_path}")
    print(f"Rows used after aggregation: {len(df)}")
    print("\nGlobal model ranking by leave-one-delta-out log RMSE:")
    cols = ["model", "alpha", "C", "log_rmse", "loo_log_rmse", "loo_rel_rmse", "log_r2"]
    print(model_df[model_df["group"] == "global"][cols].sort_values("loo_log_rmse").to_string(index=False))

    pair_free = model_df[(model_df["model"] == "free_alpha") & (model_df["group"] != "global")]
    if not pair_free.empty:
        print("\nFree exponent by pair:")
        print(pair_free[["group", "alpha", "C", "log_r2"]].to_string(index=False))
        print(f"\nMean pair exponent = {pair_free['alpha'].mean():.4f}")
        print(f"SD pair exponent   = {pair_free['alpha'].std(ddof=1):.4f}" if len(pair_free) > 1 else "")

    print(f"\nOutputs written to: {outdir}")
    print("====================================================")


if __name__ == "__main__":
    main()
