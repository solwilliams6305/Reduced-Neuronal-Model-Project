#!/usr/bin/env python3
"""
layer2c_full_fhn_hit_location_bins.py

Layer 2c:
  Full autonomous stochastic FHN first-hit location classifier with fine W/W* bins.

Why:
  Layer 2b showed that full FHN any-hit probability can cross 50%, but the first hit
  is often outside the narrow normal-form window [0.5 W*, 2 W*].

This script classifies the FIRST local hit into finer bins in

      R = W_blow / W*,

where

      W_blow = (w - w_fold_left) / eps^(2/3),
      W* = lambda^(2/3).

Default bins:
    R > 4
    2 < R <= 4
    0.5 <= R <= 2
    0.1 <= R < 0.5
    0 < R < 0.1
    no hit

It also supports a broad-window hit probability:

    broad_window = broad_lo <= R <= broad_hi

Default:
    broad_window = [0.1, 4.0]

Recommended run:

python3 layer2c_full_fhn_hit_location_bins.py \
  --overwrite \
  --a 0.7 \
  --b 0.8 \
  --delta-Is 0.01 0.02 0.04 \
  --epsilons 0.01 0.02 0.04 0.08 \
  --sigma-mults 0 0.5 1.0 1.5 2.0 2.5 3.0 4.0 5.0 6.5 8.0 \
  --n 1000 \
  --outdir results/layer2c_full_fhn_hit_location_bins

Wider noise run if needed:

python3 layer2c_full_fhn_hit_location_bins.py \
  --overwrite \
  --a 0.7 \
  --b 0.8 \
  --delta-Is 0.01 0.02 0.04 \
  --epsilons 0.01 0.02 0.04 0.08 \
  --sigma-mults 0 1 2 3 4 5 6.5 8 10 12 15 \
  --n 1000 \
  --outdir results/layer2c_full_fhn_hit_location_bins_wide

Outputs:
  layer2c_path_events.csv
  layer2c_probability_summary.csv
  layer2c_sigma50_summary.csv
  layer2c_hit_bin_stack_by_scaled_noise.png
  layer2c_hit_bin_stack_by_sigma_mult.png
  layer2c_broad_vs_narrow_collapse.png
  layer2c_median_hit_location.png
  layer2c_sigma50_any_broad_vs_half_scale.png
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


DEFAULT_BIN_EDGES = [0.0, 0.1, 0.5, 2.0, 4.0, np.inf]
DEFAULT_BIN_LABELS = [
    "R_0_to_0p1",
    "R_0p1_to_0p5",
    "R_0p5_to_2",
    "R_2_to_4",
    "R_gt_4",
]


def cubic_roots_for_w_I(w: float, I: float) -> np.ndarray:
    coeffs = [1.0, 0.0, -3.0, 3.0 * (w - I)]
    roots = np.roots(coeffs)
    real = np.real(roots[np.abs(np.imag(roots)) < 1e-8])
    return np.sort(real)


def left_attracting_v(w: float, I: float) -> float:
    roots = cubic_roots_for_w_I(w, I)
    if len(roots) == 0:
        raise RuntimeError("Cubic unexpectedly had no real roots.")
    return float(roots[0])


def interp_x_at_y(x: np.ndarray, y: np.ndarray, target: float = 0.5) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    order = np.argsort(x)
    x, y = x[order], y[order]
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]

    if len(x) < 2 or np.nanmax(y) < target:
        return np.nan
    if np.nanmin(y) >= target:
        return float(x[0])

    for i in range(len(x) - 1):
        if y[i] == target:
            return float(x[i])
        if (y[i] - target) * (y[i + 1] - target) <= 0 and y[i] != y[i + 1]:
            a = (target - y[i]) / (y[i + 1] - y[i])
            return float(x[i] + a * (x[i + 1] - x[i]))
    return np.nan


def classify_bin(R: float, edges: list[float], labels: list[str]) -> str:
    if not np.isfinite(R) or R <= 0:
        return "no_hit"
    for lo, hi, label in zip(edges[:-1], edges[1:], labels):
        # Use (lo, hi] for positive bins, except first positive includes lo=0 as R>0.
        if R > lo and R <= hi:
            return label
    return labels[-1]


def parse_edges(s: str) -> list[float]:
    vals = []
    for part in s.split(","):
        part = part.strip()
        if part.lower() in {"inf", "infty", "infinity"}:
            vals.append(np.inf)
        else:
            vals.append(float(part))
    if len(vals) < 3:
        raise argparse.ArgumentTypeError("Need at least 3 bin edges.")
    if vals[0] != 0:
        raise argparse.ArgumentTypeError("First bin edge should be 0.")
    if not all(vals[i] < vals[i+1] for i in range(len(vals)-1)):
        raise argparse.ArgumentTypeError("Bin edges must be strictly increasing.")
    return vals


def labels_from_edges(edges: list[float]) -> list[str]:
    def fmt(x):
        if np.isinf(x):
            return "inf"
        return str(x).replace(".", "p").replace("-", "m")
    return [f"R_{fmt(lo)}_to_{fmt(hi)}" for lo, hi in zip(edges[:-1], edges[1:])]


def simulate_event_batch(
    *,
    a: float,
    b: float,
    epsilon: float,
    delta_I: float,
    sigma_mult: float,
    n: int,
    dt: float,
    seed: int,
    q: float,
    narrow_lo: float,
    narrow_hi: float,
    broad_lo: float,
    broad_hi: float,
    bin_edges: list[float],
    bin_labels: list[str],
    W_start_factor: float,
    W_end_factor: float,
    T_factor: float,
    V_abs_clip: float,
    local_chart_V_abs: float,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)

    I_fold_left = 2.0 / 3.0 + (a - 1.0) / b
    I = I_fold_left + delta_I
    w_fold_left = I - 2.0 / 3.0
    lam = b * delta_I
    W_star = lam ** (2.0 / 3.0)

    sigma = sigma_mult * math.sqrt(epsilon) * math.sqrt(lam)

    W_start_blow = W_start_factor * W_star
    W_end_blow = W_end_factor * W_star

    w_start = w_fold_left + epsilon ** (2.0 / 3.0) * W_start_blow
    v_start = left_attracting_v(w_start, I)

    v = np.full(n, v_start, dtype=float)
    w = np.full(n, w_start, dtype=float)
    active = np.ones(n, dtype=bool)

    first_hit = np.zeros(n, dtype=bool)
    first_hit_t = np.full(n, np.nan)
    first_hit_V = np.full(n, np.nan)
    first_hit_W = np.full(n, np.nan)
    first_hit_R = np.full(n, np.nan)
    first_hit_bin = np.array(["no_hit"] * n, dtype=object)
    first_hit_in_narrow = np.zeros(n, dtype=bool)
    first_hit_in_broad = np.zeros(n, dtype=bool)

    in_narrow_prev = np.zeros(n, dtype=bool)
    in_broad_prev = np.zeros(n, dtype=bool)
    narrow_entry_count = np.zeros(n, dtype=int)
    broad_entry_count = np.zeros(n, dtype=int)

    min_R = np.full(n, np.inf)
    max_R = np.full(n, -np.inf)

    slow_speed0 = abs(v_start + a - b * w_start)
    T_base = abs(epsilon ** (2.0 / 3.0) * (W_start_blow - W_end_blow)) / max(slow_speed0, 1e-8)
    T_max = T_factor * T_base
    n_steps = int(np.ceil(T_max / dt))

    sqrt_dt = math.sqrt(dt)

    for step in range(n_steps):
        t = step * dt

        W_blow = (w - w_fold_left) / (epsilon ** (2.0 / 3.0))
        V_blow = (v + 1.0) / (epsilon ** (1.0 / 3.0))
        R = W_blow / W_star

        finite = np.isfinite(R)
        min_R[finite] = np.minimum(min_R[finite], R[finite])
        max_R[finite] = np.maximum(max_R[finite], R[finite])

        in_narrow = active & (R >= narrow_lo) & (R <= narrow_hi)
        in_broad = active & (R >= broad_lo) & (R <= broad_hi)

        narrow_entry_count[in_narrow & (~in_narrow_prev)] += 1
        broad_entry_count[in_broad & (~in_broad_prev)] += 1
        in_narrow_prev = in_narrow.copy()
        in_broad_prev = in_broad.copy()

        good_W = active & (W_blow > 0) & np.isfinite(V_blow) & np.isfinite(W_blow)
        threshold = np.full(n, np.inf)
        threshold[good_W] = q * np.sqrt(W_blow[good_W])

        hit_now = good_W & (~first_hit) & (V_blow >= threshold)

        if np.any(hit_now):
            first_hit[hit_now] = True
            first_hit_t[hit_now] = t
            first_hit_V[hit_now] = V_blow[hit_now]
            first_hit_W[hit_now] = W_blow[hit_now]
            first_hit_R[hit_now] = R[hit_now]
            first_hit_in_narrow[hit_now] = (R[hit_now] >= narrow_lo) & (R[hit_now] <= narrow_hi)
            first_hit_in_broad[hit_now] = (R[hit_now] >= broad_lo) & (R[hit_now] <= broad_hi)

            for j in np.where(hit_now)[0]:
                first_hit_bin[j] = classify_bin(float(R[j]), bin_edges, bin_labels)

            active[hit_now] = False

        too_far_V = np.abs(V_blow) > local_chart_V_abs
        below_chart = R < W_end_factor
        above_chart = R > max(W_start_factor * 1.5, broad_hi * 2.0, narrow_hi * 4.0)
        inactive = active & (too_far_V | below_chart | above_chart)
        active[inactive] = False

        if not np.any(active):
            break

        idx = active
        vv = v[idx]
        ww = w[idx]

        drift_v = (vv - vv**3 / 3.0 - ww + I) / epsilon
        drift_w = vv + a - b * ww

        vv = vv + drift_v * dt + (sigma / math.sqrt(epsilon)) * sqrt_dt * rng.normal(size=np.count_nonzero(idx))
        ww = ww + drift_w * dt

        vv = np.clip(vv, -V_abs_clip, V_abs_clip)

        v[idx] = vv
        w[idx] = ww

    min_R[~np.isfinite(min_R)] = np.nan
    max_R[~np.isfinite(max_R)] = np.nan

    return pd.DataFrame({
        "a": a,
        "b": b,
        "epsilon": epsilon,
        "delta_I": delta_I,
        "lambda": lam,
        "I_fold_left": I_fold_left,
        "I": I,
        "w_fold_left": w_fold_left,
        "sigma_mult": sigma_mult,
        "sigma": sigma,
        "q": q,
        "narrow_lo": narrow_lo,
        "narrow_hi": narrow_hi,
        "broad_lo": broad_lo,
        "broad_hi": broad_hi,
        "W_star": W_star,
        "path": np.arange(n),
        "first_hit": first_hit,
        "first_hit_t": first_hit_t,
        "first_hit_V": first_hit_V,
        "first_hit_W": first_hit_W,
        "first_hit_R": first_hit_R,
        "first_hit_bin": first_hit_bin,
        "first_hit_in_narrow": first_hit_in_narrow,
        "first_hit_in_broad": first_hit_in_broad,
        "narrow_entry_count": narrow_entry_count,
        "broad_entry_count": broad_entry_count,
        "min_R": min_R,
        "max_R": max_R,
        "v_start": v_start,
        "w_start": w_start,
        "T_max": T_max,
        "n_steps": n_steps,
    })


def make_probability_summary(events: pd.DataFrame, bin_labels: list[str]) -> pd.DataFrame:
    group_cols = ["epsilon", "delta_I", "lambda", "sigma_mult", "sigma"]
    rows = []

    for keys, g in events.groupby(group_cols):
        row = dict(zip(group_cols, keys))
        row["n"] = len(g)

        row["p_any_hit"] = float(np.mean(g["first_hit"]))
        row["p_narrow_hit"] = float(np.mean(g["first_hit_in_narrow"]))
        row["p_broad_hit"] = float(np.mean(g["first_hit_in_broad"]))
        row["p_no_hit"] = float(np.mean(~g["first_hit"]))

        for label in bin_labels:
            row[f"p_{label}"] = float(np.mean(g["first_hit_bin"] == label))

        row["p_entered_narrow"] = float(np.mean(g["narrow_entry_count"] > 0))
        row["p_entered_broad"] = float(np.mean(g["broad_entry_count"] > 0))
        row["p_multi_narrow"] = float(np.mean(g["narrow_entry_count"] >= 2))
        row["p_multi_broad"] = float(np.mean(g["broad_entry_count"] >= 2))

        row["median_hit_R"] = float(np.nanmedian(g["first_hit_R"])) if np.any(g["first_hit"]) else np.nan
        row["q25_hit_R"] = float(np.nanpercentile(g["first_hit_R"], 25)) if np.any(g["first_hit"]) else np.nan
        row["q75_hit_R"] = float(np.nanpercentile(g["first_hit_R"], 75)) if np.any(g["first_hit"]) else np.nan

        rows.append(row)

    out = pd.DataFrame(rows).sort_values(group_cols)
    out["scale_half"] = out["sigma"] / (np.sqrt(out["epsilon"]) * np.sqrt(out["lambda"]))
    out["scale_twothirds"] = out["sigma"] / (np.sqrt(out["epsilon"]) * (out["lambda"] ** (2.0 / 3.0)))
    out["scale_eps34"] = out["sigma"] / (out["epsilon"] ** 0.75)
    return out


def make_sigma50_summary(prob: pd.DataFrame, bin_labels: list[str]) -> pd.DataFrame:
    rows = []
    pcols = ["p_any_hit", "p_broad_hit", "p_narrow_hit"] + [f"p_{label}" for label in bin_labels]

    for keys, g in prob.groupby(["epsilon", "delta_I", "lambda"]):
        eps, dI, lam = keys
        g = g.sort_values("sigma")
        row = {"epsilon": eps, "delta_I": dI, "lambda": lam}
        for pcol in pcols:
            row[pcol.replace("p_", "sigma50_")] = interp_x_at_y(g["sigma"].to_numpy(), g[pcol].to_numpy(), 0.5)

        row["pred_half_scale"] = math.sqrt(eps) * math.sqrt(lam)
        row["pred_twothirds_scale"] = math.sqrt(eps) * (lam ** (2.0 / 3.0))
        row["pred_eps34_scale"] = eps ** 0.75
        rows.append(row)

    out = pd.DataFrame(rows).sort_values(["epsilon", "delta_I"])
    for col in [c for c in out.columns if c.startswith("sigma50_")]:
        out[f"{col}_over_half"] = out[col] / out["pred_half_scale"]
        out[f"{col}_over_twothirds"] = out[col] / out["pred_twothirds_scale"]
        out[f"{col}_over_eps34"] = out[col] / out["pred_eps34_scale"]
    return out


def plot_outputs(prob: pd.DataFrame, sigma50: pd.DataFrame, bin_labels: list[str], outdir: Path) -> None:
    # Stacked hit-bin composition by half-law scaled noise.
    pooled = prob.groupby("scale_half")[[f"p_{b}" for b in bin_labels] + ["p_no_hit"]].mean().reset_index().sort_values("scale_half")
    x = np.arange(len(pooled))
    bottom = np.zeros(len(pooled))

    plt.figure(figsize=(13, 7))
    for col in [f"p_{b}" for b in reversed(bin_labels)] + ["p_no_hit"]:
        # reversed puts large-R bins visually first
        label = col[2:] if col.startswith("p_") else col
        plt.bar(x, pooled[col], bottom=bottom, label=label)
        bottom += pooled[col].to_numpy()
    plt.xticks(x, [f"{s:g}" for s in pooled["scale_half"]], rotation=45)
    plt.xlabel(r"$\sigma/(\sqrt{\epsilon}\lambda^{1/2})$")
    plt.ylabel("pooled fraction")
    plt.title("Layer 2c: first-hit location bins by half-law scaled noise")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outdir / "layer2c_hit_bin_stack_by_scaled_noise.png", dpi=180)
    plt.close()

    # Stacked by sigma_mult.
    pooled_m = prob.groupby("sigma_mult")[[f"p_{b}" for b in bin_labels] + ["p_no_hit"]].mean().reset_index().sort_values("sigma_mult")
    x = np.arange(len(pooled_m))
    bottom = np.zeros(len(pooled_m))

    plt.figure(figsize=(13, 7))
    for col in [f"p_{b}" for b in reversed(bin_labels)] + ["p_no_hit"]:
        label = col[2:] if col.startswith("p_") else col
        plt.bar(x, pooled_m[col], bottom=bottom, label=label)
        bottom += pooled_m[col].to_numpy()
    plt.xticks(x, [f"{s:g}" for s in pooled_m["sigma_mult"]])
    plt.xlabel("sigma multiplier")
    plt.ylabel("pooled fraction")
    plt.title("Layer 2c: first-hit location bins by sigma multiplier")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outdir / "layer2c_hit_bin_stack_by_sigma_mult.png", dpi=180)
    plt.close()

    # Narrow vs broad vs any collapse.
    plt.figure(figsize=(11, 7))
    for pcol, label in [
        ("p_any_hit", "any W>0 hit"),
        ("p_broad_hit", "broad-window hit"),
        ("p_narrow_hit", "narrow-window hit"),
    ]:
        pooled = prob.groupby("scale_half")[pcol].mean().reset_index().sort_values("scale_half")
        plt.plot(pooled["scale_half"], pooled[pcol], marker="o", label=label)
    plt.axhline(0.5, linestyle="--", color="black")
    plt.xlabel(r"$\sigma/(\sqrt{\epsilon}\lambda^{1/2})$")
    plt.ylabel("pooled probability")
    plt.title("Layer 2c: broad vs narrow full-FHN local-hit probabilities")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "layer2c_broad_vs_narrow_collapse.png", dpi=180)
    plt.close()

    # Median hit location.
    plt.figure(figsize=(11, 7))
    for keys, g in prob.groupby(["epsilon", "delta_I"]):
        eps, dI = keys
        g = g.sort_values("scale_half")
        plt.plot(g["scale_half"], g["median_hit_R"], marker="o", alpha=0.8, label=f"eps={eps:g}, dI={dI:g}")
    plt.axhspan(0.5, 2.0, alpha=0.15, label="narrow [0.5,2]")
    plt.axhspan(0.1, 4.0, alpha=0.08, label="broad [0.1,4]")
    plt.yscale("log")
    plt.xlabel(r"$\sigma/(\sqrt{\epsilon}\lambda^{1/2})$")
    plt.ylabel(r"median first-hit $R=W/W_*$")
    plt.title("Layer 2c: median first-hit location")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outdir / "layer2c_median_hit_location.png", dpi=180)
    plt.close()

    # Sigma50 any/broad vs half scale.
    clean = sigma50.dropna(subset=["sigma50_any_hit"])
    if len(clean) >= 2:
        plt.figure(figsize=(10, 7))
        plt.loglog(clean["pred_half_scale"], clean["sigma50_any_hit"], "o", label="any hit")
        clean_b = sigma50.dropna(subset=["sigma50_broad_hit"])
        if len(clean_b):
            plt.loglog(clean_b["pred_half_scale"], clean_b["sigma50_broad_hit"], "s", label="broad hit")
        x_min = min(clean["pred_half_scale"].min(), clean_b["pred_half_scale"].min() if len(clean_b) else clean["pred_half_scale"].min())
        x_max = max(clean["pred_half_scale"].max(), clean_b["pred_half_scale"].max() if len(clean_b) else clean["pred_half_scale"].max())
        xgrid = np.linspace(x_min, x_max, 100)
        C = np.nanmedian(clean["sigma50_any_hit"] / clean["pred_half_scale"])
        plt.loglog(xgrid, C * xgrid, "--", label=f"any median C={C:.3g}")
        plt.xlabel(r"$\sqrt{\epsilon}\lambda^{1/2}$")
        plt.ylabel(r"$\sigma_{50}$")
        plt.title("Layer 2c: sigma50 vs half-law scale")
        plt.legend()
        plt.tight_layout()
        plt.savefig(outdir / "layer2c_sigma50_any_broad_vs_half_scale.png", dpi=180)
        plt.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="results/layer2c_full_fhn_hit_location_bins")
    parser.add_argument("--overwrite", action="store_true")

    parser.add_argument("--a", type=float, default=0.7)
    parser.add_argument("--b", type=float, default=0.8)
    parser.add_argument("--delta-Is", type=float, nargs="+", required=True)
    parser.add_argument("--epsilons", type=float, nargs="+", required=True)
    parser.add_argument("--sigma-mults", type=float, nargs="+", required=True)

    parser.add_argument("--n", type=int, default=1000)
    parser.add_argument("--dt", type=float, default=2e-4)
    parser.add_argument("--seed", type=int, default=12345)

    parser.add_argument("--q", type=float, default=0.75)

    parser.add_argument("--narrow-lo", type=float, default=0.5)
    parser.add_argument("--narrow-hi", type=float, default=2.0)
    parser.add_argument("--broad-lo", type=float, default=0.1)
    parser.add_argument("--broad-hi", type=float, default=4.0)

    parser.add_argument("--bin-edges", type=parse_edges, default=DEFAULT_BIN_EDGES,
                        help="Comma-separated bin edges, default 0,0.1,0.5,2,4,inf")

    parser.add_argument("--W-start-factor", type=float, default=8.0)
    parser.add_argument("--W-end-factor", type=float, default=-1.0)
    parser.add_argument("--T-factor", type=float, default=4.0)

    parser.add_argument("--V-abs-clip", type=float, default=1e6)
    parser.add_argument("--local-chart-V-abs", type=float, default=50.0)

    args = parser.parse_args()

    bin_edges = args.bin_edges
    bin_labels = labels_from_edges(bin_edges)

    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()) and not args.overwrite:
        raise SystemExit(f"{outdir} exists and is nonempty. Use --overwrite.")
    outdir.mkdir(parents=True, exist_ok=True)

    print("\n=== Layer 2c: full-FHN first-hit location bins ===")
    print(f"a={args.a}, b={args.b}")
    print(f"q={args.q}")
    print(f"narrow=[{args.narrow_lo},{args.narrow_hi}], broad=[{args.broad_lo},{args.broad_hi}]")
    print(f"bin_edges={bin_edges}")
    print(f"bin_labels={bin_labels}")
    print(f"n={args.n}, dt={args.dt}")
    print("sigma = sigma_mult * sqrt(epsilon) * sqrt(lambda)\n")

    batches = []
    counter = 0
    for eps in args.epsilons:
        for dI in args.delta_Is:
            lam = args.b * dI
            for mult in args.sigma_mults:
                counter += 1
                df = simulate_event_batch(
                    a=args.a,
                    b=args.b,
                    epsilon=eps,
                    delta_I=dI,
                    sigma_mult=mult,
                    n=args.n,
                    dt=args.dt,
                    seed=args.seed + 100000 * counter,
                    q=args.q,
                    narrow_lo=args.narrow_lo,
                    narrow_hi=args.narrow_hi,
                    broad_lo=args.broad_lo,
                    broad_hi=args.broad_hi,
                    bin_edges=bin_edges,
                    bin_labels=bin_labels,
                    W_start_factor=args.W_start_factor,
                    W_end_factor=args.W_end_factor,
                    T_factor=args.T_factor,
                    V_abs_clip=args.V_abs_clip,
                    local_chart_V_abs=args.local_chart_V_abs,
                )
                batches.append(df)

                any_hit = df["first_hit"].mean()
                broad = df["first_hit_in_broad"].mean()
                narrow = df["first_hit_in_narrow"].mean()
                medR = np.nanmedian(df["first_hit_R"]) if df["first_hit"].any() else np.nan
                counts = df["first_hit_bin"].value_counts(normalize=True).to_dict()

                bin_msg = ", ".join([f"{lab}={counts.get(lab,0):.3f}" for lab in bin_labels])
                print(
                    f"eps={eps:g}, dI={dI:g}, lambda={lam:g}, mult={mult:g}, "
                    f"any={any_hit:.3f}, broad={broad:.3f}, narrow={narrow:.3f}, "
                    f"medR={medR if np.isfinite(medR) else np.nan:.3g}, "
                    f"{bin_msg}, no_hit={counts.get('no_hit',0):.3f}"
                )

    events = pd.concat(batches, ignore_index=True)
    prob = make_probability_summary(events, bin_labels)
    sigma50 = make_sigma50_summary(prob, bin_labels)

    events_path = outdir / "layer2c_path_events.csv"
    prob_path = outdir / "layer2c_probability_summary.csv"
    sigma50_path = outdir / "layer2c_sigma50_summary.csv"

    events.to_csv(events_path, index=False)
    prob.to_csv(prob_path, index=False)
    sigma50.to_csv(sigma50_path, index=False)

    plot_outputs(prob, sigma50, bin_labels, outdir)

    print("\n=== Layer 2c probability summary head ===")
    print(prob.head(20).to_string(index=False))

    print("\n=== Layer 2c sigma50 summary ===")
    print(sigma50.to_string(index=False))

    print("\nSaved:")
    print(f"  {events_path}")
    print(f"  {prob_path}")
    print(f"  {sigma50_path}")
    print(f"  {outdir / 'layer2c_hit_bin_stack_by_scaled_noise.png'}")
    print(f"  {outdir / 'layer2c_hit_bin_stack_by_sigma_mult.png'}")
    print(f"  {outdir / 'layer2c_broad_vs_narrow_collapse.png'}")
    print(f"  {outdir / 'layer2c_median_hit_location.png'}")
    print(f"  {outdir / 'layer2c_sigma50_any_broad_vs_half_scale.png'}")


if __name__ == "__main__":
    main()
