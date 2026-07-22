#!/usr/bin/env python3
"""
layer2b_full_fhn_event_classifier.py

Layer 2b:
  Full autonomous stochastic FHN event classifier.

Why this exists:
  The original Layer 2 only counted a hit if the full FHN path reached

      V_blow >= q sqrt(W_blow)

  inside a narrow normal-form window

      W_blow in [window_lo W*, window_hi W*].

  In full FHN, paths may instead:
    - hit earlier, above the window;
    - hit inside the window;
    - hit later, below the window but before W=0;
    - leave the local chart without hitting;
    - curl/re-enter the local region.

  This script classifies the FIRST local hit by where it occurs in W/W*,
  and separately tracks local window entries / returns.

Model convention:
    dv = (v - v^3/3 - w + I)/eps dt + (sigma/sqrt(eps)) dB_t
    dw = (v + a - b w) dt

This matches the convention eta = sigma/sqrt(eps) in blow-up coordinates.

Correct left-fold canard point:
    I_fold_left = 2/3 + (a - 1)/b

At I = I_fold_left + delta_I:
    lambda = b * delta_I

Blow-up coordinates:
    V = (v + 1)/eps^(1/3)
    W = (w - w_fold_left)/eps^(2/3)
    W* = lambda^(2/3)

Hit target:
    V >= q sqrt(W), while W > 0.

Classify first hit:
    pre_window:     W/W* > window_hi
    window:         window_lo <= W/W* <= window_hi
    post_window:    0 < W/W* < window_lo
    no_hit:         no target hit before leaving / stopping

Also tracks:
    entry_count: how many times path enters the local window
    exit_count:  how many times path exits the local window after entry
    first_entry_t
    first_exit_t
    min_W_ratio, max_W_ratio seen

Recommended run:

python3 layer2b_full_fhn_event_classifier.py \
  --overwrite \
  --a 0.7 \
  --b 0.8 \
  --delta-Is 0.01 0.02 0.04 \
  --epsilons 0.01 0.02 0.04 0.08 \
  --sigma-mults 0 0.5 1.0 1.5 2.0 2.5 3.0 4.0 5.0 6.5 8.0 \
  --n 1000 \
  --outdir results/layer2b_full_fhn_event_classifier

If probabilities still saturate below 0.5, run wider:

python3 layer2b_full_fhn_event_classifier.py \
  --overwrite \
  --a 0.7 \
  --b 0.8 \
  --delta-Is 0.01 0.02 0.04 \
  --epsilons 0.01 0.02 0.04 0.08 \
  --sigma-mults 0 1 2 3 4 5 6.5 8 10 12 15 \
  --n 1000 \
  --outdir results/layer2b_full_fhn_event_classifier_wide

Outputs:
  layer2b_path_events.csv
      one row per simulated path

  layer2b_probability_summary.csv
      probabilities by epsilon, delta_I, lambda, sigma_mult

  layer2b_sigma50_summary.csv
      estimated sigma50 for:
        any_hit
        pre_window
        window
        post_window
        no_hit is not fitted

  layer2b_probabilities_by_scaled_noise.png
  layer2b_event_stack_by_scaled_noise.png
  layer2b_anyhit_collapse_half.png
  layer2b_windowhit_collapse_half.png
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


EVENT_ORDER = ["pre_window", "window", "post_window", "no_hit"]


def cubic_roots_for_w_I(w: float, I: float) -> np.ndarray:
    # v - v^3/3 + I - w = 0
    # equivalent: v^3 - 3v + 3(w-I)=0
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


def classify_hit(W_ratio: float, window_lo: float, window_hi: float) -> str:
    if not np.isfinite(W_ratio) or W_ratio <= 0:
        return "no_hit"
    if W_ratio > window_hi:
        return "pre_window"
    if W_ratio >= window_lo:
        return "window"
    return "post_window"


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
    window_lo: float,
    window_hi: float,
    W_start_factor: float,
    W_end_factor: float,
    T_factor: float,
    V_abs_clip: float,
    local_chart_V_abs: float,
) -> pd.DataFrame:
    """
    Simulate n independent full-FHN paths and classify first local hit.
    """
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
    first_hit_W_ratio = np.full(n, np.nan)
    first_hit_class = np.array(["no_hit"] * n, dtype=object)

    # Window crossing / recurl-ish diagnostics.
    in_window_prev = np.zeros(n, dtype=bool)
    entry_count = np.zeros(n, dtype=int)
    exit_count = np.zeros(n, dtype=int)
    first_entry_t = np.full(n, np.nan)
    first_exit_t = np.full(n, np.nan)
    min_W_ratio = np.full(n, np.inf)
    max_W_ratio = np.full(n, -np.inf)

    # Conservative simulation length. Full FHN near a fold can slow down, so we use T_factor.
    # Estimate slow time needed to pass through local chart.
    slow_speed0 = abs(v_start + a - b * w_start)
    T_base = abs(epsilon ** (2.0 / 3.0) * (W_start_blow - W_end_blow)) / max(slow_speed0, 1e-8)
    T_max = T_factor * T_base
    n_steps = int(np.ceil(T_max / dt))

    sqrt_dt = math.sqrt(dt)

    for step in range(n_steps):
        t = step * dt

        W_blow = (w - w_fold_left) / (epsilon ** (2.0 / 3.0))
        V_blow = (v + 1.0) / (epsilon ** (1.0 / 3.0))
        W_ratio = W_blow / W_star

        finite = np.isfinite(W_ratio)
        min_W_ratio[finite] = np.minimum(min_W_ratio[finite], W_ratio[finite])
        max_W_ratio[finite] = np.maximum(max_W_ratio[finite], W_ratio[finite])

        in_window = active & (W_ratio >= window_lo) & (W_ratio <= window_hi)
        entering = in_window & (~in_window_prev)
        exiting = (~in_window) & in_window_prev

        entry_count[entering] += 1
        new_first_entry = entering & ~np.isfinite(first_entry_t)
        first_entry_t[new_first_entry] = t

        exit_count[exiting] += 1
        new_first_exit = exiting & ~np.isfinite(first_exit_t)
        first_exit_t[new_first_exit] = t

        in_window_prev = in_window.copy()

        # Hit target anywhere W > 0, not just inside the narrow window.
        good_W = active & (W_blow > 0) & np.isfinite(V_blow) & np.isfinite(W_blow)
        threshold = np.full(n, np.inf)
        threshold[good_W] = q * np.sqrt(W_blow[good_W])

        hit_now = good_W & (~first_hit) & (V_blow >= threshold)

        if np.any(hit_now):
            first_hit[hit_now] = True
            first_hit_t[hit_now] = t
            first_hit_V[hit_now] = V_blow[hit_now]
            first_hit_W[hit_now] = W_blow[hit_now]
            first_hit_W_ratio[hit_now] = W_ratio[hit_now]
            for j in np.where(hit_now)[0]:
                first_hit_class[j] = classify_hit(float(W_ratio[j]), window_lo, window_hi)
            active[hit_now] = False

        # Stop tracking paths that have clearly left the local chart without a hit.
        # This prevents huge-noise global wander from polluting a local event classifier.
        too_far_V = np.abs(V_blow) > local_chart_V_abs
        below_chart = W_ratio < W_end_factor
        above_chart = W_ratio > max(W_start_factor * 1.5, window_hi * 4.0)
        inactive_now = active & (too_far_V | below_chart | above_chart)
        active[inactive_now] = False

        if not np.any(active):
            break

        # Euler-Maruyama update for still active paths only.
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

    min_W_ratio[~np.isfinite(min_W_ratio)] = np.nan
    max_W_ratio[~np.isfinite(max_W_ratio)] = np.nan

    df = pd.DataFrame({
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
        "window_lo": window_lo,
        "window_hi": window_hi,
        "W_star": W_star,
        "path": np.arange(n),
        "first_hit": first_hit,
        "first_hit_class": first_hit_class,
        "first_hit_t": first_hit_t,
        "first_hit_V": first_hit_V,
        "first_hit_W": first_hit_W,
        "first_hit_W_ratio": first_hit_W_ratio,
        "entry_count": entry_count,
        "exit_count": exit_count,
        "first_entry_t": first_entry_t,
        "first_exit_t": first_exit_t,
        "min_W_ratio": min_W_ratio,
        "max_W_ratio": max_W_ratio,
        "v_start": v_start,
        "w_start": w_start,
        "T_max": T_max,
        "n_steps": n_steps,
    })
    return df


def make_probability_summary(events: pd.DataFrame) -> pd.DataFrame:
    group_cols = ["epsilon", "delta_I", "lambda", "sigma_mult", "sigma"]

    rows = []
    for keys, g in events.groupby(group_cols):
        row = dict(zip(group_cols, keys))
        n = len(g)
        row["n"] = n

        for ev in EVENT_ORDER:
            row[f"p_{ev}"] = float(np.mean(g["first_hit_class"] == ev))

        row["p_any_hit"] = float(np.mean(g["first_hit"]))
        row["p_entered_window"] = float(np.mean(g["entry_count"] > 0))
        row["mean_entry_count"] = float(np.mean(g["entry_count"]))
        row["p_multiple_entries"] = float(np.mean(g["entry_count"] >= 2))
        row["median_first_hit_W_ratio"] = float(np.nanmedian(g["first_hit_W_ratio"])) if np.any(g["first_hit"]) else np.nan
        rows.append(row)

    summary = pd.DataFrame(rows).sort_values(group_cols)
    summary["scale_half"] = summary["sigma"] / (np.sqrt(summary["epsilon"]) * np.sqrt(summary["lambda"]))
    summary["scale_twothirds"] = summary["sigma"] / (np.sqrt(summary["epsilon"]) * (summary["lambda"] ** (2.0 / 3.0)))
    summary["scale_eps34"] = summary["sigma"] / (summary["epsilon"] ** 0.75)
    return summary


def make_sigma50_summary(prob: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for keys, g in prob.groupby(["epsilon", "delta_I", "lambda"]):
        eps, dI, lam = keys
        g = g.sort_values("sigma")
        row = {"epsilon": eps, "delta_I": dI, "lambda": lam}
        for pcol in ["p_any_hit", "p_pre_window", "p_window", "p_post_window", "p_entered_window"]:
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


def plot_outputs(prob: pd.DataFrame, sigma50: pd.DataFrame, outdir: Path) -> None:
    # Any-hit collapse.
    plt.figure(figsize=(11, 7))
    for keys, g in prob.groupby(["epsilon", "delta_I"]):
        eps, dI = keys
        g = g.sort_values("scale_half")
        plt.plot(g["scale_half"], g["p_any_hit"], marker="o", label=f"eps={eps:g}, dI={dI:g}")
    plt.axhline(0.5, linestyle="--", color="black")
    plt.xlabel(r"$\sigma/(\sqrt{\epsilon}\lambda^{1/2})$")
    plt.ylabel("P(any target hit for W>0)")
    plt.title("Layer 2b: any-hit collapse using half-law scale")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outdir / "layer2b_anyhit_collapse_half.png", dpi=180)
    plt.close()

    # Window-hit collapse.
    plt.figure(figsize=(11, 7))
    for keys, g in prob.groupby(["epsilon", "delta_I"]):
        eps, dI = keys
        g = g.sort_values("scale_half")
        plt.plot(g["scale_half"], g["p_window"], marker="o", label=f"eps={eps:g}, dI={dI:g}")
    plt.axhline(0.5, linestyle="--", color="black")
    plt.xlabel(r"$\sigma/(\sqrt{\epsilon}\lambda^{1/2})$")
    plt.ylabel("P(first hit in window)")
    plt.title("Layer 2b: window-hit collapse using half-law scale")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outdir / "layer2b_windowhit_collapse_half.png", dpi=180)
    plt.close()

    # Probability components against scaled noise, all curves pooled.
    plt.figure(figsize=(11, 7))
    for pcol, label in [
        ("p_pre_window", "pre-window"),
        ("p_window", "window"),
        ("p_post_window", "post-window"),
        ("p_any_hit", "any hit"),
        ("p_entered_window", "entered window"),
    ]:
        pooled = prob.groupby("scale_half")[pcol].mean().reset_index().sort_values("scale_half")
        # groupby exact float is okay here because sigma_mult grid gives repeated values.
        plt.plot(pooled["scale_half"], pooled[pcol], marker="o", label=label)
    plt.xlabel(r"$\sigma/(\sqrt{\epsilon}\lambda^{1/2})$")
    plt.ylabel("pooled probability")
    plt.title("Layer 2b: event probabilities by half-law scaled noise")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "layer2b_probabilities_by_scaled_noise.png", dpi=180)
    plt.close()

    # Stacked event proportions by sigma_mult, pooled over eps/dI.
    pooled = prob.groupby("sigma_mult")[["p_pre_window", "p_window", "p_post_window", "p_no_hit"]].mean().reset_index().sort_values("sigma_mult")
    x = np.arange(len(pooled))
    bottom = np.zeros(len(pooled))

    plt.figure(figsize=(12, 7))
    for col, label in [
        ("p_pre_window", "pre-window"),
        ("p_window", "window"),
        ("p_post_window", "post-window"),
        ("p_no_hit", "no hit"),
    ]:
        plt.bar(x, pooled[col], bottom=bottom, label=label)
        bottom += pooled[col].to_numpy()
    plt.xticks(x, [f"{m:g}" for m in pooled["sigma_mult"]])
    plt.xlabel("sigma multiplier in sigma = mult sqrt(eps) sqrt(lambda)")
    plt.ylabel("pooled event fraction")
    plt.title("Layer 2b: pooled first-hit event composition")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "layer2b_event_stack_by_scaled_noise.png", dpi=180)
    plt.close()

    # sigma50 vs half scale if available.
    clean = sigma50.dropna(subset=["sigma50_any_hit"])
    if len(clean) >= 2:
        plt.figure(figsize=(9, 7))
        plt.loglog(clean["pred_half_scale"], clean["sigma50_any_hit"], "o", label="measured any-hit sigma50")
        C = np.nanmedian(clean["sigma50_any_hit"] / clean["pred_half_scale"])
        xgrid = np.linspace(clean["pred_half_scale"].min(), clean["pred_half_scale"].max(), 100)
        plt.loglog(xgrid, C * xgrid, "--", label=f"median C={C:.3g}")
        plt.xlabel(r"$\sqrt{\epsilon}\lambda^{1/2}$")
        plt.ylabel(r"$\sigma_{50}$")
        plt.title("Layer 2b: any-hit sigma50 vs half-law prediction scale")
        plt.legend()
        plt.tight_layout()
        plt.savefig(outdir / "layer2b_sigma50_anyhit_vs_half_scale.png", dpi=180)
        plt.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="results/layer2b_full_fhn_event_classifier")
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
    parser.add_argument("--window-lo", type=float, default=0.5)
    parser.add_argument("--window-hi", type=float, default=2.0)

    parser.add_argument("--W-start-factor", type=float, default=8.0)
    parser.add_argument("--W-end-factor", type=float, default=-1.0)
    parser.add_argument("--T-factor", type=float, default=4.0)

    parser.add_argument("--V-abs-clip", type=float, default=1e6)
    parser.add_argument("--local-chart-V-abs", type=float, default=50.0)

    args = parser.parse_args()

    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()) and not args.overwrite:
        raise SystemExit(f"{outdir} exists and is nonempty. Use --overwrite.")
    outdir.mkdir(parents=True, exist_ok=True)

    print("\n=== Layer 2b: full-FHN event classifier ===")
    print(f"a={args.a}, b={args.b}")
    print(f"q={args.q}, window=[{args.window_lo}W*, {args.window_hi}W*]")
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
                    window_lo=args.window_lo,
                    window_hi=args.window_hi,
                    W_start_factor=args.W_start_factor,
                    W_end_factor=args.W_end_factor,
                    T_factor=args.T_factor,
                    V_abs_clip=args.V_abs_clip,
                    local_chart_V_abs=args.local_chart_V_abs,
                )
                batches.append(df)

                counts = df["first_hit_class"].value_counts(normalize=True).to_dict()
                msg = (
                    f"eps={eps:g}, dI={dI:g}, lambda={lam:g}, mult={mult:g}, "
                    f"any={df['first_hit'].mean():.3f}, "
                    f"pre={counts.get('pre_window', 0):.3f}, "
                    f"win={counts.get('window', 0):.3f}, "
                    f"post={counts.get('post_window', 0):.3f}, "
                    f"entries>0={(df['entry_count'] > 0).mean():.3f}, "
                    f"multi={(df['entry_count'] >= 2).mean():.3f}"
                )
                print(msg)

    events = pd.concat(batches, ignore_index=True)
    prob = make_probability_summary(events)
    sigma50 = make_sigma50_summary(prob)

    events_path = outdir / "layer2b_path_events.csv"
    prob_path = outdir / "layer2b_probability_summary.csv"
    sigma50_path = outdir / "layer2b_sigma50_summary.csv"

    events.to_csv(events_path, index=False)
    prob.to_csv(prob_path, index=False)
    sigma50.to_csv(sigma50_path, index=False)

    plot_outputs(prob, sigma50, outdir)

    print("\n=== Layer 2b probability summary head ===")
    print(prob.head(20).to_string(index=False))

    print("\n=== Layer 2b sigma50 summary ===")
    print(sigma50.to_string(index=False))

    print("\nSaved:")
    print(f"  {events_path}")
    print(f"  {prob_path}")
    print(f"  {sigma50_path}")
    print(f"  {outdir / 'layer2b_anyhit_collapse_half.png'}")
    print(f"  {outdir / 'layer2b_windowhit_collapse_half.png'}")
    print(f"  {outdir / 'layer2b_probabilities_by_scaled_noise.png'}")
    print(f"  {outdir / 'layer2b_event_stack_by_scaled_noise.png'}")


if __name__ == "__main__":
    main()
