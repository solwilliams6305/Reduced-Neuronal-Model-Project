#!/usr/bin/env python3
"""
layer2_full_fhn_collapse.py

Layer 2 script:
  Full autonomous stochastic FHN validation/collapse test.

Purpose:
  Check whether full FHN escape probabilities collapse better under

      sigma / (sqrt(epsilon) * lambda^(1/2))

  than under

      sigma / (sqrt(epsilon) * lambda^(2/3))
      sigma / epsilon^(3/4)

This is not meant to be a full theory of autonomous FHN. It is a bridge:
  local blow-up law -> full model evidence.

Model convention used here:
    dV = (v - v^3/3 - w + I)/epsilon dt + (sigma/sqrt(epsilon)) dB_t
    dW = (v + a - b w) dt

The sigma convention is chosen so that the blow-up noise satisfies eta=sigma/sqrt(epsilon),
matching the normal-form work. If your existing kernel uses a different noise convention,
adjust the noise line in simulate_full_fhn_probability().

Event:
  Start near the left attracting branch at a chosen offset above the left fold.
  Count "local early hit" if, while in the mapped blow-up window
      W_blow = (w - w_fold_left) / epsilon^(2/3)
      W_blow in [window_lo * lambda^(2/3), window_hi * lambda^(2/3)],
  the path reaches
      V_blow >= q sqrt(W_blow),
  where
      V_blow = (v + 1) / epsilon^(1/3).

This intentionally mirrors the normal-form observable.

Recommended run:

python3 layer2_full_fhn_collapse.py \
  --overwrite \
  --a 0.7 \
  --b 0.8 \
  --delta-Is 0.01 0.02 0.04 \
  --epsilons 0.01 0.02 0.04 0.08 \
  --sigma-mults 0 0.35 0.5 0.7 0.9 1.1 1.3 1.6 2.0 2.5 3.0 \
  --n 1000 \
  --outdir results/layer2_full_fhn_collapse

Outputs:
  layer2_full_fhn_results.csv
  layer2_sigma50_summary.csv
  layer2_collapse_half.png
  layer2_collapse_twothirds.png
  layer2_collapse_eps34.png
  layer2_sigma50_vs_prediction.png
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def cubic_roots_for_w_I(w: float, I: float) -> np.ndarray:
    # v - v^3/3 + I - w = 0
    # equivalently v^3 - 3v + 3(w-I) = 0
    coeffs = [1.0, 0.0, -3.0, 3.0 * (w - I)]
    roots = np.roots(coeffs)
    real = np.real(roots[np.abs(np.imag(roots)) < 1e-8])
    return np.sort(real)


def left_attracting_v(w: float, I: float) -> float:
    roots = cubic_roots_for_w_I(w, I)
    if len(roots) == 0:
        raise ValueError("No real root found, impossible for cubic.")
    return float(roots[0])


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


def simulate_full_fhn_probability(
    *,
    a: float,
    b: float,
    epsilon: float,
    delta_I: float,
    sigma: float,
    n: int,
    dt: float,
    seed: int,
    q: float,
    window_lo: float,
    window_hi: float,
    W_start_factor: float,
    W_end_factor: float,
    V_abs_clip: float,
) -> dict:
    """
    Full FHN local-window hit probability.

    I is set to I_fold_left + delta_I.
    lambda = b * delta_I, matching your convention from the README.
    """
    rng = np.random.default_rng(seed)

    I_fold_left = 2.0 / 3.0 + (a - 1.0) / b  # from equilibrium at v=-1, w=I-2/3 and w=(v+a)/b
    I = I_fold_left + delta_I
    w_fold_left = I - 2.0 / 3.0
    lam = b * delta_I

    W_star = lam ** (2.0 / 3.0)
    W_start_blow = W_start_factor * W_star
    W_end_blow = W_end_factor * W_star

    w_start = w_fold_left + (epsilon ** (2.0 / 3.0)) * W_start_blow
    w_end = w_fold_left + (epsilon ** (2.0 / 3.0)) * W_end_blow

    v_start = left_attracting_v(w_start, I)

    v = np.full(n, v_start, dtype=float)
    w = np.full(n, w_start, dtype=float)
    active = np.ones(n, dtype=bool)
    hit = np.zeros(n, dtype=bool)

    # Integrate long enough to pass from W_start to W_end under slow drift.
    # Use conservative time cap. If all active paths pass below the local window, stop.
    T_max = 4.0 * abs(w_start - w_end) / max(abs(v_start + a - b * w_start), 1e-8)
    n_steps = int(np.ceil(T_max / dt))
    sqrt_dt = math.sqrt(dt)

    W_lo = window_lo * W_star
    W_hi = window_hi * W_star

    for _ in range(n_steps):
        if not np.any(active):
            break

        idx = active
        vv = v[idx]
        ww = w[idx]

        drift_v = (vv - vv**3 / 3.0 - ww + I) / epsilon
        drift_w = vv + a - b * ww

        # Noise convention matching eta=sigma/sqrt(epsilon) in blow-up coordinates.
        vv = vv + drift_v * dt + (sigma / math.sqrt(epsilon)) * sqrt_dt * rng.normal(size=np.count_nonzero(idx))
        ww = ww + drift_w * dt

        vv = np.clip(vv, -V_abs_clip, V_abs_clip)

        v[idx] = vv
        w[idx] = ww

        W_blow = (w - w_fold_left) / (epsilon ** (2.0 / 3.0))
        V_blow = (v + 1.0) / (epsilon ** (1.0 / 3.0))

        in_window = active & (W_blow > 0) & (W_blow >= W_lo) & (W_blow <= W_hi)
        threshold = np.zeros(n)
        goodW = W_blow > 0
        threshold[goodW] = q * np.sqrt(W_blow[goodW])
        hit_now = in_window & (V_blow >= threshold)
        hit |= hit_now
        active[hit_now] = False

        # Once below local window/fold region, no future hit for this local event.
        below = active & (W_blow < min(W_lo, 0.0))
        active[below] = False

        # Also stop paths that exploded far to the right/left.
        escaped_chart = active & (np.abs(V_blow) > 50)
        active[escaped_chart] = False

    return {
        "a": a,
        "b": b,
        "I_fold_left": I_fold_left,
        "I": I,
        "delta_I": delta_I,
        "lambda": lam,
        "epsilon": epsilon,
        "sigma": sigma,
        "q": q,
        "window_lo": window_lo,
        "window_hi": window_hi,
        "hit_prob": float(np.mean(hit)),
        "n": n,
        "W_star": W_star,
        "W_lo": W_lo,
        "W_hi": W_hi,
        "v_start": v_start,
        "w_start": w_start,
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", default="results/layer2_full_fhn_collapse")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--a", type=float, default=0.7)
    p.add_argument("--b", type=float, default=0.8)
    p.add_argument("--delta-Is", type=float, nargs="+", required=True)
    p.add_argument("--epsilons", type=float, nargs="+", required=True)

    # sigma = multiplier * sqrt(eps) * sqrt(lambda)
    p.add_argument("--sigma-mults", type=float, nargs="+", required=True)

    p.add_argument("--n", type=int, default=1000)
    p.add_argument("--dt", type=float, default=2e-4)
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--q", type=float, default=0.75)
    p.add_argument("--window-lo", type=float, default=0.5)
    p.add_argument("--window-hi", type=float, default=2.0)
    p.add_argument("--W-start-factor", type=float, default=8.0)
    p.add_argument("--W-end-factor", type=float, default=-1.0)
    p.add_argument("--V-abs-clip", type=float, default=1e6)
    args = p.parse_args()

    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()) and not args.overwrite:
        raise SystemExit(f"{outdir} exists and is nonempty. Use --overwrite.")
    outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    counter = 0

    print("\n=== Layer 2: full FHN collapse validation ===")
    print("sigma grid uses sigma = multiplier * sqrt(epsilon) * lambda^(1/2)\n")

    for eps in args.epsilons:
        for dI in args.delta_Is:
            lam = args.b * dI
            for m in args.sigma_mults:
                sigma = m * math.sqrt(eps) * math.sqrt(lam)
                counter += 1
                row = simulate_full_fhn_probability(
                    a=args.a,
                    b=args.b,
                    epsilon=eps,
                    delta_I=dI,
                    sigma=sigma,
                    n=args.n,
                    dt=args.dt,
                    seed=args.seed + 100000 * counter,
                    q=args.q,
                    window_lo=args.window_lo,
                    window_hi=args.window_hi,
                    W_start_factor=args.W_start_factor,
                    W_end_factor=args.W_end_factor,
                    V_abs_clip=args.V_abs_clip,
                )
                rows.append(row)
                print(f"eps={eps:g}, dI={dI:g}, lambda={lam:g}, sigma={sigma:.5g}, mult={m:g}, P={row['hit_prob']:.3f}")

    results = pd.DataFrame(rows)
    results["scale_half"] = results["sigma"] / (np.sqrt(results["epsilon"]) * np.sqrt(results["lambda"]))
    results["scale_twothirds"] = results["sigma"] / (np.sqrt(results["epsilon"]) * (results["lambda"] ** (2/3)))
    results["scale_eps34"] = results["sigma"] / (results["epsilon"] ** 0.75)

    sigma50_rows = []
    for keys, g in results.groupby(["epsilon", "delta_I", "lambda"]):
        eps, dI, lam = keys
        g = g.sort_values("sigma")
        sigma50 = interp_x_at_y(g["sigma"].to_numpy(), g["hit_prob"].to_numpy(), 0.5)
        sigma50_rows.append({
            "epsilon": eps,
            "delta_I": dI,
            "lambda": lam,
            "sigma50": sigma50,
            "pred_half_scale": math.sqrt(eps) * math.sqrt(lam),
            "pred_twothirds_scale": math.sqrt(eps) * (lam ** (2/3)),
            "pred_eps34_scale": eps ** 0.75,
            "sigma50_over_half": sigma50 / (math.sqrt(eps) * math.sqrt(lam)) if np.isfinite(sigma50) else np.nan,
            "sigma50_over_twothirds": sigma50 / (math.sqrt(eps) * (lam ** (2/3))) if np.isfinite(sigma50) else np.nan,
            "sigma50_over_eps34": sigma50 / (eps ** 0.75) if np.isfinite(sigma50) else np.nan,
        })
    summary = pd.DataFrame(sigma50_rows)

    results.to_csv(outdir / "layer2_full_fhn_results.csv", index=False)
    summary.to_csv(outdir / "layer2_sigma50_summary.csv", index=False)

    # Collapse plots.
    def plot_collapse(xcol: str, fname: str, title: str):
        plt.figure(figsize=(10, 7))
        for keys, g in results.groupby(["epsilon", "delta_I"]):
            eps, dI = keys
            g = g.sort_values(xcol)
            plt.plot(g[xcol], g["hit_prob"], marker="o", label=f"eps={eps:g}, dI={dI:g}")
        plt.axhline(0.5, linestyle="--", color="black")
        plt.xlabel(xcol)
        plt.ylabel("P(local window hit)")
        plt.title(title)
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(outdir / fname, dpi=180)
        plt.close()

    plot_collapse("scale_half", "layer2_collapse_half.png", r"Full FHN collapse using $\sigma/(\sqrt{\epsilon}\lambda^{1/2})$")
    plot_collapse("scale_twothirds", "layer2_collapse_twothirds.png", r"Full FHN collapse using $\sigma/(\sqrt{\epsilon}\lambda^{2/3})$")
    plot_collapse("scale_eps34", "layer2_collapse_eps34.png", r"Full FHN collapse using $\sigma/\epsilon^{3/4}$")

    # sigma50 predicted-scale comparison.
    plt.figure(figsize=(10, 7))
    clean = summary.dropna(subset=["sigma50"])
    plt.loglog(clean["pred_half_scale"], clean["sigma50"], "o", label=r"$\sqrt{\epsilon}\lambda^{1/2}$")
    if len(clean) >= 2:
        x = np.linspace(clean["pred_half_scale"].min(), clean["pred_half_scale"].max(), 100)
        C = np.nanmedian(clean["sigma50"] / clean["pred_half_scale"])
        plt.loglog(x, C * x, "--", label=f"median C={C:.3g}")
    plt.xlabel(r"$\sqrt{\epsilon}\lambda^{1/2}$")
    plt.ylabel(r"$\sigma_{50}$")
    plt.title("Full FHN sigma50 against blow-up prediction scale")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "layer2_sigma50_vs_prediction.png", dpi=180)
    plt.close()

    print("\n=== Layer 2 sigma50 summary ===")
    print(summary.to_string(index=False))
    print(f"\nSaved outputs to {outdir}")


if __name__ == "__main__":
    main()
