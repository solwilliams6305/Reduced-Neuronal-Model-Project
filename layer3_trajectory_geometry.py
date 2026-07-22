#!/usr/bin/env python3
"""
layer3_trajectory_geometry.py

Layer 3 script:
  Representative trajectories showing how full FHN paths enter/leave the local canard window.

It plots both:
  1. physical (v,w) phase plane with cubic nullcline and local window;
  2. blown-up coordinates (V,W), with attracting/repelling branches
       V = ±sqrt(W)
     and target
       V = q sqrt(W).

Recommended run:

python3 layer3_trajectory_geometry.py \
  --overwrite \
  --a 0.7 \
  --b 0.8 \
  --epsilon 0.02 \
  --delta-I 0.02 \
  --sigma-mults 0 1.0 2.2 3.0 \
  --n-paths 8 \
  --outdir results/layer3_trajectory_geometry

sigma is defined as:
  sigma = sigma_mult * sqrt(epsilon) * lambda^(1/2)

So sigma_mult around 2.2 is approximately the normal-form q=0.75 eta50 prefactor.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def cubic_roots_for_w_I(w: float, I: float) -> np.ndarray:
    coeffs = [1.0, 0.0, -3.0, 3.0 * (w - I)]
    roots = np.roots(coeffs)
    real = np.real(roots[np.abs(np.imag(roots)) < 1e-8])
    return np.sort(real)


def left_attracting_v(w: float, I: float) -> float:
    roots = cubic_roots_for_w_I(w, I)
    return float(roots[0])


def simulate_paths(
    *,
    a: float,
    b: float,
    epsilon: float,
    delta_I: float,
    sigma_mult: float,
    n_paths: int,
    dt: float,
    seed: int,
    q: float,
    window_lo: float,
    window_hi: float,
    W_start_factor: float,
    W_end_factor: float,
):
    rng = np.random.default_rng(seed)

    I_fold_left = 2.0 / 3.0 + (a - 1.0) / b
    I = I_fold_left + delta_I
    w_fold = I - 2.0 / 3.0
    lam = b * delta_I
    W_star = lam ** (2.0 / 3.0)
    sigma = sigma_mult * math.sqrt(epsilon) * math.sqrt(lam)

    W_start_blow = W_start_factor * W_star
    W_end_blow = W_end_factor * W_star
    w_start = w_fold + epsilon ** (2.0 / 3.0) * W_start_blow
    w_end = w_fold + epsilon ** (2.0 / 3.0) * W_end_blow
    v_start = left_attracting_v(w_start, I)

    T_max = 4.0 * abs(w_start - w_end) / max(abs(v_start + a - b * w_start), 1e-8)
    n_steps = int(np.ceil(T_max / dt))
    stride = max(1, n_steps // 2500)

    v = np.full(n_paths, v_start)
    w = np.full(n_paths, w_start)
    hit = np.zeros(n_paths, dtype=bool)
    hit_t = np.full(n_paths, np.nan)
    hit_V = np.full(n_paths, np.nan)
    hit_W = np.full(n_paths, np.nan)

    rows = []

    for step in range(n_steps):
        t = step * dt

        W_blow = (w - w_fold) / (epsilon ** (2.0 / 3.0))
        V_blow = (v + 1.0) / (epsilon ** (1.0 / 3.0))

        if step % stride == 0:
            for j in range(n_paths):
                rows.append({
                    "path": j,
                    "t": t,
                    "v": v[j],
                    "w": w[j],
                    "V_blow": V_blow[j],
                    "W_blow": W_blow[j],
                    "hit": bool(hit[j]),
                    "sigma_mult": sigma_mult,
                    "sigma": sigma,
                    "lambda": lam,
                    "epsilon": epsilon,
                    "I": I,
                    "w_fold": w_fold,
                })

        in_window = (W_blow > 0) & (W_blow >= window_lo * W_star) & (W_blow <= window_hi * W_star)
        threshold = np.zeros_like(W_blow)
        good = W_blow > 0
        threshold[good] = q * np.sqrt(W_blow[good])
        hit_now = (~hit) & in_window & (V_blow >= threshold)
        hit[hit_now] = True
        hit_t[hit_now] = t
        hit_V[hit_now] = V_blow[hit_now]
        hit_W[hit_now] = W_blow[hit_now]

        # Euler-Maruyama
        drift_v = (v - v**3 / 3.0 - w + I) / epsilon
        drift_w = v + a - b * w
        v = v + drift_v * dt + (sigma / math.sqrt(epsilon)) * math.sqrt(dt) * rng.normal(size=n_paths)
        w = w + drift_w * dt

        # stop once all paths are below the local chart
        if np.all(W_blow < min(0.0, window_lo * W_star)):
            break

    traj = pd.DataFrame(rows)
    hits = pd.DataFrame({
        "path": np.arange(n_paths),
        "hit": hit,
        "hit_t": hit_t,
        "hit_V": hit_V,
        "hit_W": hit_W,
        "sigma_mult": sigma_mult,
        "sigma": sigma,
        "lambda": lam,
        "epsilon": epsilon,
        "I": I,
        "w_fold": w_fold,
        "W_star": W_star,
        "window_lo_W": window_lo * W_star,
        "window_hi_W": window_hi * W_star,
        "q": q,
    })
    return traj, hits


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", default="results/layer3_trajectory_geometry")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--a", type=float, default=0.7)
    p.add_argument("--b", type=float, default=0.8)
    p.add_argument("--epsilon", type=float, default=0.02)
    p.add_argument("--delta-I", type=float, default=0.02)
    p.add_argument("--sigma-mults", type=float, nargs="+", default=[0.0, 1.0, 2.2, 3.0])
    p.add_argument("--n-paths", type=int, default=8)
    p.add_argument("--dt", type=float, default=2e-4)
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--q", type=float, default=0.75)
    p.add_argument("--window-lo", type=float, default=0.5)
    p.add_argument("--window-hi", type=float, default=2.0)
    p.add_argument("--W-start-factor", type=float, default=8.0)
    p.add_argument("--W-end-factor", type=float, default=-1.0)
    args = p.parse_args()

    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()) and not args.overwrite:
        raise SystemExit(f"{outdir} exists and is nonempty. Use --overwrite.")
    outdir.mkdir(parents=True, exist_ok=True)

    all_traj = []
    all_hits = []
    for i, m in enumerate(args.sigma_mults):
        traj, hits = simulate_paths(
            a=args.a,
            b=args.b,
            epsilon=args.epsilon,
            delta_I=args.delta_I,
            sigma_mult=m,
            n_paths=args.n_paths,
            dt=args.dt,
            seed=args.seed + 1000 * i,
            q=args.q,
            window_lo=args.window_lo,
            window_hi=args.window_hi,
            W_start_factor=args.W_start_factor,
            W_end_factor=args.W_end_factor,
        )
        all_traj.append(traj)
        all_hits.append(hits)
        print(f"sigma_mult={m:g}, hit fraction={hits['hit'].mean():.3f}")

    traj_df = pd.concat(all_traj, ignore_index=True)
    hits_df = pd.concat(all_hits, ignore_index=True)

    traj_df.to_csv(outdir / "layer3_trajectories.csv", index=False)
    hits_df.to_csv(outdir / "layer3_hit_summary.csv", index=False)

    # Common constants.
    first = hits_df.iloc[0]
    I = first["I"]
    w_fold = first["w_fold"]
    W_star = first["W_star"]
    eps = first["epsilon"]
    q = first["q"]
    Wlo = first["window_lo_W"]
    Whi = first["window_hi_W"]

    # Physical phase plane plots, one per sigma multiplier.
    for m, g in traj_df.groupby("sigma_mult"):
        plt.figure(figsize=(9, 7))

        # Cubic nullcline.
        v_grid = np.linspace(-1.8, 0.4, 500)
        w_cubic = v_grid - v_grid**3 / 3.0 + I
        plt.plot(v_grid, w_cubic, color="black", linewidth=2, label=r"$w=v-v^3/3+I$")

        # local window in physical w coordinates
        wlo_phys = w_fold + eps ** (2/3) * Wlo
        whi_phys = w_fold + eps ** (2/3) * Whi
        plt.axhspan(wlo_phys, whi_phys, alpha=0.15, label="local canard window")
        plt.axhline(w_fold, linestyle="--", color="black", linewidth=1, label="left fold w")

        for path, gp in g.groupby("path"):
            plt.plot(gp["v"], gp["w"], linewidth=1.2, alpha=0.85)
        hs = hits_df[(hits_df["sigma_mult"] == m) & (hits_df["hit"])]
        if len(hs):
            vhit = -1 + eps ** (1/3) * hs["hit_V"]
            whit = w_fold + eps ** (2/3) * hs["hit_W"]
            plt.scatter(vhit, whit, marker="x", s=80, label="first local hit")

        plt.xlabel("v")
        plt.ylabel("w")
        plt.title(f"Full FHN trajectories, sigma_mult={m:g}")
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(outdir / f"layer3_phase_sigma_mult_{m:g}.png", dpi=180)
        plt.close()

    # Blow-up coordinate plots, one per sigma multiplier.
    for m, g in traj_df.groupby("sigma_mult"):
        plt.figure(figsize=(9, 7))

        Wmax = max(Whi * 1.3, g["W_blow"].max())
        W_grid = np.linspace(1e-8, Wmax, 500)
        plt.plot(np.sqrt(W_grid), W_grid, color="black", linewidth=2, label=r"$V=+\sqrt{W}$")
        plt.plot(-np.sqrt(W_grid), W_grid, color="black", linewidth=2, label=r"$V=-\sqrt{W}$")
        plt.plot(q * np.sqrt(W_grid), W_grid, linestyle="--", color="black", label=rf"$V={q:g}\sqrt{{W}}$ target")
        plt.axhspan(Wlo, Whi, alpha=0.15, label="window")

        for path, gp in g.groupby("path"):
            local = gp[(gp["W_blow"] > -0.02) & (gp["W_blow"] < Wmax)]
            plt.plot(local["V_blow"], local["W_blow"], linewidth=1.2, alpha=0.85)

        hs = hits_df[(hits_df["sigma_mult"] == m) & (hits_df["hit"])]
        if len(hs):
            plt.scatter(hs["hit_V"], hs["hit_W"], marker="x", s=80, label="first local hit")

        plt.xlabel(r"$V=(v+1)/\epsilon^{1/3}$")
        plt.ylabel(r"$W=(w-w_f)/\epsilon^{2/3}$")
        plt.title(f"Blown-up local geometry, sigma_mult={m:g}")
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(outdir / f"layer3_blowup_sigma_mult_{m:g}.png", dpi=180)
        plt.close()

    print(f"\nSaved outputs to {outdir}")


if __name__ == "__main__":
    main()
