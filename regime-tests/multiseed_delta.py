#!/usr/bin/env python3
"""
multiseed_delta.py

Multi-seed barrier measurement at the canonical (a,b,eps) over a delta grid, to
put Monte-Carlo error bars on B(delta) and test whether the exponent crossover
survives. Appends one row per (delta, seed) to a CSV so it can be run in chunks.
"""
from __future__ import annotations
import argparse, csv, math, os
from pathlib import Path
import numpy as np
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_resonator_phase_difference import (
    FHNParams, choose_stable_spiral_fp, make_phase_basis, simulate_many,
)
from arrhenius_escape import survival_rate


def hopf_I(a, b, eps):
    v = -math.sqrt(1 - b * eps); w = (v + a) / b
    return w - v + v ** 3 / 3.0


def measure_B(p, basis, sigmas, N, T, dt, r0, seed, S_hi, S_lo):
    xs, ys = [], []
    for s in sigmas:
        sim = simulate_many(p=p, basis=basis, sigma=s, N=N, T=T, dt=dt,
                            r0=r0, phi0=0.0, spike_threshold=1.0, seed=seed)
        res = survival_rate(sim["spike_time"], sim["spiked"], sim["t"], S_hi, S_lo)
        if res is None or res["k"] <= 0:
            continue
        xs.append(1.0 / s ** 2); ys.append(math.log(res["k"]))
    if len(xs) < 3:
        return None
    x = np.array(xs); y = np.array(ys)
    A = np.vstack([x, np.ones_like(x)]).T
    slope, icpt = np.linalg.lstsq(A, y, rcond=None)[0]
    pred = A @ np.array([slope, icpt])
    r2 = 1 - ((y - pred) ** 2).sum() / (((y - y.mean()) ** 2).sum() + 1e-30)
    return {"B": -float(slope), "lnA": float(icpt), "r2": float(r2), "npts": len(xs)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", type=float, default=0.7)
    ap.add_argument("--b", type=float, default=0.8)
    ap.add_argument("--eps", type=float, default=0.08)
    ap.add_argument("--deltas", type=float, nargs="+", required=True)
    ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3])
    ap.add_argument("--sigmas", type=float, nargs="+",
                    default=[0.04, 0.05, 0.06, 0.07, 0.08])
    ap.add_argument("--N", type=int, default=300)
    ap.add_argument("--T", type=float, default=180.0)
    ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--r0", type=float, default=0.06)
    ap.add_argument("--S-hi", type=float, default=0.85)
    ap.add_argument("--S-lo", type=float, default=0.10)
    ap.add_argument("--csv", type=str, default="results/multiseed_delta/barriers.csv")
    args = ap.parse_args()

    out = Path(args.csv); out.parent.mkdir(parents=True, exist_ok=True)
    fields = ["a", "b", "eps", "delta", "I", "seed", "B", "lnA", "r2", "npts"]
    done = set()
    if out.exists():
        with open(out) as f:
            for r in csv.DictReader(f):
                done.add((float(r["delta"]), int(r["seed"])))
    new = not out.exists()
    fh = open(out, "a", newline="")
    wr = csv.DictWriter(fh, fieldnames=fields)
    if new:
        wr.writeheader()
    for d in args.deltas:
        I = hopf_I(args.a, args.b, args.eps) - d
        p = FHNParams(I=I, a=args.a, b=args.b, eps=args.eps)
        v_fp, w_fp, cls = choose_stable_spiral_fp(p)
        basis = make_phase_basis(v_fp, w_fp, p)
        for sd in args.seeds:
            if (d, sd) in done:
                print(f"  skip delta={d} seed={sd} (done)"); continue
            res = measure_B(p, basis, args.sigmas, args.N, args.T, args.dt,
                            args.r0, sd, args.S_hi, args.S_lo)
            if res is None:
                print(f"  delta={d} seed={sd}: insufficient, skip"); continue
            row = {"a": args.a, "b": args.b, "eps": args.eps, "delta": d,
                   "I": round(I, 6), "seed": sd, "B": res["B"], "lnA": res["lnA"],
                   "r2": res["r2"], "npts": res["npts"]}
            wr.writerow(row); fh.flush()
            print(f"  delta={d} seed={sd}: B={res['B']:.4e} r2={res['r2']:.3f}")
    fh.close()
    print(f"-> {out}")


if __name__ == "__main__":
    main()
