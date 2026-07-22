#!/usr/bin/env python3
"""
canard_escape_autonomous.py
---------------------------
Autonomous canard-escape probe for the stochastic FHN project.

Purpose
-------
This is NOT the dynamic-ramp canard test.  It fixes I and asks whether noise
causes earlier departure from the local canard / repelling slow-manifold
neighbourhood.  It is designed to test the blow-up/FW prediction

    sigma_crit ~ [b (I - I_fold_L)]^(2/3) * sqrt(2 eps)

by extracting a crude sigma_50 boundary and fitting

    sigma_50 ~ eps^alpha                 expected alpha ~= 1/2
    sigma_50 / sqrt(eps) ~ lambda^beta   expected beta  ~= 2/3.

The observable is intentionally diagnostic:
  - start from a Krupa-Szmolyan blow-up initial condition near the left fold;
  - simulate the full FHN SDE at fixed I;
  - measure when each trajectory leaves a tube around the repelling branch
    V = +sqrt(W) in blow-up coordinates;
  - define 'early peel-off' relative to the deterministic sigma=0 baseline.

Because autonomous canard observables are delicate, the CSV includes raw metrics
so you can inspect whether the boundary is meaningful before trusting exponents.

Usage
-----
Quick smoke test:
    python3 canard_escape_autonomous.py --quick

Full-ish first pass:
    python3 canard_escape_autonomous.py

Fit an existing CSV only:
    python3 canard_escape_autonomous.py --fit-only --csv results/canard_escape_autonomous/results.csv

Outputs
-------
    results/canard_escape_autonomous/results.csv
    results/canard_escape_autonomous/sigma50.csv
    results/canard_escape_autonomous/*.png
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from pathlib import Path
from dataclasses import dataclass

import numpy as np

# Import local kernel.py from the project root / script directory.
# Put this script in the same directory as kernel.py, or run it from the project root.
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(Path.cwd()))

try:
    from kernel import FHN2D
except Exception as e:
    raise RuntimeError(
        "Could not import FHN2D from kernel.py. Put this script in the project root "
        "beside kernel.py, or run from that directory. Original error: " + repr(e)
    )


A_DEFAULT = 0.7
B_DEFAULT = 0.8


@dataclass
class PointResult:
    deltaI: float
    lam: float
    I: float
    eps: float
    sigma: float
    n: int
    frac_fired: float
    frac_peeled: float
    frac_early: float
    mean_t_peel: float
    mean_W_peel: float
    mean_t_fire: float
    det_t_peel: float
    det_W_peel: float
    pred_sigma: float
    T_phys: float


def blowup_initial_condition(I: float, eps: float, lam: float, W_factor: float) -> tuple[float, float, float, float]:
    """
    Start on the attracting branch of the local blow-up parabola W = V^2:

        V0 = -sqrt(W0),  W0 = W_factor * lambda^(2/3).

    Physical coordinates:
        v = -1 + eps^(1/3) V,
        w = w_f + eps^(2/3) W,  w_f = I - 2/3.
    """
    W0 = W_factor * max(lam, 1e-12) ** (2.0 / 3.0)
    V0 = -math.sqrt(W0)
    v0 = -1.0 + eps ** (1.0 / 3.0) * V0
    wf = I - 2.0 / 3.0
    w0 = wf + eps ** (2.0 / 3.0) * W0
    return v0, w0, V0, W0


def simulate_peeloff(
    *,
    I: float,
    eps: float,
    sigma: float,
    a: float,
    b: float,
    lam: float,
    n: int,
    dt: float,
    T_blow: float,
    W_factor: float,
    tube_tol: float,
    min_repelling_time: float,
    threshold: float,
    seed: int,
) -> dict:
    """
    Vectorized Euler-Maruyama simulation from a local blow-up IC.

    Returns peel/fire summaries.  Peel-off is detected after entering the
    repelling side (V>0,W>0), when the trajectory leaves the tube around
    V=sqrt(W), or when it fires.
    """
    rng = np.random.default_rng(seed)
    v0, w0, V0, W0 = blowup_initial_condition(I, eps, lam, W_factor)

    # t = eps^{-1/3} T, so integrate O(1) time in blow-up coordinates.
    T_phys = float(T_blow * eps ** (-1.0 / 3.0))
    n_steps = int(max(1, math.ceil(T_phys / dt)))
    sqrt_dt = math.sqrt(dt)
    wf = I - 2.0 / 3.0
    eps13 = eps ** (1.0 / 3.0)
    eps23 = eps ** (2.0 / 3.0)

    v = np.full(n, v0, dtype=float)
    w = np.full(n, w0, dtype=float)

    entered_rep = np.zeros(n, dtype=bool)
    t_enter_rep = np.full(n, np.nan)
    peeled = np.zeros(n, dtype=bool)
    t_peel = np.full(n, np.nan)
    W_peel = np.full(n, np.nan)
    fired = np.zeros(n, dtype=bool)
    t_fire = np.full(n, np.nan)

    # A small absolute floor prevents absurdly narrow tube near W=0.
    # In blow-up coordinates, this is still a local test.
    abs_floor = 0.08

    for k in range(n_steps):
        t = (k + 1) * dt
        noise = rng.standard_normal(n) if sigma > 0 else 0.0
        v += (v - v**3 / 3.0 - w + I) * dt + sigma * sqrt_dt * noise
        w += eps * (v + a - b * w) * dt

        V = (v + 1.0) / eps13
        W = (w - wf) / eps23
        sqrtW = np.sqrt(np.maximum(W, 0.0))
        dist_rep = np.abs(V - sqrtW)
        tube_width = tube_tol * np.maximum(sqrtW, abs_floor)

        newly_fire = (~fired) & (v >= threshold)
        if np.any(newly_fire):
            fired[newly_fire] = True
            t_fire[newly_fire] = t

        newly_enter = (~entered_rep) & (V > 0.0) & (W > 0.0)
        if np.any(newly_enter):
            entered_rep[newly_enter] = True
            t_enter_rep[newly_enter] = t

        eligible = entered_rep & (~peeled) & ((t - t_enter_rep) >= min_repelling_time)
        # Either leaves tube around repelling branch while W is still local, or fires.
        newly_peel = eligible & ((dist_rep > tube_width) | fired)
        if np.any(newly_peel):
            peeled[newly_peel] = True
            t_peel[newly_peel] = t
            W_peel[newly_peel] = W[newly_peel]

        if np.all(peeled | fired):
            # If fired trajectories have not been marked as peeled because they never
            # formally entered the repelling tube, mark them separately at firing.
            fired_unpeeled = fired & (~peeled)
            if np.any(fired_unpeeled):
                peeled[fired_unpeeled] = True
                t_peel[fired_unpeeled] = t_fire[fired_unpeeled]
                V_now = (v[fired_unpeeled] + 1.0) / eps13
                W_now = (w[fired_unpeeled] - wf) / eps23
                W_peel[fired_unpeeled] = W_now
            break

    return {
        "frac_fired": float(np.mean(fired)),
        "frac_peeled": float(np.mean(peeled)),
        "t_peel": t_peel,
        "W_peel": W_peel,
        "t_fire": t_fire,
        "T_phys": T_phys,
        "v0": v0,
        "w0": w0,
        "V0": V0,
        "W0": W0,
    }


def safe_nanmean(x: np.ndarray) -> float:
    return float(np.nanmean(x)) if np.any(~np.isnan(x)) else float("nan")


def predicted_sigma(eps: float, lam: float) -> float:
    # sigma_crit = lambda^(2/3) * sqrt(2 eps)
    return float((lam ** (2.0 / 3.0)) * math.sqrt(2.0 * eps))


def run_grid(args) -> list[PointResult]:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = Path(args.csv)
    csv_path.parent.mkdir(parents=True, exist_ok=True)

    ref = FHN2D(I=0.0, a=args.a, b=args.b)
    I_fold_L = ref.I_snic  # backward-compatible alias in kernel.py

    print("\n=== Autonomous canard-escape sweep ===")
    print(f"a={args.a}, b={args.b}")
    print(f"I_fold_L = {I_fold_L:.8f}")
    print(f"deltaI values: {args.deltaI}")
    print(f"eps values   : {args.eps}")
    print(f"sigma values : {args.sigmas}")
    print(f"n={args.n}, dt={args.dt}, T_blow={args.T_blow}\n")

    rows: list[PointResult] = []

    write_header = not csv_path.exists() or args.overwrite
    mode = "w" if args.overwrite else "a"
    with open(csv_path, mode, newline="") as fh:
        wr = csv.writer(fh)
        if write_header:
            wr.writerow([
                "deltaI", "lambda", "I", "eps", "sigma", "n",
                "frac_fired", "frac_peeled", "frac_early",
                "mean_t_peel", "mean_W_peel", "mean_t_fire",
                "det_t_peel", "det_W_peel", "pred_sigma", "T_phys",
            ])

        for dI in args.deltaI:
            I = I_fold_L + dI
            lam = args.b * dI
            for eps in args.eps:
                # Deterministic baseline for this (I, eps).
                det = simulate_peeloff(
                    I=I, eps=eps, sigma=0.0, a=args.a, b=args.b, lam=lam,
                    n=1, dt=args.dt, T_blow=args.T_blow, W_factor=args.W_factor,
                    tube_tol=args.tube_tol, min_repelling_time=args.min_repelling_time,
                    threshold=args.threshold, seed=args.seed,
                )
                det_t = safe_nanmean(det["t_peel"])
                det_W = safe_nanmean(det["W_peel"])
                pred = predicted_sigma(eps, lam)

                print(f"\n--- deltaI={dI:.5f}, lambda={lam:.5f}, I={I:.6f}, eps={eps:.4f} ---")
                print(f"pred sigma={pred:.5f}; deterministic t_peel={det_t:.5g}, W_peel={det_W:.5g}")

                for sigma in args.sigmas:
                    # Reuse deterministic result if sigma=0 and n=1 baseline is enough? No:
                    # run n trajectories at sigma=0 too so CSV has matching row.
                    res = simulate_peeloff(
                        I=I, eps=eps, sigma=sigma, a=args.a, b=args.b, lam=lam,
                        n=args.n, dt=args.dt, T_blow=args.T_blow, W_factor=args.W_factor,
                        tube_tol=args.tube_tol, min_repelling_time=args.min_repelling_time,
                        threshold=args.threshold, seed=args.seed + int(1e6*dI) + int(1e4*eps) + int(1e5*sigma),
                    )
                    Wp = res["W_peel"]
                    tp = res["t_peel"]
                    tf = res["t_fire"]

                    # Early means: peel occurs meaningfully before the deterministic peel-off.
                    # Use relative tolerances so tiny stochastic jitter around the deterministic
                    # peel-off is not counted as a real canard-destroying event.

                    t_margin = max(args.early_t_margin, args.early_frac * abs(det_t)) if np.isfinite(det_t) else args.early_t_margin
                    W_margin = max(args.early_W_margin, args.early_W_frac * abs(det_W)) if np.isfinite(det_W) else args.early_W_margin

                    early_by_W = (~np.isnan(Wp)) & np.isfinite(det_W) & (Wp > det_W + W_margin)
                    early_by_t = (~np.isnan(tp)) & np.isfinite(det_t) & (tp < det_t - t_margin)

                    early = early_by_W | early_by_t
                    frac_early = float(np.mean(early))

                    pr = PointResult(
                        deltaI=float(dI), lam=float(lam), I=float(I), eps=float(eps), sigma=float(sigma), n=int(args.n),
                        frac_fired=float(res["frac_fired"]), frac_peeled=float(res["frac_peeled"]), frac_early=frac_early,
                        mean_t_peel=safe_nanmean(tp), mean_W_peel=safe_nanmean(Wp), mean_t_fire=safe_nanmean(tf),
                        det_t_peel=float(det_t), det_W_peel=float(det_W), pred_sigma=pred, T_phys=float(res["T_phys"]),
                    )
                    rows.append(pr)
                    wr.writerow([
                        pr.deltaI, pr.lam, pr.I, pr.eps, pr.sigma, pr.n,
                        pr.frac_fired, pr.frac_peeled, pr.frac_early,
                        pr.mean_t_peel, pr.mean_W_peel, pr.mean_t_fire,
                        pr.det_t_peel, pr.det_W_peel, pr.pred_sigma, pr.T_phys,
                    ])
                    fh.flush()
                    print(
                        f" sigma={sigma:.5f} | early={frac_early:.3f} "
                        f"peeled={pr.frac_peeled:.3f} fired={pr.frac_fired:.3f} "
                        f"meanW={pr.mean_W_peel:.4g} meant={pr.mean_t_peel:.4g}"
                    )

    return rows


def read_csv(csv_path: str | Path) -> list[PointResult]:
    rows = []
    with open(csv_path) as fh:
        for r in csv.DictReader(fh):
            rows.append(PointResult(
                deltaI=float(r["deltaI"]), lam=float(r["lambda"]), I=float(r["I"]), eps=float(r["eps"]),
                sigma=float(r["sigma"]), n=int(float(r["n"])), frac_fired=float(r["frac_fired"]),
                frac_peeled=float(r["frac_peeled"]), frac_early=float(r["frac_early"]),
                mean_t_peel=float(r["mean_t_peel"]), mean_W_peel=float(r["mean_W_peel"]), mean_t_fire=float(r["mean_t_fire"]),
                det_t_peel=float(r["det_t_peel"]), det_W_peel=float(r["det_W_peel"]),
                pred_sigma=float(r["pred_sigma"]), T_phys=float(r["T_phys"]),
            ))
    return rows


def interpolate_sigma50(sigmas: np.ndarray, probs: np.ndarray, target: float = 0.5) -> float:
    order = np.argsort(sigmas)
    s = sigmas[order]
    p = probs[order]
    mask = np.isfinite(s) & np.isfinite(p)
    s, p = s[mask], p[mask]
    if len(s) == 0:
        return float("nan")
    above = np.where(p >= target)[0]
    if len(above) == 0:
        return float("nan")
    k = above[0]
    if k == 0:
        return float(s[0])
    if p[k] == p[k-1]:
        return float(s[k])
    return float(s[k-1] + (target - p[k-1]) * (s[k] - s[k-1]) / (p[k] - p[k-1]))


def fit_power(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    mask = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return float("nan"), float("nan"), float("nan")
    lx, ly = np.log(x[mask]), np.log(y[mask])
    p, logC = np.polyfit(lx, ly, 1)
    pred = p * lx + logC
    ss_res = float(np.sum((ly - pred) ** 2))
    ss_tot = float(np.sum((ly - np.mean(ly)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(p), float(math.exp(logC)), r2


def analyse(rows: list[PointResult], args) -> None:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    keys = sorted({(r.deltaI, r.lam, r.I, r.eps) for r in rows})
    sigma50_rows = []
    print("\n=== sigma_50 extraction from frac_early ===")
    for dI, lam, I, eps in keys:
        subset = [r for r in rows if r.deltaI == dI and r.eps == eps]
        sig = np.array([r.sigma for r in subset], dtype=float)
        prob = np.array([r.frac_early for r in subset], dtype=float)
        sc50 = interpolate_sigma50(sig, prob, target=args.target)
        pred = predicted_sigma(eps, lam)
        sigma50_rows.append((dI, lam, I, eps, sc50, pred, len(subset)))
        print(f"deltaI={dI:.5f} eps={eps:.4f}: sigma{int(args.target*100)}={sc50:.5g}, pred={pred:.5g}, ratio={sc50/pred if pred>0 and np.isfinite(sc50) else np.nan:.3g}")

    sigma50_path = outdir / "sigma50.csv"
    with open(sigma50_path, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["deltaI", "lambda", "I", "eps", f"sigma{int(args.target*100)}", "pred_sigma", "n_sigma_points"])
        wr.writerows(sigma50_rows)
    print(f"\n-> wrote {sigma50_path}")

    # Fits by deltaI: sigma50 ~ eps^alpha.
    print("\n=== Fit by fixed lambda: sigma50 ~ C eps^alpha ===")
    for dI in sorted({r[0] for r in sigma50_rows}):
        arr = np.array([[eps, sc] for dd, lam, I, eps, sc, pred, npt in sigma50_rows if dd == dI], dtype=float)
        if len(arr):
            alpha, C, r2 = fit_power(arr[:, 0], arr[:, 1])
            print(f"deltaI={dI:.5f}: alpha={alpha:.3f}, C={C:.4g}, R2={r2:.3f}   (target alpha=0.5)")

    # Fits by eps: sigma50/sqrt(eps) ~ lambda^beta.
    print("\n=== Fit by fixed eps: sigma50/sqrt(eps) ~ C lambda^beta ===")
    for eps in sorted({r[3] for r in sigma50_rows}):
        arr = np.array([[lam, sc / math.sqrt(eps)] for dd, lam, I, ee, sc, pred, npt in sigma50_rows if ee == eps], dtype=float)
        if len(arr):
            beta, C, r2 = fit_power(arr[:, 0], arr[:, 1])
            print(f"eps={eps:.4f}: beta={beta:.3f}, C={C:.4g}, R2={r2:.3f}   (target beta=2/3)")

    # Overall collapse: sigma50 / pred_sigma should be O(1), preferably flat.
    valid = np.array([[eps, lam, sc, pred] for dI, lam, I, eps, sc, pred, npt in sigma50_rows if np.isfinite(sc) and sc > 0 and pred > 0])
    if valid.shape[0] >= 3:
        ratio = valid[:, 2] / valid[:, 3]
        print("\n=== Collapse check ===")
        print(f"median sigma50/pred = {np.median(ratio):.3g}")
        print(f"range  sigma50/pred = [{np.min(ratio):.3g}, {np.max(ratio):.3g}]")

    # Plots are optional; skip cleanly if matplotlib unavailable.
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # Plot probability curves.
        for dI in sorted({r.deltaI for r in rows}):
            plt.figure(figsize=(7.2, 5.0))
            for eps in sorted({r.eps for r in rows if r.deltaI == dI}):
                subset = sorted([r for r in rows if r.deltaI == dI and r.eps == eps], key=lambda z: z.sigma)
                plt.plot([r.sigma for r in subset], [r.frac_early for r in subset], "o-", label=f"eps={eps:g}")
            plt.axhline(args.target, color="k", ls=":", lw=1)
            plt.xlabel("sigma")
            plt.ylabel("P(early peel-off)")
            plt.title(f"Autonomous canard escape curves, deltaI={dI:g}")
            plt.legend(fontsize=8)
            plt.tight_layout()
            p = outdir / f"early_probability_deltaI_{dI:g}.png"
            plt.savefig(p, dpi=160)
            plt.close()
            print(f"-> wrote {p}")

        # sigma50 vs eps by deltaI.
        plt.figure(figsize=(7.2, 5.0))
        for dI in sorted({r[0] for r in sigma50_rows}):
            arr = np.array([[eps, sc, pred] for dd, lam, I, eps, sc, pred, npt in sigma50_rows if dd == dI and np.isfinite(sc) and sc > 0])
            if len(arr):
                plt.loglog(arr[:, 0], arr[:, 1], "o-", label=f"measured dI={dI:g}")
                plt.loglog(arr[:, 0], arr[:, 2], "--", alpha=0.6, label=f"pred dI={dI:g}")
        plt.xlabel("eps")
        plt.ylabel(f"sigma{int(args.target*100)}")
        plt.title("Autonomous canard escape: measured vs lambda^(2/3) sqrt(2 eps)")
        plt.legend(fontsize=7)
        plt.tight_layout()
        p = outdir / "sigma50_vs_eps.png"
        plt.savefig(p, dpi=160)
        plt.close()
        print(f"-> wrote {p}")

        # Collapse plot.
        if valid.shape[0] >= 3:
            plt.figure(figsize=(6.4, 4.6))
            plt.semilogx(valid[:, 3], valid[:, 2] / valid[:, 3], "o")
            plt.axhline(np.median(valid[:, 2] / valid[:, 3]), color="k", ls="--", lw=1)
            plt.xlabel("predicted sigma")
            plt.ylabel("measured / predicted")
            plt.title("Prediction collapse check")
            plt.tight_layout()
            p = outdir / "collapse_ratio.png"
            plt.savefig(p, dpi=160)
            plt.close()
            print(f"-> wrote {p}")
    except Exception as e:
        print(f"[plotting skipped: {e!r}]")


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="small smoke-test grid")
    ap.add_argument("--fit-only", action="store_true", help="do not run simulations; analyse existing CSV")
    ap.add_argument("--overwrite", action="store_true", help="overwrite CSV instead of appending")
    ap.add_argument("--a", type=float, default=A_DEFAULT)
    ap.add_argument("--b", type=float, default=B_DEFAULT)
    ap.add_argument("--deltaI", type=float, nargs="+", default=None)
    ap.add_argument("--eps", type=float, nargs="+", default=None)
    ap.add_argument("--sigmas", type=float, nargs="+", default=None)
    ap.add_argument("--n", type=int, default=None)
    ap.add_argument("--dt", type=float, default=2.5e-3)
    ap.add_argument("--T-blow", type=float, default=10.0, help="integration horizon in blow-up T units")
    ap.add_argument("--W-factor", type=float, default=3.0, help="start W0 = W_factor * lambda^(2/3)")
    ap.add_argument("--tube-tol", type=float, default=1.0, help="relative tube width around V=sqrt(W)")
    ap.add_argument("--min-repelling-time", type=float, default=0.02)
    ap.add_argument("--early-W-margin", type=float, default=0.02, help="early if W_peel > W_det + margin")
    ap.add_argument("--early-t-margin", type=float, default=0.02, help="early if t_peel < t_det - margin")
    ap.add_argument("--threshold", type=float, default=1.0)
    ap.add_argument("--target", type=float, default=0.5, help="probability level for sigma_target")
    ap.add_argument("--seed", type=int, default=123)
    ap.add_argument("--outdir", type=str, default="results/canard_escape_autonomous")
    ap.add_argument("--csv", type=str, default="results/canard_escape_autonomous/results.csv")
    ap.add_argument("--early-frac", type=float, default=0.10,
                help="relative time margin: early if t_peel < t_det*(1-early_frac)")
    ap.add_argument("--early-W-frac", type=float, default=0.10,
                help="relative W margin: early if W_peel > W_det + early_W_frac*abs(W_det)")
    return ap


def finalize_args(args):
    if args.quick:
        if args.deltaI is None:
            args.deltaI = [0.02]
        if args.eps is None:
            args.eps = [0.02, 0.04, 0.08]
        if args.sigmas is None:
            args.sigmas = [0.0, 0.006, 0.010, 0.014, 0.018, 0.024, 0.032, 0.045, 0.060]
        if args.n is None:
            args.n = 250
        # Keep quick run shorter.
        args.T_blow = min(args.T_blow, 8.0)
    else:
        if args.deltaI is None:
            args.deltaI = [0.01, 0.02, 0.04]
        if args.eps is None:
            args.eps = [0.01, 0.015, 0.02, 0.03, 0.04, 0.06, 0.08]
        if args.sigmas is None:
            args.sigmas = [0.0, 0.003, 0.005, 0.007, 0.010, 0.014, 0.018, 0.024, 0.032, 0.045, 0.060, 0.080]
        if args.n is None:
            args.n = 800
    args.deltaI = [float(x) for x in args.deltaI]
    args.eps = [float(x) for x in args.eps]
    args.sigmas = [float(x) for x in args.sigmas]
    args.n = int(args.n)
    return args


def main():
    args = finalize_args(build_parser().parse_args())
    if args.fit_only:
        rows = read_csv(args.csv)
    else:
        rows = run_grid(args)
    analyse(rows, args)


if __name__ == "__main__":
    main()
