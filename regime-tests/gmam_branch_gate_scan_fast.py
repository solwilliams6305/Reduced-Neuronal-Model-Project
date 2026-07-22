#!/usr/bin/env python3
"""
gmam_branch_gate_scan_fast.py

Faster/safer version of gmam_branch_gate_scan.py.

Why this exists
---------------
The full branch scan can be very slow because it runs every (a,b,delta,init)
serially. Worse, if the gMAM optimiser has not converged, the resulting actions
can look non-monotone in delta, which is a red flag.

This version:
    1. parallelises independent (a,b,delta,init) jobs;
    2. defaults to a smaller branch set: straight and 2-turn;
    3. writes each branch result as it completes;
    4. summarises best branch/gate/action after the run;
    5. includes convergence diagnostics:
          final action, initial action, relative action drop,
          final action / previous action change;
    6. includes monotonicity warnings for best action vs delta.

Use this as a quick diagnostic, not the final high-accuracy instanton solver.
"""

from __future__ import annotations

import argparse
import math
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
import traceback

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


@dataclass
class FHNParams:
    I: float
    a: float = 0.7
    b: float = 0.8
    eps: float = 0.08


@dataclass
class SpiralBasis:
    fp: np.ndarray
    B: np.ndarray
    Binv: np.ndarray
    alpha: float
    omega: float


def hopf_current(a: float, b: float, eps: float) -> float:
    arg = 1.0 - b * eps
    if arg <= 0:
        return float("nan")
    v = -math.sqrt(arg)
    w = (v + a) / b
    return w - v + v ** 3 / 3.0


def fixed_point_roots(I: float, a: float, b: float) -> np.ndarray:
    coeff = np.array([-b, 0.0, 3.0 * b - 3.0, 3.0 * b * I - 3.0 * a], dtype=float)
    roots = np.roots(coeff)
    return np.sort(roots[np.abs(roots.imag) < 1e-9].real)


def jacobian(v: float, p: FHNParams) -> np.ndarray:
    return np.array([[1.0 - v * v, -1.0],
                     [p.eps,       -p.eps * p.b]], dtype=float)


def choose_stable_focus(p: FHNParams) -> tuple[float, float, np.ndarray]:
    infos = []
    for v in fixed_point_roots(p.I, p.a, p.b):
        w = (v + p.a) / p.b
        eig = np.linalg.eigvals(jacobian(v, p))
        infos.append((float(np.max(eig.real)), float(v), float(w), eig))
    if not infos:
        raise RuntimeError(f"No fixed point at I={p.I}")
    infos.sort(key=lambda x: x[0])
    max_real, v, w, eig = infos[0]
    if max_real >= 0:
        raise RuntimeError(f"Selected fixed point is not stable at I={p.I}: max Re={max_real}")
    if np.max(np.abs(eig.imag)) < 1e-9:
        raise RuntimeError(f"Selected fixed point is not a focus at I={p.I}: eig={eig}")
    return v, w, eig


def make_spiral_basis(v_fp: float, w_fp: float, p: FHNParams) -> SpiralBasis:
    J = jacobian(v_fp, p)
    vals, vecs = np.linalg.eig(J)
    idxs = np.where(vals.imag > 1e-9)[0]
    idx = int(idxs[0]) if len(idxs) else int(np.argmax(np.abs(vals.imag)))
    evec = vecs[:, idx]
    B = np.column_stack([evec.real, evec.imag]).astype(float)
    if abs(np.linalg.det(B)) < 1e-12:
        raise RuntimeError("Degenerate eigenbasis")
    return SpiralBasis(
        fp=np.array([v_fp, w_fp], dtype=float),
        B=B,
        Binv=np.linalg.inv(B),
        alpha=float(vals[idx].real),
        omega=float(abs(vals[idx].imag)),
    )


def phase_radius(X: np.ndarray, basis: SpiralBasis) -> tuple[np.ndarray, np.ndarray]:
    Z = (basis.Binv @ (X - basis.fp).T).T
    ph = np.arctan2(Z[:, 1], Z[:, 0])
    rad = np.sqrt(Z[:, 0] ** 2 + Z[:, 1] ** 2)
    return ph, rad


def drift(X: np.ndarray, p: FHNParams) -> np.ndarray:
    v = X[:, 0]; w = X[:, 1]
    return np.column_stack([
        v - v ** 3 / 3.0 - w + p.I,
        p.eps * (v + p.a - p.b * w),
    ])


def geo_action_and_grad(X: np.ndarray, p: FHNParams, lam: float) -> tuple[float, np.ndarray]:
    M = np.array([1.0, lam])
    b = drift(X, p)
    dX = X[1:] - X[:-1]
    bi = b[:-1]
    nd = np.sqrt(np.sum(dX * dX * M, axis=1))
    nb = np.sqrt(np.sum(bi * bi * M, axis=1))
    dot = np.sum(dX * bi * M, axis=1)
    eps_n = 1e-12
    terms = nd * nb - dot
    S = float(np.sum(terms))

    dterm_dd = (nb / (nd + eps_n))[:, None] * (dX * M) - (bi * M)
    dterm_db = (nd / (nb + eps_n))[:, None] * (bi * M) - (dX * M)

    grad = np.zeros_like(X)
    N = X.shape[0] - 1
    grad[:-1] += -dterm_dd
    grad[1:] += dterm_dd

    for i in range(N):
        v = X[i, 0]
        J = np.array([[1.0 - v * v, -1.0],
                      [p.eps,       -p.eps * p.b]], dtype=float)
        grad[i] += J.T @ dterm_db[i]

    grad[0] = 0.0
    grad[-1] = 0.0
    return S, grad


def reparametrize(X: np.ndarray, lam: float) -> np.ndarray:
    M = np.array([1.0, lam])
    seg = np.sqrt(np.sum((X[1:] - X[:-1]) ** 2 * M, axis=1))
    s = np.concatenate([[0.0], np.cumsum(seg)])
    if s[-1] <= 1e-14:
        return X
    s /= s[-1]
    snew = np.linspace(0.0, 1.0, len(X))
    Y = np.column_stack([
        np.interp(snew, s, X[:, 0]),
        np.interp(snew, s, X[:, 1]),
    ])
    Y[0] = X[0]; Y[-1] = X[-1]
    return Y


def gmam_minimise(X0: np.ndarray, p: FHNParams, lam: float, iters: int, lr: float,
                  reparam_every: int = 10) -> tuple[np.ndarray, list[float]]:
    X = reparametrize(X0.copy(), lam)
    step = lr
    hist = []
    last_improve_iter = 0

    for it in range(iters):
        S, g = geo_action_and_grad(X, p, lam)
        hist.append(S)
        gn = float(np.sqrt(np.sum(g * g)))
        if gn < 1e-10:
            break

        accepted = False
        for _ in range(20):
            Xtry = X - step * g
            Xtry[0] = X[0]; Xtry[-1] = X[-1]
            Stry, _ = geo_action_and_grad(Xtry, p, lam)
            if np.isfinite(Stry) and Stry < S:
                X = Xtry
                step *= 1.03
                accepted = True
                last_improve_iter = it
                break
            step *= 0.5

        if not accepted and step < 1e-14:
            break

        if (it + 1) % reparam_every == 0:
            X = reparametrize(X, lam)

        # If it hasn't improved for ages, it is probably stuck.
        if it - last_improve_iter > 500:
            break

    X = reparametrize(X, lam)
    hist.append(geo_action_and_grad(X, p, lam)[0])
    return X, hist


def target_point(p: FHNParams, v_end: float) -> np.ndarray:
    w_end = v_end - v_end ** 3 / 3.0 + p.I
    return np.array([v_end, w_end], dtype=float)


def initialise_path(start: np.ndarray, target: np.ndarray, basis: SpiralBasis, M: int,
                    n_turns: float) -> np.ndarray:
    s = np.linspace(0.0, 1.0, M + 1)
    if n_turns == 0:
        X = (1 - s)[:, None] * start + s[:, None] * target
        X[0] = start; X[-1] = target
        return X

    target_local = basis.Binv @ (target - basis.fp)
    R_target = float(np.sqrt(np.sum(target_local ** 2)))
    ph_target = float(math.atan2(target_local[1], target_local[0]))

    rr = s * R_target
    ang = ph_target + 2.0 * math.pi * n_turns * (s - 1.0)
    local = np.column_stack([rr * np.cos(ang), rr * np.sin(ang)])
    X = basis.fp[None, :] + local @ basis.B.T
    X[0] = start; X[-1] = target
    return X


def count_windings(X: np.ndarray, basis: SpiralBasis) -> float:
    ph, _ = phase_radius(X, basis)
    return float(np.sum(np.abs(np.diff(np.unwrap(ph)))) / (2.0 * math.pi))


def crossing_at_radius(X: np.ndarray, basis: SpiralBasis, R_mid: float) -> tuple[float, float]:
    ph, rad = phase_radius(X, basis)
    above = np.flatnonzero(rad >= R_mid)
    if len(above) == 0:
        return float("nan"), float("nan")
    below = np.flatnonzero(rad < R_mid)
    if len(below) == 0:
        k = int(above[0])
    else:
        k = min(int(below[-1]) + 1, len(rad) - 1)
    return float(ph[k]), float(rad[k])


def one_branch_job(job: dict) -> dict:
    try:
        a = job["a"]; b = job["b"]; eps = job["eps"]; delta = job["delta"]
        I_H = hopf_current(a, b, eps)
        I = I_H - delta
        p = FHNParams(I=I, a=a, b=b, eps=eps)
        v_fp, w_fp, eig = choose_stable_focus(p)
        basis = make_spiral_basis(v_fp, w_fp, p)
        start = basis.fp
        target = target_point(p, job["v_end"])

        _, R_target_arr = phase_radius(target[None, :], basis)
        R_target = float(R_target_arr[0])
        R_mid = job["mid_frac"] * R_target

        X0 = initialise_path(start, target, basis, job["M"], job["turns"])
        X, hist = gmam_minimise(X0, p, lam=job["lam"], iters=job["iters"], lr=job["lr"])

        S0 = float(hist[0])
        Sf = float(hist[-1])
        prev = float(hist[-2]) if len(hist) >= 2 else np.nan
        last_rel_change = abs(prev - Sf) / max(abs(Sf), 1e-300) if np.isfinite(prev) else np.nan
        rel_drop = (S0 - Sf) / max(abs(S0), 1e-300)

        ph_mid, R_cross = crossing_at_radius(X, basis, R_mid)
        ph_t, _ = phase_radius(target[None, :], basis)

        return {
            "ok": True,
            "a": a, "b": b, "eps": eps, "delta": delta,
            "I_hopf": I_H, "I": I,
            "init": job["tag"],
            "turns": job["turns"],
            "v_fp": v_fp, "w_fp": w_fp,
            "alpha": basis.alpha, "omega": basis.omega,
            "kappa": abs(basis.alpha) / basis.omega if basis.omega else np.nan,
            "S_geo": Sf,
            "B_gmam_half": 0.5 * Sf,
            "S_initial": S0,
            "rel_action_drop": rel_drop,
            "last_rel_change": last_rel_change,
            "n_hist": len(hist),
            "windings": count_windings(X, basis),
            "R_target": R_target,
            "R_mid": R_mid,
            "mid_phase_last": ph_mid,
            "mid_radius_crossed": R_cross,
            "target_phase": float(ph_t[0]),
            "error": "",
        }
    except Exception as exc:
        return {
            "ok": False,
            "a": job.get("a"), "b": job.get("b"), "eps": job.get("eps"),
            "delta": job.get("delta"), "init": job.get("tag"), "turns": job.get("turns"),
            "error": f"{exc}\n{traceback.format_exc()}",
        }


def parse_pairs(pair_strings, a_values, b_values):
    if pair_strings:
        out = []
        for s in pair_strings:
            a, b = s.split(":")
            out.append((float(a), float(b)))
        return out
    return [(a, b) for a in a_values for b in b_values]


def summarise(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    good = df[df["ok"] == True].copy()
    if good.empty:
        return pd.DataFrame()
    for key, g in good.groupby(["a", "b", "eps", "delta", "I"]):
        gs = g.sort_values("S_geo")
        best = gs.iloc[0]
        second = gs.iloc[1] if len(gs) > 1 else None
        gap = float(second["S_geo"] - best["S_geo"]) if second is not None else np.nan
        rows.append({
            "a": key[0], "b": key[1], "eps": key[2], "delta": key[3], "I": key[4],
            "best_init": best["init"],
            "best_B_gmam_half": best["B_gmam_half"],
            "second_B_gmam_half": float(0.5 * second["S_geo"]) if second is not None else np.nan,
            "relative_action_gap": gap / max(abs(float(best["S_geo"])), 1e-300) if np.isfinite(gap) else np.nan,
            "best_windings": best["windings"],
            "best_mid_phase_last": best["mid_phase_last"],
            "best_last_rel_change": best["last_rel_change"],
            "best_rel_action_drop": best["rel_action_drop"],
            "n_branches": len(gs),
        })
    out = pd.DataFrame(rows).sort_values(["a", "b", "delta"])
    # monotonicity flag per pair: count decreases in best action as delta increases.
    out["nonmonotone_drop_from_previous"] = False
    for (a, b), idx in out.groupby(["a", "b"]).groups.items():
        ids = list(idx)
        sub = out.loc[ids].sort_values("delta")
        vals = sub["best_B_gmam_half"].to_numpy(float)
        drops = np.r_[False, np.diff(vals) < -1e-12]
        out.loc[sub.index, "nonmonotone_drop_from_previous"] = drops
    return out


def fit_scaling(summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (a, b), g in summary.groupby(["a", "b"]):
        g = g.sort_values("delta")
        if len(g) < 2:
            continue
        x = np.log(g["delta"].to_numpy(float))
        y = np.log(g["best_B_gmam_half"].to_numpy(float))
        X = np.vstack([x, np.ones_like(x)]).T
        alpha, logC = np.linalg.lstsq(X, y, rcond=None)[0]
        pred = X @ np.array([alpha, logC])
        ss_res = float(np.sum((y - pred) ** 2))
        ss_tot = float(np.sum((y - y.mean()) ** 2)) + 1e-300
        rows.append({"a": a, "b": b, "alpha_gmam": float(alpha),
                     "C_gmam": float(math.exp(logC)), "log_r2": 1 - ss_res / ss_tot,
                     "n": len(g)})
    return pd.DataFrame(rows)


def make_plots(summary: pd.DataFrame, all_df: pd.DataFrame, scaling: pd.DataFrame, outdir: Path):
    good = all_df[all_df["ok"] == True].copy()
    if summary.empty or good.empty:
        return

    fig, ax = plt.subplots(figsize=(8,5))
    for (a,b), g in summary.groupby(["a","b"]):
        g = g.sort_values("delta")
        ax.plot(g["delta"], g["best_B_gmam_half"], marker="o", label=f"a={a:g}, b={b:g}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("δ")
    ax.set_ylabel("best gMAM S/2")
    ax.set_title("Best gMAM barrier vs Hopf deficit")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "best_gmam_barrier_vs_delta.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8,5))
    for (a,b), g in summary.groupby(["a","b"]):
        g = g.sort_values("delta")
        ax.plot(g["delta"], g["best_mid_phase_last"], marker="o", label=f"a={a:g}, b={b:g}")
    ax.set_xscale("log")
    ax.set_xlabel("δ")
    ax.set_ylabel("best-path mid phase")
    ax.set_title("Gate phase drift")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "gate_phase_vs_delta.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8,5))
    for (a,b), g in summary.groupby(["a","b"]):
        g = g.sort_values("delta")
        ax.plot(g["delta"], g["relative_action_gap"], marker="o", label=f"a={a:g}, b={b:g}")
    ax.set_xscale("log")
    ax.set_xlabel("δ")
    ax.set_ylabel("relative branch action gap")
    ax.set_title("Branch competition")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "relative_branch_gap.png", dpi=180)
    plt.close(fig)

    if not scaling.empty:
        fig, ax = plt.subplots(figsize=(7,4.5))
        labels = [f"a={r.a:g},b={r.b:g}" for _, r in scaling.iterrows()]
        ax.bar(np.arange(len(scaling)), scaling["alpha_gmam"])
        ax.axhline(2/3, linestyle="--", label="2/3")
        ax.axhline(1.0, linestyle=":", label="1")
        ax.set_xticks(np.arange(len(scaling)))
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_ylabel("gMAM exponent")
        ax.set_title("Best-path gMAM scaling exponent")
        ax.legend()
        fig.tight_layout()
        fig.savefig(outdir / "gmam_exponent_by_pair.png", dpi=180)
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a-values", type=float, nargs="*", default=[0.6,0.7,0.8])
    ap.add_argument("--b-values", type=float, nargs="*", default=[0.7,0.8,0.9])
    ap.add_argument("--only-pairs", type=str, nargs="*", default=None)
    ap.add_argument("--deltas", type=float, nargs="*", default=[0.015,0.03,0.05])
    ap.add_argument("--eps", type=float, default=0.08)
    ap.add_argument("--v-end", type=float, default=0.9)
    ap.add_argument("--M", type=int, default=200)
    ap.add_argument("--iters", type=int, default=2000)
    ap.add_argument("--lr", type=float, default=2e-4)
    ap.add_argument("--lam", type=float, default=2000.0)
    ap.add_argument("--turns", type=float, nargs="*", default=[0.0, 2.0],
                    help="Initial winding counts. 0 means straight.")
    ap.add_argument("--mid-frac", type=float, default=0.55)
    ap.add_argument("--n-workers", type=int, default=4)
    ap.add_argument("--outdir", type=str, default="results/gmam_branch_gate_scan_fast")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    pairs = parse_pairs(args.only_pairs, args.a_values, args.b_values)

    jobs = []
    for a,b in pairs:
        for delta in args.deltas:
            for turn in args.turns:
                tag = "straight" if abs(turn) < 1e-12 else f"spiral_{turn:g}"
                jobs.append({
                    "a": a, "b": b, "eps": args.eps, "delta": delta,
                    "v_end": args.v_end, "M": args.M, "iters": args.iters,
                    "lr": args.lr, "lam": args.lam, "turns": turn, "tag": tag,
                    "mid_frac": args.mid_frac,
                })

    print(f"Running {len(jobs)} branch jobs with {args.n_workers} workers...")
    rows = []
    if args.n_workers == 1:
        for j in jobs:
            r = one_branch_job(j)
            rows.append(r)
            status = "OK" if r["ok"] else "FAIL"
            if r["ok"]:
                print(f"{status} a={r['a']:g} b={r['b']:g} d={r['delta']:g} {r['init']} "
                      f"S/2={r['B_gmam_half']:.4e} wind={r['windings']:.2f} "
                      f"phase={r['mid_phase_last']:+.3f}")
            else:
                print(f"{status} a={r['a']} b={r['b']} d={r['delta']} {r['init']}: {r['error'].splitlines()[0]}")
    else:
        with ProcessPoolExecutor(max_workers=args.n_workers) as ex:
            futs = [ex.submit(one_branch_job, j) for j in jobs]
            for fut in as_completed(futs):
                r = fut.result()
                rows.append(r)
                status = "OK" if r["ok"] else "FAIL"
                if r["ok"]:
                    print(f"{status} a={r['a']:g} b={r['b']:g} d={r['delta']:g} {r['init']} "
                          f"S/2={r['B_gmam_half']:.4e} wind={r['windings']:.2f} "
                          f"phase={r['mid_phase_last']:+.3f}")
                else:
                    print(f"{status} a={r['a']} b={r['b']} d={r['delta']} {r['init']}: {r['error'].splitlines()[0]}")

    all_df = pd.DataFrame(rows)
    all_df.to_csv(outdir / "gmam_fast_all_branches.csv", index=False)

    summary = summarise(all_df)
    summary.to_csv(outdir / "gmam_fast_branch_summary.csv", index=False)

    scaling = fit_scaling(summary) if not summary.empty else pd.DataFrame()
    scaling.to_csv(outdir / "gmam_fast_action_scaling.csv", index=False)

    make_plots(summary, all_df, scaling, outdir)

    print("\n================ FAST GMAM SUMMARY ================")
    if not summary.empty:
        print(summary[["a","b","delta","best_init","best_B_gmam_half","relative_action_gap",
                       "best_windings","best_mid_phase_last","nonmonotone_drop_from_previous"]].to_string(index=False))
        n_bad = int(summary["nonmonotone_drop_from_previous"].sum())
        if n_bad:
            print(f"\nWARNING: {n_bad} non-monotone drops in best S/2 vs delta. "
                  "Treat gMAM action scaling as not converged/reliable yet.")
    if not scaling.empty:
        print("\nScaling:")
        print(scaling.to_string(index=False))
    print(f"\nOutputs written to: {outdir}")
    print("===================================================")


if __name__ == "__main__":
    main()
