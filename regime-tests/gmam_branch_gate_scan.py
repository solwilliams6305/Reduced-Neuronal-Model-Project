#!/usr/bin/env python3
"""
gmam_branch_gate_scan.py

Immediate Test 2: is the apparent fractional barrier law connected to nongradient
escape-path branch competition / moving-gate geometry?

This script is intentionally self-contained. It scans Hopf deficits delta for one
or more (a,b) pairs and computes geometric minimum-action paths (gMAM-style) from
multiple initialisations:

    straight, 1-turn, 2-turn, 4-turn spirals

For each case it records:
    - action from each initialisation
    - action gap between best and second-best branch
    - number of windings in the minimised path
    - committed crossing phase at a mid-shell
    - gate phase drift versus delta
    - approximate barrier B_gMAM = S_geo / 2

Why this test matters
---------------------
If a 2/3-ish law is real for nongradient reasons, one plausible mechanism is a
Maier--Stein-type most-probable escape path (MPEP) selection layer: competing
instanton branches, caustics, or a moving exit gate. This script tests that.

Strong evidence for a nongradient branch/crossover mechanism:
    1. different initialisations converge to distinct local minima with close actions;
    2. the best branch changes as delta varies;
    3. the gate phase/radius drifts systematically with delta;
    4. S_geo/2 and the measured Arrhenius barrier have similar scaling.

Caveat
------
This is a practical diagnostic, not a polished theorem-prover. gMAM is numerical
and local-minimum-sensitive by design; that sensitivity is exactly what we are
probing here.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

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
    """Lower-branch Hopf current for standard FHN fixed point."""
    arg = 1.0 - b * eps
    if arg <= 0:
        return float("nan")
    v = -math.sqrt(arg)
    w = (v + a) / b
    return w - v + v ** 3 / 3.0


def fixed_point_roots(I: float, a: float, b: float) -> np.ndarray:
    # Fixed point: w=(v+a)/b and v - v^3/3 - w + I = 0.
    # -b v^3 + (3b-3)v + (3bI-3a)=0
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
    """
    Geometric action:
        S = sum_i ( ||dX_i||_M ||b_i||_M - dX_i^T M b_i )
    with M=diag(1,lam), a finite penalty approximation for voltage-only noise.
    """
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
                  reparam_every: int = 10, tol_grad: float = 1e-10) -> tuple[np.ndarray, list[float]]:
    X = reparametrize(X0.copy(), lam)
    step = lr
    hist = []

    for it in range(iters):
        S, g = geo_action_and_grad(X, p, lam)
        hist.append(S)
        gn = float(np.sqrt(np.sum(g * g)))
        if gn < tol_grad:
            break

        accepted = False
        # Basic backtracking line search.
        for _ in range(25):
            Xtry = X - step * g
            Xtry[0] = X[0]; Xtry[-1] = X[-1]
            Stry, _ = geo_action_and_grad(Xtry, p, lam)
            if np.isfinite(Stry) and Stry < S:
                X = Xtry
                step *= 1.05
                accepted = True
                break
            step *= 0.5
        if not accepted and step < 1e-14:
            break

        if (it + 1) % reparam_every == 0:
            X = reparametrize(X, lam)

    X = reparametrize(X, lam)
    hist.append(geo_action_and_grad(X, p, lam)[0])
    return X, hist


def target_point(p: FHNParams, v_end: float) -> np.ndarray:
    # Put target on v-nullcline on spiking side.
    w_end = v_end - v_end ** 3 / 3.0 + p.I
    return np.array([v_end, w_end], dtype=float)


def initialise_path(start: np.ndarray, target: np.ndarray, basis: SpiralBasis, M: int,
                    n_turns: float, spiral: bool) -> np.ndarray:
    s = np.linspace(0.0, 1.0, M + 1)
    if not spiral or n_turns == 0:
        X = (1 - s)[:, None] * start + s[:, None] * target
        X[0] = start; X[-1] = target
        return X

    # Spiral in local Hopf coordinates from centre outward to target phase.
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


def crossing_at_radius(X: np.ndarray, basis: SpiralBasis, R_mid: float, which: str = "last") -> tuple[float, float]:
    ph, rad = phase_radius(X, basis)
    above = np.flatnonzero(rad >= R_mid)
    if len(above) == 0:
        return float("nan"), float("nan")
    if which == "first":
        k = int(above[0])
    else:
        below = np.flatnonzero(rad < R_mid)
        if len(below) == 0:
            k = int(above[0])
        else:
            k = min(int(below[-1]) + 1, len(rad) - 1)
    return float(ph[k]), float(rad[k])


def wrapped(x: float) -> float:
    return (x + math.pi) % (2 * math.pi) - math.pi


def run_case(a: float, b: float, eps: float, delta: float, v_end: float, M: int,
             lam: float, iters: int, lr: float, mid_frac: float,
             turns: list[float], out_paths: bool, outdir: Path) -> list[dict]:
    I_H = hopf_current(a, b, eps)
    I = I_H - delta
    p = FHNParams(I=I, a=a, b=b, eps=eps)
    v_fp, w_fp, eig = choose_stable_focus(p)
    basis = make_spiral_basis(v_fp, w_fp, p)

    start = basis.fp
    target = target_point(p, v_end)
    _, R_target_arr = phase_radius(target[None, :], basis)
    R_target = float(R_target_arr[0])
    R_mid = mid_frac * R_target

    rows = []
    branch_specs = [("straight", 0.0, False)] + [(f"spiral_{t:g}", float(t), True) for t in turns]

    path_dir = outdir / "paths"
    if out_paths:
        path_dir.mkdir(parents=True, exist_ok=True)

    for tag, nturns, spiral in branch_specs:
        X0 = initialise_path(start, target, basis, M, nturns, spiral)
        X, hist = gmam_minimise(X0, p, lam=lam, iters=iters, lr=lr)

        S = float(hist[-1])
        wind = count_windings(X, basis)
        ph_mid, R_cross = crossing_at_radius(X, basis, R_mid, which="last")
        ph_first, _ = crossing_at_radius(X, basis, R_mid, which="first")
        ph_t, _ = phase_radius(target[None, :], basis)
        target_phase = float(ph_t[0])

        row = {
            "a": a, "b": b, "eps": eps, "delta": delta,
            "I_hopf": I_H, "I": I,
            "v_fp": v_fp, "w_fp": w_fp,
            "alpha": basis.alpha, "omega": basis.omega,
            "kappa": abs(basis.alpha) / basis.omega if basis.omega else np.nan,
            "init": tag,
            "S_geo": S,
            "B_gmam_half": 0.5 * S,
            "windings": wind,
            "R_target": R_target,
            "R_mid": R_mid,
            "mid_phase_last": ph_mid,
            "mid_phase_first": ph_first,
            "mid_radius_crossed": R_cross,
            "target_phase": target_phase,
            "lead_mid_to_target": wrapped(ph_mid - target_phase) if np.isfinite(ph_mid) else np.nan,
            "n_iters_recorded": len(hist),
            "S_initial": float(hist[0]) if hist else np.nan,
            "S_final": S,
        }
        rows.append(row)

        if out_paths:
            fname = path_dir / f"path_a{a:g}_b{b:g}_delta{delta:g}_{tag}.csv"
            pd.DataFrame({"v": X[:, 0], "w": X[:, 1]}).to_csv(fname, index=False)

    return rows


def summarise_branches(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for key, g in df.groupby(["a", "b", "eps", "delta", "I"]):
        gs = g.sort_values("S_geo").reset_index(drop=True)
        best = gs.iloc[0]
        second = gs.iloc[1] if len(gs) > 1 else None
        gap = float(second["S_geo"] - best["S_geo"]) if second is not None else np.nan
        rel_gap = gap / max(abs(float(best["S_geo"])), 1e-300) if np.isfinite(gap) else np.nan
        rows.append({
            "a": key[0], "b": key[1], "eps": key[2], "delta": key[3], "I": key[4],
            "best_init": best["init"],
            "best_S_geo": best["S_geo"],
            "best_B_gmam_half": best["B_gmam_half"],
            "second_S_geo": float(second["S_geo"]) if second is not None else np.nan,
            "action_gap": gap,
            "relative_action_gap": rel_gap,
            "best_windings": best["windings"],
            "best_mid_phase_last": best["mid_phase_last"],
            "best_R_mid": best["R_mid"],
            "best_R_target": best["R_target"],
            "alpha": best["alpha"],
            "omega": best["omega"],
            "kappa": best["kappa"],
            "n_branches": len(gs),
            "phase_spread_across_inits": float(np.abs(np.mean(np.exp(1j * g["mid_phase_last"].dropna().to_numpy()))))
                                        if g["mid_phase_last"].notna().any() else np.nan,
            # Circular spread as 1-R; large means init-dependent phase.
            "circular_phase_spread_1_minus_R": 1.0 - float(np.abs(np.mean(np.exp(1j * g["mid_phase_last"].dropna().to_numpy()))))
                                        if g["mid_phase_last"].notna().any() else np.nan,
        })
    return pd.DataFrame(rows).sort_values(["a", "b", "delta"])


def fit_action_scaling(summary: pd.DataFrame) -> pd.DataFrame:
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
        rows.append({
            "a": a, "b": b,
            "alpha_gmam": float(alpha),
            "C_gmam": float(math.exp(logC)),
            "log_r2": 1.0 - ss_res / ss_tot,
            "n": len(g),
        })
    return pd.DataFrame(rows)


def make_plots(all_df: pd.DataFrame, summary: pd.DataFrame, scaling: pd.DataFrame,
               outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    # Action branches by delta.
    for (a, b), g in all_df.groupby(["a", "b"]):
        fig, ax = plt.subplots(figsize=(8, 5.2))
        for init, gi in g.groupby("init"):
            gi = gi.sort_values("delta")
            ax.plot(gi["delta"], gi["B_gmam_half"], marker="o", label=init)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("Hopf deficit δ")
        ax.set_ylabel("gMAM barrier estimate S_geo/2")
        ax.set_title(f"Action branches from different initialisations, a={a:g}, b={b:g}")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(outdir / f"gmam_branches_a{a:g}_b{b:g}.png", dpi=180)
        plt.close(fig)

    # Relative action gap.
    fig, ax = plt.subplots(figsize=(8, 5.2))
    for (a, b), g in summary.groupby(["a", "b"]):
        g = g.sort_values("delta")
        ax.plot(g["delta"], g["relative_action_gap"], marker="o", label=f"a={a:g}, b={b:g}")
    ax.set_xscale("log")
    ax.set_xlabel("δ")
    ax.set_ylabel("(S_second - S_best) / S_best")
    ax.set_title("Branch competition diagnostic: small gaps imply near-degenerate MPEPs")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "branch_action_gap.png", dpi=180)
    plt.close(fig)

    # Gate phase movement.
    fig, ax = plt.subplots(figsize=(8, 5.2))
    for (a, b), g in summary.groupby(["a", "b"]):
        g = g.sort_values("delta")
        ax.plot(g["delta"], g["best_mid_phase_last"], marker="o", label=f"a={a:g}, b={b:g}")
    ax.set_xscale("log")
    ax.set_xlabel("δ")
    ax.set_ylabel("best-path mid-shell phase")
    ax.set_title("Moving-gate diagnostic: instanton crossing phase vs δ")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "gate_phase_vs_delta.png", dpi=180)
    plt.close(fig)

    # Windings of best path.
    fig, ax = plt.subplots(figsize=(8, 5.2))
    for (a, b), g in summary.groupby(["a", "b"]):
        g = g.sort_values("delta")
        ax.plot(g["delta"], g["best_windings"], marker="o", label=f"a={a:g}, b={b:g}")
    ax.set_xscale("log")
    ax.set_xlabel("δ")
    ax.set_ylabel("windings in best path")
    ax.set_title("Does the preferred MPEP branch change with δ?")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "best_windings_vs_delta.png", dpi=180)
    plt.close(fig)

    # gMAM exponent by pair.
    if not scaling.empty:
        fig, ax = plt.subplots(figsize=(7.5, 4.8))
        labels = [f"a={r.a:g},b={r.b:g}" for _, r in scaling.iterrows()]
        ax.bar(np.arange(len(scaling)), scaling["alpha_gmam"])
        ax.axhline(2/3, linestyle="--", linewidth=1, label="2/3")
        ax.axhline(1.0, linestyle=":", linewidth=1, label="1")
        ax.set_xticks(np.arange(len(scaling)))
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_ylabel("gMAM barrier exponent")
        ax.set_title("Scaling of best gMAM action")
        ax.legend()
        fig.tight_layout()
        fig.savefig(outdir / "gmam_exponent_by_pair.png", dpi=180)
        plt.close(fig)


def parse_pairs(pair_strings: list[str] | None, a_values: list[float], b_values: list[float]) -> list[tuple[float, float]]:
    if pair_strings:
        pairs = []
        for s in pair_strings:
            if ":" not in s:
                raise ValueError(f"Pair must be formatted a:b, got {s}")
            a, b = s.split(":", 1)
            pairs.append((float(a), float(b)))
        return pairs
    return [(a, b) for a in a_values for b in b_values]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a-values", type=float, nargs="*", default=[0.6, 0.7, 0.8])
    ap.add_argument("--b-values", type=float, nargs="*", default=[0.7, 0.8, 0.9])
    ap.add_argument("--only-pairs", type=str, nargs="*", default=None,
                    help="Optional pair list like 0.6:0.7 0.7:0.8 0.8:0.9")
    ap.add_argument("--deltas", type=float, nargs="*", default=[0.015, 0.03, 0.05])
    ap.add_argument("--eps", type=float, default=0.08)
    ap.add_argument("--v-end", type=float, default=0.9)
    ap.add_argument("--M", type=int, default=300)
    ap.add_argument("--lam", type=float, default=2000.0)
    ap.add_argument("--iters", type=int, default=5000)
    ap.add_argument("--lr", type=float, default=2e-4)
    ap.add_argument("--turns", type=float, nargs="*", default=[1.0, 2.0, 4.0])
    ap.add_argument("--mid-frac", type=float, default=0.55)
    ap.add_argument("--save-paths", action="store_true")
    ap.add_argument("--outdir", type=str, default="results/gmam_branch_gate_scan")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    pairs = parse_pairs(args.only_pairs, args.a_values, args.b_values)

    all_rows = []
    for a, b in pairs:
        I_H = hopf_current(a, b, args.eps)
        print(f"\n=== a={a:g}, b={b:g}, eps={args.eps:g}, I_H={I_H:.6f} ===")
        for delta in args.deltas:
            try:
                rows = run_case(
                    a=a, b=b, eps=args.eps, delta=delta, v_end=args.v_end,
                    M=args.M, lam=args.lam, iters=args.iters, lr=args.lr,
                    mid_frac=args.mid_frac, turns=args.turns,
                    out_paths=args.save_paths, outdir=outdir,
                )
                all_rows.extend(rows)
                best = min(rows, key=lambda r: r["S_geo"])
                print(f"  delta={delta:g}: best={best['init']:>9s} "
                      f"S/2={best['B_gmam_half']:.4e} "
                      f"wind={best['windings']:.2f} "
                      f"phase={best['mid_phase_last']:+.3f}")
            except Exception as exc:
                print(f"  delta={delta:g}: FAILED: {exc}")

    if not all_rows:
        print("No successful gMAM cases.")
        return

    all_df = pd.DataFrame(all_rows).sort_values(["a", "b", "delta", "S_geo"])
    summary = summarise_branches(all_df)
    scaling = fit_action_scaling(summary)

    all_df.to_csv(outdir / "gmam_all_branches.csv", index=False)
    summary.to_csv(outdir / "gmam_branch_summary.csv", index=False)
    scaling.to_csv(outdir / "gmam_action_scaling.csv", index=False)

    make_plots(all_df, summary, scaling, outdir)

    print("\n================ GMAM BRANCH/GATE DIAGNOSTIC ================")
    print("\nBest branch summary:")
    print(summary[["a", "b", "delta", "best_init", "best_B_gmam_half",
                   "relative_action_gap", "best_windings", "best_mid_phase_last"]].to_string(index=False))

    if not scaling.empty:
        print("\nBest-path gMAM action scaling:")
        print(scaling.to_string(index=False))

    print(f"\nOutputs written to: {outdir}")
    print("=============================================================")


if __name__ == "__main__":
    main()
