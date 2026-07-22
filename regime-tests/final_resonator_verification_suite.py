#!/usr/bin/env python3
"""
Final resonator verification suite for stochastic FitzHugh--Nagumo.

Purpose
-------
This script runs the last sanity checks needed before making a near-Hopf
phase-gated stochastic commitment claim:

1. Deterministic fixed-point / eigenvalue scan in I.
2. Numerical Hopf localisation from Re(lambda)=0.
3. Stable-side near-Hopf sigma sweeps at one or more horizons T.
4. Sigma50 extraction from excursion probabilities.
5. Commitment-vs-threshold phase concentration checks.
6. Rotations-before-commitment summaries.

Model
-----
    dv = (v - v^3/3 - w + I) dt + sigma dW
    dw = eps (v + a - b w) dt

Noise acts only in v.

Commitment definition
---------------------
For each trajectory, we define a local spiral coordinate around the stable
fixed point using the real and imaginary parts of the complex eigenvector.

A candidate commitment event occurs when the local spiral radius crosses
`commit_radius` outward. That candidate is accepted if the trajectory later
hits the large-excursion threshold `v_excursion` before returning inside
`reset_radius`. This is intentionally earlier than the large-excursion
threshold crossing, so it tests whether the phase gate is upstream of the
spike upstroke rather than a threshold-crossing artefact.

Outputs
-------
Inside --outdir:
    deterministic_scan.csv
    deterministic_scan.png
    hopf_estimate.txt
    all_sigma_sweeps.csv
    sigma50_summary.csv
    sigma50_vs_I_by_T.png
    escape_probability_curves_T*.png
    phase_R_commit_vs_escape.png
    rotations_before_commitment.png

Example
-------
python final_resonator_verification_suite.py \
  --I-min 0.320 --I-max 0.345 --I-step 0.0005 \
  --sweep-I-values 0.324 0.326 0.328 0.330 0.332 0.334 0.336 0.338 \
  --sigmas 0.002 0.004 0.006 0.008 0.010 0.012 0.014 0.016 0.018 0.020 0.024 \
  --T-values 250 500 \
  --N 500 --dt 0.005 \
  --commit-radius 1.0 --reset-radius 0.55 \
  --outdir results/final_resonator_verification
"""

from __future__ import annotations

import argparse
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    from scipy.optimize import brentq, curve_fit
except Exception as exc:  # pragma: no cover
    raise RuntimeError("This script needs scipy installed.") from exc


# -----------------------------
# Deterministic FHN utilities
# -----------------------------

def fixed_point_v_roots(I: float, a: float, b: float) -> np.ndarray:
    """Return real roots v of fixed-point equation.

    Fixed point satisfies w=(v+a)/b and
        v - v^3/3 - (v+a)/b + I = 0.
    """
    # Multiply by 3b:
    # 3 b v - b v^3 - 3(v+a) + 3 b I = 0
    # -b v^3 + (3b-3)v + (3bI-3a)=0
    coeff = np.array([-b, 0.0, 3.0 * b - 3.0, 3.0 * b * I - 3.0 * a], dtype=float)
    roots = np.roots(coeff)
    real_roots = roots[np.abs(roots.imag) < 1e-10].real
    return np.sort(real_roots)


def jacobian(v: float, eps: float, b: float) -> np.ndarray:
    return np.array([[1.0 - v * v, -1.0], [eps, -eps * b]], dtype=float)


@dataclass
class FixedPointInfo:
    I: float
    v: float
    w: float
    eig1_real: float
    eig1_imag: float
    eig2_real: float
    eig2_imag: float
    max_real: float
    trace: float
    determinant: float
    stable: bool
    spiral: bool


def select_fixed_point(I: float, a: float, b: float, eps: float) -> FixedPointInfo:
    """Select the biologically relevant fixed point.

    For the standard parameter values in this project, there is usually one
    real root in the resonator window. If there are multiple, choose the one
    with smallest maximum eigenvalue real part; this avoids accidentally
    selecting an unstable saddle when scanning more widely.
    """
    infos = []
    for v in fixed_point_v_roots(I, a, b):
        w = (v + a) / b
        J = jacobian(v, eps, b)
        eig = np.linalg.eigvals(J)
        info = FixedPointInfo(
            I=I,
            v=float(v),
            w=float(w),
            eig1_real=float(eig[0].real),
            eig1_imag=float(eig[0].imag),
            eig2_real=float(eig[1].real),
            eig2_imag=float(eig[1].imag),
            max_real=float(np.max(eig.real)),
            trace=float(np.trace(J)),
            determinant=float(np.linalg.det(J)),
            stable=bool(np.max(eig.real) < 0),
            spiral=bool(np.any(np.abs(eig.imag) > 1e-9)),
        )
        infos.append(info)
    if not infos:
        raise RuntimeError(f"No real fixed point found at I={I}")
    return sorted(infos, key=lambda x: x.max_real)[0]


def deterministic_scan(I_values: np.ndarray, a: float, b: float, eps: float) -> pd.DataFrame:
    rows = [select_fixed_point(float(I), a, b, eps).__dict__ for I in I_values]
    return pd.DataFrame(rows)


def estimate_hopf_from_scan(df: pd.DataFrame, a: float, b: float, eps: float) -> Optional[float]:
    """Estimate Hopf point by bracketing a sign change of max_real."""
    I = df["I"].to_numpy()
    y = df["max_real"].to_numpy()
    for k in range(len(I) - 1):
        if y[k] == 0:
            return float(I[k])
        if y[k] * y[k + 1] < 0:
            lo, hi = float(I[k]), float(I[k + 1])
            f = lambda x: select_fixed_point(x, a, b, eps).max_real
            return float(brentq(f, lo, hi, xtol=1e-12, rtol=1e-12, maxiter=100))
    return None


def plot_deterministic_scan(df: pd.DataFrame, hopf_I: Optional[float], outpath: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(df["I"], df["max_real"], marker="o", markersize=3)
    ax.axhline(0, linestyle="--", linewidth=1)
    if hopf_I is not None:
        ax.axvline(hopf_I, linestyle="--", linewidth=1)
        ax.text(hopf_I, ax.get_ylim()[1] * 0.85, f"Hopf ~ {hopf_I:.6f}", rotation=90,
                va="top", ha="right")
    ax.set_xlabel("I")
    ax.set_ylabel("max Re(lambda)")
    ax.set_title("Deterministic fixed-point stability scan")
    fig.tight_layout()
    fig.savefig(outpath, dpi=180)
    plt.close(fig)


# -----------------------------
# Phase / circular statistics
# -----------------------------

def local_spiral_basis(v_star: float, eps: float, b: float) -> tuple[np.ndarray, np.ndarray]:
    """Return fixed-point Jacobian eigenbasis P and inverse Pinv.

    P columns are Re(eigenvector), Im(eigenvector) for the complex eigenvalue
    with positive imaginary part. Coordinates z = Pinv @ ([v,w]-[v*,w*]).
    """
    J = jacobian(v_star, eps, b)
    vals, vecs = np.linalg.eig(J)
    idx_candidates = np.where(vals.imag > 1e-9)[0]
    if len(idx_candidates) == 0:
        # Fall back to eigenvalue with largest imaginary magnitude.
        idx = int(np.argmax(np.abs(vals.imag)))
    else:
        idx = int(idx_candidates[0])
    evec = vecs[:, idx]
    P = np.column_stack([evec.real, evec.imag]).astype(float)
    if abs(np.linalg.det(P)) < 1e-12:
        raise RuntimeError("Degenerate local spiral eigenbasis. Are you outside the spiral regime?")
    return P, np.linalg.inv(P)


def circular_R(phases: np.ndarray) -> float:
    phases = np.asarray(phases, dtype=float)
    phases = phases[np.isfinite(phases)]
    if phases.size == 0:
        return float("nan")
    return float(np.abs(np.mean(np.exp(1j * phases))))


def circular_mean(phases: np.ndarray) -> float:
    phases = np.asarray(phases, dtype=float)
    phases = phases[np.isfinite(phases)]
    if phases.size == 0:
        return float("nan")
    return float(np.angle(np.mean(np.exp(1j * phases))))


# -----------------------------
# Stochastic simulation
# -----------------------------

@dataclass
class SweepResult:
    I: float
    sigma: float
    T: float
    N: int
    dt: float
    stable: bool
    spiral: bool
    max_real: float
    omega: float
    linear_period: float
    excursion_fraction: float
    n_excursions: int
    median_escape_time: float
    mean_escape_time: float
    median_commit_time: float
    mean_commit_time: float
    median_rotations_before_commitment: float
    mean_rotations_before_commitment: float
    R_commit_phase: float
    mean_commit_phase: float
    R_escape_phase: float
    mean_escape_phase: float
    n_commit_phases: int


def simulate_one_setting(
    I: float,
    sigma: float,
    T: float,
    N: int,
    dt: float,
    a: float,
    b: float,
    eps: float,
    v_excursion: float,
    commit_radius: float,
    reset_radius: float,
    seed: int,
) -> SweepResult:
    """Vectorised Euler-Maruyama simulation for one (I, sigma, T)."""
    info = select_fixed_point(I, a, b, eps)
    if not info.spiral:
        # Still run, but phase basis may be meaningless. Raise early for safety.
        raise RuntimeError(f"I={I}: fixed point is not a spiral; cannot define local spiral phase.")

    P, Pinv = local_spiral_basis(info.v, eps, b)
    omega = float(abs(info.eig1_imag) if abs(info.eig1_imag) > 1e-12 else abs(info.eig2_imag))
    linear_period = float(2.0 * math.pi / omega) if omega > 0 else float("nan")

    rng = np.random.default_rng(seed)
    steps = int(round(T / dt))
    sqrt_dt = math.sqrt(dt)

    v = np.full(N, info.v, dtype=float)
    w = np.full(N, info.w, dtype=float)

    # Phase tracking.
    coords = Pinv @ np.vstack([v - info.v, w - info.w])
    prev_phase = np.arctan2(coords[1], coords[0])
    cum_phase = np.zeros(N, dtype=float)

    # Candidate commitment bookkeeping.
    candidate_active = np.zeros(N, dtype=bool)
    candidate_time = np.full(N, np.nan)
    candidate_phase = np.full(N, np.nan)
    candidate_rotations = np.full(N, np.nan)

    committed = np.zeros(N, dtype=bool)
    commit_time = np.full(N, np.nan)
    commit_phase = np.full(N, np.nan)
    commit_rotations = np.full(N, np.nan)

    escaped = np.zeros(N, dtype=bool)
    escape_time = np.full(N, np.nan)
    escape_phase = np.full(N, np.nan)

    prev_radius = np.zeros(N, dtype=float)

    for step in range(1, steps + 1):
        t = step * dt
        active = ~escaped
        if not np.any(active):
            break

        # Euler-Maruyama update.
        va = v[active]
        wa = w[active]
        noise = sigma * sqrt_dt * rng.standard_normal(np.sum(active))
        dv_det = va - va**3 / 3.0 - wa + I
        dw_det = eps * (va + a - b * wa)
        v[active] = va + dv_det * dt + noise
        w[active] = wa + dw_det * dt

        # Local phase/radius for all paths, including already escaped for indexing ease.
        coords = Pinv @ np.vstack([v - info.v, w - info.w])
        phase = np.arctan2(coords[1], coords[0])
        radius = np.sqrt(coords[0] ** 2 + coords[1] ** 2)

        dphi = np.angle(np.exp(1j * (phase - prev_phase)))
        cum_phase += dphi
        prev_phase = phase

        # New outward crossings of commitment radius.
        new_candidate = active & (~candidate_active) & (prev_radius < commit_radius) & (radius >= commit_radius)
        candidate_active[new_candidate] = True
        candidate_time[new_candidate] = t
        candidate_phase[new_candidate] = phase[new_candidate]
        candidate_rotations[new_candidate] = np.abs(cum_phase[new_candidate]) / (2.0 * math.pi)

        # If a candidate falls back into the local basin before excursion, clear it.
        fell_back = active & candidate_active & (radius <= reset_radius)
        candidate_active[fell_back] = False
        candidate_time[fell_back] = np.nan
        candidate_phase[fell_back] = np.nan
        candidate_rotations[fell_back] = np.nan

        # Detect large excursion.
        newly_escaped = active & (v >= v_excursion)
        if np.any(newly_escaped):
            escaped[newly_escaped] = True
            escape_time[newly_escaped] = t
            escape_phase[newly_escaped] = phase[newly_escaped]

            # Accept the most recent candidate. If a path reached v_excursion without
            # crossing commit_radius, use the threshold crossing itself as a fallback.
            has_candidate = newly_escaped & candidate_active
            no_candidate = newly_escaped & (~candidate_active)

            committed[has_candidate] = True
            commit_time[has_candidate] = candidate_time[has_candidate]
            commit_phase[has_candidate] = candidate_phase[has_candidate]
            commit_rotations[has_candidate] = candidate_rotations[has_candidate]

            committed[no_candidate] = True
            commit_time[no_candidate] = t
            commit_phase[no_candidate] = phase[no_candidate]
            commit_rotations[no_candidate] = np.abs(cum_phase[no_candidate]) / (2.0 * math.pi)

        prev_radius = radius

    escaped_idx = escaped
    committed_idx = committed & escaped

    def nanmedian(x):
        x = np.asarray(x, dtype=float)
        return float(np.nanmedian(x)) if np.any(np.isfinite(x)) else float("nan")

    def nanmean(x):
        x = np.asarray(x, dtype=float)
        return float(np.nanmean(x)) if np.any(np.isfinite(x)) else float("nan")

    return SweepResult(
        I=float(I),
        sigma=float(sigma),
        T=float(T),
        N=int(N),
        dt=float(dt),
        stable=info.stable,
        spiral=info.spiral,
        max_real=info.max_real,
        omega=omega,
        linear_period=linear_period,
        excursion_fraction=float(np.mean(escaped_idx)),
        n_excursions=int(np.sum(escaped_idx)),
        median_escape_time=nanmedian(escape_time[escaped_idx]),
        mean_escape_time=nanmean(escape_time[escaped_idx]),
        median_commit_time=nanmedian(commit_time[committed_idx]),
        mean_commit_time=nanmean(commit_time[committed_idx]),
        median_rotations_before_commitment=nanmedian(commit_rotations[committed_idx]),
        mean_rotations_before_commitment=nanmean(commit_rotations[committed_idx]),
        R_commit_phase=circular_R(commit_phase[committed_idx]),
        mean_commit_phase=circular_mean(commit_phase[committed_idx]),
        R_escape_phase=circular_R(escape_phase[escaped_idx]),
        mean_escape_phase=circular_mean(escape_phase[escaped_idx]),
        n_commit_phases=int(np.sum(committed_idx)),
    )


# -----------------------------
# Sigma50 fitting
# -----------------------------

def logistic(sigma: np.ndarray, sigma50: float, width: float) -> np.ndarray:
    width = max(float(width), 1e-12)
    return 1.0 / (1.0 + np.exp(-(sigma - sigma50) / width))


def interpolate_sigma_at_p(sigmas: np.ndarray, probs: np.ndarray, p: float) -> float:
    order = np.argsort(sigmas)
    x = np.asarray(sigmas, dtype=float)[order]
    y = np.asarray(probs, dtype=float)[order]

    # Enforce monotonicity gently for interpolation.
    y_mono = np.maximum.accumulate(y)
    if p < y_mono[0] or p > y_mono[-1]:
        return float("nan")
    return float(np.interp(p, y_mono, x))


def summarise_sigma50(all_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (T, I), g in all_df.groupby(["T", "I"]):
        g = g.sort_values("sigma")
        sig = g["sigma"].to_numpy(float)
        prob = g["excursion_fraction"].to_numpy(float)

        sigma50_interp = interpolate_sigma_at_p(sig, prob, 0.5)
        sigma10_interp = interpolate_sigma_at_p(sig, prob, 0.1)
        sigma90_interp = interpolate_sigma_at_p(sig, prob, 0.9)

        sigma50_fit = float("nan")
        width_fit = float("nan")
        try:
            if np.nanmin(prob) < 0.5 < np.nanmax(prob):
                p0 = [sigma50_interp if np.isfinite(sigma50_interp) else np.median(sig), 0.002]
                bounds = ([min(sig) - 0.02, 1e-5], [max(sig) + 0.02, 0.1])
                popt, _ = curve_fit(logistic, sig, prob, p0=p0, bounds=bounds, maxfev=10000)
                sigma50_fit, width_fit = float(popt[0]), float(popt[1])
        except Exception:
            pass

        # Interpolate diagnostics near sigma50 using observed summaries.
        target = sigma50_interp if np.isfinite(sigma50_interp) else sigma50_fit
        med_rot_at_50 = float("nan")
        med_time_at_50 = float("nan")
        if np.isfinite(target):
            med_rot_at_50 = float(np.interp(target, sig, g["median_rotations_before_commitment"].to_numpy(float)))
            med_time_at_50 = float(np.interp(target, sig, g["median_escape_time"].to_numpy(float)))

        rows.append({
            "T": T,
            "I": I,
            "sigma50_interp": sigma50_interp,
            "sigma50_fit": sigma50_fit,
            "sigma10_interp": sigma10_interp,
            "sigma90_interp": sigma90_interp,
            "transition_width_10_90": sigma90_interp - sigma10_interp if np.isfinite(sigma10_interp) and np.isfinite(sigma90_interp) else float("nan"),
            "logistic_width_fit": width_fit,
            "median_rotations_at_sigma50": med_rot_at_50,
            "median_escape_time_at_sigma50": med_time_at_50,
            "min_R_commit_phase": float(g["R_commit_phase"].min(skipna=True)),
            "median_R_commit_phase": float(g["R_commit_phase"].median(skipna=True)),
            "min_R_escape_phase": float(g["R_escape_phase"].min(skipna=True)),
            "median_R_escape_phase": float(g["R_escape_phase"].median(skipna=True)),
            "stable_all_sigmas": bool(g["stable"].all()),
        })
    return pd.DataFrame(rows).sort_values(["T", "I"])


# -----------------------------
# Plotting
# -----------------------------

def plot_probability_curves(all_df: pd.DataFrame, outdir: Path) -> None:
    for T, gT in all_df.groupby("T"):
        fig, ax = plt.subplots(figsize=(9, 5))
        for I, g in gT.groupby("I"):
            g = g.sort_values("sigma")
            ax.plot(g["sigma"], g["excursion_fraction"], marker="o", label=f"I={I:.3f}")
        ax.axhline(0.5, linestyle="--", linewidth=1)
        ax.set_xlabel("sigma")
        ax.set_ylabel("P(large excursion before T)")
        ax.set_title(f"Resonator excursion probability curves, T={T:g}")
        ax.legend()
        fig.tight_layout()
        fig.savefig(outdir / f"escape_probability_curves_T{T:g}.png", dpi=180)
        plt.close(fig)


def plot_sigma50(summary: pd.DataFrame, hopf_I: Optional[float], outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 5))
    for T, g in summary.groupby("T"):
        ax.plot(g["I"], g["sigma50_interp"], marker="o", label=f"interp, T={T:g}")
        if g["sigma50_fit"].notna().any():
            ax.plot(g["I"], g["sigma50_fit"], marker="x", linestyle="--", label=f"logistic fit, T={T:g}")
    if hopf_I is not None:
        ax.axvline(hopf_I, linestyle="--", linewidth=1, label=f"Hopf ~ {hopf_I:.5f}")
    ax.set_xlabel("I")
    ax.set_ylabel("sigma50")
    ax.set_title("Finite-time resonator failure boundary")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "sigma50_vs_I_by_T.png", dpi=180)
    plt.close(fig)


def plot_phase_R(all_df: pd.DataFrame, outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 5))
    for (T, I), g in all_df.groupby(["T", "I"]):
        g = g.sort_values("sigma")
        ax.plot(g["sigma"], g["R_commit_phase"], marker="o", label=f"commit I={I:.3f}, T={T:g}")
        ax.plot(g["sigma"], g["R_escape_phase"], marker="x", linestyle="--", label=f"escape I={I:.3f}, T={T:g}")
    ax.set_xlabel("sigma")
    ax.set_ylabel("circular concentration R")
    ax.set_title("Commitment-phase vs threshold-crossing phase concentration")
    ax.set_ylim(-0.05, 1.05)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(outdir / "phase_R_commit_vs_escape.png", dpi=180)
    plt.close(fig)


def plot_rotations(all_df: pd.DataFrame, outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 5))
    for (T, I), g in all_df.groupby(["T", "I"]):
        g = g.sort_values("sigma")
        ax.plot(g["sigma"], g["median_rotations_before_commitment"], marker="o", label=f"I={I:.3f}, T={T:g}")
    ax.set_xlabel("sigma")
    ax.set_ylabel("median rotations before commitment")
    ax.set_title("How many resonator rotations before commitment?")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "rotations_before_commitment.png", dpi=180)
    plt.close(fig)



# -----------------------------
# Parallel worker
# -----------------------------

def run_stochastic_job(job: dict) -> dict:
    """Top-level worker so multiprocessing works cleanly on macOS/Windows."""
    res = simulate_one_setting(
        I=job["I"],
        sigma=job["sigma"],
        T=job["T"],
        N=job["N"],
        dt=job["dt"],
        a=job["a"],
        b=job["b"],
        eps=job["eps"],
        v_excursion=job["v_excursion"],
        commit_radius=job["commit_radius"],
        reset_radius=job["reset_radius"],
        seed=job["seed"],
    )
    out = res.__dict__
    out["job_index"] = job["job_index"]
    return out

# -----------------------------
# CLI
# -----------------------------

def parse_float_list(xs: Optional[list[str]]) -> list[float]:
    if xs is None:
        return []
    return [float(x) for x in xs]


def main() -> None:
    parser = argparse.ArgumentParser(description="Final resonator verification suite")
    parser.add_argument("--a", type=float, default=0.7)
    parser.add_argument("--b", type=float, default=0.8)
    parser.add_argument("--eps", type=float, default=0.08)

    parser.add_argument("--I-min", type=float, default=0.320)
    parser.add_argument("--I-max", type=float, default=0.345)
    parser.add_argument("--I-step", type=float, default=0.0005)
    parser.add_argument("--sweep-I-values", nargs="*", default=None)
    parser.add_argument("--sigmas", nargs="*", default=None)
    parser.add_argument("--T-values", nargs="*", default=["250"])

    parser.add_argument("--N", type=int, default=500)
    parser.add_argument("--dt", type=float, default=0.005)
    parser.add_argument("--v-excursion", type=float, default=1.0)
    parser.add_argument("--commit-radius", type=float, default=1.0)
    parser.add_argument("--reset-radius", type=float, default=0.55)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--outdir", type=str, default="results/final_resonator_verification")
    parser.add_argument("--skip-stochastic", action="store_true")
    parser.add_argument("--allow-unstable", action="store_true",
                        help="If omitted, stochastic sweeps are skipped for I values past Hopf/unstable fixed point.")
    parser.add_argument("--workers", type=int, default=1,
                        help="Number of parallel worker processes for stochastic sweeps. Use 1 for serial.")
    parser.add_argument("--chunksize", type=int, default=1,
                        help="Task chunk size for multiprocessing. Usually leave at 1 because jobs are heavy.")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    I_scan = np.arange(args.I_min, args.I_max + 0.5 * args.I_step, args.I_step)
    det_df = deterministic_scan(I_scan, args.a, args.b, args.eps)
    det_df.to_csv(outdir / "deterministic_scan.csv", index=False)
    hopf_I = estimate_hopf_from_scan(det_df, args.a, args.b, args.eps)
    plot_deterministic_scan(det_df, hopf_I, outdir / "deterministic_scan.png")

    with open(outdir / "hopf_estimate.txt", "w") as f:
        if hopf_I is None:
            f.write("No Hopf sign change found inside deterministic scan range.\n")
        else:
            f.write(f"Estimated Hopf I: {hopf_I:.12f}\n")

    print("Deterministic scan saved.")
    print(f"Estimated Hopf I: {hopf_I}" if hopf_I is not None else "No Hopf found in scan interval.")

    if args.skip_stochastic:
        return

    sweep_I_values = parse_float_list(args.sweep_I_values)
    if not sweep_I_values:
        if hopf_I is not None:
            sweep_I_values = [float(x) for x in np.arange(max(args.I_min, hopf_I - 0.016), hopf_I - 0.001, 0.002)]
        else:
            sweep_I_values = [0.324, 0.326, 0.328, 0.330]

    sigmas = parse_float_list(args.sigmas)
    if not sigmas:
        sigmas = [0.002, 0.004, 0.006, 0.008, 0.010, 0.012, 0.014, 0.016, 0.018, 0.020, 0.024]

    T_values = parse_float_list(args.T_values)

    # Build stochastic job list first. This makes parallel execution simple and reproducible.
    jobs = []
    skipped = []
    job_index = 0
    for T in T_values:
        for I in sweep_I_values:
            info = select_fixed_point(I, args.a, args.b, args.eps)
            if (not args.allow_unstable) and (not info.stable):
                msg = f"Skipping I={I:.6f}: deterministic fixed point is not stable (max_real={info.max_real:.6g})."
                print(msg)
                skipped.append({"I": I, "T": T, "max_real": info.max_real, "reason": "unstable fixed point"})
                continue
            for sigma in sigmas:
                job_index += 1
                seed = args.seed + int(round(1_000_000 * I)) + int(round(100_000 * sigma)) + int(round(10 * T))
                jobs.append({
                    "job_index": job_index,
                    "I": float(I),
                    "sigma": float(sigma),
                    "T": float(T),
                    "N": int(args.N),
                    "dt": float(args.dt),
                    "a": float(args.a),
                    "b": float(args.b),
                    "eps": float(args.eps),
                    "v_excursion": float(args.v_excursion),
                    "commit_radius": float(args.commit_radius),
                    "reset_radius": float(args.reset_radius),
                    "seed": int(seed),
                })

    if skipped:
        pd.DataFrame(skipped).to_csv(outdir / "skipped_unstable_I_values.csv", index=False)

    rows = []
    total_jobs = len(jobs)
    workers = max(1, int(args.workers))
    print(f"Prepared {total_jobs} stochastic jobs. Using workers={workers}.")

    if workers == 1:
        for k, job in enumerate(jobs, start=1):
            print(f"[{k}/{total_jobs}] I={job['I']:.6f}, sigma={job['sigma']:.5f}, T={job['T']:g}")
            rows.append(run_stochastic_job(job))
            pd.DataFrame(rows).sort_values(["T", "I", "sigma"]).to_csv(outdir / "all_sigma_sweeps.csv", index=False)
    else:
        # cap workers to job count to avoid spawning silly numbers of processes
        workers = min(workers, total_jobs) if total_jobs else 1
        completed = 0
        with ProcessPoolExecutor(max_workers=workers) as ex:
            future_to_job = {ex.submit(run_stochastic_job, job): job for job in jobs}
            for fut in as_completed(future_to_job):
                job = future_to_job[fut]
                completed += 1
                try:
                    row = fut.result()
                    rows.append(row)
                    print(f"[{completed}/{total_jobs}] done I={job['I']:.6f}, sigma={job['sigma']:.5f}, T={job['T']:g}")
                except Exception as exc:
                    print(f"[{completed}/{total_jobs}] FAILED I={job['I']:.6f}, sigma={job['sigma']:.5f}, T={job['T']:g}: {exc}")
                # checkpoint after every completed job
                if rows:
                    pd.DataFrame(rows).sort_values(["T", "I", "sigma"]).to_csv(outdir / "all_sigma_sweeps.csv", index=False)

    all_df = pd.DataFrame(rows).sort_values(["T", "I", "sigma"]) if rows else pd.DataFrame(rows)
    if all_df.empty:
        print("No stochastic results generated.")
        return

    all_df.to_csv(outdir / "all_sigma_sweeps.csv", index=False)
    summary = summarise_sigma50(all_df)
    summary.to_csv(outdir / "sigma50_summary.csv", index=False)

    plot_probability_curves(all_df, outdir)
    plot_sigma50(summary, hopf_I, outdir)
    plot_phase_R(all_df, outdir)
    plot_rotations(all_df, outdir)

    print("\nDone.")
    print(f"Outputs saved in: {outdir}")
    print("\nSigma50 summary:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
