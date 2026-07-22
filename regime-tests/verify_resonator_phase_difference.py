#!/usr/bin/env python3
"""
verify_resonator_phase_difference.py

Numerical verification script for the FitzHugh--Nagumo resonator regime.

Goal
----
In the resonator window the deterministic system has a stable spiral fixed point.
There is no deterministic tonic limit cycle, so the useful object is not an ISI
period. Instead, this script quantifies how noise perturbs the *spiral phase*
around the stable focus before either

  (i) the oscillation damps back to the fixed point, or
 (ii) the trajectory makes a large excursion/spike.

It computes:
  - fixed point and Jacobian eigenvalues;
  - linear resonator frequency omega = Im(lambda);
  - stochastic phase around the stable spiral in the eigenbasis;
  - phase difference against the deterministic sigma=0 reference;
  - phase variance / circular concentration CONDITIONAL on not having escaped;
  - survival curve S(t) = P(no large excursion by t);
  - probability of a large excursion before T;
  - escape/commitment phase, radius, location, and number of local rotations before escape.

Default parameters target the resonator window of the project:
  a=0.7, b=0.8, eps=0.08, I=0.30.

Usage examples
--------------
  python verify_resonator_phase_difference.py

  python verify_resonator_phase_difference.py --I 0.30 --eps 0.08 --sigma 0.04 --N 500

  python verify_resonator_phase_difference.py --sweep-sigma 0.0 0.02 0.04 0.06 0.08

Outputs
-------
Creates an output directory containing:
  resonator_phase_summary.csv
  phase_timeseries_mean.csv
  phase_difference_timeseries.png
  phase_diff_final_hist.png
  circular_concentration.png
  survival_curve.png
  escape_time_hist.png
  escape_phase_hist.png
  escape_events.csv
  sample_phase_paths.png

Notes
-----
The phase is defined locally around the stable spiral. It is therefore meaningful
only while the amplitude remains away from zero and before a large excursion.
The script masks times where the local amplitude is too small.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# Model
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class FHNParams:
    I: float = 0.30
    a: float = 0.7
    b: float = 0.8
    eps: float = 0.08


def fhn_rhs(v: np.ndarray, w: np.ndarray, p: FHNParams) -> Tuple[np.ndarray, np.ndarray]:
    """Deterministic FHN vector field."""
    dv = v - (v ** 3) / 3.0 - w + p.I
    dw = p.eps * (v + p.a - p.b * w)
    return dv, dw


def jacobian_at(v: float, p: FHNParams) -> np.ndarray:
    """Jacobian at a fixed point with voltage coordinate v."""
    return np.array(
        [[1.0 - v * v, -1.0], [p.eps, -p.eps * p.b]], dtype=float
    )


def fixed_points(p: FHNParams) -> List[Tuple[float, float]]:
    """
    Solve fixed point equation:
      w = (v+a)/b and w = v - v^3/3 + I.
    So:
      -v^3/3 + v - (v+a)/b + I = 0.
    """
    # polynomial: -v^3/3 + (1 - 1/b)v + (I - a/b) = 0
    coeff = [-1.0 / 3.0, 0.0, 1.0 - 1.0 / p.b, p.I - p.a / p.b]
    roots = np.roots(coeff)
    real_roots = sorted(float(r.real) for r in roots if abs(r.imag) < 1e-9)
    return [(v, (v + p.a) / p.b) for v in real_roots]


def classify_fp(v: float, p: FHNParams) -> Dict[str, float | str]:
    J = jacobian_at(v, p)
    eig = np.linalg.eigvals(J)
    tr = float(np.trace(J))
    det = float(np.linalg.det(J))
    disc = float(tr * tr - 4.0 * det)
    if det < 0:
        kind = "saddle"
    elif tr < 0:
        kind = "stable spiral" if disc < 0 else "stable node"
    elif tr > 0:
        kind = "unstable spiral" if disc < 0 else "unstable node"
    else:
        kind = "center/Hopf threshold"
    return {
        "trace": tr,
        "det": det,
        "disc": disc,
        "eig_real_1": float(eig[0].real),
        "eig_imag_1": float(eig[0].imag),
        "eig_real_2": float(eig[1].real),
        "eig_imag_2": float(eig[1].imag),
        "kind": kind,
    }


def choose_stable_spiral_fp(p: FHNParams) -> Tuple[float, float, Dict[str, float | str]]:
    fps = fixed_points(p)
    classified = [(v, w, classify_fp(v, p)) for v, w in fps]
    spirals = [(v, w, c) for v, w, c in classified if c["kind"] == "stable spiral"]
    if spirals:
        # In the default resonator window there is one relevant stable spiral.
        return spirals[0]
    stable = [(v, w, c) for v, w, c in classified if str(c["kind"]).startswith("stable")]
    if stable:
        return stable[0]
    raise RuntimeError(
        f"No stable fixed point found for p={p}. Fixed points: {classified}"
    )


# -----------------------------------------------------------------------------
# Local resonator phase coordinate
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class PhaseBasis:
    fp: np.ndarray          # shape (2,)
    B: np.ndarray           # columns: real(eigenvector), imag(eigenvector)
    Binv: np.ndarray
    alpha: float            # decay rate Re(lambda)
    omega: float            # angular frequency abs(Im(lambda))


def make_phase_basis(v_fp: float, w_fp: float, p: FHNParams) -> PhaseBasis:
    """Build local coordinates from the complex eigenvector of the stable spiral."""
    J = jacobian_at(v_fp, p)
    eigvals, eigvecs = np.linalg.eig(J)
    idx = int(np.argmax(np.abs(eigvals.imag)))
    lam = eigvals[idx]
    if abs(lam.imag) < 1e-10:
        raise RuntimeError(
            "Chosen fixed point is not a spiral: eigenvalues are not complex. "
            "Phase coordinate is not meaningful here."
        )
    z = eigvecs[:, idx]
    e1 = np.real(z)
    e2 = np.imag(z)
    B = np.column_stack([e1, e2])
    if abs(np.linalg.det(B)) < 1e-12:
        raise RuntimeError("Degenerate eigenbasis; cannot define phase coordinates.")
    return PhaseBasis(
        fp=np.array([v_fp, w_fp], dtype=float),
        B=B,
        Binv=np.linalg.inv(B),
        alpha=float(lam.real),
        omega=float(abs(lam.imag)),
    )


def phase_and_radius(v: np.ndarray, w: np.ndarray, basis: PhaseBasis) -> Tuple[np.ndarray, np.ndarray]:
    """Return local spiral phase and local amplitude in the eigenbasis."""
    X = np.vstack([v - basis.fp[0], w - basis.fp[1]])
    C = basis.Binv @ X
    x = C[0]
    y = C[1]
    phase = np.arctan2(y, x)
    radius = np.sqrt(x * x + y * y)
    return phase, radius


def initial_condition_from_phase(basis: PhaseBasis, radius: float, phi: float) -> np.ndarray:
    local = np.array([radius * math.cos(phi), radius * math.sin(phi)])
    return basis.fp + basis.B @ local


# -----------------------------------------------------------------------------
# Simulation
# -----------------------------------------------------------------------------

def simulate_many(
    p: FHNParams,
    basis: PhaseBasis,
    sigma: float,
    N: int,
    T: float,
    dt: float,
    r0: float,
    phi0: float,
    spike_threshold: float,
    seed: int,
) -> Dict[str, np.ndarray]:
    """Vectorised Euler--Maruyama simulation for N trajectories."""
    rng = np.random.default_rng(seed)
    n_steps = int(round(T / dt))
    t = np.linspace(0.0, n_steps * dt, n_steps + 1)

    x0 = initial_condition_from_phase(basis, r0, phi0)
    v = np.full(N, x0[0], dtype=float)
    w = np.full(N, x0[1], dtype=float)

    V = np.empty((n_steps + 1, N), dtype=float)
    W = np.empty((n_steps + 1, N), dtype=float)
    V[0] = v
    W[0] = w

    spike_time = np.full(N, np.nan)
    spiked = np.zeros(N, dtype=bool)

    sqrt_dt = math.sqrt(dt)
    for k in range(1, n_steps + 1):
        dv_det, dw_det = fhn_rhs(v, w, p)
        v = v + dv_det * dt + sigma * sqrt_dt * rng.standard_normal(N)
        w = w + dw_det * dt

        newly_spiked = (~spiked) & (v >= spike_threshold)
        spike_time[newly_spiked] = t[k]
        spiked[newly_spiked] = True

        V[k] = v
        W[k] = w

    phase = np.empty_like(V)
    radius = np.empty_like(V)
    for k in range(n_steps + 1):
        ph, rr = phase_and_radius(V[k], W[k], basis)
        phase[k] = ph
        radius[k] = rr

    phase_unwrapped = np.unwrap(phase, axis=0)

    return {
        "t": t,
        "V": V,
        "W": W,
        "phase": phase,
        "phase_unwrapped": phase_unwrapped,
        "radius": radius,
        "spiked": spiked,
        "spike_time": spike_time,
    }


def deterministic_reference(
    p: FHNParams,
    basis: PhaseBasis,
    T: float,
    dt: float,
    r0: float,
    phi0: float,
    spike_threshold: float,
) -> Dict[str, np.ndarray]:
    return simulate_many(
        p=p,
        basis=basis,
        sigma=0.0,
        N=1,
        T=T,
        dt=dt,
        r0=r0,
        phi0=phi0,
        spike_threshold=spike_threshold,
        seed=123,
    )


# -----------------------------------------------------------------------------
# Analysis helpers
# -----------------------------------------------------------------------------

def circular_resultant(phases: np.ndarray) -> Tuple[float, float]:
    """Return mean angle and resultant length R in [0,1]."""
    z = np.exp(1j * phases)
    m = np.mean(z)
    return float(np.angle(m)), float(abs(m))


def summarise_phase_difference(
    sim: Dict[str, np.ndarray],
    ref: Dict[str, np.ndarray],
    min_radius: float,
) -> Dict[str, np.ndarray | float]:
    """
    Summarise local resonator phase before escape.

    Important convention:
    - phase statistics are conditional on trajectories that have not yet made a
      large excursion and whose local spiral radius is still above min_radius;
    - once a trajectory crosses the spike/excursion threshold, local spiral phase
      is treated as invalid and the trajectory is removed from R(t).
    """
    t = sim["t"]
    ph = sim["phase_unwrapped"]
    ph_wrapped_abs = sim["phase"]
    rr = sim["radius"]
    V = sim["V"]
    W = sim["W"]
    ref_phase = ref["phase_unwrapped"][:, 0]

    # Difference relative to the deterministic nonlinear reference trajectory.
    dphi = ph - ref_phase[:, None]
    dphi_wrapped = np.angle(np.exp(1j * dphi))

    n_time, n_traj = ph.shape
    escaped_by_t = np.zeros((n_time, n_traj), dtype=bool)
    active = rr > min_radius

    escape_phase = np.full(n_traj, np.nan)
    escape_phase_diff = np.full(n_traj, np.nan)
    escape_radius = np.full(n_traj, np.nan)
    escape_v = np.full(n_traj, np.nan)
    escape_w = np.full(n_traj, np.nan)
    rotations_before_escape = np.full(n_traj, np.nan)
    escape_index = np.full(n_traj, -1, dtype=int)

    # Do not trust local phase after a large excursion, but DO record the phase
    # immediately before threshold crossing. This is the phase-commitment object.
    for j, did_spike in enumerate(sim["spiked"]):
        if did_spike:
            k_spike = int(np.searchsorted(t, sim["spike_time"][j]))
            k0 = max(0, k_spike - 1)
            escape_index[j] = k_spike
            escape_phase[j] = ph_wrapped_abs[k0, j]
            escape_phase_diff[j] = dphi_wrapped[k0, j]
            escape_radius[j] = rr[k0, j]
            escape_v[j] = V[k0, j]
            escape_w[j] = W[k0, j]
            rotations_before_escape[j] = abs((ph[k0, j] - ph[0, j]) / (2.0 * math.pi))
            active[k_spike:, j] = False
            escaped_by_t[k_spike:, j] = True

    mean_wrapped = np.full(len(t), np.nan)
    circ_R = np.full(len(t), np.nan)
    circ_var = np.full(len(t), np.nan)
    std_unwrapped = np.full(len(t), np.nan)
    n_active = np.zeros(len(t), dtype=int)
    n_unescaped = np.zeros(len(t), dtype=int)
    survival = np.full(len(t), np.nan)

    for k in range(len(t)):
        mask = active[k]
        n_active[k] = int(mask.sum())
        n_unescaped[k] = int((~escaped_by_t[k]).sum())
        survival[k] = n_unescaped[k] / n_traj
        if mask.sum() >= 3:
            mean_angle, R = circular_resultant(dphi_wrapped[k, mask])
            mean_wrapped[k] = mean_angle
            circ_R[k] = R
            circ_var[k] = 1.0 - R
            std_unwrapped[k] = float(np.std(dphi[k, mask]))

    # Final active phase difference at last reliable time for each trajectory.
    final_dphi = []
    final_t = []
    for j in range(ph.shape[1]):
        idx = np.flatnonzero(active[:, j])
        if len(idx) > 0:
            k = idx[-1]
            final_dphi.append(dphi_wrapped[k, j])
            final_t.append(t[k])

    final_dphi = np.array(final_dphi, dtype=float)
    final_t = np.array(final_t, dtype=float)

    spiked = sim["spiked"]
    spike_time = sim["spike_time"]
    if np.any(spiked):
        mean_spike_time = float(np.nanmean(spike_time))
        median_spike_time = float(np.nanmedian(spike_time[spiked]))
        mean_rotations = float(np.nanmean(rotations_before_escape[spiked]))
        median_rotations = float(np.nanmedian(rotations_before_escape[spiked]))
        mean_escape_phase_R = circular_resultant(escape_phase[spiked])[1]
        mean_escape_phase_diff_R = circular_resultant(escape_phase_diff[spiked])[1]
    else:
        mean_spike_time = float("nan")
        median_spike_time = float("nan")
        mean_rotations = float("nan")
        median_rotations = float("nan")
        mean_escape_phase_R = float("nan")
        mean_escape_phase_diff_R = float("nan")

    return {
        "t": t,
        "dphi": dphi,
        "dphi_wrapped": dphi_wrapped,
        "active": active,
        "escaped_by_t": escaped_by_t,
        "survival": survival,
        "mean_wrapped": mean_wrapped,
        "circ_R": circ_R,
        "circ_var": circ_var,
        "std_unwrapped": std_unwrapped,
        "n_active": n_active,
        "n_unescaped": n_unescaped,
        "final_dphi": final_dphi,
        "final_t": final_t,
        "escape_index": escape_index,
        "escape_phase": escape_phase,
        "escape_phase_diff": escape_phase_diff,
        "escape_radius": escape_radius,
        "escape_v": escape_v,
        "escape_w": escape_w,
        "rotations_before_escape": rotations_before_escape,
        "excursion_fraction": float(np.mean(spiked)),
        "mean_spike_time": mean_spike_time,
        "median_spike_time": median_spike_time,
        "mean_rotations_before_escape": mean_rotations,
        "median_rotations_before_escape": median_rotations,
        "escape_phase_R": float(mean_escape_phase_R),
        "escape_phase_diff_R": float(mean_escape_phase_diff_R),
    }


# -----------------------------------------------------------------------------
# IO and plotting
# -----------------------------------------------------------------------------

def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def save_timeseries_csv(outdir: Path, summary: Dict[str, np.ndarray | float]) -> None:
    t = summary["t"]
    with open(outdir / "phase_timeseries_mean.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "t",
            "mean_phase_diff_wrapped",
            "circular_R",
            "circular_variance",
            "std_unwrapped",
            "n_active",
            "n_unescaped",
            "survival_no_excursion",
        ])
        for k in range(len(t)):
            writer.writerow([
                float(t[k]),
                float(summary["mean_wrapped"][k]) if np.isfinite(summary["mean_wrapped"][k]) else "",
                float(summary["circ_R"][k]) if np.isfinite(summary["circ_R"][k]) else "",
                float(summary["circ_var"][k]) if np.isfinite(summary["circ_var"][k]) else "",
                float(summary["std_unwrapped"][k]) if np.isfinite(summary["std_unwrapped"][k]) else "",
                int(summary["n_active"][k]),
                int(summary["n_unescaped"][k]),
                float(summary["survival"][k]) if np.isfinite(summary["survival"][k]) else "",
            ])



def save_escape_events_csv(
    outdir: Path,
    sim: Dict[str, np.ndarray],
    summary: Dict[str, np.ndarray | float],
) -> None:
    """Save one row per trajectory, including escape/commitment phase if it escaped."""
    with open(outdir / "escape_events.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "trajectory",
            "escaped",
            "escape_time",
            "escape_phase_abs_wrapped",
            "escape_phase_diff_wrapped",
            "escape_radius",
            "escape_v",
            "escape_w",
            "rotations_before_escape",
        ])
        n = len(sim["spiked"])
        for j in range(n):
            escaped = bool(sim["spiked"][j])
            writer.writerow([
                j,
                int(escaped),
                float(sim["spike_time"][j]) if escaped else "",
                float(summary["escape_phase"][j]) if escaped and np.isfinite(summary["escape_phase"][j]) else "",
                float(summary["escape_phase_diff"][j]) if escaped and np.isfinite(summary["escape_phase_diff"][j]) else "",
                float(summary["escape_radius"][j]) if escaped and np.isfinite(summary["escape_radius"][j]) else "",
                float(summary["escape_v"][j]) if escaped and np.isfinite(summary["escape_v"][j]) else "",
                float(summary["escape_w"][j]) if escaped and np.isfinite(summary["escape_w"][j]) else "",
                float(summary["rotations_before_escape"][j]) if escaped and np.isfinite(summary["rotations_before_escape"][j]) else "",
            ])


def plot_outputs(
    outdir: Path,
    sim: Dict[str, np.ndarray],
    ref: Dict[str, np.ndarray],
    summary: Dict[str, np.ndarray | float],
    max_paths: int = 25,
) -> None:
    t = sim["t"]
    dphi = summary["dphi"]
    active = summary["active"]

    # Phase difference time series: sample paths + circular mean.
    # These are deliberately blanked after escape, because local spiral phase is
    # no longer meaningful once the trajectory has committed to a large loop.
    plt.figure(figsize=(8, 4.5))
    n = dphi.shape[1]
    for j in range(min(max_paths, n)):
        y = np.where(active[:, j], dphi[:, j], np.nan)
        plt.plot(t, y, alpha=0.35, linewidth=0.8)
    plt.plot(t, summary["mean_wrapped"], linewidth=2.0, label="circular mean phase diff")
    plt.xlabel("t")
    plt.ylabel("phase difference vs deterministic ref")
    plt.title("Resonator phase difference paths before escape")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "phase_difference_timeseries.png", dpi=180)
    plt.close()

    # Circular concentration conditional on non-escaped, still-active paths.
    plt.figure(figsize=(8, 4.5))
    plt.plot(t, summary["circ_R"], label="R among non-escaped active paths")
    plt.plot(t, summary["circ_var"], label="1 - R")
    plt.xlabel("t")
    plt.ylabel("circular statistic")
    plt.ylim(-0.02, 1.02)
    plt.title("Conditional phase coherence before escape")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "circular_concentration.png", dpi=180)
    plt.close()

    # Survival curve: the most important plot once excursions dominate.
    plt.figure(figsize=(8, 4.5))
    plt.plot(t, summary["survival"], linewidth=2.0, label="S(t) = P(no excursion by t)")
    plt.xlabel("t")
    plt.ylabel("survival probability")
    plt.ylim(-0.02, 1.02)
    plt.title("Resonator survival curve")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "survival_curve.png", dpi=180)
    plt.close()

    # Active counts: separates 'phase lost because escaped' from 'phase lost
    # because the spiral radius decayed too close to the fixed point'.
    plt.figure(figsize=(8, 4.5))
    plt.plot(t, summary["n_unescaped"], label="not yet escaped")
    plt.plot(t, summary["n_active"], label="phase-active")
    plt.xlabel("t")
    plt.ylabel("number of trajectories")
    plt.title("Remaining trajectories for local phase statistics")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "active_counts.png", dpi=180)
    plt.close()

    # Final reliable phase difference histogram.
    final_dphi = summary["final_dphi"]
    if len(final_dphi) > 0:
        plt.figure(figsize=(7, 4.5))
        plt.hist(final_dphi, bins=40)
        plt.xlabel("final reliable wrapped phase difference")
        plt.ylabel("count")
        plt.title("Distribution of final pre-escape/pre-decay phase differences")
        plt.tight_layout()
        plt.savefig(outdir / "phase_diff_final_hist.png", dpi=180)
        plt.close()

    # Escape-time histogram.
    escaped = sim["spiked"]
    if np.any(escaped):
        spike_times = sim["spike_time"][escaped]
        plt.figure(figsize=(7, 4.5))
        plt.hist(spike_times, bins=40)
        plt.xlabel("first large-excursion time")
        plt.ylabel("count")
        plt.title("Distribution of resonator escape times")
        plt.tight_layout()
        plt.savefig(outdir / "escape_time_hist.png", dpi=180)
        plt.close()

        # Absolute local phase at escape: asks whether there are dangerous phases.
        esc_phase = summary["escape_phase"][escaped]
        esc_phase = esc_phase[np.isfinite(esc_phase)]
        if len(esc_phase) > 0:
            plt.figure(figsize=(7, 4.5))
            plt.hist(esc_phase, bins=np.linspace(-math.pi, math.pi, 41))
            plt.xlabel("local spiral phase immediately before excursion")
            plt.ylabel("count")
            plt.title("Escape phase histogram: dangerous phases?")
            plt.tight_layout()
            plt.savefig(outdir / "escape_phase_hist.png", dpi=180)
            plt.close()

        # Phase difference at escape relative to deterministic reference.
        esc_phase_diff = summary["escape_phase_diff"][escaped]
        esc_phase_diff = esc_phase_diff[np.isfinite(esc_phase_diff)]
        if len(esc_phase_diff) > 0:
            plt.figure(figsize=(7, 4.5))
            plt.hist(esc_phase_diff, bins=np.linspace(-math.pi, math.pi, 41))
            plt.xlabel("phase difference from deterministic ref immediately before excursion")
            plt.ylabel("count")
            plt.title("Relative phase at escape")
            plt.tight_layout()
            plt.savefig(outdir / "escape_phase_diff_hist.png", dpi=180)
            plt.close()

        rotations = summary["rotations_before_escape"][escaped]
        rotations = rotations[np.isfinite(rotations)]
        if len(rotations) > 0:
            plt.figure(figsize=(7, 4.5))
            plt.hist(rotations, bins=40)
            plt.xlabel("local rotations before excursion")
            plt.ylabel("count")
            plt.title("How many resonator rotations before commitment?")
            plt.tight_layout()
            plt.savefig(outdir / "rotations_before_escape_hist.png", dpi=180)
            plt.close()

    # Sample phase-plane paths.
    plt.figure(figsize=(6, 5.5))
    for j in range(min(max_paths, sim["V"].shape[1])):
        plt.plot(sim["V"][:, j], sim["W"][:, j], alpha=0.35, linewidth=0.8)
    plt.plot(ref["V"][:, 0], ref["W"][:, 0], color="black", linewidth=2.0, label="deterministic ref")
    plt.xlabel("v")
    plt.ylabel("w")
    plt.title("Sample resonator trajectories")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "sample_phase_paths.png", dpi=180)
    plt.close()


def write_summary_csv(
    outdir: Path,
    p: FHNParams,
    sigma: float,
    N: int,
    T: float,
    dt: float,
    r0: float,
    phi0: float,
    min_radius: float,
    fp_info: Tuple[float, float, Dict[str, float | str]],
    basis: PhaseBasis,
    summary: Dict[str, np.ndarray | float],
) -> None:
    v_fp, w_fp, cls = fp_info
    with open(outdir / "resonator_phase_summary.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["quantity", "value"])
        rows = [
            ("I", p.I),
            ("a", p.a),
            ("b", p.b),
            ("eps", p.eps),
            ("sigma", sigma),
            ("N", N),
            ("T", T),
            ("dt", dt),
            ("r0", r0),
            ("phi0", phi0),
            ("min_radius", min_radius),
            ("v_fp", v_fp),
            ("w_fp", w_fp),
            ("fp_kind", cls["kind"]),
            ("trace", cls["trace"]),
            ("det", cls["det"]),
            ("eig_real", basis.alpha),
            ("omega", basis.omega),
            ("linear_period_2pi_over_omega", 2.0 * math.pi / basis.omega),
            ("excursion_fraction", summary["excursion_fraction"]),
            ("mean_spike_time", summary["mean_spike_time"]),
            ("median_spike_time", summary["median_spike_time"]),
            ("mean_rotations_before_escape", summary["mean_rotations_before_escape"]),
            ("median_rotations_before_escape", summary["median_rotations_before_escape"]),
            ("escape_phase_R", summary["escape_phase_R"]),
            ("escape_phase_diff_R", summary["escape_phase_diff_R"]),
        ]
        for row in rows:
            writer.writerow(row)


def run_single(args: argparse.Namespace, sigma: float, outdir: Path) -> Dict[str, float]:
    p = FHNParams(I=args.I, a=args.a, b=args.b, eps=args.eps)
    fp_info = choose_stable_spiral_fp(p)
    v_fp, w_fp, cls = fp_info
    if cls["kind"] != "stable spiral":
        print(f"WARNING: chosen fixed point is {cls['kind']}, not stable spiral. Phase may be less meaningful.")

    basis = make_phase_basis(v_fp, w_fp, p)
    ref = deterministic_reference(
        p=p,
        basis=basis,
        T=args.T,
        dt=args.dt,
        r0=args.r0,
        phi0=args.phi0,
        spike_threshold=args.spike_threshold,
    )
    sim = simulate_many(
        p=p,
        basis=basis,
        sigma=sigma,
        N=args.N,
        T=args.T,
        dt=args.dt,
        r0=args.r0,
        phi0=args.phi0,
        spike_threshold=args.spike_threshold,
        seed=args.seed,
    )
    summary = summarise_phase_difference(sim, ref, min_radius=args.min_radius)

    write_summary_csv(
        outdir=outdir,
        p=p,
        sigma=sigma,
        N=args.N,
        T=args.T,
        dt=args.dt,
        r0=args.r0,
        phi0=args.phi0,
        min_radius=args.min_radius,
        fp_info=fp_info,
        basis=basis,
        summary=summary,
    )
    save_timeseries_csv(outdir, summary)
    save_escape_events_csv(outdir, sim, summary)
    plot_outputs(outdir, sim, ref, summary, max_paths=args.max_paths_plot)

    return {
        "sigma": sigma,
        "excursion_fraction": float(summary["excursion_fraction"]),
        "mean_spike_time": float(summary["mean_spike_time"]),
        "median_spike_time": float(summary["median_spike_time"]),
        "mean_rotations_before_escape": float(summary["mean_rotations_before_escape"]),
        "median_rotations_before_escape": float(summary["median_rotations_before_escape"]),
        "escape_phase_R": float(summary["escape_phase_R"]),
        "escape_phase_diff_R": float(summary["escape_phase_diff_R"]),
        "omega": basis.omega,
        "period_linear": 2.0 * math.pi / basis.omega,
        "final_phase_std": float(np.std(summary["final_dphi"])) if len(summary["final_dphi"]) else float("nan"),
        "final_phase_R": circular_resultant(summary["final_dphi"])[1] if len(summary["final_dphi"]) else float("nan"),
    }


def run_sigma_sweep(args: argparse.Namespace) -> None:
    root = ensure_dir(args.outdir)
    rows = []
    for sigma in args.sweep_sigma:
        label = f"sigma_{sigma:.5g}".replace(".", "p")
        outdir = ensure_dir(root / label)
        print(f"Running sigma={sigma} -> {outdir}")
        rows.append(run_single(args, sigma=sigma, outdir=outdir))

    with open(root / "sigma_sweep_summary.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    sig = np.array([r["sigma"] for r in rows], dtype=float)
    exc = np.array([r["excursion_fraction"] for r in rows], dtype=float)
    std = np.array([r["final_phase_std"] for r in rows], dtype=float)
    R = np.array([r["final_phase_R"] for r in rows], dtype=float)
    med_t = np.array([r["median_spike_time"] for r in rows], dtype=float)
    med_rot = np.array([r["median_rotations_before_escape"] for r in rows], dtype=float)

    plt.figure(figsize=(7, 4.5))
    plt.plot(sig, exc, marker="o")
    plt.xlabel("sigma")
    plt.ylabel("P(large excursion before T)")
    plt.title("Resonator excursion probability vs noise")
    plt.tight_layout()
    plt.savefig(root / "sigma_sweep_excursion_probability.png", dpi=180)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.plot(sig, std, marker="o", label="std(final phase diff)")
    plt.plot(sig, 1.0 - R, marker="o", label="1 - R(final phase diff)")
    plt.xlabel("sigma")
    plt.ylabel("phase dispersion")
    plt.title("Phase dispersion vs noise")
    plt.legend()
    plt.tight_layout()
    plt.savefig(root / "sigma_sweep_phase_dispersion.png", dpi=180)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.plot(sig, med_t, marker="o", label="median escape time")
    plt.xlabel("sigma")
    plt.ylabel("time")
    plt.title("Median resonator escape time vs noise")
    plt.tight_layout()
    plt.savefig(root / "sigma_sweep_median_escape_time.png", dpi=180)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.plot(sig, med_rot, marker="o", label="median rotations before escape")
    plt.xlabel("sigma")
    plt.ylabel("local rotations")
    plt.title("Median rotations before commitment vs noise")
    plt.tight_layout()
    plt.savefig(root / "sigma_sweep_rotations_before_escape.png", dpi=180)
    plt.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Verify resonator phase differences in stochastic FHN.")
    parser.add_argument("--I", type=float, default=0.30)
    parser.add_argument("--a", type=float, default=0.7)
    parser.add_argument("--b", type=float, default=0.8)
    parser.add_argument("--eps", type=float, default=0.08)
    parser.add_argument("--sigma", type=float, default=0.04)
    parser.add_argument("--sweep-sigma", type=float, nargs="*", default=None,
                        help="Run a sigma sweep, e.g. --sweep-sigma 0 0.02 0.04 0.06")
    parser.add_argument("--N", type=int, default=300, help="Number of trajectories.")
    parser.add_argument("--T", type=float, default=250.0, help="Simulation horizon.")
    parser.add_argument("--dt", type=float, default=0.005)
    parser.add_argument("--r0", type=float, default=0.06, help="Initial radius in local eigenbasis.")
    parser.add_argument("--phi0", type=float, default=0.0, help="Initial phase in local eigenbasis.")
    parser.add_argument("--min-radius", type=float, default=0.004,
                        help="Ignore phase once local amplitude decays below this.")
    parser.add_argument("--spike-threshold", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--outdir", type=str, default="results/resonator_phase_difference")
    parser.add_argument("--max-paths-plot", type=int, default=30)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.sweep_sigma is not None and len(args.sweep_sigma) > 0:
        run_sigma_sweep(args)
    else:
        outdir = ensure_dir(args.outdir)
        print(f"Running single resonator phase experiment -> {outdir}")
        row = run_single(args, sigma=args.sigma, outdir=outdir)
        print("Summary:")
        for k, v in row.items():
            print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
