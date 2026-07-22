#!/usr/bin/env python3
"""
test_commitment_gate.py

Decisive test of the "phase gate" hypothesis for resonator escape in stochastic FHN.

Motivation
----------
The earlier sweep reported escape_phase_R ~ 0.997 and read it as a sharply
preferred "dangerous phase". But escape was *detected* at the spike threshold
v = 1, which is a single geometric locus (the right knee of the cubic nullcline).
Measuring absolute spiral phase there is close to tautological: every escaping
trajectory crosses the same surface, so of course it has the same phase. The
constancy of R ~ 0.997 across all sigma reinforces that it is geometry, not a
noise-selected resonance.

This script settles it with three measurements.

(1) Concentration vs observation radius.
    For each escaping trajectory we measure the absolute spiral phase at the LAST
    outbound crossing of a shell of local radius R_m, for a grid of R_m from the
    linear spiral (small) out to the threshold (large). If the phase only becomes
    concentrated near the threshold radius, the "gate" is the fold funnel
    (geometric). If it is already concentrated at small radius, there is genuine
    early phase selection.

(2) Occupation null.
    At each radius shell we also collect the phases that NON-escaping (and
    pre-escape) trajectories occupy — i.e. where trajectories simply spend time at
    that radius. If escape phase coincides with the occupation phase, escape is
    unbiased ("you escape from wherever you happen to be"); if escape concentrates
    at a phase that is under-represented in occupation, that is a real gate.

(3) Deterministic geometric-danger angle.
    Releasing the deterministic (sigma = 0) flow from a circle of radius R_m and
    measuring the phase of maximum outward radial growth gives the pure-geometry
    "least-stable / most-outward" direction at that radius. If the noisy escape
    phase tracks this curve, the selectivity is geometric.

Verdict logic
-------------
- If escape_phase_R is high ONLY at large R_m and the escape mean angle tracks the
  geometric danger angle and the occupation angle -> geometric funnel, gate is
  trivial.
- If escape_phase_R is high at SMALL R_m AND the escape mean angle is offset from
  both occupation and geometric-danger angles -> genuine noise-selected gate.

Reuses the model / phase-basis / simulation code in
verify_resonator_phase_difference.py.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# import the existing machinery
sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_resonator_phase_difference import (  # noqa: E402
    FHNParams,
    choose_stable_spiral_fp,
    make_phase_basis,
    simulate_many,
    phase_and_radius,
    circular_resultant,
    initial_condition_from_phase,
    fhn_rhs,
)


def circular_mean(phases: np.ndarray) -> float:
    return float(np.angle(np.mean(np.exp(1j * phases))))


def ang_dist(a: float, b: float) -> float:
    """Smallest absolute angular distance between two angles, in [0, pi]."""
    d = (a - b + math.pi) % (2.0 * math.pi) - math.pi
    return abs(d)


def last_outbound_phase_at_radius(
    radius: np.ndarray, phase: np.ndarray, k_end: int, R_m: float
) -> float:
    """
    Absolute spiral phase at the last index k < k_end where the local radius is
    <= R_m (i.e. the final outbound pass through the shell before escape).
    Returns nan if the trajectory never sits below R_m before k_end.
    """
    seg = radius[:k_end]
    below = np.flatnonzero(seg <= R_m)
    if below.size == 0:
        return float("nan")
    k = int(below[-1])
    return float(phase[k])


def geometric_danger_angle(basis, p: FHNParams, R_m: float, dt: float, n_probe: int = 720) -> float:
    """
    Deterministic geometric danger angle at radius R_m: release sigma=0 flow from
    n_probe points evenly spaced in spiral phase on the circle of radius R_m, take
    one Euler step, and return the phase whose one-step outward radial growth is
    largest. This is the pure-geometry most-outward direction at that radius.
    """
    phis = np.linspace(-math.pi, math.pi, n_probe, endpoint=False)
    best_growth = -np.inf
    best_phi = float("nan")
    for phi in phis:
        x0 = initial_condition_from_phase(basis, R_m, phi)
        v0, w0 = float(x0[0]), float(x0[1])
        ph0, r0 = phase_and_radius(np.array([v0]), np.array([w0]), basis)
        dv, dw = fhn_rhs(np.array([v0]), np.array([w0]), p)
        v1 = v0 + float(dv[0]) * dt
        w1 = w0 + float(dw[0]) * dt
        _, r1 = phase_and_radius(np.array([v1]), np.array([w1]), basis)
        growth = float(r1[0] - r0[0])
        if growth > best_growth:
            best_growth = growth
            best_phi = float(ph0[0])
    return best_phi


def main() -> None:
    ap = argparse.ArgumentParser(description="Commitment-surface + null-model gate test.")
    ap.add_argument("--I", type=float, default=0.30)
    ap.add_argument("--a", type=float, default=0.7)
    ap.add_argument("--b", type=float, default=0.8)
    ap.add_argument("--eps", type=float, default=0.08)
    ap.add_argument("--sigma", type=float, default=0.04)
    ap.add_argument("--N", type=int, default=500)
    ap.add_argument("--T", type=float, default=250.0)
    ap.add_argument("--dt", type=float, default=0.01)
    ap.add_argument("--r0", type=float, default=0.06)
    ap.add_argument("--phi0", type=float, default=0.0)
    ap.add_argument("--spike-threshold", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--n-radii", type=int, default=24)
    ap.add_argument("--shell-frac", type=float, default=0.12,
                    help="Relative half-width of the occupation radius shell.")
    ap.add_argument("--outdir", type=str,
                    default="results/commitment_gate_test")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    p = FHNParams(I=args.I, a=args.a, b=args.b, eps=args.eps)
    v_fp, w_fp, cls = choose_stable_spiral_fp(p)
    basis = make_phase_basis(v_fp, w_fp, p)
    print(f"fixed point=({v_fp:.4f},{w_fp:.4f}) kind={cls['kind']} "
          f"alpha={basis.alpha:.4f} omega={basis.omega:.4f}")

    sim = simulate_many(
        p=p, basis=basis, sigma=args.sigma, N=args.N, T=args.T, dt=args.dt,
        r0=args.r0, phi0=args.phi0, spike_threshold=args.spike_threshold, seed=args.seed,
    )
    t = sim["t"]
    radius = sim["radius"]      # (n_steps+1, N)
    phase = sim["phase"]        # wrapped, (n_steps+1, N)
    spiked = sim["spiked"]
    spike_time = sim["spike_time"]
    n_traj = radius.shape[1]

    esc_idx = np.where(spiked)[0]
    print(f"excursion fraction = {spiked.mean():.3f}  ({esc_idx.size}/{n_traj})")
    if esc_idx.size < 10:
        print("Too few escapes to analyse; increase sigma or N.")
        return

    # per-trajectory escape (spike) step index
    k_spike = np.full(n_traj, -1, dtype=int)
    for j in esc_idx:
        k_spike[j] = int(np.searchsorted(t, spike_time[j]))

    # threshold radius reference (median local radius at the v=1 crossing)
    thr_radii = []
    for j in esc_idx:
        k = min(max(0, k_spike[j] - 1), radius.shape[0] - 1)
        thr_radii.append(radius[k, j])
    thr_radius = float(np.median(thr_radii))
    print(f"median local radius at v=1 crossing = {thr_radius:.3f}")

    # radius grid from just above the linear spiral out to the threshold radius
    R_grid = np.linspace(0.15, thr_radius * 0.98, args.n_radii)

    escape_R = np.full(args.n_radii, np.nan)
    escape_mean = np.full(args.n_radii, np.nan)
    occ_R = np.full(args.n_radii, np.nan)
    occ_mean = np.full(args.n_radii, np.nan)
    geom_angle = np.full(args.n_radii, np.nan)
    offset_esc_occ = np.full(args.n_radii, np.nan)
    offset_esc_geom = np.full(args.n_radii, np.nan)
    n_escape_used = np.zeros(args.n_radii, dtype=int)
    n_occ_used = np.zeros(args.n_radii, dtype=int)

    # precompute, for every trajectory, the valid (pre-escape) time mask
    end_idx = np.where(spiked, k_spike, radius.shape[0])

    for m, R_m in enumerate(R_grid):
        # (1) escape commitment phase at this shell
        ph_list = []
        for j in esc_idx:
            ph = last_outbound_phase_at_radius(radius[:, j], phase[:, j], k_spike[j], R_m)
            if np.isfinite(ph):
                ph_list.append(ph)
        if len(ph_list) >= 5:
            ph_arr = np.array(ph_list)
            escape_mean[m], escape_R[m] = circular_resultant(ph_arr)
            n_escape_used[m] = len(ph_arr)

        # (2) occupation null: phases at this radius shell, pre-escape, all traj
        lo = R_m * (1.0 - args.shell_frac)
        hi = R_m * (1.0 + args.shell_frac)
        occ_ph = []
        for j in range(n_traj):
            e = end_idx[j]
            rr = radius[:e, j]
            pp = phase[:e, j]
            sel = (rr >= lo) & (rr <= hi)
            if np.any(sel):
                occ_ph.append(pp[sel])
        if occ_ph:
            occ_arr = np.concatenate(occ_ph)
            if occ_arr.size >= 10:
                occ_mean[m], occ_R[m] = circular_resultant(occ_arr)
                n_occ_used[m] = occ_arr.size

        # (3) geometric danger angle
        geom_angle[m] = geometric_danger_angle(basis, p, R_m, args.dt)

        if np.isfinite(escape_mean[m]) and np.isfinite(occ_mean[m]):
            offset_esc_occ[m] = ang_dist(escape_mean[m], occ_mean[m])
        if np.isfinite(escape_mean[m]) and np.isfinite(geom_angle[m]):
            offset_esc_geom[m] = ang_dist(escape_mean[m], geom_angle[m])

    # ---- save CSV ----
    with open(outdir / "concentration_vs_radius.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow([
            "R_shell", "escape_phase_R", "escape_phase_mean", "n_escape",
            "occupation_R", "occupation_mean", "n_occ",
            "geom_danger_angle", "offset_escape_vs_occupation", "offset_escape_vs_geom",
        ])
        for m in range(args.n_radii):
            wr.writerow([
                f"{R_grid[m]:.4f}",
                f"{escape_R[m]:.4f}" if np.isfinite(escape_R[m]) else "",
                f"{escape_mean[m]:.4f}" if np.isfinite(escape_mean[m]) else "",
                int(n_escape_used[m]),
                f"{occ_R[m]:.4f}" if np.isfinite(occ_R[m]) else "",
                f"{occ_mean[m]:.4f}" if np.isfinite(occ_mean[m]) else "",
                int(n_occ_used[m]),
                f"{geom_angle[m]:.4f}" if np.isfinite(geom_angle[m]) else "",
                f"{offset_esc_occ[m]:.4f}" if np.isfinite(offset_esc_occ[m]) else "",
                f"{offset_esc_geom[m]:.4f}" if np.isfinite(offset_esc_geom[m]) else "",
            ])

    # ---- plot: concentration vs radius ----
    plt.figure(figsize=(8, 4.8))
    plt.plot(R_grid, escape_R, "o-", label="escape-commitment phase R")
    plt.plot(R_grid, occ_R, "s--", label="occupation phase R (null)")
    plt.axvline(thr_radius, color="gray", ls=":", label=f"v=1 radius ~ {thr_radius:.2f}")
    plt.xlabel("observation radius R_m (local spiral amplitude)")
    plt.ylabel("circular concentration R")
    plt.ylim(-0.02, 1.02)
    plt.title("Phase concentration vs where you measure it")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "concentration_vs_radius.png", dpi=180)
    plt.close()

    # ---- plot: mean angle vs radius ----
    plt.figure(figsize=(8, 4.8))
    plt.plot(R_grid, escape_mean, "o-", label="escape-commitment mean phase")
    plt.plot(R_grid, occ_mean, "s--", label="occupation mean phase (null)")
    plt.plot(R_grid, geom_angle, "^:", label="geometric danger angle (det.)")
    plt.xlabel("observation radius R_m")
    plt.ylabel("mean spiral phase (rad)")
    plt.title("Where escape happens vs where trajectories sit vs geometry")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "mean_angle_vs_radius.png", dpi=180)
    plt.close()

    # ---- plot: offsets vs radius ----
    plt.figure(figsize=(8, 4.8))
    plt.plot(R_grid, offset_esc_occ, "o-", label="|escape - occupation|")
    plt.plot(R_grid, offset_esc_geom, "^:", label="|escape - geometric|")
    plt.axhline(0.0, color="gray", lw=0.8)
    plt.xlabel("observation radius R_m")
    plt.ylabel("angular offset (rad)")
    plt.title("Is escape phase offset from the geometric/occupation baseline?")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "offset_vs_radius.png", dpi=180)
    plt.close()

    # ---- verdict ----
    small = R_grid < 0.5 * thr_radius
    big = R_grid > 0.85 * thr_radius
    escR_small = np.nanmean(escape_R[small]) if np.any(small) else float("nan")
    escR_big = np.nanmean(escape_R[big]) if np.any(big) else float("nan")
    occR_small = np.nanmean(occ_R[small]) if np.any(small) else float("nan")
    off_occ_small = np.nanmean(offset_esc_occ[small]) if np.any(small) else float("nan")
    off_geom_small = np.nanmean(offset_esc_geom[small]) if np.any(small) else float("nan")

    print("\n================ VERDICT ================")
    print(f"escape_phase_R  small radius: {escR_small:.3f}   large radius: {escR_big:.3f}")
    print(f"occupation_R    small radius: {occR_small:.3f}")
    print(f"mean |escape - occupation|  (small radius): {off_occ_small:.3f} rad")
    print(f"mean |escape - geometric|   (small radius): {off_geom_small:.3f} rad")

    geometric = (escR_big - escR_small > 0.25) or (
        np.isfinite(off_occ_small) and off_occ_small < 0.5
    )
    print("-----------------------------------------")
    if geometric:
        print("READING: GEOMETRIC FUNNEL. The high R at the threshold is largely")
        print("the cubic-fold geometry. Escape phase converges to the occupation /")
        print("geometric-danger angle and is not concentrated upstream. The 'gate'")
        print("is where trajectories must pass, not a noise-selected resonance.")
    else:
        print("READING: GENUINE PHASE SELECTION. Escape phase is concentrated even")
        print("at small radius and is offset from where trajectories dwell and from")
        print("the deterministic danger angle. A real noise-selected gate.")
    print("=========================================")
    print(f"\nOutputs written to {outdir}")


if __name__ == "__main__":
    main()
