#!/usr/bin/env python3
"""
resonance_scaling.py

Does the resonator commitment-phase "gate" behave like a real resonance?

Background
----------
test_commitment_gate.py showed that, near threshold noise, the phase at which a
trajectory makes its final outbound commitment sits a noise-dependent ~0.5 rad
*ahead* of the deterministic steepest-outward ("geometric danger") direction, and
that this lead is sharp deep inside the spiral (not just at the cubic fold). That
is the signature of a noise-selected gate rather than pure fold geometry.

If the mechanism is genuinely resonant, the lead should be predictable from the
local spiral eigenstructure -- specifically the damping-to-frequency ratio
    kappa = alpha / omega = Re(lambda) / |Im(lambda)|
of the stable focus. A weakly damped focus (small |alpha|/omega) rings many times
and should commit close to the geometric direction; a strongly damped focus should
lead it more. This script measures the lead across operating points that
independently move alpha (distance below the Hopf) and omega (via eps), and plots
lead vs alpha/omega.

Design
------
For each eps we locate the deterministic Hopf current I_H analytically and place
operating points at fixed deficits delta = I_H - I (so alpha ~ kappa_lin * delta).
Varying eps moves omega ~ sqrt(eps (1 - eps b^2)). At each point we:

  1. calibrate a noise sigma near the finite-time threshold (escape fraction ~0.5),
  2. run a stochastic ensemble,
  3. measure the escape-commitment mean phase at a mid radius (a fixed fraction of
     the threshold radius), and the deterministic geometric danger angle there,
  4. record the signed wrapped lead = escape_mean - geom_angle,
     alpha, omega, alpha/omega.

Outputs: resonance_scaling_summary.csv, lead_vs_kappa.png, lead_vs_eps_I.png.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_resonator_phase_difference import (  # noqa: E402
    FHNParams, choose_stable_spiral_fp, make_phase_basis, simulate_many,
    circular_resultant,
)
from test_commitment_gate import (  # noqa: E402
    last_outbound_phase_at_radius, geometric_danger_angle,
)


def hopf_current(a: float, b: float, eps: float) -> float:
    """Deterministic Hopf I for the lower (tonic-birth) branch: v* = -sqrt(1-b eps)."""
    arg = 1.0 - b * eps
    if arg <= 0:
        return float("nan")
    v = -math.sqrt(arg)
    w = (v + a) / b
    return w - v + v ** 3 / 3.0


def wrapped(d: float) -> float:
    return (d + math.pi) % (2.0 * math.pi) - math.pi


def measure_point(p: FHNParams, sigma: float, N: int, T: float, dt: float,
                  r0: float, mid_frac: float, seed: int):
    """Run one ensemble, return escape fraction and the mid-radius lead + eigendata."""
    v_fp, w_fp, cls = choose_stable_spiral_fp(p)
    if not str(cls["kind"]).startswith("stable"):
        return None
    basis = make_phase_basis(v_fp, w_fp, p)
    sim = simulate_many(p=p, basis=basis, sigma=sigma, N=N, T=T, dt=dt,
                        r0=r0, phi0=0.0, spike_threshold=1.0, seed=seed)
    t = sim["t"]; radius = sim["radius"]; phase = sim["phase"]
    spiked = sim["spiked"]; spike_time = sim["spike_time"]
    esc = np.where(spiked)[0]
    frac = float(spiked.mean())
    if esc.size < 15:
        return {"frac": frac, "ok": False, "alpha": basis.alpha, "omega": basis.omega}

    k_spike = {j: int(np.searchsorted(t, spike_time[j])) for j in esc}
    thr_radius = float(np.median([radius[max(0, k_spike[j] - 1), j] for j in esc]))

    # measure lead at two mid radii and average (avoid relying on a single shell)
    leads = []
    esc_R_vals = []
    for frac_m in (mid_frac, mid_frac * 1.25):
        R_m = frac_m * thr_radius
        phs = [last_outbound_phase_at_radius(radius[:, j], phase[:, j], k_spike[j], R_m)
               for j in esc]
        phs = np.array([x for x in phs if np.isfinite(x)])
        if phs.size < 10:
            continue
        esc_mean, esc_R = circular_resultant(phs)
        geom = geometric_danger_angle(basis, p, R_m, dt)
        leads.append(wrapped(esc_mean - geom))
        esc_R_vals.append(esc_R)
    if not leads:
        return {"frac": frac, "ok": False, "alpha": basis.alpha, "omega": basis.omega}

    # average the (circular) leads
    lead = float(np.angle(np.mean(np.exp(1j * np.array(leads)))))
    return {
        "frac": frac, "ok": True,
        "alpha": basis.alpha, "omega": basis.omega,
        "kappa": abs(basis.alpha) / basis.omega,
        "lead_signed": lead, "lead_mag": abs(lead),
        "escape_R": float(np.mean(esc_R_vals)), "thr_radius": thr_radius,
    }


def calibrate_sigma(p, sigmas, N, T, dt, r0, seed):
    """Pick sigma whose escape fraction is closest to 0.5 (near finite-time threshold)."""
    best = None
    for s in sigmas:
        v_fp, w_fp, cls = choose_stable_spiral_fp(p)
        basis = make_phase_basis(v_fp, w_fp, p)
        sim = simulate_many(p=p, basis=basis, sigma=s, N=N, T=T, dt=dt,
                            r0=r0, phi0=0.0, spike_threshold=1.0, seed=seed)
        frac = float(sim["spiked"].mean())
        score = abs(frac - 0.5)
        if best is None or score < best[1]:
            best = (s, score, frac)
        if frac > 0.9:   # no need to keep raising sigma
            break
    return best[0], best[2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", type=float, default=0.7)
    ap.add_argument("--b", type=float, default=0.8)
    ap.add_argument("--eps-values", type=float, nargs="*", default=[0.04, 0.08, 0.12])
    ap.add_argument("--deltas", type=float, nargs="*", default=[0.030, 0.020, 0.012, 0.006],
                    help="Deficits below the Hopf current: I = I_H - delta.")
    ap.add_argument("--sigma-grid", type=float, nargs="*",
                    default=[0.01, 0.015, 0.02, 0.03, 0.04, 0.06, 0.08])
    ap.add_argument("--N-cal", type=int, default=200)
    ap.add_argument("--N", type=int, default=500)
    ap.add_argument("--T", type=float, default=220.0)
    ap.add_argument("--dt", type=float, default=0.01)
    ap.add_argument("--r0", type=float, default=0.06)
    ap.add_argument("--mid-frac", type=float, default=0.25)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--outdir", type=str, default="results/resonance_scaling")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for eps in args.eps_values:
        I_H = hopf_current(args.a, args.b, eps)
        print(f"\n=== eps={eps}  I_Hopf={I_H:.4f} ===")
        for delta in args.deltas:
            I = I_H - delta
            p = FHNParams(I=I, a=args.a, b=args.b, eps=eps)
            sigma, cal_frac = calibrate_sigma(p, args.sigma_grid, args.N_cal,
                                              args.T, args.dt, args.r0, args.seed)
            res = measure_point(p, sigma, args.N, args.T, args.dt,
                                args.r0, args.mid_frac, args.seed)
            if res is None:
                print(f"  delta={delta:.3f} I={I:.4f}: no stable spiral, skip")
                continue
            if not res["ok"]:
                print(f"  delta={delta:.3f} I={I:.4f}: too few escapes "
                      f"(sigma={sigma}, frac={res['frac']:.2f}), skip")
                continue
            row = {"eps": eps, "I": I, "I_Hopf": I_H, "delta": delta,
                   "sigma": sigma, "escape_frac": res["frac"],
                   "alpha": res["alpha"], "omega": res["omega"],
                   "kappa": res["kappa"], "lead_signed": res["lead_signed"],
                   "lead_mag": res["lead_mag"], "escape_R": res["escape_R"],
                   "thr_radius": res["thr_radius"]}
            rows.append(row)
            print(f"  delta={delta:.3f} I={I:.4f} sigma={sigma:.3f} "
                  f"frac={res['frac']:.2f} alpha={res['alpha']:.4f} omega={res['omega']:.4f} "
                  f"kappa={res['kappa']:.4f} lead={res['lead_signed']:+.3f} rad")

    if not rows:
        print("No usable points.")
        return

    with open(outdir / "resonance_scaling_summary.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wr.writeheader(); wr.writerows(rows)

    kappa = np.array([r["kappa"] for r in rows])
    lead = np.array([r["lead_signed"] for r in rows])
    epsv = np.array([r["eps"] for r in rows])

    # lead vs kappa, colour by eps
    plt.figure(figsize=(7.5, 5))
    for e in sorted(set(epsv)):
        m = epsv == e
        order = np.argsort(kappa[m])
        plt.plot(kappa[m][order], lead[m][order], "o-", label=f"eps={e}")
    plt.xlabel(r"$\kappa = |\alpha|/\omega$  (damping / frequency)")
    plt.ylabel("commitment-phase lead over geometric direction (rad)")
    plt.title("Does the gate lead scale with spiral damping/frequency?")
    plt.axhline(0, color="gray", lw=0.8)
    plt.legend(); plt.tight_layout()
    plt.savefig(outdir / "lead_vs_kappa.png", dpi=180); plt.close()

    # lead vs delta (distance below Hopf), colour by eps
    plt.figure(figsize=(7.5, 5))
    for e in sorted(set(epsv)):
        m = epsv == e
        dd = np.array([r["delta"] for r in rows])[m]
        order = np.argsort(dd)
        plt.plot(dd[order], lead[m][order], "o-", label=f"eps={e}")
    plt.xlabel(r"$\delta = I_{Hopf} - I$  (distance below Hopf)")
    plt.ylabel("commitment-phase lead (rad)")
    plt.title("Lead vs distance below the Hopf")
    plt.axhline(0, color="gray", lw=0.8)
    plt.legend(); plt.tight_layout()
    plt.savefig(outdir / "lead_vs_eps_I.png", dpi=180); plt.close()

    # correlation report
    if len(kappa) >= 3:
        r_pear = float(np.corrcoef(kappa, lead)[0, 1])
        A = np.vstack([kappa, np.ones_like(kappa)]).T
        slope, intercept = np.linalg.lstsq(A, lead, rcond=None)[0]
        print("\n================ SCALING VERDICT ================")
        print(f"points used: {len(kappa)}")
        print(f"Pearson corr(lead, kappa) = {r_pear:+.3f}")
        print(f"linear fit: lead ~ {slope:+.3f} * kappa + {intercept:+.3f}")
        if abs(r_pear) > 0.6:
            print("=> lead tracks alpha/omega: consistent with a resonance mechanism")
            print("   predictable from the local spiral eigenstructure.")
        else:
            print("=> lead does NOT track alpha/omega cleanly: the gate is not simply")
            print("   set by the linear damping/frequency ratio (geometry or higher-order).")
        print("=================================================")
    print(f"\nOutputs -> {outdir}")


if __name__ == "__main__":
    main()
