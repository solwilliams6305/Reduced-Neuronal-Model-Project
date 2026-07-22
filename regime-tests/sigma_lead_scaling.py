#!/usr/bin/env python3
"""
sigma_lead_scaling.py

Does the measured commitment lead climb toward the sigma->0 instanton value
(-1.06 rad, from instanton_exit.py gMAM) as the noise shrinks?

For the canonical resonator point (I=0.30, a=0.7, b=0.8, eps=0.08) we run
stochastic ensembles at a ladder of decreasing sigma, and at each one measure the
escape-commitment lead over the geometric-danger direction at the mid shell
(0.25 * threshold radius), exactly as in test_commitment_gate.py. Smaller sigma
escapes are rarer, so we lengthen T to keep enough escapers.

Output: sigma_lead_scaling.csv, sigma_lead_scaling.png.
"""
from __future__ import annotations
import argparse, csv, math
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_resonator_phase_difference import (
    FHNParams, choose_stable_spiral_fp, make_phase_basis, simulate_many,
    circular_resultant,
)
from test_commitment_gate import last_outbound_phase_at_radius, geometric_danger_angle


def wrapped(d): return (d + math.pi) % (2.0 * math.pi) - math.pi


def measure_lead(p, basis, sigma, N, T, dt, r0, mid_frac, seed):
    sim = simulate_many(p=p, basis=basis, sigma=sigma, N=N, T=T, dt=dt,
                        r0=r0, phi0=0.0, spike_threshold=1.0, seed=seed)
    t = sim["t"]; radius = sim["radius"]; phase = sim["phase"]
    spiked = sim["spiked"]; spike_time = sim["spike_time"]
    esc = np.where(spiked)[0]; frac = float(spiked.mean())
    if esc.size < 15:
        return {"frac": frac, "n_esc": int(esc.size), "ok": False}
    k_spike = {j: int(np.searchsorted(t, spike_time[j])) for j in esc}
    thr_radius = float(np.median([radius[max(0, k_spike[j]-1), j] for j in esc]))
    leads, Rs = [], []
    for fm in (mid_frac, mid_frac*1.25):
        R_m = fm * thr_radius
        phs = [last_outbound_phase_at_radius(radius[:, j], phase[:, j], k_spike[j], R_m) for j in esc]
        phs = np.array([x for x in phs if np.isfinite(x)])
        if phs.size < 10: continue
        esc_mean, R = circular_resultant(phs)
        geom = geometric_danger_angle(basis, p, R_m, dt)
        leads.append(wrapped(esc_mean - geom)); Rs.append(R)
    if not leads:
        return {"frac": frac, "n_esc": int(esc.size), "ok": False}
    lead = float(np.angle(np.mean(np.exp(1j*np.array(leads)))))
    return {"frac": frac, "n_esc": int(esc.size), "ok": True,
            "lead": lead, "R": float(np.mean(Rs)), "thr_radius": thr_radius}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--I", type=float, default=0.30)
    ap.add_argument("--a", type=float, default=0.7)
    ap.add_argument("--b", type=float, default=0.8)
    ap.add_argument("--eps", type=float, default=0.08)
    ap.add_argument("--sigmas", type=float, nargs="*",
                    default=[0.06, 0.05, 0.04, 0.03, 0.025])
    ap.add_argument("--N", type=int, default=400)
    ap.add_argument("--T", type=float, default=400.0)
    ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--r0", type=float, default=0.06)
    ap.add_argument("--mid-frac", type=float, default=0.25)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--instanton-lead", type=float, default=-1.06)
    ap.add_argument("--outdir", type=str, default="results/sigma_lead_scaling")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    p = FHNParams(I=args.I, a=args.a, b=args.b, eps=args.eps)
    v_fp, w_fp, cls = choose_stable_spiral_fp(p)
    basis = make_phase_basis(v_fp, w_fp, p)
    print(f"canonical point I={p.I} eps={p.eps}  alpha={basis.alpha:.4f} omega={basis.omega:.4f}")
    print(f"instanton (sigma->0) lead target = {args.instanton_lead:+.3f} rad\n")

    rows = []
    print("  sigma   n_esc  frac    lead(rad)   R")
    for s in sorted(args.sigmas, reverse=True):
        res = measure_lead(p, basis, s, args.N, args.T, args.dt, args.r0, args.mid_frac, args.seed)
        if not res["ok"]:
            print(f"  {s:.3f}   {res['n_esc']:4d}  {res['frac']:.2f}   too few escapes, skip")
            continue
        rows.append({"sigma": s, "n_esc": res["n_esc"], "frac": res["frac"],
                     "lead": res["lead"], "R": res["R"]})
        print(f"  {s:.3f}   {res['n_esc']:4d}  {res['frac']:.2f}   {res['lead']:+.3f}    {res['R']:.3f}")

    if len(rows) < 2:
        print("\nNot enough points."); return
    with open(outdir/"sigma_lead_scaling.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)

    sig = np.array([r["sigma"] for r in rows]); lead = np.array([r["lead"] for r in rows])
    # trend: is lead becoming more negative (toward instanton) as sigma decreases?
    slope = float(np.polyfit(sig, lead, 1)[0])
    print("\n================ SIGMA->0 VERDICT ================")
    print(f"points: {len(rows)}")
    print(f"lead vs sigma slope = {slope:+.3f} rad per unit sigma")
    print(f"smallest-sigma lead = {lead[np.argmin(sig)]:+.3f} rad  (sigma={sig.min():.3f})")
    print(f"instanton target    = {args.instanton_lead:+.3f} rad")
    if slope > 0.5:
        print("=> lead grows (more negative) as sigma shrinks: trending toward the")
        print("   instanton ceiling. Finite-noise correction confirmed.")
    else:
        print("=> lead does not clearly trend toward the instanton as sigma shrinks.")
    print("==================================================")

    plt.figure(figsize=(7.5, 5))
    plt.plot(sig, lead, "o-", label="measured commitment lead")
    plt.axhline(args.instanton_lead, color="r", ls="--",
                label=f"instanton (sigma->0) = {args.instanton_lead:+.2f}")
    plt.axhline(0, color="gray", lw=0.8)
    plt.xlabel("noise sigma"); plt.ylabel("commitment lead over geometric direction (rad)")
    plt.title("Does the lead climb toward the instanton as noise shrinks?")
    plt.gca().invert_xaxis(); plt.legend(); plt.tight_layout()
    plt.savefig(outdir/"sigma_lead_scaling.png", dpi=170); plt.close()
    print(f"\nOutputs -> {outdir}")


if __name__ == "__main__":
    main()
