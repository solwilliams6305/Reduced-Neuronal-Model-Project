#!/usr/bin/env python3
"""
canard_ramp_extrap.py

Pin the canard-strip exponent by removing the protocol bias. The dynamic-passage
sigma_crit(eps) is ramp-dependent (0.94 at ramp=0.5eps, 1.60 at 0.25eps) because the
observable convolves the canard strip width (the eps^{3/4} physics we want) with the
ramp-driven bifurcation delay D0(eps, ramp). The clean exponent lives in the autonomous
ramp->0 limit.

Method: for each eps, measure sigma_crit (delay halves) at several ramp fractions, then
extrapolate sigma_crit(ramp_frac) -> ramp_frac=0 (linear intercept). Fit the resulting
ramp-free sigma_crit*(eps) ~ C eps^p and compare p to 3/4 (canard) vs 1/2 (fold).

Run once per eps (appends rows), then --fit.
"""
from __future__ import annotations
import argparse, csv, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import canard_sigma_crit as csc   # sets up scipy shim + kernel import; provides mean_delay


def sigma_crit_for(eps, ramp_frac, sigmas, I_start_below, I_overshoot, dt, n, seed):
    """Return (D0, sigma_crit) where mean delay falls to D0/2; nan if not bracketed."""
    ramp = ramp_frac * eps
    sg = np.array(sorted(sigmas))
    D = np.empty_like(sg)
    for i, s in enumerate(sg):
        _, I_sp, _ = csc.mean_delay(eps, s, ramp, I_start_below, I_overshoot, dt, n, seed)
        D[i] = I_sp  # I_spike; convert to delay below using I_H1
    # I_H1 identical across sigma; recover delay via first call's structure:
    I_H1, _, _ = csc.mean_delay(eps, 0.0, ramp, I_start_below, I_overshoot, dt, n, seed)
    delay = D - I_H1
    D0 = delay[sg == 0][0] if (sg == 0).any() else delay[0]
    target = 0.5 * D0
    below = np.where(delay <= target)[0]
    if D0 <= 0 or len(below) == 0:
        return D0, float("nan")
    k = below[0]
    if k == 0:
        return D0, float(sg[0])
    sc = sg[k-1] + (target - delay[k-1]) * (sg[k] - sg[k-1]) / (delay[k] - delay[k-1])
    return D0, float(sc)


def run(args):
    csc.A, csc.B = args.a, args.b
    out = Path(args.csv); out.parent.mkdir(parents=True, exist_ok=True)
    new = not out.exists()
    fh = open(out, "a", newline=""); wr = csv.writer(fh)
    if new:
        wr.writerow(["eps", "ramp_frac", "D0", "sigma_crit"])
    print(f"eps={args.eps}  (a,b)=({args.a},{args.b})")
    print("  ramp_frac   D0        sigma_crit")
    for rf in args.ramp_fracs:
        D0, sc = sigma_crit_for(args.eps, rf, args.sigmas, args.I_start_below,
                                args.I_overshoot, args.dt, args.n, args.seed)
        wr.writerow([args.eps, rf, D0, sc]); fh.flush()
        print(f"  {rf:.3f}      {D0:.5f}   {sc:.4f}")
    fh.close(); print(f"-> {out}")


def fit(args):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    rows = {}
    with open(args.csv) as f:
        for r in csv.DictReader(f):
            sc = float(r["sigma_crit"])
            if np.isnan(sc):
                continue
            rows.setdefault(float(r["eps"]), []).append((float(r["ramp_frac"]), sc))
    eps_list, sc0_list = [], []
    print("  eps     ramp->0 sigma_crit*   (from N ramp points)")
    for eps in sorted(rows):
        pts = sorted(rows[eps])
        rf = np.array([p[0] for p in pts]); sc = np.array([p[1] for p in pts])
        if len(rf) < 2:
            print(f"  {eps:.3f}  need >=2 ramp points"); continue
        # linear extrapolation sigma_crit(rf) -> rf=0
        slope, intercept = np.polyfit(rf, sc, 1)
        eps_list.append(eps); sc0_list.append(intercept)
        print(f"  {eps:.3f}  {intercept:.4f}             ({len(rf)} pts, dsc/drf={slope:+.3f})")
    if len(eps_list) < 3:
        print("\nNeed >=3 eps to fit exponent."); return
    e = np.array(eps_list); s = np.array(sc0_list)
    ok = s > 0
    p, lnC = np.polyfit(np.log(e[ok]), np.log(s[ok]), 1)
    print(f"\n========= RAMP-FREE CANARD EXPONENT =========")
    print(f"  sigma_crit* ~ C eps^p   (ramp -> 0 extrapolated)")
    print(f"  fitted p = {p:.3f}")
    print(f"  canard (BG) = 0.75   fold = 0.50")
    print(f"  => {'3/4 (canard)' if abs(p-0.75)<abs(p-0.5) else '1/2 (fold)'}")
    print(f"=============================================")
    plt.figure(figsize=(7,5.3))
    plt.loglog(e[ok], s[ok], "o", ms=8, label=f"ramp->0 extrapolated (p={p:.2f})")
    xs = np.linspace(e.min(), e.max(), 50)
    plt.loglog(xs, s[ok][0]*(xs/e[ok][0])**0.75, "r--", label="eps^{3/4} (canard)")
    plt.loglog(xs, s[ok][0]*(xs/e[ok][0])**0.5,  "b:", label="eps^{1/2} (fold)")
    plt.xlabel("eps"); plt.ylabel("ramp-free sigma_crit*")
    plt.title("Canard exponent, protocol bias removed (ramp->0)"); plt.legend()
    outp = Path(args.csv).parent / "canard_ramp_extrap.png"
    plt.tight_layout(); plt.savefig(outp, dpi=160); plt.close()
    print(f"-> {outp}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eps", type=float, default=0.02)
    ap.add_argument("--ramp-fracs", type=float, nargs="+", default=[0.5, 0.35, 0.25, 0.15])
    ap.add_argument("--sigmas", type=float, nargs="+",
                    default=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9, 1.1])
    ap.add_argument("--a", type=float, default=0.7)
    ap.add_argument("--b", type=float, default=0.8)
    ap.add_argument("--I-start-below", type=float, default=0.04)
    ap.add_argument("--I-overshoot", type=float, default=0.5)
    ap.add_argument("--dt", type=float, default=2.5e-3)
    ap.add_argument("--n", type=int, default=140)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--csv", type=str, default="results/canard_sigma_crit/ramp_extrap.csv")
    ap.add_argument("--fit", action="store_true")
    args = ap.parse_args()
    if args.fit: fit(args)
    else: run(args)


if __name__ == "__main__":
    main()
