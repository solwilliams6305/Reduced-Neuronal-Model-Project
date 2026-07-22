#!/usr/bin/env python3
"""
canard_sigma_crit.py

Validate the canard-strip critical-noise scaling  sigma_crit(eps) ~ eps^{3/4}
(Berglund-Gentz), versus the fold value eps^{1/2}.

Method (dynamic passage, BG-standard): slowly ramp the current I(t) upward through
the lower Hopf I_H1(eps); the deterministic trajectory tracks the repelling slow
manifold past I_H1 by a bifurcation delay D0 = <I_spike> - I_H1. Noise erodes the
delay. sigma_crit is defined (user choice) as the noise level where the mean delay
falls to HALF its deterministic value:  D(sigma_crit) = 1/2 D0.

Run once per eps (cheap, resumable); rows appended to CSV. Then run with --fit to
interpolate sigma_crit per eps and fit the exponent.
"""
from __future__ import annotations
import argparse, csv, sys, types
from pathlib import Path
import numpy as np

# --- scipy.optimize.brentq shim (sandbox has no scipy) ---
def _brentq(f, lo, hi, xtol=1e-12, maxiter=200):
    flo, fhi = f(lo), f(hi)
    if flo == 0: return lo
    if fhi == 0: return hi
    if flo*fhi > 0: raise ValueError("no bracket")
    for _ in range(maxiter):
        mid = 0.5*(lo+hi); fm = f(mid)
        if abs(hi-lo) < xtol or fm == 0: return mid
        if flo*fm < 0: hi, fhi = mid, fm
        else: lo, flo = mid, fm
    return 0.5*(lo+hi)
_sp = types.ModuleType("scipy"); _op = types.ModuleType("scipy.optimize")
_op.brentq = _brentq; _sp.optimize = _op
sys.modules.setdefault("scipy", _sp); sys.modules.setdefault("scipy.optimize", _op)
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from kernel import FHN2D

A, B = 0.7, 0.8


def mean_delay(eps, sigma, ramp, I_start_below, I_overshoot, dt, n, seed):
    """Ensemble-mean bifurcation delay D = <I_spike> - I_H1 for a slow I-ramp."""
    m = FHN2D(I=0.0, a=A, b=B)
    I_H1 = m.I_hopf_lower_at(eps)
    I0 = I_H1 - I_start_below
    m0 = FHN2D(I=I0, a=A, b=B)
    rng = np.random.default_rng(seed)
    v = np.full(n, m0.V_FP); w = np.full(n, m0.W_FP)
    I_spike = np.full(n, np.nan); fired = np.zeros(n, bool)
    sdt = np.sqrt(dt); t = 0.0; I = I0
    I_max = I_H1 + I_overshoot
    max_steps = int((I_max - I0) / ramp / dt) + 10
    for _ in range(max_steps):
        noise = rng.standard_normal(n) if sigma > 0 else 0.0
        v = v + (v - v**3/3.0 - w + I)*dt + sigma*sdt*noise
        w = w + eps*(v + A - B*w)*dt
        t += dt; I = I0 + ramp*t
        newly = (~fired) & (v >= 1.0)
        I_spike[newly] = I; fired |= newly
        if fired.all() or I > I_max:
            break
    return I_H1, float(np.nanmean(I_spike)), float(fired.mean())


def run_eps(args):
    out = Path(args.csv); out.parent.mkdir(parents=True, exist_ok=True)
    new = not out.exists()
    fh = open(out, "a", newline=""); wr = csv.writer(fh)
    if new:
        wr.writerow(["eps", "ramp", "sigma", "I_H1", "I_spike", "delay", "frac_fired"])
    eps = args.eps; ramp = args.ramp_frac * eps
    print(f"eps={eps}  ramp={ramp:.5f} (={args.ramp_frac}*eps)")
    print("  sigma     delay      frac")
    for sg in args.sigmas:
        I_H1, I_sp, frac = mean_delay(eps, sg, ramp, args.I_start_below,
                                      args.I_overshoot, args.dt, args.n, args.seed)
        wr.writerow([eps, ramp, sg, I_H1, I_sp, I_sp - I_H1, frac]); fh.flush()
        print(f"  {sg:.3f}   {I_sp-I_H1:+.5f}   {frac:.2f}")
    fh.close(); print(f"-> appended to {out}")


def fit(args):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    rows = {}
    with open(args.csv) as f:
        for r in csv.DictReader(f):
            rows.setdefault(float(r["eps"]), []).append((float(r["sigma"]), float(r["delay"])))
    eps_list, sc_list = [], []
    print("  eps      D0(det)    sigma_crit (delay=D0/2)")
    for eps in sorted(rows):
        pts = sorted(rows[eps]); sg = np.array([p[0] for p in pts]); D = np.array([p[1] for p in pts])
        D0 = D[sg == 0][0] if (sg == 0).any() else D[0]
        target = 0.5 * D0
        # delay decreases with sigma; find first crossing below target
        below = np.where(D <= target)[0]
        if len(below) == 0 or D0 <= 0:
            print(f"  {eps:.3f}   {D0:.5f}   not reached (need larger sigma)")
            continue
        k = below[0]
        if k == 0:
            sc = sg[0]
        else:
            # linear interp in sigma
            sc = sg[k-1] + (target - D[k-1]) * (sg[k] - sg[k-1]) / (D[k] - D[k-1])
        eps_list.append(eps); sc_list.append(sc)
        print(f"  {eps:.3f}   {D0:.5f}   {sc:.4f}")
    if len(eps_list) < 3:
        print("\nNeed >=3 eps with a crossing to fit. Extend sigma grid."); return
    e = np.array(eps_list); s = np.array(sc_list)
    p, lnC = np.polyfit(np.log(e), np.log(s), 1)
    print(f"\n================ CANARD-STRIP VERDICT ================")
    print(f"  sigma_crit ~ C * eps^p")
    print(f"  fitted exponent p = {p:.3f}")
    print(f"  Berglund-Gentz canard prediction p = 0.75")
    print(f"  fold (excitable)   prediction p = 0.50")
    print(f"  => closer to {'3/4 (canard)' if abs(p-0.75)<abs(p-0.5) else '1/2 (fold)'}")
    print(f"======================================================")
    plt.figure(figsize=(7,5.3))
    plt.loglog(e, s, "o", ms=8, label=f"measured (p={p:.2f})")
    xs = np.linspace(e.min(), e.max(), 50)
    plt.loglog(xs, s[0]*(xs/e[0])**0.75, "r--", label="eps^{3/4} (canard, BG)")
    plt.loglog(xs, s[0]*(xs/e[0])**0.5,  "b:", label="eps^{1/2} (fold)")
    plt.xlabel("eps"); plt.ylabel("sigma_crit (delay halves)")
    plt.title("Canard-strip critical noise: 3/4 vs 1/2"); plt.legend()
    outp = Path(args.csv).parent / "canard_sigma_crit.png"
    plt.tight_layout(); plt.savefig(outp, dpi=160); plt.close()
    print(f"-> {outp}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eps", type=float, default=0.02)
    ap.add_argument("--sigmas", type=float, nargs="+",
                    default=[0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4])
    ap.add_argument("--ramp-frac", type=float, default=0.5)
    ap.add_argument("--I-start-below", type=float, default=0.04)
    ap.add_argument("--I-overshoot", type=float, default=0.20)
    ap.add_argument("--dt", type=float, default=2.5e-3)
    ap.add_argument("--n", type=int, default=250)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--csv", type=str, default="results/canard_sigma_crit/delays.csv")
    ap.add_argument("--a", type=float, default=0.7)
    ap.add_argument("--b", type=float, default=0.8)
    ap.add_argument("--fit", action="store_true")
    args = ap.parse_args()
    global A, B
    A, B = args.a, args.b
    if args.fit: fit(args)
    else: run_eps(args)


if __name__ == "__main__":
    main()
