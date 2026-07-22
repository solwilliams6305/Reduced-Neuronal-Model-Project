"""
Numerical analyst B (independent): real-pole-subtracted Borel phase extraction.

Model (EXACTLY as _borel_analysis.darboux_two_sing; order index k=n+1):
  v_n = 2 C zeta^{-(n+1)} Gamma(n+1+alpha) cos((n+1) theta - phi)   [complex pair]
      + Cr zr^{-(n+1)} Gamma(n+1+ar)                                [real instanton, PINNED zr, ar]

Procedure: profile theta over 38..60 deg; at each fixed (zr, ar, theta) fit
remaining params (C, zeta, phi, alpha, Cr) on RELATIVE residuals. Multi-start.
Also alpha=0-pinned variant and no-real-pole baselines.
"""
import numpy as np
from math import pi
from scipy.special import gamma as Gamma
from scipy.optimize import least_squares

V_BASE = [0.134, 0.111, 0.104, -0.030, -0.451, -1.19, -1.90]

def model(C, zeta, theta, phi, alpha, Cr, zr, ar, n):
    pair = 2*C * zeta**(-(n+1.0)) * Gamma(n+1.0+alpha) * np.cos((n+1.0)*theta - phi)
    real = Cr * zr**(-(n+1.0)) * Gamma(n+1.0+ar)
    return pair + real

# ---- fit remaining params at FIXED theta, zr, ar ----
# free params order: [C, zeta, phi, alpha, Cr]   (alpha dropped if alpha_fixed given)
def fit_fixed(v, theta, zr, ar, alpha_fixed=None, use_real=True):
    ns = np.arange(len(v), dtype=float)
    y = np.asarray(v, float)
    scale = np.abs(y) + 1e-3

    def unpack(p):
        p = list(p)
        Cr = p.pop() if use_real else 0.0
        if alpha_fixed is None:
            C, zeta, phi, alpha = p
        else:
            C, zeta, phi = p
            alpha = alpha_fixed
        return C, zeta, phi, alpha, Cr

    def resid(p):
        C, zeta, phi, alpha, Cr = unpack(p)
        return (model(C, zeta, theta, phi, alpha, Cr, zr, ar, ns) - y) / scale

    # bounds
    if alpha_fixed is None:
        lo = [1e-6, 0.3, -2*pi, -1.4, -10.0]
        hi = [10.0, 5.0,  2*pi,  1.4,  10.0]
    else:
        lo = [1e-6, 0.3, -2*pi, -10.0]
        hi = [10.0, 5.0,  2*pi,  10.0]
    if not use_real:  # drop Cr bound slot
        lo = lo[:-1]; hi = hi[:-1]

    best = None
    for phi0 in [0.3, 1.0, 1.8, 2.6]:
        for z0 in [1.5, 2.0, 2.5]:
            for cr0 in ([0.05, -0.05] if use_real else [0.0]):
                if alpha_fixed is None:
                    p0 = [0.1, z0, phi0, 0.0]
                else:
                    p0 = [0.1, z0, phi0]
                if use_real:
                    p0 = p0 + [cr0]
                p0 = np.array(p0, float)
                p0 = np.clip(p0, np.array(lo)+1e-9, np.array(hi)-1e-9)
                try:
                    r = least_squares(resid, p0, bounds=(lo, hi),
                                      xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=4000)
                except ValueError:
                    continue
                if best is None or r.cost < best.cost:
                    best = r
    p = best.x
    C, zeta, phi, alpha, Cr = unpack(p)
    maxres = np.max(np.abs(best.fun))
    # railing check
    lo_a = np.array(lo); hi_a = np.array(hi)
    span = hi_a - lo_a
    railed = np.any(np.abs(p - lo_a) < 0.01*span) or np.any(np.abs(p - hi_a) < 0.01*span)
    return dict(C=C, zeta=zeta, phi=phi, alpha=alpha, Cr=Cr, maxres=maxres, railed=bool(railed))

# ---- profile over theta ----
THETAS = np.arange(38, 61, 1)  # deg

def profile(v, zr, ar, alpha_fixed=None, use_real=True):
    recs = []
    for thd in THETAS:
        th = np.radians(thd)
        rec = fit_fixed(v, th, zr, ar, alpha_fixed=alpha_fixed, use_real=use_real)
        rec['theta'] = thd
        recs.append(rec)
    res = np.array([r['maxres'] for r in recs])
    imin = int(np.argmin(res))
    rmin = res[imin]
    argmin = THETAS[imin]
    band = THETAS[res <= 1.25*rmin]
    band_lo, band_hi = int(band.min()), int(band.max())
    best = recs[imin]
    incl45 = (band_lo <= 45 <= band_hi)
    return dict(argmin=int(argmin), rmin=float(rmin), band=(band_lo, band_hi),
                incl45=bool(incl45), Cr=best['Cr'], alpha=best['alpha'],
                railed=best['railed'], res=res.tolist(), recs=recs)

def fmt_row(label, p):
    r45 = "yes" if p['incl45'] else "no"
    return (f"{label:26s} argmin={p['argmin']:2d}  band=[{p['band'][0]:2d},{p['band'][1]:2d}]  "
            f"45in={r45:3s}  Cr={p['Cr']:+.4f}  alpha={p['alpha']:+.3f}  "
            f"rmin={p['rmin']:.3f}  railed={'Y' if p['railed'] else 'n'}")

if __name__ == "__main__":
    import sys
    v = V_BASE

    ZR_GRID = [1.4, 1.6, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6]
    AR_GRID = [-0.5, 0.0, 0.5]

    print("="*90)
    print("BASELINE (no real pole)")
    print("="*90)
    b_af = profile(v, 1.5, 0.0, alpha_fixed=None, use_real=False)
    b_a0 = profile(v, 1.5, 0.0, alpha_fixed=0.0, use_real=False)
    print(fmt_row("baseline alpha-free", b_af))
    print(fmt_row("baseline alpha=0", b_a0))

    print()
    print("="*90)
    print("ALPHA-FREE variant (6 eff params / 1 dof)")
    print("="*90)
    grid_af = {}
    for zr in ZR_GRID:
        for ar in AR_GRID:
            p = profile(v, zr, ar, alpha_fixed=None, use_real=True)
            grid_af[(zr, ar)] = p
            print(fmt_row(f"zr={zr} ar={ar:+.1f}", p))

    print()
    print("="*90)
    print("ALPHA=0 PINNED variant (5 eff params / 2 dof)")
    print("="*90)
    grid_a0 = {}
    for zr in ZR_GRID:
        for ar in AR_GRID:
            p = profile(v, zr, ar, alpha_fixed=0.0, use_real=True)
            grid_a0[(zr, ar)] = p
            print(fmt_row(f"zr={zr} ar={ar:+.1f}", p))

    # diagnostics: argmin vs zr at ar=0
    print()
    print("="*90)
    print("DIAGNOSTIC: argmin-theta(zr) at ar=0 (alpha-free / alpha=0)")
    print("="*90)
    for zr in ZR_GRID:
        af = grid_af[(zr, 0.0)]; a0 = grid_a0[(zr, 0.0)]
        print(f"  zr={zr}: af argmin={af['argmin']} Cr={af['Cr']:+.4f} 45in={af['incl45']} | "
              f"a0 argmin={a0['argmin']} Cr={a0['Cr']:+.4f} 45in={a0['incl45']}")

    import pickle
    with open("/Users/solomonwilliams/Reduced Neuronal Model Project/coupled-atlas/_tmp_subB_grid.pkl", "wb") as f:
        pickle.dump(dict(grid_af=grid_af, grid_a0=grid_a0, baseline_af=b_af, baseline_a0=b_a0,
                         ZR_GRID=ZR_GRID, AR_GRID=AR_GRID, THETAS=THETAS.tolist()), f)
    print("\nsaved _tmp_subB_grid.pkl")
