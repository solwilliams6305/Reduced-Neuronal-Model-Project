"""Sub-agent A: real-pole-subtracted Borel phase extraction.
Pinned real instanton (zr,ar on grid); theta FIXED on scan; fit remaining params by LS on
relative residuals (as in _borel_analysis). Profile max-rel-resid over theta.

Conventions EXACTLY as _borel_analysis.darboux_two_sing (order index k=n+1):
  v_n = 2C zeta^-(n+1) Gamma(n+1+alpha) cos((n+1)theta - phi)  +  Cr zr^-(n+1) Gamma(n+1+ar)
"""
import numpy as np
from math import pi
from scipy.special import gamma as Gamma
from scipy.optimize import least_squares

V = np.array([0.134, 0.111, 0.104, -0.030, -0.451, -1.19, -1.90])

def pair_term(C, zeta, theta, phi, alpha, n):
    return 2*C*zeta**(-(n+1.0))*Gamma(n+1.0+alpha)*np.cos((n+1.0)*theta-phi)

def real_term(Cr, zr, ar, n):
    return Cr*zr**(-(n+1.0))*Gamma(n+1.0+ar)

# bounds for the fitted params (theta held fixed outside)
# order: C, zeta, phi, alpha, Cr
LO = np.array([1e-6, 0.3, -2*pi, -1.4, -10.0])
HI = np.array([10.0, 5.0,  2*pi,  1.4,  10.0])

def _railed(p, lo, hi, tol=0.01):
    """flag any param within 1% of a bound (relative to bound span)."""
    span = hi - lo
    return bool(np.any((p - lo) < tol*span) or np.any((hi - p) < tol*span))

def fit_fixed_theta(v, theta, zr=None, ar=None, alpha_fixed=None):
    """Fit at fixed theta. If zr is None -> baseline (no real pole).
    Returns dict with params, maxrelresid, railed flag."""
    ns = np.arange(len(v), dtype=float)
    y = np.asarray(v, float)
    scale = np.abs(y) + 1e-3
    have_real = zr is not None
    free_alpha = alpha_fixed is None

    # active param indices among [C, zeta, phi, alpha, Cr]
    idx = [0, 1, 2]                       # C, zeta, phi always
    if free_alpha: idx.append(3)         # alpha
    if have_real:  idx.append(4)         # Cr
    idx = np.array(idx)
    lo = LO[idx]; hi = HI[idx]

    def expand(p):
        full = np.array([0.1, 1.5, 0.5, (0.0 if free_alpha else alpha_fixed), 0.0])
        full[idx] = p
        return full

    def resid(p):
        C, zeta, phi, alpha, Cr = expand(p)
        m = pair_term(C, zeta, theta, phi, alpha, ns)
        if have_real:
            m = m + real_term(Cr, zr, ar, ns)
        return (m - y) / scale

    best = None
    # multistart: phi x zeta x Cr (+ alpha starts)
    for phi0 in [0.3, 1.0, 1.8, 2.6]:
        for z0 in [1.5, 2.0, 2.5]:
            crs = [0.05, -0.05] if have_real else [0.0]
            for cr0 in crs:
                a_starts = [-0.5, 0.0, 0.5] if free_alpha else [alpha_fixed]
                for a0 in a_starts:
                    full0 = np.array([0.1, z0, phi0, a0, cr0])
                    p0 = np.clip(full0[idx], lo, hi)
                    try:
                        r = least_squares(resid, p0, bounds=(lo, hi),
                                          xtol=1e-15, ftol=1e-15, max_nfev=8000)
                    except ValueError:
                        continue
                    if best is None or r.cost < best.cost:
                        best = r
    C, zeta, phi, alpha, Cr = expand(best.x)
    return dict(C=C, zeta=zeta, phi=phi, alpha=alpha, Cr=Cr,
                res=float(np.max(np.abs(best.fun))),
                railed=_railed(best.x, lo, hi))

THETAS = np.arange(38, 61, 1)  # degrees

def profile(v, zr=None, ar=None, alpha_fixed=None):
    """theta scan -> list of per-theta fit dicts."""
    out = []
    for thd in THETAS:
        d = fit_fixed_theta(v, np.radians(thd), zr=zr, ar=ar, alpha_fixed=alpha_fixed)
        d['theta'] = thd
        out.append(d)
    return out

def summarize(prof):
    res = np.array([d['res'] for d in prof])
    imin = int(np.argmin(res))
    rmin = res[imin]
    band_mask = res <= 1.25*rmin
    band_th = THETAS[band_mask]
    dmin = prof[imin]
    inband45 = bool(45 in band_th)
    return dict(argmin=int(THETAS[imin]), rmin=float(rmin),
                band=(int(band_th.min()), int(band_th.max())),
                band_list=[int(x) for x in band_th],
                Cr=dmin['Cr'], alpha=dmin['alpha'], zeta=dmin['zeta'],
                railed=dmin['railed'], inband45=inband45, res=res)

def cell_line(tag, s):
    band = f"{s['band'][0]}-{s['band'][1]}"
    r = "RAILED" if s['railed'] else "ok"
    return (f"{tag:24s} argmin={s['argmin']:>3}  band[{band:>7}]  "
            f"Cr={s['Cr']:+.4f}  alpha={s['alpha']:+.3f}  zeta={s['zeta']:.3f}  "
            f"45inband={str(s['inband45']):5s}  rmin={s['rmin']:.4f}  {r}")

if __name__ == "__main__":
    import sys
    ZR_GRID = [1.4, 1.6, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6]
    AR_GRID = [-0.5, 0.0, 0.5]

    print("="*100)
    print("BASELINE (no real pole)")
    print("="*100)
    for af, lab in [(None, "alpha-free (5p:C,z,th,phi,a)"), (0.0, "alpha=0  (4p:C,z,th,phi)")]:
        s = summarize(profile(V, alpha_fixed=af))
        print(cell_line("baseline "+lab, s))

    print("\n"+"="*100)
    print("ALPHA-FREE variant (6 free params total incl theta; 1 dof)")
    print("="*100)
    resultsA = {}
    for zr in ZR_GRID:
        for ar in AR_GRID:
            s = summarize(profile(V, zr=zr, ar=ar, alpha_fixed=None))
            resultsA[(zr, ar)] = s
            print(cell_line(f"zr={zr:.1f} ar={ar:+.1f}", s))

    print("\n"+"="*100)
    print("ALPHA=0 PINNED variant (5 free params total incl theta; 2 dof)")
    print("="*100)
    resultsB = {}
    for zr in ZR_GRID:
        for ar in AR_GRID:
            s = summarize(profile(V, zr=zr, ar=ar, alpha_fixed=0.0))
            resultsB[(zr, ar)] = s
            print(cell_line(f"zr={zr:.1f} ar={ar:+.1f}", s))

    # central cell: zr=2.0, ar=0.0 (mid of Pade 1.9-2.4 window, neutral exponent)
    print("\n"+"="*100)
    print("ROBUSTNESS on central cell zr=2.0 ar=0.0 (v6 and v3 variants)")
    print("="*100)
    for v6 in [-1.75, -1.90, -2.05]:
        for v3 in [-0.020, -0.030]:
            vv = V.copy(); vv[3] = v3; vv[6] = v6
            sA = summarize(profile(vv, zr=2.0, ar=0.0, alpha_fixed=None))
            sB = summarize(profile(vv, zr=2.0, ar=0.0, alpha_fixed=0.0))
            print(f"v3={v3:+.3f} v6={v6:+.2f} | A(af): "+cell_line("", sA))
            print(f"                    | B(a0): "+cell_line("", sB))

    # diagnostic: argmin_theta(zr) trend at ar=0 (both variants)
    print("\n"+"="*100)
    print("DIAGNOSTIC: argmin-theta vs zr at ar=0")
    print("="*100)
    print("  zr :   A(alpha-free) argmin / Cr     |  B(alpha=0) argmin / Cr")
    for zr in ZR_GRID:
        a = resultsA[(zr,0.0)]; b = resultsB[(zr,0.0)]
        print(f"  {zr:.1f}:   {a['argmin']:>3} / {a['Cr']:+.4f}  (45inb={a['inband45']})   |  "
              f"{b['argmin']:>3} / {b['Cr']:+.4f}  (45inb={b['inband45']})")
