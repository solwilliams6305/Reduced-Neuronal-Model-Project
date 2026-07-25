"""
JOINT (variance + skew) large-order fit: two observables, ONE Borel geometry.

The point.  Resurgence says the Borel singularity LOCATIONS are properties of the action, not of the
observable.  So Var and kappa_3 -- expanded in the same variable x=eta^2 -- share (|zeta|, theta)
while carrying their own amplitudes, phases and exponents:

    v_k ~ 2 C_v |zeta|^-(k+1) Gamma(k+1+a_v) cos((k+1) theta - phi_v)
    t_j ~ 2 C_t |zeta|^-(j+1) Gamma(j+1+a_t) cos((j+1) theta - phi_t)

Fitting them TOGETHER puts 13 data points (7 v's + 6 t's) against 8 parameters, of which only 2 are
shared, instead of 7 points against 5.  Unlike the conformal map and the direct non-perturbative
extraction -- both of which needed (A, theta) as INPUT and so were circular -- this route only needs
the coefficients we can actually compute.

Structural corroboration already visible without any fitting: within each catastrophe BOTH ladders
first turn negative at the SAME rung (cusp: v_3 and t_3; swallowtail: v_5 and t_5), which is what a
shared cos(k theta - phi) envelope demands, and the later swallowtail node implies a SMALLER theta.

VALIDATION GATE: the cusp, where theta = 50 +- 2 deg and |zeta| ~ 1.9 are published.

Run:  python3 swtl_joint.py
"""
import json, os, glob
import numpy as np
from math import pi
from scipy.special import gamma as Gamma
from scipy.optimize import least_squares

import swtl_borel as S

HERE = os.path.dirname(os.path.abspath(__file__))

CUSP_V = S.CUSP_V
CUSP_T = [0.0587778562375058, 0.15256611939205722, 0.19098597184122607,
          -0.1867840040972798, -2.0721239879114206, -4.731359476811322]      # _w_kappa3_n20.json


def _ladder_model(C, zeta, theta, phi, alpha, n):
    return 2 * C * zeta ** (-(n + 1.0)) * Gamma(n + 1.0 + alpha) * np.cos((n + 1.0) * theta - phi)


def joint_fit(v, t, alpha_v=0.0, alpha_t=0.0, verbose=True, label=""):
    """Shared (|zeta|, theta); per-ladder (C, phi).  alphas fixed for over-fit protection."""
    nv = np.arange(len(v), dtype=float)
    nt = np.arange(len(t), dtype=float)
    yv = np.asarray(v, float); yt = np.asarray(t, float)
    sv = np.abs(yv) + 1e-3; st = np.abs(yt) + 1e-3

    def resid(p):
        zeta, theta, Cv, phiv, Ct, phit = p
        rv = (_ladder_model(Cv, zeta, theta, phiv, alpha_v, nv) - yv) / sv
        rt = (_ladder_model(Ct, zeta, theta, phit, alpha_t, nt) - yt) / st
        return np.concatenate([rv, rt])

    best = None
    for th0 in np.radians([20, 28, 35, 42, 50, 58, 66]):
        for z0 in (0.8, 1.1, 1.5, 1.9, 2.5):
            p0 = [z0, th0, 0.05, 0.5, 0.05, 0.5]
            lo = [0.3, np.radians(5), -50, -2 * pi, -50, -2 * pi]
            hi = [6.0, np.radians(120), 50, 2 * pi, 50, 2 * pi]
            try:
                r = least_squares(resid, p0, bounds=(lo, hi), xtol=1e-14, ftol=1e-14)
            except ValueError:
                continue
            if best is None or r.cost < best.cost:
                best = r
    zeta, theta, Cv, phiv, Ct, phit = best.x
    res = float(np.max(np.abs(best.fun)))
    if verbose:
        print(f"  {label:>26}: |zeta|={zeta:6.3f}  theta={np.degrees(theta):6.1f} deg   "
              f"maxrelresid={res:6.3f}")
    return zeta, np.degrees(theta), res


def variance_only(v, alpha=0.0, label=""):
    n = np.arange(len(v), dtype=float)
    y = np.asarray(v, float); sc = np.abs(y) + 1e-3

    def resid(p):
        C, zeta, theta, phi = p
        return (_ladder_model(C, zeta, theta, phi, alpha, n) - y) / sc

    best = None
    for th0 in np.radians([20, 28, 35, 42, 50, 58, 66]):
        for z0 in (0.8, 1.1, 1.5, 1.9, 2.5):
            try:
                r = least_squares(resid, [0.05, z0, th0, 0.5],
                                  bounds=([-50, 0.3, np.radians(5), -2 * pi],
                                          [50, 6.0, np.radians(120), 2 * pi]),
                                  xtol=1e-14, ftol=1e-14)
            except ValueError:
                continue
            if best is None or r.cost < best.cost:
                best = r
    C, zeta, theta, phi = best.x
    print(f"  {label:>26}: |zeta|={zeta:6.3f}  theta={np.degrees(theta):6.1f} deg   "
          f"maxrelresid={float(np.max(np.abs(best.fun))):6.3f}")
    return zeta, np.degrees(theta)


def load_skew(verbose=False):
    """Extrapolated skew ladder t_j(n->inf), plus the raw per-n ladders.

    Each t_j is extrapolated across the grid with the SAME screened estimator used for the
    variance ladder (wrong-side / bad-fit / sign-flip / plateau detection), because the skew
    coefficients show the same behaviour: t_5 runs -7.99, -6.40, -5.65, -5.39 over n=12..24,
    decelerating faster than 1/n -- the plateau signature, where raw beats Richardson.
    Grids with n < 2+2j are dropped (the resolution rule: n points cannot resolve Y_idx, idx>n).
    """
    files = sorted(glob.glob(os.path.join(HERE, 'swtl_kappa3_n*.json')),
                   key=lambda f: int(f.split('_n')[-1].split('.')[0]))
    if not files:
        return None, {}
    per_n = {}
    for f in files:
        d = json.load(open(f))
        per_n[d['n']] = d['t']
    jmax = max(len(v) for v in per_n.values()) - 1
    # Use only grids that computed the FULL ladder.  Partial (smoke-test) runs are not merely
    # shorter, they are coarser: the n=8 jmax=2 probe gives t_0=0.0265 against 0.0175-0.0178 for
    # n>=12, a 49% outlier, even though n=8 clears the index rule for t_0.  The index rule is
    # necessary, not sufficient -- so require membership of the production grid.
    dropped = {n for n, tv in per_n.items() if len(tv) != jmax + 1}
    if dropped and verbose:
        print(f"    [dropping partial/smoke grids n={sorted(dropped)}: not full-ladder runs]")
    per_n = {n: tv for n, tv in per_n.items() if len(tv) == jmax + 1}
    out, bands = [], []
    for j in range(jmax + 1):
        need = 2 + 2 * j                      # max chaos index entering t_j
        pts = sorted((n, tv[j]) for n, tv in per_n.items() if len(tv) > j and n >= need)
        if not pts:
            break
        if len(pts) == 1:
            out.append(pts[0][1]); bands.append((pts[0][1], pts[0][1]))
            continue
        best, lo, hi, ests, note = S.robust_extrap([p[0] for p in pts], [p[1] for p in pts])
        out.append(best); bands.append((lo, hi))
        if verbose:
            print(f"    t_{j}: {best:+.5f}  band [{lo:+.5f},{hi:+.5f}]  "
                  f"grid {[(n, round(x,5)) for n, x in pts]}" + (f"   <<< {note}" if note else ""))
    return out, per_n


ALPHAS = (-0.5, 0.0, 0.5)

if __name__ == "__main__":
    print("=" * 78)
    print("GATE: cusp (published theta = 50 +- 2 deg, |zeta| ~ 1.9)")
    print("=" * 78)
    print("  variance ladder ALONE (7 points, 4 params):")
    for a in ALPHAS:
        variance_only(CUSP_V, alpha=a, label=f"alpha={a:+.1f}")
    print("  JOINT variance + skew (13 points, 6 params, 2 shared):")
    tj = []
    for a in ALPHAS:
        z, th, r = joint_fit(CUSP_V, CUSP_T, alpha_v=a, alpha_t=a, label=f"alpha={a:+.1f}")
        tj.append(th)
    print(f"\n  cusp joint theta spread over alpha: {min(tj):.1f}-{max(tj):.1f} deg "
          f"(median {np.median(tj):.1f})   published 50+-2")

    print("\n" + "=" * 78)
    print("SWALLOWTAIL q=3")
    print("=" * 78)
    v, band, cont = S.load_ladder(verbose=False)
    print("  skew ladder, extrapolated per coefficient:")
    t, per_n = load_skew(verbose=True)
    if t is None:
        print("  no skew ladder yet (swtl_kappa3_n*.json missing) -- run swtl_kappa3.py first")
        raise SystemExit(0)
    print(f"  skew grids available: n={sorted(per_n)}")
    print(f"  v = {[round(x,4) for x in cont]}")
    print(f"  t = {[round(x,4) for x in t]}")
    print("\n  variance ladder ALONE:")
    for a in ALPHAS:
        variance_only(cont, alpha=a, label=f"alpha={a:+.1f}")
    print("  JOINT variance + skew:")
    ths = []
    for a in ALPHAS:
        z, th, r = joint_fit(cont, t, alpha_v=a, alpha_t=a, label=f"alpha={a:+.1f}")
        ths.append(th)
    print(f"\n  swallowtail joint theta spread over alpha: {min(ths):.1f}-{max(ths):.1f} deg "
          f"(median {np.median(ths):.1f})")
    print("\n  stability across the skew grid (theta at alpha=0, per skew n):")
    for nn in sorted(per_n):
        z, th, r = joint_fit(cont, per_n[nn], alpha_v=0.0, alpha_t=0.0, verbose=False)
        print(f"    skew n={nn:>3}: theta={th:6.1f} deg  |zeta|={z:.3f}")
