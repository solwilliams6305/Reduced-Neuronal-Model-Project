"""Robustness of theta=48+-3deg to ladder uncertainties.
Theta-fixed single complex-pair Darboux fit; scan theta 38..66 step 1deg.
For each combo of (v2,v3,v4,v5,v6) record argmin-theta and 1.25x-min interval."""
import sys, itertools
import numpy as np
from math import pi, factorial
from scipy.special import gamma as Gamma
from scipy.optimize import least_squares
from numpy.polynomial import polynomial as P

V0, V1 = 0.134, 0.111

def darboux_model(C, zeta, theta, phi, alpha, n):
    return 2*C*zeta**(-(n+1.0))*Gamma(n+1.0+alpha)*np.cos((n+1.0)*theta - phi)

def fit_fixed_theta(v, theta_rad):
    """Fit C,zeta,phi,alpha at fixed theta. Multi-start phi,zeta. Return max-rel-resid."""
    ns = np.arange(len(v), dtype=float)
    y = np.asarray(v, float)
    scale = np.abs(y) + 1e-3
    def resid(p):
        C, zeta, phi, alpha = p
        return (darboux_model(C, zeta, theta_rad, phi, alpha, ns) - y)/scale
    best = None
    for phi0 in np.radians([0, 60, 120, 180, 240, 300]):
        for z0 in [1.0, 1.4, 1.8, 2.2, 2.8]:
            p0 = [0.1, z0, phi0, 0.0]
            lo = [1e-6, 0.3, -2*pi, -1.5]
            hi = [10.0, 5.0, 2*pi, 1.5]
            try:
                r = least_squares(resid, p0, bounds=(lo, hi), xtol=1e-14, ftol=1e-14)
            except ValueError:
                continue
            if best is None or r.cost < best.cost:
                best = r
    return np.max(np.abs(best.fun))

def theta_scan(v, thetas_deg):
    res = {}
    for td in thetas_deg:
        res[td] = fit_fixed_theta(v, np.radians(td))
    return res

def pade(c, L, M):
    c = np.asarray(c, float); N = L+M
    A = np.zeros((M, M)); b = np.zeros(M)
    for i in range(1, M+1):
        for j in range(1, M+1):
            k = L+i-j
            A[i-1, j-1] = c[k] if k >= 0 else 0.0
        b[i-1] = -c[L+i]
    qtail = np.linalg.solve(A, b)
    q = np.concatenate([[1.0], qtail])
    return q

def borel_pade_poles(v, L, M):
    b = [v[n]/factorial(n) for n in range(len(v))]
    q = pade(b, L, M)
    roots = P.polyroots(q)
    return sorted(roots, key=abs)

# ---- grids ----
V6 = [-1.75, -1.90, -2.05]
V5 = [-1.18, -1.19, -1.20]
V4 = [-0.448, -0.451, -0.454]
V3 = [-0.020, -0.030]
V2 = [0.100, 0.104]

thetas = list(range(38, 67))
combos = list(itertools.product(V2, V3, V4, V5, V6))
print(f"# combos = {len(combos)}; thetas {thetas[0]}..{thetas[-1]}")

records = []
for (v2, v3, v4, v5, v6) in combos:
    v = [V0, V1, v2, v3, v4, v5, v6]
    sc = theta_scan(v, thetas)
    rmin = min(sc.values())
    argmin = min(sc, key=sc.get)
    band = sorted([t for t in thetas if sc[t] <= 1.25*rmin])
    records.append(dict(v2=v2, v3=v3, v4=v4, v5=v5, v6=v6,
                        argmin=argmin, rmin=rmin,
                        band_lo=band[0], band_hi=band[-1],
                        r45=sc[45], r53=sc[53], r56=sc[56],
                        contig=(band == list(range(band[0], band[-1]+1))),
                        sc=sc))

# ---- (a) argmin range ----
argmins = [r['argmin'] for r in records]
print(f"\n(a) argmin range: {min(argmins)}..{max(argmins)}  (median {int(np.median(argmins))})")
from collections import Counter
print("    argmin distribution:", dict(sorted(Counter(argmins).items())))

# ---- (b) union band worst case ----
lo = min(r['band_lo'] for r in records)
hi = max(r['band_hi'] for r in records)
print(f"\n(b) UNION 1.25x band (worst case): {lo}..{hi} deg")
print(f"    per-combo band widths: min {min(r['band_hi']-r['band_lo'] for r in records)}, "
      f"max {max(r['band_hi']-r['band_lo'] for r in records)} deg")

# ---- (c) sensitivity ranking (one-at-a-time from central combo) ----
central = dict(v2=0.104, v3=-0.030, v4=-0.451, v5=-1.19, v6=-1.90)
def argmin_for(vv):
    v = [V0, V1, vv['v2'], vv['v3'], vv['v4'], vv['v5'], vv['v6']]
    sc = theta_scan(v, thetas)
    rmin = min(sc.values()); am = min(sc, key=sc.get)
    band = sorted([t for t in thetas if sc[t] <= 1.25*rmin])
    return am, band[0], band[-1]
print("\n(c) One-at-a-time sensitivity (central v6=-1.90,v5=-1.19,v4=-0.451,v3=-0.030,v2=0.104):")
base_am, base_lo, base_hi = argmin_for(central)
print(f"    CENTRAL: argmin={base_am} band {base_lo}..{base_hi}")
sens = {}
for key, vals in [('v6', V6), ('v5', V5), ('v4', V4), ('v3', V3), ('v2', V2)]:
    ams = []
    for val in vals:
        c2 = dict(central); c2[key] = val
        am, blo, bhi = argmin_for(c2)
        ams.append((val, am, blo, bhi))
    spread = max(a[1] for a in ams) - min(a[1] for a in ams)
    band_spread = max(a[3] for a in ams) - min(a[2] for a in ams)
    sens[key] = (spread, band_spread, ams)
    print(f"    {key}: argmin spread {spread}deg, band-union spread {band_spread}deg  "
          + " | ".join(f"{v:+.4g}->am{am}[{blo}-{bhi}]" for v, am, blo, bhi in ams))
print("\n    Sensitivity ranking (by argmin spread):",
      sorted(sens, key=lambda k: -sens[k][0]))

# ---- (d) >=53 exclusion ----
never53 = all(r['band_hi'] < 53 for r in records)
max_bandhi = max(r['band_hi'] for r in records)
n_reach53 = sum(1 for r in records if r['band_hi'] >= 53)
print(f"\n(d) >=53deg exclusion survives every combo? {never53}  "
      f"(max band_hi over all combos = {max_bandhi}; #combos whose band reaches>=53: {n_reach53}/{len(records)})")

# ---- (e) 45 inside band ----
n45 = sum(1 for r in records if r['band_lo'] <= 45 <= r['band_hi'])
print(f"\n(e) 45deg inside 1.25x band: {n45}/{len(records)} combos "
      f"({'ALL' if n45==len(records) else f'{100*n45/len(records):.0f}%'})")
r45_ratio = [r['r45']/r['rmin'] for r in records]
print(f"    r(45)/rmin: min {min(r45_ratio):.3f}, max {max(r45_ratio):.3f}, "
      f"median {np.median(r45_ratio):.3f}")

# ---- Borel-Pade poles for extreme combos ----
print("\n--- Borel-Pade pair-angle for extreme combos ---")
extremes = {
    'central': [V0,V1,0.104,-0.030,-0.451,-1.19,-1.90],
    'v6=-1.75': [V0,V1,0.104,-0.030,-0.451,-1.19,-1.75],
    'v6=-2.05': [V0,V1,0.104,-0.030,-0.451,-1.19,-2.05],
    'v3=-0.020': [V0,V1,0.104,-0.020,-0.451,-1.19,-1.90],
    'lo-corner': [V0,V1,0.100,-0.020,-0.448,-1.18,-1.75],
    'hi-corner': [V0,V1,0.104,-0.030,-0.454,-1.20,-2.05],
}
for name, v in extremes.items():
    K = len(v)-1
    angs = []
    for M in range(2, K):
        L = K-M
        if L < 0: continue
        try:
            roots = borel_pade_poles(v, L, M)
        except np.linalg.LinAlgError:
            continue
        for r in roots:
            ang = np.degrees(np.angle(r))
            if abs(ang) > 5 and abs(ang) < 175:  # complex pair only
                angs.append(abs(ang))
    if angs:
        print(f"  {name:12s}: complex-pair angles {sorted(round(a,1) for a in angs)}  "
              f"(min {min(angs):.1f}, max {max(angs):.1f})")
    else:
        print(f"  {name:12s}: no clear complex pair")

# save summary
import pickle
with open('_tmp_theta_robust_records.pkl','wb') as f:
    pickle.dump(records, f)
print("\nsaved _tmp_theta_robust_records.pkl")
