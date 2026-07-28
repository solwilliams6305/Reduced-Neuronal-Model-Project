"""Build the junction explorable (module delta1).

PREREQ_GRAPH.md finding F2: the geometry trunk and the asymptotics trunk share no hard edge
until they meet, and they meet at exactly one node -- the coefficient ladder.  Everything
geometric produces it; everything analytic consumes it.  Until this module exists the course is
two disconnected halves rather than one argument.

The quantity that actually passes between them is the INSTANTON ACTION.  On the geometry side it
is an elementary area.  On the analysis side the same number sets the tail exponent AND the
location of the Borel singularity, hence the divergence rate of the v_k.

Everything here is computed:
  * the action law, by direct quadrature against the closed form, at every rung;
  * the cusp coefficient ladder, taken from the published values used in W_make_figures.py;
  * the Domb-Sykes style ratio |v_k|/|v_{k+1}|, which is what the coefficients say about where
    the nearest Borel singularity is.

Run:  python3 build_junction.py   ->  junction.html
"""
import json
import math
import pathlib
import numpy as np
from scipy.integrate import quad

HERE = pathlib.Path(__file__).parent
CAPS = json.loads((HERE.parent / 'captions.json').read_text())
REQUIRED = ('headline', 'setup', 'geometry_side', 'analysis_side', 'the_family',
            'distortion_factor2', 'distortion_leading', 'distortion_onerung',
            'ratio_setup', 'ratio_finding', 'ratio_limit')
missing = [k for k in REQUIRED if k not in CAPS.get('delta1', {})]
if missing:
    raise SystemExit(f'captions.json is missing delta1 keys: {missing}')

QS = [1, 2, 3, 4, 5]
NAMES = {1: 'fold', 2: 'cusp', 3: 'swallowtail', 4: 'butterfly', 5: 'rung 5'}

# ---- the action law, verified by quadrature rather than quoted -----------------------------
family, worst = [], 0.0
for q in QS:
    coef = 1.0 / (2 * (2 * q + 1))
    for s in (0.7, 2.0, 6.0):
        num = 0.5 * quad(lambda t: (t ** q) ** 2, 0.0, s)[0]
        closed = coef * s ** (2 * q + 1)
        worst = max(worst, abs(num - closed) / max(1.0, abs(closed)))
    family.append({'q': q, 'name': NAMES[q], 'exp': 2 * q + 1,
                   'den': 2 * (2 * q + 1), 'coef': coef})

# BVP-measured ratios recorded in coupled-atlas/instanton_action_q.py (s=10).  Quoted, not
# recomputed here: those are boundary-value solves and belong on the Python side of the split.
bvp = {1: 0.839, 2: 0.930, 3: 0.981, 4: 0.997, 5: 0.9995}
for f in family:
    f['bvpRatio'] = bvp[f['q']]

# ---- the cusp ladder, and what it says about the singularity -------------------------------
# The published cusp coefficients (same values as W_make_figures.py).
v = np.array([0.134, 0.111, 0.104, -0.030, -0.451, -1.19, -1.90])
kk = np.arange(len(v))
fact = np.array([float(math.factorial(int(i))) for i in kk])
b = v / fact                                   # Borel coefficients v_k / k!

# Signed ratios.  A single singularity on the positive real axis at distance A would drive these
# to 1/A.  They do not settle -- and that is the finding, not a defect.
ratios = [{'k': int(k), 'r': float(b[k + 1] / b[k])} for k in range(len(b) - 1)]

# Does a complex-conjugate pair account for the swing?  Model b_k ~ rho^-k cos(k theta - phi)
# with rho and theta FIXED to the reported Borel data; only the phase and an overall scale are
# fitted, so this is a two-parameter check of a stated result, not a free fit.
RHO, THETA = 1.9, np.deg2rad(50.0)
best = None
for phi in np.linspace(0, 2 * np.pi, 2001):
    m = RHO ** (-kk.astype(float)) * np.cos(kk * THETA - phi)
    sc = float(np.dot(m, b) / np.dot(m, m))
    res = float(np.sum((sc * m - b) ** 2) / np.sum(b ** 2))
    if best is None or res < best[0]:
        best = (res, float(phi), sc, m)
mres, mphi, msc, mvec = best
modelRatios = [{'k': int(k), 'r': float((msc * mvec)[k + 1] / (msc * mvec)[k])}
               for k in range(len(b) - 1)]

# Root test: what the modulus would look like if you tried to read it off directly.
rootTest = [{'k': int(i), 'val': float(abs(b[i]) ** (-1.0 / i))} for i in range(1, len(b))]

ref = {
    'family': family,
    'v': [float(x) for x in v],
    'b': [float(x) for x in b],
    'ratios': ratios,
    'modelRatios': modelRatios,
    'rootTest': rootTest,
    'model': {'rho': RHO, 'theta': 50.0, 'phi': float(np.degrees(mphi)), 'residual': mres},
    'invA': 1.0 / RHO,
    'captions': CAPS['delta1'],
}

# fail loudly if the headline claim stops holding
assert mres < 0.02, f'conjugate-pair model no longer fits the ladder: residual {mres:.4f}'

tpl = (HERE / 'junction.template.html').read_text()
(HERE / 'junction.html').write_text(tpl.replace('/*__REFERENCE__*/', json.dumps(ref)))

print('wrote junction.html')
print(f'  action law verified by quadrature at {len(QS)} rungs x 3 depths;'
      f' worst relative error {worst:.2e}')
for f in family:
    print(f'    q={f["q"]:d} {f["name"]:11s} I(s) = s^{f["exp"]:d}/{f["den"]:d}'
          f'   BVP/closed at s=10: {f["bvpRatio"]}')
print(f'  cusp ladder: {len(v)} coefficients, {len(ratios)} consecutive Borel ratios')
print(f'  conjugate-pair model rho={RHO}, theta=50deg, phi={np.degrees(mphi):.1f}deg'
      f'  -> relative residual {mres:.4f}')
print(f'  data ratios  : {[round(r["r"], 3) for r in ratios]}')
print(f'  model ratios : {[round(r["r"], 3) for r in modelRatios]}')
print(f'  root test |b_k|^(-1/k) at k=4,5,6 : '
      f'{[round(r["val"], 2) for r in rootTest[-3:]]}  (true |zeta| ~ {RHO}) -> NOT converged')
print('  NB the swallowtail ladder is in flux (v7 revised 2026-07-28) and is NOT used here')
