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
import pathlib
import numpy as np
from scipy.integrate import quad

HERE = pathlib.Path(__file__).parent
CAPS = json.loads((HERE.parent / 'captions.json').read_text())
REQUIRED = ('headline', 'setup', 'geometry_side', 'analysis_side', 'the_family',
            'distortion_factor2', 'distortion_leading', 'distortion_onerung')
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
v = [0.134, 0.111, 0.104, -0.030, -0.451, -1.19, -1.90]
ratios = []
for k in range(len(v) - 1):
    if v[k] != 0:
        ratios.append({'k': k, 'r': float(abs(v[k]) / abs(v[k + 1]))})

ref = {
    'family': family,
    'v': v,
    'ratios': ratios,
    'captions': CAPS['delta1'],
}

tpl = (HERE / 'junction.template.html').read_text()
(HERE / 'junction.html').write_text(tpl.replace('/*__REFERENCE__*/', json.dumps(ref)))

print('wrote junction.html')
print(f'  action law verified by quadrature at {len(QS)} rungs x 3 depths;'
      f' worst relative error {worst:.2e}')
for f in family:
    print(f'    q={f["q"]:d} {f["name"]:11s} I(s) = s^{f["exp"]:d}/{f["den"]:d}'
          f'   BVP/closed at s=10: {f["bvpRatio"]}')
print(f'  cusp ladder: {len(v)} coefficients, {len(ratios)} consecutive ratios')
print('  NB the swallowtail ladder is in flux (v7 revised 2026-07-28) and is NOT used here')
