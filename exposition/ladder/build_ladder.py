"""Build the ladder explorable (module delta2) -- the closing argument.

Two claims, pulling in opposite directions, and the module exists to hold them together:

  THE SKELETON CARRIES OVER EXACTLY.  u'' = Y^q u is solved by Bessel functions of order
  1/(q+2) -- 1/3 at the fold, 1/4 at the cusp, 1/5 at the swallowtail.  Verified below at six
  rungs against direct integration.  The family is real, not a coincidence of the first two.

  THE SOLVABILITY DOES NOT.  Fold: Airy, exact.  Cusp: Weber, closed-form and Gamma-type.
  Swallowtail: the cubic oscillator -- exactly characterised by a T-Q system but NOT Gamma-type,
  so there is no closed phase to compare against.  Above that, nothing known.

Hence the closing line of the analysis thread: resurgence is not a stylistic preference, it is
what is left.

NUMERICAL NOTE.  The decaying solution is seeded from WKB at a radius chosen so the action is
about 18.  An earlier draft fixed the radius at 9 for every rung, which at q=4 makes the action
243 and the seed exp(-243) -- underflow, and the verification 'failed' at 63% relative error for
purely arithmetic reasons.  The radius is now set per rung from the action.

Run:  python3 build_ladder.py   ->  ladder.html
"""
import json
import pathlib
import numpy as np
from scipy.special import kv
from scipy.integrate import solve_ivp

HERE = pathlib.Path(__file__).parent
CAPS = json.loads((HERE.parent / 'captions.json').read_text())
REQUIRED = ('headline', 'setup', 'skeleton', 'degrade', 'refuge',
            'distortion_dlmf', 'distortion_absence', 'distortion_notthom')
missing = [k for k in REQUIRED if k not in CAPS.get('delta2', {})]
if missing:
    raise SystemExit(f'captions.json is missing delta2 keys: {missing}')

QS = [1, 2, 3, 4, 5, 6]
NAME = {1: 'fold', 2: 'cusp', 3: 'swallowtail', 4: 'butterfly', 5: 'rung 5', 6: 'rung 6'}
# solvability status, and what it costs.  This is the project's own knowledge, not a computation.
SOLV = {
    1: ('Airy', 'exact', 'The connection problem is Airy. Everything is closed form, and '
                         'Tracy-Widom at this rung is a theorem.'),
    2: ('Weber', 'closed', 'Parabolic cylinder. The connection datum is Gamma-type, so there is '
                           'a closed phase to compare a numerical answer against.'),
    3: ('cubic oscillator', 'characterised', 'Exactly characterised by a T-Q / TBA system — '
                            'equivalently cluster coordinates on a Z/5Z pentagon — but NOT of '
                            'Gamma-type. Exact, and no closed phase.'),
    4: ('—', 'open', 'No reduction known.'),
    5: ('—', 'open', 'No reduction known.'),
    6: ('—', 'open', 'No reduction known.'),
}

CURVE_N = 220
rows, curves = [], {}
for q in QS:
    nu = 1.0 / (q + 2)
    p = (q + 2) / 2.0
    S = lambda Y: (2.0 / (q + 2)) * Y ** p
    # start where the WKB action is ~18: big enough to be deep in the decaying regime, small
    # enough that exp(-S) does not underflow
    Ymax = (18.0 * (q + 2) / 2.0) ** (1.0 / p)
    Ylo = 0.25 * Ymax
    u0 = Ymax ** (-q / 4.0) * np.exp(-S(Ymax))
    du0 = -(Ymax ** (q / 2.0)) * u0
    grid = np.linspace(Ymax, Ylo, CURVE_N)
    sol = solve_ivp(lambda Y, y: [y[1], Y ** q * y[0]], [Ymax, Ylo], [u0, du0],
                    t_eval=grid, rtol=1e-12, atol=1e-300)
    Yv, uv = sol.t[::-1], sol.y[0, ::-1]
    bess = np.sqrt(Yv) * kv(nu, S(Yv))
    # compare SHAPES: both solve the same equation, so they agree up to one overall constant
    scale = uv[len(uv) // 2] / bess[len(bess) // 2]
    rel = float(np.max(np.abs(uv - scale * bess) / np.max(np.abs(uv))))
    assert rel < 1e-6, f'q={q}: Bessel-1/{q+2} skeleton does not match, rel {rel:.2e}'

    fn, status, why = SOLV[q]
    rows.append({'q': q, 'name': NAME[q], 'nu': f'1/{q+2}', 'rel': rel,
                 'fn': fn, 'status': status, 'why': why,
                 'tailExp': 2 * q + 1})
    nrm = float(np.max(np.abs(uv)))
    curves[q] = {'Y': [float(y) for y in Yv],
                 'u': [float(x / nrm) for x in uv],
                 'b': [float(scale * x / nrm) for x in bess],
                 'Ylo': float(Ylo), 'Yhi': float(Ymax)}

ref = {'rows': rows, 'curves': curves, 'qs': QS, 'captions': CAPS['delta2']}
tpl = (HERE / 'ladder.template.html').read_text()
(HERE / 'ladder.html').write_text(tpl.replace('/*__REFERENCE__*/', json.dumps(ref)))

print('wrote ladder.html')
print('  u\'\' = Y^q u  vs  sqrt(Y) K_{1/(q+2)}( 2 Y^{(q+2)/2} / (q+2) ):')
for r in rows:
    print(f"    q={r['q']}  {r['name']:12s}  order {r['nu']:4s}  max rel deviation "
          f"{r['rel']:.2e}   tail exponent {r['tailExp']:2d}   [{r['status']}]")
print(f"  worst deviation across all {len(QS)} rungs: "
      f"{max(r['rel'] for r in rows):.2e}  -- the skeleton is exact")
