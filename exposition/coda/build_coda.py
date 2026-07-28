"""Build the coda explorable (module delta3): three families, one word.

Why this module exists.  ARCHITECTURE.md sequences delta3 as "the two-ladders split (swept vs
genuine Berry-Upstill)".  That is now one family short.  The swallowtail paper's positioning
section distinguishes THREE constructions that share the word "multicritical":

  1. swept multicritical      u'' = (sign(Y)|Y|^q - eta xi) u        <- this course.  NOT integrable.
  2. genuine Berry-Upstill    u^(n) = (Y - eta xi) u                 <- n-th order equation.
  3. RMT multicritical edge   higher-order Tracy-Widom               <- Painleve I/II hierarchies.

Family 3 is the dangerous one and is the reason this module was upgraded: it collides on BOTH
words ("multicritical", and "higher-order analogue of Tracy-Widom"), and it is integrable by
construction -- Fredholm determinants of higher Airy-type kernels -- whereas family 1 is provably
not.  A reader who arrives expecting "multicritical edge => Painleve hierarchy" will expect exactly
the wrong answer about everything this course has built.

Everything on the swept column is COMPUTED here, not quoted:
  * the recessive solution phi_q of phi'' = sign(Y)|Y|^q phi, integrated from a WKB start;
  * its first node swept from +infinity -- the escape level -- by bracketed root-find;
  * the leading edge-variance v_0^(q) = int phi^4 / phi'(node)^4 by quadrature.
The q=1 node must come out as the first Airy zero (-2.338107) and v_0^(3) as 0.04953; both are
asserted below and the build fails if they drift.

The skew values are MEASURED elsewhere (Monte-Carlo for swept, numerical integration for
Berry-Upstill) and are carried as data with that provenance stated on-screen -- see the
distortion caption `distortion_skew`.

Run:  python3 build_coda.py   ->  coda.html
"""
import json
import pathlib

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

HERE = pathlib.Path(__file__).parent
CAPS = json.loads((HERE.parent / 'captions.json').read_text())
REQUIRED = ('headline', 'setup', 'swept', 'genuine', 'rmt', 'divergence', 'falsify',
            'distortion_slice', 'distortion_thirdcol', 'distortion_skew')
missing = [k for k in REQUIRED if k not in CAPS.get('delta3', {})]
if missing:
    raise SystemExit(f'captions.json is missing delta3 keys: {missing}')

QS = [1, 2, 3, 4, 5]
NAMES = {1: 'fold  A2', 2: 'cusp  A3', 3: 'swallowtail  A4', 4: 'butterfly  A5', 5: 'A6'}


def swept_rung(q, Y0=6.0):
    """Recessive solution, its first node, and v_0^(q) -- all by integration/quadrature."""
    def rhs(Y, y):
        return [y[1], np.sign(Y) * abs(Y) ** q * y[0]]

    S = 2.0 / (q + 2) * Y0 ** ((q + 2) / 2)          # WKB action at the start
    y0 = [np.exp(-S), -Y0 ** (q / 2) * np.exp(-S)]   # recessive branch
    sol = solve_ivp(rhs, [Y0, -6.0], y0, rtol=1e-11, atol=1e-30,
                    dense_output=True, max_step=0.01)
    f = lambda Y: sol.sol(Y)[0]
    grid = np.linspace(Y0, -6.0, 20000)
    vals = f(grid)
    cross = np.where(np.sign(vals[:-1]) != np.sign(vals[1:]))[0]
    node = brentq(f, grid[cross[0]], grid[cross[0] + 1])
    dphi = sol.sol(node)[1]
    v0 = quad(lambda Y: f(Y) ** 4, node, Y0, limit=400)[0] / dphi ** 4
    return float(node), float(v0)


swept = []
for q in QS:
    node, v0 = swept_rung(q)
    swept.append({
        'q': q,
        'name': NAMES[q],
        'node': node,
        'v0': v0,
        'besselOrder': f'1/{q + 2}',
        'instExp': 2 * q + 1,
        'instDen': 2 * (2 * q + 1),
    })

# --- assertions: the build fails rather than shipping drifted numbers ---
airy1 = next(s for s in swept if s['q'] == 1)['node']
assert abs(airy1 - (-2.338107)) < 2e-6, f'q=1 node {airy1} is not the first Airy zero'
v0_swtl = next(s for s in swept if s['q'] == 3)['v0']
assert abs(v0_swtl - 0.04953) < 5e-5, f'v_0^(3) = {v0_swtl}, expected 0.04953'
v0_cusp = next(s for s in swept if s['q'] == 2)['v0']
assert abs(v0_cusp - 0.13429) < 5e-5, f'v_0^(2) = {v0_cusp}, expected 0.13429'

# --- measured skew trends: DATA, not computed here (provenance stated on-screen) ---
skew = {
    'swept':   [{'i': 1, 'v': 0.20}, {'i': 2, 'v': 0.61}, {'i': 3, 'v': 0.95}],
    'genuine': [{'i': 1, 'v': 0.225}, {'i': 2, 'v': 0.120}, {'i': 3, 'v': 0.058}],
    'note': 'swept: 2.8e6 Monte-Carlo escapes per rung at beta=2; the q=3 value is corroborated '
            'MC-free by the Fokker-Planck solve at +0.945. genuine: numerical integration.',
}

families = [
    {'key': 'swept', 'title': 'swept multicritical',
     'eqn': "u'' = ( sign(Y)|Y|^q  -  ηξ ) u", 'index': 'q  (turning order)',
     'members': ['fold → Tracy–Widom', 'cusp → \U0001d4b2', 'swallowtail → \U0001d4a2', 'butterfly', '…'],
     'integrable': False,
     'why': 'no Fredholm/soft-edge determinant (operator unbounded below); '
            'no Painlevé σ-reduction (skeleton is an isomonodromy fixed point)'},
    {'key': 'genuine', 'title': 'genuine diffraction catastrophe',
     'eqn': "u^(n) = ( Y  -  ηξ ) u", 'index': 'n  (order of the equation)',
     'members': ['Airy', 'Pearcey', 'swallowtail function', '…'],
     'integrable': True,
     'why': 'classical special functions of the unfolding'},
    {'key': 'rmt', 'title': 'random-matrix multicritical edge',
     'eqn': 'higher-order Tracy–Widom  (Fredholm determinant)', 'index': 'k  (order of vanishing)',
     'members': ['Airy kernel → TW', 'Painlevé I/II hierarchy kernels', '…'],
     'integrable': True,
     'why': 'Fredholm determinants of higher Airy-type kernels; Lax pairs throughout'},
]

ref = {'swept': swept, 'skew': skew, 'families': families, 'captions': CAPS['delta3']}

tpl = (HERE / 'coda.template.html').read_text()
(HERE / 'coda.html').write_text(tpl.replace('/*__REFERENCE__*/', json.dumps(ref)))

print('wrote coda.html')
print(f'  swept column computed at {len(QS)} rungs (ODE + root-find + quadrature)')
for s in swept:
    print(f"    q={s['q']}  node {s['node']:+.6f}  v0 {s['v0']:.5f}  "
          f"Bessel-{s['besselOrder']}  I(s)=s^{s['instExp']}/{s['instDen']}")
print('  assertions passed: q=1 node = first Airy zero; v_0^(2), v_0^(3) match the capstone')
print('  skew trends carried as MEASURED data, provenance on-screen')
