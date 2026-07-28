"""Build the divergence explorable (module gamma1) -- the opening of the analysis thread.

PREREQ_GRAPH.md F1: this node depends on NOTHING else in the course, so the analysis thread can
start at episode 1 rather than waiting eight episodes for the dynamics.  This is that start.

The example is Euler's series, chosen for one reason: the function is exactly computable, so the
divergent series can be watched AGAINST THE TRUTH.  Almost nothing else in this course can be.

    E(x) = int_0^inf e^-t / (1 + x t) dt        exact, by quadrature
         ~ sum_k (-1)^k k! x^k                  asymptotic, radius of convergence ZERO

Class B (here): the Gauss-Laguerre RULE (nodes and weights), computed once.  Class A (browser):
applies the rule, and forms the partial sums.

An earlier version shipped E(x) tabulated on a grid and interpolated it in the browser.  That was
wrong in a way that mattered: linear interpolation carried ~5e-8 error, which is LARGER than the
optimal-truncation floor it is supposed to resolve, and the browser reported the wrong best-N at
small x.  The quadrature rule instead is exact to ~4e-13 and has no such failure mode.  The golden
vectors now cover E(x) as well as the partial sums, because the original set only checked the
sums and so missed this entirely.

Run:  python3 build_divergence.py   ->  divergence.html
"""
import json
import math
import pathlib
import numpy as np
from scipy.integrate import quad

HERE = pathlib.Path(__file__).parent
CAPS = json.loads((HERE.parent / 'captions.json').read_text())
REQUIRED = ('headline', 'setup', 'the_shock', 'the_floor', 'the_hook',
            'distortion_alternating', 'distortion_optimal')
missing = [k for k in REQUIRED if k not in CAPS.get('gamma1', {})]
if missing:
    raise SystemExit(f'captions.json is missing gamma1 keys: {missing}')

NMAX = 40


def E(x):
    return quad(lambda t: math.exp(-t) / (1.0 + x * t), 0.0, np.inf, limit=400)[0]


# The rule the browser will apply:  int_0^inf e^-t f(t) dt  ~  sum w_i f(t_i)
GL_N = 64
GL_T, GL_W = np.polynomial.laguerre.laggauss(GL_N)
E_gl = lambda x: float(np.sum(GL_W / (1.0 + x * GL_T)))

gl_worst = max(abs(E_gl(x) - E(x)) / E(x)
               for x in (0.02, 0.05, 0.12, 0.25, 0.45, 0.6))
assert gl_worst < 1e-10, f'Gauss-Laguerre rule not accurate enough: {gl_worst:.2e}'

# where the series is at its best, and how good that is
report = []
for x in (0.05, 0.12, 0.25, 0.45):
    terms = np.array([(-1.0) ** k * math.factorial(k) * x ** k for k in range(NMAX + 1)])
    partial = np.cumsum(terms)
    err = np.abs(partial - E(x))
    n_best = int(np.argmin(err))
    report.append({'x': x, 'nBest': n_best, 'minErr': float(err[n_best]),
                   'expMinusInvX': float(math.exp(-1.0 / x)),
                   'ratio': float(err[n_best] / math.exp(-1.0 / x))})

# golden vectors: the browser must reproduce these partial sums exactly
golden = []
for x in (0.08, 0.3):
    acc = 0.0
    for k in range(0, 13):
        acc += (-1.0) ** k * math.factorial(k) * x ** k
        if k in (0, 4, 8, 12):
            golden.append({'x': x, 'n': k, 'S': acc})

# golden vectors must cover E(x) too -- the original set checked only partial sums, which is
# exactly why the interpolation defect got through.
goldenE = [{'x': float(x), 'E': float(E(float(x)))}
           for x in (0.02, 0.05, 0.12, 0.25, 0.45, 0.6)]

ref = {'glT': [float(v) for v in GL_T], 'glW': [float(v) for v in GL_W],
       'NMAX': NMAX, 'golden': golden, 'goldenE': goldenE, 'captions': CAPS['gamma1']}

tpl = (HERE / 'divergence.template.html').read_text()
(HERE / 'divergence.html').write_text(tpl.replace('/*__REFERENCE__*/', json.dumps(ref)))

print('wrote divergence.html')
print(f'  Gauss-Laguerre rule shipped: {GL_N} nodes; worst rel error vs quad {gl_worst:.2e}')
print('  optimal truncation, and the floor against exp(-1/x):')
for r in report:
    print(f'    x={r["x"]:.2f}: best at N={r["nBest"]:2d}, min error {r["minErr"]:.3e}, '
          f'exp(-1/x)={r["expMinusInvX"]:.3e}, ratio {r["ratio"]:.2f}')
print(f'  golden: {len(golden)} partial sums + {len(goldenE)} values of E(x)')
# the headline claim: the floor tracks exp(-1/x) up to an algebraic factor
assert all(0.5 < r['ratio'] < 12 for r in report), 'floor no longer tracks exp(-1/x)'
