"""Build the fast-slow explorable (modules alpha1 + alpha2 merged).

Merged deliberately.  PREREQ_GRAPH F5 says blow-up must OPEN on the failure of normal
hyperbolicity, so something has to CLOSE on it -- and a module whose whole job is "here is where
the tools stop working" only makes sense once the tools have been shown working.  So: act one is
the geometry, act two is the number that kills it.

FitzHugh-Nagumo, planar:
    v' = v - v^3/3 - w + I        (fast)
    w' = eps (v + a - b w)        (slow)

Critical manifold (fast nullcline):  w = v - v^3/3 + I, a cubic.
Attraction rate along it:            lambda(v) = d/dv (v - v^3/3 - w + I) = 1 - v^2
  |v| > 1  ->  lambda < 0, attracting (the outer branches)
  |v| < 1  ->  lambda > 0, repelling  (the middle branch)
  v = +-1  ->  lambda = 0 EXACTLY.  Fenichel needs it bounded away from zero; here it is not.

That last line is the entire point of the module and it is exact, not asymptotic.

Run:  python3 build_fastslow.py   ->  fastslow.html
"""
import json
import pathlib
import numpy as np

HERE = pathlib.Path(__file__).parent
CAPS = json.loads((HERE.parent / 'captions.json').read_text())
REQUIRED = ('headline', 'setup', 'manifold', 'attraction', 'breakdown',
            'distortion_planar', 'distortion_rate')
missing = [k for k in REQUIRED if k not in CAPS.get('alpha12', {})]
if missing:
    raise SystemExit(f'captions.json is missing alpha12 keys: {missing}')

A, B, I = 0.7, 0.8, 0.5


def rhs(v, w, eps):
    return v - v ** 3 / 3.0 - w + I, eps * (v + A - B * w)


def integrate(eps, v0, w0, T, dt):
    """RK4.  Must stay in lockstep with the browser's integrator."""
    n = int(T / dt)
    V = np.empty(n + 1); W = np.empty(n + 1)
    V[0], W[0] = v0, w0
    v, w = v0, w0
    for i in range(n):
        k1 = rhs(v, w, eps)
        k2 = rhs(v + dt / 2 * k1[0], w + dt / 2 * k1[1], eps)
        k3 = rhs(v + dt / 2 * k2[0], w + dt / 2 * k2[1], eps)
        k4 = rhs(v + dt * k3[0], w + dt * k3[1], eps)
        v += dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        w += dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        V[i + 1], W[i + 1] = v, w
    return V, W


# does it actually oscillate at these parameters?  check, do not assume.
V, W = integrate(0.08, -1.6, -0.4, 400.0, 0.01)
tail = V[len(V) // 2:]
amp = float(tail.max() - tail.min())
assert amp > 1.0, f'FHN is not oscillating at I={I}: tail amplitude {amp:.3f}'

# the folds, and the attraction rate there
folds = [-1.0, 1.0]
lam = lambda v: 1.0 - v * v
assert all(abs(lam(f)) < 1e-15 for f in folds), 'fold is not where the attraction rate vanishes'

# golden vectors from a short deterministic run
gv, gw = integrate(0.08, -1.6, -0.4, 20.0, 0.01)
idx = np.linspace(0, len(gv) - 1, 21).astype(int)
golden = [{'i': int(i), 'v': float(gv[i]), 'w': float(gw[i])} for i in idx]

ref = {
    'A': A, 'B': B, 'I': I,
    'folds': folds,
    'golden': {'eps': 0.08, 'v0': -1.6, 'w0': -0.4, 'dt': 0.01, 'pts': golden},
    'captions': CAPS['alpha12'],
}

tpl = (HERE / 'fastslow.template.html').read_text()
(HERE / 'fastslow.html').write_text(tpl.replace('/*__REFERENCE__*/', json.dumps(ref)))

print('wrote fastslow.html')
print(f'  FHN  a={A} b={B} I={I}  -> oscillates, tail amplitude {amp:.3f}')
print(f'  critical manifold w = v - v^3/3 + {I};  attraction rate lambda(v) = 1 - v^2')
print(f'  folds at v = {folds},  lambda there = '
      f'{[float(lam(f)) for f in folds]}  <- exactly zero, which is the point')
for v in (-2.0, -1.5, -1.0, -0.5, 0.0):
    print(f'    lambda({v:+.1f}) = {lam(v):+.2f}  '
          f'{"attracting" if lam(v) < 0 else "repelling" if lam(v) > 0 else "NEITHER"}')
print(f'  golden vectors: {len(golden)}')
