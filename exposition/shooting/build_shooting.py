"""Build the shooting explorable (module alpha4c of DESIGN_inner_chain.md).

Closes the inner chain.  alpha3 showed blow-up producing the Riccati; alpha4a traded the
Riccati's blow-ups for the zeros of u; this shows those zeros ARE the spectrum, so the escape
levels and the eigenvalues are one set of numbers counted twice.

Architecture split, as in build_scrubber.py:
  Class B (here, offline) -- the MISS CURVE over a range of trial energies.  400 shoots is far
    too much to redo on every slider move, and it never changes, so it ships as data.
  Class A (browser, live) -- ONE shoot at the current lambda, a few milliseconds.
  Golden vectors bind the two: same RK4 scheme, checked to agree, not checked to be exact.

Run:  python3 build_shooting.py   ->  shooting.html
"""
import json
import pathlib
import numpy as np
from scipy.special import ai_zeros

HERE = pathlib.Path(__file__).parent
CAPS = json.loads((HERE.parent / 'captions.json').read_text())
REQUIRED = ('headline', 'setup', 'the_miss', 'the_join',
            'distortion_scale', 'distortion_halfline')
missing = [k for k in REQUIRED if k not in CAPS.get('alpha4c', {})]
if missing:
    raise SystemExit(f'captions.json is missing alpha4c keys: {missing}')

XMAX, N = 20.0, 4000
x = np.linspace(0.0, XMAX, N)
h = -(x[1] - x[0])                      # negative: we integrate DOWN from the decaying end


def shoot(lam):
    """u'' = (x - lam) u, seeded at XMAX with the WKB decaying ratio u'/u = -sqrt(x-lam).

    Integrated downward, so the recessive solution grows -- which is the numerically stable
    direction for it.  Overall scale is arbitrary: an eigenfunction is defined up to scale.
    """
    u = np.empty(N); v = np.empty(N)
    u[-1] = 1.0
    v[-1] = -np.sqrt(XMAX - lam)
    for i in range(N - 1, 0, -1):
        xi = x[i]
        k1u, k1v = v[i],               (xi - lam) * u[i]
        k2u, k2v = v[i] + h/2*k1v, (xi + h/2 - lam) * (u[i] + h/2*k1u)
        k3u, k3v = v[i] + h/2*k2v, (xi + h/2 - lam) * (u[i] + h/2*k2u)
        k4u, k4v = v[i] + h*k3v,   (xi + h   - lam) * (u[i] + h*k3u)
        u[i-1] = u[i] + h/6*(k1u + 2*k2u + 2*k3u + k4u)
        v[i-1] = v[i] + h/6*(k1v + 2*k2v + 2*k3v + k4v)
    return u, v


def miss(lam):
    u, _ = shoot(lam)
    return u[0] / np.max(np.abs(u))     # normalised, so the curve is comparable across lambda


LAMS = np.linspace(0.05, 14.0, 400)
MISS = np.array([miss(l) for l in LAMS])

# true eigenvalues: u(0)=0 with the decaying solution means Ai(-lam)=0, i.e. lam = |a_n|
eig = np.sort(-ai_zeros(8)[0])
eig = eig[eig < 14.0]

# how well does the shooting scheme actually land on them?
resid = [float(abs(miss(float(e)))) for e in eig]

idx = np.linspace(0, N - 1, 25).astype(int)
u_ref, v_ref = shoot(6.0)               # a lambda that is NOT an eigenvalue, so nothing is special
golden = [{'i': int(i), 'u': float(u_ref[i]), 'v': float(v_ref[i])} for i in idx]

ref = {
    'XMAX': XMAX, 'N': N,
    'lams': [float(v) for v in LAMS],
    'miss': [float(v) for v in MISS],
    'eig': [float(v) for v in eig],
    'goldenLam': 6.0,
    'golden': golden,
    'captions': CAPS['alpha4c'],
}

tpl = (HERE / 'shooting.template.html').read_text()
(HERE / 'shooting.html').write_text(tpl.replace('/*__REFERENCE__*/', json.dumps(ref)))

print('wrote shooting.html')
print(f'  eigenvalues (|a_n|) below 14 : {np.round(eig, 4)}')
print(f'  |miss| at each eigenvalue    : {[f"{r:.2e}" for r in resid]}')
print(f'  miss curve: {len(LAMS)} shoots precomputed (Class B); browser does 1 per slider move')
print(f'  golden vectors: {len(golden)} at lambda={ref["goldenLam"]}')
