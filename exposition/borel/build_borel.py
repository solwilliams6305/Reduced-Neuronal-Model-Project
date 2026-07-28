"""Build the Borel explorable (modules gamma2 + gamma3) -- the linchpin, promoted from stills.

This was the original pilot.  It passed gate 3 as a storyboard and then never got promoted,
so the module carrying the single highest-value beat in the course was the one without a live
artifact.  This closes that.

PREREQ_GRAPH F6: the chain is  Borel transform has a cut -> the sum depends which side you pass
-> THAT AMBIGUITY IS EXACTLY e^{-A/x}.  Everything downstream (trans-series, Stokes, resurgence)
is bookkeeping on that one identity, so it is the beat most worth making tactile rather than
asserted.

What makes an interactive possible here: for the model series m_k = k!/A^(k+1) the Borel
transform is EXACTLY

    B(zeta) = sum_k m_k zeta^k / k! = sum_k zeta^k / A^(k+1) = 1 / (A - zeta)

-- a simple pole at zeta = A.  So the Laplace contour integral and its ambiguity are both exactly
computable, and the viewer can drag the contour around the pole and watch the answer change.  The
storyboard's panels 7-8 were schematic; here they are computed.

Two facts the widget lets you discover by hand:
  * the value does NOT depend on the detour height -- Cauchy;
  * it DOES depend on the side, and the difference is exactly 2 pi i exp(-A/x).

The real cusp coefficients are kept for the bridge panel, because the honest motivation for going
to the Borel plane at all is that seven coefficients cannot show the least-term law directly.

Run:  python3 build_borel.py   ->  borel.html
"""
import json
import math
import pathlib
import numpy as np

HERE = pathlib.Path(__file__).parent
CAPS = json.loads((HERE.parent / 'captions.json').read_text())
REQUIRED = ('headline', 'setup', 'cut', 'twochoices', 'identity', 'sowhat',
            'distortion_model', 'distortion_pole', 'distortion_quadrature')
missing = [k for k in REQUIRED if k not in CAPS.get('gamma23', {})]
if missing:
    raise SystemExit(f'captions.json is missing gamma23 keys: {missing}')

A = 1.9                    # the real instanton, S = s^5/10, as in the cusp paper
T, NQ, WID = 45.0, 40000, 0.55    # the browser uses exactly these


def contour_value(x, h, n=NQ, T=T, w=WID):
    """Laplace integral of B along a contour that detours over/under the pole by height h."""
    t = np.linspace(1e-9, T, n)
    bump = np.exp(-((t - A) / w) ** 2)
    z = t + 1j * h * bump
    dz = 1 + 1j * h * bump * (-2 * (t - A) / w ** 2)
    return complex(np.trapz(np.exp(-z / x) / (A - z) * dz, t))


# --- the physics, checked at HIGH resolution (independent of what the browser does) ----------
phys = []
for x in (0.25, 0.40, 0.60):
    up = contour_value(x, +0.9, n=400000)
    dn = contour_value(x, -0.9, n=400000)
    pred = 2j * math.pi * math.exp(-A / x)
    rel = abs((up - dn) - pred) / abs(pred)
    phys.append({'x': x, 'rel': rel})
    assert rel < 5e-3, f'ambiguity at x={x} misses the residue prediction: {rel:.2e}'

# Cauchy: the answer must not care how high the detour goes
base = contour_value(0.40, 0.6, n=400000)
for h in (0.9, 1.4, 2.0):
    d = abs(contour_value(0.40, h, n=400000) - base) / abs(base)
    assert d < 2e-4, f'value depends on detour height h={h}: {d:.2e}'

# --- golden vectors at the resolution the BROWSER uses, so the port test is exact ------------
golden = [{'x': x, 'h': h,
           're': contour_value(x, h).real, 'im': contour_value(x, h).imag}
          for x in (0.30, 0.50) for h in (-1.1, -0.7, 0.7, 1.1)]

# the real cusp ladder, for the bridge panel
v = [0.134, 0.111, 0.104, -0.030, -0.451, -1.19, -1.90]
b = [v[k] / math.factorial(k) for k in range(len(v))]

ref = {
    'A': A, 'T': T, 'NQ': NQ, 'WID': WID,
    'golden': golden, 'phys': phys,
    'v': v, 'borelCoef': b,
    'captions': CAPS['gamma23'],
}

tpl = (HERE / 'borel.template.html').read_text()
(HERE / 'borel.html').write_text(tpl.replace('/*__REFERENCE__*/', json.dumps(ref)))

print('wrote borel.html')
print(f'  B(zeta) = 1/({A} - zeta) exactly, for the model series m_k = k!/A^(k+1)')
print('  ambiguity (above - below) vs 2 pi i exp(-A/x), at high resolution:')
for p in phys:
    print(f'    x={p["x"]:.2f}: relative error {p["rel"]:.2e}')
print(f'  Cauchy check: value independent of detour height to < 2e-4')
print(f'  golden vectors at browser resolution (n={NQ}): {len(golden)}')
