"""Build the Stokes explorable (module gamma8).

Prior-art verdict for this topic: GAP.  Nothing found that animates a Stokes crossing or shows
the error-function smoothing; the best public material is prose with static figures.

The claim under test is the one the distortion ledger has carried since the first day of this
project: a Stokes line drawn as a sharp boundary implies the coefficient jumps discontinuously,
and it does not.  So the module is built around the contrast -- textbook step against Berry's
error function -- with the radius in the viewer's hand so they can watch the ramp sharpen and
see WHY the step gets drawn.

Two things are exact here and are computed, not asserted:
  * the ray angles.  zeta = (2/3) z^(3/2); Im zeta = 0 at arg z = 0, +-2pi/3 (Stokes, Berry's
    convention) and Re zeta = 0 at arg z = +-pi/3, pi (anti-Stokes).
  * the Stokes constant.  Ai(z) + w Ai(wz) + wbar Ai(wbar z) = 0 with w = exp(2 pi i/3) holds
    identically; it is what fixes the constant to i.  Verified below to machine precision.

What is NOT exact is the erf profile itself: it is the leading term of Berry's smoothing,
asymptotic in |zeta|.  Said on the page.

Run:  python3 build_stokes.py   ->  stokes.html
"""
import json
import pathlib
import numpy as np
from scipy.special import airy, erf

HERE = pathlib.Path(__file__).parent
CAPS = json.loads((HERE.parent / 'captions.json').read_text())
REQUIRED = ('headline', 'setup', 'textbook', 'berry', 'what_switches', 'constant',
            'distortion_convention', 'distortion_erf')
missing = [k for k in REQUIRED if k not in CAPS.get('gamma8', {})]
if missing:
    raise SystemExit(f'captions.json is missing gamma8 keys: {missing}')

# ---- the connection identity that fixes the Stokes constant -------------------------------
w = np.exp(2j * np.pi / 3)
probe = [1.4 + 0.6j, -2.1 + 1.3j, 0.3 - 2.2j, 3.0 + 0.0j, -1.0 - 1.0j]
resid = []
for z in probe:
    a0 = airy(z)[0]
    a1 = airy(w * z)[0]
    a2 = airy(np.conj(w) * z)[0]
    resid.append(abs(a0 + w * a1 + np.conj(w) * a2))
worst = max(resid)

# ---- exact ray angles ----------------------------------------------------------------------
zeta = lambda z: (2.0 / 3.0) * z ** 1.5
stokes_rays = [0.0, 2 * np.pi / 3, -2 * np.pi / 3]          # Im zeta = 0
anti_rays = [np.pi / 3, -np.pi / 3, np.pi]                  # Re zeta = 0

# sanity: check the classification numerically rather than trusting the algebra
chk = {}
for th in stokes_rays:
    chk[f'Im zeta at {np.degrees(th):+.0f} deg'] = float(abs(zeta(2.0 * np.exp(1j * th)).imag))
for th in anti_rays[:2]:
    chk[f'Re zeta at {np.degrees(th):+.0f} deg'] = float(abs(zeta(2.0 * np.exp(1j * th)).real))


def sigma(r, th):
    """Berry's natural variable across the Stokes line at arg z = 0."""
    zt = zeta(r * np.exp(1j * th))
    return zt.imag / np.sqrt(2.0 * zt.real)


# golden samples so the browser's sigma/erf agree with these
golden = []
for r in (1.5, 8.0, 60.0):
    for th in (-0.5, -0.12, 0.0, 0.12, 0.5):
        s = float(sigma(r, th))
        golden.append({'r': r, 'th': th, 'sigma': s, 'S': float(0.5 * (1 + erf(s)))})

# angular width where |sigma| <= 1, as a function of radius
widths = []
for r in np.geomspace(1.0, 300.0, 60):
    th = np.linspace(1e-6, np.pi / 3 - 1e-3, 4000)
    s = np.array([sigma(r, t) for t in th])
    k = np.argmin(np.abs(s - 1.0))
    widths.append({'r': float(r), 'halfWidth': float(th[k])})

ref = {
    'stokesRays': [float(t) for t in stokes_rays],
    'antiRays': [float(t) for t in anti_rays],
    'golden': golden,
    'widths': widths,
    'connectionResidual': float(worst),
    'captions': CAPS['gamma8'],
}

tpl = (HERE / 'stokes.template.html').read_text()
(HERE / 'stokes.html').write_text(tpl.replace('/*__REFERENCE__*/', json.dumps(ref)))

print('wrote stokes.html')
print(f'  connection identity  Ai(z) + w Ai(wz) + wbar Ai(wbar z) = 0')
print(f'    worst |residual| over {len(probe)} probe points : {worst:.3e}   '
      f'=> Stokes constant is exactly i')
for k, v in chk.items():
    print(f'  {k:28s} : {v:.2e}')
print(f'  golden sigma/erf samples : {len(golden)}')
print(f'  half-width at r=1 : {widths[0]["halfWidth"]:.4f} rad'
      f'   at r=300 : {widths[-1]["halfWidth"]:.4f} rad'
      f'   (ratio {widths[0]["halfWidth"]/widths[-1]["halfWidth"]:.1f}x)')
