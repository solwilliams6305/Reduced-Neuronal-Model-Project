"""Build the twin-panel scrubber  (gate 4: first web artifact).

Demonstrates the architecture from exposition/ARCHITECTURE.md in miniature:

  Class B (offline, Python)  ->  a versioned data artifact  ->  Class A (live, JS)

Python owns nothing the browser needs at interactive rate.  It emits three things:
  * the SEED -- Ai and Ai' at the right-hand endpoint, so the browser can integrate inward
    without needing a special-function library;
  * the Airy ZEROS in the window, for marking;
  * GOLDEN VECTORS -- u(Y) at check points from the *same* Euler scheme the JS will run.

The golden vectors test that the two implementations AGREE, not that either is exact.  That is
the distinction that matters: the drift we are guarding against is the hand-ported JS silently
diverging from the Python source of truth, which is the failure already present elsewhere in
this repo (see SURVEY_BRIEF.md).  Discretisation error against true Ai is a separate question and
is reported separately below.

Run:  python3 build_scrubber.py   ->  twin_scrubber.html
"""
import json
import numpy as np
from scipy.special import airy, ai_zeros

Y_LO, Y_HI, N = -9.0, 2.0, 30000
Y = np.linspace(Y_LO, Y_HI, N)
h = Y[1] - Y[0]
Ai, Aip, _, _ = airy(Y)

az = ai_zeros(12)[0]
zeros = np.sort(az[(az > Y_LO) & (az < Y_HI)])


def euler_deterministic():
    """The EXACT scheme the browser runs, at eta = 0. Must stay in lockstep with the JS."""
    u = np.empty(N); p = np.empty(N)
    u[-1], p[-1] = Ai[-1], Aip[-1]
    for i in range(N - 1, 0, -1):
        p[i-1] = p[i] - h * Y[i] * u[i]
        u[i-1] = u[i] - h * p[i]
    return u, p


u0, p0 = euler_deterministic()

# golden vectors at 25 evenly spaced check indices
idx = np.linspace(0, N - 1, 25).astype(int)
golden = [{'i': int(i), 'Y': float(Y[i]), 'u': float(u0[i]), 'p': float(p0[i])} for i in idx]

# honest report of the SCHEME's error against the true Airy function
scheme_err = float(np.max(np.abs(u0 - Ai)))
rel = scheme_err / float(np.max(np.abs(Ai)))

ref = {
    'Y_LO': Y_LO, 'Y_HI': Y_HI, 'N': N,
    'seed': {'u': float(Ai[-1]), 'p': float(Aip[-1])},
    'zeros': [float(z) for z in zeros],
    'golden': golden,
    'schemeAbsErrVsAiry': scheme_err,
    'schemeRelErrVsAiry': rel,
}

tpl = open('twin_scrubber.template.html').read()
out = tpl.replace('/*__REFERENCE__*/', json.dumps(ref))
open('twin_scrubber.html', 'w').write(out)

print(f'wrote twin_scrubber.html')
print(f'  zeros in window : {np.round(zeros, 4)}')
print(f'  seed at Y={Y_HI}: u={Ai[-1]:.6e}  p={Aip[-1]:.6e}')
print(f'  golden vectors  : {len(golden)}')
print(f'  scheme vs true Ai: abs {scheme_err:.3e}  rel {rel:.3e}   '
      f'({"first-order Euler, as expected" if rel < 5e-2 else "TOO LARGE — investigate"})')
