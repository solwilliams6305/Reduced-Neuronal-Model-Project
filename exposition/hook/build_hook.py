"""Build the opening hook (module ep0).

Its only job is to make two questions unavoidable, and then get out of the way.  It explains
nothing: the whole course is the answer.

FitzHugh-Nagumo in the EXCITABLE regime -- I chosen so the deterministic system is SILENT.  That
matters: every spike you see is caused by the noise, so the randomness is not decoration on top
of an oscillation, it is the mechanism.  Verified here rather than assumed.

Class B (here): the long reference run -- 10^5-ish spikes for stable moments, far too slow to do
in a browser.  Class A (browser): a live run the viewer watches accumulate toward it.

Run:  python3 build_hook.py   ->  hook.html
"""
import json
import pathlib
import numpy as np

HERE = pathlib.Path(__file__).parent
CAPS = json.loads((HERE.parent / 'captions.json').read_text())
REQUIRED = ('headline', 'setup', 'q1', 'q2', 'promise',
            'distortion_notreal', 'distortion_isi')
missing = [k for k in REQUIRED if k not in CAPS.get('ep0', {})]
if missing:
    raise SystemExit(f'captions.json is missing ep0 keys: {missing}')

A, B, I, EPS, DT, VTHRESH, VRESET = 0.7, 0.8, 0.30, 0.08, 0.01, 1.0, -1.0
# VRESET gives the detector HYSTERESIS (a Schmitt trigger).  Without it the noisy voltage
# wobbles across the threshold and each wobble counts: at sigma=0.10 that inflated the spike
# count by ~3x and put a large spurious spike at gap ~ 0 in the histogram.  Measured: 1611
# crossings of which 658 had gaps under one time unit; with hysteresis, 563 spikes and zero.


def run(sigma, T, seed, v0=-1.2, w0=-0.6):
    rng = np.random.default_rng(seed)
    n = int(T / DT); sq = np.sqrt(DT)
    v, w, prev = v0, w0, v0
    times = []
    armed = True
    for i in range(n):
        v += (v - v ** 3 / 3 - w + I) * DT + sigma * sq * rng.standard_normal()
        w += EPS * (v + A - B * w) * DT
        if armed and prev < VTHRESH <= v:
            times.append(i * DT)
            armed = False
        if v < VRESET:
            armed = True
        prev = v
    return np.array(times)


# EXCITABLE: silent without noise.  (One crossing is the initial transient, not a spike train.)
quiet = run(0.0, 4000.0, seed=1)
assert len(quiet) <= 1, f'not excitable at I={I}: {len(quiet)} spikes with no noise'
loud = run(0.10, 4000.0, seed=1)
# NB this bound was 100 when the detector had no hysteresis and was counting threshold jitter.
# With real spikes only, mean ISI ~53, so T=4000 yields ~75.  The assertion was calibrated
# against an inflated count -- corrected here rather than relaxed silently.
assert len(loud) > 50, f'noise does not drive spiking at I={I}: only {len(loud)}'

# the long reference the browser's live histogram should approach
ref_isi = []
for s in range(6):
    t = run(0.10, 30000.0, seed=100 + s)
    ref_isi.append(np.diff(t))
ISI = np.concatenate(ref_isi)
m, sd = float(ISI.mean()), float(ISI.std())
skew = float(((ISI - m) ** 3).mean() / sd ** 3)
HMAX = 260.0
hist, edges = np.histogram(ISI, bins=44, range=(0.0, HMAX))

# golden vectors: the deterministic (sigma = 0) trajectory
gv, gw, out = -1.2, -0.6, []
for i in range(2001):
    if i % 200 == 0:
        out.append({'i': i, 'v': float(gv), 'w': float(gw)})
    gv2 = gv + (gv - gv ** 3 / 3 - gw + I) * DT
    gw = gw + EPS * (gv + A - B * gw) * DT
    gv = gv2

ref = {
    'A': A, 'B': B, 'I': I, 'eps': EPS, 'dt': DT, 'vth': VTHRESH, 'vreset': VRESET,
    'refHist': [int(h) for h in hist], 'refEdges': [float(e) for e in edges],
    'refStats': {'n': int(ISI.size), 'mean': m, 'cv': sd / m, 'skew': skew},
    'golden': out,
    'captions': CAPS['ep0'],
}

tpl = (HERE / 'hook.template.html').read_text()
(HERE / 'hook.html').write_text(tpl.replace('/*__REFERENCE__*/', json.dumps(ref)))

print('wrote hook.html')
print(f'  regime check at I={I}: sigma=0 -> {len(quiet)} crossing(s) [EXCITABLE], '
      f'sigma=0.10 -> {len(loud)} spikes')
print(f'  reference ISI: {ISI.size} intervals from 6 runs of T=30000')
print(f'    mean {m:.2f}   CV {sd/m:.3f}   skew {skew:+.3f}')
print(f'  golden vectors: {len(out)}')
assert skew > 0.2, f'reference ISI is not visibly right-skewed (skew {skew:.3f})'
# no spurious near-zero gaps may survive the detector
assert (ISI < 1.0).sum() == 0, f'{(ISI < 1.0).sum()} sub-unit gaps: hysteresis is not working'
print(f'  detector: hysteresis on, {int((ISI < 1.0).sum())} gaps below 1.0 time unit')
