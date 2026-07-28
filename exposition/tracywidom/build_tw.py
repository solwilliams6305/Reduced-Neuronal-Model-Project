"""Build the Tracy-Widom explorable (module alpha5) -- the payoff, not a prerequisite.

PREREQ_GRAPH F4: Tracy-Widom is reachable with ZERO random matrix theory, via the stochastic
Airy operator.  So RMT is not a gateway that has to be taught first -- it is a REVEAL, and this
module is where it lands.

Prior-art verdict for this topic was PARTIAL: Quanta covers why TW matters better than we would,
but nothing found visualises the edge-scaling MECHANISM.  That is the gap this fills.

The mechanism is COLLAPSE, not slow convergence.  Raw largest-eigenvalue distributions at
different N look nothing alike, because the scale shrinks like N^(-2/3).  Rescale by that and
every N lands on one curve.  Notably the collapse is already excellent at N=25 -- so the story is
"a fixed shape hiding behind a shrinking scale", not "wait for N to get big".

Class B (here): the sampling.  Diagonalising thousands of matrices is not a browser activity.
Class A (browser): rendering, and the raw/rescaled toggle that is the whole point.

Normalisation note.  GUE with density ~ exp(-N tr H^2 / 2): off-diagonal (x+iy)/sqrt(2N) with the
symmetrisation dividing by sqrt(2), diagonal x/sqrt(N).  This puts the semicircle on [-2,2] and
makes N^(2/3)(lmax - 2) converge to TW_2.  An earlier draft got the factor wrong, the spectrum sat
on the wrong support, and the rescaled mean DIVERGED with N instead of converging -- so the
support is asserted below rather than assumed.

Run:  python3 build_tw.py   ->  tracywidom.html
"""
import json
import pathlib
import numpy as np

HERE = pathlib.Path(__file__).parent
CAPS = json.loads((HERE.parent / 'captions.json').read_text())
REQUIRED = ('headline', 'setup', 'raw', 'collapse', 'reveal', 'why',
            'distortion_beta2', 'distortion_cited', 'distortion_finite')
missing = [k for k in REQUIRED if k not in CAPS.get('alpha5', {})]
if missing:
    raise SystemExit(f'captions.json is missing alpha5 keys: {missing}')

# Tracy-Widom beta=2 reference moments (Bornemann, Tracy-Widom tables)
TW2 = {'mean': -1.771087, 'sd': 0.901773, 'skew': 0.224084}
rng = np.random.default_rng(20260728)


def gue(N):
    X = rng.standard_normal((N, N)); Y = rng.standard_normal((N, N))
    H = (X + 1j * Y) / np.sqrt(2 * N)
    H = (H + H.conj().T) / np.sqrt(2)
    np.fill_diagonal(H, rng.standard_normal(N) / np.sqrt(N))
    return H


# the semicircle support must be [-2, 2] or the rescaling below is meaningless
ev = np.linalg.eigvalsh(gue(400))
assert 1.90 < ev[-1] < 2.05 and -2.05 < ev[0] < -1.90, \
    f'GUE normalisation wrong: spectrum on [{ev[0]:.3f}, {ev[-1]:.3f}], want about [-2, 2]'

SIZES = [20, 40, 80, 160]
REPS = {20: 4000, 40: 3000, 80: 1600, 160: 700}

raw, scaled, stats = {}, {}, []
for N in SIZES:
    lm = np.array([np.linalg.eigvalsh(gue(N))[-1] for _ in range(REPS[N])])
    s = N ** (2 / 3) * (lm - 2.0)
    raw[N] = lm; scaled[N] = s
    sk = float(((s - s.mean()) ** 3).mean() / s.std() ** 3)
    stats.append({'N': N, 'reps': REPS[N],
                  'rawMean': float(lm.mean()), 'rawSd': float(lm.std()),
                  'mean': float(s.mean()), 'sd': float(s.std()), 'skew': sk})

# every size must agree with TW2 on mean and width -- that IS the collapse
for st in stats:
    assert abs(st['mean'] - TW2['mean']) < 0.12, \
        f"N={st['N']} rescaled mean {st['mean']:.3f} vs TW2 {TW2['mean']:.3f}"
    assert abs(st['sd'] - TW2['sd']) < 0.09, \
        f"N={st['N']} rescaled sd {st['sd']:.3f} vs TW2 {TW2['sd']:.3f}"

NB = 40
RAW_LO, RAW_HI = 1.0, 2.3
SC_LO, SC_HI = -5.0, 2.5
hraw = {N: np.histogram(raw[N], bins=NB, range=(RAW_LO, RAW_HI))[0].tolist() for N in SIZES}
hsc = {N: np.histogram(scaled[N], bins=NB, range=(SC_LO, SC_HI))[0].tolist() for N in SIZES}

# the semicircle, averaged over several matrices -- one matrix's spectrum is too noisy at 60
# bins to read as a semicircle at all, which defeats the panel's only job.
SPEC_N, SPEC_REPS = 400, 24
spec = np.concatenate([np.linalg.eigvalsh(gue(SPEC_N)) for _ in range(SPEC_REPS)])
sh, se = np.histogram(spec, bins=60, range=(-2.4, 2.4), density=True)
# exact Wigner semicircle for overlay: rho(x) = sqrt(4 - x^2) / (2 pi) on [-2, 2]
mid = 0.5 * (se[:-1] + se[1:])
theory = np.where(np.abs(mid) < 2.0, np.sqrt(np.maximum(4 - mid**2, 0)) / (2*np.pi), 0.0)
# the empirical density must track the exact one across the bulk
bulk = np.abs(mid) < 1.8
rel = float(np.max(np.abs(sh[bulk] - theory[bulk]) / theory[bulk]))
assert rel < 0.10, f'semicircle does not match theory in the bulk: {rel:.3f} relative'

ref = {
    'sizes': SIZES, 'stats': stats, 'tw2': TW2,
    'rawHist': hraw, 'scHist': hsc, 'NB': NB,
    'rawRange': [RAW_LO, RAW_HI], 'scRange': [SC_LO, SC_HI],
    'specHist': sh.tolist(), 'specEdges': se.tolist(), 'specTheory': theory.tolist(),
    'specN': SPEC_N, 'specReps': SPEC_REPS,
    'captions': CAPS['alpha5'],
}

tpl = (HERE / 'tracywidom.template.html').read_text()
(HERE / 'tracywidom.html').write_text(tpl.replace('/*__REFERENCE__*/', json.dumps(ref)))

print('wrote tracywidom.html')
print(f'  GUE support check: [{ev[0]:+.3f}, {ev[-1]:+.3f}]  (semicircle on [-2,2])')
print(f'  semicircle: {SPEC_REPS} x N={SPEC_N}; max bulk deviation from exact rho {rel:.3%}')
print(f"  TW2 target:  mean {TW2['mean']:+.4f}  sd {TW2['sd']:.4f}  skew {TW2['skew']:+.4f}")
for st in stats:
    print(f"    N={st['N']:4d} ({st['reps']:5d} reps)  raw sd {st['rawSd']:.5f}"
          f"   ->  rescaled mean {st['mean']:+.4f}  sd {st['sd']:.4f}  skew {st['skew']:+.4f}")
r = [st['rawSd'] for st in stats]
print(f'  raw width shrinks {r[0]/r[-1]:.1f}x from N={SIZES[0]} to N={SIZES[-1]};'
      f' rescaled widths agree to '
      f'{max(abs(st["sd"] - TW2["sd"]) for st in stats):.3f} of TW2 -- that is the collapse')
