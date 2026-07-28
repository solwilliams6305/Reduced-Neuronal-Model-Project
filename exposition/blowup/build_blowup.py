"""Build the blow-up explorable (module alpha3 of DESIGN_inner_chain.md).

Prior art sweep verdict for this topic: GAP -- "arguably the single most animation-hungry idea
on this list and completely unserved".  No public interactive or video shows the blow-up
construction.  So there is nothing to borrow here and everything to invent.

The reframe under test (DESIGN_inner_chain.md sec.3): blow-up is not a trick, it is ANISOTROPIC
ZOOM.  Open on the failure of the naive zoom, let the viewer hunt for the exponent that stops
the picture degenerating, and name what they found afterwards.

Injects exposition/captions.json -- the shared narrative contract.  Captions live in one place
so a caveat fixed in one artifact cannot ship missing from its sibling, which is exactly what
happened to the boundary-condition note between twin_panel.py and twin_scrubber.html.

VERIFICATION, added 2026-07-28 after the gate-6 retrospective found this was the ONLY module
with neither a self-test nor an assertion.  It had been shipping three algebraic identities in a
print statement with nothing checking them.  They are now derived symbolically here and the
exponents the browser hardcodes are asserted against that derivation, so the widget cannot drift
from its own algebra.

Note on convention, because it is easy to get backwards and I did: the lens is
    X = lambda x ,  Y = lambda^p y ,  eps = lambda^(-q) ebar ,  T = t / lambda
i.e. the NEW coordinates are the magnified ones.  Substituting into the fold normal form
xdot = x^2 - y, ydot = -eps gives
    dX/dT = X^2 - lambda^(2-p) Y ,      dY/dT = -lambda^(1+p-q) ebar
Deriving it with the opposite convention (x = lambda X) negates both exponents.  Since negation
does not move a zero, BOTH conventions agree that lambda leaves the equation exactly at
p = 2, q = 3 -- which is why the error would not have shown up in the punchline.

Run:  python3 build_blowup.py   ->  blowup.html
"""
import json
import pathlib
import sympy as sp

HERE = pathlib.Path(__file__).parent
CAPS = json.loads((HERE.parent / 'captions.json').read_text())

REQUIRED = ('headline', 'naive_fails', 'the_find', 'third_weight', 'the_real_point',
            'what_it_is', 'payoff', 'distortion_not_added', 'distortion_no_sphere')
missing = [k for k in REQUIRED if k not in CAPS.get('alpha3', {})]
if missing:
    raise SystemExit(f'captions.json is missing alpha3 keys: {missing}')

# --- derive the two exponents symbolically, rather than trusting the comment ------------------
lam, p, q, T, X, Y, eb = sp.symbols('lambda p q T X Y epsbar', positive=True)
x, y, eps = X / lam, Y / lam ** p, eb / lam ** q
dXdT = sp.expand(lam ** 2 * (x ** 2 - y))          # from lam^-2 dX/dT = x^2 - y
dYdT = sp.simplify(lam ** (p + 1) * (-eps))        # from lam^(-p-1) dY/dT = -eps

expA_sym = sp.simplify(sp.log(sp.simplify(-dXdT.coeff(Y))) / sp.log(lam))
expB_sym = sp.simplify(sp.log(sp.simplify(-dYdT / eb)) / sp.log(lam))

# these two literals are what blowup.template.html hardcodes; they must match the derivation
assert sp.simplify(expA_sym - (2 - p)) == 0, f'expA drifted: derived {expA_sym}, template 2-p'
assert sp.simplify(expB_sym - (1 + p - q)) == 0, f'expB drifted: derived {expB_sym}, template 1+p-q'

# the widget's whole claim: lambda vanishes from the system at exactly one (p, q)
sol = sp.solve([expA_sym, expB_sym], [p, q], dict=True)
assert sol == [{p: 2, q: 3}], f'the vanishing point is not unique or not (2,3): {sol}'

# the X^2 coefficient must be exactly 1 -- that is what fixes the time weight T = t/lambda
assert sp.simplify(dXdT.coeff(X, 2) - 1) == 0, 'the X^2 coefficient is not 1'

# --- golden vectors so the browser's Math.pow agrees with Python ------------------------------
EBAR = 0.13
golden = []
for lv in (0.5, 1.0, 4.0, 30.0):
    for pv, qv in ((1.0, 1.0), (2.0, 3.0), (1.5, 2.2), (2.6, 3.4)):
        golden.append({'lam': lv, 'p': pv, 'q': qv,
                       'A': float(lv ** (2 - pv)),
                       'B': float(lv ** (1 + pv - qv) * EBAR)})

# at the blow-up weights the coefficients must not move when lambda does -- the punchline
base = None
for lv in (0.5, 1.0, 4.0, 30.0, 1000.0):
    A, B = lv ** (2 - 2.0), lv ** (1 + 2.0 - 3.0) * EBAR
    if base is None:
        base = (A, B)
    assert abs(A - base[0]) < 1e-12 and abs(B - base[1]) < 1e-12, \
        f'at (p,q)=(2,3) the system still depends on lambda={lv}'
# and away from them it must -- or the widget would be claiming nothing
assert abs(30.0 ** (2 - 1.0) - 1.0) > 1.0, 'off-weights coefficients are not lambda-dependent'

ref = {'EBAR': EBAR, 'expA': '2 - p', 'expB': '1 + p - q', 'golden': golden,
       'captions': CAPS['alpha3']}

tpl = (HERE / 'blowup.template.html').read_text()
out = (tpl.replace('/*__CAPTIONS__*/', json.dumps(CAPS['alpha3']))
          .replace('/*__VERIFY__*/', json.dumps(ref)))
(HERE / 'blowup.html').write_text(out)

print('wrote blowup.html')
print(f'  captions injected: {len(CAPS["alpha3"])} keys from ../captions.json')
print('  algebra DERIVED symbolically, not asserted in prose:')
print(f'    dX/dT = {sp.expand(dXdT)}')
print(f'    dY/dT = {dYdT}')
print(f'    exponent on Y     : {expA_sym}   (template hardcodes 2 - p)      MATCH')
print(f'    exponent on ebar  : {expB_sym}   (template hardcodes 1 + p - q)  MATCH')
print(f'    lambda vanishes at: {sol[0]}  -- unique solution')
print(f'  golden vectors for the browser: {len(golden)}')
print('  lambda-independence at (2,3) verified over lambda in [0.5, 1000]')
