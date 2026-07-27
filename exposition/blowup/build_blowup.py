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

Run:  python3 build_blowup.py   ->  blowup.html
"""
import json
import pathlib

HERE = pathlib.Path(__file__).parent
CAPS = json.loads((HERE.parent / 'captions.json').read_text())

REQUIRED = ('headline', 'naive_fails', 'the_find', 'third_weight', 'the_real_point',
            'what_it_is', 'payoff', 'distortion_not_added', 'distortion_no_sphere')
missing = [k for k in REQUIRED if k not in CAPS.get('alpha3', {})]
if missing:
    raise SystemExit(f'captions.json is missing alpha3 keys: {missing}')

tpl = (HERE / 'blowup.template.html').read_text()
out = tpl.replace('/*__CAPTIONS__*/', json.dumps(CAPS['alpha3']))
(HERE / 'blowup.html').write_text(out)

print('wrote blowup.html')
print(f'  captions injected: {len(CAPS["alpha3"])} keys from ../captions.json')
print('  exact identities driving the widget:')
print('    manifold   Y = lambda^(p-2) X^2')
print('    dynamics   dX/dT = X^2 - lambda^(2-p) Y ,  dY/dT = -lambda^(1+p-q) epsbar')
print('    lambda leaves the equation iff  p = 2  and  q = 3')
