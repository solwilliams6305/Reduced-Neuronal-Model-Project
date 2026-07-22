"""Test that the FULL engine (quad + dense route + numpy DP) reproduces the pure-python
reference on the ACTUAL moments of v_k -- including cross-term and boundary moments, which the
earlier pure-bulk batteries did not cover. Any disagreement localizes the routing bug."""
import os, sys, time, pickle, random
os.environ.setdefault('CT_NUMBA', '0')
import numpy as np, sympy as sp
import chaos_diagram as CD
import chaos_transfer as CT

k = int(sys.argv[1]) if len(sys.argv) > 1 else 6
n = int(sys.argv[2]) if len(sys.argv) > 2 else 10
budget = float(sys.argv[3]) if len(sys.argv) > 3 else 38.0
t0 = time.time()
CD.setup(n_grid=n, MAXORD=CT.DENSE_MAXORD); CT.clear_all()
u0ps = CD.GEO['u0ps']; Vst = CD.GEO['Vst']
from _v6_driver import vk_terms
Us = [[sp.symbols(f'U{a}_{m}') for m in range(16)] for a in range(16)]
s = [sp.symbols(f's{j}') for j in range(14)]; Vs = [sp.symbols(f'Vd{j}') for j in range(14)]
Y = {idx: CD.Ybase(idx, Us, s, Vs, None, u0ps, Vst)
     for idx in sorted(set(p for pr in vk_terms(k) for p in pr[:2]))}
# load the needed-moment list if present, else enumerate
MOMF = f'_ckpt_{k}_{n}_moments.pkl'
if os.path.exists(MOMF):
    needed = pickle.load(open(MOMF, 'rb'))
else:
    needed = set()
    for a, b, _ in vk_terms(k):
        for c, at in Y[a]:
            if c != 0.0: needed.add(tuple(sorted(at)))
        for c, at in Y[b]:
            if c != 0.0: needed.add(tuple(sorted(at)))
        for ca, aa in Y[a]:
            for cb, ab in Y[b]:
                if ca and cb: needed.add(tuple(sorted(aa + ab)))
    needed = list(needed)
print(f'{len(needed)} needed; sampling & comparing full-engine vs pure-python', flush=True)

def pure(atoms):
    # reference = numpy-DP path only (validated == pure-python, but fast): quad OFF, dense route OFF
    _q, _d = CT.QUAD_ENABLE, CT.DENSE_ROUTE
    CT.QUAD_ENABLE = False; CT.DENSE_ROUTE = False
    CT._HMEMO.clear(); CT._TMEMO.clear(); CT._RMEMO.clear(); CD._MCACHE.clear()
    v = CT.moment_hybrid(list(atoms))
    CT.QUAD_ENABLE, CT.DENSE_ROUTE = _q, _d
    CT._HMEMO.clear(); CT._TMEMO.clear(); CT._RMEMO.clear(); CD._MCACHE.clear()
    return v
def full(atoms):
    CT._HMEMO.clear(); CT._TMEMO.clear(); CT._RMEMO.clear(); CD._MCACHE.clear()
    return CT.moment_hybrid(list(atoms))

random.seed(1); random.shuffle(needed)
maxre = 0.0; nchk = 0; worst = None
for m in needed:
    if time.time() - t0 > budget: break
    vf = full(m); vp = pure(m)
    d = abs(vf - vp); sc = max(abs(vf), abs(vp), 1e-300); re = d / sc
    nchk += 1
    if re > maxre:
        maxre = re; worst = (m, vf, vp, re)
print(f'[{time.time()-t0:.0f}s] checked {nchk} moments; max relerr {maxre:.2e}', flush=True)
if worst and maxre > 1e-9:
    m, vf, vp, re = worst
    bulk = tuple(sorted(a[1] for a in m if a[0]=='U')); nb = sum(1 for a in m if a[0]=='s')
    print(f'  WORST: bulk={bulk} nb={nb}  full={vf:+.6e} pure={vp:+.6e} re={re:.2e}', flush=True)
    print(f'  atoms={sorted(m)}', flush=True)
else:
    print('  ALL SAMPLED AGREE (<=1e-9) -> routing is faithful; discrepancy is elsewhere', flush=True)
