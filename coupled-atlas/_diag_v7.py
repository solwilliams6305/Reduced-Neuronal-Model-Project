"""Diagnostic: reproduce + localize the v7 wall.  Categorize Var(Y8) monomial-pair
moment classes and time representative hard moments through the current engine."""
import os, sys, time
os.environ.setdefault('CT_NUMBA', '1')
import numpy as np
import chaos_diagram as CD
import chaos_transfer as CT
from collections import Counter

t0 = time.time()
n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
CD.setup(n_grid=n, MAXORD=CT.DENSE_MAXORD); CT.clear_all()
print(f'[{time.time()-t0:.1f}s] setup done, n={CD.GEO["n"]}', flush=True)
print(f'  numba available: {CT._get_numba_bulk() is not None}', flush=True)
u0ps = CD.GEO['u0ps']; Vst = CD.GEO['Vst']
# ybase pkls are all cached -> skip the slow 34MB sympify parse entirely (pass Ye=None).
import sympy as sp
Us = [[sp.symbols(f'U{k}_{m}') for m in range(16)] for k in range(16)]
s = [sp.symbols(f's{j}') for j in range(14)]; Vs = [sp.symbols(f'Vd{j}') for j in range(14)]; Ye = None
print(f'[{time.time()-t0:.1f}s] symbols ready (skipped file parse)', flush=True)
Y8 = CD.Ybase(8, Us, s, Vs, Ye, u0ps, Vst)
print(f'[{time.time()-t0:.1f}s] Y8 built: {len(Y8)} monomials', flush=True)

# Enumerate the DISTINCT product-moment atom-multisets for Var(Y8) = <Y8 Y8> - <Y8>^2
# and categorize by (sorted bulk orders, #boundary).
def sig(atoms):
    bulk = tuple(sorted(a[1] for a in atoms if a[0] == 'U'))
    nb = sum(1 for a in atoms if a[0] == 's')
    return (bulk, nb)

classcnt = Counter()
distinct = {}
for ca, aa in Y8:
    for cb, ab in Y8:
        if ca == 0.0 or cb == 0.0:
            continue
        atoms = aa + ab
        key = tuple(sorted(atoms))
        if key in distinct:
            continue
        distinct[key] = atoms
        classcnt[sig(atoms)] += 1
print(f'[{time.time()-t0:.1f}s] Var(Y8): {len(distinct)} distinct moments, {len(classcnt)} classes', flush=True)

# Characterize the "hard" classes: pure bulk, >=3 chains (order>=2), with a chain order>=3,
# NOT all-order<=2 (so quad path does NOT apply).
def is_hard(bulk, nb):
    if nb > 0:
        return False
    chains = [o for o in bulk if o >= 2]
    if all(o <= 2 for o in bulk):   # quad-eligible
        return False
    if len(chains) >= 3:
        return True
    return False

hard = [(sg, c) for sg, c in classcnt.items() if is_hard(*sg)]
hard.sort(key=lambda kv: (-len([o for o in kv[0][0] if o >= 2]), -sum(kv[0][0])))
print(f'[{time.time()-t0:.1f}s] HARD classes (pure-bulk, >=3 chains, not all<=2): {len(hard)}', flush=True)
tot_hard = sum(c for _, c in hard)
print(f'  total distinct hard moments: {tot_hard}', flush=True)
for sg, c in hard[:20]:
    nchain = len([o for o in sg[0] if o >= 2])
    print(f'   {sg}  nchains={nchain} x{c}', flush=True)

# also count quad-eligible pure-bulk (all order<=2) and 2-chain
quad_cnt = sum(c for sg, c in classcnt.items() if sg[1] == 0 and all(o <= 2 for o in sg[0]))
twochain = sum(c for sg, c in classcnt.items() if sg[1] == 0 and len([o for o in sg[0] if o >= 2]) == 2 and not all(o <= 2 for o in sg[0]))
bnd_cnt = sum(c for sg, c in classcnt.items() if sg[1] > 0)
print(f'  quad-eligible(all<=2): {quad_cnt}   two-chain(>2): {twochain}   boundary: {bnd_cnt}', flush=True)
