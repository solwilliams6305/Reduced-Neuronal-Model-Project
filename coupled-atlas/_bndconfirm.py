"""Confirm the v6(n=10) anchor = the DENSE boundary convention. Reuse the BND='transfer' cache
(all 9722 moments) for everything EXCEPT the maxk<=5 boundary moments, which we recompute with
BND='dense' (beval). Bulk + high-order-boundary moments are identical between conventions, so this
isolates the convention shift cheaply. Checkpointed dense-boundary cache."""
import os, sys, time, pickle
os.environ.setdefault('CT_NUMBA', '0')
import numpy as np, sympy as sp
import chaos_diagram as CD
import chaos_transfer as CT
from _v6_driver import vk_terms

k, n = 6, 10
budget = float(sys.argv[1]) if len(sys.argv) > 1 else 36.0
t0 = time.time()
CD.setup(n_grid=n, MAXORD=CT.DENSE_MAXORD); CT.clear_all()
u0ps = CD.GEO['u0ps']; Vst = CD.GEO['Vst']
Us = [[sp.symbols(f'U{a}_{m}') for m in range(16)] for a in range(16)]
s = [sp.symbols(f's{j}') for j in range(14)]; Vs = [sp.symbols(f'Vd{j}') for j in range(14)]
Y = {idx: CD.Ybase(idx, Us, s, Vs, None, u0ps, Vst)
     for idx in sorted(set(p for pr in vk_terms(k) for p in pr[:2]))}
needed = pickle.load(open(f'_ckpt_{k}_{n}_moments.pkl', 'rb'))
tcache = pickle.load(open(f'_ckpt_{k}_{n}_cache.pkl', 'rb'))     # BND=transfer values (bulk identical)

def is_lowbnd(m):
    # only j>=2 boundary legs at order<=5 differ between transfer and dense; everything else is
    # identical between conventions and already correct in the transfer cache.
    js = [a[1] for a in m if a[0] == 's']
    if not js or max(js) < 2:
        return False
    maxk = max((a[1] for a in m if a[0] == 'U'), default=0)
    return maxk <= CT.DENSE_MAXORD

lowbnd = [m for m in needed if is_lowbnd(m)]
DCF = f'_bndconfirm_{k}_{n}_cache.pkl'
dcache = pickle.load(open(DCF, 'rb')) if os.path.exists(DCF) else {}
print(f'[{time.time()-t0:.0f}s] {len(lowbnd)} maxk<=5 boundary moments; {len(dcache)} recomputed', flush=True)

CT.BND_ENGINE = 'dense'
# Do NOT clear caches per moment -- but do NOT pre-seed with the transfer cache either (those hold
# the WRONG value for the j>=2 moments we are recomputing).  Each j>=2 moment is computed fresh via
# CD.moment (dense/beval); sub-caches accumulate harmlessly.
todo = [m for m in lowbnd if m not in dcache]
i = 0
for m in todo:
    t1 = time.time()
    dcache[m] = CT.moment_hybrid(list(m))     # dense (beval) -- must be exact, no fallback
    dt = time.time() - t1
    if dt > 3:
        bulk = tuple(sorted(a[1] for a in m if a[0]=='U')); js = sorted(a[1] for a in m if a[0]=='s')
        print(f'    slow {dt:.1f}s: bulk={bulk} s={js}', flush=True)
    i += 1
    if time.time() - t0 > budget: break
    # keep memory bounded on the 3.8 GB box
    if len(CT._HMEMO) > 400000:
        CT._HMEMO.clear(); CT._TMEMO.clear(); CT._RMEMO.clear(); CD._MCACHE.clear(); CT._HMEMO.update(tcache)
pickle.dump(dcache, open(DCF, 'wb'))
rem = sum(1 for m in lowbnd if m not in dcache)
print(f'[{time.time()-t0:.0f}s] recomputed {i} this run; {len(dcache)}/{len(lowbnd)}; {rem} remaining', flush=True)

if rem == 0:
    def mom(at):
        key = tuple(sorted(at))
        return dcache[key] if key in dcache else tcache[key]   # dense for low-bnd, transfer else
    def pm(A, B):
        tot = 0.0
        for ca, aa in A:
            for cb, ab in B:
                if ca == 0.0 or cb == 0.0: continue
                tot += ca * cb * mom(aa + ab)
        return tot
    def mean(A): return sum(c * mom(at) for c, at in A)
    tot = 0.0
    for a, b, coef in vk_terms(k):
        if a == b:
            mm = mean(Y[a]); val = pm(Y[a], Y[a]) - mm * mm
        else:
            val = pm(Y[a], Y[b]) - mean(Y[a]) * mean(Y[b])
        tot += coef * val
        print(f'  ({a},{b})x{coef}: {val:+.5f}', flush=True)
    print(f'\nv6(n=10) DENSE-boundary = {tot:+.6f}   [anchor +80.94197]  |diff|={abs(tot-80.94197):.3f}', flush=True)
else:
    print(f'  RESUME ({rem} left)', flush=True)
