"""Compute v_k with the FIXED engine, a heartbeat, and a slow-moment logger (prints any single
moment that takes > SLOW s, with its atom signature). Bypasses the 34MB parse via cached ybase."""
import os, sys, time
os.environ.setdefault('CT_NUMBA', '1')       # numba absent -> numpy DP auto
import numpy as np, sympy as sp
import chaos_diagram as CD
import chaos_transfer as CT
from _v6_driver import vk_terms

k = int(sys.argv[1]) if len(sys.argv) > 1 else 6
n = int(sys.argv[2]) if len(sys.argv) > 2 else 10
SLOW = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0
HB = 400
t0 = time.time()
CD.setup(n_grid=n, MAXORD=CT.DENSE_MAXORD); CT.clear_all()
u0ps = CD.GEO['u0ps']; Vst = CD.GEO['Vst']
Us = [[sp.symbols(f'U{a}_{m}') for m in range(16)] for a in range(16)]
s = [sp.symbols(f's{j}') for j in range(14)]; Vs = [sp.symbols(f'Vd{j}') for j in range(14)]
Y = {idx: CD.Ybase(idx, Us, s, Vs, None, u0ps, Vst)
     for idx in sorted(set(p for pr in vk_terms(k) for p in pr[:2]))}
print(f'[{time.time()-t0:.0f}s] Y built', flush=True)

# trace any slow bulk-engine call (catches whichever route is slow, regardless of atom count)
if os.environ.get('TRACE_BULK'):
    _obg = CT._moment_bulk_gen
    def _traced(ga):
        import time as _t
        t1 = _t.time(); v = _obg(ga); dt = _t.time() - t1
        if dt > 2:
            orders = sorted(o for o, _ in ga); nch = sum(1 for o, _ in ga if o >= 2)
            nlt = len(set(h.tobytes() for o, h in ga if o == 1)); mk = max((o for o, _ in ga), default=0)
            print(f'    BULK_SLOW {dt:.1f}s orders={orders} nchain={nch} maxk={mk} nlegtypes={nlt}', flush=True)
        return v
    CT._moment_bulk_gen = _traced

_calls = [0]; _base = CT.moment_hybrid
def mom(atoms):
    _calls[0] += 1
    if _calls[0] % HB == 0:
        print(f'    .. {_calls[0]} calls, {time.time()-t0:.0f}s', flush=True)
    # pre-log EVERY uncached moment so the last line before a stall = the exact culprit
    key = tuple(sorted(atoms))
    if os.environ.get('TRACE_ALL') and key not in CT._HMEMO:
        bulk = tuple(sorted(a[1] for a in atoms if a[0]=='U')); nb = sum(1 for a in atoms if a[0]=='s')
        print(f'      ->#{_calls[0]} bulk={bulk} nb={nb}', flush=True)
    t1 = time.time(); v = _base(atoms); dt = time.time() - t1
    if dt > SLOW:
        bulk = tuple(sorted(a[1] for a in atoms if a[0]=='U')); nb = sum(1 for a in atoms if a[0]=='s')
        print(f'    SLOW {dt:.1f}s: bulk={bulk} nb={nb}', flush=True)
    return v
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
    t1 = time.time()
    if a == b:
        mm = mean(Y[a]); val = pm(Y[a], Y[a]) - mm * mm
    else:
        val = pm(Y[a], Y[b]) - mean(Y[a]) * mean(Y[b])
    tot += coef * val
    print(f'  ({a},{b})x{coef}: {val:+.5f}  [{time.time()-t0:.0f}s, term {time.time()-t1:.0f}s]', flush=True)
anchors = {6: (10, 80.94197), 7: (10, -178.60840)}
print(f'\n[{time.time()-t0:.0f}s] v{k}(n={n}) = {tot:+.6f}', flush=True)
if k in anchors and anchors[k][0] == n:
    ac = anchors[k][1]
    print(f'  anchor {ac:+.5f}  |diff|={abs(tot-ac):.2e}  {"MATCH" if abs(tot-ac)<1e-3 else "MISMATCH"}', flush=True)
