"""Checkpointed v_k computation: the workspace kills background jobs when the launching bash
call ends, so we compute within <=~38s calls and persist the moment cache to disk between calls.

Usage:  python3 _ckpt.py <k> <n> [budget_s]
  - enumerates the distinct moments needed for v_k (cached to disk once),
  - computes uncached ones until the time budget, saving the cache,
  - when all are cached, assembles v_k and prints it + the anchor check.
Re-run repeatedly until it prints 'v{k}(n=..) ='.  State: _ckpt_{k}_{n}_*.pkl
"""
import os, sys, time, pickle
os.environ.setdefault('CT_NUMBA', '1')      # numba absent -> numpy DP auto
import numpy as np, sympy as sp
import chaos_diagram as CD
import chaos_transfer as CT
from _v6_driver import vk_terms

k = int(sys.argv[1]); n = int(sys.argv[2])
budget = float(sys.argv[3]) if len(sys.argv) > 3 else 38.0
t0 = time.time()
tag = ''
if os.environ.get('NO_DENSE'):
    CT.DENSE_ROUTE = False; CT.DENSE_NCHAIN = 999   # disable dense routing everywhere
    tag = '_nodense'
if os.environ.get('NO_QUAD'):
    CT.QUAD_ENABLE = False; tag += '_noquad'
if os.environ.get('BND'):
    CT.BND_ENGINE = os.environ['BND']; tag += '_bnd' + os.environ['BND']
MOMF = f'_ckpt_{k}_{n}_moments.pkl'
CACHEF = f'_ckpt_{k}_{n}{tag}_cache.pkl'

CD.setup(n_grid=n, MAXORD=CT.DENSE_MAXORD); CT.clear_all()
u0ps = CD.GEO['u0ps']; Vst = CD.GEO['Vst']
Us = [[sp.symbols(f'U{a}_{m}') for m in range(16)] for a in range(16)]
s = [sp.symbols(f's{j}') for j in range(14)]; Vs = [sp.symbols(f'Vd{j}') for j in range(14)]
Y = {idx: CD.Ybase(idx, Us, s, Vs, None, u0ps, Vst)
     for idx in sorted(set(p for pr in vk_terms(k) for p in pr[:2]))}

# ---- enumerate distinct moments (cache to disk once) ----
if os.path.exists(MOMF):
    with open(MOMF, 'rb') as f:
        needed = pickle.load(f)
else:
    needed = set()
    for a, b, _ in vk_terms(k):
        for c, at in Y[a]:
            if c != 0.0:
                needed.add(tuple(sorted(at)))            # <Y_a> mean moment
        if b != a:
            for c, at in Y[b]:
                if c != 0.0:
                    needed.add(tuple(sorted(at)))
        for ca, aa in Y[a]:
            if ca == 0.0: continue
            for cb, ab in Y[b]:
                if cb == 0.0: continue
                needed.add(tuple(sorted(aa + ab)))
    needed = list(needed)
    with open(MOMF, 'wb') as f:
        pickle.dump(needed, f)
print(f'[{time.time()-t0:.0f}s] {len(needed)} distinct moments needed for v{k}(n={n})', flush=True)

# ---- load cache ----
cache = {}
if os.path.exists(CACHEF):
    with open(CACHEF, 'rb') as f:
        cache = pickle.load(f)
CT._HMEMO.update(cache)
print(f'[{time.time()-t0:.0f}s] cache loaded: {len(cache)}/{len(needed)} done', flush=True)

# ---- compute uncached, time-budgeted ----
todo = [m for m in needed if m not in CT._HMEMO]
done_this = 0
_trace = os.environ.get('TRACE')
SAVE_EVERY = 2000            # periodic checkpoint + heartbeat, so a long uninterrupted run is
                            # crash-safe and visibly progressing (matters for the n=12/14 grids)
def _flush_cache():
    for mm in needed:
        if mm in CT._HMEMO:
            cache[mm] = CT._HMEMO[mm]
    with open(CACHEF, 'wb') as f:
        pickle.dump(cache, f)
for i, m in enumerate(todo):
    t1 = time.time()
    CT.moment_hybrid(list(m))
    dt = time.time() - t1
    if _trace and dt > 1.0:
        bulk = tuple(sorted(a[1] for a in m if a[0]=='U')); nb = sum(1 for a in m if a[0]=='s')
        print(f'    SLOW {dt:.1f}s bulk={bulk} nb={nb}', flush=True)
    done_this += 1
    if done_this % SAVE_EVERY == 0:
        _flush_cache()
        ncached = sum(1 for mm in needed if mm in cache)
        print(f'    .. heartbeat: {ncached}/{len(needed)} cached, {time.time()-t0:.0f}s elapsed', flush=True)
    if (time.time() - t0) > budget:
        break
_flush_cache()
remaining = sum(1 for m in needed if m not in cache)
print(f'[{time.time()-t0:.0f}s] computed {done_this} this run; {len(cache)}/{len(needed)} cached; {remaining} remaining', flush=True)

# ---- assemble if complete ----
if remaining == 0:
    def mom(at): return cache[tuple(sorted(at))]
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
    anchors = {6: (10, 80.94197), 7: (10, -178.60840)}
    print(f'\nv{k}(n={n}) = {tot:+.6f}', flush=True)
    if k in anchors and anchors[k][0] == n:
        ac = anchors[k][1]
        print(f'  anchor {ac:+.5f}  |diff|={abs(tot-ac):.2e}  {"MATCH" if abs(tot-ac)<1e-3 else "MISMATCH"}', flush=True)
else:
    print(f'  RESUME: re-run to continue ({remaining} left)', flush=True)
