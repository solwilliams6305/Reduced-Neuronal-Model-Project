"""Measure the pure-python transfer DP on the maxk>=6 multi-chain 'killer' moments
(the ones dense cannot build). Also instrument peak reachable-state count."""
import os, sys, time, signal
os.environ.setdefault('CT_NUMBA', '1')
import numpy as np, sympy as sp
import chaos_diagram as CD
import chaos_transfer as CT

n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
t0 = time.time()
CD.setup(n_grid=n, MAXORD=CT.DENSE_MAXORD); CT.clear_all()
u0ps = CD.GEO['u0ps']; Vst = CD.GEO['Vst']
Us = [[sp.symbols(f'U{k}_{m}') for m in range(16)] for k in range(16)]
s = [sp.symbols(f's{j}') for j in range(14)]; Vs = [sp.symbols(f'Vd{j}') for j in range(14)]
Y8 = CD.Ybase(8, Us, s, Vs, None, u0ps, Vst)

def sig(atoms):
    bulk = tuple(sorted(a[1] for a in atoms if a[0] == 'U'))
    nb = sum(1 for a in atoms if a[0] == 's')
    return (bulk, nb)
reps = {}
for ca, aa in Y8:
    for cb, ab in Y8:
        if ca == 0.0 or cb == 0.0: continue
        reps.setdefault(sig(aa+ab), aa+ab)
print(f'[{time.time()-t0:.1f}s] n={n} ready', flush=True)

# instrument peak states: monkeypatch _moment_bulk_gen to count. Easiest: wrap the DP by
# temporarily counting via a global in a patched copy. Instead, just time via moment_transfer.
killers = [(2,2,2,2,8),(2,2,2,2,2,6),(2,2,2,6),(2,2,8),(3,3,4,6),(2,2,3,3,6),
           (2,2,2,4,6),(1,2,2,2,2,7),(1,1,2,2,2,2,6),(2,2,2,3,3,4)]
class TO(Exception): pass
signal.signal(signal.SIGALRM, lambda *a: (_ for _ in ()).throw(TO()))
for tg in killers:
    sg = (tuple(sorted(tg)),0)
    if sg not in reps:
        print(f'  {tg}: no rep', flush=True); continue
    atoms = reps[sg]
    nchain = len([o for o in sg[0] if o>=2]); maxk=max(sg[0])
    CT._HMEMO.clear(); CT._TMEMO.clear(); CT._RMEMO.clear()
    signal.alarm(70); t1=time.time()
    try:
        v = CT.moment_hybrid(atoms); dt=time.time()-t1; signal.alarm(0)
        print(f'  {tg} nchain={nchain} maxk={maxk}: {v:+.4e} in {dt:.2f}s', flush=True)
    except TO:
        signal.alarm(0)
        print(f'  {tg} nchain={nchain} maxk={maxk}: >70s TIMEOUT', flush=True)
print(f'[{time.time()-t0:.1f}s] done', flush=True)
