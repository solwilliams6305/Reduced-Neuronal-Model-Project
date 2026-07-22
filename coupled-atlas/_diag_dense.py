"""Compare dense engine vs transfer on representative hard classes, to see which route
is viable per class. Dense builds tensors only up to the max single-atom order, so
low-max-order many-chain moments may be cheap in dense even when transfer explodes."""
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
print(f'[{time.time()-t0:.1f}s] Y8 {len(Y8)} monomials, n={n}', flush=True)

def sig(atoms):
    bulk = tuple(sorted(a[1] for a in atoms if a[0] == 'U'))
    nb = sum(1 for a in atoms if a[0] == 's')
    return (bulk, nb)

reps = {}
for ca, aa in Y8:
    for cb, ab in Y8:
        if ca == 0.0 or cb == 0.0: continue
        atoms = aa + ab
        sg = sig(atoms)
        if sg[1] != 0: continue
        reps.setdefault(sg, atoms)

targets = [(2,2,2,2,2,3,3),(2,2,2,2,2,2,4),(2,2,2,2,4,4),(2,2,2,2,3,5),
           (2,2,3,3,3,3),(2,2,2,2,8),(2,2,2,6),(2,2,8),(3,3,4,6),(2,2,2,2,2,6)]
class TO(Exception): pass
signal.signal(signal.SIGALRM, lambda *a: (_ for _ in ()).throw(TO()))

for tg in targets:
    sg = (tuple(sorted(tg)), 0)
    if sg not in reps:
        print(f'  {tg}: (no rep)', flush=True); continue
    atoms = reps[sg]
    maxk = max(a[1] for a in atoms if a[0]=='U')
    # DENSE
    CD._MCACHE.clear()
    signal.alarm(40); t1=time.time()
    try:
        vd = CD.moment(atoms); dtd = time.time()-t1; signal.alarm(0)
        ds = f'{vd:+.4e} in {dtd:.2f}s'
    except TO:
        signal.alarm(0); ds = '>40s TIMEOUT'
    except MemoryError:
        signal.alarm(0); ds = 'OOM'
    # TRANSFER
    CT._HMEMO.clear(); CT._TMEMO.clear(); CT._RMEMO.clear()
    signal.alarm(40); t1=time.time()
    try:
        vt = CT.moment_hybrid(atoms); dtt=time.time()-t1; signal.alarm(0)
        ts = f'{vt:+.4e} in {dtt:.2f}s'
    except TO:
        signal.alarm(0); ts = '>40s TIMEOUT'
    print(f'  {tg} maxk={maxk}:  dense[{ds}]   transfer[{ts}]', flush=True)
print(f'[{time.time()-t0:.1f}s] done', flush=True)
