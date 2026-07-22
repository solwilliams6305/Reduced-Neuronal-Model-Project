"""Dense-only sweep over ALL pure-bulk hard classes of Var(Y8): time each, flag max-order
and any that are slow/OOM. Tells us exactly which classes dense can absorb."""
import os, sys, time, signal, tracemalloc
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
        reps.setdefault(sg, atoms)

# pure-bulk classes with max order >=3 (the ones with no current fast path)
purebulk = [sg for sg in reps if sg[1] == 0 and sg[0] and max(sg[0]) >= 3]
purebulk.sort(key=lambda sg: (max(sg[0]) if sg[0] else 0, len(sg[0])))
class TO(Exception): pass
signal.signal(signal.SIGALRM, lambda *a: (_ for _ in ()).throw(TO()))

worst = 0.0; total_dense = 0.0; slowcnt = 0; oomcnt = 0
bymax = {}
for sg in purebulk:
    atoms = reps[sg]
    maxk = max(a[1] for a in atoms if a[0]=='U') if any(a[0]=='U' for a in atoms) else 0
    CD._MCACHE.clear()
    signal.alarm(30); t1=time.time()
    try:
        vd = CD.moment(atoms); dt = time.time()-t1; signal.alarm(0)
    except TO:
        signal.alarm(0); dt = float('inf'); slowcnt += 1
        print(f'  SLOW>30s: {sg[0]} maxk={maxk}', flush=True); continue
    except MemoryError:
        signal.alarm(0); oomcnt += 1
        print(f'  OOM: {sg[0]} maxk={maxk}', flush=True); continue
    total_dense += dt
    bymax.setdefault(maxk, [0,0.0]); bymax[maxk][0]+=1; bymax[maxk][1]+=dt
    if dt > worst:
        worst = dt
    if dt > 2.0:
        print(f'  slowish {dt:.2f}s: {sg[0]} maxk={maxk}', flush=True)
print(f'[{time.time()-t0:.1f}s] dense sweep done over {len(purebulk)} classes', flush=True)
print(f'  worst single class: {worst:.2f}s; sum over 1-rep-per-class: {total_dense:.1f}s; slow>30s={slowcnt}; oom={oomcnt}', flush=True)
for mk in sorted(bymax):
    print(f'  maxk={mk}: {bymax[mk][0]} classes, {bymax[mk][1]:.2f}s total (rep)', flush=True)
