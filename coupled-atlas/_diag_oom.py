"""Find the OOM culprit in Var(Y7) at n=10 with the NEW routing. Print each distinct moment's
route + signature + current RSS right BEFORE computing it, flushed, so the last line before a
SIGKILL identifies the offending moment."""
import os, sys, time, resource
os.environ.setdefault('CT_NUMBA', '1')
import numpy as np, sympy as sp
import chaos_diagram as CD
import chaos_transfer as CT

n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
CD.setup(n_grid=n, MAXORD=CT.DENSE_MAXORD); CT.clear_all()
u0ps = CD.GEO['u0ps']; Vst = CD.GEO['Vst']
Us = [[sp.symbols(f'U{k}_{m}') for m in range(16)] for k in range(16)]
s = [sp.symbols(f's{j}') for j in range(14)]; Vs = [sp.symbols(f'Vd{j}') for j in range(14)]
Y7 = CD.Ybase(7, Us, s, Vs, None, u0ps, Vst)
print(f'Y7 {len(Y7)} monomials, n={n}', flush=True)

def rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024   # linux: KB -> MB

def sig(atoms):
    return (tuple(sorted(a[1] for a in atoms if a[0]=='U')), sum(1 for a in atoms if a[0]=='s'))
def route(atoms):
    if all(a[0]=='s' or a[1]==1 for a in atoms): return 'hafnian'
    if any(a[0]=='s' for a in atoms): return 'boundary'
    maxk = max((a[1] for a in atoms if a[0]=='U'), default=0)
    nchain = sum(1 for a in atoms if a[0]=='U' and a[1]>=2)
    if all(a[1]<=2 for a in atoms if a[0]=='U'): return 'quad'
    if 3<=maxk<=CT.DENSE_MAXORD and nchain>=CT.DENSE_NCHAIN: return 'DENSE'
    return 'transfer'

# enumerate distinct moments of Var(Y7)=<Y7 Y7>, compute each, track RSS
distinct = {}
for ca, aa in Y7:
    for cb, ab in Y7:
        if ca==0.0 or cb==0.0: continue
        key = tuple(sorted(aa+ab))
        if key not in distinct: distinct[key] = aa+ab
print(f'{len(distinct)} distinct moments in Var(Y7)', flush=True)
order = sorted(distinct.values(), key=lambda at: (route(at), -sig(at)[0][-1] if sig(at)[0] else 0))
t0 = time.time(); i = 0; peak = 0
import signal
class TO(Exception): pass
signal.signal(signal.SIGALRM, lambda *a: (_ for _ in ()).throw(TO()))
slow = []
for atoms in order:
    i += 1
    r = route(atoms); sg = sig(atoms)
    m0 = rss_mb()
    if m0 > peak: peak = m0
    if i % 200 == 0:
        print(f'  [{i}/{len(distinct)}] {r:9s} {sg}  RSS={m0:.0f}MB t={time.time()-t0:.0f}s', flush=True)
    if r in ('transfer', 'DENSE'):
        print(f'   ->[{i}] {r} {sg}', flush=True)   # print BEFORE compute: last line = culprit
    t1 = time.time()
    signal.alarm(15)
    try:
        CT.moment_hybrid(atoms)
        signal.alarm(0)
        dt = time.time() - t1
        if dt > 2:
            print(f'   SLOW {dt:.1f}s: {r} {sg}  RSS={rss_mb():.0f}MB', flush=True)
            slow.append((dt, r, sg))
    except TO:
        signal.alarm(0)
        print(f'   >>>> HANG(>15s): {r} {sg}  RSS={rss_mb():.0f}MB  atoms={sorted(atoms)}', flush=True)
        slow.append((999, r, sg))
print(f'DONE {i} moments, peak RSS {peak:.0f}MB, {time.time()-t0:.0f}s, slow={len(slow)}', flush=True)
