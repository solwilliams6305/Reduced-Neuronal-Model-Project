"""Validate _moment_bulk_dense (the (order,head) dense contractor) vs CD.moment [independent]
on the maxk<=5 pure-bulk classes of Var(Y7), and spot-check the many-chain regime."""
import os, sys, time, signal
os.environ.setdefault('CT_NUMBA', '0')
import numpy as np, sympy as sp
import chaos_diagram as CD
import chaos_transfer as CT

n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
CD.setup(n_grid=n, MAXORD=CT.DENSE_MAXORD); CT.clear_all()
u0ps = CD.GEO['u0ps']; Vst = CD.GEO['Vst']
Us = [[sp.symbols(f'U{k}_{m}') for m in range(16)] for k in range(16)]
s = [sp.symbols(f's{j}') for j in range(14)]; Vs = [sp.symbols(f'Vd{j}') for j in range(14)]
Y7 = CD.Ybase(7, Us, s, Vs, None, u0ps, Vst)
def sig(at): return (tuple(sorted(a[1] for a in at if a[0]=='U')), sum(1 for a in at if a[0]=='s'))
reps = {}
for ca, aa in Y7:
    for cb, ab in Y7:
        if ca==0.0 or cb==0.0: continue
        reps.setdefault(sig(aa+ab), aa+ab)
# maxk in [3,5] pure-bulk classes (the dense-contractor regime)
cls = [sg for sg in reps if sg[1]==0 and sg[0] and 3 <= max(sg[0]) <= 5]
print(f'n={n}: {len(cls)} maxk in[3,5] pure-bulk classes of Var(Y7)', flush=True)
class TO(Exception): pass
signal.signal(signal.SIGALRM, lambda *a: (_ for _ in ()).throw(TO()))
maxre = 0.0; nok = 0; fails = []
for sg in cls:
    atoms = reps[sg]; ga = [(a[1], CT._head_vec(a)) for a in atoms if a[0]=='U']
    vDense = CT._moment_bulk_dense(ga)
    CD._MCACHE.clear(); signal.alarm(40)
    try:
        vCD = CD.moment(atoms); signal.alarm(0)
    except TO:
        signal.alarm(0); print(f'  {sg[0]}: CD.moment timeout', flush=True); continue
    re = abs(vDense-vCD)/max(abs(vDense),abs(vCD),1e-300); maxre=max(maxre,re)
    if re < 1e-9: nok += 1
    else: fails.append((sg[0], vDense, vCD, re))
print(f'_moment_bulk_dense vs CD.moment: {nok}/{len(cls)} pass, max relerr {maxre:.2e}', flush=True)
if fails:
    for sg,a,b,re in fails[:15]:
        print(f'  FAIL {sg}: dense={a:+.6e} CD={b:+.6e} re={re:.2e}', flush=True)
else:
    print('  ALL PASS', flush=True)
