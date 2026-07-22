"""Measure dense-contractor vs numpy-DP on every maxk in [3,5] pure-bulk class of Var(Y7),
to set the router. Reports per class: nchain, nlegs, dense time, numpy time (capped), winner,
and cross-checks they AGREE."""
import os, sys, time, signal
os.environ.setdefault('CT_NUMBA', '0')
import numpy as np, sympy as sp
import chaos_diagram as CD
import chaos_transfer as CT
import chaos_dp_numpy as DPN

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
cls = sorted([sg for sg in reps if sg[1]==0 and sg[0] and 3<=max(sg[0])<=5],
             key=lambda sg: (len([o for o in sg[0] if o>=2]), sg[0].count(1)))
print(f'n={n}: {len(cls)} maxk in[3,5] pure-bulk classes', flush=True)
class TO(Exception): pass
signal.signal(signal.SIGALRM, lambda *a:(_ for _ in ()).throw(TO()))
def timed(fn, ga, cap):
    signal.alarm(cap)
    try:
        t=time.time(); v=fn(ga); signal.alarm(0); return v, time.time()-t
    except TO:
        signal.alarm(0); return None, float('inf')
for sg in cls:
    atoms=reps[sg]; ga=[(a[1],CT._head_vec(a)) for a in atoms if a[0]=='U']
    nch=len([o for o in sg[0] if o>=2]); nl=sg[0].count(1)
    vD,dD = timed(CT._moment_bulk_dense, ga, 25)
    vN,dN = timed(DPN.moment_bulk_numpy, ga, 25)
    re = abs(vD-vN)/max(abs(vD),abs(vN),1e-300) if (vD is not None and vN is not None) else -1
    win = 'DENSE' if dD < dN else 'numpy'
    flag = '' if re<1e-9 or re<0 else '  MISMATCH!'
    print(f'  {sg[0]} nch={nch} nl={nl}: dense={dD:6.2f}s numpy={dN:6.2f}s -> {win}  re={re:.0e}{flag}', flush=True)
print('done', flush=True)
