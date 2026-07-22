"""Focused validation of chaos_dp_numpy over the numpy DP's PRODUCTION scope: every pure-bulk
maxk>=6 class of Var(Y8) (an order>=6 chain forces <=6 chains + heavy pruning). Reference =
pure-python DP. Prints one line per class, flushed, so partial progress survives a crash."""
import os, sys, time, signal
os.environ.setdefault('CT_NUMBA', '0')
import numpy as np, sympy as sp
import chaos_diagram as CD
import chaos_transfer as CT
import chaos_dp_numpy as DPN

n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
t0 = time.time()
CD.setup(n_grid=n, MAXORD=CT.DENSE_MAXORD); CT.clear_all()
u0ps = CD.GEO['u0ps']; Vst = CD.GEO['Vst']
Us = [[sp.symbols(f'U{k}_{m}') for m in range(16)] for k in range(16)]
s = [sp.symbols(f's{j}') for j in range(14)]; Vs = [sp.symbols(f'Vd{j}') for j in range(14)]
Y8 = CD.Ybase(8, Us, s, Vs, None, u0ps, Vst)
def sig(at): return (tuple(sorted(a[1] for a in at if a[0]=='U')), sum(1 for a in at if a[0]=='s'))
reps = {}
for ca, aa in Y8:
    for cb, ab in Y8:
        if ca == 0.0 or cb == 0.0: continue
        reps.setdefault(sig(aa+ab), aa+ab)
hi = sorted([sg for sg in reps if sg[1]==0 and sg[0] and max(sg[0])>=6],
            key=lambda sg:(len([o for o in sg[0] if o>=2]), max(sg[0])))
print(f'[{time.time()-t0:.1f}s] n={n}: {len(hi)} maxk>=6 pure-bulk classes', flush=True)
CT.QUAD_ENABLE = False; CT.NUMBA_DP = False
maxre = 0.0; npass = 0; fails = []
for sg in hi:
    ga = [(a[1], CT._head_vec(a)) for a in reps[sg] if a[0]=='U']
    nch = len([o for o in sg[0] if o>=2])
    t1 = time.time(); vN = DPN.moment_bulk_numpy(ga); dN = time.time()-t1; pk = DPN.LAST_PEAK
    t1 = time.time(); vP = CT._moment_bulk_gen(ga); dP = time.time()-t1
    if vN is None:
        print(f'  {sg[0]} nch={nch}: numpy None; pyDP={vP:+.4e} in {dP:.1f}s', flush=True); continue
    re = abs(vN-vP)/max(abs(vN),abs(vP),1e-300); maxre = max(maxre, re)
    ok = 'OK' if re < 1e-9 else 'FAIL'
    if re >= 1e-9: fails.append((sg[0], vN, vP, re))
    else: npass += 1
    print(f'  {sg[0]} nch={nch}: numpy={vN:+.5e}[{dN:.2f}s,pk{pk}] pyDP[{dP:.2f}s] re={re:.1e} {ok}', flush=True)
print(f'[{time.time()-t0:.1f}s] DONE {npass} pass, max relerr {maxre:.2e}', flush=True)
print('ALL PASS' if not fails else f'{len(fails)} FAILURES', flush=True)
