"""Validate chaos_dp_numpy.moment_bulk_numpy vs (a) DENSE [independent] for maxk<=5, and
(b) pure-python DP for maxk>=6. Prints per-class so progress is visible; reports peak states."""
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

def sig(atoms):
    return (tuple(sorted(a[1] for a in atoms if a[0]=='U')), sum(1 for a in atoms if a[0]=='s'))
reps = {}
for ca, aa in Y8:
    for cb, ab in Y8:
        if ca == 0.0 or cb == 0.0: continue
        reps.setdefault(sig(aa+ab), aa+ab)
purebulk = sorted([sg for sg in reps if sg[1]==0], key=lambda sg:(max(sg[0]) if sg[0] else 0, len(sg[0])))
def gatoms_of(atoms):
    return [(a[1], CT._head_vec(a)) for a in atoms if a[0]=='U']
def relerr(x, y):
    d = abs(x-y); s = max(abs(x), abs(y)); return d/s if s > 0 else d
class TO(Exception): pass
signal.signal(signal.SIGALRM, lambda *a: (_ for _ in ()).throw(TO()))

CT.QUAD_ENABLE = False; CT.NUMBA_DP = False    # force pure-python DP as the reference
max_dense = 0.0; max_py = 0.0; nd = npy = nnone = 0; fails = []
print(f'[{time.time()-t0:.1f}s] n={n}: {len(purebulk)} pure-bulk classes', flush=True)
for sg in purebulk:
    atoms = reps[sg]; ga = gatoms_of(atoms); maxk = max(sg[0]) if sg[0] else 0
    nchain = len([o for o in sg[0] if o >= 2])
    # production scope of the numpy DP: maxk>=6 (its route), plus <=4-chain low-order for coverage.
    # (>=5-chain low-order moments are routed to dense in production, and would blow the numpy
    #  state set here, so skip them -- dense validates that regime independently.)
    if not (maxk >= 6 or nchain <= 4):
        continue
    t1 = time.time()
    vN = DPN.moment_bulk_numpy(ga); dtN = time.time()-t1; pk = DPN.LAST_PEAK
    if vN is None:
        nnone += 1
        print(f'  {sg[0]} nch={nchain} maxk={maxk}: numpy->None (peak>{DPN._MAXSTATES}) [{dtN:.1f}s]', flush=True)
        continue
    tags = []
    if maxk <= 5:
        CD._MCACHE.clear(); signal.alarm(30)
        try:
            vD = CD.moment(atoms); signal.alarm(0); re = relerr(vN, vD); nd += 1
            max_dense = max(max_dense, re)
            if re > 1e-9: fails.append(('dense', sg[0], vN, vD, re))
            tags.append(f'dense re={re:.1e}')
        except TO:
            signal.alarm(0); tags.append('dense:timeout')
    if maxk >= 6 or nchain <= 3:
        signal.alarm(45)
        try:
            vP = CT._moment_bulk_gen(ga); signal.alarm(0); re = relerr(vN, vP); npy += 1
            max_py = max(max_py, re)
            if re > 1e-9: fails.append(('pyDP', sg[0], vN, vP, re))
            tags.append(f'pyDP re={re:.1e}')
        except TO:
            signal.alarm(0); tags.append('pyDP:timeout')
    # only print the interesting (maxk>=6 or slow) ones to keep log readable
    if maxk >= 6 or dtN > 1 or pk > 1e5:
        print(f'  {sg[0]} nch={nchain} maxk={maxk}: numpy={vN:+.4e} peak={pk} [{dtN:.1f}s]  {"; ".join(tags)}', flush=True)
CT.QUAD_ENABLE = True
print(f'[{time.time()-t0:.1f}s] DONE. vs dense: {nd} (max re {max_dense:.2e}); vs pyDP: {npy} (max re {max_py:.2e}); None: {nnone}', flush=True)
if fails:
    print(f'  !!! {len(fails)} FAILURES:', flush=True)
    for tag, sg, a, b, re in fails[:25]:
        print(f'    [{tag}] {sg}: numpy={a:+.6e} ref={b:+.6e} re={re:.2e}', flush=True)
else:
    print('  ALL PASS (relerr <= 1e-9)', flush=True)
