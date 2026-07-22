"""Time representative hard moments through the CURRENT engine (pure-python DP, no numba).
Uses SIGALRM to cap each moment so one killer doesn't stall the probe."""
import os, sys, time, signal
os.environ.setdefault('CT_NUMBA', '1')
import numpy as np, sympy as sp
import chaos_diagram as CD
import chaos_transfer as CT

n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
CAP = int(sys.argv[2]) if len(sys.argv) > 2 else 30
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

# collect one representative atom-list per pure-bulk class (nb==0)
reps = {}
for ca, aa in Y8:
    for cb, ab in Y8:
        if ca == 0.0 or cb == 0.0:
            continue
        atoms = aa + ab
        sg = sig(atoms)
        if sg[1] != 0:
            continue
        if sg not in reps:
            reps[sg] = atoms

def is_hard(bulk):
    chains = [o for o in bulk if o >= 2]
    return (not all(o <= 2 for o in bulk)) and len(chains) >= 3

hard = sorted([sg for sg in reps if is_hard(sg[0])],
              key=lambda sg: (-len([o for o in sg[0] if o >= 2]), -sum(sg[0])))

class TO(Exception):
    pass
def _h(sig_, frm):
    raise TO()
signal.signal(signal.SIGALRM, _h)

print(f'timing {len(hard)} hard pure-bulk classes (cap {CAP}s each):', flush=True)
slow = []
for sg in hard:
    atoms = reps[sg]
    CT._HMEMO.clear(); CT._TMEMO.clear(); CT._RMEMO.clear()
    nchain = len([o for o in sg[0] if o >= 2])
    signal.alarm(CAP)
    t1 = time.time()
    try:
        v = CT.moment_hybrid(atoms)
        dt = time.time() - t1
        signal.alarm(0)
        tag = '' if dt < 1 else ('  <<< SLOW' if dt < CAP else '')
        if dt > 1:
            print(f'  {sg[0]} nchain={nchain}: {dt:.2f}s{tag}', flush=True)
        if dt > 1:
            slow.append((sg[0], nchain, dt))
    except TO:
        signal.alarm(0)
        print(f'  {sg[0]} nchain={nchain}: >{CAP}s  TIMEOUT <<<<', flush=True)
        slow.append((sg[0], nchain, float('inf')))
print(f'[{time.time()-t0:.1f}s] done. slow(>1s) classes: {len(slow)}', flush=True)
