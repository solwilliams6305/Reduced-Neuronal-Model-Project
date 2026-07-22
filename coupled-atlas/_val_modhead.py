"""Validate _moment_bulk_dense vs numpy DP vs pure-python DP on MODIFIED-HEAD (boundary-reduced)
chains -- the case CD.moment cannot represent, so it was never cross-checked. Builds reduced
moments the way _mt_recurse does (reduce_head), then compares all three bulk engines."""
import os, sys, time
os.environ.setdefault('CT_NUMBA', '0')
import numpy as np, sympy as sp
import chaos_diagram as CD
import chaos_transfer as CT
import chaos_dp_numpy as DPN

n = 10
CD.setup(n_grid=n, MAXORD=CT.DENSE_MAXORD); CT.clear_all()
u0ps = CD.GEO['u0ps']; Vst = CD.GEO['Vst']

# build some modified heads via reduce_head/reduce_scalar (as _mt_recurse does)
def head_std(k, m): return CT._head_vec(('U', k, m))
mods = []
for (k, m) in [(3,0),(4,0),(2,1),(5,0),(3,1)]:
    for j in [0,1,2]:
        mods.append(CT._reduce_head(head_std(k, m), j))   # node-pinned head over next slot

# construct test moments mixing standard chains, modified-head chains, and legs
rng = np.random.default_rng(0)
def relerr(x,y):
    d=abs(x-y); s=max(abs(x),abs(y),1e-300); return d/s
CT.QUAD_ENABLE=False; CT.NUMBA_DP=False
tests = []
# (a) modified-head chains only, various nchain
tests.append([(2, mods[0]),(2, mods[1]),(3, mods[3]),(2, mods[5])])           # 4 chains
tests.append([(2, mods[0]),(2, mods[1]),(2, mods[2]),(2, mods[6]),(2, mods[8])]) # 5 chains
tests.append([(3, mods[3]),(2, mods[0]),(2, mods[7]),(1, head_std(1,0)),(1, head_std(1,0))]) # legs
tests.append([(2, mods[0]),(2, mods[1]),(2, mods[4]),(2, mods[7]),(1, head_std(1,1)),(1, head_std(1,1))])
tests.append([(4, mods[9]),(2, mods[0]),(2, mods[1]),(1, head_std(1,0))])
tests.append([(2, head_std(2,0)),(2, mods[0]),(3, mods[3]),(2, mods[6])])       # mix std + mod
maxre = 0.0
for i, ga in enumerate(tests):
    orders = sorted(o for o,_ in ga)
    if sum(orders) % 2 != 0:
        ga = ga + [(1, head_std(1,0))]; orders = sorted(o for o,_ in ga)
    vD = CT._moment_bulk_dense(ga)
    vN = DPN.moment_bulk_numpy(ga)
    vP = CT._moment_bulk_gen(ga) if False else None   # pure DP via forcing below
    # force pure-python DP: disable numpy+dense route
    _sr, _dr = CT.NUMPY_DP, CT.DENSE_ROUTE
    CT.NUMPY_DP=False; CT.DENSE_ROUTE=False
    vP = CT._moment_bulk_gen(ga)
    CT.NUMPY_DP, CT.DENSE_ROUTE = _sr, _dr
    reDN = relerr(vD, vN); reDP = relerr(vD, vP); rePN = relerr(vP, vN)
    maxre = max(maxre, reDN, reDP, rePN)
    tag = 'OK' if max(reDN,reDP,rePN) < 1e-9 else 'FAIL <<<'
    print(f'  test{i} orders={orders}: dense={vD:+.5e} numpy={vN:+.5e} pyDP={vP:+.5e}  reDN={reDN:.1e} reDP={reDP:.1e} {tag}', flush=True)
print(f'max relerr {maxre:.2e}  {"ALL PASS" if maxre<1e-9 else "MISMATCH FOUND"}', flush=True)
