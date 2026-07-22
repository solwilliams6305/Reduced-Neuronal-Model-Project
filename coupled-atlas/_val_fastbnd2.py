"""Validate the integrated fast dense-boundary path (moment_hybrid BND='dense') vs CD.moment
(pure dense) on actual j>=2 boundary moments of v6, and time both."""
import os, sys, time, pickle
os.environ.setdefault('CT_NUMBA', '0')
import numpy as np
import chaos_diagram as CD
import chaos_transfer as CT

n = 10
CD.setup(n_grid=n, MAXORD=CT.DENSE_MAXORD); CT.clear_all()
needed = pickle.load(open('_ckpt_6_10_moments.pkl', 'rb'))
def is_j2(m):
    js = [a[1] for a in m if a[0] == 's']
    if not js or max(js) < 2: return False
    return max((a[1] for a in m if a[0] == 'U'), default=0) <= CT.DENSE_MAXORD
j2 = [m for m in needed if is_j2(m)]
import random; random.seed(3); random.shuffle(j2)
print(f'{len(j2)} j>=2 boundary moments; sampling', flush=True)
import signal
class TO(Exception): pass
signal.signal(signal.SIGALRM, lambda *a: (_ for _ in ()).throw(TO()))
CT.BND_ENGINE = 'dense'
maxre = 0.0; t_fast = 0.0; nchk = 0; nfall = 0; nord2 = 0
t0 = time.time()
for m in j2:
    if time.time() - t0 > 38 or nchk >= 40: break
    atoms = list(m)
    has3 = any(a[0]=='U' and a[1]>=3 for a in atoms)
    if not has3: nord2 += 1     # all-order<=2: the fast-path class
    CT._HMEMO.clear(); CT._TMEMO.clear(); CT._RMEMO.clear()
    t1 = time.time(); vfast = CT.moment_hybrid(atoms); dtf = time.time()-t1; t_fast += dtf
    CD._MCACHE.clear()
    signal.alarm(9)
    try:
        vcd = CD.moment(atoms); signal.alarm(0)
    except TO:
        signal.alarm(0); continue     # reference too slow; skip (fast path already computed)
    re = abs(vfast-vcd)/max(abs(vfast),abs(vcd),1e-300); maxre = max(maxre, re); nchk += 1
    if has3: nfall += 1
    if re > 1e-9:
        bulk=tuple(sorted(a[1] for a in atoms if a[0]=='U')); js=sorted(a[1] for a in atoms if a[0]=='s')
        print(f'  DISAGREE bulk={bulk} s={js}: fast={vfast:+.6e} cd={vcd:+.6e} re={re:.1e} [ford2={not has3}]', flush=True)
print(f'checked {nchk} vs CD.moment ({nord2} all-order<=2, {nfall} order>=3 fallbacks); max relerr {maxre:.2e}', flush=True)
print(f'  fast-path total {t_fast:.2f}s over {nchk} moments (~{t_fast/max(nchk,1)*1000:.0f} ms/moment)', flush=True)
print('  ' + ('PASS' if maxre < 1e-9 else 'FAIL'), flush=True)
