"""Compute low v_k (transfer convention) across grids and save to a results pickle for
Richardson extrapolation.  Usage: python3 _ladder_low.py <k> <n1,n2,...>"""
import os, sys, time, pickle
os.environ.setdefault('CT_NUMBA', '0')
import _val_ladder as VL

k = int(sys.argv[1])
grids = [int(x) for x in sys.argv[2].split(',')]
RES = '_ladder_transfer_results.pkl'
def _load():
    try:
        return pickle.load(open(RES, 'rb'))
    except Exception:
        return {}
def _save(d):
    tmp = RES + '.tmp'
    with open(tmp, 'wb') as f:
        pickle.dump(d, f)
    os.replace(tmp, RES)            # atomic -> never leaves a truncated file
res = _load()
for n in grids:
    if (k, n) in res:
        print(f'  v{k}(n={n}) = {res[(k,n)]:+.6f}  [cached]', flush=True); continue
    t0 = time.time()
    v = VL.compute_vk(k, n, verbose=False)
    res = _load(); res[(k, n)] = v; _save(res)
    print(f'  v{k}(n={n}) = {v:+.6f}   [{time.time()-t0:.0f}s]', flush=True)
