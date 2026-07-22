"""End-to-end correctness gates with the FIXED engine (numpy DP + dense router).
Computes v3, v4, v6 and checks against the documented values. Bypasses the slow 34MB
Y-expr parse by using the cached _ybase_*.pkl (all present)."""
import os, sys, time
os.environ.setdefault('CT_NUMBA', '1')      # numba absent -> auto numpy DP
import numpy as np, sympy as sp
import chaos_diagram as CD
import chaos_transfer as CT
from _v6_driver import vk_terms

def compute_vk(k, n, verbose=True):
    setup_ord = min(k + 1, CT.DENSE_MAXORD)
    CD.setup(n_grid=n, MAXORD=setup_ord); CT.clear_all()
    u0ps = CD.GEO['u0ps']; Vst = CD.GEO['Vst']
    Us = [[sp.symbols(f'U{a}_{m}') for m in range(16)] for a in range(16)]
    s = [sp.symbols(f's{j}') for j in range(14)]; Vs = [sp.symbols(f'Vd{j}') for j in range(14)]
    Y = {idx: CD.Ybase(idx, Us, s, Vs, None, u0ps, Vst)
         for idx in sorted(set(p for pr in vk_terms(k) for p in pr[:2]))}
    mom = CT.moment_hybrid
    def pm(A, B):
        t = 0.0
        for ca, aa in A:
            for cb, ab in B:
                if ca == 0.0 or cb == 0.0: continue
                t += ca * cb * mom(aa + ab)
        return t
    def mean(A): return sum(c * mom(at) for c, at in A)
    tot = 0.0; t0 = time.time()
    for a, b, coef in vk_terms(k):
        if a == b:
            mm = mean(Y[a]); val = pm(Y[a], Y[a]) - mm * mm
        else:
            val = pm(Y[a], Y[b]) - mean(Y[a]) * mean(Y[b])
        tot += coef * val
        if verbose:
            print(f'    ({a},{b})x{coef}: {val:+.5f}  [{time.time()-t0:.0f}s]', flush=True)
    return tot

if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    t0 = time.time()
    jobs = {'v3': (3, 20, -0.060), 'v4': (4, 24, -0.454), 'v6': (6, 10, +80.94197)}
    todo = [which] if which in jobs else ['v3', 'v4', 'v6']
    for name in todo:
        k, n, expect = jobs[name]
        v = compute_vk(k, n, verbose=True)
        tag = 'OK' if abs(v - expect) < max(1e-3, 0.02 * abs(expect)) else 'CHECK'
        print(f'[{time.time()-t0:.0f}s] {name}(n={n}) = {v:+.6f}   expect ~{expect:+.5f}   {tag}', flush=True)
