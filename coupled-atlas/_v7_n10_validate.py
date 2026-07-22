"""v7 at n=10 with a moment-level heartbeat, to validate the engine against the anchor
v7(n=10) = -178.60840 and make progress VISIBLE (no more guessing from RSS).
Prints a heartbeat every HB moments and each (a,b) term as it completes."""
import os, sys, time
os.environ.setdefault('CT_NUMBA', '1')
import chaos_diagram as CD
import chaos_transfer as CT
from _v6_driver import vk_terms, _yfile

HB = 500
n = 10; k = 7
t0 = time.time()
CD.setup(n_grid=n, MAXORD=CT.DENSE_MAXORD); CT.clear_all()
u0ps = CD.GEO['u0ps']; Vst = CD.GEO['Vst']
Us, s, Vs, Ye = CD.load_Y(fname=_yfile(2 * k + 1))
print(f"loaded [{time.time()-t0:.0f}s]", flush=True)
Y = {idx: CD.Ybase(idx, Us, s, Vs, Ye, u0ps, Vst)
     for idx in sorted(set(p for pr in vk_terms(k) for p in pr[:2]))}
print(f"Y built [{time.time()-t0:.0f}s]", flush=True)

# heartbeat wrapper around the engine
_calls = [0]; _last = [time.time()]
_base = CT.moment_hybrid
def mom(atoms):
    _calls[0] += 1
    if _calls[0] % HB == 0:
        now = time.time()
        print(f"    .. {_calls[0]} moment-calls, {now-t0:.0f}s (+{now-_last[0]:.0f}s/{HB})", flush=True)
        _last[0] = now
    return _base(atoms)

def prod_moment(A, B):
    tot = 0.0
    for ca, aa in A:
        for cb, ab in B:
            if ca == 0.0 or cb == 0.0:
                continue
            tot += ca * cb * mom(aa + ab)
    return tot
def mean(A):
    return sum(c * mom(at) for c, at in A)

tot = 0.0
for a, b, coef in vk_terms(k):
    t1 = time.time()
    if a == b:
        mm = mean(Y[a]); val = prod_moment(Y[a], Y[a]) - mm * mm
    else:
        val = prod_moment(Y[a], Y[b]) - mean(Y[a]) * mean(Y[b])
    tot += coef * val
    print(f"  ({a},{b}) x{coef}: {val:+.5f}   [{time.time()-t0:.0f}s, term {time.time()-t1:.0f}s]", flush=True)

anchor = -178.60840
print(f"\nv7(n=10) = {tot:+.6f}   anchor {anchor:+.5f}   |diff| = {abs(tot-anchor):.2e}", flush=True)
print("MATCH" if abs(tot - anchor) < 1e-3 else "MISMATCH", flush=True)
