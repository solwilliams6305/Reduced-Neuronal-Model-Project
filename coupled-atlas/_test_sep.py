"""Validate _moment_sep (rank-2 separable DP) against the generic transfer DP, incl. the residual
v6 slow class (>=3 chains, one of order >=3) and boundary moments; then time the speedup."""
import time
import chaos_diagram as CD
import chaos_transfer as CT

def both(atoms, n):
    CD.setup(n_grid=n, MAXORD=2)
    CT.SEP_ENABLE = True; CT.QUAD_ENABLE = True; CT.clear_all()
    t0 = time.time(); vs = CT.moment_transfer(atoms); ts = time.time() - t0
    CT.SEP_ENABLE = False; CT.QUAD_ENABLE = False; CT.clear_all()
    t0 = time.time(); vd = CT.moment_transfer(atoms); td = time.time() - t0
    CT.SEP_ENABLE = True; CT.QUAD_ENABLE = True
    rel = abs(vs - vd) / max(abs(vd), 1e-300)
    return vs, vd, rel, ts, td

cases = [
    [('U', 3, 0), ('U', 2, 0), ('U', 1, 0)],
    [('U', 3, 0), ('U', 3, 1)],                                   # vs _moment_two? (2 chains -> two-path; force sep? routed same)
    [('U', 3, 0), ('U', 2, 1), ('U', 2, 0), ('U', 1, 1)],
    [('U', 3, 0), ('U', 3, 1), ('U', 2, 0)],
    [('U', 4, 0), ('U', 3, 1), ('U', 3, 0)],
    [('U', 5, 0), ('U', 2, 1), ('U', 1, 0)],
    [('U', 3, 0), ('U', 2, 0), ('U', 2, 1), ('U', 1, 0), ('U', 1, 1), ('U', 1, 0), ('U', 1, 1)],
    [('U', 4, 1), ('U', 2, 0), ('U', 2, 0), ('U', 1, 0), ('U', 1, 0)],
    [('U', 3, 0), ('U', 3, 0), ('U', 2, 0), ('U', 1, 0)] ,
    [('U', 7, 0), ('U', 3, 1), ('U', 2, 0), ('U', 1, 0), ('U', 1, 1)],
    # boundary moments (bottom out through sep after reduce_head)
    [('U', 3, 0), ('U', 2, 0), ('U', 2, 1), ('s', 0), ('U', 1, 0), ('U', 1, 1)],
    [('U', 4, 0), ('U', 3, 1), ('U', 2, 0), ('s', 1), ('s', 0), ('U', 1, 0)],
]

maxrel = 0.0
for atoms in cases:
    vs, vd, rel, ts, td = both(atoms, 8)
    maxrel = max(maxrel, rel)
    print(f"{'OK ' if rel < 1e-9 else 'FAIL'} {atoms}: sep={vs:+.10e} dp={vd:+.10e} rel={rel:.2e} ({ts*1e3:.0f}ms vs {td*1e3:.0f}ms)", flush=True)
print(f"max rel deviation (n=8): {maxrel:.2e}")

# precision + timing at production grids on the residual slow class
slow = [('U', 3, 0), ('U', 3, 1), ('U', 2, 0), ('U', 1, 0), ('U', 1, 1), ('U', 1, 0), ('U', 1, 1)]
for n in [12, 16, 20]:
    vs, vd, rel, ts, td = both(slow, n)
    print(f"n={n}: slow-class sep {ts*1e3:.0f}ms vs DP {td*1e3:.0f}ms ({td/max(ts,1e-9):.0f}x) rel={rel:.1e}", flush=True)
