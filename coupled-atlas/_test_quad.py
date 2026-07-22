"""Validate _moment_quad against the transfer DP (machine precision) on all-order<=2 bulk
moments, incl. the v6 bottleneck class (3 chains + 8 legs), then time the speedup."""
import time
import numpy as np
import chaos_diagram as CD
import chaos_transfer as CT

CD.setup(n_grid=8, MAXORD=2)

cases = [
    [('U', 2, 0), ('U', 2, 0)],
    [('U', 2, 0), ('U', 2, 1)],
    [('U', 2, 1), ('U', 1, 0), ('U', 1, 1)],
    [('U', 2, 0), ('U', 2, 1), ('U', 2, 0), ('U', 2, 1)],
    [('U', 2, 0)] * 3 + [('U', 1, 0)] * 2,
    [('U', 2, 1)] * 2 + [('U', 1, 0)] * 4,
    [('U', 2, 0), ('U', 2, 1)] + [('U', 1, 0), ('U', 1, 1)] * 2,
    [('U', 1, 0)] * 6,
    [('U', 2, 0)] * 4,
    [('U', 2, 0)] * 5 + [('U', 1, 1)] * 2,
    # the v6 bottleneck class: 3 order-2 chains + 8 order-1 legs (order-14 moment)
    [('U', 2, 0)] * 3 + [('U', 1, 0)] * 4 + [('U', 1, 1)] * 4,
    [('U', 2, 1)] * 3 + [('U', 1, 0)] * 8,
]

maxrel = 0.0
for atoms in cases:
    CT.QUAD_ENABLE = True; CT.clear_all()
    t0 = time.time(); vq = CT.moment_transfer(atoms); tq = time.time() - t0
    CT.QUAD_ENABLE = False; CT.clear_all()
    t0 = time.time(); vd = CT.moment_transfer(atoms); td = time.time() - t0
    rel = abs(vq - vd) / max(abs(vd), 1e-300)
    maxrel = max(maxrel, rel)
    tag = "OK " if rel < 1e-10 else "FAIL"
    print(f"{tag} {atoms}: quad={vq:+.12e} dp={vd:+.12e} rel={rel:.2e}  ({tq*1e3:.1f}ms vs {td*1e3:.1f}ms)", flush=True)
CT.QUAD_ENABLE = True

# boundary moment exercising _mt_recurse -> quad bottom-out
for atoms in [[('U', 2, 0)] * 2 + [('U', 1, 0)] * 3 + [('s', 0)],
              [('U', 2, 1)] * 3 + [('U', 1, 0)] * 2 + [('s', 1), ('s', 0)]]:
    CT.QUAD_ENABLE = True; CT.clear_all(); vq = CT.moment_transfer(atoms)
    CT.QUAD_ENABLE = False; CT.clear_all(); vd = CT.moment_transfer(atoms)
    rel = abs(vq - vd) / max(abs(vd), 1e-300)
    maxrel = max(maxrel, rel)
    print(f"{'OK ' if rel < 1e-10 else 'FAIL'} bnd {atoms}: rel={rel:.2e}", flush=True)
CT.QUAD_ENABLE = True

print(f"max rel deviation: {maxrel:.2e}")

# timing at production grid for the bottleneck class
for n in [14, 16, 18]:
    CD.setup(n_grid=n, MAXORD=2)
    atoms = [('U', 2, 0)] * 3 + [('U', 1, 0)] * 4 + [('U', 1, 1)] * 4
    CT.QUAD_ENABLE = True; CT.clear_all()
    t0 = time.time(); vq = CT.moment_transfer(atoms); tq = time.time() - t0
    CT.QUAD_ENABLE = False; CT.clear_all()
    t0 = time.time(); vd = CT.moment_transfer(atoms); td = time.time() - t0
    rel = abs(vq - vd) / max(abs(vd), 1e-300)
    print(f"n={n}: bottleneck moment quad {tq*1e3:.0f}ms vs DP {td*1e3:.0f}ms  ({td/max(tq,1e-9):.0f}x), rel={rel:.1e}", flush=True)
CT.QUAD_ENABLE = True
