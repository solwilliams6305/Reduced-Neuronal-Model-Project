"""Quad-accelerated validation ladder: v4 and v5 at larger grids, to pin their continuum values
(calibrates the grid extrapolation used for v6)."""
import time
import chaos_transfer as CT
from _v6_driver import compute_vk

for k in [4, 5]:
    for n in [16, 20, 24, 28]:
        t0 = time.time()
        v = compute_vk(k, n, CT.moment_hybrid, verbose=False)
        print(f"v{k}(n={n}) = {v:+.6f}   [{time.time()-t0:.0f}s]", flush=True)
