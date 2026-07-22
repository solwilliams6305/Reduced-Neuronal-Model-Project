"""Validation sweep: v3, v4, v5 at the same small grids used for v6 (n=10,12,14),
so the 1/n Richardson extrapolation can be validated against the known ladder
(v3~-0.03, v4~-0.45, v5~-1.1) before trusting it for v6."""
import time
import chaos_transfer as CT
from _v6_driver import compute_vk

for k in [3, 4, 5]:
    for n in [10, 12, 14]:
        t0 = time.time()
        v = compute_vk(k, n, CT.moment_hybrid, verbose=False)
        print(f"v{k}(n={n}) = {v:+.6f}   [{time.time()-t0:.0f}s]", flush=True)
