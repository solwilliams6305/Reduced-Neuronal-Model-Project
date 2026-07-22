"""Ladder extension: v4 and v5 at n=32 (v4 plateau confirmation; v5 continuum anchor)."""
import time
import chaos_transfer as CT
from _v6_driver import compute_vk

for k, n in [(4, 32), (5, 32)]:
    t0 = time.time()
    v = compute_vk(k, n, CT.moment_hybrid, verbose=False)
    print(f"v{k}(n={n}) = {v:+.6f}   [{time.time()-t0:.0f}s]", flush=True)
