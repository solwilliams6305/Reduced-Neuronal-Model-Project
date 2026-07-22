"""v6 across grids with the quad-accelerated hybrid engine (n<=14 overlap cross-checks the
old-engine runs; n>=16 is the regime the v3/v4/v5 validation showed is needed)."""
import time
import chaos_transfer as CT
from _v6_driver import compute_vk

for n in [10, 12, 14, 16, 18, 20]:
    t0 = time.time()
    print(f"--- v6 n={n} ---", flush=True)
    v = compute_vk(6, n, CT.moment_hybrid, verbose=True)
    print(f"v6(n={n}) = {v:+.6f}   [{time.time()-t0:.0f}s]", flush=True)
