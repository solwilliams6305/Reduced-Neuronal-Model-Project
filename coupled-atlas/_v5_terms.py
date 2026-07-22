"""Per-term v5 grid sequences (validation set for term-wise extrapolation: the assembled
continuum v5 = -1.19 +/- 0.01 is known from the n<=32 ladder)."""
import time
import chaos_transfer as CT
from _v6_driver import compute_vk

for n in [10, 12, 14, 16, 18, 20]:
    t0 = time.time()
    print(f"--- v5 n={n} ---", flush=True)
    v = compute_vk(5, n, CT.moment_hybrid, verbose=True)
    print(f"v5(n={n}) = {v:+.6f}   [{time.time()-t0:.0f}s]", flush=True)
