"""Extension of the validation sweep: v3, v4, v5 at n=16,18,20 (the small-grid n<=14 values are
non-monotonic; find where the smooth regime starts; v3(n=20) has a documented value -0.060)."""
import time
import chaos_transfer as CT
from _v6_driver import compute_vk

for k in [3, 4, 5]:
    for n in [16, 18, 20]:
        t0 = time.time()
        v = compute_vk(k, n, CT.moment_hybrid, verbose=False)
        print(f"v{k}(n={n}) = {v:+.6f}   [{time.time()-t0:.0f}s]", flush=True)
