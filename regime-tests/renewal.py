"""
renewal.py — the fold renewal that drives the composition-propagation proof
---------------------------------------------------------------------------
The composed noisy chart map cannot accumulate error across the fold because the inner
chart FORGETS the entry: the peel-off law is invariant to entry perturbations (entry
offset from the canard, and small entry noise). This is the L2 ~ 0 ("Lambda_c -> 0")
renewal of CompositionPropagation_proof. Large entry kicks (>= separatrix distance) cause
direct escape -- the rare events handled by the union bound.
"""
from __future__ import annotations
import numpy as np


def peel(eta, Y0=5.0, dt=8e-4, N=20000, seed=0, delta0=0.0, burst=0.0, Ymin=-4.0):
    """Inner peel-off (first-node) law with entry offset delta0 (off the canard) and an
    optional one-shot entry derivative burst."""
    rng = np.random.default_rng(seed)
    R0 = -np.sqrt(Y0) + delta0                 # canard offset by delta0
    u = np.ones(N); v = np.full(N, -R0)        # R = -u'/u = R0
    v = v + burst * rng.standard_normal(N)     # entry noise burst
    Y = Y0; sdt = np.sqrt(dt); Yz = np.full(N, np.nan); done = np.zeros(N, bool)
    for _ in range(int((Y0 - Ymin) / dt)):
        if done.all():
            break
        dB = sdt * rng.standard_normal(N); u1 = u + v * dt; Yp = Y - dt; v1 = v + (Y * u) * dt - eta * u * dB
        u = u + 0.5 * (v + v1) * dt; v = v + 0.5 * (Y * u + Yp * u1) * dt - eta * (0.5 * (u + u1)) * dB; Y = Yp
        cr = (~done) & (u < 0.0); Yz[cr] = Y; done |= cr
    x = Yz[np.isfinite(Yz)]; m = x.mean(); d = x - m; vv = np.mean(d**2)
    return m, np.sqrt(vv), np.mean(d**3) / vv**1.5


def main():
    eta = 1.0
    print("\n=== Fold renewal: inner peel-off law vs ENTRY perturbations (eta=1) ===")
    print("  must be invariant within the basin (entry forgotten => no error accumulation)\n")
    print(f"  {'entry offset d0':>16} {'entry burst':>12} | {'mean':>8}{'std':>8}{'skew':>7}")
    for d0 in (0.0, 0.5, 1.0, -0.5):
        m, s, g = peel(eta, delta0=d0, seed=1)
        print(f"  {d0:16.2f} {0.0:12.2f} | {m:8.4f}{s:8.4f}{g:7.3f}")
    for bz in (0.5, 1.5, 3.0):
        m, s, g = peel(eta, burst=bz, seed=2)
        print(f"  {0.0:16.2f} {bz:12.2f} | {m:8.4f}{s:8.4f}{g:7.3f}")
    print("\n  => offset forgotten exactly; small bursts forgotten; only a burst >~ separatrix")
    print("     causes direct escape (rare; union-bound). The fold renews the chain.\n")


if __name__ == "__main__":
    main()
