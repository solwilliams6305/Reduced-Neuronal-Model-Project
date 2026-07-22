"""
blowdown.py — verify the blow-down noise scaling (BlowDown_proof)
----------------------------------------------------------------
PHYSICAL normal form near the fold:  dr=(r^2-y)dt+sigma dW,  dy=-eps2 dt.
Fold-weighted blow-up r=eps2^{1/3}R, y=eps2^{2/3}Y, T=eps2^{1/3}t maps this EXACTLY to the
inner Riccati dR=(R^2-Y)dT+(sigma/sqrt eps2)dB.  Hence the early-escape law depends on
(sigma,eps2) ONLY through eta=sigma/sqrt(eps2), and the peel-off rescales as
Y=y_esc/eps2^{2/3}.  VERIFY the collapse directly in the physical variables.
"""
from __future__ import annotations
import numpy as np


def passage(sigma, eps2, N=20000, seed=0, rcross=2.0):
    Y0 = 4.0; y0 = eps2**(2 / 3) * Y0
    r = -np.sqrt(y0) * np.ones(N); y = y0
    dT = 2e-3; dt = dT * eps2**(1 / 3)
    sdt = np.sqrt(dt); rng = np.random.default_rng(seed)
    yesc = np.full(N, np.nan); done = np.zeros(N, bool)
    n = int((y0 + eps2**(2 / 3) * 4.0) / eps2 / dt)
    rc = eps2**(1 / 3) * rcross
    for _ in range(n):
        al = ~done
        if not al.any():
            break
        r = r + (r**2 - y) * dt + sigma * sdt * rng.standard_normal(N)
        r = np.minimum(r, 1e3); y = y - eps2 * dt
        cr = al & (r >= rc); yesc[cr] = y; done |= cr
    Y = yesc[np.isfinite(yesc)] / eps2**(2 / 3)
    return np.mean(yesc > 0), Y.mean(), Y.std()


def main():
    print("\n=== Blow-down: physical normal form collapses in eta = sigma/sqrt(eps2) ===")
    print(f"  {'eta':>5}{'eps2':>8}{'sigma':>9} | {'P_early':>9}{'<Y_peel>':>10}{'std Y':>8}")
    for eta in (1.5, 1.0):
        for eps2 in (0.05, 0.20):
            sigma = eta * np.sqrt(eps2)
            P, mY, sY = passage(sigma, eps2, seed=1)
            print(f"  {eta:5.1f}{eps2:8.3f}{sigma:9.4f} | {P:9.4f}{mY:10.4f}{sY:8.4f}")
        print(f"     (rows at fixed eta match => collapse; the exact eta=sigma/sqrt(eps2) scaling)")
    print()


if __name__ == "__main__":
    main()
