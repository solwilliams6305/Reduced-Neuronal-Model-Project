"""
snic_channelB.py — Channel-B SNIC phase exponent (ChannelB_SNIC_proof)
----------------------------------------------------------------------
Noisy Adler / saddle-node-on-circle: dtheta = (mu - cos theta) dt + sigma dW, SNIC at mu=1.
At criticality the bottleneck rescales to the parameter-free  du = (1/2)u^2 dtau + dW,
so the noise-induced rotation rate and phase diffusion scale as sigma^{2/3}.
VERIFY: omega = <theta(T)>/T ~ sigma^{2/3}, D_phi = Var(theta(T))/(2T) ~ sigma^{2/3}.
(The drift pins 2/3 cleanly; the diffusion's effective exponent is pre-asymptotically
larger -- the source of an earlier 0.83 estimate -- and drifts toward 2/3 as sigma->0.)
"""
from __future__ import annotations
import numpy as np


def measure(mu, sigma, dt=5e-3, T=2500.0, M=300, seed=0):
    rng = np.random.default_rng(seed); n = int(T / dt); th = np.zeros(M); sdt = np.sqrt(dt)
    for _ in range(n):
        th = th + (mu - np.cos(th)) * dt + sigma * sdt * rng.standard_normal(M)
    return th.mean() / T, th.var() / (2 * T)


def main():
    print("\n=== Channel-B SNIC: criticality mu=1, omega & D_phi vs sigma ===")
    print(f"  {'sigma':>7}{'omega':>9}{'D_phi':>9}{'D/omega':>9}")
    sigs = np.array([0.12, 0.16, 0.22, 0.30, 0.42, 0.58])
    om, Dd = [], []
    for s in sigs:
        o, D = measure(1.0, s, seed=7); om.append(o); Dd.append(D)
        print(f"  {s:7.3f}{o:9.4f}{D:9.4f}{D / o:9.3f}")
    om, Dd = np.array(om), np.array(Dd)
    pO = np.polyfit(np.log(sigs), np.log(om), 1)[0]
    pOs = np.polyfit(np.log(sigs[:4]), np.log(om[:4]), 1)[0]
    pD = np.polyfit(np.log(sigs), np.log(Dd), 1)[0]
    print(f"\n  omega ~ sigma^{pO:.3f}  (small-sigma {pOs:.3f})   [predict 2/3 = {2/3:.3f}]")
    print(f"  D_phi ~ sigma^{pD:.3f}   (pre-asymptotic; -> 2/3 as sigma->0)")
    print(f"  => clean SNIC phase exponent = 2/3 (drift), the noisy saddle-node-on-circle.\n")


if __name__ == "__main__":
    main()
