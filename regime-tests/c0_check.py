#!/usr/bin/env python3
"""
c0_check.py — numerical confirmation of the BGK canard-spacing constant c0.

Tests the leading-order result of MMO_C0_PINNED.md by root-finding the entry
point z0 of the k-th secondary canard from the EXACT twist condition (full
omega, no constant approximation) and computing c0 = z0^2/((2k+1)^2 mu^2):

    twist:  ∫_{z0}^0 omega(z) dz = pi*mu*(2k+1)/4 ,   omega(z) = sqrt(1+mu - z^2)
    dist  ~ exp(-z0^2/mu) = exp(-c0 (2k+1)^2 mu)

Confirms: (i) fixed-k limit c0 -> pi^2/16 as mu->0  (BGK Theorem 4.4 constant);
          (ii) the closed-form depth function c0(phi0); (iii) deep-canard c0 -> 1.

Reproduce: python3 regime-tests/c0_check.py
"""
import numpy as np
PI2_16 = np.pi**2/16.0

def _trap(y, x):                       # version-proof trapezoid
    return float(np.sum(0.5*(y[:-1]+y[1:])*np.diff(x)))
def omega(z, mu): return np.sqrt(np.maximum(1+mu - z*z, 0.0))
def twist_int(z0, mu, n=40000):
    zs = np.linspace(z0, 0.0, n); return _trap(omega(zs, mu), zs)
def find_z0(k, mu):
    target = np.pi*mu*(2*k+1)/4.0
    zmax = -np.sqrt(1+mu)*(1-1e-12)
    if twist_int(zmax, mu) < target: return None   # k beyond s_max for this mu
    lo, hi = zmax, 0.0
    for _ in range(200):
        mid = 0.5*(lo+hi)
        lo, hi = (mid, hi) if twist_int(mid, mu) > target else (lo, mid)
    return 0.5*(lo+hi)
def c0_closed(z0, mu):                 # closed-form depth function
    phi0 = np.arcsin(min(-z0/np.sqrt(1+mu), 1.0))
    return np.pi**2*np.sin(phi0)**2/((1+mu)*(2*phi0+np.sin(2*phi0))**2)

def main():
    print(f"pi^2/16 = {PI2_16:.5f}   (target shallow / fixed-k limit)\n")
    mus=[0.2,0.1,0.05,0.02,0.01,0.005]; ks=[0,1,2,3]
    print("c0(k,mu) from EXACT twist integral  ->  pi^2/16 = 0.6169 as mu->0 (fixed k):")
    print("   mu   " + "".join(f"    k={k}   " for k in ks))
    for mu in mus:
        row=f"{mu:7.3f}"
        for k in ks:
            z0=find_z0(k,mu)
            row += "    --    " if z0 is None else f"{z0**2/((2*k+1)**2*mu**2):10.5f}"
        print(row)
    print("\nclosed-form depth function c0(phi0) vs direct (k=1) — exact match:")
    for mu in mus:
        z0=find_z0(1,mu)
        if z0: print(f"  mu={mu:.3f}: direct={z0**2/(9*mu**2):.5f}  closed={c0_closed(z0,mu):.5f}")
    mu=0.02; smax=(1-mu)/(2*mu)
    print(f"\ndeep-canard (mu={mu}, s_max={smax:.0f}): c0 rises pi^2/16 -> 1 as k -> s_max")
    for k in [0,1,3,6,10,15,20,23]:
        z0=find_z0(k,mu)
        if z0 is None: print(f"  k={k:2d}: -- (beyond s_max)"); continue
        c0=z0**2/((2*k+1)**2*mu**2); phi=np.degrees(np.arcsin(-z0/np.sqrt(1+mu)))
        print(f"  k={k:2d}: c0={c0:.4f}   depth phi0={phi:5.1f} deg")

if __name__=="__main__": main()
