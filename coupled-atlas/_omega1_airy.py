"""Omega_1 for TW_beta via the Airy Green's function (sum-over-states), CORRECT normalization.

Omega_1 = 4 sum_{n>=1} <psi_0^2|psi_n^2>/(Lambda_0-Lambda_n),  psi_n = Ai(x-Lambda_n)/||.||,
Lambda_n = |n+1-th zero of Ai|.  Closed-form norm:  ||Ai(.-Lambda_n)||^2_{[0,inf)} = Ai'(-Lambda_n)^2
(identity int_z^inf Ai^2 = Ai'(z)^2 - z Ai(z)^2, and Ai(-Lambda_n)=0).  Terms ~ n^{-4/3} -> converge.

This is the O(1/beta) noise-averaged shift of the Airy connection root = -4 int psi_0^2 G_red(x,x),
the direct transpose of the cusp Omega_2 (= E[connection-root shift] at O(eta^2)).
"""
import numpy as np
from scipy.special import airy, ai_zeros

def compute(M=300, X=8.0, dx=1e-3):
    a, ap, ai_at_a, aip_at_ap = ai_zeros(M)   # a = zeros of Ai (negative)
    Lam = -a                                   # Lambda_n = |zero|  (Lambda_0=2.3381,...)
    norm2 = np.array([airy(-Lam[n])[1]**2 for n in range(M)])  # ||psi_n||^2 = Ai'(-Lambda_n)^2
    x = np.arange(dx/2, X, dx)
    # unnormalized psi_n on the overlap window (psi_0 lives in [0,6]; window [0,8] captures it)
    Aivals = np.array([airy(x - Lam[n])[0] for n in range(M)])
    p0sq = Aivals[0]**2
    ov = np.array([np.sum(p0sq * Aivals[n]**2)*dx for n in range(M)])   # int_0^X psi0u^2 psinu^2
    O = ov / (norm2[0]*norm2)                                            # normalized overlaps
    t = np.array([4.0*O[n]/(Lam[0]-Lam[n]) for n in range(1, M)])
    ns = np.arange(1, M)
    partial = np.cumsum(t)
    # tail: t_n ~ C n^{-p} (expect p~4/3); fit on upper half, sum_{M}^inf ~ C M^{1-p}/(p-1)
    m = ns >= M//2
    pf, lC = np.polyfit(np.log(ns[m]), np.log(-t[m]), 1)
    p = -pf; C = np.exp(lC); tail = -C*M**(1-p)/(p-1)
    return Lam, ns, t, partial, p, tail

if __name__ == "__main__":
    Lam, ns, t, partial, p, tail = compute(M=300)
    print(f"Lambda_0={Lam[0]:.5f}, Lambda_1={Lam[1]:.5f}")
    print("partial sums:", {c: round(float(partial[c-1]),5) for c in [10,50,100,200,299]})
    print(f"term decay exponent p = {p:.3f}  (expect ~4/3=1.33 -> convergent)")
    print(f"tail (n>=300) = {tail:+.5f}")
    Om1 = partial[-1] + tail
    print(f"\n==> Omega_1 = {partial[-1]:+.5f} + ({tail:+.5f}) = {Om1:+.4f}")
    print(f"    E[Lambda_0](beta) = {Lam[0]:.4f} + ({Om1:+.3f})/beta + O(1/beta^2)")
    # first few terms + robustness of tail window
    print("\nfirst terms:", [round(float(t[k]),4) for k in range(6)])
    for frac in (0.4,0.6,0.8):
        lo=int(300*frac); mm=ns>=lo; pf2,lC2=np.polyfit(np.log(ns[mm]),np.log(-t[mm]),1)
        pp=-pf2; CC=np.exp(lC2); tl=-CC*300**(1-pp)/(pp-1)
        print(f"   [window>={lo}: p={pp:.3f} tail={tl:+.4f} Omega_1={partial[-1]+tl:+.4f}]")
