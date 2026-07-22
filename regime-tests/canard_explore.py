#!/usr/bin/env python3
"""
canard_explore.py  — feasibility probe for the canard-strip sigma_crit(eps)~eps^{3/4}.

Two things to establish before any expensive sweep:
  (1) where the lower Hopf I_H1(eps) sits, and the Hopf frequency;
  (2) that we can actually SEE a canard / bifurcation delay numerically, and that
      noise shortens it -- the observable whose collapse defines sigma_crit.

We use the DYNAMIC-passage (slow current ramp) setup, which is the Berglund-Gentz
standard and avoids having to sit inside the exponentially-thin autonomous strip:
ramp I(t) slowly upward through I_H1 and record the delay = how far past I_H1 the
trajectory tracks the (now repelling) slow manifold before jumping to spike.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np

# --- shim scipy.optimize.brentq (sandbox has no scipy) ---
import types
def _brentq(f, lo, hi, xtol=1e-12, maxiter=200):
    flo, fhi = f(lo), f(hi)
    if flo == 0: return lo
    if fhi == 0: return hi
    if flo*fhi > 0: raise ValueError("f(lo) and f(hi) must bracket a root")
    for _ in range(maxiter):
        mid = 0.5*(lo+hi); fm = f(mid)
        if abs(hi-lo) < xtol or fm == 0: return mid
        if flo*fm < 0: hi, fhi = mid, fm
        else: lo, flo = mid, fm
    return 0.5*(lo+hi)
_scipy = types.ModuleType("scipy"); _opt = types.ModuleType("scipy.optimize")
_opt.brentq = _brentq; _scipy.optimize = _opt
sys.modules.setdefault("scipy", _scipy); sys.modules.setdefault("scipy.optimize", _opt)

from kernel import FHN2D

A, B = 0.7, 0.8


def passage_delay(eps, sigma, ramp, I_start_below, dt=2e-3, seed=0, n=400):
    """
    Slow upward ramp I(t)=I0 + ramp*t through I_H1. Return ensemble-mean current
    at the moment of spike (v>=1). Delay = I_spike - I_H1.
    """
    m = FHN2D(I=0.0, a=A, b=B)
    I_H1 = m.I_hopf_lower_at(eps)
    I0 = I_H1 - I_start_below
    rng = np.random.default_rng(seed)
    # start on the lower stable branch fixed point at I0
    m0 = FHN2D(I=I0, a=A, b=B)
    v = np.full(n, m0.V_FP); w = np.full(n, m0.W_FP)
    I = I0; t = 0.0; sdt = np.sqrt(dt)
    I_spike = np.full(n, np.nan); fired = np.zeros(n, bool)
    while I < I_H1 + 0.15 and not fired.all():
        noise = rng.standard_normal(n)
        v = v + (v - v**3/3.0 - w + I)*dt + sigma*sdt*noise
        w = w + eps*(v + A - B*w)*dt
        I = I0 + ramp*t; t += dt
        newly = (~fired) & (v >= 1.0)
        I_spike[newly] = I; fired |= newly
    return I_H1, float(np.nanmean(I_spike)), float(np.nanmean(fired))


def main():
    print("Lower-Hopf landmarks (a=0.7,b=0.8):")
    print("  eps     I_H1        omega      nascent_period")
    for eps in (0.01, 0.02, 0.05, 0.08):
        m = FHN2D(I=0.0, a=A, b=B)
        h = m.hopf_at(eps)
        print(f"  {eps:.3f}  {h['I_hopf_lower']:.6f}  {h['omega']:.5f}   {h['nascent_period']:.2f}")

    print("\nDynamic passage: delay = I_spike - I_H1 (deterministic vs noisy), eps=0.02")
    eps = 0.02; ramp = eps*0.5   # ramp slow relative to eps (slow passage)
    print(f"  ramp dI/dt = {ramp:.4f}")
    print("  sigma     I_H1      <I_spike>   delay      frac_fired")
    for sigma in (0.0, 0.005, 0.01, 0.02, 0.04, 0.08):
        I_H1, I_sp, frac = passage_delay(eps, sigma, ramp, I_start_below=0.04, seed=1)
        print(f"  {sigma:.3f}   {I_H1:.5f}   {I_sp:.5f}   {I_sp-I_H1:+.5f}    {frac:.2f}")


if __name__ == "__main__":
    main()
