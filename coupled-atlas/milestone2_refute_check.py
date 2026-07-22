"""Adversarial re-derivation of Mechanism B. Test whether -54 is robust or tunable."""
import numpy as np
from scipy.optimize import brentq

lam0 = 0.88956 - 0.88957j
C    = 0.43561 - 0.18049j
arg0 = np.angle(lam0)*180/np.pi

def zeta_saddle(x, n, steps=600):
    """x = eta2. saddle with continuous branch from x=0."""
    prev = None; sq = None
    for e in np.linspace(0.0, x, steps):
        disc = lam0**2 - 4.0*n*e*C
        s = np.sqrt(disc)
        if prev is not None and abs(s - prev) > abs(-s - prev):
            s = -s
        prev = s; sq = s
    return lam0 + (-lam0 + sq)/2.0

# The saddle depends ONLY on the product p = n*eta2. Verify:
print("=== Does zeta_saddle depend only on product p=n*eta2? ===")
for p in (0.3, 0.5, 0.514, 0.7, 1.0):
    vals = []
    for n in (3,4,5,6,7,8):
        z = zeta_saddle(p/n, n)
        vals.append(np.angle(z)*180/np.pi)
    print(f"  p={p:.3f}: arg across n=3..8 = {[f'{v:.2f}' for v in vals]}  spread={max(vals)-min(vals):.3f}")

# So arg is a MONOTONE function of p. Find p that gives various targets:
print("\n=== arg(zeta_eff) as function of p=n*eta2 (uses n=5 rep since p-invariant) ===")
def argp(p, n=5): return np.angle(zeta_saddle(p/n, n))*180/np.pi
for p in (0.0, 0.2, 0.4, 0.514, 0.6, 0.8, 1.0, 1.5):
    z = zeta_saddle(p/5, 5)
    print(f"  p={p:.3f}: arg={np.angle(z)*180/np.pi:+.2f} deg  |zeta|={abs(z):.4f}")

# Critical: the target -54 is reached at some p. But the ACTUAL Borel order n is unknown,
# and eta2 is a free parameter of the physics (the theory should predict arg AT a given small eta2).
# The claim "eta2=0.103 at n=5" is just: pick n=5, then p=0.514 forces eta2=0.103.
# Q: what does arg do as a function of eta2 at FIXED n, for the SMALL-eta2 limit that physics wants?
print("\n=== arg vs eta2 at fixed n (small-eta2 physical regime) ===")
for n in (5,6):
    print(f"  n={n}:")
    for eta2 in (0.01,0.02,0.05,0.10,0.20):
        z = zeta_saddle(eta2,n)
        print(f"    eta2={eta2:.2f}: arg={np.angle(z)*180/np.pi:+.2f}  |zeta|={abs(z):.4f}")

# The measured Borel target is -54 to -63 deg. Show arg sweeps THROUGH the whole band
# continuously as p grows -- so "hitting -54" is not a prediction, it's choosing where to stop.
print("\n=== p needed to hit each angle in the measured band -54..-63 ===")
for target in (-50,-54,-58,-63):
    try:
        p = brentq(lambda p: np.angle(zeta_saddle(p/5,5))*180/np.pi - target, 1e-3, 3.0)
        z = zeta_saddle(p/5,5)
        print(f"  arg={target}: p=n*eta2={p:.3f}  (eta2={p/5:.3f} at n=5)  |zeta|={abs(z):.4f}")
    except ValueError as e:
        print(f"  arg={target}: no root ({e})")
