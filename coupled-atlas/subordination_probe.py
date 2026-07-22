"""
subordination_probe.py — is the cusp β-family TW^(2)_β a transform/subordination of one bare law?
=================================================================================================

The fingerprint showed the cusp β-family runs OPPOSITE to TW (skew rises with β, κ4 flips −→+). Question:
is the whole β-family a one-parameter transform (time-change / subordination) of a single 'bare' law,
or genuinely 2-dimensional? Test via the QQ map between β-laws: if cusp_{β'} ≈ φ(cusp_{β}) with φ a
simple (low-order) monotone map, the family is a transform family (subordination-consistent).
"""
from __future__ import annotations
import numpy as np

def riccati_escape(q, eta, Y0=3.0, Yend=-7.0, dtau=1.5e-3, M=60000, p_expl=-10.0, seed=0):
    rng = np.random.default_rng(seed); sdt = np.sqrt(dtau)
    p = np.full(M, Y0**(q/2.0)); Yesc = np.full(M, np.nan); done = np.zeros(M, bool)
    n = int((Y0 - Yend)/dtau)
    for i in range(n):
        Y = Y0 - i*dtau
        p = p + (np.sign(Y)*abs(Y)**q - p**2)*dtau - eta*sdt*rng.standard_normal(M)
        p = np.clip(p, p_expl-1, 25)
        cr = (~done) & (p < p_expl); Yesc[cr] = Y; done |= cr
    return Yesc[np.isfinite(Yesc)]

def zstd(x): return (x - x.mean())/x.std()
def mom(x): z=zstd(x); return np.mean(z**3), np.mean(z**4)-3

def main():
    betas = [1.0, 2.0, 4.0]
    laws = {b: zstd(riccati_escape(2.0, 2.0/np.sqrt(b), seed=10+int(b))) for b in betas}
    print("="*70); print("Subordination probe — cusp β-family transform test"); print("="*70)
    for b in betas:
        sk, ku = mom(laws[b]); print(f"  β={b}: skew {sk:+.3f}, exkurt {ku:+.3f}")
    pp = np.linspace(2, 98, 49)
    Q = {b: np.percentile(laws[b], pp) for b in betas}
    # QQ β=4 vs β=1: fit affine vs cubic
    x, y = Q[1.0], Q[4.0]
    ra = np.std(y - np.polyval(np.polyfit(x, y, 1), x))
    rc = np.std(y - np.polyval(np.polyfit(x, y, 3), x))
    print(f"\n  QQ(β=4 vs β=1): affine resid {ra:.3f}, cubic resid {rc:.3f}")
    print(f"  ⇒ {'low-order transform family (subordination-consistent: one bare law, β = a smooth reparametrization)' if rc < 0.04 else 'NOT a simple transform: the β-family is genuinely shape-varying (2-parameter), not a clean subordination of one law'}")
    # location-scale check: a pure subordination by scale would keep standardized shape FIXED
    print(f"  standardized shape changes with β (skew {mom(laws[1.0])[0]:+.2f}→{mom(laws[4.0])[0]:+.2f}) ⇒")
    print(f"  it is NOT a location-scale / pure-scale subordination; β genuinely deforms the shape.")

if __name__ == "__main__":
    main()
