"""
split_gnd_test.py — do the two derived FW tails (exponents 5 left, 3 right) PIN the cusp body?
=============================================================================================

Cusp tails: F~exp(−|s|⁵/20) (left, α=5), S~exp(−(4/3)|s|³) (right, α=3). The two-sided
"split generalized-normal" with these exponents,  f(z) ∝ exp(−|m−z|⁵/sL⁵) (z<m), exp(−|z−m|³/sR³)
(z≥m),  has ONE shape parameter after standardising (the scale ratio sR/sL; m and overall scale drop
out). Trace its (skew, exkurt) locus and ask: at the cusp skew +0.61, is the implied exkurt ≈ −0.24
(measured)? If yes, the two FW tails essentially determine the whole law. Instant (no MC).
"""
from __future__ import annotations
import numpy as np

z = np.linspace(-10, 10, 8000); dz = z[1] - z[0]


def std_moments(sL, sR):
    d = np.where(z < 0, np.exp(-np.abs(z)**5/sL**5), np.exp(-np.abs(z)**3/sR**3))
    d = d/(d.sum()*dz)
    mu = np.sum(z*d)*dz; var = np.sum((z - mu)**2*d)*dz; sd = np.sqrt(var)
    zz = (z - mu)/sd
    return np.sum(zz**3*d)*dz, np.sum(zz**4*d)*dz - 3.0


def main():
    print("=" * 72)
    print("Does the (5,3)-split generalized-normal reproduce the cusp law?")
    print("=" * 72)
    print("  (5,3) locus — vary the right/left scale ratio sR/sL (sL=1):")
    print(f"     {'sR':>5} {'skew':>8} {'exkurt':>8}")
    rows = []
    for sR in np.linspace(0.5, 2.4, 24):
        sk, ku = std_moments(1.0, sR); rows.append((sR, sk, ku))
        if abs(sR - round(sR, 1)) < 0.05 or True:
            pass
    for sR, sk, ku in rows[::3]:
        print(f"     {sR:5.2f} {sk:+8.3f} {ku:+8.3f}")

    # find the ratio giving cusp skew +0.611, read kurtosis
    sks = np.array([r[1] for r in rows]); kus = np.array([r[2] for r in rows])
    srs = np.array([r[0] for r in rows])
    j = int(np.argmin(np.abs(sks - 0.611)))
    print(f"\n  at cusp skew +0.611 (sR≈{srs[j]:.2f}):  (5,3)-split exkurt = {kus[j]:+.3f}   "
          f"(measured cusp κ4 = −0.24)")
    print(f"  ⇒ {'MATCH — the two FW tails essentially PIN the law' if abs(kus[j]+0.24) < 0.12 else 'partial — tails constrain but body has residual freedom (the PIV refinement)'}")


if __name__ == "__main__":
    main()
