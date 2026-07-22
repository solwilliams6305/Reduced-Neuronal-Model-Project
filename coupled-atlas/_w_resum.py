"""Median/lateral Borel-Pade resummation of the weak-noise series
      f(x) = Var(Y*)/eta^2 = sum_{n>=0} v_n x^n,   x = eta^2,
and comparison to the MC-free ground truth (fp_cusp) at each eta, in particular at
beta=2 (x=2) which lies OUTSIDE the perturbative radius (x_c ~ 1.2) where the naive
partial sums diverge.

Method:
  Borel transform B(u) = sum v_n u^n / n!  (finite radius |zeta|~1.9).
  [L/M] Pade of the truncated Borel series encodes the singularities (complex pair
  + real instanton) from the 7 coefficients.
  Borel sum   f(x) = ∫_0^∞ e^{-t} B(x t) dt.
  The real instanton pole sits on the positive t-axis -> the integral is ambiguous;
  LATERAL sums rotate the contour to angle ±phi, and the MEDIAN 1/2(f_+ + f_-) is
  the real (physical) value (imaginary parts cancel).
"""
import numpy as np
from scipy.interpolate import pade as scipy_pade

V = np.array([0.134, 0.111, 0.104, -0.030, -0.451, -1.19, -1.90])   # v_0..v_6

def borel_pade(v, L, M):
    """Return callable B(u) = [L/M] Pade of sum v_n u^n / n!."""
    from math import factorial
    b = np.array([v[n]/factorial(n) for n in range(len(v))])
    p, q = scipy_pade(b, M, L)          # scipy: pade(taylor, m[=den deg]) -> (num,den); order args
    return lambda u: p(u)/q(u), (p, q)

def borel_sum_lateral(v, x, phi_deg, L=3, M=3, R=60.0, Nt=6000):
    """Lateral Borel sum along contour t = r e^{i phi}, r in [0,R]."""
    B, _ = borel_pade(v, L, M)
    phi = np.deg2rad(phi_deg)
    r = np.linspace(1e-6, R, Nt)
    t = r*np.exp(1j*phi)
    integ = np.exp(-t) * B(x*t) * np.exp(1j*phi)     # dt = e^{i phi} dr
    return np.trapz(integ, r)

def median_resum(v, x, phi_deg=25.0, L=3, M=3, **kw):
    fp = borel_sum_lateral(v, x, +phi_deg, L, M, **kw)
    fm = borel_sum_lateral(v, x, -phi_deg, L, M, **kw)
    return 0.5*(fp+fm), fp, fm

def naive_partial(v, x):
    return np.cumsum([v[n]*x**n for n in range(len(v))])

if __name__ == "__main__":
    import json, sys, os
    print("Borel-Pade pole check (|zeta|, theta) for a few [L/M]:")
    from math import factorial
    b = np.array([V[n]/factorial(n) for n in range(len(V))])
    for (L, M) in [(3,3),(2,4),(4,2)]:
        try:
            p, q = scipy_pade(b, M, L)
            roots = q.roots
            for z in roots:
                if abs(z) < 20:
                    print(f"  [{L}/{M}]  |zeta|={abs(z):.3f}  theta={np.degrees(np.angle(z)):+.1f}")
        except Exception as e:
            print(f"  [{L}/{M}] failed: {e}")

    gt = {}
    if os.path.exists('_w_groundtruth.json'):
        for row in json.load(open('_w_groundtruth.json')):
            gt[round(row['eta2'],4)] = row['f']

    print("\n x=eta^2  beta   naive(7-term)   median-resum   ground-truth   (phi-stability)")
    for x in [0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 2.25]:
        npart = naive_partial(V, x)[-1]
        med, fp, fm = median_resum(V, x, phi_deg=25.0)
        med2, _, _ = median_resum(V, x, phi_deg=15.0)
        gtv = gt.get(round(x,4), None)
        gts = f"{gtv:.4f}" if gtv is not None else "   -   "
        beta = 4.0/x
        print(f" {x:5.2f}  {beta:5.2f}  {npart:+13.3f}  {med.real:+13.4f}  {gts:>12}   "
              f"(dphi: {abs(med.real-med2.real):.4f})", flush=True)
