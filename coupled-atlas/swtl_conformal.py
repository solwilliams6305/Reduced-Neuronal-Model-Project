"""
Conformal-map Borel resummation for the swallowtail (q=3), with analytic instanton subtraction.

Why: plain Borel-Pade evaluates f(x)=sum v_n x^n at x=2 while the nearest Borel singularity sits at
|zeta| ~ 1.1, i.e. well outside the disc, and the ladder has a near-node at v4 that leaves a [3/3]
Pade badly conditioned.  A conformal map of the cut Borel plane onto the unit disc is the standard
fix: it pushes every singularity onto |w|=1 and turns a divergent/ill-conditioned extrapolation into
an interpolation inside the disc.  This buys accuracy WITHOUT new coefficients -- the relevant point
given v7 is compute-walled.

The map.  For a conjugate pair of branch points at zeta = A e^{+-i theta}, set
    sigma(u) = (1 - u/zeta)(1 - u/zeta_bar) = 1 - 2 (u/A) cos(theta) + (u/A)^2
    w(u)     = (1 - sqrt(sigma)) / (1 + sqrt(sigma))
so u=0 -> w=0 and both branch points (sigma=0) -> w=1.  The doubly-cut u-plane maps to |w|<1.
Inverting, with sigma(w) = ((1-w)/(1+w))^2,
    u(w) = A [ cos(theta) - sqrt(sigma(w) - sin^2(theta)) ]
(the branch that vanishes at w=0).  We re-expand B(u) as a series in w, Pade in w, and evaluate the
lateral Borel integrals through w(u).

Instanton subtraction.  The real Freidlin-Wentzell singularity sits on R+ at u = z_r (with
I(s)=s^7/14, so z_r <-> escape depth s=(14 z_r)^(1/7)).  Its contribution is removed analytically
before mapping, so the map only has to handle the complex pair.

VALIDATION FIRST: everything is checked on the cusp ladder, where the answer is published
(f(2)=0.235, truth 0.237), before being believed at q=3.

Run:  python3 swtl_conformal.py
"""
import numpy as np
from math import factorial, pi
from numpy.polynomial import polynomial as P

import swtl_borel as S

VAR_TRUTH_BETA2 = S.VAR_TRUTH_BETA2
F_TRUTH = VAR_TRUTH_BETA2 / 2.0
CUSP_TRUTH_F = {1.0: 0.232, 1.5: 0.239, 2.0: 0.237, 2.25: 0.234}


# ---------------------------------------------------------------- conformal map
def u_of_w_series(A, theta, N):
    """Power-series coefficients of u(w) = A[cos - sqrt(sigma(w) - sin^2)] about w=0, to order N."""
    # sigma(w) = ((1-w)/(1+w))^2 ; expand numerically via polynomial arithmetic on truncated series
    w = np.zeros(N + 1); w[1] = 1.0
    one = np.zeros(N + 1); one[0] = 1.0
    num = one - w              # 1 - w
    den = one + w              # 1 + w
    inv_den = _series_inv(den, N)
    r = _series_mul(num, inv_den, N)          # (1-w)/(1+w)
    sigma = _series_mul(r, r, N)
    c, s2 = np.cos(theta), np.sin(theta) ** 2
    arg = sigma.copy(); arg[0] -= s2          # sigma - sin^2  (constant term = cos^2)
    root = _series_sqrt(arg, N)               # sqrt -> constant term cos(theta)
    u = -root
    u[0] += c
    return A * u


def _series_inv(a, N):
    out = np.zeros(N + 1); out[0] = 1.0 / a[0]
    for n in range(1, N + 1):
        out[n] = -sum(a[k] * out[n - k] for k in range(1, n + 1)) / a[0]
    return out


def _series_mul(a, b, N):
    out = np.zeros(N + 1)
    for n in range(N + 1):
        out[n] = sum(a[k] * b[n - k] for k in range(n + 1))
    return out


def _series_sqrt(a, N):
    out = np.zeros(N + 1); out[0] = np.sqrt(a[0])
    for n in range(1, N + 1):
        out[n] = (a[n] - sum(out[k] * out[n - k] for k in range(1, n))) / (2 * out[0])
    return out


def _series_compose(b, u, N):
    """Given B(u)=sum b_n u^n and u(w)=sum u_m w^m (u[0]=0), return B(u(w)) to order N in w."""
    out = np.zeros(N + 1); out[0] = b[0]
    pw = np.zeros(N + 1); pw[0] = 1.0          # u(w)^0
    for n in range(1, len(b)):
        pw = _series_mul(pw, u, N)
        out += b[n] * pw
    return out


def w_of_u(u, A, theta):
    sigma = 1.0 - 2.0 * (u / A) * np.cos(theta) + (u / A) ** 2
    r = np.sqrt(sigma + 0j)
    return (1.0 - r) / (1.0 + r)


# ---------------------------------------------------------------- resummation
def conformal_borel(v, x, A, theta, L=None, M=None, phi_deg=35.0, zr=None, Cr=0.0, gamma=1.0,
                    Rmax=60.0, npts=8000):
    """Median (lateral) Borel sum of f(x) using a Pade in the conformal variable w.

    If zr is given, the real instanton pole C_r/(1-u/zr)^gamma is subtracted from B before mapping
    and added back analytically (its Borel integral is done on the same rotated contour).
    """
    K = len(v) - 1
    b = np.array([v[n] / factorial(n) for n in range(K + 1)])
    if zr is not None and Cr != 0.0:
        # subtract the instanton's Taylor coefficients: Cr*(1-u/zr)^-gamma = Cr*sum C(g+n-1,n) (u/zr)^n
        sub = np.array([Cr * _binom_rise(gamma, n) / zr ** n for n in range(K + 1)])
        b = b - sub
    N = K
    uw = u_of_w_series(A, theta, N)
    c = _series_compose(b, uw, N)
    if M is None:
        M = N // 2
    if L is None:
        L = N - M
    p, q = S.pade(c, L, M)

    phi = np.radians(phi_deg)
    r = np.linspace(1e-9, Rmax, npts)
    t = r * np.exp(1j * phi)
    u = x * t
    wv = w_of_u(u, A, theta)
    Bv = np.polyval(p[::-1], wv) / np.polyval(q[::-1], wv)
    if zr is not None and Cr != 0.0:
        Bv = Bv + Cr * (1.0 - u / zr) ** (-gamma)
    integ = np.exp(-t) * Bv * np.exp(1j * phi)
    return float(np.real(np.trapz(integ, r)))


def _binom_rise(g, n):
    """C(g+n-1, n) = g(g+1)...(g+n-1)/n!"""
    num = 1.0
    for j in range(n):
        num *= (g + j)
    return num / factorial(n)


# ---------------------------------------------------------------- drivers
def validate_cusp():
    print("=" * 78)
    print("VALIDATION: conformal Borel on the CUSP (published f(2)=0.235, truth 0.237)")
    print("=" * 78)
    print("  plain diagonal [3/3] (no map):")
    for x in (1.0, 1.5, 2.0, 2.25):
        f, _ = S.median_borel(S.CUSP_V, x, L=3, M=3, phi_deg=35.0)
        print(f"    x={x:.2f}: {f:.4f}   truth {CUSP_TRUTH_F[x]:.3f}")
    print("\n  conformal map with the paper's Borel data (A=1.9, theta=50 deg):")
    for x in (1.0, 1.5, 2.0, 2.25):
        f = conformal_borel(S.CUSP_V, x, A=1.9, theta=np.radians(50.0), phi_deg=35.0)
        t = CUSP_TRUTH_F[x]
        print(f"    x={x:.2f}: {f:.4f}   truth {t:.3f}   err {abs(f/t-1)*100:5.1f}%")
    print("\n  sensitivity of the map parameters at x=2 (truth 0.237):")
    print(f"    {'A':>5} " + " ".join(f"th={th:>3.0f}" for th in (40, 45, 50, 55, 60)))
    for A in (1.5, 1.7, 1.9, 2.1, 2.5):
        row = []
        for th in (40, 45, 50, 55, 60):
            try:
                row.append(conformal_borel(S.CUSP_V, 2.0, A=A, theta=np.radians(th), phi_deg=35.0))
            except Exception:
                row.append(np.nan)
        print(f"    {A:>5.1f} " + " ".join(f"{r:6.3f}" for r in row))


def apply_swallowtail(v):
    print("\n" + "=" * 78)
    print("SWALLOWTAIL q=3: conformal Borel (truth f(2)=%.4f)" % F_TRUTH)
    print("=" * 78)
    f_plain, _ = S.median_borel(v, 2.0, L=3, M=3, phi_deg=35.0)
    print(f"  plain diagonal [3/3], no map: f(2) = {f_plain:.4f}   (truth {F_TRUTH:.4f})")
    print(f"\n  conformal, scanning the (A, theta) the ladder actually supports")
    print(f"  (|zeta| ~ 1.1-1.5 from Borel-Pade; theta UNIDENTIFIED, so it is scanned):")
    print(f"    {'A':>5} " + " ".join(f"th={th:>3.0f}" for th in (35, 40, 45, 50, 55, 60)))
    grid = {}
    for A in (1.0, 1.1, 1.2, 1.3, 1.5):
        row = []
        for th in (35, 40, 45, 50, 55, 60):
            try:
                f = conformal_borel(v, 2.0, A=A, theta=np.radians(th), phi_deg=35.0)
            except Exception:
                f = np.nan
            row.append(f); grid[(A, th)] = f
        print(f"    {A:>5.1f} " + " ".join(f"{r:6.3f}" for r in row))
    vals = np.array([x for x in grid.values() if np.isfinite(x)])
    print(f"\n  across the scan: median {np.median(vals):.4f}, range [{vals.min():.4f}, {vals.max():.4f}]"
          f"   truth {F_TRUTH:.4f}")
    print(f"  spread/median = {(vals.max()-vals.min())/abs(np.median(vals))*100:.0f}%")
    return grid


if __name__ == "__main__":
    validate_cusp()
    v, band, cont = S.load_ladder(verbose=False)
    apply_swallowtail(cont)
