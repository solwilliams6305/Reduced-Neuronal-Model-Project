"""Borel/Darboux large-order analysis of the weak-noise ladder v_0..v_k.

Tools (PROGRAM2_CONSOLIDATION.md §3ter):
  1. Richardson extrapolation of grid values v_k(n) in 1/n.
  2. Borel-Pade: poles of Pade approximants to B(t) = sum v_n t^n / n!.
  3. Darboux/Dingle late-terms fit: v_n ~ 2|C| |zeta|^{-n} Gamma(n+alpha) cos(n*theta - phi)
     (complex-conjugate Borel pair), optionally + real-instanton term C_r zr^{-n} Gamma(n+alpha).
  4. Hypergeometric (Mera-Pedersen-Nikolic style) ratio approximant to the Borel coefficients.

Run:  python3 _borel_analysis.py            # analysis on the baseline ladder (v6 placeholder excluded)
      python3 _borel_analysis.py <v6> [v7]  # include the new coefficient(s)
"""
import sys
import numpy as np
from math import factorial, pi
from scipy.special import gamma as Gamma
from scipy.optimize import least_squares
from numpy.polynomial import polynomial as P

# ---------------- Richardson ----------------
def richardson(ns, vals, order=1):
    """Fit vals(n) = v_inf + a/n (+ b/n^2 if order=2, needs >=3 pts). Returns v_inf and fit residual."""
    ns = np.asarray(ns, float); vals = np.asarray(vals, float)
    cols = [np.ones_like(ns)] + [1.0 / ns**p for p in range(1, order + 1)]
    A = np.vstack(cols).T
    coef, res, *_ = np.linalg.lstsq(A, vals, rcond=None)
    fit = A @ coef
    return coef[0], coef, np.max(np.abs(fit - vals))

# ---------------- Borel-Pade ----------------
def pade(c, L, M):
    """Pade [L/M] to series with coefficients c[0..L+M]. Returns (p,q) coefficient arrays, q[0]=1."""
    c = np.asarray(c, float)
    N = L + M
    if len(c) < N + 1:
        raise ValueError("need L+M+1 coefficients")
    # solve for q: sum_{j=1..M} q_j c_{L+i-j} = -c_{L+i}, i=1..M
    A = np.zeros((M, M)); b = np.zeros(M)
    for i in range(1, M + 1):
        for j in range(1, M + 1):
            k = L + i - j
            A[i - 1, j - 1] = c[k] if k >= 0 else 0.0
        b[i - 1] = -c[L + i]
    qtail = np.linalg.solve(A, b)
    q = np.concatenate([[1.0], qtail])
    p = np.array([sum(c[k - j] * q[j] for j in range(0, min(k, M) + 1)) for k in range(L + 1)])
    return p, q

def borel_pade_poles(v, L, M):
    """Poles of the [L/M] Pade approximant to the Borel transform of v."""
    b = [v[n] / factorial(n) for n in range(len(v))]
    p, q = pade(b, L, M)
    roots = P.polyroots(q)
    return sorted(roots, key=abs)

def report_borel_pade(v, label=""):
    print(f"--- Borel-Pade poles {label} (|zeta|, theta deg) ---")
    K = len(v) - 1
    for M in range(2, K):
        L = K - M
        if L < 0:
            continue
        try:
            roots = borel_pade_poles(v, L, M)
        except np.linalg.LinAlgError:
            continue
        s = ", ".join(f"{abs(r):.3f}@{np.degrees(np.angle(r)):+.1f}deg" for r in roots)
        print(f"  [{L}/{M}]: {s}")

# ---------------- Darboux/Dingle fit ----------------
def darboux_model(params, n):
    C, zeta, theta, phi, alpha = params
    return 2 * C * zeta ** (-(n + 1.0)) * Gamma(n + 1.0 + alpha) * np.cos((n + 1.0) * theta - phi)

def fit_darboux(v, nmin=0, alpha_fixed=None, verbose=True):
    """Fit v_n ~ 2C zeta^{-(n+1)} Gamma(n+1+alpha) cos((n+1)theta - phi) for n>=nmin."""
    ns = np.arange(len(v), dtype=float)[nmin:]
    y = np.asarray(v, float)[nmin:]
    scale = np.abs(y) + 1e-3   # relative residuals, envelope grows factorially

    def resid(p):
        if alpha_fixed is not None:
            p = np.concatenate([p[:4], [alpha_fixed]])
        return (darboux_model(p, ns) - y) / scale

    best = None
    for th0 in np.radians([40, 50, 55, 60, 65, 75]):
        for z0 in [1.0, 1.2, 1.5, 2.0]:
            p0 = [0.1, z0, th0, 0.5, 0.0]
            if alpha_fixed is not None:
                p0 = p0[:4]
            lo = [1e-6, 0.3, np.radians(5), -2 * pi] + ([] if alpha_fixed is not None else [-1.5])
            hi = [10.0, 5.0, np.radians(120), 2 * pi] + ([] if alpha_fixed is not None else [1.5])
            try:
                r = least_squares(resid, p0, bounds=(lo, hi), xtol=1e-14, ftol=1e-14)
            except ValueError:
                continue
            if best is None or r.cost < best.cost:
                best = r
    p = best.x if alpha_fixed is None else np.concatenate([best.x[:4], [alpha_fixed]])
    C, zeta, theta, phi, alpha = p
    res = np.max(np.abs(best.fun))
    if verbose:
        print(f"  Darboux pair fit (n>={nmin}, alpha {'free' if alpha_fixed is None else alpha_fixed}): "
              f"|zeta|={zeta:.3f} theta={np.degrees(theta):.1f}deg C={C:.4f} phi={np.degrees(phi):.1f}deg "
              f"alpha={alpha:+.3f}  max-rel-resid={res:.3f}")
        pred = darboux_model(p, np.arange(len(v) + 2, dtype=float))
        print(f"    model vs data: " + " ".join(f"{m:+.3f}/{d:+.3f}" for m, d in zip(pred, v)))
        print(f"    predicts v{len(v)}={pred[len(v)]:+.3f}, v{len(v)+1}={pred[len(v)+1]:+.3f}")
    return p, res

def darboux_two_sing(params, n):
    """Complex pair + real instanton: shared Gamma growth, alpha per singularity."""
    C, zeta, theta, phi, alpha, Cr, zr, ar = params
    pair = 2 * C * zeta ** (-(n + 1.0)) * Gamma(n + 1.0 + alpha) * np.cos((n + 1.0) * theta - phi)
    real = Cr * zr ** (-(n + 1.0)) * Gamma(n + 1.0 + ar)
    return pair + real

def fit_two_sing(v, nmin=0, verbose=True):
    ns = np.arange(len(v), dtype=float)[nmin:]
    y = np.asarray(v, float)[nmin:]
    scale = np.abs(y) + 1e-3

    def resid(p):
        return (darboux_two_sing(p, ns) - y) / scale

    best = None
    for th0 in np.radians([50, 55, 60, 65]):
        for sr in [+1, -1]:
            p0 = [0.1, 1.3, th0, 0.5, 0.0, sr * 0.05, 1.5, 0.0]
            lo = [1e-6, 0.3, np.radians(5), -2 * pi, -1.5, -10, 0.3, -1.5]
            hi = [10.0, 5.0, np.radians(120), 2 * pi, 1.5, 10, 5.0, 1.5]
            try:
                r = least_squares(resid, p0, bounds=(lo, hi), xtol=1e-14, ftol=1e-14)
            except ValueError:
                continue
            if best is None or r.cost < best.cost:
                best = r
    C, zeta, theta, phi, alpha, Cr, zr, ar = best.x
    res = np.max(np.abs(best.fun))
    if verbose:
        print(f"  Two-singularity fit (n>={nmin}): pair |zeta|={zeta:.3f} theta={np.degrees(theta):.1f}deg "
              f"alpha={alpha:+.2f} | real zr={zr:.3f} Cr={Cr:+.4f} ar={ar:+.2f}  max-rel-resid={res:.3f}")
    return best.x, res

# ---------------- hypergeometric ratio approximant ----------------
def hypergeom_ratio(v, verbose=True):
    """Fit Borel-coefficient ratios b_{n+1}/b_n as rational [1/1] in n (2F1-type ansatz):
    r_n = (n + a)/(z (n + b)) => Borel singularity at t = z. Complex pair shows as failure of the
    real fit; report the effective |zeta| from |b_n|^{-1/n} and oscillation period from sign pattern."""
    b = np.array([v[n] / factorial(n) for n in range(len(v))])
    r = b[1:] / b[:-1]
    ns = np.arange(len(r), dtype=float)
    # linear system: r_n * z * (n + bb) = n + a  ->  unknowns (a, z, z*bb)
    A = np.vstack([np.ones_like(ns), -r * ns, -r]).T
    rhs = -ns
    (a, z, zbb), *_ = np.linalg.lstsq(A, rhs, rcond=None)
    if verbose:
        env = np.abs(b[np.abs(b) > 1e-12])
        radius = np.abs(b[-1]) ** (-1.0 / (len(b) - 1))
        print(f"  Hypergeom [1/1] ratio fit: singularity t={1.0/z if z!=0 else np.inf:.3f} "
              f"(real-axis ansatz; complex pair -> poor fit expected); crude radius from |b_n|^(-1/n): {radius:.3f}")
    return a, z, zbb

# ---------------- main ----------------
BASE = [0.134, 0.111, 0.104, -0.030, -0.449, -1.1]   # v0..v5 (continuum/best estimates)

if __name__ == "__main__":
    v = list(BASE)
    args = [float(x) for x in sys.argv[1:]]
    v += args
    K = len(v) - 1
    print(f"ladder v0..v{K}: {[f'{x:+.3f}' for x in v]}")
    report_borel_pade(v)
    print("--- Darboux/Dingle single complex pair ---")
    fit_darboux(v, nmin=0, alpha_fixed=0.0)
    fit_darboux(v, nmin=0, alpha_fixed=None)
    fit_darboux(v, nmin=1, alpha_fixed=0.0)
    if K >= 7:
        print("--- pair + real instanton (needs >=8 data) ---")
        fit_two_sing(v, nmin=0)
    print("--- hypergeometric ratio diagnostic ---")
    hypergeom_ratio(v)
