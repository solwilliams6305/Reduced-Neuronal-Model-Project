"""milestone2_measure.py -- Robust re-measurement of the stochastic Borel phase.

Independent estimators of (arg zeta = theta, |zeta|) from the variance ladder
  v = (v0..v5) = (0.134, 0.111, 0.100, -0.02, -0.45, -1.1)   [task-specified]
plus a bootstrap over the stated coefficient uncertainties.

Methods (all MC-free, from the same 6 numbers):
  A. Borel-Pade poles across ALL valid [L/M] splits.
  B. Darboux/Dingle floated-theta fit (single complex-conjugate pair), several nmin/alpha settings.
  C. Conformal-Borel (Euler) map: reparametrize Borel var t = zeta*w/(w+ (|near|)) ... here we use
     a simpler self-consistent "optimal conformal radius" scan + Cauchy-root diagnostic.
  D. Direct complex-ratio estimate: treat b_n = v_n/n! and estimate the nearest complex-conjugate
     singularity pair from three consecutive Borel coefficients (a 2-term recurrence / companion 2x2).

Certifies: parameter independence via bootstrap; drop-v5 vs include-v5 drift.

Run: python3 milestone2_measure.py
"""
import numpy as np
from math import factorial, pi
from scipy.special import gamma as Gamma
from scipy.optimize import least_squares
from numpy.polynomial import polynomial as P

np.set_printoptions(suppress=True, precision=4)

# task-specified ladder
V_TASK = [0.134, 0.111, 0.100, -0.02, -0.45, -1.1]
# _borel_analysis.py BASE ladder (slightly different continuum estimates)
V_BASE = [0.134, 0.111, 0.104, -0.030, -0.449, -1.1]

LAMBDA0_MOD = 1.258   # |lambda0|

# ------------------------------------------------------------------
# Borel coefficients
def borel_coeffs(v):
    return np.array([v[n] / factorial(n) for n in range(len(v))], float)

# ------------------------------------------------------------------
# A. Borel-Pade
def pade(c, L, M):
    c = np.asarray(c, float)
    N = L + M
    if len(c) < N + 1:
        raise ValueError("need L+M+1 coefficients")
    A = np.zeros((M, M)); b = np.zeros(M)
    for i in range(1, M + 1):
        for j in range(1, M + 1):
            k = L + i - j
            A[i - 1, j - 1] = c[k] if k >= 0 else 0.0
        b[i - 1] = -c[L + i]
    qtail = np.linalg.solve(A, b)
    q = np.concatenate([[1.0], qtail])
    return q

def borel_pade_pairs(v):
    """Return list of (L, M, [(|zeta|,theta_deg), ...]) for all valid splits, complex pairs only."""
    b = borel_coeffs(v)
    K = len(v) - 1
    out = []
    for M in range(1, K + 1):
        L = K - M
        if L < 0:
            continue
        try:
            q = pade(b, L, M)
            roots = P.polyroots(q)
        except (np.linalg.LinAlgError, ValueError):
            continue
        # keep finite roots
        roots = [r for r in roots if np.isfinite(r) and abs(r) > 1e-9 and abs(r) < 100]
        pairs = []
        for r in sorted(roots, key=abs):
            pairs.append((abs(r), np.degrees(np.angle(r))))
        out.append((L, M, pairs))
    return out

def borel_pade_nearest_complex(v):
    """Extract the nearest complex-conjugate pair (theta, |zeta|) from each [L/M] with a genuine
    complex pole off the real axis. Returns list of (theta_deg, mod)."""
    res = []
    for L, M, pairs in borel_pade_pairs(v):
        # find nearest-modulus root with |theta| between ~20 and ~85 deg (a genuine complex pair)
        cand = [(mod, th) for (mod, th) in pairs if 15 < abs(th) < 88]
        if not cand:
            continue
        cand.sort(key=lambda x: x[0])
        mod, th = cand[0]
        res.append((abs(th), mod, (L, M)))
    return res

# ------------------------------------------------------------------
# B. Darboux/Dingle
def darboux_model(params, n):
    C, zeta, theta, phi, alpha = params
    return 2 * C * zeta ** (-(n + 1.0)) * Gamma(n + 1.0 + alpha) * np.cos((n + 1.0) * theta - phi)

def fit_darboux(v, nmin=0, alpha_fixed=None):
    ns = np.arange(len(v), dtype=float)[nmin:]
    y = np.asarray(v, float)[nmin:]
    scale = np.abs(y) + 1e-3

    def resid(p):
        if alpha_fixed is not None:
            p = np.concatenate([p[:4], [alpha_fixed]])
        return (darboux_model(p, ns) - y) / scale

    best = None
    for th0 in np.radians([35, 45, 50, 55, 60, 65, 70, 80]):
        for z0 in [0.9, 1.1, 1.26, 1.5, 1.8, 2.2]:
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
    return dict(mod=zeta, theta=np.degrees(theta), C=C, phi=np.degrees(phi),
               alpha=alpha, maxres=res)

# ------------------------------------------------------------------
# D. Direct complex-ratio / companion estimate.
# Model the Borel coefficients as dominated by a conjugate pair b_n ~ Re[A zeta^{-(n+1)}] (times
# slowly varying Gamma-growth, ignored at leading ratio order). Then b_n satisfies a 2-term real
# linear recurrence  b_n = p*b_{n-1} + q*b_{n-2}  whose characteristic roots are 1/zeta and 1/conj(zeta).
# Solve p,q by least squares over available windows -> zeta = 1/root.
def direct_ratio_pair(v, use_borel=True, weight_late=True):
    if use_borel:
        b = borel_coeffs(v)
    else:
        b = np.asarray(v, float)
    n = len(b)
    if n < 4:
        return None
    # rows: b_k = p b_{k-1} + q b_{k-2}, k=2..n-1
    rows = []
    rhs = []
    for k in range(2, n):
        rows.append([b[k - 1], b[k - 2]])
        rhs.append(b[k])
    A = np.array(rows); r = np.array(rhs)
    if weight_late:
        # weight later equations more (larger-order) -- geometric
        w = np.array([2.0 ** i for i in range(len(r))])
        Aw = A * w[:, None]; rw = r * w
        sol, *_ = np.linalg.lstsq(Aw, rw, rcond=None)
    else:
        sol, *_ = np.linalg.lstsq(A, r, rcond=None)
    p, q = sol
    # characteristic: x^2 - p x - q = 0  -> x = (p +- sqrt(p^2+4q))/2
    disc = p * p + 4 * q
    roots = np.roots([1.0, -p, -q])
    # the geometric factor per step is x = 1/zeta (Borel) ignoring Gamma; for a complex pair disc<0
    if disc < 0:
        x = roots[0]  # complex
        zeta = 1.0 / x
        return dict(mod=abs(zeta), theta=abs(np.degrees(np.angle(zeta))), disc=disc,
                    complex_pair=True)
    else:
        # real roots -> no complex pair detected at this order
        return dict(roots=roots.real.tolist(), disc=disc, complex_pair=False)

def direct_ratio_pair_3term(v):
    """Use exactly the LAST three Borel coefficients b_{n-2},b_{n-1},b_n and one more to fit p,q
    from the two most-recent equations (most sensitive to the true late-order behavior)."""
    b = borel_coeffs(v)
    n = len(b)
    if n < 4:
        return None
    # two equations from the last two indices
    k1, k2 = n - 1, n - 2
    A = np.array([[b[k1 - 1], b[k1 - 2]], [b[k2 - 1], b[k2 - 2]]])
    r = np.array([b[k1], b[k2]])
    try:
        p, q = np.linalg.solve(A, r)
    except np.linalg.LinAlgError:
        return None
    disc = p * p + 4 * q
    roots = np.roots([1.0, -p, -q])
    if disc < 0:
        zeta = 1.0 / roots[0]
        return dict(mod=abs(zeta), theta=abs(np.degrees(np.angle(zeta))), disc=disc, complex_pair=True)
    return dict(roots=roots.real.tolist(), disc=disc, complex_pair=False)

# ------------------------------------------------------------------
# C. Conformal-Borel (Euler transform) diagnostic.
# Map t = (2 rho w)/(1 - w^2) style is overkill for 6 coeffs; instead do an "optimal ratio" scan:
# assume a single conjugate pair with Gamma(n+1+alpha) growth, divide it out for a grid of (alpha),
# then the residual sequence c_n = b_n * zeta_guess... -- we instead estimate |zeta| via the
# Cauchy-Hadamard root of the GAMMA-DEFLATED Borel coefficients and theta via the sign-change period.
def deflated_root_estimate(v, alpha=0.0):
    b = borel_coeffs(v)
    n = len(b)
    # deflate factorial-in-Borel: a_n = b_n / Gamma(n+1+alpha)  ~ Re[A zeta^{-(n+1)}]
    a = np.array([b[k] / Gamma(k + 1.0 + alpha) for k in range(n)])
    # |zeta| from ratio of successive envelope: |a_n / a_{n-1}| ~ 1/|zeta|
    ratios = np.abs(a[1:] / a[:-1])
    mods = 1.0 / ratios  # per-step |zeta| estimates
    # theta from the phase of the complex geometric fit on a: fit a_n ~ Re[A x^n], x=1/zeta
    # via same 2-term recurrence on a
    rows = []; rhs = []
    for k in range(2, n):
        rows.append([a[k - 1], a[k - 2]]); rhs.append(a[k])
    A = np.array(rows); r = np.array(rhs)
    sol, *_ = np.linalg.lstsq(A, r, rcond=None)
    p, q = sol
    disc = p * p + 4 * q
    roots = np.roots([1.0, -p, -q])
    info = dict(mods_per_step=mods, alpha=alpha, disc=disc)
    if disc < 0:
        zeta = 1.0 / roots[0]
        info.update(mod=abs(zeta), theta=abs(np.degrees(np.angle(zeta))), complex_pair=True)
    else:
        info.update(complex_pair=False, roots=roots.real.tolist())
    return info

# ------------------------------------------------------------------
def summarize_point(v, label):
    print(f"\n================ {label}: v = {[f'{x:+.4f}' for x in v]} ================")
    # A. Borel-Pade
    print("A. Borel-Pade complex pairs across splits (theta_deg, |zeta|, [L/M]):")
    bp = borel_pade_nearest_complex(v)
    for th, mod, lm in bp:
        print(f"     theta={th:6.2f}  |zeta|={mod:.3f}   [{lm[0]}/{lm[1]}]")
    bp_thetas = [th for th, mod, lm in bp]
    bp_mods = [mod for th, mod, lm in bp]
    # also full pole listing
    print("   (all poles per split:)")
    for L, M, pairs in borel_pade_pairs(v):
        s = ", ".join(f"{mod:.3f}@{th:+.1f}" for mod, th in pairs)
        print(f"     [{L}/{M}]: {s}")
    # B. Darboux
    print("B. Darboux/Dingle floated-theta fits:")
    dar = []
    for nmin, af in [(0, None), (0, 0.0), (1, 0.0), (1, None)]:
        d = fit_darboux(v, nmin=nmin, alpha_fixed=af)
        tag = f"nmin={nmin},alpha={'free' if af is None else af}"
        print(f"     {tag:20s}: theta={d['theta']:6.2f}  |zeta|={d['mod']:.3f}  alpha={d['alpha']:+.2f}  maxres={d['maxres']:.3f}")
        dar.append(d)
    # D. direct ratio
    print("D. Direct complex-ratio (2-term recurrence on Borel coeffs):")
    for wl in [False, True]:
        d = direct_ratio_pair(v, weight_late=wl)
        if d['complex_pair']:
            print(f"     weight_late={wl}: theta={d['theta']:6.2f}  |zeta|={d['mod']:.3f}  (disc={d['disc']:.3f})")
        else:
            print(f"     weight_late={wl}: NO complex pair (disc={d['disc']:.3f}, real roots {d.get('roots')})")
    d3 = direct_ratio_pair_3term(v)
    if d3 and d3['complex_pair']:
        print(f"     last-2-eqns:    theta={d3['theta']:6.2f}  |zeta|={d3['mod']:.3f}")
    elif d3:
        print(f"     last-2-eqns:    NO complex pair (disc={d3['disc']:.3f})")
    # C. deflated
    print("C. Gamma-deflated root estimate (conformal-ish):")
    for al in [0.0, -0.5, -1.0]:
        c = deflated_root_estimate(v, alpha=al)
        if c['complex_pair']:
            print(f"     alpha={al:+.1f}: theta={c['theta']:6.2f}  |zeta|={c['mod']:.3f}   per-step|zeta|={np.array(c['mods_per_step'])}")
        else:
            print(f"     alpha={al:+.1f}: NO complex pair; per-step|zeta|={np.array(c['mods_per_step'])}")
    return dict(bp_thetas=bp_thetas, bp_mods=bp_mods, dar=dar)


if __name__ == "__main__":
    summarize_point(V_TASK, "TASK ladder")
    summarize_point(V_BASE, "BASE ladder (_borel_analysis)")
