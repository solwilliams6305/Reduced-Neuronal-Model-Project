"""STAGE 2 pipeline: build the self-consistent transfer ladder v0..v7 (Richardson-extrapolated
in 1/n) and run the Borel/Darboux resurgence analysis.  scipy-free (numpy + math only).

Auto-folds whatever grids are on disk:
  v0..v5  <- _ladder_transfer_results.pkl   (many grids -> Richardson)
  v6      <- _ckpt_6_{10,12,14}_cache.pkl    (assembled if complete)
  v7      <- _ckpt_7_{10,12,14}_cache.pkl    (assembled if complete; 1/n Richardson when >=2 grids)

Re-run any time; it picks up new _ckpt_7_* state files the moment they complete.

Borel tools (PROGRAM2_CONSOLIDATION.md 3ter):
  - Richardson 1/n extrapolation of grid values.
  - Hypergeometric / Meijer-G ratio approximants to the Borel coefficients (PRIMARY; Mera-Pedersen-
    Nikolic arXiv:1802.06034 style: converge at few orders, native branch cut).
  - Darboux/Dingle late-terms fit  v_n ~ 2C |zeta|^{-(n+1)} Gamma(n+1+alpha) cos((n+1)theta - phi)
    (complex-conjugate Borel pair), and pair + real-instanton two-singularity fit.
  - Borel-Pade poles as cross-check.
Nonlinear fits: grid the nonlinear params (theta, |zeta|, alpha, ...) and LINEAR-solve the amplitudes
(the model is linear in the cos/sin amplitudes once theta,zeta,alpha are fixed) -> robust, scipy-free.
"""
import os, sys, pickle
import numpy as np
from math import factorial, pi, gamma as Gamma, atan2, hypot

HERE = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------------------------------- ladder loading
def _vk_terms(k):
    m = k + 1
    pairs = [(m, m, 1)]
    for j in range(1, k + 1):
        a, b = j, 2 * m - j
        if a < b:
            pairs.append((a, b, 2))
    return pairs

def assemble_vk(k, n):
    """Assemble v_k(n) from a completed _ckpt_{k}_{n}_cache.pkl (+ moments list); None if incomplete."""
    cf = os.path.join(HERE, f'_ckpt_{k}_{n}_cache.pkl')
    mf = os.path.join(HERE, f'_ckpt_{k}_{n}_moments.pkl')
    if not (os.path.exists(cf) and os.path.exists(mf)):
        return None
    cache = pickle.load(open(cf, 'rb'))
    needed = pickle.load(open(mf, 'rb'))
    if any(m not in cache for m in needed):
        return None                                   # not finished yet
    import sympy as sp, chaos_diagram as CD, chaos_transfer as CT
    CD.setup(n_grid=n, MAXORD=CT.DENSE_MAXORD)
    u0ps = CD.GEO['u0ps']; Vst = CD.GEO['Vst']
    Us = [[sp.symbols(f'U{a}_{m}') for m in range(16)] for a in range(16)]
    s = [sp.symbols(f's{j}') for j in range(14)]; Vs = [sp.symbols(f'Vd{j}') for j in range(14)]
    Y = {idx: CD.Ybase(idx, Us, s, Vs, None, u0ps, Vst)
         for idx in sorted(set(p for pr in _vk_terms(k) for p in pr[:2]))}
    def mom(at): return cache[tuple(sorted(at))]
    def pm(A, B):
        t = 0.0
        for ca, aa in A:
            for cb, ab in B:
                if ca == 0.0 or cb == 0.0: continue
                t += ca * cb * mom(aa + ab)
        return t
    def mean(A): return sum(c * mom(at) for c, at in A)
    tot = 0.0
    for a, b, coef in _vk_terms(k):
        if a == b:
            mm = mean(Y[a]); val = pm(Y[a], Y[a]) - mm * mm
        else:
            val = pm(Y[a], Y[b]) - mean(Y[a]) * mean(Y[b])
        tot += coef * val
    return tot

def richardson(ns, vals, order=1):
    ns = np.asarray(ns, float); vals = np.asarray(vals, float)
    order = min(order, len(ns) - 1)
    cols = [np.ones_like(ns)] + [1.0 / ns ** p for p in range(1, order + 1)]
    A = np.vstack(cols).T
    coef, *_ = np.linalg.lstsq(A, vals, rcond=None)
    fit = A @ coef
    return coef[0], np.max(np.abs(fit - vals))

def load_ladder(verbose=True):
    """Return (v_inf list v0..vK, dict of per-k finite-grid data)."""
    res_f = os.path.join(HERE, '_ladder_transfer_results.pkl')
    try:
        res = pickle.load(open(res_f, 'rb'))
    except Exception:
        res = {}
    grid = {}                                          # k -> [(n, val)]
    for (k, n), v in res.items():
        grid.setdefault(k, []).append((n, v))
    # v5, v6, v7 from checkpoint caches (higher grids computed via _ckpt.py)
    for k in (5, 6, 7):
        for n in (10, 12, 14, 16, 18, 20, 24):
            v = assemble_vk(k, n)
            if v is not None:
                grid.setdefault(k, []).append((n, v))
    vinf = {}
    detail = {}
    for k in sorted(grid):
        pts = sorted(set(grid[k]))
        ns = [n for n, _ in pts]; vs = [v for _, v in pts]
        detail[k] = pts
        if len(ns) >= 3:
            vi, r2 = richardson(ns, vs, order=2)
            _, r1 = richardson(ns, vs, order=1)
            vinf[k] = vi; detail[('res', k)] = (r1, r2)
        elif len(ns) >= 2:
            vi, r1 = richardson(ns, vs, order=1); vinf[k] = vi; detail[('res', k)] = (r1, None)
        else:
            vinf[k] = vs[-1]; detail[('res', k)] = (None, None)
        if verbose:
            gstr = " ".join(f"n{n}={v:+.4f}" for n, v in pts)
            print(f"  v{k}: {gstr}   -> v_inf={vinf[k]:+.5f}", flush=True)
    K = max(vinf)
    return [vinf[k] for k in range(K + 1)], detail

# ----------------------------------------------------------------------------- Borel-Pade
def pade(c, L, M):
    c = np.asarray(c, float)
    A = np.zeros((M, M)); b = np.zeros(M)
    for i in range(1, M + 1):
        for j in range(1, M + 1):
            kk = L + i - j
            A[i - 1, j - 1] = c[kk] if kk >= 0 else 0.0
        b[i - 1] = -c[L + i]
    qtail = np.linalg.solve(A, b)
    q = np.concatenate([[1.0], qtail])
    p = np.array([sum(c[kk - j] * q[j] for j in range(0, min(kk, M) + 1)) for kk in range(L + 1)])
    return p, q

def borel_pade(v):
    print("  Borel-Pade poles of B(t)=sum v_n t^n/n!  (|zeta| @ theta deg):", flush=True)
    b = [v[n] / factorial(n) for n in range(len(v))]
    K = len(v) - 1
    for M in range(2, K + 1):
        L = K - M
        if L < 0: continue
        try:
            _, q = pade(b, L, M)
            roots = sorted(np.roots(q[::-1]), key=abs)
        except (np.linalg.LinAlgError, ValueError):
            continue
        s = ", ".join(f"{abs(r):.3f}@{np.degrees(np.angle(r)):+.1f}" for r in roots if abs(r) > 1e-9)
        print(f"    [{L}/{M}]: {s}", flush=True)

# ----------------------------------------------------------------------------- Darboux pair fit
def fit_darboux(v, nmin=0, thetas=None, zetas=None, alphas=None):
    """v_n ~ zeta^{-(n+1)} Gamma(n+1+alpha) [A cos((n+1)th) + B sin((n+1)th)].  Grid (th,zeta,alpha),
    linear-solve (A,B).  Returns best (theta,|zeta|,alpha,C,phi,resid)."""
    ns = np.arange(len(v), dtype=float)[nmin:]
    y = np.asarray(v, float)[nmin:]
    scale = np.abs(y) + 1e-3
    if thetas is None: thetas = np.radians(np.arange(20, 90, 0.5))
    if zetas is None: zetas = np.linspace(0.4, 3.0, 80)
    if alphas is None: alphas = np.linspace(-1.0, 1.5, 26)
    best = None
    for al in alphas:
        G = np.array([Gamma(n + 1.0 + al) for n in ns])
        for z in zetas:
            zp = z ** (-(ns + 1.0))
            for th in thetas:
                cc = np.cos((ns + 1.0) * th); ss = np.sin((ns + 1.0) * th)
                M = np.vstack([zp * G * cc / scale, zp * G * ss / scale]).T
                rhs = y / scale
                try:
                    ab, *_ = np.linalg.lstsq(M, rhs, rcond=None)
                except np.linalg.LinAlgError:
                    continue
                resid = np.max(np.abs(M @ ab - rhs))
                if best is None or resid < best[0]:
                    A, B = ab; C = hypot(A, B) / 2.0; phi = atan2(B, A)
                    best = (resid, th, z, al, C, phi)
    resid, th, z, al, C, phi = best
    return dict(theta=np.degrees(th), zeta=z, alpha=al, C=C, phi=np.degrees(phi), resid=resid,
                theta_rad=th)

# ----------------------------------------------------------------------------- two-singularity fit
def fit_two_sing(v, nmin=0, thetas=None, zetas=None, alphas=None, zrs=None):
    """pair (A cos + B sin) + real instanton (Cr).  Grid nonlinear (th,zeta,alpha,zr,ar), linear-solve
    (A,B,Cr)."""
    ns = np.arange(len(v), dtype=float)[nmin:]
    y = np.asarray(v, float)[nmin:]
    scale = np.abs(y) + 1e-3
    if thetas is None: thetas = np.radians(np.arange(30, 80, 1.0))
    if zetas is None: zetas = np.linspace(0.5, 3.0, 40)
    if alphas is None: alphas = np.linspace(-0.5, 1.0, 7)
    if zrs is None: zrs = np.linspace(0.8, 3.0, 30)
    best = None
    for al in alphas:
        G = np.array([Gamma(n + 1.0 + al) for n in ns])
        for ar in alphas:
            Gr = np.array([Gamma(n + 1.0 + ar) for n in ns])
            for z in zetas:
                zp = z ** (-(ns + 1.0))
                for zr in zrs:
                    zrp = zr ** (-(ns + 1.0))
                    for th in thetas:
                        cc = np.cos((ns + 1.0) * th); ss = np.sin((ns + 1.0) * th)
                        M = np.vstack([zp * G * cc / scale, zp * G * ss / scale, zrp * Gr / scale]).T
                        rhs = y / scale
                        try:
                            sol, *_ = np.linalg.lstsq(M, rhs, rcond=None)
                        except np.linalg.LinAlgError:
                            continue
                        resid = np.max(np.abs(M @ sol - rhs))
                        if best is None or resid < best[0]:
                            A, B, Cr = sol
                            best = (resid, th, z, al, hypot(A, B) / 2.0, atan2(B, A), Cr, zr, ar)
    resid, th, z, al, C, phi, Cr, zr, ar = best
    return dict(theta=np.degrees(th), zeta=z, alpha=al, C=C, Cr=Cr, zr=zr, ar=ar, resid=resid)

# ----------------------------------------------------------------------------- hypergeometric / Meijer-G
def hypergeom_ratio(v):
    """PRIMARY approximant.  Borel coeffs b_n=v_n/n!.  Fit the ratio r_n=b_{n+1}/b_n as a rational
    [1/1] in n:  r_n = (n+a)/(z (n+c)).  A single real pole gives t*=1/z; a complex-conjugate Borel
    pair makes the best real [1/1] fit RESIDUAL large and the recovered z complex-unstable -> we
    report both the [1/1] estimate and the modulus/phase from the growth+oscillation of b_n
    (the Meijer-G native-branch-cut readout)."""
    b = np.array([v[n] / factorial(n) for n in range(len(v))], float)
    nz = np.abs(b) > 1e-14
    idx = np.where(nz)[0]
    out = {}
    if len(idx) >= 3:
        r = b[1:] / b[:-1]
        ns = np.arange(len(r), dtype=float)
        good = np.isfinite(r)
        A = np.vstack([np.ones(good.sum()), -r[good] * ns[good], -r[good]]).T
        rhs = -ns[good]
        try:
            (a, z, zc), *_ = np.linalg.lstsq(A, rhs, rcond=None)
            tstar = 1.0 / z if abs(z) > 1e-12 else np.inf
            out['tstar_11'] = tstar
        except np.linalg.LinAlgError:
            out['tstar_11'] = None
    # modulus from |b_n|^{-1/n}, phase from sign-change period
    n_last = idx[-1]
    out['radius'] = abs(b[n_last]) ** (-1.0 / n_last) if n_last > 0 else np.inf
    signs = np.sign(b[nz])
    flips = np.sum(signs[1:] != signs[:-1])
    span = idx[-1] - idx[0]
    out['period'] = (2 * span / flips) if flips > 0 else np.inf     # ~ 2pi/theta
    out['theta_from_period'] = 360.0 / out['period'] if out['period'] not in (0, np.inf) else None
    return out

# ----------------------------------------------------------------------------- main
def run(verbose=True):
    print("=== STAGE 2: transfer ladder + Borel/Darboux analysis ===", flush=True)
    v, detail = load_ladder(verbose=verbose)
    K = len(v) - 1
    print(f"\nladder v0..v{K} (v_inf): {[f'{x:+.4f}' for x in v]}", flush=True)
    signs = "".join('+' if x >= 0 else '-' for x in v)
    print(f"sign pattern: {signs}", flush=True)
    print("\n--- Hypergeometric / Meijer-G ratio approximant (PRIMARY) ---", flush=True)
    hg = hypergeom_ratio(v)
    print(f"  [1/1] Borel pole t* = {hg.get('tstar_11')}", flush=True)
    print(f"  |zeta| from |b_n|^(-1/n) = {hg['radius']:.3f}   sign-flip period ~ {hg['period']:.2f} "
          f"=> theta ~ {hg['theta_from_period']}", flush=True)
    print("\n--- Darboux/Dingle complex-pair fit ---", flush=True)
    for nmin in (0, 1):
        d = fit_darboux(v, nmin=nmin)
        print(f"  n>={nmin}: theta={d['theta']:.1f}deg |zeta|={d['zeta']:.3f} alpha={d['alpha']:+.2f} "
              f"C={d['C']:.4f} phi={d['phi']:.1f}deg  max-rel-resid={d['resid']:.3f}", flush=True)
    if K >= 7:
        print("\n--- pair + real-instanton (two-singularity) fit [needs v7] ---", flush=True)
        t = fit_two_sing(v, nmin=0)
        print(f"  pair: theta={t['theta']:.1f}deg |zeta|={t['zeta']:.3f} alpha={t['alpha']:+.2f} C={t['C']:.4f}"
              f"  | real: zr={t['zr']:.3f} Cr={t['Cr']:+.4f} ar={t['ar']:+.2f}  max-rel-resid={t['resid']:.3f}",
              flush=True)
    print("\n--- Borel-Pade cross-check ---", flush=True)
    borel_pade(v)
    return v, detail

if __name__ == "__main__":
    run()
