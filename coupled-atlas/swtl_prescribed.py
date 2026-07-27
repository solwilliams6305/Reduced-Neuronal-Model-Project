"""
Prescribed-singularity Borel resummation, using the imported deterministic anchor.

Why not just retry the conformal map.  The conjugate-pair conformal map failed for a STRUCTURAL
reason (swtl_conformal.py): the Borel plane has THREE singularities and the pair map handles two,
leaving the real Freidlin-Wentzell instanton strictly inside the unit disc (|w|=0.01-0.22), where
the mapped series cannot converge.  Knowing theta does not repair that -- the third singularity is
still inside.  What the new information DOES enable is a different and better-posed method.

The method.  Instead of letting a Pade approximant discover the singularities (which, with 8
coefficients and a near-node at v_4, it does badly), BUILD THEM IN.  Write

    B(u) = A_r (1-u/z_r)^{-g_r}  +  2 Re[ A_p (1-u/zeta)^{-g_p} ]  +  (polynomial remainder)

with the pair angle PRESCRIBED at theta = pi/(q+2) -- the deterministic anchor phase, 45 deg at
q=2 (where it is exactly known) and 36 deg at q=3 -- and only the amplitudes, exponents and the
remainder fitted to the known Taylor coefficients b_k = v_k/k!.  The Borel integral is then done on
a rotated contour, median of the two laterals, exactly as before.

This uses the anchor for what it is: an input that fixes WHERE the singularities are, which is
precisely the information the ladder could not supply and which both previously refuted methods
needed.

VALIDATION GATE: the cusp, where the published answer is f(2)=0.235 against truth 0.237, and where
plain diagonal [3/3] already achieves 0.2347.  A method that cannot at least match that is not worth
carrying to q=3.

Run:  python3 swtl_prescribed.py
"""
import numpy as np
from math import factorial, pi
from scipy.optimize import least_squares

import swtl_borel as S

CUSP_TRUTH = {1.0: 0.232, 1.5: 0.239, 2.0: 0.237, 2.25: 0.234}
SWTL_TRUTH = {1.0: 0.145188, 1.5: 0.160, 2.0: 0.164683}


def anchor_theta(q):
    """Deterministic resonance-root phase = pi/(q+2) = phi/4 with Voros's phi = 4pi/(q+2).
    Exact at q=2 (the cusp lambda_0 phase is exactly -pi/4)."""
    return pi / (q + 2)


def _rise(g, n):
    """Taylor coefficient of (1-z)^{-g}: C(g+n-1,n)."""
    r = 1.0
    for j in range(n):
        r *= (g + j) / (j + 1)
    return r


def build(params, K, theta, npoly):
    """Taylor coefficients (to order K) of the prescribed-singularity ansatz."""
    Ar, zr, gr, Apr, Api, zmod, gp = params[:7]
    poly = params[7:7 + npoly]
    zeta = zmod * np.exp(1j * theta)
    Ap = Apr + 1j * Api
    out = np.zeros(K + 1)
    for n in range(K + 1):
        val = Ar * _rise(gr, n) / zr ** n
        val += 2.0 * np.real(Ap * _rise(gp, n) / zeta ** n)
        if n < npoly:
            val += poly[n]
        out[n] = val
    return out


def fit(v, theta, npoly=3, zr0=2.0, zmod0=1.8):
    K = len(v) - 1
    b = np.array([v[n] / factorial(n) for n in range(K + 1)])
    sc = np.abs(b) + 1e-6

    def resid(p):
        return (build(p, K, theta, npoly) - b) / sc

    best = None
    for zr in (0.8 * zr0, zr0, 1.5 * zr0):
        for zm in (0.7 * zmod0, zmod0, 1.4 * zmod0):
            p0 = np.concatenate([[0.01, zr, 1.0, 0.01, 0.0, zm, 1.0], np.zeros(npoly)])
            lo = np.concatenate([[-50, 0.3, 0.05, -50, -50, 0.3, 0.05], -50 * np.ones(npoly)])
            hi = np.concatenate([[50, 20.0, 4.0, 50, 50, 20.0, 4.0], 50 * np.ones(npoly)])
            try:
                r = least_squares(resid, p0, bounds=(lo, hi), xtol=1e-15, ftol=1e-15, max_nfev=20000)
            except ValueError:
                continue
            if best is None or r.cost < best.cost:
                best = r
    return best.x, float(np.max(np.abs(best.fun)))


def resum(params, x, theta, npoly, phi_deg=35.0, Rmax=80.0, npts=12000):
    Ar, zr, gr, Apr, Api, zmod, gp = params[:7]
    poly = params[7:7 + npoly]
    zeta = zmod * np.exp(1j * theta)
    Ap = Apr + 1j * Api
    phi = np.radians(phi_deg)
    r = np.linspace(1e-9, Rmax, npts)
    t = r * np.exp(1j * phi)
    u = x * t
    B = Ar * (1.0 - u / zr) ** (-gr)
    B = B + Ap * (1.0 - u / zeta) ** (-gp) + np.conj(Ap) * (1.0 - u / np.conj(zeta)) ** (-gp)
    for n, c in enumerate(poly):
        B = B + c * u ** n
    return float(np.real(np.trapz(np.exp(-t) * B * np.exp(1j * phi), r)))


def run(name, v, q, truth, theta_deg=None):
    th = anchor_theta(q) if theta_deg is None else np.radians(theta_deg)
    print(f"\n{name}: prescribed pair angle theta = {np.degrees(th):.1f} deg "
          f"(= pi/(q+2), q={q})")
    for npoly in (2, 3, 4):
        p, res = fit(v, th, npoly=npoly)
        zr, zmod = p[1], p[5]
        line = f"  npoly={npoly}  fitresid={res:7.4f}  z_r={zr:6.3f} |zeta|={zmod:6.3f} :"
        for x in sorted(truth):
            f = resum(p, x, th, npoly)
            line += f"  f({x})={f:+.4f}[{(f/truth[x]-1)*100:+.0f}%]"
        print(line)


if __name__ == "__main__":
    print("=" * 78)
    print("GATE: cusp.  Plain diagonal [3/3] gives f(2)=0.2347 vs truth 0.237 (-1.0%).")
    print("=" * 78)
    run("cusp", S.CUSP_V, 2, CUSP_TRUTH)
    print("\n  (cusp anchor 45 deg is EXACT; measured stochastic pair is ~50 deg, so also try 50)")
    run("cusp", S.CUSP_V, 2, CUSP_TRUTH, theta_deg=50.0)

    v, band, cont = S.load_ladder(verbose=False)
    print("\n" + "=" * 78)
    print("SWALLOWTAIL (8 coefficients).  Plain diagonal [3/4] gives f(2)=0.2177 vs 0.1647 (+32%).")
    print("=" * 78)
    run("swallowtail", cont, 3, SWTL_TRUTH)
    print("\n  (and at the measured stochastic phase 39.4 deg)")
    run("swallowtail", cont, 3, SWTL_TRUTH, theta_deg=39.4)
