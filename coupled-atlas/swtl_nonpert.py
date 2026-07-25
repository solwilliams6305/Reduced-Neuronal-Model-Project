"""
ROUTE 1: read the NON-PERTURBATIVE sector straight off the observable, with no new rungs.

Context.  v7 is compute-expensive and v8 is doubly walled (symbolic + assembly), so the Borel data
(A, theta) cannot be pinned by pushing the ladder further.  But the trans-series says

    f(x) = [perturbative sector]  +  [non-perturbative sectors],
    f_np(x) ~ C_r x^{b_r} e^{-z_r/x}                                   (real FW instanton)
            + 2|C| x^{b} e^{-A cos(theta)/x} cos( A sin(theta)/x - phi )   (conjugate pair)

so the remainder R(x) = f_exact(x) - f_pert(x) DECAYS like e^{-A cos(theta)/x} and OSCILLATES in
1/x with frequency A sin(theta).  Fitting R over a range of x therefore yields A and theta directly
-- measured from the physics rather than inferred from seven large-order coefficients.

f_pert is the median (lateral) Borel sum of sum v_n x^n, which is reliable at small/moderate x where
the Borel-Pade integral is dominated by small u.  f_exact comes from the FP-PDE.

THE ACCURACY PROBLEM, AND THE FIX.  R is only a few percent of f, so f_exact must be very accurate.
The existing solver is first-order upwind, whose numerical diffusion |drift|*dp/2 is INDEPENDENT of
eta and swamps the physical D=eta^2/2 at small eta (this is the "small-eta bias" seen earlier: the
error halves as dp halves, 59%->30%->15%->8%->4%).  Brute-force dp refinement costs 4x per halving
(cells x steps).  Instead this module uses a SECOND-ORDER TVD (MUSCL + van Leer limiter)
reconstruction of the advective flux, which removes the O(dp) diffusion, with dp-Richardson kept as
a safety net.

VALIDATION GATE: at q=2, eta=0.3 the trusted cusp series is converged at f = 0.144770.  The solver
must reproduce that before any q=3 extraction is believed.

Run:  python3 swtl_nonpert.py [--gate] [--remainder]
"""
import sys
import numpy as np

import swtl_borel as S

CUSP_V = S.CUSP_V
CUSP_SERIES_ETA03 = 0.144770          # converged partial sum at x=0.09


def _van_leer(r):
    """van Leer limiter phi(r) = (r+|r|)/(1+|r|), phi(r<=0)=0."""
    a = np.abs(r)
    return np.where(np.isfinite(r), (r + a) / (1.0 + a), 2.0)


def solve_fp(q, eta, dp=0.01, dt=None, Y0=3.0, pmin=-9.0, pmax=9.0, tau_max=7.5, order=2):
    """Escape-location CDF via the Riccati FP equation.

    order=1 -> first-order upwind (the legacy scheme); order=2 -> MUSCL/van Leer TVD.
    Returns (y ascending, F(y)).
    """
    D = eta * eta / 2.0
    pc = np.arange(pmin + dp / 2, pmax, dp)
    pf = np.arange(pmin, pmax + dp / 2, dp)
    pf2 = pf * pf
    if dt is None:                       # keep advective CFL and diffusive number fixed
        dt = min(0.35 * dp / max(abs(pmin), abs(pmax)) ** 2, 0.20 * dp * dp / max(D, 1e-12))
    p0 = np.sqrt(max(np.sign(Y0) * abs(Y0) ** q, 1e-9))
    rho = np.exp(-0.5 * ((pc - p0) / 0.30) ** 2)
    rho /= rho.sum() * dp
    nst = int(round(tau_max / dt))
    taus = np.empty(nst + 1); Sv = np.empty(nst + 1)
    taus[0] = 0.0; Sv[0] = rho.sum() * dp
    eps = 1e-300
    for k in range(nst):
        tau = k * dt
        Y = Y0 - tau
        V = np.sign(Y) * abs(Y) ** q
        drift = V - pf2                                  # at faces, size N+1
        # ghost-padded cell values for reconstruction
        g = np.concatenate(([0.0, 0.0], rho, [0.0, 0.0]))   # index shift +2
        # for face i (between cell i-1 and cell i):  cells i-2,i-1,i,i+1 -> g[i], g[i+1], g[i+2], g[i+3]
        cm1, c0, c1, c2 = g[0:-3], g[1:-2], g[2:-1], g[3:]
        if order >= 2:
            # left state: reconstruct at right edge of cell i-1 (values c0 with neighbours cm1, c1)
            dL_up = c0 - cm1
            dL_dn = c1 - c0
            rL = np.where(np.abs(dL_up) > eps, dL_dn / np.where(np.abs(dL_up) > eps, dL_up, 1.0), 0.0)
            rho_L = c0 + 0.5 * _van_leer(rL) * dL_up
            # right state: reconstruct at left edge of cell i (values c1 with neighbours c0, c2)
            dR_up = c2 - c1
            dR_dn = c1 - c0
            rR = np.where(np.abs(dR_up) > eps, dR_dn / np.where(np.abs(dR_up) > eps, dR_up, 1.0), 0.0)
            rho_R = c1 - 0.5 * _van_leer(rR) * dR_up
            rho_L = np.maximum(rho_L, 0.0); rho_R = np.maximum(rho_R, 0.0)
        else:
            rho_L, rho_R = c0, c1
        adv = np.where(drift > 0, rho_L, rho_R) * drift
        diff = -D * (c1 - c0) / dp
        J = adv + diff
        J[-1] = 0.0                                      # reflecting on the right
        rho = rho - dt * (J[1:] - J[:-1]) / dp
        np.maximum(rho, 0.0, out=rho)
        taus[k + 1] = tau + dt
        Sv[k + 1] = rho.sum() * dp
    y = Y0 - taus
    return y[::-1], Sv[::-1]


def f_of_x(q, x, dp=0.01, order=2, **kw):
    """f(x) = Var(Y*)/eta^2 at x = eta^2."""
    import fp_cusp as FP
    eta = np.sqrt(x)
    y, F = solve_fp(q, eta, dp=dp, order=order, **kw)
    m1, sd, sk, ek, c5, c6 = FP.cumulants_from_F(y, F)
    return sd * sd / x


def f_richardson(q, x, dps=(0.02, 0.01), order=2, **kw):
    """Richardson-extrapolate f in dp.  order=2 -> assume O(dp^2); order=1 -> O(dp)."""
    vals = [f_of_x(q, x, dp=d, order=order, **kw) for d in dps]
    d0, d1 = dps[-2], dps[-1]
    p = 2.0 if order >= 2 else 1.0
    r = (d0 / d1) ** p
    return vals[-1] + (vals[-1] - vals[-2]) / (r - 1.0), vals


def gate():
    print("=" * 78)
    print("VALIDATION GATE: q=2, eta=0.3 (x=0.09).  Trusted cusp series = %.6f" % CUSP_SERIES_ETA03)
    print("=" * 78)
    x = 0.09
    print(f"  {'scheme':>18} {'dp':>8} {'f':>11} {'err':>9}")
    for order, tag in ((1, 'upwind (legacy)'), (2, 'MUSCL/van Leer')):
        for dp in (0.04, 0.02, 0.01):
            f = f_of_x(2, x, dp=dp, order=order)
            print(f"  {tag:>18} {dp:>8.3f} {f:>11.6f} {(f/CUSP_SERIES_ETA03-1)*100:>8.2f}%")
    fr, vals = f_richardson(2, x, dps=(0.02, 0.01), order=2)
    print(f"\n  2nd-order + Richardson(dp=0.02,0.01): f = {fr:.6f}   "
          f"err {(fr/CUSP_SERIES_ETA03-1)*100:+.2f}%")
    ok = abs(fr / CUSP_SERIES_ETA03 - 1) < 0.01
    print(f"  => {'PASS (<1%)' if ok else 'FAIL -- do not trust the extraction'}")
    return ok


if __name__ == "__main__":
    gate()
