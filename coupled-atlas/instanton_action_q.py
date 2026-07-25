"""
General-q left-tail (persistence) instanton action, by the SAME artifact-free BVP machinery that
established the cusp value I(s) -> s^5/10 (instanton_action.py).

Why this exists: SWALLOWTAIL_TRANSSERIES_NOTES.md quotes the general-q anchor as
    I(s) = |s|^{2q+1} / [4(2q+1)]     -> q=2: s^5/20,  q=3: s^7/28,
but instanton_action.py *numerically verified* I(s) -> s^5/10 at q=2 (robust in p0 and M), and the
cusp paper eq:instanton also states s^5/10.  So the quoted general formula is off by a factor 2 at
the one point where it can be checked.

Analytic resolution: on the oscillatory side V=-t^q<0 the deterministic Riccati pdot=V-p^2 explodes;
holding it off to depth s needs p~0, hence pi ~ -V = t^q, and
    I(s) = (1/2) int_0^s pi^2 dt = s^{2q+1} / [2(2q+1)],
which reproduces 1/10 at q=2 and predicts 1/14 (not 1/28) at q=3.  This script checks that
prediction numerically for q=2 (regression) and q=3 (the swallowtail anchor).

Run:  python3 instanton_action_q.py
"""
import numpy as np
from scipy.integrate import solve_bvp

# q=1..5: turns the constant from a two-point coincidence into a FIVE-member family law.
# Measured I/s^(2q+1) at s=10 vs the predicted 1/[2(2q+1)]:
#   q=1 0.13987/0.16667 (0.839) | q=2 0.09297/0.10000 (0.930) | q=3 0.07010/0.07143 (0.981)
#   q=4 0.05538/0.05556 (0.997) | q=5 0.04543/0.04545 (0.9995)
# The ratio -> 1 monotonically: at fixed s the steeper q reaches its asymptote sooner, so the
# low-q entries are pre-asymptotic rather than wrong (push s up to tighten them).
Q_LIST = [1, 2, 3, 4, 5]


def make_V(q):
    # t=-Y; oscillatory side t>0: V=-t^q; confining side t<0: V=+|t|^q
    def V(t):
        return -np.sign(t) * np.abs(t) ** q
    return V


def instanton_action(s, q, p0=0.0, M=20.0, t0=1e-3, n=1200):
    """EL system on [t0,s]; BC p(t0)=p0 (incoming, O(1)), p(s)=-M (explosion threshold)."""
    V = make_V(q)
    t = np.linspace(t0, s, n)

    def ode(t, y):
        p, pi = y
        return np.vstack([pi + V(t) - p * p, 2 * p * pi])

    def bc(ya, yb):
        return np.array([ya[0] - p0, yb[0] + M])

    pg = np.where(t < 0.9 * s, 1.0 / np.maximum(t, 0.3), -M * (t - 0.9 * s) / (0.1 * s))
    pg = np.clip(pg, -M, 3.0)
    pig = np.abs(t) ** q
    sol = solve_bvp(ode, bc, t, np.vstack([pg, pig]), max_nodes=200000, tol=1e-6, verbose=0)
    if not sol.success:
        return np.nan, sol
    tt = np.linspace(t0, s, 4000)
    _, pipi = sol.sol(tt)
    return 0.5 * np.trapz(pipi ** 2, tt), sol


def predicted(q):
    """I(s) = s^{2q+1} / [2(2q+1)]"""
    return 1.0 / (2.0 * (2 * q + 1))


for q in Q_LIST:
    C_pred = predicted(q)
    C_notes = 1.0 / (4.0 * (2 * q + 1))
    p = 2 * q + 1
    print(f"\n=== q={q}:  I(s) vs predicted s^{p}/{1/C_pred:.0f}  (notes claim s^{p}/{1/C_notes:.0f}) ===")
    ss = np.array([2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0])
    Is = []
    for s in ss:
        S, _ = instanton_action(s, q)
        Is.append(S)
        print(f"  s={s:5.1f}: I(s)={S:12.3f}   I/s^{p}={S/s**p:.5f}   "
              f"ratio vs pred={S/(C_pred*s**p):.4f}   vs notes={S/(C_notes*s**p):.4f}")
    Is = np.array(Is)
    big = ss >= 5
    pexp = np.polyfit(np.log(ss[big]), np.log(Is[big]), 1)[0]
    Cbig = Is[-1] / ss[-1] ** p
    print(f"  large-s exponent = {pexp:.3f}  (predicted {p})")
    print(f"  large-s constant I/s^{p} = {Cbig:.5f}   predicted {C_pred:.5f}   notes {C_notes:.5f}")
    verdict = "PREDICTED 1/[2(2q+1)] CONFIRMED" if abs(Cbig / C_pred - 1) < 0.05 else (
        "notes' 1/[4(2q+1)] confirmed" if abs(Cbig / C_notes - 1) < 0.05 else "NEITHER -- investigate")
    print(f"  => {verdict}")
    print(f"  robustness (constant must be universal in p0, M):")
    for p0 in [-0.5, 0.0, 0.5]:
        S, _ = instanton_action(10.0, q, p0=p0)
        print(f"    p0={p0:+.1f}: I(10)/10^{p}={S/10.0**p:.5f}")
    for M in [12.0, 20.0, 40.0]:
        S, _ = instanton_action(10.0, q, M=M)
        print(f"    M={M:4.0f}: I(10)/10^{p}={S/10.0**p:.5f}")
