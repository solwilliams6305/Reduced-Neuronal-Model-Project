"""
_tw_stokes_a_dependence.py -- the a-dependence of the TW_beta tail Stokes constant S(a).

RESULT (this note): S(a) has NO new transcendental a-dependence. Under the exact scaling
p=sqrt(a) P of the backward-Kolmogorov PDE (verified symbolically), the tail's dominant
boundary-escape reduces to an a-INDEPENDENT cubic-barrier problem in the single variable
    g := beta a^{3/2}   (eps_eff^2 = 1/g),
so all a-dependence of S(a) is the elementary Jacobian; the transcendental content is the
a-independent reduced Stokes constant S0 computed here.

Reduced problem: dP = -(P^2-1) dtau + 2 eps_eff dW, escape well(+1)->barrier(-1)->-inf.
Exact escape rate (mean-first-passage; the two Laplace peaks at P=-1,+1 decouple to all orders):
    R(g) = (1/pi) e^{-(2/3)g} / Qtilde(2/g)^2 ,
    Qtilde(D) = sum_{j>=0} qt_j D^j ,   qt_j = (6j)! / ( 576^j (3j)! (2j)! )   (exact rationals).
The reduced fluctuation series F(g)=Qtilde(2/g)^{-2}=sum_m d_m g^{-m} has factorial growth
    d_m ~ S0 * Gamma(m+gamma') * (3/2)^m         (Borel singularity at g=2/3=Phi; two-instanton).
We compute d_m exactly, extract (S0, gamma'), and compare S0 to the HM backbone
    C = sqrt(2/(3 pi^3)) = 0.146632271193848478...
"""
from fractions import Fraction as Fr
import mpmath as mp
import os

mp.mp.dps = 80
HERE = os.path.dirname(os.path.abspath(__file__))
C_HM = mp.sqrt(mp.mpf(2) / (3 * mp.pi**3))


# ---------- Qtilde coefficients as mpf via a running product (no giant integers) ----------
# qt_j = (6j)!/(576^j (3j)!(2j)!);  qt_j/qt_{j-1} = [prod_{i=0}^{5}(6j-i)] / (576*(3j)(3j-1)(3j-2)*(2j)(2j-1))
def qt_coeffs(J):
    out = [mp.mpf(1)]
    for j in range(1, J + 1):
        num = mp.mpf(1)
        for i in range(6):
            num *= (6 * j - i)
        den = mp.mpf(576) * (3*j) * (3*j - 1) * (3*j - 2) * (2*j) * (2*j - 1)
        out.append(out[-1] * num / den)
    return out


# ---------- series reciprocal-square: Qtilde^{-2} (mpf) ----------
def recip_square(q, M):
    s = [mp.mpf(0)] * (M + 1)          # reciprocal of Qtilde
    s[0] = mp.mpf(1)
    for m in range(1, M + 1):
        s[m] = -mp.fsum(q[k] * s[m - k] for k in range(1, m + 1))
    e = [mp.mpf(0)] * (M + 1)          # square of the reciprocal
    for m in range(M + 1):
        e[m] = mp.fsum(s[k] * s[m - k] for k in range(m + 1))
    return e


def neville_to_zero(xs, ys):
    n = len(xs)
    P = [mp.mpf(v) for v in ys]
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            P[i] = (xs[i] * P[i - 1] - xs[i - j] * P[i]) / (xs[i] - xs[i - j])
    return P[n - 1]


def to_mpf(x):
    return x   # already mpf


# ---------- cross-check R(g) against the direct MFPT double integral ----------
def R_direct(g):
    # tau = (1/D) int_{-inf}^{1} dy e^{V(y)/D} int_{y}^{inf} dz e^{-V(z)/D},  V=P^3/3 - P, D=2/g
    D = mp.mpf(2) / g
    V = lambda t: t**3 / 3 - t
    def inner(y):
        return mp.quad(lambda z: mp.e**(-V(z) / D), [y, 1, mp.inf])
    tau = (1 / D) * mp.quad(lambda y: mp.e**(V(y) / D) * inner(y), [-mp.inf, -1, 1])
    return 1 / tau


def R_formula(g, qt):
    D = mp.mpf(2) / g
    Q = mp.fsum(to_mpf(qt[j]) * D**j for j in range(len(qt)))   # asymptotic; truncate sensibly
    return (1 / mp.pi) * mp.e**(-mp.mpf(2) / 3 * g) / Q**2


def main():
    J = 170
    qt = qt_coeffs(J)
    print("# qt_j = (6j)!/(576^j (3j)!(2j)!)")
    print("  qt_0..3:", [mp.nstr(qt[j], 10) for j in range(4)], " (qt_1 should be 5/48=0.10416..)")
    assert abs(qt[1] - Fr(5, 48)) < mp.mpf(10)**(-40)

    # ---- (2) reduced fluctuation coefficients d_m and their large-order growth ----
    M = J
    e = recip_square(qt, M)                 # Qtilde^{-2} = sum e_m D^m
    dm = [e[m] * mp.mpf(2)**m for m in range(M + 1)]   # F(g)=sum d_m g^{-m}
    print("\n# (2) reduced fluctuation series F(g)=Qtilde(2/g)^{-2}=sum d_m g^{-m}")
    print("  d_0..4:", [mp.nstr(dm[m], 8) for m in range(5)])
    print("  signs d_1..12:", "".join("+" if dm[m] > 0 else "-" for m in range(1, 13)))

    rho = mp.mpf(3) / 2                        # analytic: q~_j Borel sing at D=4/3 => d_m base 3/2
    ratio = lambda m: dm[m] / dm[m - 1]
    # base check: ratio(m)-ratio(m-1) -> rho
    rho_est = lambda m: ratio(m) - ratio(m - 1)
    # gamma' estimator: ratio(m)/rho - (m-1) -> gamma'   (d_m ~ B Gamma(m+gamma') rho^m)
    gam_est = lambda m: ratio(m) / rho - (m - 1)
    sign = 1 if dm[M] > 0 else -1
    S0_est = lambda m, gam: sign * dm[m] / (mp.gamma(m + gam) * rho**m)

    # STABLE extrapolation: moderate Neville windows only (large windows are ill-conditioned).
    def extrap(fn, windows=(15, 20, 25, 30)):
        out = {}
        for W in windows:
            ks = list(range(M - W + 1, M + 1))
            xs = [mp.mpf(1) / k for k in ks]
            out[W] = neville_to_zero(xs, [fn(k) for k in ks])
        return out

    rvals = extrap(rho_est)
    print("  rho extrapolants (expect 3/2):", {W: mp.nstr(v, 12) for W, v in rvals.items()})
    gvals = extrap(gam_est)
    print("  gamma' extrapolants (expect 0):", {W: mp.nstr(v, 12) for W, v in gvals.items()})
    gstar = gvals[25]

    for gam_label, gam in (("gamma'=0", mp.mpf(0)), ("gamma'=gstar", gstar)):
        Svals = extrap(lambda m: S0_est(m, gam))
        S0 = Svals[25]
        print(f"\n  [{gam_label}] S0 extrapolants:", {W: mp.nstr(v, 14) for W, v in Svals.items()})
        print(f"  S0 = {mp.nstr(S0, 22)}   (sign {sign:+d})")
        print(f"    1/pi                      = {mp.nstr(1/mp.pi, 22)}   S0-1/pi = {mp.nstr(S0-1/mp.pi,3)}")
        print(f"    C_HM = sqrt(2/(3pi^3))     = {mp.nstr(C_HM, 16)}")
        print(f"    S0/C_HM   = {mp.nstr(S0/C_HM, 12)}   S0*pi = {mp.nstr(S0*mp.pi, 12)}   "
              f"S0*pi^2 = {mp.nstr(S0*mp.pi**2, 12)}")
        print(f"    S0*2pi = {mp.nstr(S0*2*mp.pi, 12)}   S0*4 = {mp.nstr(S0*4, 12)}   "
              f"S0*2pi^3/... C_HM/S0 = {mp.nstr(C_HM/S0, 12)}")
        for probe, val in (("S0", S0), ("S0*pi", S0*mp.pi), ("S0*pi^2", S0*mp.pi**2),
                           ("S0*2pi", S0*2*mp.pi)):
            try:
                idn = mp.identify(val, ['pi', 'sqrt(pi)', 'sqrt(2)', 'sqrt(3)', 'sqrt(6)'])
                if idn: print(f"    identify({probe}) = {idn}")
            except Exception:
                pass

    # ---- (3) sanity: escape-rate formula vs direct MFPT integral (low precision; MFPT=1 to run) ----
    if not os.environ.get("MFPT"):
        return
    print("\n# (3) R(g): formula vs direct MFPT double integral")
    with mp.workdps(30):
        print(f"   {'g':>5} {'R_direct':>16} {'R_formula':>16} {'rel diff':>10} {'K*':>4}")
        for g in [mp.mpf(v) for v in (10, 16)]:
            Rd = R_direct(g)
            D = mp.mpf(2) / g
            terms = [abs(qt[j] * D**j) for j in range(min(60, len(qt)))]
            Kmin = min(range(1, len(terms)), key=lambda k: terms[k])
            Q = mp.fsum(qt[j] * D**j for j in range(Kmin + 1))
            Rf = (1 / mp.pi) * mp.e**(-mp.mpf(2)/3 * g) / Q**2
            print(f"   {int(g):>5} {mp.nstr(Rd,10):>16} {mp.nstr(Rf,10):>16} "
                  f"{mp.nstr(abs(Rd-Rf)/Rd,3):>10} {Kmin:>4}")


if __name__ == "__main__":
    main()
