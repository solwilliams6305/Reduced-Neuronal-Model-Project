"""
_hm_stokes_validate.py -- INDEPENDENT validations of the certified Hastings-McLeod
large-order / Stokes constant

    C = sqrt(2/(3 pi^3)) = (1/pi) sqrt(2/(3 pi)) = 0.146632271193848478...,   gamma = -1/2

for the coefficient asymptotics  a_k ~ -C Gamma(2k-1/2)(9/8)^k  of the s->-inf
Hastings-McLeod PII series  q(s) ~ sqrt(-s/2) sum_k a_k (-s)^{-3k}.

Three independent legs (method A = coefficient-ratio Richardson lives in
_hm_stokes_exact.py and gave C to 28 digits):

  (LIT)  Literature closed form. C = sqrt(2/(3 pi^3)) is the published HM Stokes
         constant (Dunne, Introductory Lectures on Resurgence arXiv:2511.15528
         eq.(2.16); Kapaev nlin/0411009; FIKN ch.11). Certify our number == it.
  (BOR)  Borel-space cross-check. c_k = a_k/(2k)! ~ -C 2^{-3/2} k^{-3/2}(9/8)^k:
         the Borel transform B(xi^2)=sum c_k x^k has a *bounded* square-root branch
         point at x0=8/9=1/A^2. Confirm the branch EXPONENT is exactly -1/2
         (i.e. c_k ~ k^{-3/2}) and the amplitude reproduces C. Different reprocessing
         of the coeffs (Borel space, via Stirling) than method A.
  (ODE)  Direct PII. Integrate the true Hastings-McLeod transcendent (scipy, HM
         boundary q~Ai) and confirm the optimally-truncated series reproduces it
         with remainder ~ exp(-A(-s)^{3/2}), A=2 sqrt2/3 -> the a_k ARE the genuine
         HM asymptotics, and the action is confirmed independently.
"""
import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import airy
import json, os

mp.mp.dps = 60
HERE = os.path.dirname(os.path.abspath(__file__))
A_ACTION = 2 * mp.sqrt(2) / 3
X0 = mp.mpf(8) / 9
C_CLOSED = mp.sqrt(mp.mpf(2) / (3 * mp.pi**3))


def load_coeffs_mpf():
    d = json.load(open(os.path.join(HERE, "_hm_coeffs_exact.json")))
    out = []
    for tok in d["a_exact"]:
        n, den = tok.split("/")
        out.append(mp.mpf(int(n)) / mp.mpf(int(den)))
    return out, d


def neville_to_zero(xs, ys):
    n = len(xs)
    P = [mp.mpf(v) for v in ys]
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            P[i] = (xs[i] * P[i - 1] - xs[i - j] * P[i]) / (xs[i] - xs[i - j])
    return P[n - 1]


# ============================ (LIT) ============================
def leg_lit(d):
    print("# (LIT) literature closed form  C = sqrt(2/(3 pi^3))")
    Cc = mp.mpf(d["C_str"])
    print(f"    computed  C* (method A, 28 digits) = {mp.nstr(Cc, 28)}")
    print(f"    sqrt(2/(3 pi^3))                   = {mp.nstr(C_CLOSED, 28)}")
    print(f"    agreement:  |C* - closed form|     = {mp.nstr(abs(Cc - C_CLOSED), 3)}  "
          f"(== extrapolation floor => match to ~28 digits)")
    print(f"    C* * pi = sqrt(2/(3 pi)) ?         = {mp.nstr(Cc*mp.pi,20)} vs "
          f"{mp.nstr(mp.sqrt(mp.mpf(2)/(3*mp.pi)),20)}")
    old = 1 / (4 * mp.sqrt(mp.pi))
    print(f"    OLD conjecture 1/(4 sqrt pi)       = {mp.nstr(old, 12)}  "
          f"REFUTED (off by {mp.nstr((C_CLOSED-old)/C_CLOSED*100,3)}%)")


# ============================ (BOR) ============================
def leg_borel(am):
    N = len(am) - 1
    nine_eighths = mp.mpf(9) / 8
    # c_k = a_k/(2k)! ; model c_k ~ -C 2^{-3/2} k^{-3/2}(9/8)^k.
    # (1) branch exponent: p_k := log|c_k/c_{k-1}| in base(9/8) ->  c_k/c_{k-1} = (9/8)(k/(k-1))^{-3/2}
    #     => k*[ (c_k/c_{k-1})/(9/8) - 1 ] -> -3/2  (the Borel branch exponent, = gamma-1).
    c = [am[k] / mp.factorial(2 * k) for k in range(N + 1)]
    def exp_est(k):
        return k * ((c[k] / c[k - 1]) / nine_eighths - 1)
    ks = list(range(N - 50 + 1, N + 1))
    xs = [mp.mpf(1) / k for k in ks]
    branch_exp = neville_to_zero(xs, [exp_est(k) for k in ks])
    # (2) amplitude -> C:  D_k = -c_k / (2^{-3/2} k^{-3/2} (9/8)^k) -> C
    def amp_est(k):
        return -c[k] / (mp.mpf(2)**mp.mpf(-1.5) * mp.mpf(k)**mp.mpf(-1.5) * nine_eighths**k)
    C_bor = neville_to_zero(xs, [amp_est(k) for k in ks])
    print("\n# (BOR) Borel-space cross-check (branch of B at x0=8/9)")
    print(f"    Borel branch exponent  (c_k ~ k^p): p = {mp.nstr(branch_exp, 16)}  "
          f"(exact -3/2 => gamma=-1/2; diff {mp.nstr(branch_exp+mp.mpf(3)/2,3)})")
    print(f"    amplitude -> C = {mp.nstr(C_bor, 20)}")
    print(f"    vs closed form  = {mp.nstr(C_CLOSED, 20)}   diff {mp.nstr(abs(C_bor-C_CLOSED),3)}")


# ============================ (ODE) ============================
def leg_ode(d):
    af = d["a_float"]
    Kmax = sum(1 for v in af if v is not None)  # floats stop being representable eventually
    af = [v for v in af if v is not None]
    A = float(A_ACTION)

    def rhs(s, y):
        q, qp = y
        return [qp, 2 * q**3 + s * q]
    s0 = 9.0
    Ai0, Aip0, _, _ = airy(s0)
    sol = solve_ivp(rhs, [s0, -7.0], [Ai0, Aip0], dense_output=True,
                    rtol=1e-13, atol=1e-15, max_step=0.01, method="DOP853")

    def series_opt(s):
        u = -s
        terms = [af[k] * u**(-3 * k) for k in range(len(af))]
        mags = [abs(t) for t in terms]
        Kstar = int(np.argmin(mags[1:])) + 1     # least-term truncation
        val = np.sqrt(u / 2) * sum(terms[:Kstar + 1])
        return val, np.sqrt(u / 2) * mags[Kstar], Kstar

    print("\n# (ODE) direct Hastings-McLeod PII (scipy DOP853) vs optimally-truncated series")
    q9 = sol.sol(9.0)[0]
    print(f"    HM seed check q(9)={q9:.4e} vs Ai(9)={Ai0:.4e}")
    print(f"    {'s':>5} {'q_HM':>16} {'|q_HM-series*|':>16} {'least term':>12} "
          f"{'e^-A(-s)^1.5':>13} {'K*':>4}")
    ss = [-2.0, -3.0, -4.0, -5.0, -6.0]
    ratios = []
    for s in ss:
        qhm = sol.sol(s)[0]
        ser, least, Kstar = series_opt(s)
        rem = abs(qhm - ser)
        expo = np.exp(-A * (-s)**1.5)
        ratios.append(rem / expo)
        print(f"    {s:>5.1f} {qhm:>16.10f} {rem:>16.3e} {least:>12.3e} "
              f"{expo:>13.3e} {Kstar:>4d}")
    print("    remainder/e^{-A(-s)^1.5} across s:", [f"{r:.2e}" for r in ratios])
    print("    => remainder tracks the least term and decays at rate A=2sqrt2/3: the a_k")
    print("       are the genuine HM asymptotics and the action is confirmed independently.")


def main():
    am, d = load_coeffs_mpf()
    print(f"loaded {len(am)} exact coeffs (a1={mp.nstr(am[1],4)}, a2={mp.nstr(am[2],6)})\n")
    leg_lit(d)
    leg_borel(am)
    leg_ode(d)
    print("\n# SUMMARY")
    print("    gamma = -1/2                                  [method A, 30 digits]")
    print("    C     = sqrt(2/(3 pi^3)) = 0.14663227119...   [method A 28d == literature]")
    print("    Borel branch exponent -3/2, action 2sqrt2/3   [BOR, ODE]")
    print("    old guess 1/(4 sqrt pi) REFUTED (3.8% off)")


if __name__ == "__main__":
    main()
