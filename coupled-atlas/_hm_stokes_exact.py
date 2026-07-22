"""
_hm_stokes_exact.py  -- EXACT large-order / Stokes constant of the Hastings-McLeod
Painlevé II trans-series backbone (TW_beta halo, Theorem T2c).

Upgrades _hm_stokes.py (which took N=10 sympy coeffs and reported the *guess*
C ~ 0.139 ~ 1/(4 sqrt pi), gamma ~ -1/2) to a CERTIFIED, high-order result.

Model (Hastings-McLeod, alpha=0):   q'' = 2 q^3 + s q,
    q(s) ~ sqrt(-s/2) * sum_{k>=0} a_k (-s)^{-3k}     as s -> -inf.
The tail law inherits these coeffs by Theorem T1 (Lax-pair map to PII).

Claim to certify:   a_k ~ -C * Gamma(2k + gamma) * (9/8)^k ,
    with 9/8 = 1/A^2, action A = 2 sqrt(2)/3  (Borel singularity at xi = A in
    the variable xi = (-s)^{3/2}); gamma =? -1/2 ; C =? 1/(4 sqrt pi).

Method A  : exact rational recursion -> Neville/Richardson extrapolation of
            (gamma, C) with stability error bars.
(Methods B/ODE live in _hm_stokes_validate.py.)
"""
from fractions import Fraction as Fr
import mpmath as mp
import json, os

mp.mp.dps = 80
HERE = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------------------
# Exact coefficients via the closed-form rational recursion.
#
# Substituting q = 2^{-1/2} sum_k a_k u^{1/2-3k}  (u=-s) into  q_uu = 2q^3 - u q
# and collecting the coefficient of u^{3/2-3m} gives, for m>=1,
#
#     a_m = [ a_{m-1} * p * (p-1)  -  T(m) ] / 2 ,   p = 1/2 - 3(m-1),
#     T(m) = sum_{i+j+l=m, 0<=i,j,l<=m-1} a_i a_j a_l    (cubic self-convolution,
#                                                          m-th index excluded).
# a_0 = 1. Derivation + hand-checks (a_1=-1/8, a_2=-73/128) in the notes.
# ---------------------------------------------------------------------------
def hm_coeffs_exact(N):
    a = [Fr(0)] * (N + 1)
    a[0] = Fr(1)
    for m in range(1, N + 1):
        p = Fr(1, 2) - 3 * (m - 1)
        drive = a[m - 1] * p * (p - 1)
        T = Fr(0)
        for i in range(m):            # i = 0..m-1
            ai = a[i]
            if ai == 0:
                continue
            for j in range(m):        # j = 0..m-1
                l = m - i - j
                if 0 <= l <= m - 1:
                    T += ai * a[j] * a[l]
        a[m] = (drive - T) / 2
    return a


def to_mpf(fr):
    return mp.mpf(fr.numerator) / mp.mpf(fr.denominator)


# ---------------------------------------------------------------------------
# Neville extrapolation of a sequence f_k (nodes x_k = 1/k) to x = 0.
# Works for f_k = L + c1/k + c2/k^2 + ...  (integer-power asymptotic series).
# ---------------------------------------------------------------------------
def neville_to_zero(xs, ys):
    n = len(xs)
    P = [mp.mpf(v) for v in ys]
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            # extrapolate to x=0:
            P[i] = (xs[i] * P[i - 1] - xs[i - j] * P[i]) / (xs[i] - xs[i - j])
    return P[n - 1]


def main():
    N = 120
    print(f"# Hastings-McLeod PII s->-inf coefficients, exact recursion, N={N}")
    a = hm_coeffs_exact(N)

    # ---- certification 1: known low-order values ----
    assert a[0] == Fr(1),          f"a0={a[0]}"
    assert a[1] == Fr(-1, 8),      f"a1={a[1]}"
    assert a[2] == Fr(-73, 128),   f"a2={a[2]}"
    print("cert: a0=1, a1=-1/8, a2=-73/128  OK")

    # ---- certification 2: reproduce existing _hm_coeffs.json (from N=10 sympy) ----
    jpath = os.path.join(HERE, "_hm_coeffs.json")
    if os.path.exists(jpath):
        js = json.load(open(jpath))
        maxerr = max(abs(float(a[k]) - js[k]) / max(1.0, abs(js[k]))
                     for k in range(len(js)))
        print(f"cert: matches _hm_coeffs.json (k<{len(js)}) to rel {maxerr:.1e}  "
              f"{'OK' if maxerr < 1e-9 else 'FAIL'}")

    # ---- signs (Cleri-Dunne non-alternating signature) ----
    signs = "".join("+" if a[k] > 0 else "-" for k in range(1, 11))
    print(f"signs a1..a10: {signs}  (all '-' => positive-axis Borel singularity)")

    am = [to_mpf(x) for x in a]
    nine_eighths = mp.mpf(9) / 8

    # ======================= METHOD A ==============================
    # gamma estimator from consecutive ratios:
    #   rho_k = (a_k/a_{k-1})/(9/8) = (2k+gamma-1)(2k+gamma-2)(1+O(1/k))
    #   => 2k+gamma = [3 + sqrt(1+4 rho_k)]/2  =>  gamma_k -> gamma
    def gamma_est(k):
        rho = (am[k] / am[k - 1]) / nine_eighths
        return (3 + mp.sqrt(1 + 4 * rho)) / 2 - 2 * k

    print("\n# METHOD A -- Richardson/Neville extrapolation")
    print("gamma_k (raw, tail):",
          [mp.nstr(gamma_est(k), 8) for k in range(N - 4, N + 1)])

    # extrapolate gamma using several window sizes for a stability bar
    gamma_vals = {}
    for M in (20, 30, 40, 50, 60):
        ks = list(range(N - M + 1, N + 1))
        xs = [mp.mpf(1) / k for k in ks]
        ys = [gamma_est(k) for k in ks]
        gamma_vals[M] = neville_to_zero(xs, ys)
    print("gamma extrapolants (window M -> value):")
    for M, v in gamma_vals.items():
        print(f"   M={M:3d}:  gamma* = {mp.nstr(v, 20)}")
    gstar = gamma_vals[40]
    spread_g = max(abs(gamma_vals[M] - gstar) for M in gamma_vals)
    print(f"gamma*  = {mp.nstr(gstar, 16)}   (window spread ~ {mp.nstr(spread_g,3)})")
    print(f"vs -1/2 : diff = {mp.nstr(gstar + mp.mpf(1)/2, 5)}")

    # Now extract C.  Do it BOTH with the fitted gamma* and with the clean
    # hypothesis gamma = -1/2, so we can see which gives a stable C.
    def C_series(gamma):
        # C_k = -a_k / (Gamma(2k+gamma) (9/8)^k) -> C
        out = {}
        for M in (20, 30, 40, 50, 60):
            ks = list(range(N - M + 1, N + 1))
            xs = [mp.mpf(1) / k for k in ks]
            ys = [-am[k] / (mp.gamma(2 * k + gamma) * nine_eighths**k) for k in ks]
            out[M] = neville_to_zero(xs, ys)
        return out

    for label, gamma in (("gamma=gamma*", gstar), ("gamma=-1/2", mp.mpf(-1)/2)):
        Cvals = C_series(gamma)
        Cstar = Cvals[40]
        spread = max(abs(Cvals[M] - Cstar) for M in Cvals)
        print(f"\nC extrapolants [{label}]:")
        for M, v in Cvals.items():
            print(f"   M={M:3d}:  C* = {mp.nstr(v, 20)}")
        print(f"C*  = {mp.nstr(Cstar, 18)}   (window spread ~ {mp.nstr(spread,3)})")
        # closed-form probes
        inv4sqrtpi = 1 / (4 * mp.sqrt(mp.pi))
        print(f"   C* * 4 sqrt(pi)      = {mp.nstr(Cstar*4*mp.sqrt(mp.pi), 16)}  (=1 if C=1/(4 sqrt pi))")
        print(f"   1/(4 sqrt pi)        = {mp.nstr(inv4sqrtpi, 16)}")
        print(f"   C* - 1/(4 sqrt pi)   = {mp.nstr(Cstar - inv4sqrtpi, 6)}")
        try:
            ident = mp.identify(Cstar, ['pi', 'sqrt(pi)', 'sqrt(2)', 'sqrt(3)'])
            print(f"   mp.identify(C*)      = {ident}")
        except Exception as e:
            print(f"   mp.identify failed: {e}")

    # ---- closed-form identification of C via PSLQ (integer relation) ----
    print("\n# CLOSED-FORM SEARCH for C")
    Cstar = C_series(mp.mpf(-1)/2)[60]
    print(f"C* (gamma=-1/2, M=60) = {mp.nstr(Cstar, 30)}")
    A = 2 * mp.sqrt(2) / 3
    probes = {
        "1/(4 sqrt pi)":       1/(4*mp.sqrt(mp.pi)),
        "C*/sqrt(pi)":         Cstar/mp.sqrt(mp.pi),
        "C* * 2 pi":           Cstar*2*mp.pi,
        "C* * 2 pi / sqrt(A)": Cstar*2*mp.pi/mp.sqrt(A),
        "C* * sqrt(2 pi)":     Cstar*mp.sqrt(2*mp.pi),
        "C* * pi":             Cstar*mp.pi,
    }
    for k, v in probes.items():
        print(f"   {k:22s} = {mp.nstr(v, 18)}")
    # PSLQ: is  C = 2^p1 * 3^p2 * pi^p3  ?  -> relation among logs
    try:
        rel = mp.pslq([mp.log(Cstar), mp.log(2), mp.log(3), mp.log(mp.pi)],
                      maxcoeff=10**6, maxsteps=10**5)
        print(f"   pslq[lnC, ln2, ln3, lnpi] = {rel}")
        if rel and rel[0] != 0:
            n0, n2, n3, npi = rel
            print(f"     => C^{n0} = 2^{-n2} 3^{-n3} pi^{-npi}  "
                  f"(i.e. C = 2^{mp.mpf(-n2)/n0} 3^{mp.mpf(-n3)/n0} pi^{mp.mpf(-npi)/n0})")
    except Exception as e:
        print(f"   pslq(logs) failed: {e}")
    for label, target in (("C", Cstar), ("C*sqrt(pi)", Cstar*mp.sqrt(mp.pi)),
                          ("C*2pi", Cstar*2*mp.pi)):
        try:
            ident = mp.identify(target, ['pi', 'sqrt(pi)', 'sqrt(2)', 'sqrt(3)',
                                          'sqrt(6)', 'log(2)', 'gamma(1/4)'])
            print(f"   identify({label}) = {ident}")
        except Exception as e:
            print(f"   identify({label}) failed: {e}")

    # save extended coeffs (exact strings; floats only where representable)
    def safe_float(x):
        try:
            return float(x)
        except OverflowError:
            return None
    out = {"N": N, "gamma": "-1/2",
           "C_str": mp.nstr(Cstar, 30),
           "a_exact": [f"{x.numerator}/{x.denominator}" for x in a],
           "a_float": [safe_float(x) for x in a]}
    json.dump(out, open(os.path.join(HERE, "_hm_coeffs_exact.json"), "w"))
    print(f"\nsaved {N+1} exact coeffs + C -> _hm_coeffs_exact.json")


if __name__ == "__main__":
    main()
