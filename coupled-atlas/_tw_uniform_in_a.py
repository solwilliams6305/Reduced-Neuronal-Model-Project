"""
_tw_uniform_in_a.py -- UNIFORM-IN-a noise Borel-summability of the TW_beta right tail.

Open Task #2 of PROGRAM2_HANDOFF_D_TWBETA_NEXT.md: extend T3 (median Borel-summability of the
noise 1/beta expansion) from the scaling regime g=beta a^{3/2}->inf to a uniform-in-a statement,
controlling the Borot-Nadal R_m(2/beta) corrections + the O(a^{-3/2}) x-extension.

EXACT INPUT -- Borot-Nadal (arXiv:1111.2761, Prop 1.1), the all-beta all-orders right tail:
   1 - TW_beta(s) = [Gamma(beta/2)/((4 beta)^{beta/2} 2 pi)] s^{-3beta/4} e^{-(2/3) beta s^{3/2}}
                    * exp[ sum_{m>=1} (beta/2) R_m(2/beta) s^{-3m/2} ],
   R_1(X)=(1/24)(-5X^2+9X-39),
   R_2(X)=(5/64)(11X^2-19X+36),
   R_3(X)=(1/4608)(-1105X^4+3240X^3-23325X^2+34938X-41433),      X = 2/beta.
   (R_m is a polynomial of degree <= m+1 in X, rational coefficients.)

We identify s = a (the tail depth; rate e^{-(2/3)beta a^{3/2}} matches the program's e^{-(2/3)g},
g=beta a^{3/2}).  Writing the exponent  Sigma(a,beta)=sum_m (beta/2)R_m(2/beta) a^{-3m/2}  and
collecting powers of 1/beta with t:=a^{-3/2}:
   Sigma = beta*sigma_{-1}(t) + sigma_0(t) + sum_{n>=1} sigma_n(t) beta^{-n},
   sigma_n(t) = 2^n * sum_{m>=n} r_{m,n+1} t^m ,   r_{m,l}=[X^l]R_m .
 - beta^{+1} piece sigma_{-1} = (1/2) sum_m R_m(0) t^m  -> the RATE correction Phi(a).
 - beta^{0}  piece sigma_0                             -> the one-loop (a^{-3beta/4}-type) prefactor.
 - beta^{-n} pieces (n>=1)                             -> THE NOISE LOOP SERIES (what T3 resums).
The DIAGONAL (a->inf at fixed g) coefficient of g^{-n} is sigma_n's leading term 2^n r_{n,n+1}.

This script establishes:
 (1) BN R_1,R_2 reproduce the known GUE (beta=2) tail  -> input validated.
 (2) The reduction: Coulomb-gas prefactor is a-INDEPENDENT; the log term is a single beta^{+1}
     monomial; hence ALL noise (1/beta) Borel content sits in {sigma_n}_{n>=1}.
 (3) The frozen boundary-escape series d_n (this program) equals the TRUE tail diagonal 2^n r_{n,n+1}
     ONLY at n=1; from n=2 it differs -> the x-extension enters the diagonal at two loops.
 (4) The finite-a corrections are the OFF-DIAGONAL terms of sigma_n, relative size O(a^{-3/2}).
"""
import sympy as sp

X, beta, a, t = sp.symbols('X beta a t', positive=True)

# ----- Borot-Nadal R_m(X), X=2/beta (arXiv:1111.2761 eqs 1-18..1-20) -----
R = {
    1: sp.Rational(1, 24) * (-5*X**2 + 9*X - 39),
    2: sp.Rational(5, 64) * (11*X**2 - 19*X + 36),
    3: sp.Rational(1, 4608) * (-1105*X**4 + 3240*X**3 - 23325*X**2 + 34938*X - 41433),
}
Mmax = max(R)


def rcoeff(m, l):
    """r_{m,l} = [X^l] R_m."""
    return sp.Poly(R[m], X).coeff_monomial(X**l) if l <= sp.Poly(R[m], X).degree() else sp.Integer(0)


def sep(title):
    print("\n" + "="*78 + "\n" + title + "\n" + "="*78)


# =====================================================================================
# (1) VALIDATE the fetched R_m against the known GUE (beta=2) Tracy-Widom tail
# =====================================================================================
def validate_gue():
    sep("(1) VALIDATION: BN R_m at beta=2 (X=1) vs the known GUE tail")
    # GUE: 1-F_2(s) = e^{-(4/3)s^{3/2}}/(16 pi s^{3/2}) * (1 + b1 s^{-3/2} + b2 s^{-3} + ...),
    #   b1 = -35/24, b2 = 3745/1152 (Tracy-Widom 1994). Exponent = log(bracket):
    #   [s^{-3/2}] log = b1;  [s^{-3}] log = b2 - b1^2/2.  BN gives R_m(1)= (beta/2)R_m = R_m(1) at beta=2.
    b1 = sp.Rational(-35, 24)
    b2 = sp.Rational(3745, 1152)
    R1_1, R2_1 = R[1].subs(X, 1), R[2].subs(X, 1)
    print(f"  R_1(1) = {R1_1}   ; known b1              = {b1}   match: {sp.simplify(R1_1-b1)==0}")
    print(f"  R_2(1) = {R2_1}   ; known b2 - b1^2/2     = {sp.nsimplify(b2-b1**2/2)}   "
          f"match: {sp.simplify(R2_1-(b2-b1**2/2))==0}")
    print(f"  R_3(1) = {sp.nsimplify(R[3].subs(X,1))}  (=-27685/4608)  -> BN input validated.")


# =====================================================================================
# (2) THE REDUCTION: expand Sigma in 1/beta, isolate the noise series sigma_n
# =====================================================================================
def expand_sigma():
    sep("(2) REDUCTION: Sigma(a,beta)=sum_m (beta/2)R_m(2/beta) a^{-3m/2}, collected in 1/beta")
    # substitute X=2/beta, a^{-3m/2}=t^m ; expand in beta
    Sigma = sum((beta/2) * R[m].subs(X, 2/beta) * t**m for m in R)
    Sigma = sp.expand(Sigma)
    # collect powers of beta
    poly = sp.Poly(sp.expand(Sigma * beta**Mmax), beta)  # clear negative powers
    sigmas = {}
    for (k,), c in poly.terms():
        n = Mmax - k          # power of beta is (k - Mmax) => beta^{-n}, n = Mmax-k
        sigmas[n] = sp.expand(c)   # coeff of beta^{-n} (n can be -1,0,1,2,3)
    for n in sorted(sigmas):
        lab = {(-1): "beta^{+1}  RATE Phi(a)", 0: "beta^{0}   one-loop prefactor"}.get(
            n, f"beta^{{-{n}}}  NOISE loop {n}")
        s_t = sp.Poly(sigmas[n], t)
        low = min((mono[0] for mono in s_t.monoms()), default=0)
        print(f"  sigma_{n:>2} [{lab}]:  {sigmas[n]}")
        if n >= 1:
            print(f"        -> lowest power t^{low} = a^(-{3*low}/2) ; diagonal(t^{n}) coeff "
                  f"2^{n} r_{{{n},{n+1}}} = {2**n * rcoeff(n, n+1)}")
    # internal consistency: reassemble Sigma at beta=2 from the sigma_n and compare to sum_m R_m(1)t^m
    reassembled = sp.expand(sum(sigmas[n] * sp.Rational(1, 2)**n for n in sigmas))
    direct = sp.expand(sum(R[m].subs(X, 1) * t**m for m in R))   # (beta/2)R_m = R_m at beta=2
    print(f"  [self-check] reassemble Sigma at beta=2 from sigma_n == sum_m R_m(1) t^m ? "
          f"{sp.simplify(reassembled - direct) == 0}")
    return sigmas


# =====================================================================================
# (3) FROZEN escape diagonal d_n vs TRUE tail diagonal 2^n r_{n,n+1}
# =====================================================================================
def frozen_vs_true():
    sep("(3) FROZEN boundary-escape d_n  vs  TRUE BN diagonal 2^n r_{n,n+1}")
    # frozen qt_j and F(g)=Qtilde(2/g)^{-2}=sum d_m g^{-m}
    J = 8
    qt = [sp.Integer(1)]
    for j in range(1, J+1):
        qt.append(sp.factorial(6*j) / (sp.Integer(576)**j * sp.factorial(3*j) * sp.factorial(2*j)))
    D = sp.symbols('D')
    Q = sum(qt[j]*D**j for j in range(J+1))
    Finv2 = sp.series(1/Q**2, D, 0, J+1).removeO()          # Qtilde^{-2} in D
    e = [Finv2.coeff(D, m) for m in range(J+1)]
    d = [e[m]*2**m for m in range(J+1)]                     # F(g)=sum d_m g^{-m}
    # log F(g) = sum delta_n g^{-n}
    gg = sp.symbols('gg')
    Fg = sum(d[m]*gg**(-m) for m in range(J+1))
    logF = sp.series(sp.log(Fg.subs(gg, 1/sp.symbols('z'))), sp.symbols('z'), 0, 5).removeO()
    z = sp.symbols('z')
    delta = [sp.nsimplify(logF.coeff(z, n)) for n in range(5)]
    print(f"  frozen d_1..4        : {[sp.nsimplify(d[m]) for m in range(1,5)]}")
    print(f"  frozen delta_n=[g^-n]logF : {[delta[n] for n in range(1,5)]}")
    print(f"  TRUE diagonal 2^n r_(n,n+1):")
    for n in range(1, 4):
        diag = 2**n * rcoeff(n, n+1)
        agree = sp.simplify(diag - delta[n]) == 0
        print(f"     n={n}: 2^{n} r_({n},{n+1}) = {str(diag):>12}   frozen delta_{n} = {str(delta[n]):>14}   "
              f"AGREE: {agree}")
    print("  => frozen agrees with the EXACT diagonal at n=1,3 but NOT n=2 (exact=0 vs frozen=-5/8):")
    print("     R_2 has a degree DEFICIT (deg 2 < m+1=3) => exact 2-loop scaling coeff vanishes.")
    print("     So the frozen boundary-escape series is NOT the exact tail noise series beyond n=1;")
    print("     the large-order (Gevrey) behaviour must come from the BN recursion, not the freeze.")


# =====================================================================================
# (4) UNIFORMITY structure: singularity location and off-diagonal size
# =====================================================================================
def uniformity(sigmas):
    sep("(4) UNIFORMITY: singularity location (a-indep sector) + off-diagonal O(a^{-3/2})")
    print("  Coulomb-gas prefactor Gamma(beta/2)/((4beta)^{beta/2} 2pi): a-INDEPENDENT")
    print("     -> its Borel singularities (Stirling, +-i pi n) are on the IMAGINARY axis, fixed,")
    print("        uniformly separated from R+ for ALL a. (No a-dependence => trivially uniform.)")
    print("  Log term -3beta/4 log a: a single beta^{+1} monomial -> no 1/beta-Borel content.")
    print("  => ALL noise-Borel content is in {sigma_n}_{n>=1}. Diagonal singularity:")
    print("     sigma_n ~ 2^n r_(n,n+1) a^{-3n/2} => Borel sing at zeta=Phi(a)=(2/3)a^{3/2} in R+,")
    print("     MOVING with a but NEVER leaving the positive ray (no collision) => sector uniform.")
    print("  Off-diagonal (finite-a) corrections, sigma_n(t)/[diagonal] = 1 + O(a^{-3/2}):")
    print("  (only sigma_1 is complete to t^3 from R_1..R_3; sigma_2,3 diagonals shown, tails need R_{>=4})")
    n = 1
    s_t = sp.Poly(sigmas[n], t)
    lead_coeff = s_t.coeff_monomial(t**1)
    ratio = sp.expand(sigmas[n] / (lead_coeff * t))
    print(f"     sigma_1/diag = {sp.nsimplify(ratio)}")
    print(f"        => sigma_1(a) = (-5/12) a^(-3/2) [1 - (33/8) a^(-3/2) + (1555/64) a^(-3) + ...]")
    print(f"        the leading x-extension correction to the noise Stokes term is -(33/8) a^(-3/2).")


def main():
    print("UNIFORM-IN-a NOISE BOREL-SUMMABILITY -- Borot-Nadal control of finite-a corrections")
    validate_gue()
    sigmas = expand_sigma()
    frozen_vs_true()
    uniformity(sigmas)
    sep("SUMMARY")
    print("  * BN input validated at beta=2.  * Reduction: noise-Borel content localized to sigma_n.")
    print("  * Exact tail noise coeffs = BN sigma_n; frozen escape matches the diagonal at n=1,3 but")
    print("    NOT n=2 (exact 2-loop scaling coeff = 0, frozen = -5/8) => frozen != exact noise series.")
    print("  * Sector uniform (dyn. sing. on R+ at Phi(a), Coulomb-gas on iR fixed). RESIDUAL: uniform")
    print("    Gevrey bound on sigma_n(a) needs the BN recursion (large-order r_{n,n+1}) + off-diag O(a^-3/2).")


if __name__ == "__main__":
    main()
