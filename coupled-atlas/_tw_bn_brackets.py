"""
_tw_bn_brackets.py -- exact Borot-Nadal (arXiv:1111.2761) tail data + cross-validation.

Page 5 of 1111.2761 tabulates the RIGHT-tail brackets for beta=1,2,4 to O(s^{-21/2}):
   1 - TW_beta(s) = [prefactor] * ( 1 + sum_{m=1}^{6} b_m^{(beta)} s^{-3m/2} + O(s^{-21/2}) ).
By Prop 1.1,  log(bracket) = sum_m (beta/2) R_m(2/beta) s^{-3m/2},  so
   R_m(2/beta) = (2/beta) * [s^{-3m/2}] log(bracket_beta).
With X=2/beta:  beta=1 -> X=2,  beta=2 -> X=1,  beta=4 -> X=1/2.

This script:
 (1) extracts R_m(2), R_m(1), R_m(1/2) for m=1..6 from the brackets (EXACT rationals);
 (2) VALIDATES the paper's R_1,R_2,R_3 polynomials against those three points;
 (3) tabulates the new m=4,5,6 data (R_4,R_5,R_6 at X=2,1,1/2) for downstream use.

These three X-points do NOT fix the full R_m (deg<=m+1) for m>=3, but they (a) independently
confirm the whole framework and (b) give exact finite-beta anchors for the uniformity analysis.
"""
import sympy as sp

X = sp.symbols('X')

# --- paper R_1,R_2,R_3 (eqs 1-18..1-20) ---
Rpoly = {
    1: sp.Rational(1, 24)*(-5*X**2 + 9*X - 39),
    2: sp.Rational(5, 64)*(11*X**2 - 19*X + 36),
    3: sp.Rational(1, 4608)*(-1105*X**4 + 3240*X**3 - 23325*X**2 + 34938*X - 41433),
}

# --- page-5 brackets: b_m^{(beta)} as EXACT rationals, m=1..6 (coeff of s^{-3m/2}) ---
# beta=2 (GUE):
b2 = [sp.Rational(-35, 2**3*3), sp.Rational(3745, 2**7*3**2), sp.Rational(-805805, 2**10*3**4),
      sp.Rational(289554265, 2**15*3**5), sp.Rational(-31241084875, 2**18*3**6),
      sp.Rational(23604769513325, 2**22*3**8)]
# beta=1 (GOE):
b1 = [sp.Rational(-41, 2**4*3), sp.Rational(9241, 2**9*3**2), sp.Rational(-5075225, 2**13*3**4),
      sp.Rational(5153008945, 2**19*3**5), sp.Rational(-1674966309205, 2**23*3**6),
      sp.Rational(3985569631633205, 2**28*3**8)]
# beta=4 (GSE):
b4 = [sp.Rational(-143, 2**4*3), sp.Rational(41509, 2**9*3**2), sp.Rational(-20443229, 2**13*3**4),
      sp.Rational(15418569025, 2**19*3**5), sp.Rational(-3330409204735, 2**23*3**6),
      sp.Rational(4908974519795465, 2**28*3**8)]

def log_bracket_coeffs(bcoeffs):
    """coefficients of log(1 + sum b_m u^m) in u (u = s^{-3/2}), m=1..len."""
    u = sp.symbols('u')
    bracket = 1 + sum(bcoeffs[m-1]*u**m for m in range(1, len(bcoeffs)+1))
    ser = sp.series(sp.log(bracket), u, 0, len(bcoeffs)+1).removeO()
    return [sp.nsimplify(ser.coeff(u, m)) for m in range(1, len(bcoeffs)+1)]

# R_m(2/beta) = (2/beta) * [s^{-3m/2}] log(bracket)
Lb1 = log_bracket_coeffs(b1)   # beta=1 -> gives (1/2) R_m(2)
Lb2 = log_bracket_coeffs(b2)   # beta=2 -> gives       R_m(1)
Lb4 = log_bracket_coeffs(b4)   # beta=4 -> gives  2    R_m(1/2)

Rm_at = {}   # Rm_at[m] = (R_m(2), R_m(1), R_m(1/2))
for m in range(1, 7):
    Rm2   = 2*Lb1[m-1]                 # beta=1
    Rm1   =   Lb2[m-1]                 # beta=2
    Rmhalf= sp.Rational(1,2)*Lb4[m-1]  # 2 R_m(1/2) = Lb4  -> R_m(1/2)=Lb4/2 ... wait: 2 R_m(1/2)=Lb4
    Rmhalf= sp.Rational(1,2)*Lb4[m-1]
    Rm_at[m] = (Rm2, Rm1, Rmhalf)

print("="*74)
print("(1) VALIDATE paper R_1,R_2,R_3 against beta=1,2,4 brackets (X=2,1,1/2)")
print("="*74)
ok = True
for m in (1, 2, 3):
    R2p, R1p, Rhp = Rpoly[m].subs(X, 2), Rpoly[m].subs(X, 1), Rpoly[m].subs(X, sp.Rational(1,2))
    R2b, R1b, Rhb = Rm_at[m]
    c2 = sp.simplify(R2p - R2b) == 0
    c1 = sp.simplify(R1p - R1b) == 0
    ch = sp.simplify(Rhp - Rhb) == 0
    ok = ok and c2 and c1 and ch
    print(f"  R_{m}: X=2 {R2p}=={R2b}?{c2}   X=1 {R1p}=={R1b}?{c1}   X=1/2 {Rhp}=={Rhb}?{ch}")
print(f"\n  ALL R_1,R_2,R_3 match the brackets at X=2,1,1/2:  {ok}")

print("\n" + "="*74)
print("(2) NEW exact data from the brackets: R_m at X=2,1,1/2 for m=4,5,6")
print("="*74)
for m in (4, 5, 6):
    R2b, R1b, Rhb = Rm_at[m]
    print(f"  R_{m}(2)   = {R2b}")
    print(f"  R_{m}(1)   = {R1b}")
    print(f"  R_{m}(1/2) = {Rhb}")

print("\n" + "="*74)
print("(3) NOTE on completeness")
print("="*74)
print("  R_m has degree <= m+1 in X (m+2 unknown coeffs). Three X-points fix R_1,R_2 (deg 2)")
print("  and check R_3 (deg 4, over-constrained here only with the paper value).")
print("  For m>=4 the full polynomial needs the BN topological recursion (2-19/2-20); these")
print("  brackets give finite-beta anchors (X=2,1,1/2) but not the diagonal r_{m,m+1} (X->inf).")
