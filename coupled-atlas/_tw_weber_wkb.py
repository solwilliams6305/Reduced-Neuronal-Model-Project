"""
_tw_weber_wkb.py -- Weber-curve WKB for the TW_beta tail diagonal (Kidwai-Osuga 2204.12431, eq 4.23).

The TW_beta tail lives on the WEBER refined spectral curve (Gaussian-beta-ensemble semicircle):
    eps^2 psi'' + [C + (4 m_inf - x^2)/4] psi = 0,   eps = eps_1 = hbar sqrt(beta), edge x=2 (m_inf=1).
Pure WKB (C=0) of the Riccati W=psi'/psi, W = (1/eps) sum_k eps^k W_k, W_k ~ w_k (x-2)^{(1-3k)/2} at the
edge, gives a CLOSED number recursion for the leading edge coefficients w_k:
    w_0 = -1,   w_k = (1/2)[ sum_{j=1}^{k-1} w_j w_{k-j} + ((4-3k)/2) w_{k-1} ].
Edge double-scaling (a=2+eps^{2/3}s, lam=beta from the rate) maps the k=m+1 leading term to the
TOP-DEGREE-in-X coefficient of R_m:   r_hat_{m} := -4 w_{m+1}/(3m).

PARTIAL RESULT: r_hat_m matches the EXACT top coefficient r_{m,m+1}=[X^{m+1}]R_m of Borot-Nadal at the
NON-deficit orders  R_1: -5/24,  R_3: -1105/4608  (exact), and differs at R_2 (r_hat_2=5/32 vs exact
r_{2,3}=0 -- the degree-DEFICIT / 2-loop cancellation also seen in PROGRAM2_TWBETA_UNIFORM_IN_A_NOTES).

*** WHY THIS DOES NOT GIVE R_4 (structural obstruction, established here) ***
The pure-WKB r_hat_m ALTERNATE in sign (base ->~3/2), whereas the TRUE diagonal 2^n r_{n,n+1}
(= -5/12, 0, -1105/576, ...) is NON-alternating (median-summable). The sign flip and the deficits are
produced by the beta-dependence (the refinement), which the pure WKB drops. And a full edge extraction
of the Weber quantum curve (eq 4.23, PURE Schrodinger form -- NO first-derivative term, beta entering
ONLY through the constant K) shows the K-dependence CANCELS between the edge-scale and the coefficient:
with K=1+q1*X one finds R_1 = -15/8 identically, independent of q1. So eq 4.23 as written cannot inject
the X-dependence of R_m -- the beta-content must enter via the refined-TR<->Gaussian-betaE parameter
DICTIONARY (N-scaling of m_inf/mu/nu), which is not pinned here.

CONCLUSION: this route (like the loop recursion and the trivial refined-Airy curve) does NOT yield a
validated R_4 in this session. What IS delivered: the clean WKB number-recursion for w_k, the exact
top-coefficient match at non-deficit orders, and the structural proof that the constant-K Weber ODE is
insufficient. The blocker is the same everywhere: the exact beta<->spectral-curve parameter map.
"""
import mpmath as mp
from fractions import Fraction as Fr

def w_coeffs(K):
    w = [Fr(-1)]
    for k in range(1, K + 1):
        conv = sum(w[j]*w[k-j] for j in range(1, k))
        w.append(Fr(1, 2)*(conv + Fr(4 - 3*k, 2)*w[k-1]))
    return w

def main():
    K = 40
    w = w_coeffs(K)
    # r_hat_m = -4 w_{m+1}/(3 m)  = leading-WKB top-degree coeff of R_m
    print("# WKB diagonal r_hat_m = -4 w_{m+1}/(3m)  vs exact r_{m,m+1} (top-degree BN coeff)")
    exact = {1: Fr(-5,24), 2: Fr(0), 3: Fr(-1105,4608)}   # exact [X^{m+1}] R_m
    for m in range(1, 7):
        rh = Fr(-4,1)*w[m+1]/(3*m)
        tag = f"  exact r_{{{m},{m+1}}} = {exact[m]}" if m in exact else ""
        print(f"   m={m}: r_hat = {rh}{tag}")

    # large-order growth of r_hat_m -> Borel singularity of the WKB-diagonal series sum 2^m r_hat_m g^{-m}
    mp.mp.dps = 30
    rh = [None] + [mp.mpf(int((Fr(-4,1)*w[m+1]/(3*m)).numerator)) /
                   int((Fr(-4,1)*w[m+1]/(3*m)).denominator) for m in range(1, K)]
    print("\n# growth: r_hat_m / r_hat_{m-1}  (times 1/(2) for the 2^m in g-series) -> 1/(Borel radius)")
    print("#   Borel series sum 2^m r_hat_m g^{-m}: ratio 2 r_hat_m/r_hat_{m-1} -> 1/zeta_c")
    for m in (10, 20, 30, 38):
        ratio = 2*rh[m]/rh[m-1]
        print(f"   m={m:>3}: 2 r_hat_m/r_hat_(m-1) = {mp.nstr(ratio, 10)}   (|.| ~ m for factorial growth)")
    # Domb-Sykes on |r_hat_m|^{1/m} and ratio/m to see factorial base
    print("\n# factorial base: (2 r_hat_m/r_hat_(m-1))/m -> base of Gamma(m) growth")
    for m in (20, 30, 38):
        print(f"   m={m:>3}: ratio/m = {mp.nstr(2*rh[m]/rh[m-1]/m, 10)}")
    print("\n# sign pattern of r_hat_m (median vs alternating):",
          "".join("+" if rh[m] > 0 else "-" for m in range(1, 20)))

if __name__ == "__main__":
    main()
