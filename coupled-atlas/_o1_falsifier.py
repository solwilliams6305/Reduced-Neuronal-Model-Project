"""O1 falsifier, pass 1 (foundation / necessary-condition checks).

Two clean, self-contained numerical checks that must hold if the O1 reduction is correct:

(A) PAINLEVE-II BACKBONE -> DUMAZ-VIRAG TAIL at beta=2.
    Solve Hastings-McLeod PII  q'' = 2q^3 + s q,  q(s)~Ai(s) (s->+inf).
    Build the TW_2 tail g(s) = int_s^inf (x-s) q(x)^2 dx  (= -log F_2 for the GUE edge).
    Check its right tail:  -log g(s) = (4/3) s^{3/2} + (3/2) log s + const,
    i.e. rate 4/3 = (2/3)*beta at beta=2  AND prefactor exponent -3/2 = -3 beta/4 at beta=2.
    If the PII backbone did NOT reproduce the Dumaz-Virag exponent+prefactor, O1's
    "the tail dresses PII" premise is falsified.

(B) WEAK-NOISE VARIANCE (first-order perturbation theory) for the stochastic Airy operator.
    Ground state psi0(x)=Ai(x - a1) of -d^2+x on R_+ Dirichlet (a1 = 2.3381, first Airy zero).
    First-order eigenvalue shift dLambda = (2/sqrt(beta)) int psi0^2 b'  =>  Var ~ s2/beta,
    s2 = 4 int psi0^4 / (int psi0^2)^2.  Reports s2 (the leading TW_beta variance coefficient
    the CLDS Painleve-II prediction must reproduce).
"""
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.special import airy

# ---------------- (A) Hastings-McLeod PII and the TW_2 tail ----------------
def solve_hm(s_max=8.0, s_min=1.5, rtol=1e-11, atol=1e-13):
    Ai_max, Aip_max, _, _ = airy(s_max)
    def rhs(s, y):
        q, qp = y
        return [qp, 2*q**3 + s*q]
    sol = solve_ivp(rhs, [s_max, s_min], [Ai_max, Aip_max],
                    dense_output=True, rtol=rtol, atol=atol, max_step=0.005)
    return sol

def tw2_tail_g(s, sol, x_hi=12.0):
    # g(s) = int_s^inf (x-s) q(x)^2 dx ; q from HM on [s,x_hi], Ai^2 beyond (q->Ai)
    def q_of(x):
        if x <= sol.t[0]:   # within integrated range (sol.t descending: t[0]=s_max)
            return sol.sol(x)[0]
        Ai_x, _, _, _ = airy(x)
        return Ai_x
    val, _ = quad(lambda x: (x - s) * q_of(x)**2, s, x_hi, limit=200)
    return val

def check_A():
    print("=== (A) PII (Hastings-McLeod) -> TW_2 tail: rate & prefactor ===")
    print("    model: -log g(s) = (4/3) s^{3/2} + (3/2) log s + C")
    print("    i.e. rate 4/3 = (2/3)beta  and tail prefactor s^{-3/2} = s^{-3beta/4}  (beta=2)")
    sol = solve_hm(s_max=10.0, s_min=1.5)
    q8 = sol.sol(8.0)[0]; Ai8,_,_,_ = airy(8.0)
    print(f"  HM sanity: q(8)={q8:.3e} vs Ai(8)={Ai8:.3e}  rel={abs(q8-Ai8)/abs(Ai8):.1e}")
    ss = np.array([3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0])
    G = np.array([-np.log(tw2_tail_g(s, sol, x_hi=16.0)) for s in ss])
    # full 3-param regression
    M = np.vstack([ss**1.5, np.log(ss), np.ones_like(ss)]).T
    coef, *_ = np.linalg.lstsq(M, G, rcond=None)
    # rate-fixed extraction of prefactor: r(s) = -log g - (4/3)s^1.5 = (3/2)log s + C
    r = G - (4.0/3.0)*ss**1.5
    Mp = np.vstack([np.log(ss), np.ones_like(ss)]).T
    bpref, *_ = np.linalg.lstsq(Mp, r, rcond=None)
    print(f"  {'s':>4} {'-log g':>11} {'r=-logg-(4/3)s^1.5':>20}")
    for s, g, rr in zip(ss, G, r):
        print(f"  {s:>4.1f} {g:>11.4f} {rr:>20.4f}")
    print(f"  free 3-param:   rate A={coef[0]:.4f} (exp 1.3333),  prefactor B={coef[1]:.4f} (exp +1.5)")
    print(f"  rate-fixed:     prefactor B={bpref[0]:.4f} (exp +1.5),  C={bpref[1]:.4f}")
    ok_rate = abs(coef[0]-4/3) < 0.03
    ok_pref = abs(bpref[0]-1.5) < 0.2
    print(f"  VERDICT (A): rate {'PASS' if ok_rate else 'FAIL'}, "
          f"prefactor {'PASS' if ok_pref else 'FAIL'}")
    return ok_rate and ok_pref

# ---------------- (B) first-order weak-noise variance ----------------
def check_B():
    print("\n=== (B) weak-noise variance s2 = 4 int psi0^4 (first-order PT) ===")
    a1 = 2.338107410459767   # |first zero of Ai|
    def psi0(x):
        Ai_v, _, _, _ = airy(x - a1)   # Ai(x - a1): =0 at x=0, nodeless ground state on R+
        return Ai_v
    norm2, _ = quad(lambda x: psi0(x)**2, 0, np.inf, limit=200)
    int4, _  = quad(lambda x: psi0(x)**4, 0, np.inf, limit=200)
    s2 = 4*int4/norm2**2
    print(f"  int psi0^2 = {norm2:.6f},  int psi0^4 = {int4:.6f}")
    print(f"  s2 = 4 int psi0^4/(int psi0^2)^2 = {s2:.5f}")
    print(f"  => Var(Lambda_0) ~ {s2:.4f}/beta   (the leading TW_beta variance coeff;")
    print(f"     CLDS's Painleve-II 'Var(a_i) = s2_tilde/beta' must reproduce this number)")
    return s2

if __name__ == "__main__":
    okA = check_A()
    s2 = check_B()
    print("\n=== SUMMARY ===")
    print(f"(A) PII backbone reproduces DV tail (rate+prefactor): {'PASS' if okA else 'FAIL'}")
    print(f"(B) leading weak-noise variance coefficient s2 = {s2:.4f}")
