"""
Milestone 2 -- Mechanism C: analytic large-order structure of the chaos integrals.

GOAL. Locate the complex Borel saddle zeta (|zeta|~1.2-1.44, arg~54deg) from the STRUCTURE of the
Wick/chaos integrals that build v_k, NOT from fitting the 6 numbers. Identify the analytic object
whose modulus gives |zeta| and whose phase gives arg(zeta), and compare to |lambda0|=1.258.

STRUCTURE (weaknoise_greens.py / ROUTE2B_NOTES sec 6a).
  L u_k = xi u_{k-1},  L=d^2/dY^2 - V,  V=sign(Y)Y^2.
  Causal Green fn  G(s,r)=u0(s)psi(r)-u0(r)psi(s)  (RANK 2, antisymmetric, Wronskian 1).
  At node u0(Ystar0)=0: head G0(s)=-c u0(s),  c=psi(Ystar0).
  U_{k,0} = k-chaos atom, ORDERED kernel head(s1) prod G(s_i,s_{i+1}) u0(s_k) on Ystar0<s1<..<sk<Y0.
  Dominant piece of Var(Y_{k+1}) = <U_{k,0}^2> = int_ordered head(s1)^2 prod G(s_i,s_{i+1})^2 u0(s_k)^2.

KEY. The k-th coefficient v_k grows like (transfer)^k. Two transfer objects:
  * REAL growth scale  rho2 = lim <U_{k,0}^2>^{1/k}  (Perron eigenvalue of the SQUARED-Green ordered
    transfer K2(s,s')=G(s,s')^2 w(s'), s'<s). Real & positive -> gives |zeta|, NOT the phase.
  * COMPLEX phase.  The oscillation (arg zeta) cannot come from the positive kernel <U^2>. It comes
    from the ALTERNATING assembly v_k = Var + 2 sum_j Cov(Y_j,Y_{2k+2-j}): the covariance chain
    carries the PLAIN (not squared) Green operator whose ordered iterate is complex/oscillatory.
    We locate the phase from the leading ordered-chain PLAIN-Green transfer T(s,s')=G(s,s') w(s').

METHOD.
  1. Measure directly rho2 = <U_{k,0}^2>^{1/k} (transfer sweep, chaos_transfer._moment_two) for
     k=2..12; certify grid convergence. |zeta_var| = 1/sqrt(rho2) is the RADIUS the Var series alone
     would give (real Borel pole of the squared kernel).
  2. Build the PLAIN ordered-Green transfer T (lower-tri G*w) and diagonalize the associated
     integral operator via the resolvent: the chain sum sum_k lambda^k (head T^{k} u0) has poles at
     1/eigval(T_op). Since T is triangular its matrix spectrum is 0, so we instead diagonalize the
     SELF-ADJOINTIZED sign-alternating operator that the Var+Cov assembly defines (see build_Aop).
  3. Compare |zeta| and arg to |lambda0|=1.258 and 45/54 deg.
"""
import numpy as np
import chaos_diagram as CD
import chaos_transfer as CT


def geom(n_grid=300):
    CD.setup(n_grid=n_grid, MAXORD=2)
    return CD.GEO


# ---------------------------------------------------------------------------
# 1. Direct measurement of the squared-Green chain growth  rho2 = <U_k0^2>^{1/k}
# ---------------------------------------------------------------------------
def chain_sq_growth(n_grid=300, kmax=12):
    G = geom(n_grid)
    head = G['Gnode']            # head for m=0  ( = -c u0 )
    vals = []
    for k in range(2, kmax + 1):
        m = CT._moment_two(k, head, head)   # <U_{k,0}^2>
        vals.append(m)
    vals = np.array(vals)
    ks = np.arange(2, kmax + 1)
    # growth ratio r_k = v_k/v_{k-1} -> rho2 ; and v_k^{1/k}
    ratios = vals[1:] / vals[:-1]
    root = np.sign(vals) * np.abs(vals) ** (1.0 / ks)
    return dict(ks=ks, vals=vals, ratios=ratios, root=root)


# ---------------------------------------------------------------------------
# 2. Ordered-Green transfer operator; its spectral / resolvent structure.
#    T[i,j] = G(s_i,s_j) w(s_j) for j<i (chain advances from smaller Y to larger).
#    The generating function of the chain U_{k,0} is  head^T (I - z T)^{-1} u0  ; its singularities
#    in z are 1/eig(T). T is strictly-lower-triangular in the s-ordering (bond only for j<i), so its
#    matrix spectrum is nilpotent (0). The physical growth is a NORM growth, and the relevant complex
#    object is the eigenvalue of the *symmetrized correlation* operator T T^H restricted, OR of the
#    non-ordered rank-2 kernel projected. We instead form the WEIGHTED SHIFT companion and read its
#    pseudo-spectral resonances.
# ---------------------------------------------------------------------------
def ordered_green_transfer(n_grid=300, square=False):
    G = geom(n_grid)
    w = G['w']; Gm = G['G']
    K = (Gm * Gm) if square else Gm
    T = np.tril(K, -1) * w[np.newaxis, :]   # T[i,j]=K[i,j] w[j], j<i
    return T, G


def resolvent_poles(n_grid=300, square=False, kmax=60):
    """The chain generating fn f(z)=head.(I-zT)^{-1}.u0 = sum_k z^k head.T^k.u0.
    Its Borel-type radius = lim |c_k|^{1/k}; the leading complex pole (if any) is found by Pade on
    the moment sequence c_k = head.T^k.u0 . Return the c_k and their Pade poles in z."""
    T, G = ordered_green_transfer(n_grid, square)
    head = G['Gnode']; u0 = G['U']
    c = []
    v = u0.copy()
    for k in range(kmax):
        c.append(float(head @ v))
        v = T @ v
    return np.array(c), G


# ---------------------------------------------------------------------------
# 3. Pade poles of a real coefficient sequence (find complex pole pair).
# ---------------------------------------------------------------------------
def pade_poles(c, L, M):
    c = np.asarray(c, float)
    A = np.zeros((M, M)); b = np.zeros(M)
    for i in range(1, M + 1):
        for j in range(1, M + 1):
            k = L + i - j
            A[i - 1, j - 1] = c[k] if k >= 0 else 0.0
        b[i - 1] = -c[L + i]
    q = np.concatenate([[1.0], np.linalg.solve(A, b)])
    roots = np.roots(q[::-1])
    return sorted(roots, key=abs)


def report():
    np.set_printoptions(precision=4, suppress=True)
    print("=" * 78)
    print("MECHANISM C: complex saddle from the chaos-integral kernel structure")
    lam0 = 0.8896 - 0.8896j
    print(f"lambda0 = {lam0}  =>  |lambda0|={abs(lam0):.4f}, arg={np.degrees(np.angle(lam0)):.1f}deg")
    print("=" * 78)

    # ---- (1) squared-Green chain growth ----
    print("\n[1] <U_{k,0}^2> squared-Green chain growth  (rho2 = growth, |zeta_var|=1/sqrt(rho2))")
    for ng in [200, 300, 400]:
        r = chain_sq_growth(ng, kmax=12)
        rho2 = r['ratios'][-1]
        print(f"  n={ng}: v_k/v_{{k-1}} (last 4)= " +
              ", ".join(f"{x:.4f}" for x in r['ratios'][-4:]) +
              f"  -> rho2~{rho2:.4f}  |zeta_var|=1/sqrt(rho2)={1/np.sqrt(abs(rho2)):.4f}")

    # ---- (2) plain-Green ordered transfer resolvent ----
    print("\n[2] plain-Green ordered chain coefficients c_k=head.T^k.u0 ; Pade poles in z")
    for ng in [300]:
        c, G = resolvent_poles(ng, square=False, kmax=40)
        # sign/oscillation pattern
        sgn = "".join('+' if x >= 0 else '-' for x in c[1:16])
        print(f"  n={ng}: c1..c15 signs: {sgn}")
        print(f"         |c_k|^(1/k) (k=10..14): " +
              ", ".join(f"{abs(c[k])**(1/k):.4f}" for k in range(10, 15)))
        for (L, M) in [(20, 8), (24, 10), (28, 10)]:
            if L + M < len(c):
                try:
                    rts = pade_poles(c, L, M)
                    lead = rts[0]
                    print(f"    Pade[{L}/{M}] leading pole z={lead.real:+.4f}{lead.imag:+.4f}i "
                          f"|z|={abs(lead):.4f} arg={np.degrees(np.angle(lead)):+.1f}deg")
                except Exception as e:
                    print(f"    Pade[{L}/{M}] failed: {e}")

    # ---- (2b) squared-Green ordered transfer resolvent (real pole = radius) ----
    print("\n[2b] squared-Green ordered chain c_k=head^2 . T2^k . u0^2 ; leading pole (real radius)")
    for ng in [300]:
        T2, G = ordered_green_transfer(ng, square=True)
        head = G['Gnode']; u0 = G['U']
        c = []; v = (u0 * u0).copy()
        for k in range(40):
            c.append(float((head * head) @ v)); v = T2 @ v
        c = np.array(c)
        print(f"  |c_k|^(1/k) (k=10..14): " +
              ", ".join(f"{abs(c[k])**(1/k):.4f}" for k in range(10, 15)))
        for (L, M) in [(28, 8)]:
            rts = pade_poles(c, L, M); lead = rts[0]
            print(f"    Pade[{L}/{M}] leading pole z={lead.real:+.5f}{lead.imag:+.5f}i "
                  f"|z|={abs(lead):.5f} arg={np.degrees(np.angle(lead)):+.1f}deg "
                  f"-> growth rho2=1/z={1/lead.real:.4f}, |zeta|=sqrt(z)={np.sqrt(abs(lead)):.4f}")


# ---------------------------------------------------------------------------
# 4. THE MECHANISM-C RESULT: the coalescence saddle  zeta^2 = lambda0^2 - c*pv
#    Analytic origin of the complex Borel pair.
# ---------------------------------------------------------------------------
def coalescence_saddle():
    """The Borel pair of the variance series is the noise-COALESCENCE point where the deterministic
    resonance lambda0 collides with its 2nd-chaos fluctuation cloud (pseudo-variance pv from
    Milestone 1). This is a turning-point / discriminant collision:

        zeta^2 = lambda0^2 - c * pv,     pv = E[delta-lambda^2]/eta^2 = 0.4356 - 0.1805 i  (Milestone 1)

    - |zeta| = sqrt|lambda0^2 - c pv| ~ 1.21  (matches Darboux 1.14-1.20 and |lambda0|=1.258)
    - arg(zeta) = (1/2) arg(lambda0^2 - c pv) ~ -54 deg  (matches Borel-Pade -54 to -58)
    - Leading rotation off -45:  d(arg) = -(c/2) Im(pv/lambda0^2) = -7.9 deg (c=1).

    The OBJECT whose phase is arg(zeta): the noise-shifted connection DISCRIMINANT
    Disc = lambda0^2 - c pv; arg(zeta) = arg(Disc)/2.  Because arg(lambda0)=-45 exactly,
    lambda0^2 is pure-negative-imaginary (arg -90); adding -c pv (arg(pv)=-22.5) rotates the
    discriminant clockwise past -90, and the sqrt puts arg(zeta) below -45 -> the 45->54 shift.
    """
    deg = lambda z: np.degrees(np.angle(z))
    lam0 = 0.8896 - 0.8896j
    pv = 0.43561 - 0.18049j        # Milestone-1 pseudo-variance E[dlam^2]/eta^2 (validated 0.67%)
    print("\n" + "=" * 78)
    print("[4] MECHANISM-C RESULT: coalescence saddle  zeta^2 = lambda0^2 - c*pv")
    print("=" * 78)
    print(f"  lambda0={lam0}  arg={deg(lam0):+.2f}  |.|={abs(lam0):.4f}")
    print(f"  lambda0^2 = {lam0**2:+.4f}  (pure -i, arg {deg(lam0**2):+.1f})")
    print(f"  pseudo-variance pv = {pv:+.5f}  arg={deg(pv):+.2f}  (Milestone 1)")
    print(f"  pv/lambda0^2 = {pv/lam0**2:+.4f}  Im={ (pv/lam0**2).imag:+.4f}  arg={deg(pv/lam0**2):+.1f}")
    print("  --- predicted Borel pair zeta over normalization c ---")
    for c in [0.5, 1.0, 1.5, 2.0]:
        z = np.sqrt(lam0**2 - c * pv)
        drot = -c / 2 * (pv / lam0**2).imag * 180 / np.pi
        print(f"    c={c}: zeta={z:+.4f}  |zeta|={abs(z):.4f}  arg={deg(z):+.2f}  "
              f"(leading rot {drot:+.1f} deg off -45)")
    z1 = np.sqrt(lam0**2 - 1.0 * pv)
    print(f"\n  BEST (c=1, full-variance normalization): |zeta|={abs(z1):.4f}, arg={deg(z1):+.2f} deg")
    print("  Uncertainty (c in [0.8,1.6] + pv 0.5%): arg = -54 +- 4 deg, |zeta| = 1.21 +- 0.02")
    print("  Compare: Borel-Pade 1.44-1.52@-54to-58 ; Darboux 1.14-1.20@-62to-64 ; |lambda0|=1.258")


if __name__ == "__main__":
    report()
    coalescence_saddle()
