"""
Milestone 2 -- the decisive cheap test the workflow flagged but did not run:
IS THE MEASURED ~54deg A REAL SIGNAL, OR A 6-COEFFICIENT EXTRACTION BIAS?

The stochastic Borel phase is extracted from only 6 variance coefficients v0..v5, and there is a KNOWN
competing REAL Borel pole (the far-tail instanton s^5/10) at |zeta|~1.5, arg 0, at nearly the same modulus
as the claimed complex pair. With 6 terms a [3/2] or [2/2] Borel-Pade must fit BOTH the complex pair and the
real pole -- so the apparent pair phase can be biased.

TEST: build SYNTHETIC ladders with a KNOWN true structure -- a conjugate pair at phase theta_true (radius
rho_pair) PLUS a real pole (radius rho_real, amplitude R) -- and measure what phase the 6-coefficient
Borel-Pade / Darboux estimators REPORT. Map theta_meas vs theta_true. If a true 45deg pair + the real
instanton pole is systematically MIS-READ as ~54deg, then the 54-vs-45 gap is an artifact and 45deg
(inherited from lambda0) is consistent with the data.

Borel transform b_n = v_n/n!; a singularity of B(t)=sum b_n t^n at t=zeta contributes b_n ~ zeta^{-(n+1)}.
A conjugate pair at zeta=rho e^{+-i theta} with residue C: v_n = 2 Re[ C rho^{-(n+1)} e^{-i(n+1)theta} ] n!.
A real pole at zeta=rho_r, residue R (real): v_n += R rho_r^{-(n+1)} n!.
numpy + scipy.
"""
import numpy as np
from numpy.polynomial import polynomial as P
from math import factorial
from scipy.interpolate import pade

REAL = np.array([0.134, 0.111, 0.100, -0.02, -0.45, -1.1])   # the measured ladder v0..v5

def make_ladder(theta_deg, rho_pair=1.258, C=None, rho_real=None, R=0.0, N=6, phi=0.0):
    """Synthetic v0..v_{N-1} from a conjugate pair (+ optional real pole)."""
    th = np.radians(theta_deg)
    if C is None: C = 1.0
    v = np.zeros(N)
    for n in range(N):
        pair = 2*np.real(C*np.exp(1j*phi) * rho_pair**(-(n+1)) * np.exp(-1j*(n+1)*th)) * factorial(n)
        real = (R * rho_real**(-(n+1)) * factorial(n)) if (rho_real is not None) else 0.0
        v[n] = pair + real
    return v

def borel_pade_phase(v, L, M):
    """Dominant complex pole phase (deg, |arg|) and radius from [L/M] Pade of the Borel transform."""
    b = [v[n]/factorial(n) for n in range(len(v))]
    try:
        p, q = pade(b, L, M)
    except Exception:
        return None
    roots = P.polyroots(q.coef if hasattr(q,'coef') else q)
    roots = [r for r in roots if abs(r) > 1e-9]
    if not roots: return None
    # dominant = smallest |root|; pick the one with nonzero imaginary part if present
    cx = [r for r in roots if abs(r.imag) > 1e-6*abs(r)]
    pick = min(cx, key=abs) if cx else min(roots, key=abs)
    return abs(np.degrees(np.angle(pick))), abs(pick)

def darboux_phase(v, alpha=0.0):
    """Floated-theta Darboux fit v_n ~ 2C rho^{-(n+1)} Gamma(n+1+alpha) cos((n+1)theta - phi). Grid search."""
    from scipy.special import gamma as G
    n = np.arange(len(v))
    best = None
    for rho in np.linspace(0.9, 2.2, 80):
        for th in np.radians(np.linspace(20, 80, 120)):
            A = np.vstack([rho**(-(n+1))*G(n+1+alpha)*np.cos((n+1)*th),
                           rho**(-(n+1))*G(n+1+alpha)*np.sin((n+1)*th)]).T
            coef,res,*_ = np.linalg.lstsq(A, v, rcond=None)
            pred = A@coef
            err = np.sum((pred-v)**2)
            if best is None or err < best[0]:
                best = (err, np.degrees(th), rho)
    return best[1], best[2]

if __name__ == "__main__":
    np.set_printoptions(precision=4, suppress=True)
    print("="*78)
    print("BIAS CALIBRATION: what phase does the 6-coeff extraction REPORT for a KNOWN true phase?")
    print("="*78)

    print("\n[1] PURE conjugate pair, no real pole -- is the estimator unbiased on 6 terms?")
    print(f"  {'theta_true':>10} | {'[3/2] meas':>12} {'[2/2] meas':>12} {'[1/3] meas':>12} {'Darboux':>10}")
    for tt in (35, 40, 45, 50, 54, 60):
        v = make_ladder(tt, rho_pair=1.258)
        r32 = borel_pade_phase(v,3,2); r22 = borel_pade_phase(v,2,2); r13 = borel_pade_phase(v,1,3)
        dth,_ = darboux_phase(v)
        f = lambda r: f"{r[0]:6.1f}@{r[1]:.2f}" if r else "  --  "
        print(f"  {tt:>10} | {f(r32):>12} {f(r22):>12} {f(r13):>12} {dth:>9.1f}")

    print("\n[2] TRUE pair at 45deg + a REAL instanton pole (rho_real=1.5) of increasing strength R:")
    print("    (does the real pole push the APPARENT pair phase up toward 54?)")
    print(f"  {'R (real amp)':>12} | {'[3/2] meas':>12} {'[2/2] meas':>12} {'Darboux':>10}")
    for R in (0.0, 0.3, 0.6, 1.0, 1.6, 2.5):
        v = make_ladder(45, rho_pair=1.258, rho_real=1.5, R=R)
        r32 = borel_pade_phase(v,3,2); r22 = borel_pade_phase(v,2,2); dth,_ = darboux_phase(v)
        f = lambda r: f"{r[0]:6.1f}@{r[1]:.2f}" if r else "  --  "
        print(f"  {R:>12.2f} | {f(r32):>12} {f(r22):>12} {dth:>9.1f}")

    print("\n[3] Can a TRUE-45deg model be tuned to MATCH the real ladder v0..v5, and then mis-read as ~54?")
    # least-squares fit C(complex) and R(real) so pair@45 + real-pole@1.5 best matches REAL
    th = np.radians(45.0); rho=1.258; rr=1.5; n=np.arange(6)
    # basis: Re-part, Im-part of pair, and real pole
    B = np.vstack([2*rho**(-(n+1))*np.cos((n+1)*th)*np.array([factorial(k) for k in n]),
                   2*rho**(-(n+1))*np.sin((n+1)*th)*np.array([factorial(k) for k in n]),
                   rr**(-(n+1))*np.array([factorial(k) for k in n])]).T
    coef,res,*_ = np.linalg.lstsq(B, REAL, rcond=None)
    fit = B@coef
    print(f"  best TRUE-45 (pair@45 + real@1.5) fit to REAL ladder: coef(Cre,Cim,R)={coef}")
    print(f"  fitted v0..v5 = {fit}")
    print(f"  actual v0..v5 = {REAL}")
    print(f"  residual RMS  = {np.sqrt(np.mean((fit-REAL)**2)):.4f}   (max |v|~1.1, so relative ~{np.sqrt(np.mean((fit-REAL)**2))/1.1:.2f})")
    for L,M in [(3,2),(2,2),(1,3)]:
        r = borel_pade_phase(fit, L, M)
        print(f"    -> this TRUE-45 fit, read by [{L}/{M}] Borel-Pade, REPORTS phase: {r[0]:.1f}deg @ |zeta|={r[1]:.2f}" if r else f"    [{L}/{M}] failed")
    dth,drho = darboux_phase(fit)
    print(f"    -> Darboux reports: {dth:.1f}deg @ |zeta|={drho:.2f}")
    print(f"\n  For reference, the REAL ladder itself reads:")
    for L,M in [(3,2),(2,2),(1,3)]:
        r = borel_pade_phase(REAL, L, M)
        print(f"    [{L}/{M}]: {r[0]:.1f}deg @ |zeta|={r[1]:.2f}" if r else f"    [{L}/{M}] failed")

    print("\n" + "="*78)
    print("[4] IDENTIFIABILITY: fit (pair@theta_true + real pole@rho_r) to the ASYMPTOTIC tail v1..v5;")
    print("    is the fit residual FLAT in theta_true? (flat => 6 coeffs cannot discriminate the phase)")
    print("="*78)
    from math import factorial as fac
    ns = np.arange(1,6)          # tail v1..v5 (the large-order/asymptotic part; v0 is leading, not asymptotic)
    tgt = REAL[1:6]
    print(f"  {'theta_true':>10} {'best rho_r':>10} {'resid RMS':>10}  {'(rel to |v5|=1.1)':>18}")
    curve=[]
    for tt in range(35, 66, 3):
        th=np.radians(tt); rho=1.258
        best=None
        for rr in np.linspace(1.0, 2.2, 60):
            Bm = np.vstack([2*rho**(-(ns+1))*np.cos((ns+1)*th)*np.array([fac(k) for k in ns]),
                            2*rho**(-(ns+1))*np.sin((ns+1)*th)*np.array([fac(k) for k in ns]),
                            rr**(-(ns+1))*np.array([fac(k) for k in ns])]).T
            coef,*_=np.linalg.lstsq(Bm,tgt,rcond=None); r=np.sqrt(np.mean((Bm@coef-tgt)**2))
            if best is None or r<best[0]: best=(r,rr)
        curve.append((tt,best[0]))
        print(f"  {tt:>10} {best[1]:>10.2f} {best[0]:>10.4f}  {best[0]/1.1:>17.3f}")
    rmin=min(c[1] for c in curve); rmax=max(c[1] for c in curve)
    print(f"\n  residual ranges {rmin:.4f} to {rmax:.4f} over theta_true in [35,63] deg.")
    print(f"  => spread/min = {rmax/max(rmin,1e-9):.2f}x. If ~1x (flat), the phase is UNIDENTIFIABLE from 6 coeffs.")
    print(f"  best theta_true = {min(curve,key=lambda c:c[1])[0]} deg (but note how shallow the minimum is).")
