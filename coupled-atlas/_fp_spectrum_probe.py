"""
_fp_spectrum_probe.py — frontier move #1: complex spectrum of the noise-dressed first-explosion generator.

Deep-research Finding 5 (Dao Duc-Schuss-Holcman): non-self-adjoint Fokker-Planck first-passage generators
carry COMPLEX eigenvalues computable by semiclassical asymptotics -> oscillatory survival. In resurgence the
Borel singularities of a perturbative series are governed by the operator's complex resonances, so the complex
eigenvalues of the NOISE-DRESSED generator are the operator-level 'noise-dressed lambda0'. This probe extracts
them and asks: does the leading complex-eigenvalue PHASE sit near the deterministic lambda0 (-45 deg / |arg|=45)
or near the measured Borel pair (~54-63 deg), and how does it MOVE with noise eta?

Object: the Prufer-phase Fokker-Planck generator L(theta; Y, eta) on [0,pi], absorbing at theta=pi, reflecting
at theta=0, at frozen sweep coordinate Y (W(Y)=sign(Y)|Y|^q).  Built with the same Scharfetter-Gummel
discretization validated in _riccati_pde_probe.py; here we eigendecompose the operator instead of time-stepping.
"""
import numpy as np
from _riccati_pde_probe import _bernoulli

def fp_generator(q, eta, Y, N=400):
    """Return the (N x N) non-self-adjoint FP generator matrix L (interior nodes 0..N-1; node N absorbing)."""
    th = np.linspace(0.0, np.pi, N + 1); h = th[1] - th[0]
    thf = 0.5*(th[:-1] + th[1:]); sf = np.sin(thf); cf = np.cos(thf)
    Df = 0.5*eta**2*sf**4
    W = np.sign(Y)*abs(Y)**q
    a = cf*cf - W*sf*sf - eta**2*sf**3*cf            # effective advection drift at faces
    Dsafe = np.where(Df > 1e-300, Df, 1e-300)
    Pe = a*h/Dsafe
    alpha = (Df/h)*_bernoulli(-Pe)                    # F_{k+1/2} = alpha_k C_k - beta_k C_{k+1}
    beta  = (Df/h)*_bernoulli(Pe)
    L = np.zeros((N, N)); inv_h = 1.0/h
    for k in range(N):
        am = alpha[k-1] if k >= 1 else 0.0
        bm = beta[k-1]  if k >= 1 else 0.0
        L[k, k] = (-(alpha[k] + bm))*inv_h if k >= 1 else (-(alpha[k]))*inv_h
        if k >= 1:      L[k, k-1] = am*inv_h
        if k <= N - 2:  L[k, k+1] = beta[k]*inv_h
    return L

def leading_eigs(L, ntop=8):
    ev = np.linalg.eigvals(L)
    ev = ev[np.argsort(-ev.real)]                     # slowest-decaying first (real part closest to 0, least negative)
    return ev[:ntop]

def main():
    eta_phys = np.sqrt(2.0)
    print("="*80)
    print("Complex spectrum of the noise-dressed first-explosion generator (q=2 cusp)")
    print("  deterministic resonance lambda0 = 0.8896 - 0.8896i  (|arg| = 45 deg)")
    print("  measured Borel pair phase ~ 54-63 deg")
    print("="*80)
    q = 2.0
    for Y in [-2.19, -1.5, -1.0, -0.5]:
        print(f"\n  frozen Y = {Y}  (W = {np.sign(Y)*abs(Y)**q:+.3f}):")
        print(f"    {'eta':>6} | leading complex eigenvalues  Re + Im i   (|arg| in deg)")
        for eta in [0.3, 0.6, eta_phys, 2.0, 3.0]:
            L = fp_generator(q, eta, Y, N=400)
            ev = leading_eigs(L, ntop=6)
            comps = [z for z in ev if abs(z.imag) > 1e-6*(abs(z)+1e-30)]
            if comps:
                z = comps[0]
                ang = np.degrees(np.arctan2(abs(z.imag), abs(z.real)))
                tag = "eta=sqrt2" if abs(eta-eta_phys) < 1e-6 else ""
                print(f"    {eta:6.3f} | {z.real:+.4f} {z.imag:+.4f}i   |arg|={ang:5.1f} deg   {tag}")
            else:
                print(f"    {eta:6.3f} | (leading eigenvalues all real: {', '.join(f'{z.real:+.3f}' for z in ev[:3])})")

if __name__ == "__main__":
    main()
