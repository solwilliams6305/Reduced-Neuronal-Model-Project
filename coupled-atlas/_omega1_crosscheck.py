"""Independent cross-check of Omega_1 (=-1.12 from sum-over-states) by DIRECT DIAGONALIZATION
Monte Carlo -- no eigenstate sum, no Green's function.

Discretize H0 = -d^2/dx^2 + x on [0,L] Dirichlet (symmetric tridiagonal).  Perturb by the
white-noise multiplication operator 2 eps b'(x): on the grid, b'_i ~ N(0, 1/h) (so E[b'_i b'_j]
= delta_ij/h represents delta(x-y)).  H = H0 + diag(2 eps b') stays symmetric tridiagonal.

For each realization: lowest eigenvalue Lambda0(eps); exact first-order shift eps*Lam1 with
Lam1 = 2 sum psi0^2 b' h (psi0 from H0).  Variance-reduced 2nd-order sample:
    Lam2_hat = (Lambda0(eps) - Lambda0_det - eps*Lam1)/eps^2  = Lambda2 + O(eps),
and E[Lam2_hat] = Omega_1 + O(eps^2)  (the O(eps) term ~ E[Lambda3] = 0 by parity).
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal

def run(L=12.0, h=0.02, eps=0.05, NMC=6000, seed=0):
    x = np.arange(h, L, h)                 # interior grid (Dirichlet at 0 and L)
    n = len(x)
    diag0 = 2.0/h**2 + x
    off = -np.ones(n-1)/h**2
    # deterministic ground state
    w0, v0 = eigh_tridiagonal(diag0, off, select='i', select_range=(0, 0))
    L0 = w0[0]; psi0 = v0[:, 0]
    psi0 /= np.sqrt(np.sum(psi0**2)*h)     # normalize wrt measure h
    p0sq_h = psi0**2 * h
    rng = np.random.default_rng(seed)
    samples = []
    for _ in range(NMC):
        bp = rng.normal(0.0, np.sqrt(1.0/h), size=n)     # white noise, E[bp_i bp_j]=delta_ij/h
        d = diag0 + 2.0*eps*bp
        w, _ = eigh_tridiagonal(d, off, select='i', select_range=(0, 0))
        Lam1 = 2.0*np.sum(p0sq_h*bp)                       # exact first-order (coeff of eps)
        Lam2_hat = (w[0] - L0 - eps*Lam1)/eps**2
        samples.append(Lam2_hat)
    samples = np.array(samples)
    return L0, samples.mean(), samples.std()/np.sqrt(NMC), samples.std()

if __name__ == "__main__":
    print("Independent diagonalization-MC cross-check of Omega_1  (sum-over-states gave -1.12)")
    print(f"{'eps':>6} {'L':>5} {'h':>6} {'NMC':>6}  {'Omega_1 (MC)':>16}  {'per-sample sd':>13}")
    for eps in (0.08, 0.05, 0.03):
        L0, mean, se, sd = run(L=12.0, h=0.02, eps=eps, NMC=8000, seed=1)
        print(f"{eps:>6.3f} {12.0:>5.1f} {0.02:>6.3f} {8000:>6}  {mean:>+9.4f} +- {se:.4f}  {sd:>13.3f}")
    print("\ngrid convergence at eps=0.05:")
    for h in (0.03, 0.02, 0.015):
        L0, mean, se, sd = run(L=12.0, h=h, eps=0.05, NMC=8000, seed=2)
        print(f"  h={h:.3f}: L0={L0:.4f} (det 2.3381),  Omega_1 = {mean:+.4f} +- {se:.4f}")
