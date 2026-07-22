"""
network_field_sweep.py
=======================
Banked #4 payoff: measure the Airy_2 PROCESS covariance directly from a disordered
folded-cycle network swept through the fold (not just the inherited TW marginal).

Leading peel-off of the swept network = -2.338 - lambda_min(C)  (each eigenmode is an
independent shifted canard, peeling at Y = -2.338 - c_k). Letting the disordered coupling
C(tau) drift as a matrix Ornstein-Uhlenbeck process makes the edge a genuine process;
its increment variance is the Airy_2 covariance (linear -> 2*Var(TW)).

Findings (see FOLDED_PDE_NETWORK_SWEEP.md):
  * marginal = TW2 (var 0.79, skew +0.22);
  * increment variance linear -> saturates at 2*Var(marginal) -> 2*Var(TW2)=1.626;
  * local (nearest-neighbour) drift fails (skew +0.59, no TW);
  * the actual nonlinear Cole-Hopf field reproduces Y_node,k = -2.338 - c_k (corr 0.95).
Real (GOE) couplings give the Airy_1 / TW_1 sibling (plateau 2*Var(TW1)=3.22).
"""
import numpy as np

def gue(N, rng):
    A = rng.standard_normal((N, N)) + 1j*rng.standard_normal((N, N)); return (A + A.conj().T)/2

def ou_edge_trajectory(N=112, T=2200, corr_steps=50, seed=0):
    """Leading-edge process xi(tau) = top eigenvalue of a matrix-OU disordered coupling."""
    rng = np.random.default_rng(seed); a = np.exp(-1.0/corr_steps)
    C = gue(N, rng); lam = np.empty(T)
    for t in range(T):
        lam[t] = np.linalg.eigvalsh(C)[-1]
        C = a*C + np.sqrt(1-a*a)*gue(N, rng)            # stationary matrix-OU
    return N**(1/6.)*(lam - 2*np.sqrt(N))               # edge-rescaled -> TW2 units

def increment_variance(xi, lags):
    return np.array([(xi[L:]-xi[:-L]).var() for L in lags])

def tw2_marginal(N=112, S=2000, seed=1):
    rng = np.random.default_rng(seed)
    ed = np.array([np.linalg.eigvalsh(gue(N, rng))[-1] for _ in range(S)])
    xi = N**(1/6.)*(ed - 2*np.sqrt(N)); v = xi.var()
    return v, ((xi-xi.mean())**3).mean()/v**1.5

def nonlinear_field_sweep(Nf=64, eta=0.1, edge=0.6, Y0=4.0, Ymin=-5.0, dY=0.003, seed=2):
    """Genuine real Cole-Hopf network field swept through the fold; returns (sim, predicted)."""
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((Nf, Nf)); Cf = (A + A.T)/np.sqrt(2*Nf)*(edge/2)   # GOE, spectral edge ~ edge
    c, V = np.linalg.eigh(Cf)
    U = V@np.ones(Nf); W = V@(-np.sqrt(Y0+c))          # each eigenmode on its recessive branch
    Y = Y0; Ynode = np.full(Nf, np.nan); rng2 = np.random.default_rng(seed+7)
    for _ in range(int((Y0-Ymin)/dY)):
        up = V.T@U
        U2 = U - dY*W; W = W - dY*(Y*U + Cf@U) + eta*U*np.sqrt(dY)*rng2.standard_normal(Nf)
        U = U2; Y -= dY
        un = V.T@U; nn = np.isnan(Ynode) & (un*up < 0); Ynode[nn] = Y
        nrm = np.max(np.abs(U)) + 1e-300; U = U/nrm; W = W/nrm
    return Ynode, -2.338 - c

if __name__ == "__main__":
    lags = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256])
    xi = ou_edge_trajectory()
    V = increment_variance(xi, lags); v, s = tw2_marginal()
    print(f"marginal: var={v:.3f} skew={s:+.3f}  (TW2 0.81/+0.22)")
    print(f"increment variance: {np.round(V,2)}")
    print(f"plateau ~ {V[-2]:.2f}  vs  2*Var(TW2)=1.63  (finite-N via 2*Var(marginal)={2*v:.2f})")
    sim, pred = nonlinear_field_sweep()
    print(f"nonlinear sweep: corr(sim, -2.338-c_k) = {np.corrcoef(sim, pred)[0,1]:.3f}; "
          f"leading {np.nanmax(sim):.2f} vs {pred.max():.2f}")
