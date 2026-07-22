"""
twotime_folded_pde.py
=====================
Capstone: two-time disordered folded medium. The leading escape edge
xi(tau) = N^{1/6}(lambda_top(C(tau)) - 2 sqrt(N)) under a PHYSICAL slow drive of the
disordered coupling C(tau). Which slow drive yields the Airy_2 process?

Discriminator = small-lag exponent of the increment variance Var[xi(tau+D)-xi(tau)]:
  Brownian (slope 1) = Airy_2 ;  ballistic (slope 2) = smooth, NOT Airy_2.

Result (see FOLDED_PDE_TWOTIME.md):
  * stochastic coupling drift (noisy plasticity)  -> exponent ~1, plateau ~2 Var(TW2)=1.63 : AIRY_2
  * deterministic slow drive (rotation)           -> exponent ~2 (ballistic)               : NOT Airy_2
  * BOTH share the same instantaneous marginal = TW2 -> the marginal cannot distinguish them.
Conclusion: Airy_2 = edge of Dyson BROWNIAN motion; the slow drift must be STOCHASTIC.
"""
import numpy as np

def gue(N, rng):
    A = rng.standard_normal((N, N)) + 1j*rng.standard_normal((N, N)); return (A + A.conj().T)/2

def edge(lam, N):
    return N**(1/6.)*(lam - 2*np.sqrt(N))

def traj_stochastic(N=96, T=2200, corr=50, seed=0):
    rng = np.random.default_rng(seed); a = np.exp(-1.0/corr); C = gue(N, rng); out = np.empty(T)
    for t in range(T):
        out[t] = np.linalg.eigvalsh(C)[-1]; C = a*C + np.sqrt(1-a*a)*gue(N, rng)
    return edge(out, N)

def traj_deterministic(N=96, T=2200, period=260, seed=0):
    rng = np.random.default_rng(seed); A = gue(N, rng); B = gue(N, rng); w = 2*np.pi/period; out = np.empty(T)
    for t in range(T):
        out[t] = np.linalg.eigvalsh(np.cos(w*t)*A + np.sin(w*t)*B)[-1]   # GUE at every tau, smooth motion
    return edge(out, N)

def exponent(xi, lags=(1, 2, 4, 8)):
    lags = np.array(lags); V = np.array([(xi[L:]-xi[:-L]).var() for L in lags])
    return np.polyfit(np.log(lags), np.log(V), 1)[0], V

if __name__ == "__main__":
    xs = traj_stochastic(); xd = traj_deterministic()
    es, Vs = exponent(xs); ed_, Vd = exponent(xd)
    lags = np.array([1, 2, 4, 8, 16, 32, 64, 128])
    Vs = np.array([(xs[L:]-xs[:-L]).var() for L in lags])
    print(f"stochastic drift : exponent {es:.2f} (Brownian=1 -> Airy_2); plateau {Vs[-2]:.2f} vs 2Var(TW2)=1.63")
    print(f"deterministic    : exponent {ed_:.2f} (ballistic=2 -> NOT Airy_2)")
    print("=> Airy_2 needs a STOCHASTIC slow drift of disordered long-range couplings.")
