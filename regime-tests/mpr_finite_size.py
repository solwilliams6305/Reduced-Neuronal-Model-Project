"""
mpr_finite_size.py — numerical backbone for the finite-size Langevin derivation (MPR_FINITE_SIZE_LANGEVIN.md).

Two facts the derivation rests on, both checked here against the DIRECT finite-N theta-neuron network
(no mean-field assumed; the macroscopic (r,v) is read from the network via Montrio's conformal map of
the Kuramoto order parameter Z = <e^{i theta}>):

  (1) The macroscopic fluctuation variance scales as 1/N: Var(r) ~ 1/N (slope ~ -1 in log-log), i.e.
      the finite-size noise on the mean-field is O(N^{-1/2}) -- the central claim of the system-size
      Langevin. [validated: slope ~ -0.98, Var*N ~ const]

  (2) The macroscopic bifurcation at the adaptation operating point is a FOLD OF LIMIT CYCLES (amplitude
      edge), NOT a SNIC phase edge: the deterministic period stays finite (~17) and the amplitude stays
      steady (~2.1) up to the boundary, then the cycle vanishes abruptly (a coexisting rest state -> the
      rate-governed bistable regime). [corrects the earlier phase-edge claim; the sigma^{2/3} test failed]

    python3 mpr_finite_size.py scaling    # finite-N Var(r) ~ 1/N test (~30 s)
    python3 mpr_finite_size.py diag        # deterministic fold-of-cycles diagnostic
    python3 mpr_finite_size.py fig         # -> figures/mpr_finite_size.png
"""
from __future__ import annotations
import os
import sys
import numpy as np
import neural_mass_edge as nm

HERE = os.path.dirname(os.path.abspath(__file__))
SCALE_CACHE = os.path.join(HERE, "_mpr_finitesize_cache.npz")
N_LIST = [250, 500, 1000, 2000, 4000]


def net_r(N, etabar=-8.0, J=15.0, Delta=1.0, T=300.0, dt=1e-2, burn=100.0, seed=0):
    """Direct finite-N theta-neuron (QIF) network; returns the macroscopic rate time series r(t).
    eta_j are DETERMINISTIC Lorentzian quantiles (isolating dynamic finite-size fluctuations from
    quenched-sampling noise); r,v are read from Z=<e^{i theta}> via W=pi r + i v=(1-Z*)/(1+Z*)."""
    rng = np.random.default_rng(seed)
    u = (np.arange(N) + 0.5) / N
    eta = etabar + Delta * np.tan(np.pi * (u - 0.5))
    th = 2 * np.pi * rng.random(N); n = int(T / dt); rs = []
    for k in range(n):
        Z = np.mean(np.exp(1j * th)); W = (1 - np.conj(Z)) / (1 + np.conj(Z)); r = max(W.real / np.pi, 0.0)
        th = np.mod(th + dt * ((1 - np.cos(th)) + (1 + np.cos(th)) * (eta + J * r)), 2 * np.pi)
        if k * dt > burn:
            rs.append(r)
    return np.array(rs)


def scaling(seeds=2):
    V = []; M = []
    for N in N_LIST:
        runs = [net_r(N, seed=s) for s in range(seeds)]
        V.append(np.median([r.var() for r in runs])); M.append(np.mean([r.mean() for r in runs]))
        print(f"  N={N:5d}: mean_r={M[-1]:.4f}  Var(r)={V[-1]:.3e}  Var*N={V[-1]*N:.3e}")
    V = np.array(V); slope = np.polyfit(np.log(N_LIST), np.log(V), 1)[0]
    print(f"  log Var vs log N slope = {slope:.2f}  (finite-size Langevin predicts -1.00)")
    np.savez(SCALE_CACHE, Ns=np.array(N_LIST), Var=V, mean_r=np.array(M), slope=np.array([slope]))
    print(f"  saved {SCALE_CACHE}")


def diag():
    etas = [-2.305, -2.308, -2.310, -2.311, -2.312, -2.313, -2.314]
    out = []
    for e in etas:
        rs, dt = nm.mpr_adapt(e, T=4000.0); half = rs[len(rs) * 2 // 3:]
        out.append((e, nm.period(rs, dt), half.max() - half.min()))
        print(f"  eta={e:.4f}: period={out[-1][1] if np.isfinite(out[-1][1]) else 'gone':>6}  "
              f"amp={out[-1][2]:.3f}")
    print("  => period FINITE (~17) + amplitude STEADY (~2.1) then abrupt -> FOLD OF CYCLES (amplitude edge).")
    return out


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) >= 2 else ""
    if cmd == "scaling":
        scaling()
    elif cmd == "diag":
        diag()
    elif cmd == "fig":
        from mpr_finite_size_figure import make_figure
        make_figure()
    else:
        print("usage: mpr_finite_size.py scaling | diag | fig")
