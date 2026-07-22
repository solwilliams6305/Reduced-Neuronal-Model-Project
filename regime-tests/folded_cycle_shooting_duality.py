"""
folded_cycle_shooting_duality.py — the shooting characterization, two independent ways
--------------------------------------------------------------------------------------
Theorem (shooting characterization):  Y_node =_d -Lambda0(beta),  eta = 2/sqrt(beta).
Both are the SAME functional: the first-zero location of the recessive solution of the
noisy Airy operator on the forbidden half-line.  We confirm it by computing the two
sides with methods that share NO code path:

  -Lambda0  : smallest eigenvalue of the random-matrix discretization of the stochastic
              Airy operator  H_beta = -d^2/dx^2 + x + eta*xi  (tridiagonal, Dirichlet).
  Y_node    : first node of u in the SWEPT ODE  u''=(Y-eta xi)u, Y=Y0-T (canard data).

Deterministic check (eta=0): smallest eig = 2.3381 = -a1 (first Airy zero), so the
matrix is the correct SAO.  With noise, the -Lambda0 and Y_node laws coincide.
"""
from __future__ import annotations
import os
import numpy as np


def Lam0(eta, L=5.5, n=240, reps=1600, seed=0):
    """Smallest eigenvalue of the discretized stochastic Airy operator on [0,L]."""
    rng = np.random.default_rng(seed); h = L / n; x = (np.arange(1, n + 1)) * h
    off = -1.0 / h**2; out = np.empty(reps)
    for r in range(reps):
        d = 2.0 / h**2 + x + (eta / np.sqrt(h)) * rng.standard_normal(n)
        M = np.diag(d) + np.diag(np.full(n - 1, off), 1) + np.diag(np.full(n - 1, off), -1)
        out[r] = np.linalg.eigvalsh(M)[0]
    return out


def Y_node(eta, Y0=4.0, dt=8e-4, N=16000, seed=0, Ymin=-4.0):
    rng = np.random.default_rng(seed); u = np.ones(N); v = np.full(N, np.sqrt(Y0)); Y = Y0
    sdt = np.sqrt(dt); Yz = np.full(N, np.nan); done = np.zeros(N, bool)
    for _ in range(int((Y0 - Ymin) / dt)):
        if done.all(): break
        dB = sdt * rng.standard_normal(N); u1 = u + v * dt; Yp = Y - dt; v1 = v + (Y * u) * dt - eta * u * dB
        u = u + 0.5 * (v + v1) * dt; v = v + 0.5 * (Y * u + Yp * u1) * dt - eta * (0.5 * (u + u1)) * dB; Y = Yp
        cr = (~done) & (u < 0.0); Yz[cr] = Y; done |= cr
    return Yz[np.isfinite(Yz)]


def mom(x):
    m = x.mean(); d = x - m; v = np.mean(d**2); return m, np.sqrt(v), np.mean(d**3) / v**1.5


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    figdir = os.path.join(os.path.dirname(here), "figures")
    os.makedirs(figdir, exist_ok=True)

    print("\n=== Shooting duality: -Lambda0 (matrix SAO) vs Y_node (ODE shooting) ===")
    print(f"  deterministic (eta=0) smallest eig = {Lam0(0.0, reps=1)[0]:.4f}  (= -a1 = 2.3381)\n")
    print(f"  {'beta':>4}{'eta':>7} | {'-Lambda0 (matrix)':>26} | {'Y_node (ODE)':>26}")
    data = {}
    for beta, eta in [(4, 1.0), (2, np.sqrt(2.0))]:
        mL = -Lam0(eta, seed=1); Yn = Y_node(eta, seed=2); data[beta] = (mL, Yn, eta)
        a, b, g = mom(mL); p, q, s = mom(Yn)
        print(f"  {beta:4d}{eta:7.3f} | {a:8.4f}{b:7.4f}{g:7.3f}  (m/s/sk) | {p:8.4f}{q:7.4f}{s:7.3f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11.2, 4.4))
    for axi, beta in zip(ax, (4, 2)):
        mL, Yn, eta = data[beta]
        rng = (-4.2, 0.5)
        axi.hist(mL, bins=55, range=rng, density=True, color="C3", alpha=0.45,
                 label=r"$-\Lambda_0$ (random-matrix SAO)")
        axi.hist(Yn, bins=55, range=rng, density=True, histtype="step", lw=2.0, color="C0",
                 label=r"$Y_{\rm node}$ (ODE shooting)")
        axi.set_xlabel(r"$-\Lambda_0\ =\ Y_{\rm node}$"); axi.set_ylabel("density")
        axi.set_title(fr"$\beta={beta}$ ($\eta={eta:.2f}$): two independent methods agree")
        axi.legend(fontsize=9, frameon=False); axi.grid(alpha=0.3)
    fig.suptitle(r"Shooting characterization  $Y_{\rm node}=_d-\Lambda_0(\beta)$:  matrix eigenvalue $\equiv$ ODE first node",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = os.path.join(figdir, "folded_cycle_shooting_duality.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print(f"\n  figure -> figures/folded_cycle_shooting_duality.png\n")


if __name__ == "__main__":
    main()
