"""
eigenvalue_vs_node.py — Direction C, decisive test of Conjecture 2.

Is the finite-eta roughening (alpha<1 at beta=2) in the LEVEL-CROSSING (the node Y_node=G^{-1}(0))
rather than in the spectrum?  We compute, from the SAME correlated-noise array at beta=2:
  (a) node              Y_node(x)   -- the first-node (peel-off) location  [a moving zero];
  (b) smooth Riccati    R(Yref;x) = -u'/u at a fixed level above the node  [smooth spectral data];
  (c) eigenvalue        Lambda0(x)  -- SAO ground state by inverse iteration on the discretised
                                       operator -d^2/dY^2 + Y + eta b'(Y) (Dirichlet)  [-Lambda0 =d TW_beta].
and compare local roughness alpha (V(r)=Var(Delta) ~ r^alpha).  Prediction (Conjecture 2):
alpha(R), alpha(Lambda0) ~ 1 (clean Airy_2) while alpha(Y_node) ~ 0.75 -> the roughening is the
level-crossing.
"""
from __future__ import annotations
import os
import numpy as np


def ou_field(M, Nx, ell, rng):
    a = np.exp(-1.0 / ell); b = np.sqrt(1.0 - a * a)
    Z = np.empty((M, Nx)); Z[:, 0] = rng.standard_normal(M)
    for x in range(1, Nx):
        Z[:, x] = a * Z[:, x - 1] + b * rng.standard_normal(M)
    return Z


def node_and_riccati(eta, Nx, ell, Y0=5.0, Yend=-4.0, Yref=1.0, dt=2e-3, seed=0):
    rng = np.random.default_rng(seed)
    nsteps = int((Y0 - Yend) / dt); sdt = np.sqrt(dt)
    Z = ou_field(nsteps, Nx, ell, rng)
    r0 = np.sqrt(Y0) + 1.0 / (4 * Y0)
    u = np.ones(Nx); w = np.full(Nx, r0); yn = np.full(Nx, np.nan); al = np.ones(Nx, bool)
    Rref = np.full(Nx, np.nan); mref = int((Y0 - Yref) / dt)
    for m in range(nsteps):
        Y = Y0 - dt * m
        if m == mref:
            Rref[:] = w / u                                    # smooth Riccati at Yref (above node)
        un = u + w * dt
        wn = w + (Y * u) * dt - eta * u * sdt * Z[m]
        cr = al & (un * u <= 0.0)
        yn[cr] = Y - dt * (u[cr] / (u[cr] - un[cr])); al[cr] = False
        u, w = un, wn
    yn[np.isnan(yn)] = Yend
    return yn, Rref


def lambda0_profile(eta, Nx, ell, M=320, L=6.5, sigma=1.4, iters=40, seed=0):
    """SAO ground state Lambda0(x): -d^2/dY^2 + Y + eta b'(Y) on [0,L], Dirichlet; inverse iteration."""
    rng = np.random.default_rng(seed)
    h = L / M; Yj = (np.arange(1, M + 1)) * h
    Z = ou_field(M, Nx, ell, rng)
    diag = 2.0 / h ** 2 + Yj[:, None] + eta * Z / np.sqrt(h)    # (M, Nx)
    off = -1.0 / h ** 2
    d = diag - sigma
    v = rng.standard_normal((M, Nx)); v /= np.linalg.norm(v, axis=0)
    cp = np.empty((M, Nx)); dp = np.empty((M, Nx)); wsol = np.empty((M, Nx))
    for _ in range(iters):                                     # Thomas solve (A-sigma) wsol = v
        cp[0] = off / d[0]; dp[0] = v[0] / d[0]
        for j in range(1, M):
            mm = d[j] - off * cp[j - 1]
            cp[j] = off / mm; dp[j] = (v[j] - off * dp[j - 1]) / mm
        wsol[-1] = dp[-1]
        for j in range(M - 2, -1, -1):
            wsol[j] = dp[j] - cp[j] * wsol[j + 1]
        v = wsol / np.linalg.norm(wsol, axis=0)
    Av = diag * v; Av[:-1] += off * v[1:]; Av[1:] += off * v[:-1]
    return (v * Av).sum(axis=0)                                 # Rayleigh quotient = Lambda0(x)


def alpha_of(profiles, lo=2, hi=10):
    Nx = profiles.shape[1]
    V = np.array([((profiles[:, r:] - profiles[:, :Nx - r]) ** 2).mean() for r in range(0, hi + 2)])
    return np.polyfit(np.log(np.arange(lo, hi)), np.log(V[lo:hi]), 1)[0], V


if __name__ == "__main__":
    eta = np.sqrt(2.0); Nx, ell, R = 350, 80.0, 20
    res = [node_and_riccati(eta, Nx, ell, seed=500 + k) for k in range(R)]
    YN = np.array([r[0] for r in res])
    L0 = np.array([lambda0_profile(eta, Nx, ell, M=800, L=7.0, iters=60, seed=900 + k) for k in range(R)])

    aYN, Vyn = alpha_of(YN); aL0, Vl0 = alpha_of(L0)
    print(f"{'quantity':>22} {'mean':>8} {'std':>7} {'skew':>7} {'alpha':>7}")
    for name, P, a in [("node  Y_node(x)", YN, aYN), ("eigenvalue -Lambda0", -L0, aL0)]:
        m, s = P.mean(), P.std(); sk = ((P - m) ** 3).mean() / s ** 3
        print(f"{name:>22} {m:8.3f} {s:7.3f} {sk:7.3f} {a:7.3f}")
    print(f"\n  TW_2 marginal target: mean -1.771 / std 0.902 / skew 0.224 (eigenvalue still grid-converging)")
    verdict = ("eigenvalue is ALSO rough -> roughening is INTRINSIC to the spectrum "
               "(Conjecture 2 overturned, not a level-crossing artifact)"
               if aL0 < 0.85 else "eigenvalue ~ Airy_2 (alpha~1): roughening is the level-crossing")
    print(f"  VERDICT: {verdict}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.7))
    r = np.arange(1, 11)
    for nm, V, a, c in [(fr"node $Y_{{\rm node}}$  ($\alpha={aYN:.2f}$)", Vyn, aYN, "#dc2626"),
                        (fr"eigenvalue $\Lambda_0$  ($\alpha={aL0:.2f}$)", Vl0, aL0, "#2563eb")]:
        ax[0].loglog(r, V[1:11] / V[1], "o-", color=c, ms=4, label=nm)
    ax[0].loglog(r, r ** 1.0, "--", color="#0f172a", lw=1.4, label=r"Airy$_2$ slope $r^1$ ($\beta\to\infty$)")
    ax[0].set_xlabel(r"lag $r$"); ax[0].set_ylabel(r"$V(r)/V(1)$ (normalised)")
    ax[0].set_title(r"Roughness: node vs spectrum ($\beta=2$)")
    ax[0].legend(fontsize=8.4, framealpha=0.95, edgecolor="#cbd5e1")

    ax[1].hist(YN.ravel(), bins=70, range=(-5, 1), density=True, color="#dc2626", alpha=0.45,
               label=r"$Y_{\rm node}$")
    ax[1].hist(-L0.ravel(), bins=70, range=(-5, 1), density=True, histtype="step", lw=2.0,
               color="#2563eb", label=r"$-\Lambda_0$ (eigenvalue)")
    ax[1].axvline(-1.771, color="#0f172a", lw=1.0, ls="--", label=r"TW$_2$ mean")
    ax[1].set_xlabel(r"value"); ax[1].set_ylabel("density")
    ax[1].set_title("Both marginals $=$ TW$_2$")
    ax[1].legend(fontsize=8.6, framealpha=0.95, edgecolor="#cbd5e1")

    fig.tight_layout()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures",
                       "eigenvalue_vs_node.png")
    out = os.path.abspath(out)
    fig.savefig(out, dpi=140); print("saved", out)
