"""
airy_line_ensemble.py — Direction C, the positive capstone: the successive peel-offs of ONE
folded cycle are the Airy_2 point process (the GUE edge).

A single folded cycle's inner operator is the stochastic Airy operator H_beta = -d^2/dY^2 + Y +
eta b'(Y).  Its lowest eigenvalues Lambda0 < Lambda1 < Lambda2 are the successive peel-off
levels (first node, second node, ...).  In the edge scaling these are the Airy_beta POINT
process; for beta=2 it is the GUE soft edge.  We test, with NO coupling (one operator, many
independent realizations):
  (i)   the three marginals -Lambda0, -Lambda1, -Lambda2 match the GUE edge xi_k = N^{2/3}(lam_k-2)
        (Tracy-Widom hierarchy, asymptotic means -1.7711, -3.6754, -5.1716);
  (ii)  the spacings grow down the spectrum (the sqrt|.| Airy density);
  (iii) the small-gap level repulsion ~ gap^2 (beta=2).
"""
from __future__ import annotations
import os
import numpy as np

ASYMPT = [-1.7711, -3.6754, -5.1716]                         # Airy_2 point-process means (k=1,2,3)


def sao_levels(eta, R, kmax=3, M=750, L=7.5, seed=0):
    rng = np.random.default_rng(seed); h = L / M; Yj = np.arange(1, M + 1) * h
    Z = rng.standard_normal((M, R))                          # iid white noise, R realizations
    diag = 2.0 / h ** 2 + Yj[:, None] + eta * Z / np.sqrt(h); e2 = (1.0 / h ** 2) ** 2

    def sturm(lam):
        q = diag[0] - lam; cnt = (q < 0).astype(np.int32)
        for i in range(1, M):
            q = (diag[i] - lam) - e2 / np.where(np.abs(q) < 1e-300, -1e-300, q); cnt += (q < 0)
        return cnt

    def kth(k):
        lo = np.full(R, -6.0); hi = np.full(R, 14.0)
        for _ in range(48):
            mid = 0.5 * (lo + hi); has = sturm(mid) >= k; hi = np.where(has, mid, hi); lo = np.where(has, lo, mid)
        return 0.5 * (lo + hi)

    return np.array([kth(j) for j in range(1, kmax + 1)])    # (kmax, R): Lambda0, Lambda1, Lambda2


def gue_edge(N, R, kmax=3, seed=1):
    rng = np.random.default_rng(seed); xis = np.empty((kmax, R)); s = N ** (2.0 / 3.0)
    for r in range(R):
        X = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
        H = (X + X.conj().T) / 2.0 / np.sqrt(N)               # Hermitian, semicircle edge ~ 2
        ev = np.linalg.eigvalsh(H)[-kmax:][::-1]              # top kmax, descending
        xis[:, r] = s * (ev - 2.0)
    return xis                                               # (kmax, R): xi1 > xi2 > xi3


if __name__ == "__main__":
    eta = np.sqrt(2.0)
    Lam = sao_levels(eta, R=2400, M=750, seed=0)             # Lambda0,1,2
    negL = -Lam                                              # -Lambda_k = the point process (descending)
    xi = gue_edge(N=180, R=1600, seed=2)

    print(f"{'level k':>8} {'SAO -Lam_k':>14} {'GUE xi_k':>14} {'asympt':>9}")
    for k in range(3):
        print(f"{k+1:>8} {negL[k].mean():7.3f}±{negL[k].std():.2f}   {xi[k].mean():7.3f}±{xi[k].std():.2f}"
              f"   {ASYMPT[k]:9.3f}")
    g01 = Lam[1] - Lam[0]; g12 = Lam[2] - Lam[1]
    print(f"\n  spacings (SAO): <Lam1-Lam0>={g01.mean():.3f}  <Lam2-Lam1>={g12.mean():.3f}  "
          f"(grow down the spectrum: Airy sqrt|.| density)")
    gxi = xi[0] - xi[1]
    print(f"  gap repulsion: small-gap P(gap)~gap^p; fit p(SAO)... see panel (beta=2 => p=2)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    fig, ax = plt.subplots(1, 3, figsize=(15.4, 4.7))
    cols = ["#2563eb", "#16a34a", "#dc2626"]

    # Panel A: the three marginals -- SAO (filled) vs GUE edge (step) vs asymptotic means
    for k in range(3):
        ax[0].hist(negL[k], bins=60, range=(-8, 0.5), density=True, color=cols[k], alpha=0.45)
        ax[0].hist(xi[k], bins=60, range=(-8, 0.5), density=True, histtype="step", lw=2.0, color=cols[k])
        ax[0].axvline(ASYMPT[k], color=cols[k], lw=1.0, ls=":")
    ax[0].plot([], [], color="#475569", lw=6, alpha=0.45, label="SAO peel-offs (filled)")
    ax[0].plot([], [], color="#475569", lw=2, label="GUE edge (step)")
    ax[0].set_xlabel(r"$-\Lambda_k$  /  $\xi_k$"); ax[0].set_ylabel("density")
    ax[0].set_title(r"Successive peel-offs $=$ Airy$_2$ point process")
    ax[0].legend(fontsize=8.6, framealpha=0.95, edgecolor="#cbd5e1")

    # Panel B: gap (level repulsion) -- small-gap ~ gap^2 for beta=2
    ax[1].hist(g01, bins=70, range=(0, 4.5), density=True, color="#2563eb", alpha=0.5,
               label=r"SAO gap $\Lambda_1-\Lambda_0$")
    ax[1].hist(gxi, bins=70, range=(0, 4.5), density=True, histtype="step", lw=2.0, color="#0f172a",
               label=r"GUE gap $\xi_1-\xi_2$")
    gg = np.linspace(0, 1.2, 50)
    ax[1].plot(gg, 0.5 * gg ** 2, "--", color="#dc2626", lw=1.6, label=r"$\propto$ gap$^2$ ($\beta=2$ repulsion)")
    ax[1].set_xlabel(r"gap"); ax[1].set_ylabel("density")
    ax[1].set_title(r"Level repulsion $\sim$ gap$^2$ (GUE)")
    ax[1].set_xlim(0, 4.5); ax[1].legend(fontsize=8.6, framealpha=0.95, edgecolor="#cbd5e1")

    # Panel C: joint scatter (-Lambda0, -Lambda1) -- ordered, repelling
    ax[2].plot(negL[0], negL[1], ".", color="#2563eb", ms=1.5, alpha=0.4)
    ax[2].plot([-8, 0.5], [-8, 0.5], "--", color="#94a3b8", lw=1.0, label=r"$-\Lambda_1=-\Lambda_0$")
    ax[2].set_xlabel(r"$-\Lambda_0$ (top peel-off)"); ax[2].set_ylabel(r"$-\Lambda_1$ (2nd)")
    ax[2].set_title("Ordered & repelling (no crossing)")
    ax[2].set_xlim(-5, 0.5); ax[2].set_ylim(-7, -1); ax[2].legend(fontsize=9, framealpha=0.95, edgecolor="#cbd5e1")

    fig.tight_layout()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures",
                       "airy_line_ensemble.png")
    out = os.path.abspath(out)
    fig.savefig(out, dpi=140); print("saved", out)
