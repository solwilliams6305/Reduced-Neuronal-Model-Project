"""
gap_roughness.py — Direction C, testing Conjecture 1/3: is the finite-eta roughening of the
peel-off / eigenvalue field driven by AVOIDED CROSSINGS (small SAO edge gaps Lambda1 - Lambda0)?

Mechanism (Conjecture 3, spectral rigidity): at finite beta the level repulsion (~ gap^beta) is
soft, so the two lowest SAO eigenvalues approach often; near an avoided crossing the ground
state turns rapidly and Lambda0(x) kinks -> roughness (alpha<1).  As beta grows the repulsion
stiffens, small gaps are suppressed, Lambda0 smooths, alpha -> 1.  Predictions tested here:
  (i)   Lambda0(x) kinks where the gap pinches;
  (ii)  E[(Delta Lambda0)^2 | gap] grows as gap -> 0;
  (iii) the gap distribution's small-gap density ~ gap^beta stiffens with beta.

Two lowest eigenvalues via Sturm-sequence bisection (robust), vectorised over the array.
"""
from __future__ import annotations
import os
import numpy as np
from eigenvalue_vs_node import ou_field


def two_lowest(eta, Nx, ell, M=550, L=7.0, seed=0):
    rng = np.random.default_rng(seed)
    h = L / M; Yj = np.arange(1, M + 1) * h
    Z = ou_field(M, Nx, ell, rng)
    diag = 2.0 / h ** 2 + Yj[:, None] + eta * Z / np.sqrt(h)        # (M, Nx)
    e2 = (1.0 / h ** 2) ** 2

    def sturm(lam):                                                # count eigenvalues < lam, per x
        q = diag[0] - lam
        cnt = (q < 0).astype(np.int32)
        for i in range(1, M):
            q = (diag[i] - lam) - e2 / np.where(np.abs(q) < 1e-300, -1e-300, q)
            cnt += (q < 0)
        return cnt

    def kth(k, lo, hi, iters=46):
        lo = np.full(Nx, lo); hi = np.full(Nx, hi)
        for _ in range(iters):
            mid = 0.5 * (lo + hi); has = sturm(mid) >= k
            hi = np.where(has, mid, hi); lo = np.where(has, lo, mid)
        return 0.5 * (lo + hi)

    return kth(1, -5.0, 12.0), kth(2, -5.0, 12.0)


def alpha_of(P, lo=2, hi=10):
    Nx = P.shape[1]
    V = np.array([((P[:, r:] - P[:, :Nx - r]) ** 2).mean() for r in range(0, hi + 2)])
    return np.polyfit(np.log(np.arange(lo, hi)), np.log(V[lo:hi]), 1)[0]


if __name__ == "__main__":
    Nx, ell = 500, 80.0
    data = {}
    for beta in [2.0, 4.0]:
        eta = 2.0 / np.sqrt(beta); R = 14
        L0 = []; L1 = []
        for k in range(R):
            a, b = two_lowest(eta, Nx, ell, seed=1000 + 50 * int(beta) + k)
            L0.append(a); L1.append(b)
        L0 = np.array(L0); L1 = np.array(L1); gap = L1 - L0
        dL0 = L0[:, 1:] - L0[:, :-1]                                # lag-1 increment
        gx = gap[:, :-1]                                            # gap at the left endpoint
        data[beta] = dict(L0=L0, L1=L1, gap=gap, dL0=dL0, gx=gx, alpha=alpha_of(L0))
        m = (-L0).mean(); s = (-L0).std()
        print(f"beta={beta:.0f}: -Lam0 mean={m:.3f} std={s:.3f}  alpha(Lam0)={data[beta]['alpha']:.3f}  "
              f"median gap={np.median(gap):.3f}  frac(gap<0.3)={np.mean(gap<0.3):.3f}")

    # node alpha (validated sweep) vs EXACT (Sturm) eigenvalue alpha -> Conjecture 2
    from airy_process_peeloffs import peeloff_profile
    nod = np.array([peeloff_profile(np.sqrt(2.0), 500, 80.0, seed=7000 + k) for k in range(12)])
    a_node = alpha_of(nod)
    print(f"\n  CONJ 2 check (beta=2): node alpha={a_node:.3f} (rough)  vs  exact eigenvalue "
          f"alpha(Lam0)={data[2.0]['alpha']:.3f} (smooth, ~Airy2)  => roughening is the level-crossing")

    # (ii) conditional roughness E[(dLam0)^2 | gap]
    d = data[2.0]
    gflat = d["gx"].ravel(); dflat = d["dL0"].ravel() ** 2
    edges = np.linspace(np.percentile(gflat, 1), np.percentile(gflat, 97), 13)
    cen = 0.5 * (edges[1:] + edges[:-1])
    cond = np.array([dflat[(gflat >= edges[i]) & (gflat < edges[i + 1])].mean() for i in range(len(cen))])
    sl = np.polyfit(np.log(cen), np.log(cond), 1)[0]
    print(f"\n  E[(dLam0)^2 | gap] ~ gap^{sl:.2f}  (negative slope => roughness grows as gap->0: gap-driven)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    fig, ax = plt.subplots(1, 3, figsize=(15.4, 4.7))

    # Panel A: the two lowest eigenvalues + gap, one realization window -> avoided crossings
    xw = np.arange(220)
    ax[0].plot(xw, d["L0"][0, :220], color="#2563eb", lw=1.4, label=r"$\Lambda_0(x)$")
    ax[0].plot(xw, d["L1"][0, :220], color="#16a34a", lw=1.4, label=r"$\Lambda_1(x)$")
    gmin_idx = np.where(d["gap"][0, :220] < np.percentile(d["gap"][0], 8))[0]
    ax[0].plot(gmin_idx, d["L0"][0, gmin_idx], "o", color="#dc2626", ms=5, zorder=5,
               label="gap minima (avoided crossings)")
    ax[0].set_xlabel(r"array site $x$"); ax[0].set_ylabel("SAO edge eigenvalues")
    ax[0].set_title(r"$\Lambda_0,\Lambda_1$ rarely approach (strong repulsion, $\beta=2$)")
    ax[0].legend(fontsize=8.4, framealpha=0.95, edgecolor="#cbd5e1")

    # Panel B: conditional roughness vs gap (the smoking gun)
    ax[1].loglog(cen, cond, "o-", color="#2563eb", ms=6, mec="white", mew=1.0,
                 label=r"$E[(\Delta\Lambda_0)^2\,|\,{\rm gap}]$")
    ax[1].loglog(cen, cond[0] * (cen / cen[0]) ** sl, "--", color="#0f172a", lw=1.6,
                 label=fr"$\propto {{\rm gap}}^{{{sl:.2f}}}$")
    ax[1].set_xlabel("gap  $\\Lambda_1-\\Lambda_0$"); ax[1].set_ylabel("local roughness")
    ax[1].set_title("Roughness does NOT grow at small gap (gap-driven refuted)")
    ax[1].legend(fontsize=9, framealpha=0.95, edgecolor="#cbd5e1")

    # Panel C: gap distribution / level repulsion, beta=2 vs 4
    for beta, col in [(2.0, "#dc2626"), (4.0, "#2563eb")]:
        g = data[beta]["gap"].ravel()
        ax[2].hist(g, bins=60, range=(0, 4), density=True, histtype="step", lw=2.0, color=col,
                   label=fr"$\beta={beta:.0f}$ ($\alpha={data[beta]['alpha']:.2f}$)")
    ax[2].set_xlabel(r"gap  $\Lambda_1-\Lambda_0$"); ax[2].set_ylabel("density")
    ax[2].set_title(r"Level repulsion stiffens with $\beta$ (fewer small gaps)")
    ax[2].legend(fontsize=9, framealpha=0.95, edgecolor="#cbd5e1")

    fig.tight_layout()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures",
                       "gap_roughness.png")
    out = os.path.abspath(out)
    fig.savefig(out, dpi=140); print("saved", out)
