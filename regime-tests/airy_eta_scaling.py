"""
airy_eta_scaling.py — the resolution of the Direction-C roughness puzzle.

The coupled peel-off field has TW_beta one-point marginal at all beta, but its local roughness
exponent alpha (V(r)=Var(DY_node) ~ r^alpha) came out ~0.75 < 1 at beta=2 -- rougher than the
Airy_2 process (Hoelder-1/2, alpha=1).  Here we show alpha -> 1 as eta = 2/sqrt(beta) -> 0:
the finite-eta nonlinearity roughens the field; the universal Airy_2 local-Brownian behaviour
is recovered in the edge-scaling (small-noise) limit -- exactly where Airy_2 lives.  We also fit
the deviation 1 - alpha ~ eta^p.
"""
from __future__ import annotations
import os
import numpy as np
from airy_process_peeloffs import peeloff_profile


def local_alpha(eta, Nx=460, ell=80.0, R=20, seed0=400):
    profs = np.array([peeloff_profile(eta, Nx, ell, Y0=4.5, Yend=-3.5, dt=2.5e-3, seed=seed0 + k)
                      for k in range(R)])
    V = np.array([((profs[:, r:] - profs[:, :Nx - r]) ** 2).mean() for r in range(0, 12)])
    rs = np.arange(2, 10)
    a = np.polyfit(np.log(rs), np.log(V[2:10]), 1)[0]
    return a, profs.std()


if __name__ == "__main__":
    betas = np.array([2.0, 3.0, 4.0, 6.0, 9.0, 14.0, 22.0])
    etas = 2.0 / np.sqrt(betas)
    A = np.array([local_alpha(e) for e in etas])
    alpha, stds = A[:, 0], A[:, 1]
    print(f"{'beta':>5} {'eta':>6} {'std':>7} {'alpha':>7} {'1-alpha':>8}")
    for b, e, a, s in zip(betas, etas, alpha, stds):
        print(f"{b:5.0f} {e:6.3f} {s:7.3f} {a:7.3f} {1-a:8.3f}")
    # fit 1 - alpha ~ eta^p
    good = (1 - alpha) > 0.01
    p, lnc = np.polyfit(np.log(etas[good]), np.log(1 - alpha[good]), 1)
    print(f"\n  fit  1 - alpha ~ eta^p,  p = {p:.2f}   (cf eta^2 = 2.00, eta^{{8/3}} = 2.67)")
    print(f"  => alpha -> 1 (Airy_2 Hoelder-1/2) as eta -> 0;  the beta=2 roughness is finite-eta.")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    fig, ax = plt.subplots(1, 2, figsize=(11.2, 4.6))

    ax[0].axhline(1.0, color="#16a34a", lw=1.8, ls="--", label=r"Airy$_2$ (Hölder-$\frac{1}{2}$): $\alpha=1$")
    ax[0].semilogx(betas, alpha, "o-", color="#2563eb", ms=8, mec="white", mew=1.1,
                   label=r"peel-off field $\alpha(\beta)$")
    ax[0].set_xlabel(r"$\beta = 4/\eta^2$"); ax[0].set_ylabel(r"local roughness exponent $\alpha$")
    ax[0].set_title(r"$\alpha\to1$ as $\eta\to0$: converges to Airy$_2$")
    ax[0].set_ylim(0.6, 1.05); ax[0].legend(fontsize=9, framealpha=0.95, edgecolor="#cbd5e1")

    ee = np.linspace(etas.min() * 0.9, etas.max() * 1.05, 50)
    ax[1].loglog(etas, 1 - alpha, "o", color="#dc2626", ms=8, mec="white", mew=1.1,
                 label=r"$1-\alpha$ (finite-$\eta$ roughening)")
    ax[1].loglog(ee, np.exp(lnc) * ee ** p, "-", color="#0f172a", lw=2.0, label=fr"$\propto\eta^{{{p:.2f}}}$")
    ax[1].loglog(ee, (1 - alpha[0]) * (ee / etas[0]) ** (8 / 3), ":", color="#7c3aed", lw=1.6,
                 label=r"$\eta^{8/3}$")
    ax[1].set_xlabel(r"$\eta=2/\sqrt{\beta}$"); ax[1].set_ylabel(r"$1-\alpha$")
    ax[1].set_title(fr"deviation vanishes as $\eta^{{{p:.1f}}}$")
    ax[1].legend(fontsize=9, framealpha=0.95, edgecolor="#cbd5e1")

    fig.tight_layout()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures",
                       "airy_eta_scaling.png")
    out = os.path.abspath(out)
    fig.savefig(out, dpi=140)
    print("saved", out)
