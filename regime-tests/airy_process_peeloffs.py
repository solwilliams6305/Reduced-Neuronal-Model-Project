"""
airy_process_peeloffs.py — Direction C: the Airy process for successive / coupled peel-offs.

A SINGLE peel-off is the stochastic Airy operator (SAO) ground state, Y_node =d TW_beta (the
1-point marginal of the Airy_beta process).  A SPATIAL ARRAY of folded cycles whose inner noise
is white in the sweep variable Y but CORRELATED across the array x (the coupling) is the SAO
under stationary noise dynamics -- the defining construction of the (stationary) Airy_beta
process (the SAO is the edge of the beta-ensemble; evolving its noise as Ornstein-Uhlenbeck is
edge-of-Dyson-OU, whose limit is the Airy line ensemble).

We compute the peel-off profile  Y_node(x)  (recessive Cole-Hopf first node per site x) with
OU-correlated noise across x, and test the Airy-process signatures:
  (i)   one-point marginal  Y_node(x) =d TW_beta   (mean/std/skew);
  (ii)  Hölder-1/2 roughness:  V(r) = Var(Y_node(x+r) - Y_node(x)) ~ c r^alpha,  alpha = 1
        (locally Brownian -- the Airy signature; a SMOOTH field would give alpha = 2), with the
        linear-response slope  dV/dr -> 2 eta^2 (∫ psi0^4) / ell  (the same Airy integral
        ∫ psi0^4 = 0.414 that sets the bulk SAO width), and saturation V(inf) = 2 Var(TW_beta);
  (iii) decorrelation:  C(r) = Cov(Y_node(x), Y_node(x+r))  decays from Var over ~ell.
"""
from __future__ import annotations
import os
import numpy as np

TW2 = dict(mean=-1.7711, std=0.9018, skew=0.2241)          # one-point target (beta=2)
PSI4 = 0.414                                                # ∫ psi0^4 (Airy ground state)


def peeloff_profile(eta, Nx, ell, Y0=5.0, Yend=-4.0, dt=2e-3, seed=0):
    """Peel-off profile Y_node(x), x=0..Nx-1; noise white in Y, OU(ell)-correlated across x."""
    rng = np.random.default_rng(seed)
    nsteps = int((Y0 - Yend) / dt); sdt = np.sqrt(dt)
    a = np.exp(-1.0 / ell); b = np.sqrt(1.0 - a * a)
    # OU-correlated noise field Z[m, x]: white in m (the Y-sweep), stationary OU(ell) in x
    Z = np.empty((nsteps, Nx))
    Z[:, 0] = rng.standard_normal(nsteps)
    for x in range(1, Nx):
        Z[:, x] = a * Z[:, x - 1] + b * rng.standard_normal(nsteps)
    r0 = np.sqrt(Y0) + 1.0 / (4 * Y0)
    u = np.ones(Nx); w = np.full(Nx, r0); yn = np.full(Nx, np.nan); al = np.ones(Nx, bool)
    for m in range(nsteps):
        Y = Y0 - dt * m
        if not al.any():
            break
        un = u + w * dt
        wn = w + (Y * u) * dt - eta * u * sdt * Z[m]
        cr = al & (un * u <= 0.0)
        yn[cr] = Y - dt * (u[cr] / (u[cr] - un[cr]))        # sub-grid zero-crossing (no nugget)
        al[cr] = False
        u, w = un, wn
    yn[np.isnan(yn)] = Yend
    return yn


if __name__ == "__main__":
    eta = np.sqrt(2.0)                                      # beta = 4/eta^2 = 2
    Nx, ell, R = 600, 80.0, 40
    profiles = np.array([peeloff_profile(eta, Nx, ell, seed=100 + k) for k in range(R)])
    allv = profiles.ravel()

    # (i) one-point marginal
    m, sd = allv.mean(), allv.std()
    sk = ((allv - m) ** 3).mean() / sd ** 3
    print(f"(i) marginal  mean={m:.3f} [{TW2['mean']}]  std={sd:.3f} [{TW2['std']}]  "
          f"skew={sk:.3f} [{TW2['skew']}]  => one-point law = TW_2")

    # (ii) increment variance V(r) and (iii) covariance C(r)
    rmax = 200
    V = np.zeros(rmax + 1); C = np.zeros(rmax + 1)
    pc = profiles - profiles.mean(axis=1, keepdims=True)
    for r in range(0, rmax + 1):
        d = profiles[:, r:] - profiles[:, :Nx - r]
        V[r] = (d ** 2).mean()
        C[r] = (pc[:, r:] * pc[:, :Nx - r]).mean()
    # local roughness exponent (small r)
    rs = np.arange(2, 11)                                  # local (small-r) roughness window
    alpha, lnc = np.polyfit(np.log(rs), np.log(V[rs]), 1)
    slope_pred = 2 * eta ** 2 * PSI4 / ell                 # linear-response dV/dr
    print(f"(ii) roughness V(r) ~ r^alpha, alpha = {alpha:.2f}  (Airy/local-Brownian = 1; smooth = 2)")
    print(f"     small-r slope ~ {(V[2]-V[1]):.4f}/step   vs linear-response 2 eta^2 ∫psi0^4 /ell = {slope_pred:.4f}")
    print(f"     saturation V(inf) ~ {V[rmax]:.3f}   vs 2 Var(TW_2) = {2*TW2['std']**2:.3f}")
    print(f"(iii) covariance: C(0)={C[0]:.3f} [Var(TW_2)={TW2['std']**2:.3f}], "
          f"C decays to {C[rmax]:.3f}; correlation length ~ell={ell:.0f}")

    # ============================ figure ============================
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    C_P, C_V, C_C = "#2563eb", "#dc2626", "#16a34a"
    fig, ax = plt.subplots(1, 3, figsize=(15.2, 4.6))

    # Panel A: a few peel-off profiles (the Airy-like wandering process)
    xx = np.arange(Nx)
    for k, col in zip(range(4), ["#1d4ed8", "#0891b2", "#7c3aed", "#db2777"]):
        ax[0].plot(xx, profiles[k], color=col, lw=1.1, alpha=0.85)
    ax[0].axhline(TW2["mean"], color="#475569", lw=1.0, ls="--", label=r"$\mathbb{E}\,$TW$_2$")
    ax[0].set_xlabel(r"array site $x$ (coupling direction)")
    ax[0].set_ylabel(r"peel-off $Y_{\rm node}(x)$")
    ax[0].set_title(r"Peel-off profiles (candidate Airy$_2$ process)")
    ax[0].legend(fontsize=9, framealpha=0.95, edgecolor="#cbd5e1")

    # Panel B: one-point marginal vs TW_2 moments
    ax[1].hist(allv, bins=80, range=(-5, 1), density=True, color=C_P, alpha=0.55)
    ax[1].axvline(m, color="#0f172a", lw=1.2, ls="--")
    ax[1].set_xlabel(r"$Y_{\rm node}$"); ax[1].set_ylabel("density")
    ax[1].set_title(fr"1-point marginal $=$ TW$_2$ (skew {sk:.2f} vs 0.224)")
    ax[1].text(0.02, 0.95, f"mean {m:.2f} / {TW2['mean']}\nstd  {sd:.2f} / {TW2['std']}\n"
               f"skew {sk:.2f} / {TW2['skew']}", transform=ax[1].transAxes, va="top",
               fontsize=8.6, bbox=dict(fc="white", ec="#cbd5e1"))

    # Panel C: roughness V(r) (Hölder-1/2) + covariance
    r = np.arange(1, rmax + 1); rl = np.arange(2, 16)
    ax[2].loglog(r, V[1:], "o", color=C_V, ms=3.5, label=r"$V(r)=\mathrm{Var}(\Delta Y_{\rm node})$")
    ax[2].loglog(rl, np.exp(lnc) * rl ** alpha, "-", color="#0f172a", lw=2.0,
                 label=fr"local $\propto r^{{{alpha:.2f}}}$")
    ax[2].loglog(rl, np.exp(lnc) * rl ** 1.0 / rl[0] ** (1 - alpha) * 2 ** (1 - alpha), ":",
                 color=C_C, lw=1.6, label=r"Airy$_2$ Hölder-$\frac{1}{2}$ ($r^1$)")
    ax[2].loglog(r, np.full_like(r, 2 * sd ** 2, dtype=float), ":", color="#475569", lw=1.4,
                 label=r"$2\,\mathrm{Var}$ (saturation)")
    ax[2].loglog(rl, rl ** 2 * (V[2] / 4), "--", color="#94a3b8", lw=1.2, label=r"smooth $r^2$ (excluded)")
    ax[2].set_xlabel(r"lag $r$ (sites)"); ax[2].set_ylabel(r"increment variance $V(r)$")
    ax[2].set_title(fr"Rough field: $r^{{{alpha:.2f}}}$ (rougher than Airy$_2$ $r^1$)")
    ax[2].set_ylim(1e-3, 5); ax[2].legend(fontsize=8.2, framealpha=0.95, edgecolor="#cbd5e1")

    fig.tight_layout()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures",
                       "airy_process_peeloffs.png")
    out = os.path.abspath(out)
    fig.savefig(out, dpi=140)
    print("\nsaved", out)
