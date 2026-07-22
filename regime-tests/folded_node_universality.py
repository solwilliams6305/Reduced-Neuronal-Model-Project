"""
folded_node_universality.py
===========================
Direction D: does the noisy FOLDED NODE obey the same Tracy-Widom_beta edge law as the
folded cycle?  If so, "the folded cycle is in the TW edge class" becomes "this FAMILY is."

The folded node's canard inner equation is the WEBER (parabolic-cylinder) equation, NOT the
Airy equation of the fold -- a genuinely different global operator,
        u'' = ( z^2/4 - (nu+1/2) - eta xi ) u ,
with an oscillatory well |z| < z*, z* = 2 sqrt(nu+1/2) (the rotations / secondary canards),
and a simple turning point at z* with local slope s = sqrt(nu+1/2).  Near z* the potential is
locally LINEAR (Q ~ s(z-z*)), so by soft-edge universality the noise-broadened peel-off should
be TW_beta with the SAME local rescaling as the fold: zeta = s^{1/3}(z-z*), beta = 4 s / eta^2.

We reuse the validated machinery: recessive Cole-Hopf entry (here the recessive parabolic-
cylinder D_nu, ratio -D'/D = z/2 - nu/z), first node of u = peel-off, moments vs TW_beta.

What this script establishes (see printout):
  * the deterministic peel-off (eta->0) is Weber-SHIFTED from the Airy zero -2.338, and
    -> -2.338 as the fold sharpens (s -> inf): location is model-dependent;
  * the FLUCTUATION SHAPE converges to TW_beta as s -> inf -- skew/kurt -> TW_beta, std -> TW_beta,
    with controlled finite-curvature (Weber) corrections at moderate s;
  * at a sharp fold (large s) the folded node matches TW_beta across the noise family
    (beta = 4,2,1), coinciding with the folded cycle (Airy) -- i.e. SAME edge universality class.
"""
from __future__ import annotations
import os
import numpy as np

# exact TW_beta cumulants (std, skew, excess-kurt); mean for beta=1,2
TW = {1: dict(mean=-1.2065, std=1.2680, skew=0.2935, kurt=0.1652),
      2: dict(mean=-1.7711, std=0.9018, skew=0.2241, kurt=0.0935),
      4: dict(mean=None,    std=0.7109, skew=0.1655, kurt=0.0491)}


def first_node(eta, N, Qfun, Y0, Yend, r0, dt=2e-3, seed=0):
    """Recessive Cole-Hopf sweep u''=(Q(Y)-eta xi)u; return first-node levels (escaped)."""
    rng = np.random.default_rng(seed)
    u = np.ones(N); w = np.full(N, r0)
    yn = np.full(N, np.nan); al = np.ones(N, bool)
    nsteps = int((Y0 - Yend) / dt); sdt = np.sqrt(dt)
    for k in range(nsteps):
        Y = Y0 - dt * k
        if not al.any():
            break
        un = u + w * dt
        wn = w + (Qfun(Y) * u) * dt - eta * u * sdt * rng.standard_normal(N)
        cr = al & (un * u <= 0.0); yn[cr] = Y; al[cr] = False
        u, w = un, wn
    return yn[~np.isnan(yn)]


def moments(x):
    m, sd = x.mean(), x.std(); c = x - m
    return m, sd, float((c**3).mean() / sd**3), float((c**4).mean() / sd**4 - 3)


def weber_peeloff(eta, N, nu, seed=0):
    s = np.sqrt(nu + 0.5); zstar = 2 * s
    z0 = zstar + 2.5; zend = zstar - (6.0 * s**(-1.0/3.0) + 0.5)
    r0 = zstar / 2 - nu / zstar
    zn = first_node(eta, N, lambda z: z*z/4 - (nu + 0.5), z0, zend, r0, seed=seed)
    return (zn - zstar) * s**(1.0/3.0)            # local Airy coordinate zeta


def airy_peeloff(eta, N, seed=0):                 # the fold / folded cycle (linear potential)
    Y0 = 6.0; r0 = np.sqrt(Y0) + 1/(4*Y0)
    return first_node(eta, N, lambda Y: Y, Y0, -7.0, r0, seed=seed)


if __name__ == "__main__":
    N = 30000

    # ---- 1. convergence to the TW edge as the fold sharpens (fixed beta = 2) ----
    nus = np.array([1.5, 3.5, 8.0, 20.0, 50.0, 120.0])
    ss = np.sqrt(nus + 0.5)
    conv = []   # (s, skew, kurt, std, anchor_mean)
    print("Convergence at beta=2 (eta = sqrt(2s)) -- folded node -> TW_2 as s grows:")
    print(f"{'s':>6} {'skew':>8} {'kurt':>8} {'std':>7} {'anchor(eta->0)':>15}")
    for nu, s in zip(nus, ss):
        eta = np.sqrt(2 * s)                       # beta = 4s/eta^2 = 2
        z = weber_peeloff(eta, N, nu, seed=int(10 + nu))
        _, sd, sk, ku = moments(z)
        z0 = weber_peeloff(0.10, 6000, nu, seed=99)   # near-deterministic anchor
        anc = z0.mean()
        conv.append((s, sk, ku, sd, anc))
        print(f"{s:6.2f} {sk:8.3f} {ku:8.3f} {sd:7.3f} {anc:15.3f}")
    conv = np.array(conv)
    print(f"   TW_2 targets: skew {TW[2]['skew']}, kurt {TW[2]['kurt']}, std {TW[2]['std']};"
          f" Airy zero -2.338")

    # ---- 2. beta-family at a sharp fold (nu=50) vs the folded cycle (Airy) ----
    betas = np.array([4.0, 2.0, 1.0])
    nu_sharp = 120.0; s_sharp = np.sqrt(nu_sharp + 0.5)
    fn_sk, fn_ku, ai_sk, ai_ku = [], [], [], []
    print(f"\nbeta-family at sharp fold nu={nu_sharp} (s={s_sharp:.2f}) vs folded cycle (Airy):")
    print(f"{'beta':>5} {'node skew':>10} {'cycle skew':>11} {'TW skew':>8} | "
          f"{'node kurt':>10} {'cycle kurt':>11} {'TW kurt':>8}")
    for b in betas:
        eta_n = np.sqrt(4 * s_sharp / b)
        zn = weber_peeloff(eta_n, N, nu_sharp, seed=int(200 + b))
        _, _, skn, kun = moments(zn)
        eta_a = np.sqrt(4.0 / b)                   # Airy: s=1, beta=4/eta^2
        za = airy_peeloff(eta_a, N, seed=int(300 + b))
        _, _, ska, kua = moments(za)
        fn_sk.append(skn); fn_ku.append(kun); ai_sk.append(ska); ai_ku.append(kua)
        print(f"{b:5.0f} {skn:10.3f} {ska:11.3f} {TW[b]['skew']:8.3f} | "
              f"{kun:10.3f} {kua:11.3f} {TW[b]['kurt']:8.3f}")

    # full-moment match at beta=2, sharp fold
    z2 = weber_peeloff(np.sqrt(2 * s_sharp), 60000, nu_sharp, seed=7)
    m2, sd2, sk2, ku2 = moments(z2)
    print(f"\n  sharp folded node, beta=2: (mean,std,skew,kurt) = "
          f"({m2:.3f},{sd2:.3f},{sk2:.3f},{ku2:.3f})  vs TW_2 "
          f"({TW[2]['mean']},{TW[2]['std']},{TW[2]['skew']},{TW[2]['kurt']})")

    # ---- figure ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    C_SK, C_KU, C_TW, C_AI = "#2563eb", "#ea8a0b", "#16a34a", "#dc2626"

    fig, ax = plt.subplots(1, 2, figsize=(12.4, 5.0))

    # Panel A: convergence vs s
    ax[0].plot(conv[:, 0], conv[:, 1], "o-", color=C_SK, ms=7, mec="white", mew=1.0,
               label="folded-node skew")
    ax[0].plot(conv[:, 0], conv[:, 2], "s-", color=C_KU, ms=7, mec="white", mew=1.0,
               label="folded-node excess kurt")
    ax[0].axhline(TW[2]["skew"], ls="--", color=C_SK, lw=1.5)
    ax[0].axhline(TW[2]["kurt"], ls="--", color=C_KU, lw=1.5)
    ax[0].axhline(0.0, color="#cbd5e1", lw=0.8)
    ax[0].text(conv[-1, 0], TW[2]["skew"] + 0.012, "TW$_2$ skew 0.224", color=C_SK,
               fontsize=8.5, ha="right")
    ax[0].text(conv[-1, 0], TW[2]["kurt"] - 0.03, "TW$_2$ kurt 0.093", color=C_KU,
               fontsize=8.5, ha="right")
    ax[0].set_xlabel(r"turning-point slope  $s=\sqrt{\nu+1/2}$  (fold sharpens $\rightarrow$)")
    ax[0].set_ylabel("standardised cumulant")
    ax[0].set_title(r"Folded node $\to$ TW$_2$ edge as the fold sharpens ($\beta=2$)")
    ax[0].grid(True, color="#eef2f7", lw=0.6)
    ax[0].legend(loc="lower right", fontsize=9, framealpha=0.95, edgecolor="#cbd5e1")

    # Panel B: beta-family at sharp fold
    bx = betas
    tw_sk = [TW[b]["skew"] for b in betas]
    ax[1].plot(bx, tw_sk, "-", color=C_TW, lw=2.0, zorder=2, label=r"exact TW$_\beta$ skew")
    ax[1].plot(bx, fn_sk, "o", color=C_SK, ms=10, mec="white", mew=1.1, zorder=5,
               label=fr"folded node ($\nu={int(nu_sharp)}$, sharp)")
    ax[1].plot(bx, ai_sk, "*", color=C_AI, ms=15, mec="white", mew=1.0, zorder=4,
               label="folded cycle (Airy)")
    for b, v in zip(bx, tw_sk):
        ax[1].annotate(f"{v:.3f}", (b, v), textcoords="offset points", xytext=(0, 9),
                       fontsize=8, color=C_TW, ha="center")
    ax[1].set_xlabel(r"$\beta = 4s/\eta^2$")
    ax[1].set_ylabel("peel-off skewness")
    ax[1].set_xticks([1, 2, 4])
    ax[1].set_title("Sharp folded node matches TW$_\\beta$ across the noise family")
    ax[1].grid(True, color="#eef2f7", lw=0.6)
    ax[1].legend(loc="upper right", fontsize=9, framealpha=0.95, edgecolor="#cbd5e1")

    fig.tight_layout()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures",
                       "folded_node_universality.png")
    out = os.path.abspath(out)
    fig.savefig(out, dpi=140)
    print("\nsaved", out)
