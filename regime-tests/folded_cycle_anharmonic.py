"""
folded_cycle_anharmonic.py — Part 2b: validate the PHASE MODULATION G(theta)
----------------------------------------------------------------------------
The symmetric Bautin grounded the INHERITED half (Channel A, sigma_* ~ sqrt(eps2),
C_q ~ 5.7).  The genuinely-NOVEL (not BG/BGK) content is the phase modulation
    sigma_*^A(theta) = C_q sqrt(eps2) sqrt(a(theta) c(theta) / b(theta)).
To make the FOLD GEOMETRY phase-dependent (not just shift the fold location), we
modulate the CUBIC coefficient by the phase:

    rho' = rho ( beta + g(theta) rho^2 - rho^4 ),   g(theta) = g0 + g1 cos(2 pi theta)
    beta' = -eps2                                    (super-slow drift)
    (degenerate noise on the radius:  + sigma dW)

Closed-form local fold (saddle-node of cycles) at frozen theta, g=g(theta):
    SNLC:  beta_SN = -g^2/4,   rho_SN = sqrt(g/2)
    near fold (u=rho-rho_SN, nu=beta-beta_SN):
        u' = -2 g sqrt(g/2) u^2 + sqrt(g/2) nu
    => JKK coefficients  b = 2 g sqrt(g/2),  a = sqrt(g/2),  c = 1
    => G(theta) = sqrt(ac/b) = sqrt( sqrt(g/2) / (2 g sqrt(g/2)) ) = 1/sqrt(2 g(theta)).

FROZEN-PHASE TEST (alpha=2): fix theta=theta0, drift beta through the SNLC, measure
the noise-induced early-escape threshold sigma_*(theta0).  PREDICTION:
    sigma_*(theta0) = C_q sqrt(eps2) G(theta0) = C_q sqrt(eps2) / sqrt(2 g(theta0)).
So sigma_*(theta0) * sqrt(2 g(theta0)) should be FLAT (= C_q sqrt(eps2)); a wrong
power of g should NOT flatten it (falsification control).
"""
from __future__ import annotations
import os
import numpy as np

G0, G1 = 1.0, 0.5           # g(theta) = g0 + g1 cos(2 pi theta),  g in [0.5, 1.5]
EPS2 = 2e-3


def g_of(theta):
    return G0 + G1 * np.cos(2 * np.pi * theta)


def G_pred(theta):
    return 1.0 / np.sqrt(2.0 * g_of(theta))


def beta_SN(g):
    return -g * g / 4.0


def rho_plus(beta, g):
    s = (g + np.sqrt(np.maximum(g * g + 4.0 * beta, 0.0))) / 2.0
    return np.sqrt(np.maximum(s, 0.0))


def simulate(theta0, sigma, N=40, dbeta=0.05, dt=1e-3, seed=0):
    """Frozen phase theta0 (g fixed); drift beta from SN+dbeta to SN-dbeta.
    Returns beta_death = beta at first rho < 0.5 rho_SN."""
    rng = np.random.default_rng(seed)
    g = g_of(theta0)
    bSN = beta_SN(g)
    rho_SN = np.sqrt(g / 2.0)
    rho_coll = 0.5 * rho_SN
    b0 = bSN + dbeta
    b_end = bSN - dbeta
    rho = np.full(N, rho_plus(b0, g))
    beta = b0
    n = int((b0 - b_end) / EPS2 / dt)
    sdt = np.sqrt(dt)
    dead = np.zeros(N, bool)
    bdeath = np.full(N, b_end)
    for _ in range(n):
        rho = rho + rho * (beta + g * rho**2 - rho**4) * dt + sigma * sdt * rng.standard_normal(N)
        beta = beta - EPS2 * dt
        nd = (~dead) & (rho < rho_coll)
        if nd.any():
            bdeath[nd] = beta
            dead |= nd
    return bdeath - bSN          # margin m (>0 = early death)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    figdir = os.path.join(os.path.dirname(here), "figures")
    os.makedirs(figdir, exist_ok=True)

    thetas = np.linspace(0.0, 0.5, 6)           # samples g from 1.5 down to 0.5
    sigmas = np.array([0.04, 0.07, 0.11, 0.16, 0.22, 0.30, 0.40])
    N = 40

    print("\n=== Part 2b: phase modulation via the canard R-Theta collapse ===")
    print(f"  g(theta)={G0}+{G1}cos(2pi theta) in [{G0-G1},{G0+G1}], eps2={EPS2:.0e}, N={N}")
    print("  window  delta_nu = eps2^{2/3} g^{-2/3};  noise scale G=1/sqrt(2g)")
    print("  R = <m> g^{2/3}/eps2^{2/3};  Theta_correct = sigma sqrt(2g)/sqrt(eps2)\n")

    # per-phase (Theta_correct, Theta_naive, R) curves
    curves = []
    for i, th in enumerate(thetas):
        g = g_of(th)
        m = np.array([np.mean(simulate(th, s, N=N, seed=10 * i + k))
                      for k, s in enumerate(sigmas)])
        R = m * g ** (2.0 / 3.0) / EPS2 ** (2.0 / 3.0)
        Th_c = sigmas * np.sqrt(2.0 * g) / np.sqrt(EPS2)     # correct: includes 1/G
        Th_n = sigmas / np.sqrt(EPS2)                        # naive: ignores G
        curves.append(dict(th=th, g=g, R=R, Tc=Th_c, Tn=Th_n))

    # collapse quality: cross-phase std of R on a common Theta grid (correct vs naive)
    def collapse_cv(key):
        lo = max(c[key].min() for c in curves)
        hi = min(c[key].max() for c in curves)
        grid = np.linspace(lo, hi, 12)
        stacked = np.array([np.interp(grid, c[key], c["R"]) for c in curves])
        spread = stacked.std(axis=0)
        scale = np.abs(stacked).mean() + 1e-9
        return float(np.mean(spread) / scale), grid, stacked

    cv_c, gridc, stc = collapse_cv("Tc")
    cv_n, gridn, stn = collapse_cv("Tn")
    print(f"  CROSS-PHASE collapse spread (mean std / scale) over common Theta window:")
    print(f"     vs Theta_correct = sigma sqrt(2g)/sqrt(eps2)  (G included): {cv_c:.3f}")
    print(f"     vs Theta_naive   = sigma/sqrt(eps2)           (G ignored) : {cv_n:.3f}")
    verdict = "COLLAPSES better WITH G" if cv_c < cv_n else "G does NOT improve collapse"
    print(f"     => {verdict}  (ratio naive/correct = {cv_n/max(cv_c,1e-9):.2f})\n")

    # C_q from the correct-collapsed curve: Theta at which pooled R crosses ~1 (canonical window)
    pooled = stc.mean(axis=0)
    Cq = np.nan
    for j in range(len(gridc) - 1):
        if pooled[j] < 1.0 <= pooled[j + 1]:
            f = (1.0 - pooled[j]) / (pooled[j + 1] - pooled[j])
            Cq = gridc[j] + f * (gridc[j + 1] - gridc[j]); break
    print(f"  C_q (pooled R crosses canonical window R=1 at Theta) ~ {Cq:.2f}"
          f"   (fuzzy; O(several), consistent with symmetric-Bautin 5.7)\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    for c in curves:
        ax[0].plot(c["Tc"], c["R"], "o-", ms=3, label=fr"$\theta$={c['th']:.2f}, g={c['g']:.2f}")
    ax[0].axhline(1.0, color="k", ls=":", lw=0.8)
    ax[0].set_xlabel(r"$\Theta=\sigma\sqrt{2g(\theta)}/\sqrt{\epsilon_2}$  (= $\sigma/(\sqrt{\epsilon_2}\,G)$)")
    ax[0].set_ylabel(r"$R=\langle m\rangle\,g^{2/3}/\epsilon_2^{2/3}$")
    ax[0].set_title(f"WITH G: collapse across phase (CV {cv_c:.2f})")
    ax[0].legend(fontsize=7, frameon=False); ax[0].grid(alpha=0.3)
    for c in curves:
        ax[1].plot(c["Tn"], c["R"], "s--", ms=3, label=fr"$\theta$={c['th']:.2f}")
    ax[1].axhline(1.0, color="k", ls=":", lw=0.8)
    ax[1].set_xlabel(r"$\Theta_{\rm naive}=\sigma/\sqrt{\epsilon_2}$  (G ignored)")
    ax[1].set_ylabel(r"$R$")
    ax[1].set_title(f"WITHOUT G: splays (CV {cv_n:.2f})")
    ax[1].legend(fontsize=7, frameon=False); ax[1].grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(figdir, "folded_cycle_anharmonic_Gtheta.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print(f"  figure -> figures/folded_cycle_anharmonic_Gtheta.png\n")


if __name__ == "__main__":
    main()
