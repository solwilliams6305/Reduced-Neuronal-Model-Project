"""
folded_cycle_bautin.py — Part 2 (faithful): fold of limit cycles + noise
------------------------------------------------------------------------
The VdP-through-Hopf testbed FAILED (amplitude->0 makes it a resonator; noise
induces spikes past the fold rather than killing the cycle early). The faithful
JKK object is a fold of limit cycles at FINITE amplitude: a saddle-node of
cycles (SNLC), realized by the Bautin / generalized-Hopf normal form.

Planar fast oscillator (degenerate noise on the fast Cartesian variable u):
    du = [ u (beta + R^2 - R^4) - omega v ] dt + sigma dW,   R^2 = u^2+v^2
    dv = [ v (beta + R^2 - R^4) + omega u ] dt
    dbeta = -eps2 dt                                    (super-slow drift)

For beta in (-1/4, 0): a stable FP at rho=0, an UNSTABLE cycle rho_- and a
STABLE cycle rho_+.  They annihilate at the SNLC  beta_SN = -1/4, rho_SN = 1/sqrt2.
Channel-A escape: noise knocks the trajectory off the stable cycle rho_+ across
the unstable-cycle barrier rho_- into the FP -- the cycle dies EARLY (at beta >
beta_SN).  Finite-amplitude fold => no Hopf resonator contamination.

Local fold normal form (amplitude u = rho - rho_SN, nu = beta - beta_SN):
    u' = -sqrt(2) u^2 + (1/sqrt2) nu,   nu' = -eps2
=> JKK coefficients  a = 1/sqrt2, b = sqrt2, c = 1  =>  G = sqrt(ac/b) = 1/sqrt2.
Predicted  sigma_*^A = C_q sqrt(eps2) G,  G = 0.707, expect C_q = O(8).

TEST: sigma_* proportional to sqrt(eps2) (the drift), with C_q = O(8) inheriting
the canard/VdP constant.  (omega = phase rate decouples in the symmetric model,
so eps1-independence is structural here; phase modulation G(theta) needs the
anharmonic version -> Channel B / Part 3.)
"""
from __future__ import annotations
import os
import numpy as np

BETA_SN = -0.25
RHO_SN = 1.0 / np.sqrt(2.0)
G_PRED = np.sqrt((1/np.sqrt(2)) * 1.0 / np.sqrt(2))   # sqrt(ac/b) = 1/sqrt2


def rho_plus(beta):
    """stable-cycle radius: rho^2 = (1 + sqrt(1+4 beta))/2."""
    return np.sqrt((1.0 + np.sqrt(1.0 + 4.0 * beta)) / 2.0)


def simulate(eps2, sigma, omega=1.0, N=40, beta0=-0.18, beta_end=-0.30,
             dt=1e-3, rho_collapse=0.40, seed=0):
    """Vectorised EM; returns beta_death (beta at first rho < rho_collapse)."""
    rng = np.random.default_rng(seed)
    n = int((beta0 - beta_end) / eps2 / dt)
    r0 = rho_plus(beta0)
    u = np.full(N, r0); v = np.zeros(N)
    beta = beta0
    sdt = np.sqrt(dt)
    dead = np.zeros(N, bool)
    beta_death = np.full(N, beta_end)
    for _ in range(n):
        R2 = u * u + v * v
        g = beta + R2 - R2 * R2
        u = u + (u * g - omega * v) * dt + sigma * sdt * rng.standard_normal(N)
        v = v + (v * g + omega * u) * dt
        beta = beta - eps2 * dt
        newdead = (~dead) & (np.sqrt(u * u + v * v) < rho_collapse)
        if newdead.any():
            beta_death[newdead] = beta
            dead |= newdead
    return beta_death


def cell_curve(eps2, sigmas, omega, N, seed0):
    m_mean = np.empty(len(sigmas))
    for k, s in enumerate(sigmas):
        bd = simulate(eps2, s, omega=omega, N=N, seed=seed0 + k)
        m_mean[k] = np.mean(bd - BETA_SN)        # early-death margin (>0)
    return m_mean


def sigma_star(sigmas, m, thr):
    for i in range(len(sigmas) - 1):
        if m[i] < thr <= m[i + 1]:
            f = (thr - m[i]) / (m[i + 1] - m[i])
            return sigmas[i] * (sigmas[i + 1] / sigmas[i]) ** f
    return np.nan


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    figdir = os.path.join(os.path.dirname(here), "figures")
    os.makedirs(figdir, exist_ok=True)

    sigmas = np.array([0.03, 0.06, 0.10, 0.16, 0.24, 0.36, 0.50])
    N = 50
    eps2_list = [1e-3, 2e-3, 4e-3]
    thr = 0.01

    print("\n=== Part 2 (faithful): Bautin fold-of-cycles, Channel-A sigma_* ===")
    print(f"  SNLC at beta=-0.25, rho_SN={RHO_SN:.3f}; G=sqrt(ac/b)={G_PRED:.3f}")
    print(f"  sigma grid {sigmas}, N={N}\n")

    curves = {}
    for j, e2 in enumerate(eps2_list):
        m = cell_curve(e2, sigmas, omega=1.0, N=N, seed0=1000 * j)
        curves[e2] = m
        ss = sigma_star(sigmas, m, thr)
        print(f"  eps2={e2:.0e}: sigma_*={ss:.3f}  sigma_*/sqrt(eps2)={ss/np.sqrt(e2):6.2f}"
              f"  C_q={ss/np.sqrt(e2)/G_PRED:5.2f}")
        print(f"     <m>=[" + ", ".join(f"{x:.4f}" for x in m) + "]")

    e2a = np.array(eps2_list)
    ssa = np.array([sigma_star(sigmas, curves[e], thr) for e in e2a])
    good = ~np.isnan(ssa)
    slope = np.polyfit(np.log(e2a[good]), np.log(ssa[good]), 1)[0] if good.sum() >= 2 else np.nan
    Cq = np.nanmean(ssa / np.sqrt(e2a) / G_PRED)
    print(f"\n  fitted slope d ln sigma_*/d ln eps2 = {slope:.3f}   (predict 0.50)")
    print(f"  model C_q = mean sigma_*/(sqrt(eps2) G) = {Cq:.2f}   (expect O(8))\n")

    # quick omega-independence spot check (structural: phase decouples)
    print("  omega spot-check (eps2=2e-3, sigma=0.16):")
    for om in (0.5, 1.0, 2.0):
        bd = simulate(2e-3, 0.16, omega=om, N=80, seed=55)
        print(f"     omega={om:.1f}: <m>={np.mean(bd-BETA_SN):.4f}")

    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    ax[0].loglog(e2a, ssa, "o-", color="C0", label=f"measured (slope {slope:.2f})")
    ax[0].loglog(e2a, ssa[1] * (e2a / e2a[1]) ** 0.5, "k--", lw=1, label="slope 1/2")
    ax[0].set_xlabel(r"$\epsilon_2$ (super-slow drift)")
    ax[0].set_ylabel(r"$\sigma_*$")
    ax[0].set_title(rf"$\sigma_*\propto\sqrt{{\epsilon_2}}$,  $C_q\approx{Cq:.1f}$")
    ax[0].legend(fontsize=9, frameon=False); ax[0].grid(alpha=0.3, which="both")
    for e2 in eps2_list:
        ax[1].plot(sigmas / np.sqrt(e2), curves[e2], "o-", ms=4,
                   label=fr"$\epsilon_2$={e2:.0e}")
    ax[1].axhline(thr, color="k", ls=":", lw=0.8)
    ax[1].set_xlabel(r"$\Theta=\sigma/\sqrt{\epsilon_2}$")
    ax[1].set_ylabel(r"$\langle m\rangle=\langle\beta_{\rm death}-\beta_{\rm SN}\rangle$")
    ax[1].set_title(r"Collapse vs $\sigma/\sqrt{\epsilon_2}$ (fold of cycles)")
    ax[1].legend(fontsize=8, frameon=False); ax[1].grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(figdir, "folded_cycle_bautin_sigmastar.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print(f"\n  figure -> figures/folded_cycle_bautin_sigmastar.png\n")


if __name__ == "__main__":
    main()
