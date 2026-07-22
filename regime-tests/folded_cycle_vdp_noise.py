"""
folded_cycle_vdp_noise.py — Part 2, step 2 (the C_q / scaling measurement)
--------------------------------------------------------------------------
Noisy three-timescale forced VdP:
    dx = (x - x^3/3 - w) dt + sigma dW    (fast, degenerate noise)
    dw = eps1 (x - a) dt                  (slow, relaxation cycle)
    da = eps2 dt                          (super-slow drift through the fold a=1)

Channel-A escape: as a drifts up toward the fold a=1, noise kills the cycle
EARLY (before a=1).  Observable: a_death = a at the last full bottom excursion
(x < -1.5); early-death margin m = 1 - a_death (>0 means noise killed it before
the fold).

THE DISCRIMINATING TEST.  Folded-cycle theory predicts eta = sigma/sqrt(eps2),
hence sigma_* proportional to sqrt(eps2) (the SUPER-SLOW drift), and INDEPENDENT
of eps1 (the cycle's own slow rate) at leading order.  The plain canard has only
one eps.  So:
  - sweep eps2 (fixed eps1): sigma_*(eps2) should scale ~ eps2^{1/2};
  - sweep eps1 (fixed eps2): sigma_*(eps1) should be ~ flat (slope ~0);
  - pooled collapse of <m>(sigma) vs Theta = sigma/sqrt(eps2) (collapses) vs
    sigma/sqrt(eps1) (does NOT)  -> falsification control.
Then C_q ~ sigma_* / (sqrt(eps2) * G) with G = sqrt(ac/b) ~ O(1) (expect ~8).
"""
from __future__ import annotations
import os
import numpy as np


def simulate(eps1, eps2, sigma, N=40, a0=0.91, a_end=1.01, dt=1e-3, seed=0):
    """Vectorised Euler-Maruyama over N trajectories. Returns a_death array."""
    rng = np.random.default_rng(seed)
    T = (a_end - a0) / eps2
    n = int(T / dt)
    x = np.full(N, 0.5); w = np.zeros(N)
    a = a0
    sdt = np.sqrt(dt)
    a_death = np.full(N, a0)          # last a with a full bottom excursion
    for _ in range(n):
        x = x + (x - x**3 / 3.0 - w) * dt + sigma * sdt * rng.standard_normal(N)
        w = w + eps1 * (x - a) * dt
        a = a + eps2 * dt
        low = x < -1.5
        if low.any():
            a_death[low] = a
    return a_death


def cell_curve(eps1, eps2, sigmas, N, seed0=0):
    """<m> = <1 - a_death> vs sigma for one (eps1,eps2) cell."""
    m_mean = np.empty(len(sigmas)); m_se = np.empty(len(sigmas))
    for k, s in enumerate(sigmas):
        ad = simulate(eps1, eps2, s, N=N, seed=seed0 + k)
        m = 1.0 - ad
        m_mean[k] = m.mean(); m_se[k] = m.std() / np.sqrt(N)
    return m_mean, m_se


def sigma_star(sigmas, m_mean, m_thresh=0.02):
    """interpolate sigma at which <m> first crosses m_thresh."""
    for i in range(len(sigmas) - 1):
        if m_mean[i] < m_thresh <= m_mean[i + 1]:
            f = (m_thresh - m_mean[i]) / (m_mean[i + 1] - m_mean[i])
            return sigmas[i] * (sigmas[i + 1] / sigmas[i]) ** f
    return np.nan


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    figdir = os.path.join(os.path.dirname(here), "figures")
    os.makedirs(figdir, exist_ok=True)

    sigmas = np.array([0.05, 0.09, 0.14, 0.20, 0.30, 0.44])
    N = 40

    # cells: eps2-sweep at eps1=0.1, and eps1-sweep at eps2=2e-3
    eps2_sweep = [(0.10, 1e-3), (0.10, 2e-3), (0.10, 4e-3)]
    eps1_sweep = [(0.05, 2e-3), (0.20, 2e-3)]            # (0.10,2e-3) shared
    cells = eps2_sweep + eps1_sweep

    print("\n=== Part 2 step 2: folded-VdP Channel-A sigma_* scaling ===")
    print(f"sigma grid: {sigmas}\nN={N} trajectories/cell\n")
    results = {}
    for j, (e1, e2) in enumerate(cells):
        mm, se = cell_curve(e1, e2, sigmas, N, seed0=100 * j)
        ss = sigma_star(sigmas, mm)
        results[(e1, e2)] = dict(m=mm, se=se, sstar=ss)
        print(f"  eps1={e1:.2f} eps2={e2:.0e}:  sigma_*={ss:.3f}   "
              f"<m>(sigma)=[" + ", ".join(f"{v:.3f}" for v in mm) + "]")

    # --- scaling checks ---
    print("\n-- eps2 scaling at eps1=0.10 (predict slope 1/2) --")
    e2s = np.array([1e-3, 2e-3, 4e-3])
    ss_e2 = np.array([results[(0.10, e)]["sstar"] for e in e2s])
    p2 = np.polyfit(np.log(e2s), np.log(ss_e2), 1)[0]
    for e, s in zip(e2s, ss_e2):
        print(f"   eps2={e:.0e}: sigma_*={s:.3f}   sigma_*/sqrt(eps2)={s/np.sqrt(e):.2f}")
    print(f"   fitted slope d ln sigma_* / d ln eps2 = {p2:.3f}  (predict 0.50)")

    print("\n-- eps1 scaling at eps2=2e-3 (predict slope ~0) --")
    e1s = np.array([0.05, 0.10, 0.20])
    ss_e1 = np.array([results[(e, 2e-3)]["sstar"] for e in e1s])
    p1 = np.polyfit(np.log(e1s), np.log(ss_e1), 1)[0]
    for e, s in zip(e1s, ss_e1):
        print(f"   eps1={e:.2f}: sigma_*={s:.3f}")
    print(f"   fitted slope d ln sigma_* / d ln eps1 = {p1:.3f}  (predict ~0)")

    # C_q estimate (G = sqrt(ac/b) ~ O(1); report sigma_*/sqrt(eps2))
    Cq_proxy = np.mean([results[(0.10, e)]["sstar"] / np.sqrt(e) for e in e2s])
    print(f"\n-- C_q proxy = mean sigma_*/sqrt(eps2) (eps1=0.10) = {Cq_proxy:.2f}")
    print("   (true C_q = this / G, G=sqrt(ac/b)~O(1); expect O(8) as full FHN/VdP)\n")

    # --- figure: scaling + collapse ---
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))

    ax[0].loglog(e2s, ss_e2, "o-", color="C0", label=f"vary $\\epsilon_2$ (slope {p2:.2f})")
    ax[0].loglog(e2s, ss_e2[1] * (e2s / e2s[1]) ** 0.5, "k--", lw=1, label="slope 1/2 ref")
    ax[0].loglog(e1s, ss_e1, "s-", color="C3", label=f"vary $\\epsilon_1$ (slope {p1:.2f})")
    ax[0].set_xlabel(r"$\epsilon$ (drift $\epsilon_2$ / cycle $\epsilon_1$)")
    ax[0].set_ylabel(r"$\sigma_*$")
    ax[0].set_title(r"$\sigma_*\propto\sqrt{\epsilon_2}$, flat in $\epsilon_1$")
    ax[0].legend(fontsize=8, frameon=False); ax[0].grid(alpha=0.3, which="both")

    for (e1, e2), r in results.items():
        ax[1].plot(sigmas / np.sqrt(e2), r["m"], "o-", ms=3,
                   label=f"$\\epsilon_1$={e1}, $\\epsilon_2$={e2:.0e}")
    ax[1].axhline(0.02, color="k", ls=":", lw=0.8)
    ax[1].set_xlabel(r"$\Theta=\sigma/\sqrt{\epsilon_2}$")
    ax[1].set_ylabel(r"$\langle m\rangle=\langle 1-a_{\rm death}\rangle$")
    ax[1].set_title(r"Collapse vs $\sigma/\sqrt{\epsilon_2}$ (folded-cycle scale)")
    ax[1].legend(fontsize=7, frameon=False); ax[1].grid(alpha=0.3)

    fig.tight_layout()
    out = os.path.join(figdir, "folded_cycle_vdp_sigmastar.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print(f"  figure -> figures/folded_cycle_vdp_sigmastar.png\n")


if __name__ == "__main__":
    main()
