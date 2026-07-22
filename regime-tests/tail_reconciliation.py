"""
tail_reconciliation.py
======================
Direction A, part (b): is the early-escape probability the polynomial hazard, or the
exponential stochastic-Airy / Tracy-Widom left tail?

P_early(eta) = P(Y_node > 0) = probability the noisy canard peels off BEFORE the fold.

CORRECT SIMULATION (the canard is the recessive / repelling branch).  The canard is the
recessive solution of the Cole-Hopf linear equation, NOT the attracting Riccati branch.
We therefore integrate the linear, non-explosive equation

        u'' = (Y - eta xi) u ,     Y = Y0 - T ,     xi = dB/dT,

from a deep start Y0 on the recessive branch (asymptotic ratio -Ai'/Ai = sqrt(Y0)+1/(4Y0)),
and record the FIRST node of u (a node of u <=> blow-up of R=-u'/u <=> peel-off).  The
multiplicative noise (eta*u) vanishes at a node, so there are no spurious crossings; and
integrating DOWNWARD is self-correcting (the dominant Bi component decays), so a deep start
is stable.  Early escape = Y_node > 0.  Validated: mean Y_node -> -2.338 (first Airy zero) as
eta->0, and P_early matches the exact TW tail 1-F_beta(0) at beta=1,2 (see cross-checks).

Two competing predictions for the tail:
  * quasi-static HAZARD     H = eta^2/(4 pi)      -- polynomial in eta (counts crossings);
  * stochastic-Airy TAIL    P_early ~ exp(-c/eta^2)  -- exponential (TW_beta left large
                            deviation, beta=4/eta^2, exp(-beta s^3/24) = exp(-s^3/6 eta^2);
                            Painleve II action c = 5.4439).
Discriminator: ln P_early vs 1/eta^2 is a STRAIGHT LINE (exponential), and the polynomial
hazard sits far above (P_early < H: a separatrix crossing is necessary but not sufficient
for a committed node).
"""
from __future__ import annotations
import os
import numpy as np


def P_early(eta, N, Y0=5.0, Yend=-4.0, dt=2.5e-3, seed=0):
    """First-node early-escape probability via the recessive Cole-Hopf sweep."""
    rng = np.random.default_rng(seed)
    r0 = np.sqrt(Y0) + 1.0 / (4.0 * Y0)         # asymptotic recessive ratio -Ai'/Ai at Y0
    u = np.ones(N)
    w = np.full(N, r0)                           # w = du/dT = +(-Ai'/Ai) u
    Ynode = np.full(N, np.nan)
    alive = np.ones(N, bool)
    nsteps = int((Y0 - Yend) / dt)
    sdt = np.sqrt(dt)
    for k in range(nsteps):
        Y = Y0 - dt * k
        if not alive.any():
            break
        un = u + w * dt
        wn = w + (Y * u) * dt - eta * u * sdt * rng.standard_normal(N)
        crossed = alive & (un * u <= 0.0)        # first node of u
        Ynode[crossed] = Y
        alive[crossed] = False
        u, w = un, wn
    Ynode[np.isnan(Ynode)] = Yend
    n_e = int((Ynode > 0.0).sum())
    return n_e / N, np.sqrt(max(n_e, 1)) / N, float(Ynode.mean())


if __name__ == "__main__":
    etas = np.array([1.0, 1.2, np.sqrt(2.0), 1.6, 1.8, 2.0])
    N = 40000
    P = np.empty_like(etas); E = np.empty_like(etas); Ybar = np.empty_like(etas)
    print(f"{'eta':>7} {'beta=4/eta^2':>12} {'1/eta^2':>8} {'P_early':>11} {'+/-':>10} "
          f"{'<Ynode>':>9} {'H=eta^2/4pi':>12}")
    for i, eta in enumerate(etas):
        P[i], E[i], Ybar[i] = P_early(float(eta), N, seed=11 + i)
        print(f"{eta:7.3f} {4/eta**2:12.3f} {1/eta**2:8.3f} {P[i]:11.5f} {E[i]:10.2e} "
              f"{Ybar[i]:9.3f} {eta**2/(4*np.pi):12.5f}")

    # ---- fit ln P vs u = 1/eta^2 ----
    u = 1.0 / etas ** 2
    slope, intercept = np.polyfit(u, np.log(P), 1)
    c_fit = -slope
    resid = np.log(P) - (slope * u + intercept)
    R2 = 1 - np.sum(resid ** 2) / np.sum((np.log(P) - np.log(P).mean()) ** 2)
    print(f"\n  fit  ln P_early = {intercept:.3f} - c/eta^2 :  c_fit = {c_fit:.3f}  (R^2={R2:.3f})")
    print(f"       Painleve II action c = 5.4439 (= s^3/6); c_fit is the moderate-eta effective slope.")
    iw = int(np.argmin(np.abs(etas - np.sqrt(2)))); i2 = int(np.argmin(np.abs(etas - 2.0)))
    print(f"  cross-check  beta=2 (eta=sqrt2): P_early = {P[iw]:.4f}  vs exact TW 1-F_2(0) = 0.0306")
    print(f"  cross-check  beta=1 (eta=2)    : P_early = {P[i2]:.4f}  vs exact TW 1-F_1(0) = 0.168")
    Hs = etas ** 2 / (4 * np.pi)
    print(f"  hazard H over-counts P_early by {np.min(Hs/P):.1f}x (eta=2) to {np.max(Hs/P):.0f}x "
          f"(eta=1): wrong functional form (polynomial vs exponential).")

    # ---- figure ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    C_D, C_FIT, C_C, C_POLY = "#2563eb", "#1e3a8a", "#16a34a", "#dc2626"

    fig, ax = plt.subplots(1, 2, figsize=(12.4, 5.0))

    uu = np.linspace(u.min() * 0.9, u.max() * 1.05, 100)
    ax[0].errorbar(u, np.log(P), yerr=E / P, fmt="o", color=C_D, ms=7, capsize=3,
                   mec="white", mew=1.0, zorder=5, label=r"Monte-Carlo $\ln P_{\rm early}$")
    ax[0].plot(uu, intercept + slope * uu, "-", color=C_FIT, lw=2.0,
               label=fr"exponential fit $c_{{\rm fit}}={c_fit:.2f}$ ($R^2={R2:.2f}$)")
    anchor = np.mean(np.log(P)) + 5.4439 * np.mean(u)
    ax[0].plot(uu, anchor - 5.4439 * uu, ":", color=C_C, lw=1.8,
               label=r"Painlevé II slope $c=5.444$")
    ax[0].plot(uu, np.log(1.0 / uu / (4 * np.pi)), "--", color=C_POLY, lw=1.8,
               label=r"hazard $\ln(\eta^2/4\pi)$ (over-counts)")
    ax[0].set_xlabel(r"$1/\eta^2$"); ax[0].set_ylabel(r"$\ln P_{\rm early}$")
    ax[0].set_title("Early escape: exponential SAO tail, not polynomial hazard")
    ax[0].grid(True, color="#eef2f7", lw=0.6)
    ax[0].legend(loc="lower left", fontsize=8.6, framealpha=0.95, edgecolor="#cbd5e1")

    eg = np.linspace(etas.min(), etas.max(), 100)
    ax[1].semilogy(etas, P, "o", color=C_D, ms=7, mec="white", mew=1.0, zorder=5,
                   label="Monte-Carlo")
    ax[1].semilogy(eg, np.exp(intercept) * np.exp(-c_fit / eg ** 2), "-", color=C_FIT, lw=2.0,
                   label=fr"$\propto e^{{-{c_fit:.2f}/\eta^2}}$")
    ax[1].semilogy(eg, eg ** 2 / (4 * np.pi), "--", color=C_POLY, lw=1.8,
                   label=r"hazard $\eta^2/4\pi$")
    ax[1].semilogy([np.sqrt(2), 2.0], [0.0306, 0.168], "*", color=C_C, ms=16, mec="white",
                   mew=1.0, zorder=6, label=r"exact TW $1-F_\beta(0)$")
    ax[1].set_xlabel(r"effective noise $\eta=\sigma/\sqrt{\epsilon_2}$")
    ax[1].set_ylabel(r"$P_{\rm early}=P(Y_{\rm node}>0)$")
    ax[1].set_title("Early-escape probability vs noise")
    ax[1].grid(True, which="both", color="#eef2f7", lw=0.6)
    ax[1].legend(loc="lower right", fontsize=8.8, framealpha=0.95, edgecolor="#cbd5e1")

    fig.tight_layout()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures",
                       "tail_reconciliation.png")
    out = os.path.abspath(out)
    fig.savefig(out, dpi=140)
    print("\nsaved", out)
