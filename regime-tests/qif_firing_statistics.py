"""
qif_firing_statistics.py — Direction #1: the phase edge IS the noisy Type-I (QIF / theta) neuron.

The quadratic-integrate-and-fire neuron  dv = (v^2 + I) dt + sigma dW,  reset v->v_reset when
v -> +inf, is the SNIC normal form (the theta neuron).  It is EXACTLY the off-critical noisy
saddle-node  dR = (R^2 - Y) dT + eta dB  with  R<->v, Y<->-I, eta<->sigma.  So the phase-edge
results are the near-threshold inter-spike-interval (ISI) statistics of a noisy Type-I neuron:

  * the ISI is the first-passage time of the noisy saddle-node  =>  ISI law = the quartic FPT
    law (skewed, NOT exponential or Gaussian);
  * rescaling v=sigma^{2/3} rho, t=sigma^{-2/3} s gives the universal f-I curve
        r(I, sigma) = sigma^{2/3} * Phi(nu),   nu = I / sigma^{4/3},   Phi(nu) = 1 / J(nu),
    J(nu) = canonical mean FPT.  Limits: Phi -> sqrt(nu)/pi for nu>>1 (so r -> sqrt(I)/pi, the
    deterministic Type-I rate) and Phi -> exp(-(8/3)|nu|^{3/2}) for nu<<0 (Arrhenius noise-induced
    firing).  At threshold nu=0: r = sigma^{2/3}/J(0), the universal sigma^{2/3} scaling.
  * the ISI coefficient of variation CV crosses over: ~1 (Poisson-like, sub-threshold Kramers)
    -> ~0.6 (SNIC) -> ->0 (regular firing, supra-threshold).
"""
from __future__ import annotations
import os
import numpy as np
from offcritical_phase_edge import J_of_nu


def Phi(nu):
    nu = np.asarray(nu, dtype=float)
    if nu.ndim == 0:
        return 1.0 / J_of_nu(float(nu))
    return np.array([1.0 / J_of_nu(float(x)) for x in nu.ravel()]).reshape(nu.shape)


def qif_sim(I, sigma, N=1500, dt=2.5e-3, n_isi=16, Tcap=240.0, v_th=14.0, v_reset=-14.0, seed=0):
    rng = np.random.default_rng(seed)
    v = np.full(N, float(v_reset)); last = np.zeros(N)
    nu = I / sigma ** (4.0 / 3.0)
    Tmax = min(n_isi * sigma ** (-2.0 / 3.0) * J_of_nu(nu), Tcap)
    nsteps = int(Tmax / dt); sdt = np.sqrt(dt); isis = []
    for k in range(nsteps):
        t = k * dt
        v += (v * v + I) * dt + sigma * sdt * rng.standard_normal(N)
        fired = v > v_th
        if fired.any():
            idx = np.where(fired)[0]
            isis.extend((t - last[idx]).tolist())
            last[idx] = t; v[fired] = v_reset
    isis = np.array(isis)
    isis = isis[isis > 0]
    r = 1.0 / isis.mean(); cv = isis.std() / isis.mean()
    return nu, r, cv, isis


if __name__ == "__main__":
    nu_grid = np.linspace(-2.0, 6.0, 200)
    Phi_grid = Phi(nu_grid)

    # QIF confirmation points: two noise levels, several currents -> collapse
    sigmas = [0.5, 0.3]; nus = [-0.4, 0.0, 0.7, 1.6]
    pts = []  # (sigma, I, nu, r, cv)
    print(f"{'sigma':>6} {'I':>8} {'nu':>6} {'r(QIF)':>9} {'sigma^2/3*Phi':>13} {'CV':>6}")
    for sg in sigmas:
        for nu in nus:
            I = nu * sg ** (4.0 / 3.0)
            nu_m, r, cv, isis = qif_sim(I, sg, seed=int(100 * sg) + int(10 * nu))
            pts.append((sg, I, nu_m, r, cv))
            print(f"{sg:6.2f} {I:8.4f} {nu_m:6.2f} {r:9.4f} {sg**(2/3)*Phi(nu):13.4f} {cv:6.2f}")
    pts = np.array(pts)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    fig, ax = plt.subplots(1, 3, figsize=(15.4, 4.7))
    cm = {0.5: "#dc2626", 0.3: "#2563eb"}

    # Panel A: noise-smoothed Type-I f-I curve
    Ig = np.linspace(-0.6, 1.2, 300)
    for sg in sigmas:
        ax[0].plot(Ig, sg ** (2/3) * Phi(Ig / sg ** (4/3)), color=cm[sg], lw=2.2, label=fr"$\sigma={sg}$")
    Ip = np.linspace(1e-4, 1.2, 100)
    ax[0].plot(Ip, np.sqrt(Ip) / np.pi, "--", color="#0f172a", lw=1.6, label=r"deterministic $\sqrt{I}/\pi$")
    for sg, I, nu, r, cv in pts:
        ax[0].plot(I, r, "o", color=cm[sg], ms=7, mec="white", mew=1.0, zorder=5)
    ax[0].axvline(0, color="#94a3b8", lw=1.0, ls=":")
    ax[0].set_xlabel(r"input current $I$"); ax[0].set_ylabel(r"firing rate $r$")
    ax[0].set_title("Noise-smoothed Type-I f–I curve")
    ax[0].legend(fontsize=9, framealpha=0.95, edgecolor="#cbd5e1")

    # Panel B: universal scaling collapse  r*sigma^{-2/3}  vs  nu = I/sigma^{4/3}
    ax[1].plot(nu_grid, Phi_grid, color="#0f172a", lw=2.4, label=r"$\Phi(\nu)=1/J(\nu)$ (universal)")
    npos = nu_grid[nu_grid > 0.3]; ax[1].plot(npos, np.sqrt(npos) / np.pi, "--", color="#16a34a",
                                              lw=1.6, label=r"$\sqrt{\nu}/\pi$ (deterministic)")
    nneg = nu_grid[nu_grid < -0.5]
    ax[1].plot(nneg, np.sqrt(np.abs(nneg)) / np.pi * np.exp(-(8 / 3) * np.abs(nneg) ** 1.5),
               ":", color="#dc2626", lw=1.8,
               label=r"Arrhenius $\frac{\sqrt{|\nu|}}{\pi}e^{-\frac{8}{3}|\nu|^{3/2}}$")
    for sg in sigmas:
        sel = pts[pts[:, 0] == sg]
        ax[1].plot(sel[:, 2], sel[:, 3] * sg ** (-2/3), "o", color=cm[sg], ms=8, mec="white",
                   mew=1.0, zorder=5, label=fr"QIF $\sigma={sg}$")
    ax[1].set_xlabel(r"$\nu = I/\sigma^{4/3}$"); ax[1].set_ylabel(r"$r\,\sigma^{-2/3}=\Phi(\nu)$")
    ax[1].set_title(r"Universal collapse (SNIC $\sigma^{2/3}$ scaling)")
    ax[1].set_ylim(0, 1.0); ax[1].legend(fontsize=8.6, framealpha=0.95, edgecolor="#cbd5e1")

    # Panel C: ISI regularity (CV) crossover
    for sg in sigmas:
        sel = pts[pts[:, 0] == sg]
        ax[2].plot(sel[:, 2], sel[:, 4], "o-", color=cm[sg], ms=7, mec="white", mew=1.0,
                   label=fr"QIF $\sigma={sg}$")
    ax[2].axhline(1.0, color="#94a3b8", lw=1.0, ls=":")
    ax[2].text(-1.7, 0.92, "Poisson-like\n(sub-threshold,\nKramers)", fontsize=8, color="#dc2626")
    ax[2].text(0.05, 0.5, "SNIC", fontsize=8.5, color="#7c3aed")
    ax[2].text(1.0, 0.18, "regular\n(supra-threshold)", fontsize=8, color="#16a34a")
    ax[2].set_xlabel(r"$\nu = I/\sigma^{4/3}$"); ax[2].set_ylabel("ISI coefficient of variation")
    ax[2].set_title("Spike-train regularity crossover")
    ax[2].set_ylim(0, 1.15); ax[2].legend(fontsize=9, framealpha=0.95, edgecolor="#cbd5e1")

    fig.tight_layout()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures",
                       "qif_firing_statistics.png")
    out = os.path.abspath(out)
    fig.savefig(out, dpi=140); print("saved", out)
