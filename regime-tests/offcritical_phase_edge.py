"""
offcritical_phase_edge.py — Direction B, off-critical: the SNIC phase edge as the centre of a
one-parameter crossover.

Off criticality the inner Riccati is  dR = (mu + R^2) dT + sigma dB  (Y = -mu).  Rescaling
R = sigma^{2/3} rho,  T = sigma^{-2/3} s  gives the eta-free ONE-PARAMETER family
        drho = (nu + rho^2) ds + dW ,        nu = mu / sigma^{4/3}   (rescaled detuning).
The first-passage (rotation) generator's Schrodinger form is the TILTED QUARTIC family
        - chi'' + ( rho^4 + 2 nu rho^2 + 2 rho + nu^2 ) chi = 2 lambda chi .

Three regimes in the single parameter nu:
   nu < 0  : DOUBLE WELL (fixed points +-sqrt|nu|) -> Kramers activated escape;
             mean period ~ exp( (8/3) |nu|^{3/2} )  (the Arrhenius |.|^{3/2} barrier);
   nu = 0  : SNIC critical edge (pure quartic; mean = closed-form J; sigma^{2/3} scaling);
   nu > 0  : SINGLE WELL, no fixed point -> deterministic rotation, period -> pi/sqrt(nu),
             regular phase diffusion D_phi ~ sigma^2 (the fold-of-cycles / finite-period limit).
"""
from __future__ import annotations
import os
import math
import numpy as np


def J_of_nu(nu, L=12.0, n=24000):
    """Mean rotation period = MFPT(-inf -> +inf) for drho=(nu+rho^2)ds+dW (D=1/2), by quadrature."""
    rho = np.linspace(-L, L, n); drho = rho[1] - rho[0]
    Phi = 2.0 * (nu * rho + rho ** 3 / 3.0)
    logcum = np.logaddexp.accumulate(Phi + np.log(drho))
    integ = np.exp(logcum - Phi)
    return 2.0 * np.trapezoid(integ, rho) + 2.0 / L


def quartic_E0(nu, L=6.0, n=600):
    """Ground state of -d^2 + (rho^4 + 2 nu rho^2 + 2 rho + nu^2) = the FPT 2*lambda0."""
    rho = np.linspace(-L, L, n); h = rho[1] - rho[0]
    main = 2.0 / h ** 2 + rho ** 4 + 2.0 * nu * rho ** 2 + 2.0 * rho + nu ** 2
    off = -1.0 / h ** 2 * np.ones(n - 1)
    return float(np.linalg.eigvalsh(np.diag(main) + np.diag(off, 1) + np.diag(off, -1))[0])


def fpt_sim(nu, N=14000, dt=2.5e-3, seed=0):
    rng = np.random.default_rng(seed)
    well = np.sqrt(abs(nu)) if nu < 0 else 0.0
    rho_inj = -(well + 6.0); rho_esc = well + 8.0
    rho = np.full(N, rho_inj); T = np.full(N, np.nan); alive = np.ones(N, bool)
    Tmax = min(6.0 * J_of_nu(nu), 330.0); nsteps = int(Tmax / dt); sdt = np.sqrt(dt)
    for k in range(nsteps):
        if not alive.any():
            break
        r = rho[alive]
        rho[alive] = r + (nu + r * r) * dt + sdt * rng.standard_normal(r.size)
        esc = alive & (rho > rho_esc)
        T[np.where(esc)[0]] = k * dt; alive[esc] = False
    Te = T[~np.isnan(T)]
    m, sd = Te.mean(), Te.std()
    return Te, m, sd, float(((Te - m) ** 3).mean() / sd ** 3)


if __name__ == "__main__":
    J_dB = np.sqrt(2 * np.pi) / 3 * 6 ** (1 / 6) * math.gamma(1 / 6)      # closed form J(0), D=1/2
    nus = np.linspace(-3.2, 6.0, 80)
    Jq = np.array([J_of_nu(nu) for nu in nus])
    print(f"  J(nu=0) quadrature = {J_of_nu(0.0):.4f}   vs closed form = {J_dB:.4f}")
    print(f"  {'nu':>6} {'J(nu)':>12} {'lambda0=E0/2':>13}")
    for nu in [-2.5, -1.0, 0.0, 1.0, 3.0]:
        print(f"  {nu:6.2f} {J_of_nu(nu):12.4f} {quartic_E0(nu)/2:13.4f}")

    # sims: confirm J(nu) and get the shape crossover (rescaled FPT distributions)
    sim_nu = [-1.0, 0.0, 2.0, 4.0]; sims = {}
    print("  -- sim check (mean period) & shape --")
    for nu in sim_nu:
        Te, m, sd, sk = fpt_sim(nu, seed=3); sims[nu] = (Te, m, sd, sk)
        print(f"     nu={nu:+.1f}:  sim mean={m:7.3f} (quad {J_of_nu(nu):7.3f})  CV={sd/m:.2f}  skew={sk:.2f}")

    # asymptotes
    nn = nus[nus < -0.3]; kram = np.exp((8.0 / 3.0) * np.abs(nn) ** 1.5)
    kram *= J_of_nu(-2.0) / np.exp((8.0 / 3.0) * 2.0 ** 1.5)
    npd = nus[nus > 0.4]; det = np.pi / np.sqrt(npd)

    # ============================ figure ============================
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    C_J, C_K, C_D, C_0 = "#2563eb", "#dc2626", "#16a34a", "#7c3aed"
    fig, ax = plt.subplots(1, 3, figsize=(15.2, 4.7))

    # Panel A: J(nu) crossover
    ax[0].semilogy(nus, Jq, color=C_J, lw=2.4, label=r"$J(\nu)$ (MFPT quadrature)")
    ax[0].semilogy(nn, kram, "--", color=C_K, lw=1.8, label=r"Kramers $\propto e^{(8/3)|\nu|^{3/2}}$")
    ax[0].semilogy(npd, det, "--", color=C_D, lw=1.8, label=r"deterministic $\pi/\sqrt{\nu}$")
    ax[0].semilogy(sim_nu, [sims[nu][1] for nu in sim_nu], "o", color="#0f172a", ms=7,
                   mec="white", mew=1.0, zorder=5, label="simulation")
    ax[0].plot(0, J_of_nu(0.0), "*", color=C_0, ms=16, mec="white", mew=1.0, zorder=6)
    ax[0].text(0.15, J_of_nu(0) * 1.7, "SNIC edge\n$\\nu{=}0$", color=C_0, fontsize=8.5)
    ax[0].text(-3.05, 80, "Kramers\n(Arrhenius)", color=C_K, fontsize=8.5)
    ax[0].text(3.0, 1.7, "deterministic\nrotation", color=C_D, fontsize=8.5)
    ax[0].set_xlabel(r"rescaled detuning  $\nu=\mu/\sigma^{4/3}$")
    ax[0].set_ylabel(r"mean rotation period  $J(\nu)$")
    ax[0].set_title("Off-critical crossover (one parameter)")
    ax[0].set_ylim(0.6, 400); ax[0].legend(fontsize=8.2, framealpha=0.95, edgecolor="#cbd5e1")

    # Panel B: the quartic potential family
    rr = np.linspace(-2.7, 2.7, 400)
    for nu, col, lab in [(-2.5, C_K, r"$\nu=-2.5$ (double well)"),
                         (0.0, C_0, r"$\nu=0$ (critical)"),
                         (2.5, C_D, r"$\nu=2.5$ (single well)")]:
        ax[1].plot(rr, rr ** 4 + 2 * nu * rr ** 2 + 2 * rr + nu ** 2, color=col, lw=2.2, label=lab)
    ax[1].set_xlabel(r"$\rho$")
    ax[1].set_ylabel(r"$V(\rho)=\rho^4+2\nu\rho^2+2\rho+\nu^2$")
    ax[1].set_title("FPT-operator potential family")
    ax[1].set_ylim(-2, 26); ax[1].legend(fontsize=8.4, framealpha=0.95, edgecolor="#cbd5e1")

    # Panel C: the FPT shape crossover (rescaled to unit mean)
    cmap = {-1.0: C_K, 0.0: C_0, 2.0: "#0891b2", 4.0: C_D}
    for nu in sim_nu:
        Te, m, sd, sk = sims[nu]
        ax[2].hist(Te / m, bins=70, range=(0, 4), density=True, histtype="step", lw=2.0,
                   color=cmap[nu], label=fr"$\nu={nu:+.0f}$ (CV {sd/m:.2f}, skew {sk:.1f})")
    ax[2].set_xlabel(r"$T/\langle T\rangle$"); ax[2].set_ylabel("density")
    ax[2].set_title("FPT shape: broad/exponential ($\\nu{<}0$) $\\to$ peaked ($\\nu{>}0$)")
    ax[2].legend(fontsize=8.2, framealpha=0.95, edgecolor="#cbd5e1")

    fig.tight_layout()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures",
                       "offcritical_phase_edge.png")
    out = os.path.abspath(out)
    fig.savefig(out, dpi=140)
    print("\nsaved", out)
