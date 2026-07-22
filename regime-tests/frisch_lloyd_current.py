"""
frisch_lloyd_current.py
=======================
Direction A, first move (from FOLDED_CYCLE_NOISYAIRY_CONJECTURES.md, Conjecture 2 / §7):
the Frisch-Lloyd current of the inner Riccati, and where it departs from quasi-static Kramers.

The inner (K2) Riccati is   dR = (R^2 - Y) dT + eta dB,   Y = Y0 - T.
Frozen at a level Y, the escape rate (= rate of finite-time blow-ups R -> +inf, the
rotation number / Frisch-Lloyd current of u'' = (Y - eta xi) u) is the instantaneous
hazard.  The quasi-static (WKB/Arrhenius) approximation is

        lambda_qs(Y) = (sqrt(Y)/pi) * exp( -8 Y^{3/2} / (3 eta^2) ),     Y > 0,

whose sweep-integral reproduces the closed-form integrated hazard H = eta^2 / (4 pi).
It must break in the noisy-Airy window  Y <~ eta^{4/3}  (inner scales R_* ~ eta^{2/3},
Y_* ~ eta^{4/3}); that breakdown is the sub-exponential prefactor we are after.

KEY SIMPLIFICATION (exact rescaling).  Put  R = eta^{2/3} rho,  Y = eta^{4/3} y,
T = eta^{-2/3} tau.  The SDE becomes eta-INDEPENDENT:

        d rho = (rho^2 - y) d tau + dW .

So the physical current factorises as   lambda(Y, eta) = eta^{2/3} * J(y),  y = Y / eta^{4/3},
with ONE universal current J(y).  We compute J(y) exactly (quadrature) and by Monte-Carlo,
and compare to the rescaled Kramers form  J_qs(y) = (sqrt(y)/pi) e^{-8 y^{3/2}/3}.

Two readings of the result:
  * for y >> 1 :  J(y) -> J_qs(y)              (Kramers is correct deep in the barrier regime);
  * for y <~ 1 :  J(y) departs upward, and J(0) > 0 although J_qs(0) = 0
                  -- escape persists at and past the fold; this is the inner correction.
"""
from __future__ import annotations
import os
import numpy as np

# ----------------------------------------------------------------------------------
# 1. Exact universal current J(y) by the Frisch-Lloyd / mean-first-passage quadrature.
#    1/J(y) = T(y) = 2 * \int dz e^{-Phi(z)} \int_{-inf}^{z} e^{Phi(w)} dw,
#    Phi(rho) = 2 ( rho^3/3 - y rho )   (the eta=1 inner potential exponent).
#    Computed overflow-safe in log space via np.logaddexp.accumulate.
# ----------------------------------------------------------------------------------
def J_quadrature(y, L=30.0, n=60000):
    rho = np.linspace(-L, L, n)
    drho = rho[1] - rho[0]
    Phi = 2.0 * (rho ** 3 / 3.0 - y * rho)
    log_inner = np.logaddexp.accumulate(Phi + np.log(drho))     # log \int_{-L}^{rho} e^{Phi} drho'
    integrand = np.exp(log_inner - Phi)                          # e^{-Phi} \int_{-inf}^{rho} e^{Phi}
    T = 2.0 * np.trapezoid(integrand, rho)
    T += 2.0 / L            # analytic 1/(2 rho^2) tails on both sides (integrand factor 2)
    return 1.0 / T


def J_kramers(y):
    """Quasi-static Arrhenius current (rescaled, eta=1). Zero for y<=0 (no barrier)."""
    y = np.asarray(y, float)
    out = np.where(y > 0, np.sqrt(np.clip(y, 0, None)) / np.pi
                   * np.exp(-8.0 * np.clip(y, 0, None) ** 1.5 / 3.0), 0.0)
    return out


# ----------------------------------------------------------------------------------
# 2. Monte-Carlo current: simulate d rho = (rho^2 - y) d tau + dW with reinjection.
#    Rate = (# blow-ups past rho_esc) / (N * tau_measured).
# ----------------------------------------------------------------------------------
def J_montecarlo(y, N=1500, tau=250.0, dt=5e-3, rho_esc=8.0, rho_inj=-6.0,
                 burn_frac=0.15, seed=0):
    rng = np.random.default_rng(seed)
    rho0 = -np.sqrt(y) if y > 0 else -2.0
    rho = np.full(N, rho0)
    nsteps = int(tau / dt)
    burn = int(burn_frac * nsteps)
    sdt = np.sqrt(dt)
    escapes = 0
    for k in range(nsteps):
        rho += (rho * rho - y) * dt + sdt * rng.standard_normal(N)
        esc = rho > rho_esc
        if k >= burn:
            escapes += int(esc.sum())
        if esc.any():
            rho[esc] = rho_inj
    tau_meas = (nsteps - burn) * dt
    rate = escapes / (N * tau_meas)
    err = np.sqrt(max(escapes, 1)) / (N * tau_meas)          # Poisson 1-sigma
    return rate, err


if __name__ == "__main__":
    # ---- universal current on a fine grid (quadrature) ----
    yq = np.linspace(-1.5, 3.5, 200)
    Jq = np.array([J_quadrature(y) for y in yq])
    Jk = J_kramers(yq)

    # ---- Monte-Carlo validation at representative levels ----
    y_mc = np.array([-0.5, 0.0, 0.5, 1.0, 1.5, 2.0])
    Jmc = np.empty_like(y_mc); Jmc_err = np.empty_like(y_mc)
    print(f"{'y':>6} {'J_quad':>12} {'J_MC':>12} {'+/-':>11} {'J_Kramers':>12} {'J/Kram':>9}")
    for i, y in enumerate(y_mc):
        Jmc[i], Jmc_err[i] = J_montecarlo(y, seed=100 + i)
        jk = float(J_kramers(np.array([y]))[0])
        jq = J_quadrature(y)
        ratio = jq / jk if jk > 0 else np.inf
        print(f"{y:6.2f} {jq:12.5e} {Jmc[i]:12.5e} {Jmc_err[i]:11.3e} {jk:12.5e} {ratio:9.3f}")

    # ---- headline diagnostics ----
    J0 = J_quadrature(0.0)
    # large-y check: quadrature vs Kramers
    yhi = 3.0
    print(f"\n  J(0)        = {J0:.5e}   (quasi-static Kramers predicts exactly 0 at the fold)")
    print(f"  J({yhi})/Kramers = {J_quadrature(yhi)/float(J_kramers(np.array([yhi]))[0]):.4f}"
          f"   (-> 1 confirms Kramers deep in the barrier regime)")
    # multiplicative correction kappa = J / Kramers
    ypos = yq[yq > 0.2]
    kappa = np.array([J_quadrature(y) for y in ypos]) / J_kramers(ypos)
    kmin = kappa.min(); ykmin = ypos[kappa.argmin()]
    print(f"  Kramers OVERestimates most at y ~ {ykmin:.2f}:  J/Kramers = {kmin:.3f}"
          f"  ({(1 - kmin) * 100:.0f}% high)")
    print(f"  J/Kramers -> {kappa[-1]:.3f} at y={ypos[-1]:.2f} (Kramers recovered);"
          f"  diverges as y->0+ since Kramers->0 while J(0)={J0:.3f}>0.")
    print(f"  => correction is an O({(1 - kmin) * 100:.0f}%) prefactor for y=O(1),"
          f" and QUALITATIVE (finite rate) at/below the fold y<=0.")

    # ----------------------------------------------------------------------------------
    # 3. Figure.
    # ----------------------------------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    C_Q, C_MC, C_K, C_ACC = "#2563eb", "#dc2626", "#64748b", "#16a34a"

    fig, ax = plt.subplots(1, 2, figsize=(12.2, 5.0))

    # Panel A: the universal current
    ax[0].semilogy(yq, Jq, color=C_Q, lw=2.4, label=r"Frisch–Lloyd quadrature $\mathcal{J}(y)$")
    ax[0].semilogy(yq[yq > 0], Jk[yq > 0], "--", color=C_K, lw=1.8,
                   label=r"quasi-static Kramers $\frac{\sqrt{y}}{\pi}e^{-8y^{3/2}/3}$")
    ax[0].errorbar(y_mc, Jmc, yerr=Jmc_err, fmt="o", color=C_MC, ms=6, capsize=3,
                   mec="white", mew=1.0, zorder=5, label="Monte-Carlo (inner SDE)")
    ax[0].axvline(0.0, color="#94a3b8", lw=1.0, ls=":")
    ax[0].axvspan(-1.5, 0.0, color=C_ACC, alpha=0.07)
    ax[0].text(-0.72, Jq[0]*0.6, "early escape\n(past the fold;\nKramers $=0$ here)",
               color="#15803d", fontsize=8.5, ha="center", va="top")
    ax[0].set_xlabel(r"rescaled level  $y = Y/\eta^{4/3}$")
    ax[0].set_ylabel(r"escape current  $\mathcal{J}(y)$")
    ax[0].set_title("Universal inner escape current")
    ax[0].set_xlim(-1.5, 3.5); ax[0].set_ylim(1e-8, 5)
    ax[0].grid(True, which="both", color="#eef2f7", lw=0.6)
    ax[0].legend(loc="upper right", fontsize=9, framealpha=0.95, edgecolor="#cbd5e1")

    # Panel B: multiplicative correction kappa(y) = J / Kramers
    yb = yq[(yq > 0.1) & (yq < 3.2)]
    kappab = np.array([J_quadrature(y) for y in yb]) / J_kramers(yb)
    kmin = kappab.min(); ykmin = yb[kappab.argmin()]
    ax[1].plot(yb, kappab, color=C_Q, lw=2.4)
    ax[1].axhline(1.0, color=C_K, lw=1.6, ls="--", label="Kramers exact (ratio $=1$)")
    ax[1].plot([ykmin], [kmin], "o", color=C_Q, ms=6, mec="white", mew=1.0, zorder=5)
    ax[1].annotate(f"Kramers overestimates\nby ~{(1 - kmin) * 100:.0f}% near $y\\approx{ykmin:.1f}$",
                   xy=(ykmin, kmin), xytext=(1.9, 0.83), color="#1e3a8a", fontsize=9, ha="center",
                   arrowprops=dict(arrowstyle="-|>", color="#1e3a8a", lw=1.2))
    ax[1].annotate("$\\to\\infty$ as $y\\to0^+$\n(Kramers vanishes,\n$\\mathcal{J}$ stays finite)",
                   xy=(yb[0], kappab[0]), xytext=(0.85, 1.5), color="#15803d", fontsize=8.5, ha="center",
                   arrowprops=dict(arrowstyle="-|>", color="#15803d", lw=1.1))
    ax[1].text(2.7, 1.02, "Kramers\nrecovered", color="#475569", fontsize=8.5, ha="center", va="bottom")
    ax[1].set_xlabel(r"rescaled level  $y$")
    ax[1].set_ylabel(r"$\mathcal{J}(y)\,/\,\left[\frac{\sqrt{y}}{\pi}e^{-8y^{3/2}/3}\right]$")
    ax[1].set_title("Sub-exponential correction to the quasi-static rate")
    ax[1].set_xlim(0.1, 3.2); ax[1].set_ylim(0.7, 1.75)
    ax[1].grid(True, which="both", color="#eef2f7", lw=0.6)
    ax[1].legend(loc="upper right", fontsize=9, framealpha=0.95, edgecolor="#cbd5e1")

    fig.tight_layout()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures",
                       "frisch_lloyd_current.png")
    out = os.path.abspath(out)
    fig.savefig(out, dpi=140)
    np.savez(os.path.join(os.path.dirname(os.path.abspath(__file__)), "frisch_lloyd_current.npz"),
             yq=yq, Jq=Jq, Jk=Jk, y_mc=y_mc, Jmc=Jmc, Jmc_err=Jmc_err, J0=J0)
    print("\nsaved", out)
