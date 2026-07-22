"""
phase_edge_snic.py — Direction B: the SNIC phase channel = the second (phase) edge.

The folded limit cycle has ONE inner Riccati  dR = (R^2 - Y) dT + eta dB, Cole-Hopf
u'' = (Y - eta xi)u.  The destroying bifurcation selects which edge problem is singular:

  AMPLITUDE edge (fold of cycles):  SWEEP Y = Y0 - T; peel-off = first node during the
        sweep = SAO ground state  =>  Y_node =d TW_beta  (Painleve II, closed form).

  PHASE edge (SNIC):  hold Y at criticality (Y=0); the ROTATION = successive blow-ups of R
        = node SPACINGS of u at fixed energy.  The rotation period is the FIRST-PASSAGE TIME
        of the noisy saddle-node dR = R^2 dT + eta dB; rescaling R=sigma^{2/3}rho, T=sigma^{-2/3}s
        gives the eta-free canonical  drho = rho^2 ds + dW  (so omega, D_phi ~ sigma^{2/3}).

Results pinned here:
  (1) the rate constant J = ∬_{w<u} e^{(w^3-u^3)/3} has a CLOSED FORM.  Doing the inner
      Gaussian integral exactly, J = sqrt(pi) ∫_0^inf v^{-1/2} e^{-v^3/12} dv
        = (sqrt(pi)/3) 12^{1/6} Gamma(1/6) ~ 4.976   [D=1];  (sqrt(2pi)/3) 6^{1/6} Gamma(1/6)
        ~ 6.270  [D=1/2, dB-normalisation].  The naive 2D quadrature has only 1/u^2 tails
      (g(u) ~ 1/u^2), which the closed form sidesteps.
  (2) sigma^{2/3} scaling of omega, D_phi (verified by collapse and exponent fit).
  (3) the first-passage LAW: strongly right-skewed, exponential right tail with rate lambda0;
      lambda0 = E0/2 where E0 is the ground state of the QUARTIC operator -d^2 + rho^4 + 2rho
      (the FPT generator's Schrodinger form).  The anharmonic quartic is non-integrable =>
      NO Tracy-Widom / Painleve closed form for the distribution (only the constant J is closed).
"""
from __future__ import annotations
import os
import math
import numpy as np


# =====================================================================================
# 1. Rate constant J — closed form, plus the naive integrand g(u) ~ 1/u^2 (slow tails)
# =====================================================================================
def J_closed_forms():
    Jp = np.sqrt(np.pi) / 3 * 12 ** (1 / 6) * math.gamma(1 / 6)      # D=1   (paper)
    Jd = np.sqrt(2 * np.pi) / 3 * 6 ** (1 / 6) * math.gamma(1 / 6)   # D=1/2 (dB)
    return Jp, Jd


def g_of_u(u, Nv=5000, vmax=34.0):
    # g(u) = ∫_0^inf exp(-u^2 v + u v^2 - v^3/3) dv ; cubic v-grid resolves the ~1/u^2-wide
    # small-v peak so the 1/u^2 tail is captured at large |u|.
    V = vmax * np.linspace(0.0, 1.0, Nv) ** 3
    UU = u[:, None]; VV = V[None, :]
    integ = np.exp(-UU * UU * VV + UU * VV * VV - VV ** 3 / 3.0)
    return np.trapezoid(integ, V, axis=1)


def compute_J_numeric(U=120.0, n=6000):
    u = np.linspace(-U, U, n)
    g = g_of_u(u)
    c = float(np.median((u[u > 30] ** 2) * g[u > 30]))              # g(u) ~ c/u^2,  c -> 1
    J_num = np.trapezoid(g, u) + 2 * c / U                          # with 1/u^2 tail correction
    return J_num, c, u, g


# =====================================================================================
# 2. Canonical first-passage law  drho = rho^2 ds + dW  (D=1/2, dB-normalisation)
# =====================================================================================
def fpt_canonical(N, noise=1.0, rho_inj=-12.0, rho_esc=12.0, dt=1.5e-3, Tmax=80.0, seed=0):
    rng = np.random.default_rng(seed)
    rho = np.full(N, rho_inj); T = np.full(N, np.nan); alive = np.ones(N, bool)
    nsteps = int(Tmax / dt); sdt = np.sqrt(dt)
    for k in range(nsteps):
        if not alive.any():
            break
        r = rho[alive]
        rho[alive] = r + r * r * dt + noise * sdt * rng.standard_normal(r.size)
        esc = alive & (rho > rho_esc)
        T[np.where(esc)[0]] = k * dt
        alive[esc] = False
    return T


def moments(x):
    m, sd = x.mean(), x.std(); c = x - m
    return m, sd, float((c ** 3).mean() / sd ** 3), float((c ** 4).mean() / sd ** 4 - 3)


# =====================================================================================
# 3. The quartic operator  -chi'' + (rho^4 + 2 rho) chi = 2 lambda chi  (FPT Schrodinger form)
#    Right tail  P(T) ~ e^{-lambda0 T},  lambda0 = E0/2,  E0 = ground state of the quartic.
# =====================================================================================
def quartic_ground(L=6.0, n=1200):
    rho = np.linspace(-L, L, n); h = rho[1] - rho[0]
    main = 2.0 / h ** 2 + rho ** 4 + 2.0 * rho
    off = -1.0 / h ** 2 * np.ones(n - 1)
    E = np.linalg.eigvalsh(np.diag(main) + np.diag(off, 1) + np.diag(off, -1))
    return float(E[0])


if __name__ == "__main__":
    # ---- (1) J: closed form + slow-tail diagnostic ----
    Jp, Jd = J_closed_forms()
    J_num, c_tail, ug, gg = compute_J_numeric()
    print(f"[J] CLOSED FORM  J = (sqrt(pi)/3) 12^(1/6) Gamma(1/6) = {Jp:.4f}   [D=1, paper]")
    print(f"    numeric double-integral check = {J_num:.4f}   (g(u) ~ {c_tail:.2f}/u^2 : slow 1/U tails)")
    print(f"    dB-normalisation period (D=1/2) = (sqrt(2pi)/3) 6^(1/6) Gamma(1/6) = {Jd:.4f}")

    # ---- (2) canonical FPT law (dB-normalisation, D=1/2) ----
    T = fpt_canonical(60000, noise=1.0, dt=1.5e-3, Tmax=45.0, seed=1)
    Te = T[~np.isnan(T)]; cens = np.isnan(T).mean()
    m, sd, sk, ku = moments(Te)
    print(f"\n[FPT] canonical law: mean={m:.3f} (+cutoff {2*c_tail/12:.2f} -> Jd={Jd:.2f})  "
          f"std={sd:.3f}  skew={sk:.3f}  exkurt={ku:.3f}  censored={cens:.4f}")
    ts = np.sort(Te); S = 1.0 - np.arange(1, ts.size + 1) / ts.size
    win = (ts > m + 1.5 * sd) & (S > 2e-3)
    slope, inter = np.polyfit(ts[win], np.log(S[win]), 1); lam0_fit = -slope
    E0 = quartic_ground(); lam0_q = E0 / 2.0
    print(f"[FPT] right-tail rate lambda0(sim) = {lam0_fit:.3f}  vs  quartic E0/2 = {lam0_q:.3f}"
          f"  (E0={E0:.3f})  => anharmonic, no Painleve/TW closed form for the LAW")
    print(f"[TW?] skew {sk:.2f} (TW1/TW2 = 0.29/0.22), exkurt {ku:.2f} (0.17/0.09) => NOT Tracy-Widom")

    # ---- (3) sigma^{2/3} collapse from the physical equation ----
    print("\n[sigma^2/3] physical dR=R^2 dT + sigma dB, rescaled cutoffs +-12 sigma^{2/3}:")
    sig_list = [0.5, 0.25]; resc = {}
    for sg in sig_list:
        cut = 12.0 * sg ** (2.0 / 3.0)
        Tp = fpt_canonical(16000, noise=sg, rho_inj=-cut, rho_esc=cut, dt=2.0e-3,
                           Tmax=42.0 * sg ** (-2.0 / 3.0), seed=5)
        Tp = Tp[~np.isnan(Tp)]; resc[sg] = Tp * sg ** (2.0 / 3.0)
        print(f"   sigma={sg:4.2f}:  <T_phys>={Tp.mean():8.3f}   sigma^2/3<T>={resc[sg].mean():.3f}")
    Tm = np.array([(resc[s] / s ** (2 / 3)).mean() for s in sig_list])
    p_exp = -np.polyfit(np.log(sig_list), np.log(Tm), 1)[0]
    print(f"   fitted exponent  <T_phys> ~ sigma^(-{p_exp:.3f})  (predicted 2/3 = 0.667)")

    # =================================== figure ===================================
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    C_A, C_C, C_K = "#2563eb", "#ea8a0b", "#dc2626"
    fig, ax = plt.subplots(1, 3, figsize=(15.2, 4.6))

    # Panel A: FPT law + exponential tail (= quartic ground state) ; not TW
    ax[0].hist(Te, bins=120, range=(0, 40), density=True, color=C_A, alpha=0.5, label="canonical FPT")
    tt = np.linspace(m + 1.0 * sd, 40, 50)
    ax[0].plot(tt, np.exp(inter) * lam0_fit * np.exp(-lam0_fit * tt), "-", color=C_K, lw=2.0,
               label=fr"tail $\propto e^{{-\lambda_0 T}}$, $\lambda_0={lam0_fit:.2f}=E_0/2$")
    ax[0].axvline(m, color="#475569", lw=1.2, ls="--", label=fr"mean $\approx J$")
    ax[0].set_yscale("log"); ax[0].set_xlabel(r"rotation period $T$ (canonical units)")
    ax[0].set_ylabel("density"); ax[0].set_ylim(1e-4, 1)
    ax[0].set_title(fr"Phase-edge first-passage law (skew {sk:.1f}, not TW)")
    ax[0].legend(fontsize=8.6, framealpha=0.95, edgecolor="#cbd5e1")

    # Panel B: sigma^{2/3} collapse
    cols = {0.5: "#1d4ed8", 0.25: "#16a34a"}
    for sg in sig_list:
        ax[1].hist(resc[sg], bins=80, range=(0, 30), density=True, histtype="step", lw=2.0,
                   color=cols[sg], label=fr"$\sigma={sg}$")
    ax[1].hist(Te, bins=80, range=(0, 30), density=True, histtype="step", lw=1.3, ls="--",
               color="#334155", label="canonical")
    ax[1].set_xlabel(r"$\sigma^{2/3}\, T_{\rm phys}$"); ax[1].set_ylabel("density")
    ax[1].set_title(fr"$\sigma^{{2/3}}$ collapse (fitted exponent {p_exp:.2f})")
    ax[1].legend(fontsize=8.6, framealpha=0.95, edgecolor="#cbd5e1")

    # Panel C: slow 1/u^2 tail of the J integrand vs the closed form
    pos = ug > 2
    ax[2].loglog(ug[pos], gg[pos], color=C_C, lw=2.2, label=r"$g(u)$ (naive J integrand)")
    ax[2].loglog(ug[pos], 1.0 / ug[pos] ** 2, "--", color="#475569", lw=1.6, label=r"$\sim 1/u^2$ tail")
    ax[2].set_xlabel(r"$u$"); ax[2].set_ylabel(r"$g(u)=\int_{-\infty}^u e^{(w^3-u^3)/3}dw$")
    ax[2].set_title(fr"$J=(\sqrt{{\pi}}/3)\,12^{{1/6}}\Gamma(1/6)={Jp:.3f}$  (closed form)")
    ax[2].legend(fontsize=8.8, framealpha=0.95, edgecolor="#cbd5e1")

    fig.tight_layout()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures",
                       "phase_edge_snic.png")
    out = os.path.abspath(out)
    fig.savefig(out, dpi=140)
    print("\nsaved", out)
