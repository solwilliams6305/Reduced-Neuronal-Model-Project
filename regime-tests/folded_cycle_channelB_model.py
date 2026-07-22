"""
folded_cycle_channelB_model.py — Part 2 item #2: Channel B on the real model
----------------------------------------------------------------------------
Measure the phase-diffusion coefficient D_phi of the Bautin fold-of-cycles as the
SNLC is approached, with degenerate noise on the fast Cartesian variable u.

    du = [ u (beta + R^2 - R^4) - Omega v ] dt + sigma dW,   R^2=u^2+v^2
    dv = [ v (beta + R^2 - R^4) + Omega u ] dt
    (fixed beta; stable cycle rho_+(beta); SNLC at beta=-1/4, rho=1/sqrt2)

D_phi = Var(phi_accumulated)/T across an ensemble (deterministic rotation cancels).

THE QUESTION: does D_phi BLOW UP as beta -> beta_SN (Channel B enhancement), or stay
finite?  Leading-order prediction for a SMOOTH cycle (u-noise projected on phase):
    D_phi^lead = sigma^2 <sin^2>/rho_+^2 = sigma^2 / (2 rho_+^2)   -> FINITE at SNLC.
Any enhancement must come from the TRANSVERSE (2nd-order) correction: as the radial
Floquet rate |lambda| = |F'(rho_+)| -> 0 at the SNLC, <r^2> ~ sigma^2/|lambda| blows
up and leaks into phase, giving a sigma^4 / |lambda| term.  So expect: D_phi/sigma^2
~ const (finite) at small sigma, RISING near the fold / at larger sigma.
"""
from __future__ import annotations
import os
import numpy as np

BETA_SN = -0.25
RHO_SN = 1.0 / np.sqrt(2.0)


def rho_plus(beta):
    return np.sqrt((1.0 + np.sqrt(1.0 + 4.0 * beta)) / 2.0)


def lam_perp(beta):
    """radial Floquet/relaxation rate |F'(rho_+)|, F=rho(beta+rho^2-rho^4)."""
    rp = rho_plus(beta)
    return abs(beta + 3 * rp**2 - 5 * rp**4)


def D_phi(beta, sigma, N=150, Omega=2 * np.pi, T=30.0, dt=2e-3, seed=0):
    rng = np.random.default_rng(seed)
    rp = rho_plus(beta)
    u = np.full(N, rp); v = np.zeros(N)
    phi_prev = np.zeros(N)              # atan2(0, rp)=0
    phi_acc = np.zeros(N)
    sdt = np.sqrt(dt)
    n = int(T / dt)
    for _ in range(n):
        R2 = u * u + v * v
        g = beta + R2 - R2 * R2
        u = u + (u * g - Omega * v) * dt + sigma * sdt * rng.standard_normal(N)
        v = v + (v * g + Omega * u) * dt
        phi = np.arctan2(v, u)
        dphi = phi - phi_prev
        dphi = (dphi + np.pi) % (2 * np.pi) - np.pi    # unwrap to (-pi,pi]
        phi_acc += dphi
        phi_prev = phi
    return float(np.var(phi_acc) / T)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    figdir = os.path.join(os.path.dirname(here), "figures")
    os.makedirs(figdir, exist_ok=True)

    deltas = [0.30, 0.15, 0.07, 0.03]      # beta = beta_SN + delta (approach the SNLC)
    sigmas = np.array([0.02, 0.05, 0.10, 0.18, 0.30])
    N = 150

    print("\n=== Part 2 #2: Channel B (phase diffusion) on the Bautin fold-of-cycles ===")
    print(f"  SNLC beta=-0.25, rho_SN={RHO_SN:.3f}.  D_phi/sigma^2 leading = 1/(2 rho_+^2).\n")
    print(f"  {'delta':>6s} {'rho_+':>6s} {'|lambda|':>9s} {'lead 1/2rho^2':>13s}"
          f"  D_phi/sigma^2 vs sigma ->")
    results = {}
    for d in deltas:
        beta = BETA_SN + d
        rp = rho_plus(beta); lam = lam_perp(beta)
        lead = 1.0 / (2 * rp**2)
        row = []
        for k, s in enumerate(sigmas):
            Dp = D_phi(beta, s, N=N, seed=int(1000 * d) + k)
            row.append(Dp / s**2)
        results[d] = np.array(row)
        print(f"  {d:6.2f} {rp:6.3f} {lam:9.4f} {lead:13.3f}   "
              + "  ".join(f"{x:6.2f}" for x in row))

    print("\n  Reading:")
    print("   - small-sigma column ~ the finite leading 1/(2 rho_+^2): Channel B does")
    print("     NOT blow up at leading order on this SMOOTH fold-of-cycles.")
    print("   - rise at larger sigma / smaller delta = the transverse (2nd-order)")
    print("     enhancement as |lambda|->0; this is the only Channel-B growth here.")
    # quantify: ratio of D_phi/sigma^2 at largest sigma vs smallest, per delta
    print("\n  enhancement factor (D/sig^2 at sigma_max)/(at sigma_min) per delta:")
    for d in deltas:
        r = results[d]
        print(f"     delta={d:.2f}: {r[-1]/r[0]:.2f}   (|lambda|={lam_perp(BETA_SN+d):.3f})")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    for d in deltas:
        rp = rho_plus(BETA_SN + d)
        ax[0].plot(sigmas, results[d], "o-", ms=4,
                   label=fr"$\delta$={d}, $|\lambda|$={lam_perp(BETA_SN+d):.2f}")
        ax[0].axhline(1 / (2 * rp**2), color="gray", lw=0.5, ls=":")
    ax[0].set_xlabel(r"$\sigma$"); ax[0].set_ylabel(r"$D_\phi/\sigma^2$")
    ax[0].set_title(r"Channel B on Bautin: leading $D_\phi/\sigma^2$ finite;"
                    "\n" r"transverse rise near SNLC (small $\delta$, large $\sigma$)")
    ax[0].legend(fontsize=8, frameon=False); ax[0].grid(alpha=0.3)

    # D_phi/sigma^2 (small sigma) vs |lambda| : finite, ~1/(2rho^2), not diverging
    lams = np.array([lam_perp(BETA_SN + d) for d in deltas])
    small = np.array([results[d][0] for d in deltas])
    large = np.array([results[d][-1] for d in deltas])
    ax[1].plot(lams, small, "o-", color="C0", label=r"small $\sigma$ (leading)")
    ax[1].plot(lams, large, "s-", color="C3", label=r"large $\sigma$ (+transverse)")
    ax[1].set_xlabel(r"radial Floquet rate $|\lambda_\perp|$ (→0 at SNLC)")
    ax[1].set_ylabel(r"$D_\phi/\sigma^2$")
    ax[1].set_title(r"Leading $D_\phi$ stays finite as $|\lambda_\perp|\to0$;"
                    "\n" r"only the large-$\sigma$ (transverse) part grows")
    ax[1].legend(fontsize=9, frameon=False); ax[1].grid(alpha=0.3)
    ax[1].invert_xaxis()
    fig.tight_layout()
    out = os.path.join(figdir, "folded_cycle_channelB_model.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print(f"\n  figure -> figures/folded_cycle_channelB_model.png\n")


if __name__ == "__main__":
    main()
