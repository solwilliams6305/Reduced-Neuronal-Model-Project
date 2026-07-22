"""
folded_cycle_rate_constant.py — pinning the early-escape rate constant c
------------------------------------------------------------------------
P_early(eta) ~ exp(-c/eta^2).  c is the eta->0 large-deviation rate.  By the
optimal-fluctuation (Lifshitz/instanton) principle for the stochastic-Airy ground
state lower tail:
    P(Lambda0(beta)<0) ~ exp(-beta A*),  A* = (1/8) min_psi K[psi]^2/int psi^4 ,
    K[psi]=int_0^inf (psi'^2 + x psi^2) dx ,  psi(0)=0 ,   c = 4 A* (beta=4/eta^2).
After optimising the x-scale, the optimal SHAPE solves the nonlinear-Airy ground state
    -g'' + z g = g^3 ,  g(0)=0 , decaying,
and  c = 2 (int g'^2)(int z g^2)/(int g^4).  Pohozaev identities (multiply the EL by g
and by z g') give  P=R, S=2P  with P=int g'^2, R=int z g^2, S=int g^4, hence  c = P.

RESULT: c = 5.4439  (NOT 2*pi=6.2832).  The single-fit 6.22 over eta in [1,2] is
pre-asymptotic; the local slope of ln P_early vs 1/eta^2 descends toward 5.444 as eta->0.
"""
from __future__ import annotations
import os
import numpy as np


def instanton(dz=5e-4, zmax=12.0):
    """Solve -g''+z g=g^3 (g(0)=0) by shooting; return z,g,g' and (P,R,S,c)."""
    def shoot(s):
        n = int(zmax / dz); g = 0.0; gp = s; z = 0.0
        Z = np.empty(n); G = np.empty(n); GP = np.empty(n)
        for i in range(n):
            Z[i] = z; G[i] = g; GP[i] = gp
            k1g, k1p = gp, z * g - g**3
            k2g, k2p = gp + dz / 2 * k1p, (z + dz / 2) * (g + dz / 2 * k1g) - (g + dz / 2 * k1g)**3
            k3g, k3p = gp + dz / 2 * k2p, (z + dz / 2) * (g + dz / 2 * k2g) - (g + dz / 2 * k2g)**3
            k4g, k4p = gp + dz * k3p, (z + dz) * (g + dz * k3g) - (g + dz * k3g)**3
            g += dz / 6 * (k1g + 2 * k2g + 2 * k3g + k4g); gp += dz / 6 * (k1p + 2 * k2p + 2 * k3p + k4p); z += dz
            if g < 0: return 'node', Z[:i + 1], G[:i + 1], GP[:i + 1]
            if abs(g) > 20: return 'blow', Z[:i + 1], G[:i + 1], GP[:i + 1]
        return 'decay', Z, G, GP
    lo, hi = 1.5, 3.0
    for _ in range(58):
        s = 0.5 * (lo + hi); tag = shoot(s)[0]
        if tag == 'node': hi = s
        else: lo = s
    s = 0.5 * (lo + hi); tag, Z, G, GP = shoot(s)
    imax = int(np.argmax(G)); cut = len(G)
    for i in range(imax + 2, len(G) - 1):
        if GP[i] > 0: cut = i; break
    Z, G, GP = Z[:cut], G[:cut], GP[:cut]
    P = np.trapezoid(GP**2, Z); R = np.trapezoid(Z * G**2, Z); S = np.trapezoid(G**4, Z)
    return Z, G, GP, P, R, S, 2 * P * R / S


# measured P_early(eta) (folded_cycle_prefactor.py + c_trend, larger-N at small eta)
ETA = np.array([0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
PE  = np.array([5e-5, 3.2e-4, 1.43e-3, 9.4e-3, 2.89e-2, 6.25e-2, 1.093e-1, 1.654e-1])


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    figdir = os.path.join(os.path.dirname(here), "figures")
    os.makedirs(figdir, exist_ok=True)

    Z, G, GP, P, R, S, c = instanton()
    print("\n=== Pinning c via the nonlinear-Airy instanton ===")
    print(f"  P=int g'^2 = {P:.5f}   R=int z g^2 = {R:.5f}   (Pohozaev P=R: {abs(P-R):.1e})")
    print(f"  S=int g^4  = {S:.5f}   (Pohozaev S=2P: {abs(S-2*P):.1e})")
    print(f"  c = 2PR/S = P = {c:.5f}    [2*pi = {2*np.pi:.4f}]\n")

    inv = 1 / ETA**2; lnP = np.log(PE)
    mid = 0.5 * (ETA[:-1] + ETA[1:])
    slopes = -(lnP[:-1] - lnP[1:]) / (inv[:-1] - inv[1:])
    print("  local slope c_eff (descends toward c as eta->0):")
    for m, sl in zip(mid, slopes):
        print(f"    eta~{m:.2f}:  c_eff = {sl:.3f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11.2, 4.4))
    ax[0].plot(Z, G, "-", color="C0", lw=2, label=r"instanton $g(z)$")
    ax[0].fill_between(Z, 0, G, color="C0", alpha=0.12)
    ax[0].set_xlabel(r"$z$"); ax[0].set_ylabel(r"$g$")
    ax[0].set_title(r"Optimal fluctuation: $-g''+zg=g^3$  ($c=\!\int g'^2$)")
    ax[0].legend(fontsize=9, frameon=False); ax[0].grid(alpha=0.3)

    ax[1].plot(mid, slopes, "o-", color="C0", ms=7, label=r"measured local slope $c_{\rm eff}(\eta)$")
    ax[1].axhline(c, color="C2", lw=1.8, ls="-", label=fr"instanton $c={c:.3f}$ ($\eta\!\to\!0$)")
    ax[1].axhline(6.22, color="C1", lw=1.2, ls="--", label=r"full-range fit $6.22$ (pre-asymptotic)")
    ax[1].axhline(2 * np.pi, color="C3", lw=1.2, ls=":", label=r"$2\pi=6.28$ (rejected)")
    ax[1].set_xlabel(r"$\eta$ (smaller $\to$ asymptotic)"); ax[1].set_ylabel(r"$c_{\rm eff}=-\,d\ln P_{\rm early}/d(\eta^{-2})$")
    ax[1].set_title(r"Rate descends toward $c=5.444$ as $\eta\to0$")
    ax[1].legend(fontsize=8.5, frameon=False); ax[1].grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(figdir, "folded_cycle_rate_constant.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print(f"\n  figure -> figures/folded_cycle_rate_constant.png\n")


if __name__ == "__main__":
    main()
