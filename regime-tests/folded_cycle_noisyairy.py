"""
folded_cycle_noisyairy.py — Path A last gap: the noisy-Airy inner exit measure
-------------------------------------------------------------------------------
The inner Riccati  dR=(R^2-Y)dT+eta dB,  Y=Y0-T  has an EXACT Cole-Hopf
linearisation  R=-u'/u  =>  u'' = (Y - eta*xi) u   (linear Schrodinger /
stochastic Airy equation).  Escape (R->+inf) <=> first node of u.  So the inner
exit measure = first-node law of a LINEAR stochastic ODE, which is the eigenvalue
problem of the STOCHASTIC AIRY OPERATOR  H = -d^2/dx^2 + x + (2/sqrt(beta)) b'(x)
with the dictionary  x = -Y,  eta = 2/sqrt(beta)  (Ramirez-Rider-Virag 2011).

This script verifies the two checkable consequences of that identification:

  (BULK)  typical first node = a1 + O(eta) GAUSSIAN  (1st-order perturbation of
          the deterministic Airy ground state).  Predict std ~ eta^1, mean -> a1
          = -2.3381.  [large-beta / Gaussian edge of the SAO]

  (TAIL)  early escape (first node at Y>0 <=> SAO ground state lambda0<0) is a
          left large-deviation; the quasi-static hazard H=eta^2/4pi OVERestimates
          it -- the true (smaller) prefactor is the SAO ground-state left tail.
          [this is the open sub-exponential prefactor; here we just exhibit the
          two-scale split that the SAO picture predicts.]

Companion: FOLDED_CYCLE_NOISYAIRY_CONJECTURES.md
"""
from __future__ import annotations
import os
import numpy as np


def first_node(eta, Y0=4.0, dt=5e-4, N=8000, Ymin=-4.0, seed=0):
    """March u'' = (Y-eta xi) u from the decaying Airy branch; return Y at the
    first node of u (Heun/Stratonovich).  Start u=1, u'=+sqrt(Y0) so that
    R=-u'/u=-sqrt(Y0) sits on the attracting canard."""
    rng = np.random.default_rng(seed)
    u = np.ones(N)
    v = np.full(N, np.sqrt(Y0))
    Y = Y0
    sdt = np.sqrt(dt)
    Yz = np.full(N, np.nan)
    done = np.zeros(N, bool)
    n = int((Y0 - Ymin) / dt)
    for _ in range(n):
        al = ~done
        if not al.any():
            break
        dB = sdt * rng.standard_normal(N)
        # u' = v ;  v' = (Y - eta xi) u  ->  dv = Y u dt - eta u dB
        u1 = u + v * dt
        Yp = Y - dt
        v1 = v + (Y * u) * dt - eta * u * dB
        u = u + 0.5 * (v + v1) * dt
        v = v + 0.5 * (Y * u + Yp * u1) * dt - eta * (0.5 * (u + u1)) * dB
        Y = Yp
        cr = al & (u < 0.0)            # u started > 0 : first sign change = first node
        Yz[cr] = Y
        done |= cr
    return Yz


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    figdir = os.path.join(os.path.dirname(here), "figures")
    os.makedirs(figdir, exist_ok=True)
    a1 = -2.33810741   # first zero of Ai

    print("\n=== Path A last gap: noisy-Airy = stochastic Airy operator (eta=2/sqrt(beta)) ===")
    print("  (BULK) first node = a1 + O(eta) Gaussian: predict std~eta^1, mean->a1=-2.3381\n")
    print(f"  {'eta':>6} {'beta=4/eta^2':>12} {'mean Yz':>9} {'std Yz':>8} {'std/eta':>8} {'P(early)':>9}")
    etas = np.array([0.10, 0.15, 0.20, 0.30, 0.40])
    means, stds = [], []
    for k, eta in enumerate(etas):
        Yz = first_node(eta, seed=7 + k)
        m, s = np.nanmean(Yz), np.nanstd(Yz)
        pe = np.mean(Yz > 0.0)
        means.append(m); stds.append(s)
        print(f"  {eta:6.2f} {4/eta**2:12.1f} {m:9.4f} {s:8.4f} {s/eta:8.3f} {pe:9.4f}")
    means, stds = np.array(means), np.array(stds)
    p = np.polyfit(np.log(etas), np.log(stds), 1)[0]
    print(f"\n  fitted std ~ eta^p, p = {p:.3f}  (predict 1.000: Gaussian perturbative bulk)")
    print(f"  std/eta -> {stds[0]/etas[0]:.3f} (= sqrt(integral psi0^4), an Airy constant)")
    print(f"  mean(eta=0.1) = {means[0]:.4f}  vs  a1 = {a1:.4f}\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    # (left) bulk std ~ eta^1
    ax[0].loglog(etas, stds, "o-", color="C0", ms=7, label=f"std of first node (slope {p:.2f})")
    ax[0].loglog(etas, (stds[0] / etas[0]) * etas, "k--", lw=1, label="slope 1 (O($\\eta$) Gaussian)")
    ax[0].set_xlabel(r"$\eta\;(=2/\sqrt{\beta})$")
    ax[0].set_ylabel(r"std$(Y_{\rm node})$")
    ax[0].set_title(r"(BULK) SAO edge: $Y_{\rm node}=a_1+O(\eta)$ Gaussian")
    ax[0].legend(fontsize=9, frameon=False); ax[0].grid(alpha=0.3, which="both")
    # (right) histogram of the centred first node at one eta vs a1
    eta = 0.30
    Yz = first_node(eta, N=20000, seed=99)
    ax[1].hist(Yz, bins=60, density=True, color="C0", alpha=0.8, label=fr"$Y_{{\rm node}}$, $\eta$={eta}")
    ax[1].axvline(a1, color="C3", lw=1.5, ls="--", label=r"$a_1=-2.338$ (det. Airy zero)")
    ax[1].axvline(0.0, color="k", lw=1.0, ls=":", label=r"fold $Y=0$ (early-escape edge)")
    ax[1].set_xlabel(r"$Y_{\rm node}$ (first node $=$ SAO ground state $-\lambda_0$)")
    ax[1].set_ylabel("density")
    ax[1].set_title(r"Inner exit measure $=$ SAO ground-state law")
    ax[1].legend(fontsize=8, frameon=False); ax[1].grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(figdir, "folded_cycle_noisyairy.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print(f"  figure -> figures/folded_cycle_noisyairy.png\n")


if __name__ == "__main__":
    main()
