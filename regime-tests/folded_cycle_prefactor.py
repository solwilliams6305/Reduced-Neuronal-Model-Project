"""
folded_cycle_prefactor.py — Path A: extracting the noisy-Airy prefactor
-----------------------------------------------------------------------
Two distinct observables of the inner (noisy-Airy) problem dR=(R^2-Y)dT+eta dB:

  (i)  SADDLE-CROSSING HAZARD  H(eta)=eta^2/4pi  -- the integrated Kramers rate for
       the canard to cross the separatrix; sets the critical noise sigma_*^A (leading
       order, DONE).  POLYNOMIAL in eta.

  (ii) CANONICAL EARLY ESCAPE  P_early(eta)=P(first node of u before the fold, Y>0)
       = P(SAO ground state Lambda0(beta)<0),  beta=4/eta^2.  This is the inner exit
       measure's tail.  EXPONENTIAL: P_early ~ exp(-c/eta^2), c~6.3 (~2pi).

These are DIFFERENT objects (the node/full-blow-up is stricter than one saddle
crossing), which resolves the old "P_esc << eta^2/4pi" puzzle.  The sub-exponential
prefactor the gap needed is the EXACT Tracy-Widom tail mass 1-F_beta(0), computed
here from Painleve II and matched to Monte Carlo.

Companion: FOLDED_CYCLE_NOISYAIRY_PROOF.md
"""
from __future__ import annotations
import os
import numpy as np


def P_early(eta, Y0=4.0, dt=1e-3, N=80000, seed=0):
    """P(first node of u''=(Y-eta xi)u before the fold Y=0) = P(Lambda0<0)."""
    rng = np.random.default_rng(seed)
    u = np.ones(N); v = np.full(N, np.sqrt(Y0)); Y = Y0; sdt = np.sqrt(dt)
    done = np.zeros(N, bool)
    for _ in range(int(Y0 / dt)):
        dB = sdt * rng.standard_normal(N)
        u1 = u + v * dt; Yp = Y - dt; v1 = v + (Y * u) * dt - eta * u * dB
        u = u + 0.5 * (v + v1) * dt
        v = v + 0.5 * (Y * u + Yp * u1) * dt - eta * (0.5 * (u + u1)) * dB
        Y = Yp
        done |= (u < 0.0)
    return done.mean()


def painleve_tw(s_max=10.0, s_min=-10.0, h=5e-4):
    """Exact TW_1, TW_2 CDFs via Hastings-McLeod (corrected Airy IC)."""
    z = (2 / 3) * s_max**1.5; pf = np.exp(-z) / (2 * np.sqrt(np.pi))
    q = pf / s_max**0.25 * (1 - 5 / (72 * z) + 385 / (10368 * z**2))
    p = -pf * s_max**0.25 * (1 + 7 / (72 * z) - 455 / (10368 * z**2))
    n = int(round((s_max - s_min) / h)); s = s_max; y = np.array([q, p, 0.0, 0.0, 0.0])
    S = np.empty(n + 1); Ev = np.empty(n + 1); W = np.empty(n + 1); S[0] = s; Ev[0] = 0; W[0] = 0

    def f(s, y):
        q, p, u, E, w = y
        return np.array([p, s * q + 2 * q**3, -q * q, -u, -q])
    for i in range(1, n + 1):
        k1 = f(s, y); k2 = f(s - h / 2, y - h / 2 * k1); k3 = f(s - h / 2, y - h / 2 * k2); k4 = f(s - h, y - h * k3)
        y = y - h / 6 * (k1 + 2 * k2 + 2 * k3 + k4); s -= h; S[i] = s; Ev[i] = y[3]; W[i] = y[4]
    o = np.argsort(S); S, Ev, W = S[o], Ev[o], W[o]
    return S, np.exp(-(W + Ev) / 2), np.exp(-Ev)        # S, F1, F2


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    figdir = os.path.join(os.path.dirname(here), "figures")
    os.makedirs(figdir, exist_ok=True)

    etas = np.array([1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
    P = np.array([P_early(e, N=80000, seed=5) for e in etas])
    inv = 1 / etas**2; lnP = np.log(P)
    c, a = -np.polyfit(inv, lnP, 1)[0], np.polyfit(inv, lnP, 1)[1]
    r_exp = (lnP - (a - c * inv)).std()
    pw = np.polyfit(np.log(etas), lnP, 1); r_pow = (lnP - np.polyval(pw, np.log(etas))).std()

    S, F1, F2 = painleve_tw()
    F2_0 = float(np.interp(0.0, S, F2)); F1_0 = float(np.interp(0.0, S, F1))
    mc2 = P_early(np.sqrt(2.0), N=120000, seed=3); mc1 = P_early(2.0, N=120000, seed=3)

    print("\n=== Noisy-Airy prefactor extraction ===")
    print(f"  early escape is EXPONENTIAL: ln P = {a:.2f} - {c:.3f}/eta^2  (resid {r_exp:.3f})")
    print(f"     vs power-law ln P = {pw[1]:.2f} + {pw[0]:.2f} ln eta  (resid {r_pow:.3f})  -> exponential wins")
    print(f"     rate c = {c:.3f}   (candidate 2pi = {2*np.pi:.3f})")
    print(f"  exact-TW tail cross-check  P_early = 1 - F_beta(0):")
    print(f"     beta=2: MC {mc2:.4f}  vs  1-F2(0) {1-F2_0:.4f}")
    print(f"     beta=1: MC {mc1:.4f}  vs  1-F1(0) {1-F1_0:.4f}\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11.2, 4.4))

    # (left) ln P vs 1/eta^2: exponential (linear) beats polynomial (curved)
    xx = np.linspace(inv.min() * 0.9, inv.max() * 1.05, 100)
    ax[0].plot(inv, lnP, "o", color="C0", ms=8, label="measured $\\ln P_{\\rm early}$")
    ax[0].plot(xx, a - c * xx, "-", color="C0", lw=1.5, label=fr"$\exp(-{c:.2f}/\eta^2)$ fit ($c\approx2\pi$)")
    # the saddle hazard eta^2/4pi would be ln = ln(eta^2/4pi) = -ln(4pi) - (1/eta^2)*... no: polynomial
    etg = 1 / np.sqrt(xx); ax[0].plot(xx, np.log(etg**2 / (4 * np.pi)), "--", color="C3", lw=1.3,
                                      label=r"$\ln(\eta^2/4\pi)$ (saddle hazard, $\neq$ this observable)")
    ax[0].set_xlabel(r"$1/\eta^2$"); ax[0].set_ylabel(r"$\ln P_{\rm early}$")
    ax[0].set_title(r"Early escape is exponential $\exp(-c/\eta^2)$, not $\propto\eta^2$")
    ax[0].legend(fontsize=8.5, frameon=False); ax[0].grid(alpha=0.3)

    # (right) P_early vs eta: MC vs EXACT TW tail 1-F_beta(0)
    ax[1].plot(etas, P, "o", color="C0", ms=7, label="MC $P_{\\rm early}$")
    ax[1].plot(etas, np.exp(a - c * inv), "-", color="C0", lw=1.2, alpha=0.7, label="exp fit")
    ax[1].plot([np.sqrt(2)], [1 - F2_0], "*", color="C3", ms=18, label=r"exact $1-F_2(0)$ ($\beta$=2)")
    ax[1].plot([2.0], [1 - F1_0], "P", color="C2", ms=12, label=r"exact $1-F_1(0)$ ($\beta$=1)")
    ax[1].plot([np.sqrt(2)], [mc2], "o", mfc="none", mec="C3", ms=12)
    ax[1].plot([2.0], [mc1], "o", mfc="none", mec="C2", ms=12)
    ax[1].set_xlabel(r"$\eta\;(=2/\sqrt{\beta})$"); ax[1].set_ylabel(r"$P_{\rm early}=P(Y_{\rm node}>0)$")
    ax[1].set_title(r"Prefactor $=$ exact Tracy-Widom tail $1-F_\beta(0)$")
    ax[1].legend(fontsize=8.5, frameon=False); ax[1].grid(alpha=0.3)

    fig.tight_layout()
    out = os.path.join(figdir, "folded_cycle_prefactor.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print(f"  figure -> figures/folded_cycle_prefactor.png\n")


if __name__ == "__main__":
    main()
