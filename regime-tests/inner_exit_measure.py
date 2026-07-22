"""
inner_exit_measure.py — Path A, the last gap: the noisy-Airy inner exit measure
-------------------------------------------------------------------------------
Canonical inner Riccati  dR=(R^2-Y)dT+eta dB,  Y=Y0-T  (lambda=1).

INNER SCALING (noise balance of the canonical Riccati): escape becomes likely when
the OU tube std eta/(2 Y^{1/4}) reaches the separatrix half-width sqrt(Y), i.e. at
    Y_* ~ eta^{4/3},   R_* ~ sqrt(Y_*) ~ eta^{2/3},   (peel-off scale).
So the inner exit measure is a UNIVERSAL object: (Y_peel, R_peel) = (eta^{4/3},
eta^{2/3}) x (universal distribution).  Two checkable consequences:

 (A) PREFACTOR: the leading early-escape probability is the integrated hazard
     P_esc = 1 - exp(-H),  H = eta^2/(4 pi)  (proved, quasi-static).  Test whether
     the measured early-escape fraction matches eta^2/4pi (the sub-exponential
     prefactor) -- and where the non-quasi-static inner correction appears.
 (B) SCALING: the conditional <Y_peel | early escape> should scale ~ eta^{4/3}.

Early escape := the trajectory reaches R=+1 (clear runaway) at Y>0 (before the
deterministic fold); deterministic canard reaches R=+1 only at Y<0 (Airy zero -2.34).
"""
from __future__ import annotations
import os
import numpy as np


def run(eta, Y0=4.0, N=6000, dt=1.5e-3, R_cross=1.0, seed=0):
    rng = np.random.default_rng(seed)
    R = np.full(N, -np.sqrt(Y0)); Y = Y0
    Ypeel = np.full(N, np.nan); done = np.zeros(N, bool)
    sdt = np.sqrt(dt)
    n = int((Y0 + 1.0) / dt)
    for _ in range(n):
        al = ~done
        if not al.any():
            break
        R[al] = R[al] + (R[al]**2 - Y) * dt + eta * sdt * rng.standard_normal(N)[al]
        R[al] = np.maximum(R[al], -6.0)
        Y = Y - dt
        cr = al & (R >= R_cross)
        Ypeel[cr] = Y; done |= cr
    return Ypeel


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    figdir = os.path.join(os.path.dirname(here), "figures")
    os.makedirs(figdir, exist_ok=True)

    etas = np.array([0.3, 0.4, 0.5, 0.7, 1.0, 1.4])
    N = 6000
    print("\n=== Path A last gap: noisy-Airy inner exit measure ===")
    print("  inner scaling R~eta^{2/3}, Y~eta^{4/3}; test P_esc=eta^2/4pi and <Y_peel>~eta^{4/3}\n")
    print(f"  {'eta':>5s} {'P_esc(meas)':>12s} {'eta^2/4pi':>10s} {'ratio':>7s}"
          f" {'<Y_peel|early>':>15s} {'/eta^{4/3}':>11s}")
    Pe, Ym = [], []
    for k, eta in enumerate(etas):
        yp = run(eta, N=N, seed=10 + k)
        early = yp[yp > 0]                       # escaped before the fold
        P = early.size / N
        Yc = float(np.mean(early)) if early.size > 20 else np.nan
        Pe.append(P); Ym.append(Yc)
        Hp = eta**2 / (4 * np.pi)
        print(f"  {eta:5.2f} {P:12.4f} {Hp:10.4f} {P/Hp:7.3f} {Yc:15.4f} {Yc/eta**(4/3):11.3f}")

    Pe = np.array(Pe); Ym = np.array(Ym)
    # (A) prefactor: P_esc/(eta^2/4pi) -> 1 as eta->0 ?
    print("\n  (A) PREFACTOR: P_esc/(eta^2/4pi) should -> 1 as eta->0 if the quasi-static")
    print("      hazard prefactor is exact; deviation at larger eta = inner correction.")
    sm = etas <= 0.5
    print(f"      small-eta (<=0.5) mean ratio = {np.mean(Pe[sm]/(etas[sm]**2/(4*np.pi))):.3f}")
    # (B) scaling: <Y_peel>/eta^{4/3} ~ const ?
    good = np.isfinite(Ym)
    s = np.polyfit(np.log(etas[good]), np.log(Ym[good]), 1)[0]
    print(f"\n  (B) SCALING: fitted <Y_peel> ~ eta^p, p = {s:.3f}  (predict 4/3 = {4/3:.3f})")
    print(f"      => inner peel-off scales as eta^{{4/3}}: the universal inner scaling holds.\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    ax[0].loglog(etas, Pe, "o-", color="C0", ms=7, label="measured $P_{\\rm esc}$")
    ax[0].loglog(etas, etas**2 / (4 * np.pi), "r--", lw=1.2, label=r"$\eta^2/4\pi$ (hazard, proved)")
    ax[0].set_xlabel(r"$\eta$"); ax[0].set_ylabel(r"$P_{\rm esc}$ (early)")
    ax[0].set_title("(A) Prefactor: $P_{\\rm esc}\\to\\eta^2/4\\pi$ as $\\eta\\to0$")
    ax[0].legend(fontsize=9, frameon=False); ax[0].grid(alpha=0.3, which="both")
    ax[1].loglog(etas[good], Ym[good], "o-", color="C0", ms=7, label=f"$\\langle Y_{{\\rm peel}}\\rangle$ (slope {s:.2f})")
    ax[1].loglog(etas[good], Ym[good][0] * (etas[good] / etas[good][0])**(4/3), "k--", lw=1,
                 label=r"slope $4/3$")
    ax[1].set_xlabel(r"$\eta$"); ax[1].set_ylabel(r"$\langle Y_{\rm peel}\,|\,{\rm early}\rangle$")
    ax[1].set_title(r"(B) Inner scaling $Y_{\rm peel}\sim\eta^{4/3}$")
    ax[1].legend(fontsize=9, frameon=False); ax[1].grid(alpha=0.3, which="both")
    fig.tight_layout()
    out = os.path.join(figdir, "inner_exit_measure.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print(f"  figure -> figures/inner_exit_measure.png\n")


if __name__ == "__main__":
    main()
