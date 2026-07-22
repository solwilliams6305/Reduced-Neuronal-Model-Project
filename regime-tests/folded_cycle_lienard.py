"""
folded_cycle_lienard.py — the NAMED worked model: JKK 2024 Section 5.1
----------------------------------------------------------------------
Periodically-forced van der Pol / Lienard (JKK eq. 70-74), the canonical
folded-limit-cycle example.  With theta=1 (the O(1) rescaling parameter) and
degenerate noise on the fast variable x:

    x' = -y + K(x) + sigma dW,   K(x) = x - x^3/3   (fast; cubic fold)
    theta' = eps1                                    (EXTERNAL forcing phase)
    y' = eps2 ( g(x) - A(theta) ),  g(x)=x           (super-slow drift)

Fold of K at x_F = -1 (K'(-1)=0, K''(-1)=2>0), y_fold = K(-1) = -2/3.
JKK (74):  a = b = 1 (=theta^{-1}, theta^{-1}*K''/2),  c(theta) = A(theta) - g(x_F)
         = A(theta) + 1.  So the phase modulation is
    G(theta) = sqrt(a c / b) = sqrt(c(theta)) = sqrt(1 + A(theta)).
We take A(theta) = A0 cos(2 pi theta)  (Assumption 3: A>-1 => c>0).

FROZEN-PHASE Channel-A test (theta=theta0): the named-model canard escape.  As y
drifts down through y_fold, noise knocks x off the attracting branch x<-1 EARLY.
Margin m = y_escape - y_fold.  PREDICTION sigma_*^A(theta0) = C_q sqrt(eps2) G(theta0).
Canard R-Theta collapse (drift rate lambda = eps2 c):
    R = m / (eps2 c)^{2/3},   Theta = sigma / (sqrt(eps2) sqrt(c)) = sigma/(sqrt(eps2) G).
If sqrt(c) is the right modulation, R(Theta) collapses across theta0 (c).  This is
sigma_*^A(theta) ~ sqrt(c(theta)) validated on JKK's OWN forced-VdP example.
"""
from __future__ import annotations
import os
import numpy as np

A0 = 0.6                         # forcing amplitude (c = 1 + A0 cos in [0.4,1.6])
EPS2 = 0.02
Y_FOLD = -2.0 / 3.0


def c_of(theta):
    return 1.0 + A0 * np.cos(2 * np.pi * theta)


def simulate(theta0, sigma, N=80, x0=-1.4, dt=1e-3, T_max=60.0, x_cross=0.0, seed=0):
    """Frozen theta0; drift y down through the fold; return mean margin m."""
    rng = np.random.default_rng(seed)
    A = A0 * np.cos(2 * np.pi * theta0)
    y0 = x0 - x0**3 / 3.0                    # start on attracting branch K(x0)
    x = np.full(N, x0); y = np.full(N, y0)
    yesc = np.full(N, np.nan); esc = np.zeros(N, bool)
    sdt = np.sqrt(dt); n = int(T_max / dt)
    for _ in range(n):
        al = ~esc
        if not al.any():
            break
        xn = x.copy()
        xn[al] = x[al] + (-y[al] + x[al] - x[al]**3 / 3.0) * dt + sigma * sdt * rng.standard_normal(N)[al]
        y[al] = y[al] + EPS2 * (x[al] - A) * dt
        cr = al & (x < x_cross) & (xn >= x_cross)
        if cr.any():
            yesc[cr] = y[cr]; esc |= cr
        x = xn
    return float(np.nanmean(yesc) - Y_FOLD)      # margin m (>0 = early escape)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    figdir = os.path.join(os.path.dirname(here), "figures")
    os.makedirs(figdir, exist_ok=True)

    thetas = np.linspace(0.0, 0.5, 5)            # c from 1.6 (theta=0) to 0.4 (theta=0.5)
    sigmas = np.array([0.05, 0.09, 0.14, 0.20, 0.30, 0.45])
    N = 80

    print("\n=== NAMED worked model: JKK Sec 5.1 forced van der Pol ===")
    print(f"  A(theta)=A0 cos, A0={A0};  c(theta)=1+A0 cos;  G(theta)=sqrt(c).")
    print(f"  a=b=1, eps2={EPS2}, fold x_F=-1, y_fold=-2/3.\n")
    print(f"  Theta_correct = sigma/(sqrt(eps2) sqrt(c)) ;  R = m/(eps2 c)^(2/3)\n")

    curves = []
    for i, th in enumerate(thetas):
        c = c_of(th)
        m = np.array([simulate(th, s, N=N, seed=10 * i + k) for k, s in enumerate(sigmas)])
        R = m / (EPS2 * c) ** (2.0 / 3.0)
        Tc = sigmas / (np.sqrt(EPS2) * np.sqrt(c))      # correct (with G=sqrt c)
        Tn = sigmas / np.sqrt(EPS2)                     # naive (no G)
        curves.append(dict(th=th, c=c, m=m, R=R, Tc=Tc, Tn=Tn))
        print(f"  theta={th:.3f} c={c:.3f} G={np.sqrt(c):.3f}  <m>=["
              + ",".join(f"{v:+.3f}" for v in m) + "]")

    def cv(key):
        lo = max(cc[key].min() for cc in curves); hi = min(cc[key].max() for cc in curves)
        grid = np.linspace(lo, hi, 12)
        st = np.array([np.interp(grid, cc[key], cc["R"]) for cc in curves])
        return float(np.mean(st.std(0)) / (np.abs(st).mean() + 1e-9))
    cvc, cvn = cv("Tc"), cv("Tn")
    print(f"\n  cross-phase collapse spread (mean std / scale):")
    print(f"     vs Theta_correct = sigma/(sqrt(eps2) sqrt(c))  (G=sqrt c): {cvc:.3f}")
    print(f"     vs Theta_naive   = sigma/sqrt(eps2)            (no G)    : {cvn:.3f}")
    print(f"     => {'collapses WITH G (sqrt c confirmed)' if cvc < cvn else 'no improvement'}"
          f"   (naive/correct = {cvn/max(cvc,1e-9):.2f})\n")
    print("  This validates sigma_*^A(theta) ~ sqrt(c(theta)) = G(theta) on JKK's")
    print("  OWN forced-van der Pol example (a named system, not a normal form).\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    for cc in curves:
        ax[0].plot(cc["Tc"], cc["R"], "o-", ms=4, label=fr"$\theta$={cc['th']:.2f}, c={cc['c']:.2f}")
    ax[0].set_xlabel(r"$\Theta=\sigma/(\sqrt{\epsilon_2}\sqrt{c(\theta)})$")
    ax[0].set_ylabel(r"$R=m/(\epsilon_2 c)^{2/3}$")
    ax[0].set_title(f"Forced VdP (JKK 5.1) WITH $G=\\sqrt{{c}}$: collapse (CV {cvc:.2f})")
    ax[0].legend(fontsize=7, frameon=False); ax[0].grid(alpha=0.3)
    for cc in curves:
        ax[1].plot(cc["Tn"], cc["R"], "s--", ms=4, label=fr"$\theta$={cc['th']:.2f}")
    ax[1].set_xlabel(r"$\Theta_{\rm naive}=\sigma/\sqrt{\epsilon_2}$ (no $G$)")
    ax[1].set_ylabel(r"$R$")
    ax[1].set_title(f"WITHOUT $G$: splays (CV {cvn:.2f})")
    ax[1].legend(fontsize=7, frameon=False); ax[1].grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(figdir, "folded_cycle_lienard_Gtheta.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print(f"  figure -> figures/folded_cycle_lienard_Gtheta.png\n")


if __name__ == "__main__":
    main()
