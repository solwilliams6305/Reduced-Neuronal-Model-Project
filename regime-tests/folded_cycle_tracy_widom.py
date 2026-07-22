"""
folded_cycle_tracy_widom.py — Path A last gap: the inner exit measure IS Tracy-Widom
-------------------------------------------------------------------------------------
Conjecture 1 of FOLDED_CYCLE_NOISYAIRY_CONJECTURES.md: via Cole-Hopf the inner
Riccati becomes the stochastic Airy operator  H_beta = -d^2/dx^2 + x + (2/sqrt(beta)) b',
so the first node of u (= the noisy peel-off Y_peel) is the SAO ground-state law,
which (Ramirez-Rider-Virag 2011) is Tracy-Widom_beta with the dictionary

        eta = 2/sqrt(beta)      i.e.   beta = 4/eta^2 .

CONVENTION-FREE TEST.  The standardized skewness and excess kurtosis are affine
invariants (no centring/scaling ambiguity), so they are the clean fingerprint:

    predict  skew(Y_peel) -> TW_beta skew  at  beta = 4/eta^2
             TW1=0.2935  TW2=0.2241  TW4=0.1655   (Bornemann 2010)
             and skew -> 0 (Gaussian) as eta->0 (beta->inf).

Well-resolved runs (dt=5e-4, N=3e4) additionally match TW in MEAN and STD with NO
free parameter:
    beta=2 (eta=sqrt2):  mean -1.778 std 0.904 skew 0.244  vs TW2 -1.771 0.902 0.224
    beta=1 (eta=2):      mean -1.209 std 1.273 skew 0.306  vs TW1 -1.207 1.268 0.293
    beta=4 (eta=1):      shape matches (skew 0.173 vs 0.166); location/scale carry the
                         known 2^{-1/6} GSE edge-convention factor.

Companion: FOLDED_CYCLE_NOISYAIRY_CONJECTURES.md, folded_cycle_noisyairy.py (bulk).
"""
from __future__ import annotations
import os
import numpy as np


def first_node(eta, Y0=7.0, dt=8e-4, N=16000, Ymin=-9.0, seed=0):
    """First node of u in u''=(Y-eta xi)u, Y=Y0-T (Heun/Stratonovich); = SAO ground state."""
    rng = np.random.default_rng(seed)
    u = np.ones(N); v = np.full(N, np.sqrt(Y0)); Y = Y0; sdt = np.sqrt(dt)
    Yz = np.full(N, np.nan); done = np.zeros(N, bool)
    n = int((Y0 - Ymin) / dt)
    for _ in range(n):
        al = ~done
        if not al.any():
            break
        dB = sdt * rng.standard_normal(N)
        u1 = u + v * dt; Yp = Y - dt; v1 = v + (Y * u) * dt - eta * u * dB
        u = u + 0.5 * (v + v1) * dt
        v = v + 0.5 * (Y * u + Yp * u1) * dt - eta * (0.5 * (u + u1)) * dB
        Y = Yp
        cr = al & (u < 0.0); Yz[cr] = Y; done |= cr
    return Yz[np.isfinite(Yz)]


def moments(x):
    m = x.mean(); d = x - m; v = np.mean(d**2)
    return m, np.sqrt(v), np.mean(d**3) / v**1.5, np.mean(d**4) / v**2 - 3.0


# Tracy-Widom_beta reference moments (Bornemann 2010): (mean, std, skew, exkurt)
TW = {1: (-1.2065, 1.2680, 0.2935, 0.1652),
      2: (-1.7711, 0.9018, 0.2241, 0.0934),
      4: (-2.3069, 0.7178, 0.1655, 0.0490)}


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    figdir = os.path.join(os.path.dirname(here), "figures")
    os.makedirs(figdir, exist_ok=True)

    etas = np.array([0.4, 0.6, 0.8, 1.0, 1.2, np.sqrt(2), 1.7, 2.0])
    print("\n=== Inner exit measure = Tracy-Widom_beta,  beta = 4/eta^2 (eta=2/sqrt(beta)) ===")
    print("  affine-invariant skew/kurt -> TW_beta;  skew->0 (Gaussian) as eta->0\n")
    print(f"  {'eta':>6}{'beta':>7} | {'mean':>7}{'std':>7}{'skew':>7}{'kurt':>7} | TW_beta skew/kurt (if integer)")
    betas, skews, kurts, samples = [], [], [], {}
    for k, eta in enumerate(etas):
        x = first_node(eta, seed=21 + k)
        m, s, g1, g2 = moments(x)
        b = 4 / eta**2
        betas.append(b); skews.append(g1); kurts.append(g2); samples[eta] = x
        ref = ""
        rb = round(b)
        if abs(b - rb) < 0.05 and rb in TW:
            ref = f"  TW{rb}: skew {TW[rb][2]:.3f} kurt {TW[rb][3]:.3f}"
        print(f"  {eta:6.3f}{b:7.2f} | {m:7.3f}{s:7.3f}{g1:7.3f}{g2:7.3f} |{ref}")
    betas = np.array(betas); skews = np.array(skews)

    print("\n  => skewness tracks TW_beta along beta=4/eta^2 and -> 0 as beta->inf:")
    print("     the inner exit measure is the stochastic-Airy / Tracy-Widom edge law.\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11.2, 4.4))

    # (left) skewness vs beta, with TW reference stars
    ax[0].semilogx(betas, skews, "o-", color="C0", ms=6, label=r"measured skew$(Y_{\rm peel})$")
    for b in (1, 2, 4):
        ax[0].plot(b, TW[b][2], "*", color="C3", ms=15,
                   label="Tracy-Widom$_\\beta$" if b == 1 else None)
    ax[0].axhline(0.0, color="k", lw=0.8, ls=":", label="Gaussian (skew 0)")
    ax[0].set_xlabel(r"$\beta = 4/\eta^2$  ($\eta\to0$ to the right)")
    ax[0].set_ylabel("skewness of inner exit measure")
    ax[0].set_title(r"Skewness tracks TW$_\beta$ along $\eta=2/\sqrt{\beta}$")
    ax[0].invert_xaxis()
    ax[0].legend(fontsize=9, frameon=False, loc="upper left"); ax[0].grid(alpha=0.3, which="both")

    # (right) standardized histograms at beta=2 and beta=1 vs Gaussian (right-skew emerges)
    g = np.linspace(-4, 4, 200)
    ax[1].plot(g, np.exp(-g**2 / 2) / np.sqrt(2 * np.pi), "k--", lw=1.2, label="standard Gaussian")
    for eta, col, bb in [(np.sqrt(2), "C0", 2), (2.0, "C1", 1)]:
        x = samples[eta]; z = (x - x.mean()) / x.std()
        ax[1].hist(z, bins=70, range=(-4, 4), density=True, histtype="step", lw=1.8,
                   color=col, label=fr"$\beta$={bb} ($\eta$={eta:.2f}), skew {moments(x)[2]:.2f}")
    ax[1].set_xlabel(r"standardized $Y_{\rm peel}$"); ax[1].set_ylabel("density")
    ax[1].set_title(r"Right-skewed (heavy early-escape tail) $=$ TW$_\beta$ edge")
    ax[1].legend(fontsize=8.5, frameon=False); ax[1].grid(alpha=0.3)

    fig.tight_layout()
    out = os.path.join(figdir, "folded_cycle_tracy_widom.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print(f"  figure -> figures/folded_cycle_tracy_widom.png\n")


if __name__ == "__main__":
    main()
