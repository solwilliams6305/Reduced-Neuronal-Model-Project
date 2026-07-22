"""
peeloff_tw_validation.py — T1.1 step 0: validate the noise/peel-off machinery.
==============================================================================

Reproduces the single-unit result of the NoisyFoldedCycle paper (and
regime-tests/folded_cycle_tracy_widom.py): the inner exit measure of the noisy
folded limit cycle is Tracy–Widom_β with β = 4/η².

Inner Cole–Hopf / stochastic-Airy equation (paper Lemma 7):
    u'' = (Y - η ξ) u ,   Y = Y0 - T ,   ξ = dB/dT
A finite-time blow-up of the Riccati R = -u'/u  ⇔  a simple node (zero) of u.
Sweeping Y downward from Y0, the FIRST node level Y_node is the peel-off level,
and  Y_node  =d  -Λ0(β)  =d  TW_β   (Ramírez–Rider–Virág).

This is the baseline that the coupled/forced extension (T1.1 core, T1.2) builds on.
Method mirrors regime-tests/folded_cycle_tracy_widom.py (Heun/Stratonovich).

Output: figures/peeloff_tw_validation.png + printed moment table.
"""
from __future__ import annotations
import os, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

# Tracy–Widom_β reference moments (Bornemann 2010): mean, std, skew, excess-kurt
TW = {1: (-1.2065, 1.2680, 0.2935, 0.1652),
      2: (-1.7711, 0.9018, 0.2241, 0.0934),
      4: (-2.3069, 0.7178, 0.1655, 0.0490)}


def first_node(eta, Y0=7.0, dt=8e-4, N=16000, Ymin=-9.0, seed=0):
    """First node of u in u''=(Y-η ξ)u, Y=Y0-T (Heun/Stratonovich) = peel-off level."""
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
    m = x.mean(); d = x - m; var = np.mean(d**2)
    return m, np.sqrt(var), np.mean(d**3) / var**1.5, np.mean(d**4) / var**2 - 3.0


def main():
    t0 = time.time()
    etas = np.array([0.6, 0.8, 1.0, np.sqrt(2), 1.7, 2.0])
    print("=" * 74)
    print("T1.1 baseline — inner exit measure = Tracy–Widom_β,  β = 4/η²")
    print("=" * 74)
    print(f"  {'η':>6}{'β':>6} | {'mean':>8}{'std':>7}{'skew':>7}{'kurt':>7} |"
          f"  TW_β reference (skew/kurt)")
    betas, skews, samples = [], [], {}
    checks = []
    for k, eta in enumerate(etas):
        x = first_node(eta, seed=21 + k)
        m, s, g1, g2 = moments(x)
        b = 4 / eta**2
        betas.append(b); skews.append(g1); samples[eta] = x
        rb = round(b); ref = ""
        if abs(b - rb) < 0.05 and rb in TW:
            ref = f"  TW{rb}: {TW[rb][2]:.3f}/{TW[rb][3]:.3f}"
            dskew = abs(g1 - TW[rb][2])
            checks.append((rb, dskew, abs(m - TW[rb][0]), abs(s - TW[rb][1])))
        print(f"  {eta:6.3f}{b:6.2f} | {m:8.3f}{s:7.3f}{g1:7.3f}{g2:7.3f} |{ref}")
    betas = np.array(betas); skews = np.array(skews)

    # β=4 (GSE) carries a documented 2^{-1/6} edge-convention factor in location/scale,
    # so only its SHAPE (skew) is convention-free; β=1,2 are checked in mean/std too.
    ok = all(ds < 0.035 and (rb == 4 or (dm < 0.06 and dstd < 0.06))
             for rb, ds, dm, dstd in checks)
    print("\n  skew(β=1,2,4) + mean/std(β=1,2) vs Bornemann TW reference:",
          "PASS" if ok else "CHECK", "  [β=4 location/scale = known GSE 2^(-1/6) factor]")
    print("  skewness → 0 (Gaussian) as β → ∞ (η → 0):",
          "yes" if skews[0] < skews[-1] else "CHECK")

    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
    ax[0].semilogx(betas, skews, "o-", color="#1f3b73", ms=6, label="measured skew(Y_node)")
    for b in (1, 2, 4):
        ax[0].plot(b, TW[b][2], "*", color="#b3402b", ms=15,
                   label="Tracy–Widom_β" if b == 1 else None)
    ax[0].axhline(0.0, color="k", lw=0.8, ls=":", label="Gaussian (skew 0)")
    ax[0].set_xlabel("β = 4/η²   (η→0 to the right)")
    ax[0].set_ylabel("skewness of peel-off level")
    ax[0].set_title("skewness tracks TW_β along η = 2/√β")
    ax[0].invert_xaxis(); ax[0].legend(fontsize=9, frameon=False); ax[0].grid(alpha=0.3, which="both")

    g = np.linspace(-4, 4, 200)
    ax[1].plot(g, np.exp(-g**2 / 2) / np.sqrt(2 * np.pi), "k--", lw=1.2, label="standard Gaussian")
    for eta, col, bb in [(np.sqrt(2), "#1f3b73", 2), (2.0, "#d95f0e", 1)]:
        x = samples[eta]; z = (x - x.mean()) / x.std()
        ax[1].hist(z, bins=70, range=(-4, 4), density=True, histtype="step", lw=1.8,
                   color=col, label=f"β={bb} (η={eta:.2f}), skew {moments(x)[2]:.2f}")
    ax[1].set_xlabel("standardized Y_node"); ax[1].set_ylabel("density")
    ax[1].set_title("right-skewed early-escape tail = TW_β edge")
    ax[1].legend(fontsize=8.5, frameon=False); ax[1].grid(alpha=0.3)

    fig.suptitle("T1.1 baseline — noisy folded-cycle peel-off is Tracy–Widom_β (single unit)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fp = os.path.join(FIG, "peeloff_tw_validation.png")
    fig.savefig(fp, dpi=140, bbox_inches="tight")
    print(f"\n  Figure: {fp}")
    print(f"  TOTAL {time.time()-t0:.1f}s — DONE")


if __name__ == "__main__":
    main()
