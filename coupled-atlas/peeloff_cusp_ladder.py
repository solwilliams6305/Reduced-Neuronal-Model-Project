"""
peeloff_cusp_ladder.py — probe of the catastrophe-ladder conjecture.
====================================================================

The paper's fold has a LINEAR turning point in the inner Cole–Hopf equation,
    u'' = (Y - η ξ) u ,                         (k = 1, Airy → Tracy–Widom)
and the cusp / higher catastrophes correspond to HIGHER-ORDER turning points,
    u'' = (sign(Y)|Y|^k - η ξ) u ,              (k > 1)
which is exactly the "higher-order Airy" hierarchy that governs MULTICRITICAL
random-matrix edges (k=1 generic Airy/TW; k=2,3,... multicritical edges with
distinct, higher-order Tracy–Widom laws).

Question: does the NOISE-INDUCED peel-off law (first node of the swept recessive
solution) move off Tracy–Widom as the turning order k departs from 1?
  • if YES → a ladder of noise-induced escape universality classes exists
    (fold→Airy→TW is rung one; higher catastrophes are new rungs).
  • if NO  → TW is super-universal across turning orders (a strong negative).

This is the minimal controlled model of the cusp→Pearcey idea: k is a knob for the
singularity order. (The true coupled cusp requires the coupled blow-up; this isolates
the universality question.) Fixed η=√2 (β=2 at the fold). Heun/Stratonovich, as in
peeloff_tw_validation.py. Independent noise per k.

Output: figures/peeloff_cusp_ladder.png + printed summary.
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

TW2 = dict(skew=0.2241, kurt=0.0934)      # Bornemann reference (fold, k=1, β=2)


def peeloff_ladder(ks, eta, N=6000, Y0=5.0, dt=4e-4, Ymin=-7.0, seed=0):
    """First node of u'' = (sign(Y)|Y|^k - η ξ)u for each k. Returns (K,N) peel-off levels."""
    K = len(ks); kk = np.asarray(ks)[:, None]
    rng = np.random.default_rng(seed); sdt = np.sqrt(dt)
    Vf = lambda Y: np.sign(Y) * np.abs(Y) ** kk                  # (K,1)
    u = np.ones((K, N)); v = np.sqrt(np.maximum(Vf(Y0), 0)) * np.ones((K, N)); Y = Y0
    Yz = np.full((K, N), np.nan); done = np.zeros((K, N), bool)
    n = int((Y0 - Ymin) / dt)
    for _ in range(n):
        al = ~done
        if not al.any():
            break
        dB = sdt * rng.standard_normal((K, N))
        Vc = Vf(Y); Yp = Y - dt; Vp = Vf(Yp)
        u1 = u + v * dt
        v1 = v + (Vc * u) * dt - eta * u * dB
        u = u + 0.5 * (v + v1) * dt
        v = v + 0.5 * (Vc * u + Vp * u1) * dt - eta * 0.5 * (u + u1) * dB
        Y = Yp
        cr = al & (u < 0.0); Yz[cr] = Y; done |= cr
    return Yz


def moments(x):
    x = x[np.isfinite(x)]; m = x.mean(); d = x - m; vv = np.mean(d**2)
    return m, np.sqrt(vv), np.mean(d**3) / vv**1.5, np.mean(d**4) / vv**2 - 3.0, x.size


def main():
    t0 = time.time()
    ks = np.array([0.7, 1.0, 1.5, 2.0, 3.0])
    eta = np.sqrt(2.0)
    Yz = peeloff_ladder(ks, eta)
    N = Yz.shape[1]
    se_skew = np.sqrt(6.0 / N)                  # ~MC standard error on skewness

    print("=" * 72)
    print("Catastrophe-ladder probe — peel-off law vs turning order k (η=√2)")
    print("=" * 72)
    print(f"  fold reference (k=1): TW₂ skew {TW2['skew']:.3f}, exkurt {TW2['kurt']:.3f}"
          f"   (skew MC s.e. ≈ {se_skew:.3f})")
    print(f"  {'k':>5} | {'mean':>8}{'std':>8}{'skew':>8}{'kurt':>8} {'esc%':>6}")
    rows = []
    for i, k in enumerate(ks):
        m, s, g1, g2, cnt = moments(Yz[i])
        rows.append((k, m, s, g1, g2));
        print(f"  {k:5.2f} | {m:8.3f}{s:8.3f}{g1:8.3f}{g2:8.3f} {100*cnt/N:6.0f}")
    skews = np.array([r[3] for r in rows]); kurts = np.array([r[4] for r in rows])

    i1 = int(np.argmin(np.abs(ks - 1.0)))
    tw_ok = abs(skews[i1] - TW2['skew']) < 3 * se_skew
    moved = (abs(skews[-1] - skews[i1]) > 3 * se_skew)
    trend = np.all(np.diff(skews) > -se_skew)   # roughly monotone in k
    print(f"\n  k=1 reproduces TW₂ skew: {'PASS' if tw_ok else 'CHECK'} "
          f"({skews[i1]:.3f} vs {TW2['skew']:.3f})")
    print(f"  law MOVES off TW with k (skew k=3 vs k=1 differ > 3·s.e.): "
          f"{'YES — ladder exists' if moved else 'no — TW looks super-universal'}")
    print(f"  skew trend roughly monotone in k: {'yes' if trend else 'non-monotone'}")
    print(f"  skew range over k: {skews.min():.3f} … {skews.max():.3f}")

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
    ax[0].axhline(TW2['skew'], color="#b3402b", lw=1, ls="--", label="TW₂ (fold)")
    ax[0].fill_between([ks.min(), ks.max()], TW2['skew']-3*se_skew, TW2['skew']+3*se_skew,
                       color="#b3402b", alpha=0.12)
    ax[0].plot(ks, skews, "o-", color="#1f3b73", ms=6)
    ax[0].axvline(1.0, color="grey", lw=0.6, ls=":")
    ax[0].set_xlabel("turning order  k  (k=1 = fold)"); ax[0].set_ylabel("peel-off skewness")
    ax[0].set_title("(A) does the edge law move off TW with k?")
    ax[0].legend(fontsize=9, frameon=False)

    ax[1].axhline(TW2['kurt'], color="#b3402b", lw=1, ls="--", label="TW₂ (fold)")
    ax[1].plot(ks, kurts, "s-", color="#2c7d59", ms=6)
    ax[1].axvline(1.0, color="grey", lw=0.6, ls=":")
    ax[1].set_xlabel("turning order  k"); ax[1].set_ylabel("peel-off excess kurtosis")
    ax[1].set_title("(B) excess kurtosis vs k")
    ax[1].legend(fontsize=9, frameon=False)

    cols = ["#888780", "#b3402b", "#1f3b73", "#7a3b8f"]
    for j, k in enumerate([0.7, 1.0, 2.0, 3.0]):
        idx = int(np.argmin(np.abs(ks - k))); x = Yz[idx][np.isfinite(Yz[idx])]
        z = (x - x.mean()) / x.std()
        ax[2].hist(z, bins=70, range=(-4, 4), density=True, histtype="step", lw=1.7,
                   color=cols[j], label=f"k={k} (skew {moments(Yz[idx])[2]:.2f})")
    g = np.linspace(-4, 4, 200)
    ax[2].plot(g, np.exp(-g**2/2)/np.sqrt(2*np.pi), "k:", lw=1, label="Gaussian")
    ax[2].set_xlabel("standardized peel-off"); ax[2].set_ylabel("density")
    ax[2].set_title("(C) standardized laws across k")
    ax[2].legend(fontsize=8, frameon=False)

    fig.suptitle("Catastrophe-ladder probe — noise-induced escape law vs turning order "
                 "(fold k=1 → TW; higher k = cusp/multicritical)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fp = os.path.join(FIG, "peeloff_cusp_ladder.png")
    fig.savefig(fp, dpi=140, bbox_inches="tight")
    print(f"\n  Figure: {fp}")
    print(f"  TOTAL {time.time()-t0:.1f}s — DONE")


if __name__ == "__main__":
    main()
