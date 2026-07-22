"""
coupled_peeloff.py — T1.1 core: peel-off statistics of TWO coupled folded cycles.
=================================================================================

Inner (blow-up chart K2) amplitude equation of one folded cycle (paper eq. 1):
    dR = (R^2 - Y) dT + η dB ,   Y = Y0 - T ,
with the noise-induced peel-off level Y_node (R → +∞) distributed as TW_β, β=4/η².

Electrical (gap-junction) coupling g(v2 - v1) of the two physical fast variables
is, in the rescaling chart, a coupling of the inner AMPLITUDES R1, R2 (r = ε2^{1/3} R).
So the coupled inner system is:
    dR_i = (R_i^2 - Y + g (R_j - R_i)) dT + η dB_i ,   i=1,2 (independent noise).

Numerics / escape: R blows up at peel-off, so we integrate with an amplitude cap
R_esc and record Y_node when R_i first reaches it; the escaped unit is then FROZEN
at R_esc — which keeps it pulling its partner (the dominant spike-mediated coupling
term identified by Chow–Kopell / Lewis–Rinzel). Euler–Maruyama, vectorized over g.

NB on faithfulness: this is the inner-amplitude image of gap-junction coupling; the
exact blow-down of the constant g is a derivation TODO, and the full coupled FHR
(coupled_fhr.py) is the physical cross-check. What is robust here is the QUALITATIVE
law: how coupling correlates and shifts the two peel-off distributions — the coupled
extension of the single-unit TW_β toward the §9.4 coupled joint law.

Output: figures/coupled_peeloff.png + printed summary.
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

TW2_SKEW = 0.2241   # Bornemann reference (β=2)


def coupled_peeloff(eta, g_arr, N=3000, Y0=7.0, dt=8e-4, Ymin=-6.5,
                    R_esc=3.0, seed=0):
    """Return Yn (G,N,2): peel-off level of each unit, per coupling g."""
    rng = np.random.default_rng(seed)
    G = len(g_arr)
    g = np.asarray(g_arr, float)[:, None, None]
    R = np.full((G, N, 2), -np.sqrt(Y0))
    Y = Y0; sdt = np.sqrt(dt)
    esc = np.zeros((G, N, 2), bool)
    Yn = np.full((G, N, 2), np.nan)
    n = int((Y0 - Ymin) / dt)
    for _ in range(n):
        active = ~esc
        noise = eta * sdt * rng.standard_normal((G, N, 2))
        coup = g * (R[:, :, ::-1] - R)
        dR = (R**2 - Y + coup) * dt + noise
        R = np.where(active, R + dR, R)
        R = np.clip(R, -6.0, R_esc)
        Y -= dt
        newesc = active & (R >= R_esc)
        Yn[newesc] = Y
        esc |= newesc
        if esc.all():
            break
    return Yn


def skew(x):
    d = x - x.mean(); v = np.mean(d**2)
    return np.mean(d**3) / v**1.5 if v > 0 else np.nan


def main():
    t0 = time.time()
    eta = np.sqrt(2.0)                              # β = 2
    g_arr = np.array([0.0, 0.05, 0.1, 0.2, 0.35, 0.6, 1.0])
    Yn = coupled_peeloff(eta, g_arr, seed=3)

    print("=" * 72)
    print(f"T1.1 core — two coupled folded cycles, η=√2 (β=2)")
    print("=" * 72)
    print(f"  {'g':>6} | {'corr(Y1,Y2)':>12} {'marg.mean':>10} {'marg.std':>9} {'marg.skew':>10}")
    corrs, means, stds, skews = [], [], [], []
    for gi, gg in enumerate(g_arr):
        Y1 = Yn[gi, :, 0]; Y2 = Yn[gi, :, 1]
        ok = np.isfinite(Y1) & np.isfinite(Y2)
        c = float(np.corrcoef(Y1[ok], Y2[ok])[0, 1])
        marg = np.concatenate([Y1[ok], Y2[ok]])
        corrs.append(c); means.append(marg.mean()); stds.append(marg.std()); skews.append(skew(marg))
        print(f"  {gg:6.3f} | {c:12.3f} {marg.mean():10.3f} {marg.std():9.3f} {skews[-1]:10.3f}")
    corrs = np.array(corrs)

    g0_ok = abs(corrs[0]) < 0.10
    g0_tw = abs(skews[0] - TW2_SKEW) < 0.06
    mono = np.all(np.diff(corrs) > -0.03)
    print(f"\n  g=0 consistency: corr≈0 ({corrs[0]:+.3f}) {'PASS' if g0_ok else 'CHECK'};  "
          f"marginal skew≈TW2 ({skews[0]:.3f} vs {TW2_SKEW}) {'PASS' if g0_tw else 'CHECK'}")
    print(f"  coupling correlates peel-offs (corr rises with g): "
          f"{corrs[0]:+.2f} → {corrs[-1]:+.2f}  {'PASS' if mono and corrs[-1] > 0.3 else 'CHECK'}")

    # ---- figure ----
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))

    ax[0].plot(g_arr, corrs, "o-", color="#b3402b")
    ax[0].set_ylim(-0.05, 1.02)
    ax[0].set_xlabel("inner coupling g"); ax[0].set_ylabel("corr(Y_node¹, Y_node²)")
    ax[0].set_title("(A) coupling correlates the two peel-offs")
    ax[0].axhline(0, color="grey", lw=0.5, ls=":")

    for gi, col, lab in [(0, "#888780", f"g=0 (corr {corrs[0]:+.2f})"),
                         (len(g_arr) - 1, "#1f3b73", f"g={g_arr[-1]} (corr {corrs[-1]:+.2f})")]:
        Y1 = Yn[gi, :, 0]; Y2 = Yn[gi, :, 1]
        ok = np.isfinite(Y1) & np.isfinite(Y2)
        ax[1].plot(Y1[ok], Y2[ok], ".", ms=2.5, alpha=0.25, color=col, label=lab)
    ax[1].set_xlabel("unit 1 peel-off Y¹"); ax[1].set_ylabel("unit 2 peel-off Y²")
    ax[1].set_title("(B) joint peel-off law: g=0 vs strong g")
    ax[1].legend(fontsize=8.5, frameon=False)

    ax[2].plot(g_arr, means, "o-", color="#1f3b73", label="marginal mean")
    ax[2].plot(g_arr, stds, "s--", color="#2c7d59", label="marginal std")
    ax[2].set_xlabel("inner coupling g"); ax[2].set_ylabel("peel-off level")
    ax[2].set_title("(C) marginal shifts with coupling")
    ax[2].legend(fontsize=9, frameon=False)

    fig.suptitle("T1.1 core — two coupled folded cycles: coupling correlates & shifts "
                 "the peel-off law (β=2)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fp = os.path.join(FIG, "coupled_peeloff.png")
    fig.savefig(fp, dpi=140, bbox_inches="tight")
    print(f"\n  Figure: {fp}")
    print(f"  TOTAL {time.time()-t0:.1f}s — DONE")


if __name__ == "__main__":
    main()
