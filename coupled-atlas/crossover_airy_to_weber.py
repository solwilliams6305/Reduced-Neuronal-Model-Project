"""
crossover_airy_to_weber.py — PART 2: the analytic Airy → Weber crossover (two folds merging).
=============================================================================================

Analytic shadow of the numerical g-crossover. The symmetric-mode critical manifold near the
cusp has TWO folds at ±Δ(g), with Δ → 0 at the cusp. The inner equation is
    u'' = (V(Y) - ηξ) u ,    V(Y) ≈ Y^2 - Δ^2   (two simple turning points at ±Δ).
Classical uniform asymptotics (Olver, Asymptotics & Special Functions, ch. on coalescing
turning points):
    • an ISOLATED turning point  (Δ ≫ ℓ, ℓ = inner blow-up scale) → Airy → (stochastic) → TW;
    • TWO COALESCING turning points (Δ ≲ ℓ)                       → parabolic cylinder / Weber
                                                                   → (stochastic) → cusp edge law.
So the effective turning order k_eff runs 1 → 2 as Δ/ℓ : ∞ → 0. This panelises:
  (A) the potential morph V = Y^2 - Δ^2 (two turnings coalescing),
  (B) the inner crossover — peel-off skew vs turning order k (1 = Airy/TW → 2 = Weber),
  (C) the physical g-crossover (from crossover_fold_to_cusp.py): skew swing + ~4× spread
      at the cusp onset g ≈ -0.05 … -0.11, i.e. where Δ(g) ~ ℓ.

Output: figures/crossover_airy_to_weber.png + printed summary.  Notes: NOISY_CUSP_CROSSOVER_NOTES.md
"""
from __future__ import annotations
import os, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from peeloff_cusp_ladder import peeloff_ladder

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")

# physical g-crossover results (from crossover_fold_to_cusp.py run)
G = np.array([0.10, 0.05, 0.0, -0.04, -0.08, -0.11, -0.14, -0.17, -0.20])
G_SKEW = np.array([+0.221, +0.015, -0.056, +0.103, -0.191, -0.434, -0.629, -0.690, -0.703])
G_STD = np.array([0.0049, 0.0040, 0.0042, 0.0060, 0.0089, 0.0111, 0.0136, 0.0170, 0.0205])


def skew(x):
    x = x[np.isfinite(x)]; d = x - x.mean(); v = np.mean(d**2)
    return np.mean(d**3) / v**1.5


def main():
    t0 = time.time()
    eta = np.sqrt(2.0)
    print("=" * 72)
    print("PART 2 — analytic Airy → Weber crossover (coalescing turning points)")
    print("=" * 72)

    # (B) inner turning-order crossover k: 1 → 2
    ks = np.linspace(1.0, 2.0, 9)
    sk = np.array([skew(peeloff_ladder([k], eta, N=6000, seed=30 + i)[0])
                   for i, k in enumerate(ks)])
    print("  inner crossover (turning order k → peel-off skew):")
    for k, s in zip(ks, sk):
        print(f"    k={k:.2f}: skew {s:+.3f}")
    print(f"  Airy (k=1) skew {sk[0]:+.3f} → Weber (k=2) skew {sk[-1]:+.3f}")
    print(f"\n  physical g-crossover (orientation: full-model w-release):")
    print(f"    fold side g≥0: skew≈{G_SKEW[G>=0].mean():+.2f}, std≈{G_STD[G>=0].mean():.4f}")
    print(f"    cusp side g≤-0.11: skew≈{G_SKEW[G<=-0.11].mean():+.2f}, "
          f"std≈{G_STD[G<=-0.11].mean():.4f}  (spread ×{G_STD[G<=-0.11].mean()/G_STD[G>=0].mean():.1f})")
    print("  ⇒ same TW→Weber swing; crossover where Δ(g) ~ inner scale ℓ (g≈-0.05…-0.11).")

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))

    # (A) potential morph
    Y = np.linspace(-2.2, 2.2, 400)
    for D, col, lab in [(1.5, "#b3402b", "Δ=1.5 (two folds → Airy)"),
                        (0.9, "#d07b30", "Δ=0.9"),
                        (0.45, "#5a7d2c", "Δ=0.45"),
                        (0.0, "#2c7d59", "Δ=0 (merged → Weber)")]:
        ax[0].plot(Y, Y**2 - D**2, color=col, lw=1.8, label=lab)
        if D > 0:
            ax[0].plot([-D, D], [0, 0], "o", color=col, ms=4)
    ax[0].axhline(0, color="grey", lw=0.6, ls=":")
    ax[0].set_xlabel("inner slow variable Y"); ax[0].set_ylabel("inner potential V(Y)=Y²−Δ²")
    ax[0].set_title("(A) two turning points coalesce (Airy→Weber)")
    ax[0].legend(fontsize=8, frameon=False)

    # (B) inner turning-order crossover
    ax[1].plot(ks, sk, "o-", color="#1f3b73")
    ax[1].axvline(1, color="#b3402b", lw=0.7, ls=":"); ax[1].axvline(2, color="#2c7d59", lw=0.7, ls=":")
    ax[1].annotate("Airy / TW", (1.0, sk[0]), textcoords="offset points", xytext=(6, -12), fontsize=9, color="#b3402b")
    ax[1].annotate("Weber", (2.0, sk[-1]), textcoords="offset points", xytext=(-38, 6), fontsize=9, color="#2c7d59")
    ax[1].set_xlabel("effective turning order  k  (= 1+Δ-merging)")
    ax[1].set_ylabel("inner peel-off skew")
    ax[1].set_title("(B) inner crossover: turning order 1→2")

    # (C) physical g-crossover
    ax2 = ax[2].twinx()
    ax[2].plot(G, G_SKEW, "o-", color="#7a3b8f", label="peel-off skew")
    ax2.plot(G, G_STD, "s--", color="#2c7d59", label="peel-off std", alpha=0.8)
    ax[2].axvspan(-0.11, -0.05, color="#7a3b8f", alpha=0.10)
    ax[2].invert_xaxis()
    ax[2].set_xlabel("coupling g  (→ toward cusp)")
    ax[2].set_ylabel("skew", color="#7a3b8f"); ax2.set_ylabel("std", color="#2c7d59")
    ax[2].set_title("(C) physical g-crossover (KP FHN): cusp onset shaded")
    ax[2].legend(fontsize=8, frameon=False, loc="lower left")
    ax2.legend(fontsize=8, frameon=False, loc="upper left")

    fig.suptitle("Part 2 — Airy → Weber as two folds merge: inner crossover (turning order) ↔ "
                 "physical g-crossover", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fp = os.path.join(FIG, "crossover_airy_to_weber.png")
    fig.savefig(fp, dpi=140, bbox_inches="tight")
    print(f"\n  Figure: {fp}")
    print(f"  TOTAL {time.time()-t0:.1f}s — DONE")


if __name__ == "__main__":
    main()
