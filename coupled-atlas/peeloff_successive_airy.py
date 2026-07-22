"""
peeloff_successive_airy.py — T1.2: successive peel-offs → the Airy_β point process.
===================================================================================

The first node of the swept Cole–Hopf solution u is the peel-off level
Y_node^(1) = -Λ0(β) = TW_β.  The SUCCESSIVE nodes
    Y_node^(1) > Y_node^(2) > Y_node^(3) > ...
are the successive edge eigenvalues -Λ_k(β) of the stochastic Airy operator — i.e.
the **Airy_β point process**, the natural multi-point extension of the TW_β marginal.
This is the inner-level realization of the paper's §9.4 conjecture that "finer
KPZ-class objects (the Airy process)" govern the joint law of successive peel-offs.

Parameter-free anchor: as η→0 (β→∞) the k-th node → the k-th zero a_k of Ai
(a_1=-2.338, a_2=-4.088, a_3=-5.521, a_4=-6.787; Sturm/deterministic limit).
At finite η the levels fluctuate AND correlate (level repulsion / rigidity).

Signatures measured here:
  • mean of the k-th node vs the Airy zeros (deterministic-limit check)
  • marginal histograms of the 1st/2nd/3rd peel-off
  • consecutive-level correlation corr(Y^(1), Y^(2))  (>0 ⇒ Airy-process rigidity)

NB: successive NODES of one swept operator realize the Airy point process (eigenvalue
configuration). The full Airy₂ *process* across successive forced passages — the
time-indexed object — is the next layer (needs the forcing/coupling dimension).

Output: figures/peeloff_successive_airy.png + printed summary.
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


def successive_nodes(eta, K=4, Y0=8.0, dt=8e-4, N=12000, Ymin=-13.0, seed=0):
    """First K nodes (zeros) of u in u''=(Y-η ξ)u, Y=Y0-T (Heun/Stratonovich)."""
    rng = np.random.default_rng(seed)
    u = np.ones(N); v = np.full(N, np.sqrt(Y0)); Y = Y0; sdt = np.sqrt(dt)
    Yk = np.full((N, K), np.nan); cnt = np.zeros(N, int); idx = np.arange(N)
    n = int((Y0 - Ymin) / dt)
    for _ in range(n):
        u_prev = u
        dB = sdt * rng.standard_normal(N)
        u1 = u + v * dt; Yp = Y - dt; v1 = v + (Y * u) * dt - eta * u * dB
        u = u + 0.5 * (v + v1) * dt
        v = v + 0.5 * (Y * u + Yp * u1) * dt - eta * (0.5 * (u + u1)) * dB
        Y = Yp
        crossed = ((u_prev > 0) != (u > 0)) & (cnt < K)
        if crossed.any():
            r = idx[crossed]
            Yk[r, cnt[r]] = Y
            cnt[r] += 1
        if (cnt >= K).all():
            break
    return Yk


def main():
    t0 = time.time()
    K = 4
    # first zeros of the Airy function Ai (negative); deterministic η→0 peel-off levels
    a_zeros = np.array([-2.338107, -4.087949, -5.520560, -6.786708])[:K]
    print("=" * 74)
    print("T1.2 — successive peel-offs = Airy_β point process")
    print("=" * 74)
    print(f"  Airy zeros (deterministic η→0 limit): "
          f"{', '.join(f'{a:.3f}' for a in a_zeros)}")

    runs = {}
    for eta in (0.4, 1.0):
        Yk = successive_nodes(eta, K=K, seed=7)
        runs[eta] = Yk
        means = np.nanmean(Yk, axis=0)
        frac = np.mean(np.isfinite(Yk), axis=0)
        b = 4 / eta**2
        print(f"\n  η={eta} (β={b:.2f}) — mean of k-th peel-off (fraction resolved):")
        for k in range(K):
            print(f"     k={k+1}:  {means[k]:7.3f}   (Airy zero {a_zeros[k]:7.3f}, "
                  f"resolved {100*frac[k]:.0f}%)")

    # correlation / rigidity at eta=1.0
    Yk = runs[1.0]
    ok = np.isfinite(Yk[:, 0]) & np.isfinite(Yk[:, 1])
    rho12 = float(np.corrcoef(Yk[ok, 0], Yk[ok, 1])[0, 1])
    print(f"\n  consecutive-level correlation corr(Y^(1),Y^(2)) at η=1.0:  {rho12:+.3f}"
          f"   ({'positive ⇒ rigidity (Airy-process signature)' if rho12 > 0.1 else 'check'})")

    small = np.nanmean(runs[0.4], axis=0)
    conv = np.max(np.abs(small - a_zeros))
    print(f"  η→0 check: max|mean_node - Airy_zero| at η=0.4 = {conv:.3f}  "
          f"({'PASS' if conv < 0.25 else 'CHECK'})")

    # ---- figure ----
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
    cols = ["#1f3b73", "#d95f0e", "#2c7d59", "#7a3b8f"]

    Yk1 = runs[1.0]
    for k in range(3):
        x = Yk1[np.isfinite(Yk1[:, k]), k]
        ax[0].hist(x, bins=70, density=True, histtype="step", lw=1.7, color=cols[k],
                   label=f"peel-off #{k+1}")
        ax[0].axvline(a_zeros[k], color=cols[k], lw=0.8, ls=":")
    ax[0].set_xlabel("peel-off level Y_node"); ax[0].set_ylabel("density")
    ax[0].set_title("(A) successive peel-offs, η=1.0 (β=4)\n dotted = Airy zeros (η→0 limit)")
    ax[0].legend(fontsize=9, frameon=False)

    ax[1].plot(Yk1[ok, 0], Yk1[ok, 1], ".", ms=2.5, alpha=0.25, color="#1f3b73")
    ax[1].set_xlabel("1st peel-off Y^(1)"); ax[1].set_ylabel("2nd peel-off Y^(2)")
    ax[1].set_title(f"(B) consecutive-level coupling\n corr = {rho12:+.2f} (rigidity)")

    ks = np.arange(1, K + 1)
    for eta, col, mk in [(0.4, "#888780", "s"), (1.0, "#1f3b73", "o")]:
        ax[2].plot(ks, np.nanmean(runs[eta], axis=0), mk + "-", color=col,
                   label=f"mean node, η={eta} (β={4/eta**2:.1f})")
    ax[2].plot(ks, a_zeros, "*", color="#b3402b", ms=14, label="Airy zeros a_k")
    ax[2].set_xlabel("peel-off index k"); ax[2].set_ylabel("level")
    ax[2].set_title("(C) node means → Airy zeros as η→0")
    ax[2].set_xticks(ks); ax[2].legend(fontsize=9, frameon=False)

    fig.suptitle("T1.2 — successive peel-offs realize the Airy_β point process "
                 "(multi-point extension of TW_β)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fp = os.path.join(FIG, "peeloff_successive_airy.png")
    fig.savefig(fp, dpi=140, bbox_inches="tight")
    print(f"\n  Figure: {fp}")
    print(f"  TOTAL {time.time()-t0:.1f}s — DONE")


if __name__ == "__main__":
    main()
