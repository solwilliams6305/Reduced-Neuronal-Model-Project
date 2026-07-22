"""
crossover_fold_to_cusp.py — the unifying experiment: coupling drives a fold→cusp
crossover of the peel-off edge class (fuses T1.1 and the cusp frontier).
================================================================================

Genuine Kristiansen–Pedersen coupled FHN (gap-junction on the fast variable):
    v_i' = -v_i^3 + 3 v_i - w_i + g(v_j - v_i),   w_i' = ε(v_i - c).

Reframing (per the singularity→edge-class spine): the peel-off law is NOT universally
Tracy–Widom — it depends on the singularity the system sits at. Weak/attractive coupling
keeps two generic FOLDS (fold-class escape). Repulsive coupling MERGES the folds into a
CUSP (Weber/folded-node-class escape; Kristiansen–Pedersen). So as g sweeps from the
two-fold regime into the symmetric cusp, the peel-off edge class should CROSS OVER.

We sweep g and measure the peel-off law (slow variable w_1 at the upper-fold escape):
its shape (skew), spread (std), and mean. Fold-class → cusp-class shows up as a change
in shape and a jump in spread. (Full-model observable carries outer corrections — the
crossover is the robust signature, not the absolute skew sign.)

Output: figures/crossover_fold_to_cusp.png + printed summary.
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


def sweep(g_vals, c=0.99, eps=0.015, sigma=0.02, M=350, T=2400.0, dt=0.02,
          warmup=700.0, seed=0):
    ng = len(g_vals)
    G = np.asarray(g_vals, float)[:, None]                 # (ng,1)
    rng = np.random.default_rng(seed); sdt = np.sqrt(dt)
    v = np.tile([0.2, -0.2], (ng, M, 1)).astype(float); w = np.zeros((ng, M, 2))

    def step(v, w, noise):
        coup = G[:, :, None] * (v[:, :, ::-1] - v)
        vn = v + (-v**3 + 3*v - w + coup) * dt + noise
        wn = w + eps * (v - c) * dt
        return vn, wn

    for _ in range(int(warmup/dt)):
        v, w = step(v, w, 0.0)

    on_up = np.zeros((ng, M), bool); wmax = np.full((ng, M), -9.0)
    gidx = np.tile(np.arange(ng)[:, None], (1, M))
    pw, pg = [], []
    n = int(T/dt)
    for _ in range(n):
        noise = sigma * sdt * rng.standard_normal((ng, M, 2))
        v, w = step(v, w, noise)
        v1, w1 = v[:, :, 0], w[:, :, 0]
        on_up |= (v1 > 1.2)
        wmax = np.where(on_up, np.maximum(wmax, w1), wmax)
        esc = on_up & (v1 < 0.5)
        if esc.any():
            pw.append(wmax[esc].copy()); pg.append(gidx[esc].copy())
        on_up = np.where(esc, False, on_up)
        wmax = np.where(esc, -9.0, wmax)
    return (np.concatenate(pw) if pw else np.array([]),
            np.concatenate(pg) if pg else np.array([], int))


def moments(x):
    m = x.mean(); d = x - m; var = np.mean(d**2)
    return m, np.sqrt(var), np.mean(d**3)/var**1.5, np.mean(d**4)/var**2 - 3.0


def main():
    t0 = time.time()
    g_vals = np.array([0.10, 0.05, 0.0, -0.04, -0.08, -0.11, -0.14, -0.17, -0.20])
    pw, pg = sweep(g_vals)
    print("=" * 70)
    print("Fold → cusp crossover of the peel-off edge class (KP coupled FHN)")
    print("=" * 70)
    print(f"  {'g':>7} | {'mean':>8}{'std':>9}{'skew':>8}{'kurt':>8} {'n':>6}")
    means, stds, skews = [], [], []
    samples = {}
    for gi, g in enumerate(g_vals):
        x = pw[pg == gi]
        samples[gi] = x
        if x.size > 30:
            m, s, sk, ku = moments(x)
            means.append(m); stds.append(s); skews.append(sk)
            print(f"  {g:7.3f} | {m:8.3f}{s:9.4f}{sk:+8.3f}{ku:+8.3f} {x.size:6d}")
        else:
            means.append(np.nan); stds.append(np.nan); skews.append(np.nan)
            print(f"  {g:7.3f} | (too few escapes: {x.size})")
    means, stds, skews = map(np.array, (means, stds, skews))

    # crossover diagnostics
    foldmask = g_vals >= 0.0; cuspmask = g_vals <= -0.11
    sk_fold = np.nanmean(skews[foldmask]); sk_cusp = np.nanmean(skews[cuspmask])
    sd_fold = np.nanmean(stds[foldmask]); sd_cusp = np.nanmean(stds[cuspmask])
    print(f"\n  fold side (g≥0):    skew {sk_fold:+.3f}, std {sd_fold:.4f}")
    print(f"  cusp side (g≤-0.11): skew {sk_cusp:+.3f}, std {sd_cusp:.4f}")
    print(f"  CROSSOVER: skew shifts {sk_fold:+.2f} → {sk_cusp:+.2f}; "
          f"spread ×{sd_cusp/sd_fold:.1f}  "
          f"({'clear' if (sd_cusp > 2*sd_fold and sk_cusp < sk_fold - 0.15) else 'check'})")

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    gc = -0.105   # approx cusp onset for shading
    for a_ in ax[:2]:
        a_.axvspan(g_vals.min()-0.01, gc, color="#7a3b8f", alpha=0.07)
        a_.axvspan(gc, g_vals.max()+0.01, color="#b3402b", alpha=0.06)
    ax[0].plot(g_vals, skews, "o-", color="#1f3b73")
    ax[0].set_xlabel("coupling g  (←repulsive/cusp · attractive/fold→)")
    ax[0].set_ylabel("peel-off skewness"); ax[0].set_title("(A) shape crosses over")
    ax[0].invert_xaxis()
    ax[1].plot(g_vals, stds, "s-", color="#2c7d59")
    ax[1].set_xlabel("coupling g"); ax[1].set_ylabel("peel-off spread (std)")
    ax[1].set_title("(B) spread amplifies at the cusp"); ax[1].invert_xaxis()

    for gi, col, lab in [(1, "#b3402b", f"g={g_vals[1]:+.2f} (fold side)"),
                         (7, "#7a3b8f", f"g={g_vals[7]:+.2f} (cusp side)")]:
        x = samples[gi]
        if x.size > 30:
            z = (x - x.mean())/x.std()
            ax[2].hist(z, bins=45, range=(-4, 4), density=True, histtype="step", lw=1.8,
                       color=col, label=lab)
    gg = np.linspace(-4, 4, 200)
    ax[2].plot(gg, np.exp(-gg**2/2)/np.sqrt(2*np.pi), "k:", lw=1, label="Gaussian")
    ax[2].set_xlabel("standardized peel-off"); ax[2].set_ylabel("density")
    ax[2].set_title("(C) fold-side vs cusp-side law"); ax[2].legend(fontsize=8.5, frameon=False)

    fig.suptitle("Coupling-driven fold → cusp crossover of the peel-off edge class "
                 "(unifies T1.1 & the cusp frontier)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fp = os.path.join(FIG, "crossover_fold_to_cusp.png")
    fig.savefig(fp, dpi=140, bbox_inches="tight")
    print(f"\n  Figure: {fp}")
    print(f"  TOTAL {time.time()-t0:.1f}s — DONE")


if __name__ == "__main__":
    main()
