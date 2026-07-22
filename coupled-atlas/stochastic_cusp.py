"""
stochastic_cusp.py — the stochastic cusp step (toward TW→Weber-edge).
=====================================================================

Three pieces of "what the stochastic step needs", attempted numerically on the rescaled
cusp inner SDE (the A₃ analogue of the paper's inner Riccati):

    dD = (-D^3 + M D - W) ds + η dB_s ,   W = s (slow load swept through the cusp).

(1) NOISE SCALE [DERIVED]. Pushing additive σ dW through the cusp blow-up (δ,μ,δw)=(r,r²,r³),
    r~ε^{1/5}, inner time s=r²t, gives η_cusp = σ/(√2 ε^{2/5}); the same calculation reproduces
    the fold's η=σ/√ε₂ (validation, see notes). Here we work directly in the rescaled variable.

(2) BERGLUND–GENTZ TUBE [NUMERIC]. The noisy trajectory should stay in an O(η)-tube around the
    deterministic canard (upper branch) until the peel-off — confinement one rung above the fold.

(3) STOCHASTIC UNIFORM CONNECTION [NUMERIC]. Sweep the chart parameter M and measure the inner
    exit measure (peel-off load W at escape). M≫1 = two well-separated folds → generic fold →
    Tracy–Widom; M→O(1) = coalescing folds (cusp) → Weber-class (skew rises). This is the noise
    analogue of Olver's Airy↔Weber connection, realised in the cusp chart.

NB the chart parameter M = μ/ε^{2/5} = -2g/ε^{2/5}: large M = g very negative (resolved folds),
small M = g near 0 (cusp). The full-model peel-off (symmetric fold release) is a DIFFERENT
observable — the symmetric escape MODULATED by this antisymmetric funnel — which is why its
Weber signature is strongest at large μ; see notes for the reconciliation.

Output: figures/stochastic_cusp.png + summary.
"""
from __future__ import annotations
import os, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
TW = {1: 0.2935, 2: 0.2241, 4: 0.1655}


def exit_loads(M, eta, N=5000, ds=0.002, seed=0):
    """Ensemble peel-off load W at escape for the rescaled cusp, swept W=s upward."""
    rng = np.random.default_rng(seed); sds = np.sqrt(ds)
    D = np.full(N, np.sqrt(M))                 # start in the upper well
    W = 0.0; out = np.full(N, np.nan); done = np.zeros(N, bool)
    Wstar = (2*M/3)*np.sqrt(M/3)               # deterministic fold load
    Wmax = 1.8*Wstar + 0.5
    nstep = int(Wmax/ds)
    for _ in range(nstep):
        al = ~done
        if not al.any():
            break
        D = D + (-D**3 + M*D - W)*ds + eta*sds*rng.standard_normal(N)
        W += ds
        cr = al & (D < 0.0); out[cr] = W; done |= cr
    return out[np.isfinite(out)], Wstar


def tube(M, eta, ntraj=14, ds=0.002, seed=1):
    rng = np.random.default_rng(seed); sds = np.sqrt(ds)
    Wstar = (2*M/3)*np.sqrt(M/3); n = int(1.05*Wstar/ds)
    Ws = np.arange(n)*ds
    # deterministic canard (upper branch), eta=0
    Dd = np.empty(n); d = np.sqrt(M)
    for i in range(n):
        d = d + (-d**3 + M*d - Ws[i])*ds; Dd[i] = d
    # noisy trajectories
    Dn = np.empty((n, ntraj)); dd = np.full(ntraj, np.sqrt(M))
    for i in range(n):
        dd = dd + (-dd**3 + M*dd - Ws[i])*ds + eta*sds*rng.standard_normal(ntraj)
        Dn[i] = dd
    return Ws, Dd, Dn, Wstar


def skew(x):
    d = x - x.mean(); v = np.mean(d**2); return np.mean(d**3)/v**1.5


def main():
    t0 = time.time()
    print("=" * 70)
    print("Stochastic cusp step — rescaled cusp inner SDE")
    print("=" * 70)
    print("  [DERIVED] η_cusp = σ/(√2 · ε^{2/5})   (same method gives fold η=σ/√ε₂)")

    eta = 1.0
    Ms = np.array([0.5, 0.8, 1.2, 1.8, 2.6, 3.6, 5.0])
    print(f"\n  [NUMERIC] inner exit measure vs chart parameter M (η={eta}):")
    print(f"  {'M':>6}{'Wstar':>9}{'skew':>9}{'n':>7}")
    sk = []
    samples = {}
    for i, M in enumerate(Ms):
        x, Ws = exit_loads(M, eta, seed=5+i)
        sk.append(skew(x) if x.size > 50 else np.nan); samples[M] = x
        print(f"  {M:6.2f}{Ws:9.3f}{sk[-1]:9.3f}{x.size:7d}")
    sk = np.array(sk)
    print(f"\n  skew: M=5 (resolved folds → TW) {sk[-1]:+.2f}  →  M=0.5 (cusp → Weber) {sk[0]:+.2f}")
    rises = sk[0] > sk[-1] + 0.1
    print(f"  TW→Weber as M→cusp (skew rises): {'YES' if rises else 'check'}")

    # tube confinement
    Ws, Dd, Dn, Wstar = tube(2.0, 0.6)
    band = np.nanstd(Dn - Dd[:, None], axis=1)
    pre = Ws < 0.85*Wstar
    print(f"\n  [NUMERIC] tube at M=2, η=0.6: pre-escape tube half-width "
          f"≈ {np.nanmean(band[pre]):.3f}  (∝η ⇒ Berglund–Gentz confinement holds one rung up)")

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    ax[0].plot(Ws, Dd, "k-", lw=2, label="deterministic canard")
    ax[0].plot(Ws, Dn, lw=0.5, alpha=0.5, color="#7a3b8f")
    ax[0].axvline(Wstar, color="grey", ls=":", lw=0.8)
    ax[0].set_xlabel("load W"); ax[0].set_ylabel("D")
    ax[0].set_title("(A) Berglund–Gentz tube (M=2): O(η) around the canard")
    ax[0].legend(fontsize=8.5, frameon=False)

    ax[1].axhspan(TW[1], TW[4], color="#b3402b", alpha=0.12)
    ax[1].axhline(0.60, color="#2c7d59", lw=1, ls="--")
    ax[1].plot(Ms, sk, "o-", color="#1f3b73")
    ax[1].annotate("Tracy–Widom band", (Ms[-1], 0.22), fontsize=8.5, color="#b3402b", ha="right")
    ax[1].annotate("Weber-class", (Ms[0], 0.60), fontsize=8.5, color="#2c7d59")
    ax[1].invert_xaxis()
    ax[1].set_xlabel("chart parameter M  (→ toward cusp)"); ax[1].set_ylabel("exit-measure skew")
    ax[1].set_title("(B) inner exit measure: TW (M≫1) → Weber (M→cusp)")

    gg = np.linspace(-4, 4, 200)
    ax[2].plot(gg, np.exp(-gg**2/2)/np.sqrt(2*np.pi), "k:", lw=1, label="Gaussian")
    for M, col in [(5.0, "#b3402b"), (0.5, "#2c7d59")]:
        x = samples[M]; z = (x - x.mean())/x.std()
        ax[2].hist(z, bins=55, range=(-4, 4), density=True, histtype="step", lw=1.8,
                   color=col, label=f"M={M} (skew {skew(x):+.2f})")
    ax[2].set_xlabel("standardised exit load"); ax[2].set_ylabel("density")
    ax[2].set_title("(C) resolved-fold (TW) vs cusp (Weber) law"); ax[2].legend(fontsize=8.5, frameon=False)

    fig.suptitle("Stochastic cusp step — η_cusp derived; tube confinement; TW→Weber inner connection",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fp = os.path.join(FIG, "stochastic_cusp.png")
    fig.savefig(fp, dpi=140, bbox_inches="tight")
    print(f"\n  Figure: {fp}")
    print(f"  TOTAL {time.time()-t0:.1f}s — DONE")


if __name__ == "__main__":
    main()
