"""
forced_cusp_process.py — the cusp's successive-peel-off process (Weber analogue of Airy₂).
=========================================================================================

Completes the spine for the cusp rung: the fold rung has both a marginal (TW) and a
multi-point PROCESS (forced → Airy₂-type). Here we build the cusp counterpart.

Successive forced passages of the inner equation u'' = (V(Y) - η ξ_k) u, first node = peel-off,
with V(Y)=sign(Y)|Y|^k (k=1 fold/Airy; k=2 cusp/Weber) and OU-correlated noise across passages
    ξ_k(Y) = ρ ξ_{k-1}(Y) + √(1-ρ²) ζ_k(Y)     (ρ = forcing-set inter-passage correlation).

For each k and ρ we measure: the marginal (stationary across passages), the autocovariance C(s)
(decaying ⇒ a genuine process), and the increment variance (→ 2·Var ⇒ decorrelation). We compare
the cusp process (k=2, Weber marginal) to the fold process (k=1, TW/Airy₂) — does the cusp give a
stationary WEBER-marginal process, and how does its covariance/rigidity differ?

Output: figures/forced_cusp_process.png + tagged summary.
"""
from __future__ import annotations
import os, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
TW2_SKEW = 0.2241


def forced_process(kk, eta, rho, K=32, M=1000, Y0=5.0, dt=2.5e-3, Ymin=-6.0, seed=0):
    rng = np.random.default_rng(seed); sdt = np.sqrt(dt); n = int((Y0 - Ymin) / dt)
    Ynode = np.full((M, K), np.nan)
    path_prev = rng.standard_normal((M, n))
    for k in range(K):
        fresh = rng.standard_normal((M, n))
        path = rho * path_prev + np.sqrt(1 - rho**2) * fresh
        path_prev = path
        u = np.ones(M); v = np.full(M, Y0**(kk / 2.0)); Y = Y0
        Yz = np.full(M, np.nan); done = np.zeros(M, bool)
        for step in range(n):
            al = ~done
            if not al.any():
                break
            dB = path[:, step] * sdt
            Vc = np.sign(Y) * abs(Y)**kk; Yp = Y - dt; Vp = np.sign(Yp) * abs(Yp)**kk
            u1 = u + v * dt; v1 = v + (Vc * u) * dt - eta * u * dB
            u = u + 0.5 * (v + v1) * dt
            v = v + 0.5 * (Vc * u + Vp * u1) * dt - eta * 0.5 * (u + u1) * dB
            Y = Yp
            cr = al & (u < 0.0); Yz[cr] = Y; done |= cr
        Ynode[:, k] = Yz
    return Ynode


def acov(Y, burn=6, Smax=12):
    Z = Y[:, burn:]; mu = np.nanmean(Z); var = np.nanvar(Z); Kb = Z.shape[1]
    C = np.full(Smax + 1, np.nan); V = np.full(Smax + 1, np.nan)
    for s in range(Smax + 1):
        a = Z[:, :Kb - s] - mu; b = Z[:, s:] - mu
        ok = np.isfinite(a) & np.isfinite(b)
        C[s] = np.mean(a[ok] * b[ok]) / var
        V[s] = np.nanvar(Z[:, s:] - Z[:, :Kb - s])
    return C, V, var


def skew(x):
    x = x[np.isfinite(x)]; d = x - x.mean(); v = np.mean(d**2); return np.mean(d**3) / v**1.5


def main():
    t0 = time.time()
    eta = np.sqrt(2.0)
    print("=" * 74)
    print("Cusp successive-peel-off process (Weber analogue of Airy₂),  β=2")
    print("=" * 74)
    runs = {}
    for kk in (1, 2):
        for rho in (0.0, 0.85):
            Y = forced_process(kk, eta, rho, seed=10 + kk + int(10*rho))
            C, V, var = acov(Y)
            runs[(kk, rho)] = dict(Y=Y, C=C, V=V, var=var)
        Yp = runs[(kk, 0.85)]["Y"]
        sk = skew(Yp[:, 6:].ravel())
        pm = np.nanmean(Yp[:, 6:], axis=0)
        name = "fold (k=1, Airy/TW)" if kk == 1 else "cusp (k=2, Weber)"
        print(f"\n  {name}:")
        print(f"    marginal skew {sk:+.3f} ({'≈TW₂ '+str(TW2_SKEW) if kk==1 else 'Weber-class'}); "
              f"stationary mean {pm.min():.2f}…{pm.max():.2f}")
        print(f"    autocorr (ρ=0.85): C(1)={runs[(kk,0.85)]['C'][1]:+.3f} "
              f"C(4)={runs[(kk,0.85)]['C'][4]:+.3f} C(8)={runs[(kk,0.85)]['C'][8]:+.3f}")
        print(f"    iid check (ρ=0):   C(1)={runs[(kk,0.0)]['C'][1]:+.3f}")
        print(f"    increment V(∞)/2Var = {runs[(kk,0.85)]['V'][-1]/(2*runs[(kk,0.85)]['var']):.2f}")

    # decorrelation length comparison
    def corr_len(C):
        s = np.arange(len(C))
        below = np.where(C < 0.2)[0]
        return below[0] if below.size else len(C)
    lf = corr_len(runs[(1, 0.85)]["C"]); lc = corr_len(runs[(2, 0.85)]["C"])
    print(f"\n  decorrelation lag (C<0.2): fold {lf}  vs  cusp {lc}  "
          f"⇒ cusp process is {'distinct' if lc != lf else 'similar in covariance'} "
          f"(marginal differs: TW vs Weber)")

    # ---- figure ----
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    for kk, col, lab in [(1, "#b3402b", "fold (TW/Airy₂)"), (2, "#7a3b8f", "cusp (Weber)")]:
        x = runs[(kk, 0.85)]["Y"][:, 6:].ravel(); z = (x[np.isfinite(x)]-np.nanmean(x))/np.nanstd(x)
        ax[0].hist(z, bins=60, range=(-4, 4), density=True, histtype="step", lw=1.8,
                   color=col, label=f"{lab}, skew {skew(x):+.2f}")
    gg = np.linspace(-4, 4, 200)
    ax[0].plot(gg, np.exp(-gg**2/2)/np.sqrt(2*np.pi), "k:", lw=1, label="Gaussian")
    ax[0].set_xlabel("standardised peel-off"); ax[0].set_ylabel("density")
    ax[0].set_title("(A) marginal: fold TW vs cusp Weber"); ax[0].legend(fontsize=8.5, frameon=False)

    s = np.arange(len(runs[(1, 0.85)]["C"]))
    ax[1].plot(s, runs[(1, 0.85)]["C"], "o-", color="#b3402b", label="fold C(s) (Airy₂)")
    ax[1].plot(s, runs[(2, 0.85)]["C"], "s-", color="#7a3b8f", label="cusp C(s) (Weber)")
    ax[1].plot(s, runs[(1, 0.0)]["C"], ":", color="#b3402b", alpha=0.6, label="iid (ρ=0)")
    ax[1].plot(s, runs[(2, 0.0)]["C"], ":", color="#7a3b8f", alpha=0.6)
    ax[1].axhline(0, color="grey", lw=0.5, ls=":")
    ax[1].set_xlabel("passage lag s"); ax[1].set_ylabel("autocorrelation C(s)")
    ax[1].set_title("(B) both are processes; covariance compared"); ax[1].legend(fontsize=8.5, frameon=False)

    for kk, col in [(1, "#b3402b"), (2, "#7a3b8f")]:
        V = runs[(kk, 0.85)]["V"]; ax[2].plot(s, V/(2*runs[(kk,0.85)]["var"]), "o-", color=col,
                                              label=("fold" if kk == 1 else "cusp")+" V(s)/2Var")
    ax[2].axhline(1.0, color="grey", lw=1, ls="--", label="2·Var (decorrelated)")
    ax[2].set_xlabel("passage lag s"); ax[2].set_ylabel("increment var / 2Var")
    ax[2].set_title("(C) increments saturate at 2·Var"); ax[2].legend(fontsize=8.5, frameon=False)

    fig.suptitle("Cusp successive-peel-off process — a stationary Weber-marginal process "
                 "(the Weber analogue of Airy₂)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fp = os.path.join(FIG, "forced_cusp_process.png")
    fig.savefig(fp, dpi=140, bbox_inches="tight")
    print(f"\n  Figure: {fp}")
    print(f"  TOTAL {time.time()-t0:.1f}s — DONE")


if __name__ == "__main__":
    main()
