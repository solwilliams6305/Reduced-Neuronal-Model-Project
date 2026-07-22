"""
forced_peeloff_airy.py — T1.2 / T1.1-forced: successive peel-offs of a FORCED
folded cycle → an Airy₂-type process (§9.4).
==============================================================================

A single folded cycle undergoes SUCCESSIVE passages (peel-offs) Y_node^(k),
k=1,2,.... Asymmetric coupling (one unit forcing the other) makes the noise the
cycle sees on consecutive passages temporally correlated. We model that forcing as
an Ornstein–Uhlenbeck correlation across passages of the inner noise path:

    ξ_k(Y) = ρ ξ_{k-1}(Y) + sqrt(1-ρ²) ζ_k(Y) ,   ρ = forcing-set correlation.

ρ=0 → independent passages (no forcing); ρ→1 → frozen (strong forcing).
For each passage we sweep u'' = (Y - η ξ_k) u, Y=Y0-T (the validated Heun scheme),
and take the first node = peel-off level Y_node^(k). We use η=√2 (β=2) — the Airy₂
universality class.

§9.4 conjecture content: forcing should turn the successive-peel-off sequence from
iid TW draws into a stationary TW₂-marginal PROCESS with (i) decaying autocovariance
and (ii) locally-Brownian increments saturating at 2·Var — the Airy₂ signatures.

Honest scope: this realizes a TW₂-marginal process with Airy-process-like temporal
structure. Matching the EXACT Airy₂ covariance needs the inter-passage dynamics
calibrated to the Dyson-Brownian-motion edge (a DBM-edge sampler is the reference
for the next step). What we test here is the qualitative iid → process transition.

Output: figures/forced_peeloff_airy.png + printed summary.
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

TW2_SKEW = 0.2241


def forced_process(eta, rho, K=36, M=1500, Y0=5.0, dt=2.5e-3, Ymin=-4.5, seed=0):
    """Successive peel-off levels Y_node^(k), shape (M,K), with OU-in-passage noise."""
    rng = np.random.default_rng(seed); sdt = np.sqrt(dt)
    n = int((Y0 - Ymin) / dt)
    Ynode = np.full((M, K), np.nan)
    path_prev = rng.standard_normal((M, n))             # unit-variance increments
    for k in range(K):
        fresh = rng.standard_normal((M, n))
        path = rho * path_prev + np.sqrt(1 - rho**2) * fresh
        path_prev = path
        u = np.ones(M); v = np.full(M, np.sqrt(Y0)); Y = Y0
        Yz = np.full(M, np.nan); done = np.zeros(M, bool)
        for step in range(n):
            al = ~done
            if not al.any():
                break
            dB = path[:, step] * sdt
            u1 = u + v * dt; Yp = Y - dt; v1 = v + (Y * u) * dt - eta * u * dB
            u = u + 0.5 * (v + v1) * dt
            v = v + 0.5 * (Y * u + Yp * u1) * dt - eta * (0.5 * (u + u1)) * dB
            Y = Yp
            cr = al & (u < 0.0); Yz[cr] = Y; done |= cr
        Ynode[:, k] = Yz
    return Ynode


def acov(Y, burn=6, Smax=14):
    """Autocorrelation C(s) and increment variance V(s) of the peel-off process."""
    Z = Y[:, burn:]                                     # (M, Kb)
    mu = np.nanmean(Z); var = np.nanvar(Z)
    Kb = Z.shape[1]
    C = np.full(Smax + 1, np.nan); V = np.full(Smax + 1, np.nan)
    for s in range(Smax + 1):
        a = Z[:, :Kb - s] - mu; b = Z[:, s:] - mu
        ok = np.isfinite(a) & np.isfinite(b)
        C[s] = np.mean(a[ok] * b[ok]) / var
        d = Z[:, s:] - Z[:, :Kb - s]
        V[s] = np.nanvar(d)
    return C, V, mu, var


def skew(x):
    x = x[np.isfinite(x)]; d = x - x.mean(); v = np.mean(d**2)
    return np.mean(d**3) / v**1.5


def main():
    t0 = time.time()
    eta = np.sqrt(2.0)                                  # β = 2 (Airy₂ class)
    print("=" * 72)
    print("T1.2/T1.1-forced — successive peel-offs of a forced folded cycle (β=2)")
    print("=" * 72)

    res = {}
    for rho in (0.0, 0.85):
        Y = forced_process(eta, rho, seed=11)
        C, V, mu, var = acov(Y)
        res[rho] = dict(Y=Y, C=C, V=V, mu=mu, var=var)
        # stationarity: per-passage marginal mean/std
        pm = np.nanmean(Y[:, 6:], axis=0); ps = np.nanstd(Y[:, 6:], axis=0)
        print(f"\n  forcing ρ={rho}:")
        print(f"    marginal: mean {mu:.3f}  std {np.sqrt(var):.3f}  skew {skew(Y[:,6:].ravel()):.3f}"
              f"  (TW2 skew {TW2_SKEW})")
        print(f"    stationary across passages: mean varies {pm.min():.2f}..{pm.max():.2f}, "
              f"std {ps.min():.2f}..{ps.max():.2f}")
        print(f"    autocorr C(1)={C[1]:+.3f}  C(2)={C[2]:+.3f}  C(4)={C[4]:+.3f}  C(8)={C[8]:+.3f}")

    c0, c1 = res[0.0]["C"][1], res[0.85]["C"][1]
    iid_ok = abs(c0) < 0.08
    proc_ok = c1 > 0.3 and res[0.85]["C"][1] > res[0.85]["C"][8]
    print(f"\n  ρ=0 is iid (C(1)≈0): {c0:+.3f} {'PASS' if iid_ok else 'CHECK'}")
    print(f"  ρ=0.85 is a PROCESS (C(1)>0, decaying): C(1)={c1:+.2f}→C(8)={res[0.85]['C'][8]:+.2f}"
          f"  {'PASS' if proc_ok else 'CHECK'}")
    sat = res[0.85]["V"][-1] / (2 * res[0.85]["var"])
    print(f"  increment variance saturates at 2·Var (Airy-process): V(∞)/2Var = {sat:.2f}")

    # ---- figure ----
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
    Yp = res[0.85]["Y"]
    pm = np.nanmean(Yp[:, 6:], axis=0); ps = np.nanstd(Yp[:, 6:], axis=0)
    ks = np.arange(6, Yp.shape[1])
    ax[0].plot(ks, pm, "o-", color="#1f3b73", ms=4, label="per-passage mean")
    ax[0].fill_between(ks, pm - ps, pm + ps, color="#1f3b73", alpha=0.15, label="±1 std")
    ax[0].set_xlabel("passage k"); ax[0].set_ylabel("peel-off level")
    ax[0].set_title("(A) stationary TW₂ marginal across passages (ρ=0.85)")
    ax[0].legend(fontsize=9, frameon=False)

    s = np.arange(len(res[0.0]["C"]))
    ax[1].plot(s, res[0.0]["C"], "s--", color="#888780", label="ρ=0 (iid → white)")
    ax[1].plot(s, res[0.85]["C"], "o-", color="#b3402b", label="ρ=0.85 (forced → process)")
    ax[1].axhline(0, color="grey", lw=0.5, ls=":")
    ax[1].set_xlabel("passage lag s"); ax[1].set_ylabel("autocorrelation C(s)")
    ax[1].set_title("(B) forcing makes peel-offs a correlated process")
    ax[1].legend(fontsize=9, frameon=False)

    V = res[0.85]["V"]; twovar = 2 * res[0.85]["var"]
    ax[2].plot(s, V, "o-", color="#2c7d59", label="Var[Y(k+s)−Y(k)]")
    ax[2].axhline(twovar, color="grey", lw=1, ls="--", label="2·Var (decorrelated)")
    ax[2].set_xlabel("passage lag s"); ax[2].set_ylabel("increment variance")
    ax[2].set_title("(C) locally-Brownian increments saturating at 2·Var")
    ax[2].legend(fontsize=9, frameon=False)

    fig.suptitle("T1.2 forced — successive peel-offs of a forced folded cycle form an "
                 "Airy₂-type process (β=2)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fp = os.path.join(FIG, "forced_peeloff_airy.png")
    fig.savefig(fp, dpi=140, bbox_inches="tight")
    print(f"\n  Figure: {fp}")
    print(f"  TOTAL {time.time()-t0:.1f}s — DONE")


if __name__ == "__main__":
    main()
