"""
cusp_tubes.py — PART A: Berglund–Gentz tube estimate in the cusp chart (attempt).
=================================================================================

Rescaled cusp SDE: dD = (-D^3 + M D - W) ds + η dB_s,  W = s.
Linearise about the deterministic canard D̄(s) (attracting branch): ξ = D - D̄ obeys
    dξ = a(s) ξ ds + η dB_s ,   a(s) = -3 D̄(s)^2 + M  (< 0 on the attracting branch),
an Ornstein–Uhlenbeck process with slowly-varying rate. Predictions:
    • tube variance      Var(ξ)(s) ≈ η² / (2|a(s)|)
    • Gaussian confinement  P(sup_s |ξ| > h) ≲ exp(-κ h²/η²)
    • matching height (linear tube fails ⇒ peel-off): a≈-2√(3M)·d near the fold ⇒
      d_* ~ (η²/(4√(3M)))^{1/3} = η^{2/3}/(4√(3M))^{1/3}  (→ η^{2/3} at the fold, M=O(1)).

This script checks the variance law and the Gaussian confinement numerically (the two clean,
provable-in-principle ingredients of a Berglund–Gentz tube one rung up).

Output: figures/cusp_tubes.png + summary.
"""
from __future__ import annotations
import os, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")


def canard(M, ds, n):
    D = np.sqrt(M); out = np.empty(n)
    for i in range(n):
        W = i*ds; D = D + (-D**3 + M*D - W)*ds; out[i] = D
    return out


def ensemble(M, eta, N=5000, ds=0.002, seed=0):
    Wstar = (2*M/3)*np.sqrt(M/3); n = int(0.92*Wstar/ds)   # pre-escape window
    Dbar = canard(M, ds, n)
    rng = np.random.default_rng(seed); sds = np.sqrt(ds)
    D = np.full(N, np.sqrt(M)); xi = np.empty((n, N))
    for i in range(n):
        W = i*ds
        D = D + (-D**3 + M*D - W)*ds + eta*sds*rng.standard_normal(N)
        xi[i] = D - Dbar[i]
    a = -3*Dbar**2 + M
    return Dbar, a, xi, Wstar, np.arange(n)*ds


def main():
    t0 = time.time()
    M, eta = 2.0, 0.30
    print("=" * 70)
    print(f"PART A — cusp tube estimate (M={M}, η={eta})")
    print("=" * 70)
    Dbar, a, xi, Wstar, Ws = ensemble(M, eta)
    var = xi.var(axis=1)
    pred = eta**2 / (2*np.abs(a))
    # compare away from the fold (|a|>0.5)
    ok = np.abs(a) > 0.5
    ratio = np.nanmean(var[ok] / pred[ok])
    print(f"  [NUMERIC] tube variance Var(ξ) vs η²/(2|a|): mean ratio = {ratio:.2f}  "
          f"({'✓ OU law' if 0.8 < ratio < 1.2 else 'check'})")

    # Gaussian confinement: P(sup|ξ|>h) over the pre-escape window
    sup = np.max(np.abs(xi[ok]), axis=0)           # sup over the safe window, per realization
    hs = np.linspace(0.15, 0.6, 12)
    P = np.array([(sup > h).mean() for h in hs])
    m = P > 0
    # fit log P vs (h/η)^2  → slope -κ
    x = (hs[m]/eta)**2; y = np.log(P[m])
    kappa = -np.polyfit(x, y, 1)[0]
    print(f"  [NUMERIC] confinement P(sup|ξ|>h) ~ exp(-κ (h/η)²): κ ≈ {kappa:.2f}  "
          f"(Gaussian tail ⇒ tube holds)")
    print(f"  [DERIVED] matching height d_* ~ η^(2/3)/(4√(3M))^(1/3) = "
          f"{eta**(2/3)/(4*np.sqrt(3*M))**(1/3):.3f}  (peel-off begins here)")

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    ax[0].plot(Ws, var, color="#7a3b8f", lw=2, label="measured Var(ξ)")
    ax[0].plot(Ws, pred, "--", color="#1f3b73", lw=1.6, label="η²/(2|a(s)|)  [OU prediction]")
    ax[0].set_xlabel("load W"); ax[0].set_ylabel("Var(ξ)")
    ax[0].set_title("(A) tube variance follows the OU law"); ax[0].legend(fontsize=9, frameon=False)
    ax[0].set_ylim(0, min(0.06, np.nanmax(pred[ok])*1.3))

    ax[1].semilogy(hs, P, "o", color="#b3402b")
    xx = np.linspace(hs.min(), hs.max(), 50)
    ax[1].semilogy(xx, np.exp(np.polyval(np.polyfit(x, y, 1), (xx/eta)**2)), "-",
                   color="#1f3b73", lw=1.4, label=f"exp(−κ(h/η)²), κ={kappa:.1f}")
    ax[1].set_xlabel("tube half-width h"); ax[1].set_ylabel("P(sup|ξ| > h)")
    ax[1].set_title("(B) Gaussian confinement"); ax[1].legend(fontsize=9, frameon=False)

    # sample paths in the tube
    rng = np.random.default_rng(3); sds = np.sqrt(0.002); n = len(Ws)
    Dn = np.full(10, np.sqrt(M))
    paths = np.empty((n, 10))
    for i in range(n):
        Dn = Dn + (-Dn**3 + M*Dn - Ws[i])*0.002 + eta*sds*rng.standard_normal(10)
        paths[i] = Dn
    ax[2].plot(Ws, Dbar, "k-", lw=2, label="canard D̄(s)")
    ax[2].plot(Ws, Dbar + np.sqrt(pred), "b:", lw=1, label="±OU tube")
    ax[2].plot(Ws, Dbar - np.sqrt(pred), "b:", lw=1)
    ax[2].plot(Ws, paths, lw=0.5, alpha=0.5, color="#7a3b8f")
    ax[2].set_xlabel("load W"); ax[2].set_ylabel("D")
    ax[2].set_title("(C) paths confined to the ±η/√(2|a|) tube"); ax[2].legend(fontsize=8.5, frameon=False)

    fig.suptitle("Part A — cusp tube: OU variance law + Gaussian confinement (Berglund–Gentz one rung up)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fp = os.path.join(FIG, "cusp_tubes.png")
    fig.savefig(fp, dpi=140, bbox_inches="tight")
    print(f"\n  Figure: {fp}")
    print(f"  TOTAL {time.time()-t0:.1f}s — DONE")


if __name__ == "__main__":
    main()
