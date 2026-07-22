"""
folded_cycle_chartmatch.py — Path A bucket 1: the global noisy chart-matching
-----------------------------------------------------------------------------
The load-bearing estimate for the stochastic K1->K2->K3 matching is the
Berglund-Gentz CONCENTRATION TUBE around the deterministic canard gamma_2, made
UNIFORM in the angular variable theta (JKK's central difficulty: theta cannot be
localized).  Linearising the inner Riccati

    dr2 = ( b r2^2 - a y2 ) dT + eta dB        (a=a(theta), b=b(theta), c=1)

around the attracting canard r2 = -sqrt((a/b) y2): delta = r2 - r2_canard obeys
    d(delta) = -2 sqrt(a b y2) delta dT + eta dB        (Ornstein-Uhlenbeck),
contraction rate kappa = 2 sqrt(a b y2), so the (quasi-stationary) tube variance is

    Var(delta) ≈ eta^2 / (2 kappa) = eta^2 / (4 sqrt(a b y2)).

This is UNIFORM in theta because a(theta), b(theta) are smooth and strictly positive
on the compact circle (so kappa is bounded below by 2 sqrt(min ab y2) > 0).  Escape
occurs when the tube reaches the separatrix (repelling branch, distance
2 sqrt((a/b) y2) from the canard) -> recovers the FW barrier V(Y)=8Y^{3/2}/3.

VERIFY: simulate the inner SDE on the canard for several (a,b) (= frozen theta
values) and eta; measure Var(r2) at y2-checkpoints; check the rescaled variance
    Var(r2) * 4 sqrt(a b y2) / eta^2  ≈ 1,  COLLAPSING across (a,b) (theta-uniformity).
"""
from __future__ import annotations
import os
import numpy as np

C = 1.0


def canard(y2, a, b):
    """attracting-branch canard (leading order) r2 = -sqrt((a/b) y2)."""
    return -np.sqrt((a / b) * np.maximum(y2, 0.0))


def run(a, b, eta, y0=4.0, y_checks=(3.0, 2.0, 1.5, 1.0, 0.6), N=4000, dt=1e-3, seed=0):
    """start ON the canard at y0; drift y2 down; return Var(r2) at each checkpoint."""
    rng = np.random.default_rng(seed)
    r = np.full(N, canard(y0, a, b))
    y = y0
    sdt = np.sqrt(dt)
    checks = sorted(y_checks, reverse=True)
    out = {}
    ci = 0
    n = int((y0 - min(y_checks) + 0.2) / (C * dt))
    for _ in range(n):
        r = r + (b * r**2 - a * y) * dt + eta * sdt * rng.standard_normal(N)
        y = y - C * dt
        while ci < len(checks) and y <= checks[ci]:
            # variance about the (ensemble) mean = transverse tube width
            out[checks[ci]] = float(np.var(r))
            ci += 1
        if ci >= len(checks):
            break
    return out


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    figdir = os.path.join(os.path.dirname(here), "figures")
    os.makedirs(figdir, exist_ok=True)

    geoms = [(1.0, 1.0), (2.0, 1.0), (1.0, 2.0), (1.5, 0.7)]    # = frozen theta values
    etas = [0.08, 0.12]
    y_checks = (3.0, 2.0, 1.5, 1.0, 0.6)

    print("\n=== Path A bucket 1: noise tube around the canard, uniform in theta ===")
    print("  predict Var(r2) ≈ eta^2 / (4 sqrt(a b y2));  rescaled -> 1, collapsing in (a,b)\n")
    print(f"  {'a':>4s}{'b':>4s}{'eta':>6s}   rescaled Var*4 sqrt(ab y2)/eta^2  at y2 = "
          + " ".join(f"{y:.1f}" for y in sorted(y_checks, reverse=True)))
    allres = {}
    for (a, b) in geoms:
        for eta in etas:
            res = run(a, b, eta, y_checks=y_checks, seed=int(100 * a + 10 * b + 100 * eta))
            resc = {y: res[y] * 4 * np.sqrt(a * b * y) / eta**2 for y in res}
            allres[(a, b, eta)] = resc
            print(f"  {a:4.1f}{b:4.1f}{eta:6.2f}   "
                  + "  ".join(f"{resc[y]:.2f}" for y in sorted(res, reverse=True)))

    # collapse quality: across all (a,b,eta), the rescaled variance should ~1
    ys = sorted(y_checks, reverse=True)
    print("\n  collapse (mean +/- std of rescaled Var across all 8 (a,b,eta) cells):")
    for y in ys:
        vals = np.array([allres[k][y] for k in allres if y in allres[k]])
        print(f"    y2={y:.1f}:  {vals.mean():.2f} +/- {vals.std():.2f}")
    print("\n  => rescaled variance ~ O(1) and COLLAPSES across (a,b)=frozen-theta:")
    print("     the BG concentration tube Var=eta^2/4sqrt(ab y2) is theta-UNIFORM")
    print("     (deviation grows as y2->0 where the quasi-stationary OU approx breaks,")
    print("     i.e. the fold/inner region where the K2 escape analysis takes over).\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    for (a, b) in geoms:
        for eta in etas:
            resc = allres[(a, b, eta)]
            yy = sorted(resc, reverse=True)
            ax[0].plot(yy, [resc[y] for y in yy], "o-", ms=4, alpha=0.8,
                       label=fr"$(a,b)$=({a:g},{b:g}),$\eta$={eta}")
    ax[0].axhline(1.0, color="k", lw=1.2, ls="--", label="OU prediction = 1")
    ax[0].set_xlabel(r"$y_2$ (toward fold $\to$)"); ax[0].invert_xaxis()
    ax[0].set_ylabel(r"Var$(r_2)\cdot 4\sqrt{ab\,y_2}/\eta^2$")
    ax[0].set_title("Noise tube collapses to OU prediction (uniform in $\\theta$)")
    ax[0].legend(fontsize=7, frameon=False, ncol=2); ax[0].grid(alpha=0.3)

    # raw tube width vs the prediction for one cell
    a, b, eta = 1.0, 1.0, 0.12
    res = run(a, b, eta, y_checks=tuple(np.linspace(3.5, 0.5, 12)), seed=1)
    yy = np.array(sorted(res, reverse=True))
    ax[1].plot(yy, [res[y] for y in yy], "o", ms=5, color="C0", label="measured Var$(r_2)$")
    ax[1].plot(yy, eta**2 / (4 * np.sqrt(a * b * yy)), "-", color="C3",
               label=r"$\eta^2/4\sqrt{ab\,y_2}$ (OU tube)")
    ax[1].set_xlabel(r"$y_2$"); ax[1].invert_xaxis(); ax[1].set_ylabel(r"Var$(r_2)$")
    ax[1].set_title(r"Tube widens as $y_2\to0$ (contraction $\to0$ at the fold)")
    ax[1].legend(fontsize=9, frameon=False); ax[1].grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(figdir, "folded_cycle_chartmatch.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print(f"  figure -> figures/folded_cycle_chartmatch.png\n")


if __name__ == "__main__":
    main()
