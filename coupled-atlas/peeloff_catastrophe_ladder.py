"""
peeloff_catastrophe_ladder.py — Step 1: the genuine diffraction-catastrophe ladder.
====================================================================================

Berry–Upstill diffraction catastrophes: the A_n cuspoid has diffraction integral
∫exp(t^{n+1}) and its function solves an n-th-order linear ODE:

    fold  (A₂, Airy)       u''   = Y u        (2nd order)  → the paper → Tracy–Widom
    cusp  (A₃, Pearcey)    u'''  = Y u        (3rd order)  → ???
    swallowtail (A₄)       u'''' = Y u        (4th order)  → ???

So the genuine catastrophe ladder is  u^{(n)} = (Y - η ξ) u ,  Y = Y0 - T, with noise
multiplying u exactly as in the paper's fold case. The canard is the solution that grows
along the positive real characteristic root (u^{(j)}/u = Y^{j/n} at large Y); the first
node of u, swept downward, is the peel-off level. n=2 must reproduce TW; the question is
whether n=3 (the cusp / Pearcey) gives a DIFFERENT edge universality class.

This is the faithful "genuine cusp" test (the cusp = the Pearcey ODE), complementing the
multicritical/higher-order-Airy probe in peeloff_cusp_ladder.py. Relating it to the specific
coupled-FHN cusp (Kristiansen–Pedersen) is the remaining modelling step. Heun/Stratonovich,
generalising peeloff_tw_validation.py (which is the n=2 special case).

Output: figures/peeloff_catastrophe_ladder.png + printed summary.
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

TW2 = dict(skew=0.2241, kurt=0.0934)        # fold reference (n=2, β=2)
NAMES = {2: "fold (Airy)", 3: "cusp (Pearcey)", 4: "swallowtail"}


def ladder_peeloff(n, eta, N=10000, Y0=5.0, dt=1e-3, Ymin=-7.0, seed=0):
    """First node of u^{(n)} = (Y - η ξ)u (Heun). n=2 is the fold → TW."""
    rng = np.random.default_rng(seed); sdt = np.sqrt(dt)
    U = [np.full(N, Y0 ** (j / n)) for j in range(n)]      # canard IC: u^{(j)}/u = Y0^{j/n}
    Y = Y0
    Yz = np.full(N, np.nan); done = np.zeros(N, bool)
    nstep = int((Y0 - Ymin) / dt)
    for _ in range(nstep):
        al = ~done
        if not al.any():
            break
        dB = sdt * rng.standard_normal(N); Yp = Y - dt
        # predictor
        U1 = [None] * n
        for j in range(n - 1):
            U1[j] = U[j] + U[j + 1] * dt
        U1[n - 1] = U[n - 1] + (Y * U[0]) * dt - eta * U[0] * dB
        # corrector
        Un = [None] * n
        for j in range(n - 1):
            Un[j] = U[j] + 0.5 * (U[j + 1] + U1[j + 1]) * dt
        Un[n - 1] = U[n - 1] + 0.5 * (Y * Un[0] + Yp * U1[0]) * dt \
            - eta * 0.5 * (Un[0] + U1[0]) * dB
        U = Un; Y = Yp
        cr = al & (U[0] < 0.0); Yz[cr] = Y; done |= cr
    return Yz[np.isfinite(Yz)]


def moments(x):
    m = x.mean(); d = x - m; v = np.mean(d**2)
    return m, np.sqrt(v), np.mean(d**3) / v**1.5, np.mean(d**4) / v**2 - 3.0, x.size


def main():
    t0 = time.time()
    eta = np.sqrt(2.0)
    print("=" * 72)
    print("Step 1 — the diffraction-catastrophe ladder  u^(n) = (Y-ηξ)u")
    print("=" * 72)
    print(f"  fold (n=2) reference: TW₂ skew {TW2['skew']:.3f}, exkurt {TW2['kurt']:.3f}")
    print(f"  {'n':>3} {'catastrophe':>16} | {'mean':>8}{'std':>8}{'skew':>8}{'kurt':>8} {'N':>7}")
    samples = {}
    for n in (2, 3, 4):
        x = ladder_peeloff(n, eta, seed=4 + n)
        samples[n] = x
        m, s, g1, g2, cnt = moments(x)
        print(f"  {n:>3} {NAMES[n]:>16} | {m:8.3f}{s:8.3f}{g1:8.3f}{g2:8.3f} {cnt:7d}")

    se = np.sqrt(6.0 / min(len(v) for v in samples.values()))
    sk2 = moments(samples[2])[2]; sk3 = moments(samples[3])[2]; sk4 = moments(samples[4])[2]
    print(f"\n  (skew MC s.e. ≈ {se:.3f})")
    print(f"  n=2 reproduces TW₂: {'PASS' if abs(sk2-TW2['skew'])<3*se else 'CHECK'} "
          f"({sk2:.3f} vs {TW2['skew']:.3f})")
    print(f"  cusp (n=3) departs from the fold: "
          f"{'YES — new edge class' if abs(sk3-sk2)>3*se else 'no'}  (Δskew={sk3-sk2:+.3f})")
    print(f"  swallowtail (n=4) departs further: {'YES' if abs(sk4-sk2)>3*se else 'no'}  "
          f"(Δskew={sk4-sk2:+.3f})")

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    ns = [2, 3, 4]; sk = [moments(samples[n])[2] for n in ns]; ku = [moments(samples[n])[3] for n in ns]
    ax[0].axhline(TW2['skew'], color="#b3402b", lw=1, ls="--", label="TW₂ (fold)")
    ax[0].fill_between([1.8, 4.2], TW2['skew']-3*se, TW2['skew']+3*se, color="#b3402b", alpha=0.12)
    ax[0].plot(ns, sk, "o-", color="#1f3b73", ms=8)
    for n, y in zip(ns, sk):
        ax[0].annotate(NAMES[n].split()[0], (n, y), textcoords="offset points", xytext=(6, 6), fontsize=9)
    ax[0].set_xticks(ns); ax[0].set_xlabel("catastrophe order n  (2 fold · 3 cusp · 4 swallowtail)")
    ax[0].set_ylabel("peel-off skewness")
    ax[0].set_title("(A) escape law climbs the catastrophe ladder")
    ax[0].legend(fontsize=9, frameon=False)

    cols = {2: "#b3402b", 3: "#1f3b73", 4: "#7a3b8f"}
    gg = np.linspace(-4, 4, 200)
    ax[1].plot(gg, np.exp(-gg**2/2)/np.sqrt(2*np.pi), "k:", lw=1, label="Gaussian")
    for n in ns:
        x = samples[n]; z = (x - x.mean()) / x.std()
        ax[1].hist(z, bins=70, range=(-4, 4), density=True, histtype="step", lw=1.8,
                   color=cols[n], label=f"{NAMES[n]} (skew {moments(x)[2]:.2f})")
    ax[1].set_xlabel("standardized peel-off"); ax[1].set_ylabel("density")
    ax[1].set_title("(B) fold (TW) vs cusp vs swallowtail laws")
    ax[1].legend(fontsize=8.5, frameon=False)

    fig.suptitle("Step 1 — noise through the diffraction-catastrophe ladder "
                 "(fold→Airy→TW; cusp→Pearcey; …)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fp = os.path.join(FIG, "peeloff_catastrophe_ladder.png")
    fig.savefig(fp, dpi=140, bbox_inches="tight")
    print(f"\n  Figure: {fp}")
    print(f"  TOTAL {time.time()-t0:.1f}s — DONE")


if __name__ == "__main__":
    main()
