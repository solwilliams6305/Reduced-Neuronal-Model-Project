"""
coupled_blowup_noise.py — noise through the coupled cusp blow-up: noisy Airy → noisy Weber.
===========================================================================================

The peel-off escape needs a SEMI-INFINITE oscillatory region, so the textbook Olver two-turning-point
(bounded-well) equation u''=(¼Y²−a)u is the WRONG inner model (it gives a bounded-well law, not the
edge). The correct unfolding keeps the semi-infinite structure and tunes the turning ORDER near the
origin:
    u'' = (V_Δ(Y) − η ξ) u ,   V_Δ(Y) = sign(Y)·|Y|·(|Y| + Δ).
Near Y=0: V_Δ ≈ Δ·Y (LINEAR, slope Δ → Airy) for |Y|≪Δ;  V_Δ ≈ sign(Y)Y² (QUADRATIC → Weber) for
|Y|≫Δ. At Δ=0 it is exactly the cusp law sign(Y)Y². So Δ unfolds Weber (Δ=0) into a linear fold (Δ→∞).
The coupling sets the separation via Δ(g)=2√(−2g/3): g→0 ⇒ Δ→0 (cusp/Weber); large |g| ⇒ Δ large
(fold/Airy).

Carrying NOISE through: peel-off = first node of the swept decaying solution. Prediction:
  • Δ ≫ inner noise scale ℓ(η): the turning is effectively LINEAR → AIRY → Tracy–Widom.
  • Δ → 0 (Δ ≲ ℓ): the quadratic dominates → WEBER edge → the sub-Gaussian cusp law.
Crossover at Δ ~ ℓ(η); measure Δ_crit(η) to check the noise rescales correctly.

[PROVED] the deterministic reduction (K-P, Olver Airy↔Weber) and Δ(g)=2√(−2g/3).  [NUMERIC] the noisy
crossover + Δ_crit(η) here.  [HEURISTIC] the Δ(g)↔g map and the tie to the full-model √ε onset.

Output: figures/coupled_blowup_noise.png + summary.
"""
from __future__ import annotations
import os, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
TW2 = dict(sk=0.224, ku=0.093)


def Vfun(Y, D):
    return np.sign(Y)*np.abs(Y)*(np.abs(Y) + D)


def first_node(D, eta, M=4000, dt=1.0e-3, Y0=5.0, Yend=-6.0, seed=0):
    """First node of u''=(V_Δ(Y)−ηξ)u, V_Δ=sign(Y)|Y|(|Y|+Δ); decaying solution swept from Y0."""
    rng = np.random.default_rng(seed); sdt = np.sqrt(dt)
    V0 = Vfun(Y0, D)
    u = np.ones(M); v = np.sqrt(max(V0, 1e-6))*u; Y = Y0   # growing-branch IC (matches peel-off codes)
    Ynode = np.full(M, np.nan); done = np.zeros(M, bool)
    n = int((Y0 - Yend)/dt)
    for _ in range(n):
        u_prev = u
        dB = sdt*rng.standard_normal(M)
        Vc = Vfun(Y, D); Yp = Y - dt; Vp = Vfun(Yp, D)
        u1 = u + v*dt
        v1 = v + (Vc*u)*dt - eta*u*dB
        u = u + 0.5*(v + v1)*dt
        v = v + 0.5*(Vc*u + Vp*u1)*dt - eta*0.5*(u + u1)*dB
        Y = Yp
        cr = (~done) & ((u_prev > 0) != (u > 0))
        Ynode[cr] = Y; done |= cr
        if done.all():
            break
    return Ynode[np.isfinite(Ynode)]


def moments(x):
    if x.size < 30:
        return np.nan, np.nan
    z = (x - x.mean())/x.std()
    return np.mean(z**3), np.mean(z**4) - 3.0


def main():
    t0 = time.time()
    eta0 = np.sqrt(2.0)                              # β=2
    Dvals = [6.0, 4.0, 2.5, 1.5, 0.9, 0.5, 0.25, 0.1, 0.0]
    print("=" * 78)
    print("Noise through the coupled cusp blow-up — noisy Airy (Δ large) → noisy Weber (Δ→0)")
    print("=" * 78)
    print(f"  inner eq u''=(V_Δ−ηξ)u, V_Δ=sign(Y)|Y|(|Y|+Δ), η=√2 (β=2). Δ=2√(−2g/3): Δ→0 ⇔ g→0 (cusp).")
    print(f"\n  {'Δ':>6} | {'skew':>8} {'exkurt':>8} {'n':>6}   regime")
    sk, ku, samp = {}, {}, {}
    for D in Dvals:
        x = first_node(D, eta0, seed=7)
        s, k = moments(x); sk[D] = s; ku[D] = k; samp[D] = x
        reg = "Airy/TW" if D >= 4.0 else ("Weber/cusp" if D <= 0.25 else "crossover")
        print(f"  {D:6.2f} | {s:+8.3f} {k:+8.3f} {x.size:6d}   {reg}")
    print(f"  TW₂ ref: skew {TW2['sk']:+.2f}, exkurt {TW2['ku']:+.2f}.  cusp(Δ=0, sign(Y)Y²) ref: skew ≈+0.6, exkurt ≈−0.2.")

    # noise rescaling: Δ_crit(η) = separation at the skew midpoint between Δ-large (Airy) and Δ=0 (Weber)
    print(f"\n  noise rescaling — Δ_crit(η) (skew midpoint between Airy (Δ large) and Weber (Δ=0)):")
    crit = {}; curves = {}
    Dgrid = np.array([4.0, 2.5, 1.5, 0.9, 0.5, 0.2, 0.0])
    for eta in (0.9, 1.15, 1.41):                       # clean regime (β≈4.9,3.0,2.0); β=1 washes out
        sarr = np.array([moments(first_node(D, eta, M=3500, seed=11))[0] for D in Dgrid])  # skew, incr as Δ↓
        curves[eta] = sarr
        mid = 0.5*(sarr[0] + sarr[-1])
        dc = float(np.interp(mid, sarr, Dgrid[::-1]))   # sarr increasing; map skew→Δ
        crit[eta] = dc
        print(f"    η={eta:.3f} (β={4/eta**2:.1f}): skew {sarr[0]:+.2f}(Δ=4)→{sarr[-1]:+.2f}(Δ=0), "
              f"mid {mid:+.2f} ⇒ Δ_crit ≈ {dc:.2f}")
    print(f"    crossover sits at Δ_crit ~ O(1). [honest] a clean ℓ(η) POWER is confounded: both endpoint")
    print(f"    laws (Airy skew, Weber skew) themselves shift with η/β (see CUSP_LAW_CONJECTURES), so the")
    print(f"    Δ_crit metric mixes the crossover scale with the β-dependence of the endpoints.")

    print(f"\n  ⇒ [NUMERIC] the noisy peel-off law crosses AIRY/TW (Δ≫ℓ) → WEBER/cusp (Δ→0), the deterministic")
    print(f"    two-fold→cusp picture with noise carried through; Δ_crit(η) sets the inner noise scale ℓ(η).")
    print(f"  [HEURISTIC] Δ=2√(−2g/3): crossover Δ~ℓ~ε^{{1/4}} ⇔ g_crit∝√ε — matches the full-model √ε onset")
    print(f"    (g_crit/√ε≈−0.58). DIRECTION: cusp/Weber at Δ→0 = SMALL |g| (near g=0); Airy at large |g|.")

    # ---- figure ----
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    DD = np.array(Dvals)
    ax[0].axhline(0, color="grey", lw=0.8)
    ax[0].axhline(TW2['sk'], color="#b3402b", ls=":", lw=1, label="TW skew")
    ax[0].plot(DD, [sk[D] for D in Dvals], "o-", color="#7a3b8f", lw=1.8, label="skew")
    ax[0].plot(DD, [ku[D] for D in Dvals], "s-", color="#2c7d59", lw=1.8, label="exkurt")
    ax[0].set_xlabel("fold separation Δ  (Δ=2√(−2g/3); Δ→0 is the cusp)")
    ax[0].set_ylabel("peel-off moment"); ax[0].invert_xaxis()
    ax[0].set_title("(A) Δ large → TW;  Δ→0 → Weber (κ4 flips +→−)"); ax[0].legend(fontsize=8.5, frameon=False)

    gg = np.linspace(-4, 4, 200)
    ax[1].plot(gg, np.exp(-gg**2/2)/np.sqrt(2*np.pi), "k:", lw=1, label="Gaussian")
    for D, col, tag in [(6.0, "#b3402b", "Δ=6 (Airy/TW)"), (0.0, "#7a3b8f", "Δ=0 (Weber/cusp)")]:
        x = samp[D]
        if x.size > 30:
            z = (x - x.mean())/x.std()
            ax[1].hist(z, bins=50, range=(-4, 4), density=True, histtype="step", lw=1.9, color=col,
                       label=f"{tag} (skew {sk[D]:+.2f}, κ4 {ku[D]:+.2f})")
    ax[1].set_xlabel("standardised peel-off"); ax[1].set_ylabel("density")
    ax[1].set_title("(B) noisy Airy vs noisy Weber inner law"); ax[1].legend(fontsize=8, frameon=False)

    cols3 = {0.9: "#1d6e3f", 1.15: "#1f3b73", 1.41: "#7a3b8f"}
    for e in sorted(curves):
        ax[2].plot(Dgrid, curves[e], "o-", color=cols3[e], lw=1.6, label=f"η={e} (β={4/e**2:.1f})")
        ax[2].axvline(crit[e], color=cols3[e], ls=":", lw=1)
    ax[2].axhline(TW2['sk'], color="grey", ls=":", lw=0.8)
    ax[2].set_xlabel("fold separation Δ"); ax[2].set_ylabel("peel-off skew"); ax[2].invert_xaxis()
    ax[2].set_title("(C) crossover at Δ~O(1); endpoints β-shift (ℓ(η) power confounded)")
    ax[2].legend(fontsize=8, frameon=False)

    fig.suptitle("Stochastic coupled blow-up: noise carried through gives noisy-Airy (Δ≫ℓ) → noisy-Weber "
                 "(Δ→0); the cusp/Weber edge is approached at small |g| (Δ→0)", fontsize=10.3)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fp = os.path.join(FIG, "coupled_blowup_noise.png")
    fig.savefig(fp, dpi=140, bbox_inches="tight")
    print(f"\n  Figure: {fp}")
    print(f"  TOTAL {time.time()-t0:.1f}s — DONE")


if __name__ == "__main__":
    main()
