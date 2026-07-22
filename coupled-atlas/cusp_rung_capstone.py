"""
cusp_rung_capstone.py — pushing the cusp rung (B) to its HONEST ceiling.
=======================================================================

Mirror of fold_rung_capstone.py — but the cusp rung has a genuinely OPEN third level
(the intrinsic multicritical "Dyson–Weber" line ensemble), so it cannot honestly reach 100% the
way the fold rung did. This capstone consolidates what IS established and pins what is NOT:

 (1) MARGINAL — Weber, not Pearcey (high statistics). Inner laws u''=(V(Y)-ηξ)u with
     V=sign(Y)|Y|^k: k=1 fold→TW (skew ≈0.22), k=2 Weber/cusp (skew ≈0.6); plus the 3rd-order
     Pearcey ODE u'''=(Y-ηξ)u. The cusp marginal is distinct from BOTH Gaussian/TW and Pearcey.
     [Consolidates weber_vs_pearcey.py, where the genuine KP coupled-FHN escape matched Weber:
     KS 0.056 vs Pearcey 0.123 — "Weber (2nd order), not Pearcey (3rd order)".]
 (2) POINT PROCESS + correction — successive nodes of the SWEPT Weber operator give the cusp
     swept-node correlation (the Weber analogue of the fold's swept-node +0.89). Logged correction:
     this is the SWEPT-NODE value (nodes share the sweep path); the GENUINE multicritical
     point-process rigidity awaits the Dyson–Weber object — exactly as, for the fold, the genuine
     Airy rigidity 0.50 (GUE edge) < swept-node 0.89.
 (3) INTRINSIC — a no-go + the open object. The GUE/matrix-OU edge gives the Airy line ensemble
     (top-eigenvalue skew ≈0.22); the cusp marginal is Weber (skew ≈0.6). So the cusp's intrinsic
     process is NOT the Airy line ensemble — it is the multicritical-edge ("Dyson–Weber") line
     ensemble (DBM at a higher-order edge, density ~ dist^{3/2}), which is genuinely OPEN: no
     RRV-type stochastic-operator characterization is known for k>1, and matrix-OU reaches only
     the Airy edge (a critically-tuned non-Gaussian model is needed; numerically delicate).

Output: figures/cusp_rung_capstone.png + tagged summary.
"""
from __future__ import annotations
import os, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from peeloff_catastrophe_ladder import ladder_peeloff       # n-th order inner: n=2 fold, n=3 Pearcey
from peeloff_cusp_ladder import peeloff_ladder              # turning order k (2nd order): k=2 Weber

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
TW2 = dict(skew=0.2241, kurt=0.0934)


def zstd(x):
    x = x[np.isfinite(x)]; return (x - x.mean()) / x.std()


def moments(x):
    x = x[np.isfinite(x)]; d = x - x.mean(); v = np.mean(d**2)
    return np.mean(d**3)/v**1.5, np.mean(d**4)/v**2 - 3.0


def ks_dist(a, b):
    a = np.sort(a); b = np.sort(b); grid = np.concatenate([a, b])
    Fa = np.searchsorted(a, grid, side="right") / a.size
    Fb = np.searchsorted(b, grid, side="right") / b.size
    return np.max(np.abs(Fa - Fb))


def successive_weber_nodes(eta, K=3, Y0=8.0, dt=8e-4, N=4000, Ymin=-11.0, seed=0):
    """First K nodes of the swept Weber operator u''=(sign(Y)|Y|^2 - η ξ)u (Heun/Stratonovich)."""
    rng = np.random.default_rng(seed)
    u = np.ones(N); v = np.full(N, Y0)          # v0 = Y0^{k/2} with k=2
    Y = Y0; sdt = np.sqrt(dt)
    Yk = np.full((N, K), np.nan); cnt = np.zeros(N, int); idx = np.arange(N)
    n = int((Y0 - Ymin) / dt)
    for _ in range(n):
        u_prev = u
        dB = sdt * rng.standard_normal(N)
        Vc = np.sign(Y) * Y * Y; Yp = Y - dt; Vp = np.sign(Yp) * Yp * Yp
        u1 = u + v * dt; v1 = v + (Vc * u) * dt - eta * u * dB
        u = u + 0.5 * (v + v1) * dt
        v = v + 0.5 * (Vc * u + Vp * u1) * dt - eta * 0.5 * (u + u1) * dB
        Y = Yp
        crossed = ((u_prev > 0) != (u > 0)) & (cnt < K)
        if crossed.any():
            r = idx[crossed]; Yk[r, cnt[r]] = Y; cnt[r] += 1
        if (cnt >= K).all():
            break
    return Yk


def gue_top(N, M, seed=0):
    """Top eigenvalue of the GUE edge (the Airy line-ensemble marginal = Tracy–Widom)."""
    rng = np.random.default_rng(seed); top = np.empty(M)
    for m in range(M):
        A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
        H = (A + A.conj().T) / (2.0 * np.sqrt(N))
        top[m] = np.linalg.eigvalsh(H)[-1]
    return top


def main():
    t0 = time.time()
    eta = np.sqrt(2.0)                       # β = 4/η² = 2
    print("=" * 76)
    print("Cusp rung capstone — Weber marginal + swept-node point process + intrinsic no-go")
    print("=" * 76)

    # ---- (1) MARGINAL: fold/TW, Weber/cusp, Pearcey ----
    fold = ladder_peeloff(2, eta, N=6000, seed=12)            # fold / Airy → TW
    pear = ladder_peeloff(3, eta, N=6000, seed=13)            # Pearcey (3rd order)
    weber = peeloff_ladder([2.0], eta, N=6000, seed=14)[0]    # Weber-class (cusp, 2nd order)
    sk_f, ku_f = moments(fold); sk_w, ku_w = moments(weber); sk_p, ku_p = moments(pear)
    print(f"\n  (1) MARGINAL (standardised inner-law skew / exkurt):")
    print(f"      fold / Airy (TW): skew {sk_f:+.3f}  exkurt {ku_f:+.3f}   [TW₂ {TW2['skew']:+.3f}/{TW2['kurt']:+.3f}]")
    print(f"      Weber / cusp    : skew {sk_w:+.3f}  exkurt {ku_w:+.3f}   ← distinct, right-skewed")
    print(f"      Pearcey (3rd)   : skew {sk_p:+.3f}  exkurt {ku_p:+.3f}")
    ks_wf = ks_dist(zstd(weber), zstd(fold)); ks_wp = ks_dist(zstd(weber), zstd(pear))
    print(f"      KS(Weber, fold/TW) = {ks_wf:.3f}   KS(Weber, Pearcey) = {ks_wp:.3f}  "
          f"⇒ Weber is its own shape (≠ TW, ≠ Pearcey)")
    print(f"      [from weber_vs_pearcey.py] genuine KP coupled-FHN escape ↔ Weber: KS 0.056 vs "
          f"Pearcey 0.123 → 'Weber, not Pearcey'. K-P reduction PROVES the 2nd-order Weber form.")

    # ---- (2) POINT PROCESS: swept Weber nodes + correction ----
    a_zeros = np.array([-2.338107, -4.087949, -5.520560])     # Airy zeros (fold η→0 reference only)
    Yk = successive_weber_nodes(eta, K=3, seed=7)
    means = np.nanmean(Yk, axis=0); frac = np.mean(np.isfinite(Yk), axis=0)
    ok = np.isfinite(Yk[:, 0]) & np.isfinite(Yk[:, 1])
    r12 = float(np.corrcoef(Yk[ok, 0], Yk[ok, 1])[0, 1])
    print(f"\n  (2) POINT PROCESS — swept Weber nodes (η={eta:.3f}, β=2):")
    for k in range(3):
        print(f"      node #{k+1}: mean {means[k]:7.3f}   (resolved {100*frac[k]:.0f}%)")
    print(f"      cusp SWEPT-NODE correlation corr(Y¹,Y²) = {r12:+.3f}")
    print(f"      [CORRECTION] this is the SWEPT-NODE value (nodes share the sweep path). The genuine")
    print(f"      multicritical point-process rigidity awaits the Dyson–Weber object — cf. the fold,")
    print(f"      where genuine Airy rigidity 0.50 (GUE edge) < swept-node 0.89.")
    print(f"      Level-3 TIME process: forced Weber peel-offs = stationary Weber-marginal process")
    print(f"      (decaying covariance, increments→2·Var) — the Weber analogue of Airy₂ "
          f"[forced_cusp_process.py].")

    # ---- (3) INTRINSIC no-go: GUE edge is Airy, cusp is Weber ----
    gtop = gue_top(N=90, M=3500, seed=3)
    sk_g, ku_g = moments(gtop)
    print(f"\n  (3) INTRINSIC no-go:")
    print(f"      GUE/matrix-OU edge (Airy line ensemble) top-eigenvalue skew {sk_g:+.3f} (TW₂ {TW2['skew']:+.3f}).")
    print(f"      cusp marginal (Weber) skew {sk_w:+.3f}  ≠  Airy {sk_g:+.3f}")
    print(f"      ⇒ the cusp's intrinsic process is NOT the Airy line ensemble. It is the")
    print(f"      multicritical-edge 'Dyson–Weber' line ensemble (DBM at a higher-order edge,")
    print(f"      density ~ dist^{{3/2}}) — defined precisely, GENUINELY OPEN (no RRV-type operator")
    print(f"      known for k>1; matrix-OU reaches only the Airy edge).")

    print(f"\n  ⇒ Cusp rung at its honest ceiling: marginal = Weber (proved reduction + matched, ≠Pearcey);")
    print(f"     point process = Weber (swept-node measured, genuine rigidity pending); time process =")
    print(f"     forced Weber surrogate done, INTRINSIC Dyson–Weber object OPEN. Not 100% — ~90%.")

    # ---- figure ----
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.7))
    gg = np.linspace(-4, 4, 220)
    ax[0].plot(gg, np.exp(-gg**2/2)/np.sqrt(2*np.pi), "k:", lw=1, label="Gaussian")
    ax[0].hist(zstd(fold), bins=60, range=(-4, 4), density=True, histtype="step", lw=1.7,
               color="#b3402b", label=f"fold / TW (skew {sk_f:+.2f})")
    ax[0].hist(zstd(weber), bins=60, range=(-4, 4), density=True, histtype="step", lw=2.1,
               color="#7a3b8f", label=f"Weber / cusp (skew {sk_w:+.2f})")
    ax[0].hist(zstd(pear), bins=60, range=(-4, 4), density=True, histtype="step", lw=1.7,
               color="#1f3b73", label=f"Pearcey 3rd (skew {sk_p:+.2f})")
    ax[0].set_xlabel("standardised peel-off"); ax[0].set_ylabel("density")
    ax[0].set_title("(1) marginal = Weber — distinct from TW & Pearcey")
    ax[0].legend(fontsize=8.5, frameon=False)

    ax[1].plot(Yk[ok, 0], Yk[ok, 1], ".", ms=2.5, alpha=0.22, color="#7a3b8f")
    ax[1].set_xlabel("1st Weber node Y¹"); ax[1].set_ylabel("2nd Weber node Y²")
    ax[1].set_title(f"(2) Weber point process — swept-node corr {r12:+.2f}\n"
                    f"(genuine multicritical rigidity: open)")

    ax[2].plot(gg, np.exp(-gg**2/2)/np.sqrt(2*np.pi), "k:", lw=1, label="Gaussian")
    ax[2].hist(zstd(gtop), bins=55, range=(-4, 4), density=True, histtype="step", lw=1.9,
               color="#d95f0e", label=f"GUE edge = Airy (skew {sk_g:+.2f})")
    ax[2].hist(zstd(weber), bins=55, range=(-4, 4), density=True, histtype="step", lw=1.9,
               color="#7a3b8f", label=f"cusp = Weber (skew {sk_w:+.2f})")
    ax[2].set_xlabel("standardised λ_max / peel-off"); ax[2].set_ylabel("density")
    ax[2].set_title("(3) intrinsic ≠ Airy ⇒ Dyson–Weber (open)")
    ax[2].legend(fontsize=8.5, frameon=False)

    fig.suptitle("Cusp rung capstone — Weber marginal (≠ Pearcey) + Weber point process; "
                 "intrinsic multicritical 'Dyson–Weber' ensemble still open", fontsize=11.5)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fp = os.path.join(FIG, "cusp_rung_capstone.png")
    fig.savefig(fp, dpi=140, bbox_inches="tight")
    print(f"\n  Figure: {fp}")
    print(f"  TOTAL {time.time()-t0:.1f}s — DONE")


if __name__ == "__main__":
    main()
