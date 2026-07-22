"""
cusp_law_fingerprint.py — the full fingerprint of the cusp (Weber-TW) law, to discipline conjectures.
=====================================================================================================

Beyond skew/kurtosis + tails: extract the standardized cumulants κ3…κ6, the β-family behaviour, and the
quantile map to Tracy–Widom. Together these strongly constrain what the combined/exact law can be.

Object: first-explosion law of dp=(sign(Y)|Y|^q − p²)dY + η dW (η²=4/β). q=1 = TW (fold), q=2 = cusp.

Probes:
  (1) CUMULANTS κ3,κ4,κ5,κ6 (standardized) for q=1,2,3 — the shape fingerprint.
  (2) β-FAMILY: cusp (q=2) skew/exkurt at β=1,2,4 — does it scale like a TW_β-type family?
  (3) QQ vs TW: is the cusp law a simple deterministic transform of TW? Q_cusp(p) vs Q_TW(p).
  (4) Candidate forms ruled in/out by the above.

Output: figures/cusp_law_fingerprint.png + the fingerprint table.
"""
from __future__ import annotations
import os, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
TW2 = dict(sk=0.2241, ku=0.0934)


def riccati_escape(q, eta, Y0=3.0, Yend=-7.0, dtau=2.0e-3, M=45000, p_expl=-10.0, seed=0):
    rng = np.random.default_rng(seed); sdt = np.sqrt(dtau)
    p = np.full(M, Y0**(q/2.0)); Yesc = np.full(M, np.nan); done = np.zeros(M, bool)
    n = int((Y0 - Yend)/dtau)
    for i in range(n):
        Y = Y0 - i*dtau
        p = p + (np.sign(Y)*abs(Y)**q - p**2)*dtau - eta*sdt*rng.standard_normal(M)
        p = np.clip(p, p_expl - 1.0, 25.0)
        cr = (~done) & (p < p_expl); Yesc[cr] = Y; done |= cr
    return Yesc[np.isfinite(Yesc)]


def std_cumulants(x):
    z = (x - x.mean())/x.std()
    m3 = np.mean(z**3); m4 = np.mean(z**4); m5 = np.mean(z**5); m6 = np.mean(z**6)
    k3 = m3; k4 = m4 - 3.0; k5 = m5 - 10*m3; k6 = m6 - 15*m4 - 10*m3**2 + 30.0
    return k3, k4, k5, k6


def main():
    t0 = time.time()
    print("=" * 80)
    print("Cusp (Weber-TW) law fingerprint — cumulants, β-family, QQ vs TW")
    print("=" * 80)

    # (1) cumulants along the low ladder, β=2
    eta2 = np.sqrt(2.0)
    laws = {q: riccati_escape(q, eta2, seed=3+int(2*q)) for q in (1.0, 2.0, 3.0)}
    print(f"\n  (1) standardized cumulants (β=2):")
    print(f"      {'q':>4} {'κ3(skew)':>9} {'κ4(exk)':>9} {'κ5':>8} {'κ6':>8}  {'note':>14}")
    cum = {}
    for q in (1.0, 2.0, 3.0):
        k3, k4, k5, k6 = std_cumulants(laws[q]); cum[q] = (k3, k4, k5, k6)
        nt = "fold = TW" if q == 1 else ("cusp = Weber-TW" if q == 2 else "swallowtail")
        print(f"      {q:4.1f} {k3:+9.3f} {k4:+9.3f} {k5:+8.2f} {k6:+8.2f}  {nt:>14}")
    print(f"      TW₂ ref: κ3 {TW2['sk']:+.3f}, κ4 {TW2['ku']:+.3f} (κ5,κ6 small +).")

    # (2) β-family for the cusp
    print(f"\n  (2) cusp (q=2) β-family   (η=2/√β):")
    print(f"      {'β':>4} {'skew':>8} {'exkurt':>8}")
    betafam = {}
    for beta in (1.0, 2.0, 4.0):
        eta = 2.0/np.sqrt(beta)
        x = riccati_escape(2.0, eta, seed=50+int(beta))
        k3, k4, _, _ = std_cumulants(x); betafam[beta] = (k3, k4)
        print(f"      {beta:4.1f} {k3:+8.3f} {k4:+8.3f}")
    print(f"      (TW_β skew: β=1 →0.29, β=2 →0.22, β=4 →0.17 — decreasing. Cusp pattern above.)")

    # (3) QQ vs TW
    cusp = laws[2.0]; fold = laws[1.0]
    zc = np.sort((cusp - cusp.mean())/cusp.std()); zf = np.sort((fold - fold.mean())/fold.std())
    ps = np.array([0.5, 1, 2, 5, 10, 25, 50, 75, 90, 95, 98, 99, 99.5])
    Qc = np.percentile(zc, ps); Qf = np.percentile(zf, ps)
    print(f"\n  (3) QQ — standardized quantiles, cusp vs TW(fold):")
    print(f"      {'p%':>6} {'TW(q=1)':>9} {'cusp(q=2)':>10} {'Δ=cusp−TW':>11}")
    for i, p in enumerate(ps):
        print(f"      {p:6.1f} {Qf[i]:+9.3f} {Qc[i]:+10.3f} {Qc[i]-Qf[i]:+11.3f}")
    # is cusp ≈ affine transform of TW?  (it's not, if Δ is non-constant/non-linear)
    A = np.polyfit(Qf, Qc, 1); resid_aff = np.std(Qc - (A[0]*Qf + A[1]))
    B = np.polyfit(Qf, Qc, 3); resid_cub = np.std(Qc - np.polyval(B, Qf))
    print(f"      affine fit Q_cusp≈{A[0]:.2f}·Q_TW+{A[1]:+.2f}: resid {resid_aff:.3f};  "
          f"cubic-in-Q_TW resid {resid_cub:.3f}")
    print(f"      cubic coeffs (Q_TW³,²,¹,⁰): {B[0]:+.3f} {B[1]:+.3f} {B[2]:+.3f} {B[3]:+.3f}")

    print(f"\n  READOUT:")
    print(f"   • cusp is MORE skewed than TW (+{cum[2.0][0]:.2f} vs +{cum[1.0][0]:.2f}) AND has NEGATIVE κ4")
    print(f"     ({cum[2.0][1]:+.2f}) — so it is NOT TW smeared by Gaussian (that lowers skew, keeps κ4>0).")
    print(f"   • the QQ map cusp↔TW is a smooth monotone reparametrization (cubic resid {resid_cub:.3f}):")
    print(f"     the cusp law ≈ a fixed nonlinear transform of TW — candidate for a closed 'φ(TW)' form.")

    # ---- figure ----
    fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.7))
    # (A) cumulants along the ladder
    qs = [1.0, 2.0, 3.0]
    ax[0].axhline(0, color="grey", lw=0.8)
    for idx, lab, col in [(0, "κ3 skew", "#b3402b"), (1, "κ4 exkurt", "#7a3b8f"),
                          (2, "κ5", "#2c7d59"), (3, "κ6", "#1f3b73")]:
        ax[0].plot(qs, [cum[q][idx] for q in qs], "o-", color=col, lw=1.6, label=lab)
    ax[0].set_xlabel("turning order q"); ax[0].set_ylabel("standardized cumulant")
    ax[0].set_title("(A) cumulant fingerprint along the ladder"); ax[0].legend(fontsize=8, frameon=False)

    # (B) β-family
    bb = [1.0, 2.0, 4.0]
    ax[1].plot(bb, [betafam[b][0] for b in bb], "o-", color="#7a3b8f", lw=1.8, label="cusp skew")
    ax[1].plot(bb, [betafam[b][1] for b in bb], "s-", color="#2c7d59", lw=1.8, label="cusp exkurt")
    ax[1].plot(bb, [0.29, 0.22, 0.17], "^--", color="#b3402b", lw=1.3, label="TW_β skew (ref)")
    ax[1].axhline(0, color="grey", lw=0.8)
    ax[1].set_xlabel("β = 4/η²"); ax[1].set_ylabel("moment")
    ax[1].set_title("(B) cusp β-family"); ax[1].legend(fontsize=8, frameon=False)

    # (C) QQ cusp vs TW
    pp = np.linspace(0.5, 99.5, 80)
    qc = np.percentile(zc, pp); qf = np.percentile(zf, pp)
    ax[2].plot([-4, 4], [-4, 4], "k:", lw=1, label="identity (cusp=TW)")
    ax[2].plot(qf, qc, "-", color="#7a3b8f", lw=2.2, label="Q_cusp vs Q_TW")
    ax[2].plot(qf, np.polyval(B, qf), "--", color="#b3402b", lw=1.2, label="cubic φ(TW)")
    ax[2].set_xlabel("TW (fold) quantile"); ax[2].set_ylabel("cusp quantile")
    ax[2].set_title("(C) cusp ≈ a fixed nonlinear transform of TW"); ax[2].legend(fontsize=8.5, frameon=False)

    fig.suptitle("Cusp (Weber-TW) fingerprint: cumulants, β-family, and a smooth quantile map to TW — "
                 "the law is 'beyond TW' (more skew, κ4<0), ≈ a fixed nonlinear transform of TW", fontsize=10.2)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fp = os.path.join(FIG, "cusp_law_fingerprint.png")
    fig.savefig(fp, dpi=140, bbox_inches="tight")
    print(f"\n  Figure: {fp}")
    print(f"  TOTAL {time.time()-t0:.1f}s — DONE")


if __name__ == "__main__":
    main()
