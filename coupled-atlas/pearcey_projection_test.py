"""
pearcey_projection_test.py — §2 of RH_DIRECTION_NOVEL_ANGLES: is the cusp RH "Pearcey, reduced"?
================================================================================================

Claim (§2 route A): the cusp edge law is the slow-flow reduction of the PEARCEY process — integrate out
the third/y-direction (the extra cuspoid unfolding parameter) with the KP slow-flow weighting, and the
reduced law should be the cusp (Weber) law. If so, PIV descends from the known Pearcey parent.

CONSTRUCTION STATUS (honest):
  • [CITED structure] the Pearcey object = the A₃ diffraction-catastrophe ODE p'''=y p' + i x p, i.e. the
    swept 3rd-order inner equation u'''=(Y−ηξ)u + y·u' (deformation y = the second cuspoid parameter).
    We use the swept-inner-equation form (the catastrophe-ladder Pearcey, peeloff_catastrophe_ladder),
    NOT the cited RMT Pearcey kernel — [RECONSTRUCTED proxy], honestly flagged.
  • [RECONSTRUCTED] the slow-flow weighting w(y): the y-direction is integrated out (mixture of peel-off
    laws L_y over the passage). We scan weightings (Gaussian fast-mode invariant; uniform window) and
    the y-limits, and ask whether ANY natural reduction lands on the cusp fingerprint.

Cusp fingerprint to match: skew +0.61, exkurt −0.24, κ₅≈−2.2, κ₆≈−2.7, left-tail exponent 5
(e^{−|s|⁵/20}), cubic QQ vs TW  Q≈Q_TW+0.072Q²−0.035Q³.
   PASS = reduced-Pearcey reproduces the cusp fingerprint; FAIL = it does not.

Output: figures/pearcey_projection_test.png + pass/fail summary.
"""
from __future__ import annotations
import os, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from peeloff_catastrophe_ladder import ladder_peeloff       # n=2 fold (TW), n=3 Pearcey reference

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
CUSP = dict(skew=0.61, exk=-0.24, k5=-2.2, k6=-2.7, aL=5.0, qq2=0.072, qq3=-0.035)


def deformed_pearcey(yd, eta, N=7000, Y0=5.0, dt=1.0e-3, Ymin=-7.0, seed=0):
    """First node of the deformed swept Pearcey u'''=(Y−ηξ)u + yd·u' (Heun); yd=0 is the cusp/Pearcey."""
    rng = np.random.default_rng(seed); sdt = np.sqrt(dt)
    U = [np.full(N, Y0**(j/3.0)) for j in range(3)]
    Y = Y0; Yz = np.full(N, np.nan); done = np.zeros(N, bool)
    nstep = int((Y0 - Ymin)/dt)
    for _ in range(nstep):
        al = ~done
        if not al.any():
            break
        dB = sdt*rng.standard_normal(N); Yp = Y - dt
        U1 = [None]*3
        U1[0] = U[0] + U[1]*dt
        U1[1] = U[1] + U[2]*dt
        U1[2] = U[2] + (Y*U[0] + yd*U[1])*dt - eta*U[0]*dB
        Un = [None]*3
        Un[0] = U[0] + 0.5*(U[1] + U1[1])*dt
        Un[1] = U[1] + 0.5*(U[2] + U1[2])*dt
        Un[2] = U[2] + 0.5*((Y*U[0] + yd*U[1]) + (Yp*U1[0] + yd*U1[1]))*dt - eta*0.5*(U[0] + U1[0])*dB
        U = Un; Y = Yp
        cr = al & (U[0] < 0.0); Yz[cr] = Y; done |= cr
    return Yz[np.isfinite(Yz)]


def std_cumulants(x):
    z = (x - x.mean())/x.std()
    m3, m4, m5, m6 = [np.mean(z**k) for k in (3, 4, 5, 6)]
    return m3, m4 - 3, m5 - 10*m3, m6 - 15*m4 - 10*m3**2 + 30


def left_tail_alpha(x):
    xs = np.sort(x); n = xs.size; F = (np.arange(n) + 0.5)/n
    m = (F < 0.06) & (-np.log(F) > 1.0)
    return np.polyfit(np.log(np.abs(xs[m])), np.log(-np.log(F[m])), 1)[0] if m.sum() > 20 else np.nan


def cubic_qq(x, ref):
    """QQ-map of standardized x vs standardized ref, fit Q_x = Q_ref + a Q_ref² + b Q_ref³."""
    pp = np.linspace(2, 98, 49)
    qx = np.percentile((x - x.mean())/x.std(), pp); qr = np.percentile((ref - ref.mean())/ref.std(), pp)
    A = np.linalg.lstsq(np.column_stack([qr, qr**2, qr**3]), qx, rcond=None)[0]
    return A[1], A[2]   # quadratic, cubic coefficients


def main():
    t0 = time.time()
    eta = np.sqrt(2.0)
    print("=" * 80)
    print("§2 Pearcey-projection test — is the cusp law 'Pearcey, reduced'?")
    print("=" * 80)
    tw = ladder_peeloff(2, eta, N=8000, seed=2)          # TW reference (fold)

    yds = [-6, -4, -3, -2, -1, 0, 1, 2, 3, 4, 6]
    print(f"\n  deformed Pearcey peel-off vs deformation y (y=0 = symmetric Pearcey):")
    print(f"    {'y':>5} {'skew':>7} {'exkurt':>8} {'αL':>6}")
    samp = {}; sk = {}; ku = {}
    for yd in yds:
        x = deformed_pearcey(float(yd), eta, seed=30 + yd)
        samp[yd] = x; s, k, _, _ = std_cumulants(x); sk[yd] = s; ku[yd] = k
        print(f"    {yd:5d} {s:+7.3f} {k:+8.3f} {left_tail_alpha(x):6.2f}")
    print(f"  (Pearcey y=0: skew {sk[0]:+.2f}, exk {ku[0]:+.2f}.  Cusp target: skew +0.61, exk −0.24, αL 5.)")

    # slow-flow projection: mixture of L_y over y with weighting w(y); scan reconstructed weightings
    print(f"\n  slow-flow y-projection (mixture ∫L_y w(y)dy; w reconstructed):")
    print(f"    {'weighting':>26} {'skew':>7} {'exkurt':>8} {'κ5':>7} {'κ6':>7} {'αL':>6} {'QQ(a,b)':>14}")
    def mixture(weights):
        parts = [np.random.default_rng(7).choice(samp[yd], size=int(2000*w)) for yd, w in weights.items() if w > 0]
        return np.concatenate(parts)
    results = {}
    weightings = {
        "uniform |y|≤2": {yd: (1.0 if abs(yd) <= 2 else 0.0) for yd in yds},
        "uniform all y": {yd: 1.0 for yd in yds},
        "Gaussian σ_y=2": {yd: np.exp(-yd**2/(2*2.0**2)) for yd in yds},
        "Gaussian σ_y=4": {yd: np.exp(-yd**2/(2*4.0**2)) for yd in yds},
        "large |y| (≥3)": {yd: (1.0 if abs(yd) >= 3 else 0.0) for yd in yds},
    }
    for name, w in weightings.items():
        x = mixture(w)
        s, k, k5, k6 = std_cumulants(x); aL = left_tail_alpha(x); a2, a3 = cubic_qq(x, tw)
        results[name] = (s, k, k5, k6, aL, a2, a3, x)
        print(f"    {name:>26} {s:+7.3f} {k:+8.3f} {k5:+7.2f} {k6:+7.2f} {aL:6.2f} ({a2:+.2f},{a3:+.2f})")

    # also the best single-y (closest skew to cusp) for reference
    yd_best = min(yds, key=lambda y: abs(sk[y] - CUSP['skew']))
    print(f"  best single-y for skew: y={yd_best} (skew {sk[yd_best]:+.2f}, exk {ku[yd_best]:+.2f})")

    # verdict: does ANY reduction match the cusp fingerprint?
    def matches(s, k, k5, k6, aL, a2, a3):
        return (abs(s-CUSP['skew']) < 0.12 and abs(k-CUSP['exk']) < 0.15 and
                abs(aL-CUSP['aL']) < 1.2 and abs(a2-CUSP['qq2']) < 0.05)
    any_match = any(matches(*results[n][:7]) for n in results)
    print(f"\n  CUSP fingerprint: skew +0.61, exk −0.24, κ5 −2.2, κ6 −2.7, αL 5, QQ(+0.07,−0.04).")
    print(f"  ⇒ {'PASS' if any_match else 'FAIL'}: a slow-flow reduction of Pearcey "
          f"{'reproduces' if any_match else 'does NOT reproduce'} the cusp fingerprint.")
    if not any_match:
        print(f"    The deformed Pearcey peel-off stays near the Pearcey class (skew {sk[0]:+.2f}→~{max(sk.values()):+.2f},")
        print(f"    exkurt {ku[0]:+.2f}) and does not reach the cusp/Weber (skew +0.61, exk −0.24, αL 5) under any")
        print(f"    reconstructed y-weighting. ⇒ the cusp RH is NOT simply 'Pearcey reduced' by this projection.")
    print(f"  [CITED-structure/RECONSTRUCTED-proxy] Pearcey = catastrophe-ladder 3rd-order ODE (not the RMT")
    print(f"    Pearcey kernel); [RECONSTRUCTED] slow-flow weighting; [NUMERIC] the peel-off laws & mixtures.")

    # ---- figure ----
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    yy = np.array(yds)
    ax[0].axhline(CUSP['skew'], color="#1f9e75", ls="--", lw=1.2, label="cusp skew +0.61")
    ax[0].axhline(0.224, color="#b3402b", ls=":", lw=1, label="TW 0.22")
    ax[0].plot(yy, [sk[y] for y in yds], "o-", color="#7a3b8f", lw=1.7, label="skew(y)")
    ax[0].plot(yy, [ku[y] for y in yds], "s-", color="#2c7d59", lw=1.5, label="exkurt(y)")
    ax[0].axhline(0, color="grey", lw=0.6)
    ax[0].set_xlabel("Pearcey deformation y"); ax[0].set_ylabel("moment")
    ax[0].set_title("(A) deformed-Pearcey peel-off vs y"); ax[0].legend(fontsize=8, frameon=False)

    gg = np.linspace(-4, 4, 200)
    ax[1].plot(gg, np.exp(-gg**2/2)/np.sqrt(2*np.pi), "k:", lw=1, label="Gaussian")
    for name, col in [("uniform all y", "#7a3b8f"), ("large |y| (≥3)", "#1f3b73")]:
        x = results[name][7]; z = (x - x.mean())/x.std()
        ax[1].hist(z, bins=60, range=(-4, 4), density=True, histtype="step", lw=1.7, color=col,
                   label=f"{name} (skew {results[name][0]:+.2f})")
    zc = (samp[0] - samp[0].mean())/samp[0].std()
    ax[1].hist(zc, bins=60, range=(-4, 4), density=True, histtype="step", lw=1.5, color="#888",
               label=f"Pearcey y=0 (skew {sk[0]:+.2f})")
    ax[1].set_xlabel("standardized peel-off"); ax[1].set_ylabel("density")
    ax[1].set_title("(B) reduced-Pearcey laws"); ax[1].legend(fontsize=7.5, frameon=False)

    names = list(results.keys())
    sks = [results[n][0] for n in names]; kus = [results[n][1] for n in names]
    ax[2].scatter(sks, kus, s=70, color="#7a3b8f", zorder=4)
    for n, s, k in zip(names, sks, kus):
        ax[2].annotate(n.split()[0], (s, k), textcoords="offset points", xytext=(4, 3), fontsize=7)
    ax[2].scatter([CUSP['skew']], [CUSP['exk']], marker="*", s=180, color="#1f9e75", zorder=5, label="cusp target")
    ax[2].scatter([0.224], [0.093], marker="D", s=70, color="#b3402b", zorder=5, label="TW")
    ax[2].axhline(0, color="grey", lw=0.6)
    ax[2].set_xlabel("skew"); ax[2].set_ylabel("excess kurtosis")
    ax[2].set_title("(C) reductions vs cusp target"); ax[2].legend(fontsize=8, frameon=False)

    fig.suptitle("§2 Pearcey-projection: does the slow-flow-reduced Pearcey reproduce the cusp law?", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fp = os.path.join(FIG, "pearcey_projection_test.png")
    fig.savefig(fp, dpi=140, bbox_inches="tight")
    print(f"\n  Figure: {fp}")
    print(f"  TOTAL {time.time()-t0:.1f}s — DONE")


if __name__ == "__main__":
    main()
