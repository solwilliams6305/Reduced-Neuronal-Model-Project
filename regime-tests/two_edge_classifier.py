"""
two_edge_classifier.py — Direction D: a noise-only classifier separating the two excitability classes
by WHICH edge carries the noise and with WHAT scaling exponent.

  Type-I  (SNIC / QIF)          : the PHASE edge.  spike TIMING is the noisy order parameter.
        threshold rate  r ~ sigma^{2/3} ;  ISI = quartic FPT (skewed) ;  f-I CONTINUOUS (f -> 0).
  Type-II (Hopf / Stuart-Landau): the AMPLITUDE edge.  oscillation AMPLITUDE is the order parameter.
        threshold amp   A ~ sigma^{1/2} ;  f-I DISCONTINUOUS (f -> omega_H != 0).

Classifier: measure the threshold order-parameter scaling exponent.
   2/3 + timing + quartic ISI                       => SNIC / Type-I   (phase edge)
   1/2 + amplitude + (Tracy-Widom for fold-of-cycles) => Hopf / fold-of-cycles / Type-II (amplitude edge)

CAVEAT (made precise in TWO_EDGE_CLASSIFIER.md): a *generic* Hopf gives Rayleigh amplitude
fluctuations; the TRACY-WIDOM fine structure of the amplitude edge is specific to the FOLD-OF-CYCLES
canard peel-off (the burster) — the project's main result. Here we establish the robust, computable
exponent dichotomy (2/3 vs 1/2) and the onset dichotomy.
"""
from __future__ import annotations
import os
import numpy as np
from offcritical_phase_edge import J_of_nu

HERE = os.path.dirname(os.path.abspath(__file__))
OMEGA = 2.5                                          # Hopf angular frequency (nonzero onset freq)


def Phi(nu):
    nu = np.asarray(nu, float)
    if nu.ndim == 0:
        return 1.0 / J_of_nu(float(nu))
    return np.array([1.0 / J_of_nu(float(x)) for x in nu.ravel()]).reshape(nu.shape)


# ---- Type-I: noisy QIF threshold firing rate (phase edge) ---------------------------------------
def qif_rate(I, sigma, N=2500, dt=2.0e-3, fpt_factor=10.0, Tcap=150.0,
             v_th=30.0, v_reset=-30.0, seed=0):
    rng = np.random.default_rng(seed)
    v = np.full(N, float(v_reset)); tau = np.full(N, np.nan); alive = np.ones(N, bool)
    nu = I / sigma ** (4.0 / 3.0)
    Tmax = min(fpt_factor * sigma ** (-2.0 / 3.0) * J_of_nu(nu), Tcap)
    sdt = np.sqrt(dt)
    for k in range(int(Tmax / dt)):
        na = int(alive.sum())
        if na == 0:
            break
        v[alive] += (v[alive] ** 2 + I) * dt + sigma * sdt * rng.standard_normal(na)
        cross = alive & (v > v_th); tau[cross] = k * dt; alive[cross] = False
    tau = tau[np.isfinite(tau)]
    return 1.0 / tau.mean() if len(tau) > 20 else np.nan


# ---- Type-II: noisy Stuart-Landau threshold amplitude (amplitude edge) ---------------------------
def hopf_amp(mu, sigma, N=4000, dt=2.0e-3, T=50.0, omega=OMEGA, seed=0):
    rng = np.random.default_rng(seed)
    x = 0.05 * rng.standard_normal(N); y = 0.05 * rng.standard_normal(N)
    sdt = np.sqrt(dt); nstep = int(T / dt); acc = []
    for k in range(nstep):
        r2 = x * x + y * y
        dx = (mu * x - omega * y - r2 * x) * dt + sigma * sdt * rng.standard_normal(N)
        dy = (mu * y + omega * x - r2 * y) * dt + sigma * sdt * rng.standard_normal(N)
        x += dx; y += dy
        if k > int(0.5 * nstep):                      # average amplitude over the stationary tail
            acc.append(np.sqrt(x * x + y * y).mean())
    return float(np.mean(acc))


def fit_exponent(xs, ys):
    m = np.isfinite(ys) & (ys > 0)
    s, b = np.polyfit(np.log(np.asarray(xs)[m]), np.log(np.asarray(ys)[m]), 1)
    return s, b


if __name__ == "__main__":
    sig = np.array([0.18, 0.26, 0.38, 0.55, 0.80, 1.15])

    rI = np.array([qif_rate(0.0, s, seed=int(300 * s)) for s in sig])      # Type-I @ I=0
    aH = np.array([hopf_amp(0.0, s, seed=int(300 * s)) for s in sig])      # Type-II @ mu=0
    sI, bI = fit_exponent(sig, rI)
    sII, bII = fit_exponent(sig, aH)
    print(f"  Type-I  (SNIC) threshold rate  r ~ sigma^{sI:.3f}   (predict 2/3 = 0.667)")
    print(f"  Type-II (Hopf) threshold amp   A ~ sigma^{sII:.3f}   (predict 1/2 = 0.500)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    fig, ax = plt.subplots(1, 3, figsize=(15.8, 4.8))
    C1, C2 = "#dc2626", "#2563eb"

    # Panel A: the f-I onset dichotomy (Hodgkin)
    Ig = np.linspace(-0.12, 1.2, 200); s0 = 0.05
    fI = s0 ** (2 / 3) * Phi(Ig / s0 ** (4 / 3))                # continuous, noise-smoothed sqrt(I)
    ax[0].plot(Ig, fI, color=C1, lw=2.4, label="Type-I (SNIC): continuous")
    mug = np.linspace(-0.12, 1.2, 200); fH = np.where(mug >= 0, OMEGA / (2 * np.pi), 0.0)
    ax[0].plot(mug[mug < 0], fH[mug < 0], color=C2, lw=2.4)
    ax[0].plot(mug[mug >= 0], fH[mug >= 0], color=C2, lw=2.4, label="Type-II (Hopf): jump to $f_H$")
    ax[0].plot([0, 0], [0, OMEGA / (2 * np.pi)], color=C2, lw=2.4, ls=":")
    ax[0].plot(0, OMEGA / (2 * np.pi), "o", color=C2, ms=7, mec="white", mew=1.0)
    ax[0].plot(0, 0, "o", mfc="white", mec=C2, ms=7, mew=1.5)
    ax[0].axvline(0, color="#94a3b8", lw=1.0, ls=":")
    ax[0].set_xlabel(r"drive past threshold ($I$ or $\mu$)"); ax[0].set_ylabel("rhythm frequency")
    ax[0].set_title("Onset dichotomy (f–I)")
    ax[0].legend(fontsize=8.6, framealpha=0.95, edgecolor="#cbd5e1")

    # Panel B: the noise-scaling exponent dichotomy (the classifier core)
    ax[1].plot(sig, rI, "o", color=C1, ms=8, mec="white", mew=1.0,
               label=fr"Type-I rate $\sim\sigma^{{{sI:.2f}}}$")
    ax[1].plot(sig, np.exp(bI) * sig ** sI, "-", color=C1, lw=1.6)
    ax[1].plot(sig, aH, "s", color=C2, ms=8, mec="white", mew=1.0,
               label=fr"Type-II amp $\sim\sigma^{{{sII:.2f}}}$")
    ax[1].plot(sig, np.exp(bII) * sig ** sII, "-", color=C2, lw=1.6)
    sg = np.linspace(sig.min(), sig.max(), 50)
    ax[1].plot(sg, rI[2] * (sg / sig[2]) ** (2 / 3), "--", color="#7f1d1d", lw=1.3, label=r"$\sigma^{2/3}$")
    ax[1].plot(sg, aH[2] * (sg / sig[2]) ** (1 / 2), "--", color="#1e3a8a", lw=1.3, label=r"$\sigma^{1/2}$")
    ax[1].set_xscale("log"); ax[1].set_yscale("log")
    ax[1].set_xlabel(r"noise $\sigma$"); ax[1].set_ylabel("threshold order parameter")
    ax[1].set_title(r"Exponent classifier: $2/3$ vs $1/2$")
    ax[1].legend(fontsize=8.2, framealpha=0.95, edgecolor="#cbd5e1")

    # Panel C: the two-edge fingerprint summary
    ax[2].axis("off")
    rows = [
        ["", "Type-I (SNIC)", "Type-II (Hopf/FoC)"],
        ["edge", "phase / time", "amplitude / level"],
        ["order param.", "spike timing", "oscill. amplitude"],
        [r"noise exp.", r"$\sigma^{2/3}$", r"$\sigma^{1/2}$"],
        ["onset freq.", "0 (continuous)", r"$f_H\neq0$ (jump)"],
        ["fluct. law", "quartic FPT", "Rayleigh / TW*"],
    ]
    tb = ax[2].table(cellText=rows, cellLoc="center", loc="center", bbox=[0.0, 0.08, 1.0, 0.86])
    tb.auto_set_font_size(False); tb.set_fontsize(9.2)
    for j in range(3):
        tb[(0, j)].set_facecolor("#e2e8f0"); tb[(0, j)].set_text_props(weight="bold")
    for i in range(1, 6):
        tb[(i, 1)].set_text_props(color=C1); tb[(i, 2)].set_text_props(color=C2)
    ax[2].set_title("The two-edge fingerprint")
    ax[2].text(0.0, 0.0, "*TW is the fold-of-cycles canard fine structure (main result),\n"
               "  not a generic Hopf (which is Rayleigh).", fontsize=7.6, color="#475569")

    fig.tight_layout()
    out = os.path.abspath(os.path.join(HERE, "..", "figures", "two_edge_classifier.png"))
    fig.savefig(out, dpi=140); print("saved", out)
