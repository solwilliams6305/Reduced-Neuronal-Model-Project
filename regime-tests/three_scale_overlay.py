"""
three_scale_overlay.py
======================
The THREE critical-noise scales of the noisy folded cycle, on one (sigma, rate) plane.

    (floor)    TW dynamic scale      sigma_* ~ sqrt(eps2)            eta = sigma/sqrt(eps2) ~ O(1)
               -- noise-induced early escape turns on; peel-off becomes TW_beta distributed.

    (ceiling)  static-Kramers lid    sigma_crit^2 ~ DeltaU* / log(C/eps2)
               -- recrossing over the unstable cycle voids the single-escape (TW) picture.
                  Arrhenius in sigma; only logarithmically rate-dependent.

    (separate) Hopf canard strip     sigma_crit(eps) ~ eps^{3/4}     (Berglund-Gentz)
               -- a DIFFERENT passage (tracking the repelling slow manifold past the lower
                  Hopf). Rate here is the FHN time-scale eps, NOT the fold-of-cycles eps2,
                  so it is plotted in its own colour with its own guide slope.

Floor and ceiling share eps2 (fold of cycles) and bracket the Tracy-Widom-valid window;
the window opens (widens) as eps2 -> 0 because the floor falls like sqrt(eps2) while the
ceiling falls only like 1/sqrt(log(1/eps2)). The canard-strip law is a separate, steeper
object at a different rate parameter -- shown to make the "these are not the same scale"
point visible at a glance.

Kramers constants (DeltaU*, C) are re-fit live from tw_boundary.P_recross -- the verified
Arrhenius collapse (log(eps2 * N_ret) linear in 1/sigma^2, all ramp rates on one line,
R^2 ~ 0.94). Canard points are the measured table from CANARD_SIGMA_CRIT.md
(ramp 0.5 eps, (a,b) = (0.7, 0.8)).
"""
from __future__ import annotations
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tw_boundary import P_recross


# ----------------------------------------------------------------------------------
# 1. Re-fit the Kramers ceiling constants (the verified Arrhenius collapse).
# ----------------------------------------------------------------------------------
def fit_kramers(N=1000):
    sigmas = np.linspace(0.05, 0.16, 9)
    epss = [0.004, 0.002, 0.001]
    X, Yv = [], []
    for eps in epss:
        P = np.array([P_recross(sg, eps, N=N, seed=10 + j + int(1000 * eps))
                      for j, sg in enumerate(sigmas)])
        Nret = -np.log(np.clip(1 - P, 1e-6, 1))
        m = (P > 0.01) & (P < 0.7)
        X.append(1 / sigmas[m] ** 2)
        Yv.append(np.log(eps * Nret[m]))
    X = np.concatenate(X); Yv = np.concatenate(Yv)
    slope, intercept = np.polyfit(X, Yv, 1)
    R2 = 1 - np.sum((Yv - (slope * X + intercept)) ** 2) / np.sum((Yv - Yv.mean()) ** 2)
    return -slope, float(np.exp(intercept)), float(R2)


dU, C, R2 = fit_kramers()
print(f"Kramers fit:  DeltaU* = {dU:.4f}   C = {C:.4f}   R^2 = {R2:.3f}")


# ----------------------------------------------------------------------------------
# 2. The three scales.
# ----------------------------------------------------------------------------------
# fold-of-cycles rate, eps2  (floor & ceiling live here)
eps2 = np.logspace(np.log10(8e-4), np.log10(9e-3), 300)

floor = np.sqrt(eps2)                                   # sigma_* ~ sqrt(eps2),  eta = 1
# ceiling = recrossing onset, N_ret = 5% :  sigma^2 = dU / log( (C/eps2) / 0.05 )
ceiling = np.sqrt(dU / np.log((C / eps2) / 0.05))

eps2_marks = np.array([0.004, 0.002, 0.001])            # verified ramp rates
ceil_marks = np.sqrt(dU / np.log((C / eps2_marks) / 0.05))

# Hopf canard strip -- DIFFERENT passage, rate = eps (FHN).  CANARD_SIGMA_CRIT.md
eps_can = np.array([0.01, 0.02, 0.04, 0.06, 0.08])
sig_can = np.array([0.160, 0.336, 0.573, 1.007, 1.064])
xg = np.logspace(np.log10(7e-3), np.log10(0.105), 50)
Ag = sig_can[1] / eps_can[1] ** 0.75                    # anchor BG 3/4 guide at eps = 0.02
guide_can = Ag * xg ** 0.75


# ----------------------------------------------------------------------------------
# 3. Plot.
# ----------------------------------------------------------------------------------
plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155",
                     "axes.linewidth": 0.9})
fig, ax = plt.subplots(figsize=(9.2, 6.6))

C_FLOOR, C_CEIL, C_CAN, C_BAND = "#2563eb", "#dc2626", "#ea8a0b", "#16a34a"

# TW-valid band
band = ceiling > floor
ax.fill_between(eps2, floor, ceiling, where=band, color=C_BAND, alpha=0.13, lw=0,
                label="Tracy–Widom regime (clean for $\\eta\\lesssim1$)")

# floor
ax.plot(eps2, floor, color=C_FLOOR, lw=2.4,
        label=r"floor: $\sigma_*\sim\sqrt{\epsilon_2}$  (TW scale, $\eta\approx1$)")
# ceiling
ax.plot(eps2, ceiling, color=C_CEIL, lw=2.4,
        label=r"ceiling: $\sigma_{\rm crit}^2\sim\Delta U_*/\log(1/\epsilon_2)$  (Kramers)")
ax.plot(eps2_marks, ceil_marks, "o", color=C_CEIL, ms=7, mec="white", mew=1.1, zorder=5,
        label=fr"verified recrossing ($\Delta U_*={dU:.3f}$, $R^2={R2:.2f}$)")

# canard strip (separate passage)
ax.plot(xg, guide_can, "--", color=C_CAN, lw=1.6, alpha=0.85,
        label=r"$\propto\epsilon^{3/4}$  (Berglund–Gentz canard)")
ax.plot(eps_can, sig_can, "D", color=C_CAN, ms=8, mec="white", mew=1.1, zorder=5,
        label=r"Hopf canard strip $\sigma_{\rm crit}(\epsilon)$  (different passage)")

# annotations
ax.annotate("TW window widens\nas $\\epsilon_2\\to0$",
            xy=(1.0e-3, 0.047), xytext=(1.05e-3, 0.0125),
            color="#15803d", fontsize=10, ha="center",
            arrowprops=dict(arrowstyle="-|>", color="#15803d", lw=1.3))
ax.annotate("window pinches shut\n(rate too fast for TW)",
            xy=(8.0e-3, 0.088), xytext=(3.7e-3, 0.155),
            color="#7f1d1d", fontsize=9, ha="center",
            arrowprops=dict(arrowstyle="-|>", color="#7f1d1d", lw=1.1))
ax.text(0.052, 0.93, "different canard passage:\nrate $=\\epsilon$ (FHN), not $\\epsilon_2$;\n"
                     "measured slope $\\approx0.94$ (protocol-bound)",
        color="#9a5a08", fontsize=8.6, ha="center", va="top")

ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel(r"rate parameter   ($\epsilon_2$: fold of cycles    |    $\epsilon$: FHN canard strip)")
ax.set_ylabel(r"noise amplitude  $\sigma$")
ax.set_title("Three critical-noise scales of the noisy folded cycle", fontsize=13, pad=12)
ax.set_xlim(7e-4, 1.15e-1)
ax.set_ylim(1.0e-2, 1.6)
ax.grid(True, which="both", color="#e2e8f0", lw=0.6)
ax.legend(loc="lower right", fontsize=8.8, framealpha=0.95, edgecolor="#cbd5e1")
fig.tight_layout()

out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "figures", "three_scale_overlay.png")
out = os.path.abspath(out)
fig.savefig(out, dpi=140)
print("saved", out)
