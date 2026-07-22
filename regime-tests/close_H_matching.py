"""
close_H_matching.py
===================
Direction A, part (a): the outer-matching integration that actually closes the
integrated hazard H with the inner (sub-exponential) prefactor.

H is the integrated Kramers hazard for crossing the separatrix before the fold,

        H = int_0^inf  lambda(Y) dY ,   lambda(Y, eta) = eta^{2/3} J(Y/eta^{4/3}),

so in the rescaled level y = Y/eta^{4/3},   H = eta^2 * I+,   I+ = int_0^inf J(y) dy.
Quasi-statically  lambda_qs(Y) = (sqrt(Y)/pi) e^{-8 Y^{3/2}/3 eta^2}  gives  H_qs = eta^2/(4 pi).

THE MATCHING (resolving the C3 boundary-layer worry).  The naive worry was that
int J(y) dy diverges -- true, but only on the y -> -inf side, where J(y) ~ sqrt|y|/pi.
That branch is the CERTAIN late escape PAST the fold (Y<0); it is excluded by the
early-escape definition (Y>0).  On the early-escape side Y>0 the inner current J(y) is
uniformly valid: it MATCHES the outer quasi-static rate as y -> inf (verified earlier,
J/J_qs -> 0.98 at y=3) and corrects it near the fold.  Hence int_0^inf J(y) dy is finite
and is the uniformly-correct (composite) integrated hazard -- no separate outer piece to add.

Output: the prefactor  K = I+ / (1/4pi)  on H = eta^2/4pi, and the closed H = K eta^2/4pi.
"""
from __future__ import annotations
import os
import numpy as np
from frisch_lloyd_current import J_quadrature, J_kramers

# ---- integrate over the early-escape region y > 0 ----
y = np.linspace(1e-3, 5.0, 800)
J = np.array([J_quadrature(v) for v in y])
Jqs = J_kramers(y)

I_plus = float(np.trapezoid(J, y))
I_qs = float(np.trapezoid(Jqs, y))
I_qs_exact = 1.0 / (4.0 * np.pi)
K = I_plus / I_qs_exact

print(f"  I+   = ∫_0^∞ J(y) dy        = {I_plus:.5f}")
print(f"  I_qs = ∫_0^∞ J_qs(y) dy     = {I_qs:.5f}   (analytic 1/4π = {I_qs_exact:.5f}; check)")
print(f"  prefactor  K = I+ /(1/4π)   = {K:.4f}")
print(f"  => closed  H = K·η²/4π = {K:.3f}·η²/4π = {I_plus:.5f}·η²"
      f"   (quasi-static H_qs = {I_qs_exact:.5f}·η²)")
print(f"     i.e. the sub-exponential prefactor RAISES the integrated hazard by "
      f"{(K-1)*100:.1f}% relative to the quasi-static η²/4π "
      f"(finite rate near the fold outweighs the mid-range deficit).")

# cumulative integrals (to show convergence to I+ and I_qs)
cumJ = np.concatenate(([0.0], np.cumsum(0.5 * (J[1:] + J[:-1]) * np.diff(y))))
cumQ = np.concatenate(([0.0], np.cumsum(0.5 * (Jqs[1:] + Jqs[:-1]) * np.diff(y))))

# ---- figure ----
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
C_J, C_Q, C_EX, C_DEF = "#2563eb", "#64748b", "#16a34a", "#dc2626"

fig, ax = plt.subplots(1, 2, figsize=(12.2, 5.0))

m = y <= 4.0
ax[0].plot(y[m], J[m], color=C_J, lw=2.4, label=r"inner current $\mathcal{J}(y)$")
ax[0].plot(y[m], Jqs[m], "--", color=C_Q, lw=1.9,
           label=r"quasi-static $\frac{\sqrt{y}}{\pi}e^{-8y^{3/2}/3}$")
exc = J >= Jqs
ax[0].fill_between(y[m], Jqs[m], J[m], where=exc[m], color=C_EX, alpha=0.22,
                   label="excess (finite rate near fold)")
ax[0].fill_between(y[m], Jqs[m], J[m], where=~exc[m], color=C_DEF, alpha=0.18,
                   label="deficit (Kramers over-counts)")
ax[0].set_xlabel(r"rescaled level  $y = Y/\eta^{4/3}$")
ax[0].set_ylabel(r"escape current")
ax[0].set_title("Integrand of the early-escape hazard")
ax[0].set_xlim(0, 4); ax[0].set_ylim(0, 0.17)
ax[0].grid(True, color="#eef2f7", lw=0.6)
ax[0].legend(loc="upper right", fontsize=8.6, framealpha=0.95, edgecolor="#cbd5e1")

ax[1].plot(y[m], cumJ[m], color=C_J, lw=2.4, label=r"$\int_0^y \mathcal{J}$")
ax[1].plot(y[m], cumQ[m], "--", color=C_Q, lw=1.9, label=r"$\int_0^y \mathcal{J}_{\rm qs}$")
ax[1].axhline(I_plus, color=C_J, lw=1.0, ls=":")
ax[1].axhline(I_qs_exact, color=C_Q, lw=1.0, ls=":")
ax[1].text(2.0, I_plus + 0.002, fr"$I_+={I_plus:.4f}$", color=C_J, fontsize=9)
ax[1].text(2.0, I_qs_exact - 0.006, fr"$1/4\pi={I_qs_exact:.4f}$", color="#475569", fontsize=9)
ax[1].annotate(fr"prefactor $K=\dfrac{{I_+}}{{1/4\pi}}={K:.3f}$",
               xy=(3.4, I_plus), xytext=(1.4, 0.045), fontsize=11, color="#1e3a8a",
               arrowprops=dict(arrowstyle="-|>", color="#1e3a8a", lw=1.1))
ax[1].set_xlabel(r"rescaled level  $y$")
ax[1].set_ylabel(r"cumulative hazard / $\eta^2$")
ax[1].set_title(r"Closed integrated hazard:  $H = K\,\eta^2/4\pi$")
ax[1].set_xlim(0, 4); ax[1].set_ylim(0, 0.09)
ax[1].grid(True, color="#eef2f7", lw=0.6)
ax[1].legend(loc="lower right", fontsize=9, framealpha=0.95, edgecolor="#cbd5e1")

fig.tight_layout()
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures", "close_H_matching.png")
out = os.path.abspath(out)
fig.savefig(out, dpi=140)
print("\nsaved", out)
