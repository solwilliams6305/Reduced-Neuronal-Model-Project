"""Publication figures for W_ResurgentTransseries_paper.tex.
  Fig 1: the exact coefficient ladder v_0..v_6 (sign-coloured stems + resurgent
         late-term envelope guide) -> W_fig_ladder.pdf
  Fig 2: the Borel plane -- complex-conjugate pair (theta~50, |zeta|~1.9) + real
         s^5/10 instanton, contrasted with the deterministic root lambda_0 (+-45,
         |lambda_0|=1.258) -> W_fig_borel.pdf
No LaTeX needed; serif fonts to match the paper.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
from scipy.optimize import curve_fit
from scipy.special import gamma as Gamma

plt.rcParams.update({
    'font.family': 'serif', 'font.size': 11, 'axes.linewidth': 0.8,
    'mathtext.fontset': 'dejavuserif', 'figure.dpi': 150,
})

# ---- data ----
k   = np.arange(7)
v   = np.array([0.134, 0.111, 0.104, -0.030, -0.451, -1.19, -1.90])
verr= np.array([0.0,   0.0,   0.0,    0.0,   0.003,  0.01,  0.15])   # quoted uncertainties
BLUE, RED = '#1f5fb4', '#c0392b'

# =====================================================================
# Figure 1 -- coefficient ladder
# =====================================================================
fig, ax = plt.subplots(figsize=(6.0, 3.7))

# resurgent late-term envelope guide  v_n ~ C*Gamma(n+a)*r^{-n}*cos(n*theta - phi)
# theta, |zeta| FIXED to the reported Borel data (50 deg, 1.9); only amplitude/shape
# fitted -- the guide ILLUSTRATES the stated law, it is not an independent refit
# (a free 5-param fit on 7 points is under-determined; see paper Rem. on unidentifiability).
TH_FIX = np.deg2rad(50.0); R_FIX = 1.9
def model(n, C, a, phi):
    return C * Gamma(n + a) * R_FIX**(-n) * np.cos(n*TH_FIX - phi)
try:
    popt, _ = curve_fit(model, k, v, p0=[0.02, 1.0, 0.3],
                        bounds=([-10, 0.2, -np.pi], [10, 3.0, np.pi]), maxfev=20000)
    nn = np.linspace(0, 6.35, 400)
    ax.plot(nn, model(nn, *popt), color='0.55', lw=1.1, ls='--', zorder=1,
            label=r'late-term law $\sim\Gamma(n{+}\alpha)\,|\zeta|^{-n}\cos(n\theta-\varphi)$,'
                  '\n' r'$\theta{=}50^\circ,\,|\zeta|{=}1.9$ fixed')
except Exception as e:
    print('envelope fit skipped:', e)

# sign-coloured stems
for ki, vi, ei in zip(k, v, verr):
    c = BLUE if vi >= 0 else RED
    ax.plot([ki, ki], [0, vi], color=c, lw=2.2, zorder=2)
    ax.errorbar(ki, vi, yerr=ei, fmt='o', color=c, ms=6, capsize=3,
                elinewidth=1.2, zorder=3)
ax.axhline(0, color='k', lw=0.7)

# annotations
ax.annotate('first sign flip\n(near a node)', xy=(3, -0.03), xytext=(2.15, -0.9),
            fontsize=8.5, ha='center', color=RED,
            arrowprops=dict(arrowstyle='->', color=RED, lw=0.8))
ax.annotate('growing envelope\n(factorial divergence)', xy=(6, -1.90), xytext=(4.15, -1.55),
            fontsize=8.5, ha='center',
            arrowprops=dict(arrowstyle='->', color='0.3', lw=0.8))
ax.text(0.5, 0.55, 'sign pattern  $+,+,+,-,-,-,-$', transform=ax.transAxes,
        fontsize=9.5, ha='center',
        bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='0.7', lw=0.7))

ax.set_xlabel(r'coefficient index $k$')
ax.set_ylabel(r'$v_k$')
ax.set_title(r'Exact weak-noise coefficients  $\mathrm{Var}(Y^\star)=\eta^2\sum_k v_k\,\eta^{2k}$',
             fontsize=10.5)
ax.set_xticks(k)
ax.set_ylim(-2.35, 0.6)
ax.legend(loc='lower left', fontsize=8, framealpha=0.9)
ax.grid(True, axis='y', ls=':', lw=0.5, alpha=0.6)
fig.subplots_adjust(left=0.11, right=0.97, top=0.90, bottom=0.14)
fig.savefig('W_fig_ladder.pdf'); fig.savefig('W_fig_ladder.png', dpi=150)
print('wrote W_fig_ladder.pdf', 'envelope popt=', np.round(popt,3) if 'popt' in dir() else None)

# =====================================================================
# Figure 2 -- Borel plane
# =====================================================================
fig2, ax2 = plt.subplots(figsize=(5.4, 5.4))   # extra height for the out-of-axes legend

# stochastic conjugate pair: theta = 50 +- 2 deg, |zeta| in [1.8, 2.6], best 1.9.
# NB the modulus interval is strongly asymmetric about the best value (+0.7/-0.1), so the
# admissible REGION is the primary datum and the best-fit point is subordinate to it: the
# wedge is drawn with an edge (a stated region, not shading) and the marker is kept small,
# so the figure cannot be read as asserting precision the interval disclaims.
th0 = 50.0; zabs = 1.9
Z_LO, Z_HI = 1.8, 2.6          # |zeta| interval
TH_LO, TH_HI = 48.0, 52.0      # theta interval (deg)
def polar(rr, deg):
    a = np.deg2rad(deg); return rr*np.cos(a), rr*np.sin(a)
# admissible region, upper + lower lobe (Wedge takes increasing angles, hence two calls)
for t1, t2 in ((TH_LO, TH_HI), (-TH_HI, -TH_LO)):
    ax2.add_patch(Wedge((0, 0), Z_HI, t1, t2, width=Z_HI-Z_LO, facecolor=BLUE, alpha=0.18,
                        edgecolor=BLUE, lw=0.7, zorder=1))

# best-fit point -- deliberately subordinate to the wedge above
for sgn in (+1, -1):
    x, y = polar(zabs, sgn*th0)
    ax2.plot(x, y, '*', color=BLUE, ms=9, mec='white', mew=0.5, zorder=5,
             label=('stochastic Borel pair\n'
                    r'$|\zeta|\in[1.8,2.6]$, best $1.9$' '\n'
                    r'$\theta=50^\circ\!\pm\!2^\circ$'
                    if sgn > 0 else None))

# deterministic lambda_0 at +-45 deg, |lambda_0|=1.258 (open circles, for contrast)
labs = 1.258
for sgn in (+1, -1):
    x, y = polar(labs, sgn*45)
    ax2.plot(x, y, 'o', mfc='none', mec='0.35', mew=1.4, ms=9, zorder=4,
             label=('deterministic root\n' r'$\lambda_0=1.258\,e^{\pm i\,45^\circ}$'
                    if sgn > 0 else None))

# real s^5/10 instanton on positive real axis (competes at ~same modulus)
ax2.plot(1.9, 0.0, 's', color=RED, ms=9, zorder=5,
         label='real instanton\n' r'$S=s^5/10$')

# radial guide lines at 45 and 50 deg to show the phase gap
for deg, c, ls in ((45, '0.35', '--'), (50, BLUE, '-')):
    x, y = polar(2.75, deg)
    ax2.plot([0, x], [0, y], color=c, lw=0.9, ls=ls, alpha=0.7, zorder=2)
ax2.annotate(r'$\Delta\theta\approx5^\circ$', xy=polar(1.55, 47.5), fontsize=9,
             color=BLUE, ha='left')

# unit circle + axes
tt = np.linspace(0, 2*np.pi, 300)
ax2.plot(np.cos(tt), np.sin(tt), color='0.8', lw=0.7, zorder=0)
ax2.axhline(0, color='k', lw=0.6); ax2.axvline(0, color='k', lw=0.6)

ax2.set_xlim(-0.6, 2.9); ax2.set_ylim(-2.4, 2.4)
ax2.set_aspect('equal')
ax2.set_xlabel(r'$\mathrm{Re}\,\zeta$'); ax2.set_ylabel(r'$\mathrm{Im}\,\zeta$')
ax2.set_title(r'Borel plane of $\mathcal{W}$: conjugate pair $+$ real instanton', fontsize=10.5)
# Legend sits OUTSIDE, below the axes. The conjugate pair occupies both right-hand
# quadrants and the guide rays cross the upper-left, so every in-axes corner occludes
# something -- in particular a lower-right legend hides the lower lobe and makes the
# pair read as asymmetric.
ax2.legend(loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=3, fontsize=7.4,
           framealpha=0.92, handlelength=1.2, columnspacing=1.4, borderpad=0.5)
ax2.grid(True, ls=':', lw=0.4, alpha=0.5)
fig2.subplots_adjust(left=0.13, right=0.97, top=0.93, bottom=0.23)
fig2.savefig('W_fig_borel.pdf'); fig2.savefig('W_fig_borel.png', dpi=150)
print('wrote W_fig_borel.pdf')
