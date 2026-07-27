"""PILOT STORYBOARD -- module gamma2+gamma3: the Borel plane and the ambiguity identity.

STILLS ONLY.  Per exposition/PROBE_PROTOCOL.md the explanation is tested as a sequence of
static panels BEFORE anything is animated: if the stills do not carry it, animation will only
make the failure feel smoother.

The arc, in one line:  the terms shrink then grow -> the smallest term is ~e^{-A/x} -> Borel
explains why (a singularity at A) -> the contour ambiguity across its cut is the SAME e^{-A/x}
-> therefore the answer needs an extra exponentially small term.

Panels 1-5 are computed from the real cusp coefficients v_0..v_6.
Panels 6-8 are SCHEMATIC (labelled as such on the panel): a faithful lateral Borel sum from 7
starved coefficients would be a research computation, not a storyboard, and would test the
implementation rather than the explanation.

Run:  python3 storyboard.py   ->  storyboard_contact.png  (+ panel_N.pdf if SPLIT=True)
"""
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# (no patch primitives needed: the conjugate-pair wedges were removed from panel 6)

plt.rcParams.update({
    'font.family': 'serif', 'font.size': 9, 'axes.linewidth': 0.8,
    'mathtext.fontset': 'dejavuserif', 'figure.dpi': 140,
})

# ---------------------------------------------------------------- data (real)
# cusp weak-noise variance coefficients, as used in W_make_figures.py
v = np.array([0.134, 0.111, 0.104, -0.030, -0.451, -1.19, -1.90])
k = np.arange(len(v))
A = 1.9                      # real instanton, S = s^5/10 -- the singularity ON the contour
TH_PAIR = 50.0               # conjugate pair phase (deg)
Z_LO, Z_HI = 1.8, 2.6        # |zeta| interval for the pair
BLUE, RED, GREY = '#1f5fb4', '#c0392b', '0.45'
SCHEM = dict(fontsize=7, style='italic', color=RED)
# white backing so annotations stay legible where they cross a curve
BG = dict(bbox=dict(facecolor='white', edgecolor='none', alpha=0.88, pad=1.6))

fig, axes = plt.subplots(2, 4, figsize=(15.5, 7.4))
ax = axes.ravel()

# ---------------------------------------------------------------- 1. the numbers
a = ax[0]
for ki, vi in zip(k, v):
    c = BLUE if vi >= 0 else RED
    a.plot([ki, ki], [0, vi], color=c, lw=2.4)
    a.plot(ki, vi, 'o', color=c, ms=5)
a.axhline(0, color='k', lw=0.7)
a.annotate('sign flips here', xy=(3, -0.03), xytext=(1.5, -1.15), fontsize=7.5, color=RED,
           arrowprops=dict(arrowstyle='->', color=RED, lw=0.8))
a.set_xlabel(r'$k$'); a.set_ylabel(r'$v_k$'); a.set_xticks(k)
a.set_title('1.  Seven numbers', loc='left', fontsize=10)
a.text(0.03, 0.06, 'they grow, and the sign\nturns over once',
       transform=a.transAxes, fontsize=7.5, color=GREY)

# ---------------------------------------------------------------- 2. a MODEL series
# The real ladder cannot show optimal truncation (see panel 4), so the mechanism is taught
# first on the canonical Borel-summable prototype  m_k = k!/A^{k+1}, labelled as a model.
a = ax[1]
km = np.arange(15)
mk = np.array([math.factorial(int(i)) for i in km], dtype=float) / A**(km + 1)
xs_demo = [0.22, 0.35, 0.55]
for xv, col in zip(xs_demo, [BLUE, '#7a5cc0', RED]):
    t = mk * xv**km
    a.semilogy(km, t, 'o-', color=col, ms=3.5, lw=1.1, label=rf'$x={xv}$')
    kmin = int(np.argmin(t))
    a.plot(kmin, t[kmin], 'v', color=col, ms=9, mfc='none', mew=1.6)
a.set_xlabel(r'$k$'); a.set_ylabel(r'$|m_k|\,x^k$')
a.legend(fontsize=7, loc='upper center')
a.set_title('2.  Model series: a floor that moves', loc='left', fontsize=10)
a.text(0.03, 0.05, 'triangles = smallest term.\nit slides left as $x$ grows',
       transform=a.transAxes, fontsize=7.5, color=GREY)
a.text(0.60, 0.05, r'model  $m_k=k!/A^{k+1}$', transform=a.transAxes, **SCHEM, **BG)

# ---------------------------------------------------------------- 3. least term ~ exp(-A/x)
a = ax[2]
xs = np.linspace(0.18, 0.75, 220)
least_m = np.array([np.min(mk * xx**km) for xx in xs])
a.semilogy(xs, least_m, color=BLUE, lw=1.9, label='smallest term (model)')
a.semilogy(xs, np.exp(-A / xs), color=RED, ls='--', lw=1.5, label=rf'$e^{{-A/x}}$,  $A={A}$')
a.set_xlabel(r'$x=\eta^2$'); a.set_ylabel('size')
a.legend(fontsize=7, loc='lower right')
a.set_title(r'3.  The floor is $e^{-A/x}$', loc='left', fontsize=10)
a.text(0.04, 0.95, 'an exponential falls out\nof a pure power series.\n'
                   r'$A$ is the only input.',
       transform=a.transAxes, fontsize=7.5, color=GREY, **BG, va='top')
a.text(0.04, 0.50, 'curves differ by an\nalgebraic prefactor only',
       transform=a.transAxes, **SCHEM, **BG, va='top')

# ---------------------------------------------------------------- 4. the REAL ladder is worse
a = ax[3]
for xv, col in zip([0.25, 0.40, 0.55], [BLUE, '#7a5cc0', RED]):
    t = np.abs(v) * xv**k
    a.semilogy(k, t, 'o-', color=col, ms=4, lw=1.2, label=rf'$x={xv}$')
a.axvline(3, color=RED, lw=1.0, ls=':')
a.set_xlabel(r'$k$'); a.set_ylabel(r'$|v_k|\,x^k$'); a.set_xticks(k)
a.legend(fontsize=7, loc='lower left')
a.set_title('4.  Our seven do NOT show it', loc='left', fontsize=10)
a.text(0.06, 0.96, 'the dip is pinned at $k=3$ for every $x$ —\n'
                   "that's the sign-flip node, not a moving floor.\n"
                   'seven terms is too few to read $A$ off directly.',
       transform=a.transAxes, fontsize=7.2, color=RED, **BG, va='top')
a.text(0.06, 0.60, r'$\Rightarrow$ we need a different route to $A$',
       transform=a.transAxes, fontsize=8.5, **BG, va='top')

# ---------------------------------------------------------------- 5. divide by k!
a = ax[4]
b = v / np.array([math.factorial(int(i)) for i in k], dtype=float)
for ki, bi in zip(k, b):
    c = BLUE if bi >= 0 else RED
    a.plot([ki, ki], [1e-4, abs(bi)], color=c, lw=2.4)
    a.plot(ki, abs(bi), 'o', color=c, ms=5)
a.set_yscale('log'); a.set_ylim(1e-3, 0.4)
a.set_xlabel(r'$k$'); a.set_ylabel(r'$|v_k/k!|$'); a.set_xticks(k)
a.set_title(r'5.  Divide by $k!$', loc='left', fontsize=10)
a.text(0.05, 0.20, 'now they decay: the Borel\nseries converges, with a\n'
                   'radius set by the nearest\nsingularity',
       transform=a.transAxes, fontsize=7.5, color=GREY, **BG, va='bottom')
a.text(0.05, 0.95, "they decay, but not smoothly —\nhold that thought",
       transform=a.transAxes, fontsize=7, color=GREY, **BG, va='top')

# ---------------------------------------------------------------- 6. the Borel plane + CUT
a = ax[5]
# The complex-conjugate pair is DELIBERATELY ABSENT here.  It is off-contour, so it causes no
# ambiguity for real positive x -- but drawn as shaded wedges it is the most salient object on
# the panel and viewers misattribute the ambiguity to it.  One mechanism at a time; the pair
# belongs to the later module on reading the plane's finer structure.  See DISTORTION_LEDGER.md.
# THE CUT -- mandatory, see DISTORTION_LEDGER.md
a.plot([A, 3.4], [0, 0], color=RED, lw=3.0, solid_capstyle='butt', zorder=4)
for xc in np.linspace(A, 3.35, 11):
    a.plot([xc, xc + 0.10], [0.0, 0.14], color=RED, lw=0.8, zorder=4)
a.plot(A, 0, 's', color=RED, ms=7, zorder=6)
a.annotate('branch point\n+ its CUT', xy=(2.5, 0.06), xytext=(1.70, 0.85), fontsize=7.5,
           color=RED, arrowprops=dict(arrowstyle='->', color=RED, lw=0.8))
a.axhline(0, color='k', lw=0.6); a.axvline(0, color='k', lw=0.6)
a.set_xlim(-0.5, 3.6); a.set_ylim(-1.7, 1.7); a.set_aspect('equal')
a.set_xlabel(r'$\mathrm{Re}\,\zeta$'); a.set_ylabel(r'$\mathrm{Im}\,\zeta$')
a.set_title('6.  Where it stops converging', loc='left', fontsize=10)
a.text(0.03, 0.05, 'one obstruction, sitting exactly on the path.\n(this plane has more in it — later)',
       transform=a.transAxes, fontsize=7.5, color=GREY, **BG)

# ---------------------------------------------------------------- 7. contour + two choices
a = ax[6]
a.plot([A, 3.4], [0, 0], color=RED, lw=3.0, solid_capstyle='butt', zorder=4)
a.plot(A, 0, 's', color=RED, ms=7, zorder=6)
tt = np.linspace(0, 3.3, 240)
for sgn, col, lab in ((+1, BLUE, 'pass above'), (-1, '#d9851f', 'pass below')):
    a.plot(tt, sgn * 0.50 * np.exp(-((tt - A) / 0.70)**2), color=col, lw=2.0, zorder=5,
           label=lab)
a.annotate('the contour must\ngo round it', xy=(A, 0.45), xytext=(0.05, 1.10), fontsize=8,
           arrowprops=dict(arrowstyle='->', color='k', lw=0.8))
a.text(0.25, -1.35, r'$\int_0^\infty B(\zeta)\,e^{-\zeta/x}\,d\zeta$', fontsize=12, **BG)
a.legend(fontsize=7.5, loc='upper right')
a.axhline(0, color='k', lw=0.6); a.axvline(0, color='k', lw=0.6)
a.set_xlim(-0.5, 3.6); a.set_ylim(-1.7, 1.7); a.set_aspect('equal')
a.set_xlabel(r'$\mathrm{Re}\,\zeta$')
a.set_title('7.  Two choices, no reason to prefer either', loc='left', fontsize=9.5)
a.text(0.02, 0.03, 'schematic', transform=a.transAxes, **SCHEM, **BG)

# ---------------------------------------------------------------- 8. the identity
a = ax[7]
amb = np.exp(-A / xs)
a.semilogy(xs, amb, color=RED, ls='--', lw=1.3, label=r'$e^{-A/x}$ alone')
a.semilogy(xs, 0.8 * xs**0.5 * amb, color=RED, lw=2.2, label=r'$C\,x^{b}e^{-A/x}$')
a.legend(fontsize=7, loc='lower right')
a.set_xlabel(r'$x=\eta^2$'); a.set_ylabel('difference between the two')
a.set_ylim(1e-5, 2.0)
a.set_title('8.  Their difference — the whole point', loc='left', fontsize=9.5)
a.text(0.05, 0.94, r'$\Delta \sim C\,x^{b}\,e^{-A/x}$', transform=a.transAxes, fontsize=13.5,
       color=RED, **BG, va='top')
a.text(0.05, 0.75, 'the SAME exponential as the\nfloor in panel 3 — and $A$ is the\n'
                   'singularity we just located.\n'
                   r'$C$ (a Stokes constant) and $b$ are'
                   '\nrefinements; the exponential\nis the content.',
       transform=a.transAxes, fontsize=7.8, color=GREY, **BG, va='top')
a.text(0.05, 0.50, 'so the power series alone\ncannot be the answer:\n'
                   'it is ambiguous at exactly\nthat size.\n\n'
                   r'$\Rightarrow$ add a term $\sim e^{-A/x}$.'
                   '\n' r'$\Rightarrow$ trans-series.',
       transform=a.transAxes, fontsize=8.2, **BG, va='top')
a.text(0.72, 0.03, 'schematic', transform=a.transAxes, **SCHEM, **BG)

for a in ax:
    a.grid(True, ls=':', lw=0.4, alpha=0.5)

fig.suptitle('PILOT STORYBOARD  --  the Borel plane and the ambiguity identity   '
             '(stills only; not yet animated)', fontsize=11.5, y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig('storyboard_contact.png', dpi=140)
print('wrote storyboard_contact.png')
