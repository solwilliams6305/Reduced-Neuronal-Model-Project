"""PROTOTYPE -- the R <-> u twin panel  (module alpha4a of DESIGN_inner_chain.md)

Tests ONE claim: that Cole-Hopf stops being a trick once you SEE that
    R = -u'/u  blows up exactly where  u  crosses zero.
Blow-ups are hard to see and hard to compute with; zeros are easy.  The transformation buys
you zeros.  If that lands, "escape <=> first node" is self-evident rather than asserted.

STILLS ONLY, per PROBE_PROTOCOL.md.

Everything here is EXACT, not schematic:
  row 1   u(Y) = Ai(Y), the recessive solution of  u'' = Y u
  row 2   R(Y) = -Ai'(Y)/Ai(Y), which solves  R' = R^2 - Y
          (differentiate R = -u'/u and substitute u'' = Y u)
  row 3   the same construction with noise -- u'' = (Y + eta*xi)u -- so the first zero wanders.
          That wandering IS the escape statistic.

Shared x-axis throughout: the vertical alignment IS the argument, so the rows must never drift
out of register.  Every Airy zero inside the plotted window is marked -- an unmarked asymptote
in row 2 would silently falsify the panel's whole claim.

Run:  python3 twin_panel.py  ->  twin_panel.png
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import airy, ai_zeros
from scipy.stats import skew

plt.rcParams.update({
    'font.family': 'serif', 'font.size': 10, 'axes.linewidth': 0.8,
    'mathtext.fontset': 'dejavuserif', 'figure.dpi': 140,
})
BLUE, RED, GREY, AMB = '#1f5fb4', '#c0392b', '0.45', '#d9851f'
BG = dict(bbox=dict(facecolor='white', edgecolor='none', alpha=0.9, pad=1.8))

Y_LO, Y_HI = -9.0, 2.0
Y = np.linspace(Y_LO, Y_HI, 4000)
N = len(Y)
h = Y[1] - Y[0]
Ai, Aip, _, _ = airy(Y)

# ai_zeros returns the zeros themselves, which are NEGATIVE -- do not negate.  Take plenty and
# keep every one inside the window, so row 2 has no asymptote without a dot above it.
all_zeros = ai_zeros(12)[0]
zeros = np.sort(all_zeros[(all_zeros > Y_LO) & (all_zeros < Y_HI)])
Z_FIRST = zeros[-1]              # first zero met sweeping Y downward from the recessive end

fig, ax = plt.subplots(3, 1, figsize=(9.2, 10.0), sharex=True,
                       gridspec_kw={'height_ratios': [1.0, 1.0, 1.25], 'hspace': 0.15})

# ------------------------------------------------------------------ row 1: u has zeros
a = ax[0]
a.plot(Y, Ai, color=BLUE, lw=2.0)
a.axhline(0, color='k', lw=0.7)
a.plot(zeros, np.zeros_like(zeros), 'o', color=RED, ms=7, zorder=5)
a.set_ylabel(r'$u(Y)=\mathrm{Ai}(Y)$')
a.set_title("$u$ is smooth. It crosses zero.", loc='left', fontsize=11.5)
a.text(0.015, 0.07, f'nothing dramatic happens here —\njust {len(zeros)} ordinary crossings.\n'
                    'this one curve carries through all three rows.',
       transform=a.transAxes, fontsize=8.5, color=GREY, **BG)

# ------------------------------------------------------------------ row 2: R blows up there
a = ax[1]
R = -Aip / Ai
a.plot(Y, np.ma.masked_where(np.abs(R) > 12, R), color=AMB, lw=2.0)
a.axhline(0, color='k', lw=0.7)
for z in zeros:
    a.axvline(z, color=RED, lw=1.0, ls='--', alpha=0.85, zorder=1)
a.set_ylim(-12, 12)
a.set_ylabel(r"$R(Y)=-u'/u$")
a.set_title(r'$R$ blows up. At exactly the same $Y$.', loc='left', fontsize=11.5)
a.text(0.015, 0.07, r"$R$ solves $R'=R^2-Y$ — the inner Riccati."
                    "\nevery asymptote sits on a zero above; there are no others.",
       transform=a.transAxes, fontsize=8.5, color=GREY, **BG)
a.text(0.985, 0.90, 'blow-ups are hard to see.\nzeros are easy.\n'
                    'Cole–Hopf buys you zeros.',
       transform=a.transAxes, fontsize=10, color=RED, ha='right', va='top', **BG)

# ------------------------------------------------------------------ row 3: noise moves it
a = ax[2]
rng = np.random.default_rng(20260727)
eta, n_real = 0.55, 400
u = np.empty((n_real, N)); p = np.empty((n_real, N))
u[:, -1], p[:, -1] = Ai[-1], Aip[-1]
noise = rng.standard_normal((n_real, N)) * eta / np.sqrt(abs(h))
for i in range(N - 1, 0, -1):                      # integrate downward from the recessive end
    p[:, i-1] = p[:, i] - h * (Y[i] + noise[:, i]) * u[:, i]
    u[:, i-1] = u[:, i] - h * p[:, i]

cross = np.signbit(u[:, :-1]) != np.signbit(u[:, 1:])
has = cross.any(axis=1)
idx = (N - 2) - np.argmax(cross[:, ::-1], axis=1)  # rightmost crossing = first met sweeping down
first_zeros = Y[idx][has]

un = u / np.max(np.abs(u), axis=1, keepdims=True)
for j in range(14):                                 # the ensemble is SCATTER, hence grey
    a.plot(Y, un[j], color=GREY, lw=0.7, alpha=0.40, zorder=2)
# THE INVARIANT: the very same curve drawn in rows 1 and 2, in the same colour and weight.
# Row 3 must not introduce a new object -- it shows this one trajectory with noise scattered
# around it.  See DESIGN_inner_chain.md sec.2.
a.plot(Y, Ai / np.max(np.abs(Ai)), color=BLUE, lw=2.0, zorder=6)
a.plot(Z_FIRST, 0, 'o', color=RED, ms=7, zorder=7)
a.axhline(0, color='k', lw=0.7)
a.axvline(Z_FIRST, color=RED, lw=1.4, ls='--', zorder=4)
a.plot(first_zeros, np.zeros_like(first_zeros), '|', color=RED, ms=11, mew=0.9, alpha=0.35,
       zorder=5)
a.set_ylim(-1.15, 1.15); a.set_xlim(Y_LO, Y_HI)
a.set_xlabel(r'$Y$   (the swept parameter)')
a.set_ylabel('$u$, normalised')
a.set_title('Add noise: the first zero wanders. That wandering is the escape law.',
            loc='left', fontsize=11.5)
a.text(0.015, 0.06, "blue: the SAME $u$ as rows 1–2.  grey: 14 noisy runs.\n"
                    rf"$u''=(Y+\eta\,\xi)u$,  $\eta={eta}$,  {has.sum()} runs;"
                    "  ticks mark every first zero",
       transform=a.transAxes, fontsize=8.5, color=GREY, **BG)
a.annotate('its first zero —\nthe same red dot', xy=(Z_FIRST, 0.42), xytext=(Z_FIRST - 3.6, 0.90),
           fontsize=8.5, color=RED, arrowprops=dict(arrowstyle='->', color=RED, lw=0.9))

# histogram of the first zero -- placed clear of the tick row at y=0
sk = skew(first_zeros)
ins = a.inset_axes([0.60, 0.045, 0.375, 0.335])
ins.patch.set_facecolor('white'); ins.patch.set_alpha(0.96)
ins.hist(first_zeros, bins=26, color=RED, alpha=0.6, edgecolor='none')
ins.axvline(Z_FIRST, color=RED, lw=1.2, ls='--')
ins.set_yticks([]); ins.tick_params(labelsize=6.5)
ins.set_title(f'first zero  (skew ${sk:+.2f}$, s.e. ${np.sqrt(6/len(first_zeros)):.2f}$)',
              fontsize=7, color=GREY, pad=2.5)
# HONESTY: eta=0.55 on this grid is NOT the scaling regime, so this histogram does not
# exhibit Tracy-Widom and must never be captioned as if it did.  The measured skew is not
# significantly different from zero.  See DISTORTION_LEDGER.md.
ins.text(0.03, 0.055, 'not the scaling limit —\nthis is not TW', transform=ins.transAxes,
         fontsize=6.5, style='italic', color=RED)

for a_ in ax:
    a_.grid(True, ls=':', lw=0.4, alpha=0.5)

fig.suptitle('THE TRADE:  blow-ups for zeros', fontsize=13.5, y=0.972)
fig.subplots_adjust(left=0.11, right=0.975, top=0.932, bottom=0.058)
fig.savefig('twin_panel.png', dpi=140)
print(f'wrote twin_panel.png  |  zeros in window: {np.round(zeros, 3)}  |  '
      f'{has.sum()}/{n_real} runs crossed  |  mean first zero {first_zeros.mean():.3f} '
      f'vs det {Z_FIRST:.3f}  |  skew {sk:+.3f}')
