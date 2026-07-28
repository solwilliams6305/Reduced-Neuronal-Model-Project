"""Static render of the junction's ratio panel (module delta1, second half).

Static-first, per PROBE_PROTOCOL.md: the interactive already exists, but the still is what gets
probed and what survives a broken preview. Reads the same numbers the browser does, from the
same source, so the two cannot drift.

Run:  python3 ratio_still.py   ->  ratio_still.png
"""
import json
import math
import pathlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'axes.linewidth': 0.8,
                     'mathtext.fontset': 'dejavuserif', 'figure.dpi': 150})
BLUE, RED, PUR, GRN, GREY = '#2f6fb0', '#c0392b', '#7a5cc0', '#1f9d63', '0.45'

v = np.array([0.134, 0.111, 0.104, -0.030, -0.451, -1.19, -1.90])
k = np.arange(len(v))
b = v / np.array([float(math.factorial(int(i))) for i in k])
r = b[1:] / b[:-1]

RHO, TH = 1.9, np.deg2rad(50.0)
best = None
for phi in np.linspace(0, 2 * np.pi, 2001):
    m = RHO ** (-k.astype(float)) * np.cos(k * TH - phi)
    sc = float(np.dot(m, b) / np.dot(m, m))
    res = float(np.sum((sc * m - b) ** 2) / np.sum(b ** 2))
    if best is None or res < best[0]:
        best = (res, phi, sc * m)
res, phi, mv = best
rm = mv[1:] / mv[:-1]

fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.4, 4.3),
                             gridspec_kw={'width_ratios': [1.55, 1]})

kk = np.arange(len(r))
a1.axhline(0, color='k', lw=.7)
a1.axhline(1 / RHO, color=GRN, ls='--', lw=1.7,
           label=rf'$1/A={1/RHO:.3f}$ — where a single REAL singularity would settle')
a1.plot(kk, rm, 'o-', color=PUR, ms=5, lw=1.8, alpha=.85,
        label=rf'conjugate pair $1.9\,e^{{\pm i50^\circ}}$ (residual {100*res:.1f}%)')
a1.plot(kk, r, 'o-', color=BLUE, ms=7, lw=2.4, label='the actual cusp ladder')
a1.set_xticks(kk); a1.set_xlabel(r'$k$')
a1.set_ylabel(r'$b_{k+1}/b_k$,   $b_k=v_k/k!$')
a1.set_title('They do not settle — and that is the finding', loc='left', fontsize=11)
a1.legend(fontsize=8, loc='upper right', framealpha=.94)
a1.grid(True, ls=':', lw=.4, alpha=.6)
a1.annotate('sign flip', xy=(2, r[2]), xytext=(1.15, -0.75), fontsize=8, color=RED,
            arrowprops=dict(arrowstyle='->', color=RED, lw=.8))
a1.annotate('spike', xy=(3.1, 3.0), xytext=(3.5, 2.0), fontsize=8, color=RED,
            arrowprops=dict(arrowstyle='->', color=RED, lw=.8))

root = np.array([abs(b[i]) ** (-1.0 / i) for i in range(1, len(b))])
a2.plot(np.arange(1, len(b)), root, 'o-', color=BLUE, ms=6, lw=2)
a2.axhline(RHO, color=RED, ls='--', lw=1.6, label=rf'true $|\zeta|\approx{RHO}$')
a2.set_xlabel(r'$k$'); a2.set_ylabel(r'$|b_k|^{-1/k}$')
a2.set_title('…and the modulus is not readable yet', loc='left', fontsize=11)
a2.legend(fontsize=8.5, loc='upper right')
a2.grid(True, ls=':', lw=.4, alpha=.6)
a2.text(.05, .06, 'seven coefficients give 2.5–2.7\nand are not converging',
        transform=a2.transAxes, fontsize=8.5, color=RED, style='italic',
        bbox=dict(fc='white', ec='none', alpha=.9, pad=2))

fig.suptitle('THE COEFFICIENTS KNOW THE GEOMETRY — but seven of them cannot say it plainly',
             fontsize=12, y=.985)
fig.tight_layout(rect=[0, 0, 1, .94])
fig.savefig('ratio_still.png', dpi=150)
print(f'wrote ratio_still.png  |  model residual {res:.4f}, phi {np.degrees(phi):.1f} deg')
print(f'  data  ratios: {np.round(r, 3)}')
print(f'  model ratios: {np.round(rm, 3)}')
print(f'  root test    : {np.round(root, 2)}  vs true {RHO}')
assert res < 0.02, 'conjugate-pair model no longer fits'
