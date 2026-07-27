"""
Figures for the swallowtail trans-series paper (three, mirroring the cusp paper's set).

  swtl_fig_ladder.pdf  -- the coefficient ladder, sign-coloured, with the node/anti-node marked,
                          plus the kappa_3 ladder underneath on its own scale.
  swtl_fig_borel.pdf   -- the Borel plane: the conjugate pair with its method spread, the real
                          instanton ray, the deterministic anchor phase pi/(q+2), and the cusp
                          values for contrast.
  swtl_fig_resum.pdf   -- resummation vs ground truth across x, 7-coefficient vs 8-coefficient,
                          with the |err|<10% ranges marked.

Data are taken from the committed ladders and the validation sweep, not recomputed.

Run:  python3 swtl_figures.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

V = [0.04969, 0.06323, 0.11079, 0.17705, 0.03635, -1.49615, -5.33940, -21.53045]
T = [0.01687, 0.05662, 0.16891, 0.40452, 0.51824, -5.22509]
CUSP_V = [0.134, 0.111, 0.104, -0.030, -0.451, -1.19, -1.90]

SWEEP = [  # x, beta, f_exact, f_7coef, f_8coef
    (0.09, 44.4, 0.056595, 0.0564, 0.0564),
    (0.16, 25.0, 0.063866, 0.0635, 0.0632),
    (0.25, 16.0, 0.074659, 0.0740, 0.0736),
    (0.36, 11.1, 0.088648, 0.0831, 0.0876),
    (0.49,  8.2, 0.104324, 0.0859, 0.1043),
    (0.64,  6.3, 0.119848, 0.0825, 0.1223),
    (0.81,  4.9, 0.133752, 0.0753, 0.1405),
    (1.00,  4.0, 0.145188, 0.0666, 0.1583),
    (1.44,  2.8, 0.159783, 0.0497, 0.1913),
    (2.00,  2.0, 0.164683, 0.0354, 0.2234),
    (2.50,  1.6, 0.162948, 0.0271, 0.2470),
]

BLUE, RED, GREY = '#1f5fa9', '#c0392b', '#888888'


def fig_ladder():
    fig, ax = plt.subplots(2, 1, figsize=(7.0, 5.4), sharex=False,
                           gridspec_kw={'height_ratios': [1.25, 1]})
    k = np.arange(len(V))
    for a, vals, name, nodeidx, flipidx in (
            (ax[0], V, r'variance ladder $v_k$', 4, 5),
            (ax[1], T, r'$\kappa_3$ ladder $t_j$', None, 5)):
        kk = np.arange(len(vals))
        cols = [BLUE if v > 0 else RED for v in vals]
        a.axhline(0, color='k', lw=0.6)
        a.vlines(kk, 0, vals, color=cols, lw=2.4)
        a.plot(kk, vals, 'o', ms=5, mfc='w', mec='k', mew=0.9, zorder=3)
        a.set_yscale('symlog', linthresh=0.05)
        a.set_ylabel(name)
        a.grid(alpha=0.25, ls=':')
        if nodeidx is not None:
            a.annotate('node', xy=(nodeidx, vals[nodeidx]), xytext=(nodeidx - 0.75, 1.2),
                       arrowprops=dict(arrowstyle='->', lw=0.9), fontsize=9)
            a.annotate('anti-node', xy=(nodeidx + 1, vals[nodeidx + 1]),
                       xytext=(nodeidx + 0.25, -12), arrowprops=dict(arrowstyle='->', lw=0.9),
                       fontsize=9)
        a.axvline(flipidx - 0.5, color=GREY, ls='--', lw=1.0)
        a.text(flipidx - 0.45, a.get_ylim()[1] * 0.35, 'sign change', fontsize=8, color=GREY)
    ax[1].set_xlabel(r'rung $k$ (resp. $j$)')
    ax[0].set_title(r'Swallowtail ($q=3$) ladders: both turn negative at the same rung', fontsize=10)
    fig.tight_layout()
    fig.savefig('../swtl_fig_ladder.pdf')
    print('wrote swtl_fig_ladder.pdf')


def fig_borel():
    fig, ax = plt.subplots(figsize=(5.8, 5.4))
    ax.axhline(0, color='k', lw=0.6); ax.axvline(0, color='k', lw=0.6)
    # swallowtail pair: |zeta| ~ 1.5, theta band 33-42
    zm, tlo, thi, tbest = 1.5, 33.0, 42.0, 39.4
    th = np.radians(np.linspace(tlo, thi, 80))
    for sgn in (+1, -1):
        ax.fill(np.concatenate([[0], zm * np.cos(th)]),
                np.concatenate([[0], sgn * zm * np.sin(th)]),
                color=BLUE, alpha=0.16, lw=0)
    for sgn in (+1, -1):
        ax.plot(zm * np.cos(np.radians(tbest)), sgn * zm * np.sin(np.radians(tbest)), '*',
                ms=15, color=BLUE, zorder=5,
                label=r'swallowtail pair $\theta\approx39^\circ$' if sgn > 0 else None)
    # deterministic anchor pi/(q+2) = 36 deg
    for sgn in (+1, -1):
        ax.plot([0, 2.6 * np.cos(np.radians(36))], [0, sgn * 2.6 * np.sin(np.radians(36))],
                ls='--', lw=1.2, color='#2e7d32',
                label=r'anchor $\pi/(q{+}2)=36^\circ$' if sgn > 0 else None)
    # cusp pair for contrast
    for sgn in (+1, -1):
        ax.plot(1.9 * np.cos(np.radians(50)), sgn * 1.9 * np.sin(np.radians(50)), 'o',
                ms=8, mfc='none', mec=RED, mew=1.6,
                label=r'cusp pair $50^\circ,\,|\zeta|=1.9$' if sgn > 0 else None)
    # real instanton
    ax.plot([1.5], [0], 's', ms=9, color='k', label=r'real FW instanton (on $\mathbb{R}_+$)')
    ax.set_xlim(-0.4, 2.8); ax.set_ylim(-2.2, 2.2); ax.set_aspect('equal')
    ax.set_xlabel(r'$\mathrm{Re}\,u$'); ax.set_ylabel(r'$\mathrm{Im}\,u$')
    ax.set_title('Borel plane of the swallowtail law', fontsize=10)
    ax.legend(fontsize=8, loc='upper right', framealpha=0.95)
    ax.grid(alpha=0.25, ls=':')
    fig.tight_layout(); fig.savefig('../swtl_fig_borel.pdf')
    print('wrote swtl_fig_borel.pdf')


def fig_resum():
    x = np.array([r[0] for r in SWEEP]); fe = np.array([r[2] for r in SWEEP])
    f7 = np.array([r[3] for r in SWEEP]); f8 = np.array([r[4] for r in SWEEP])
    fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.7))
    ax[0].plot(x, fe, 'k-', lw=2, label='ground truth (FP-PDE)')
    ax[0].plot(x, f7, 'o--', color=GREY, ms=4, label='7 coefficients [3/3]')
    ax[0].plot(x, f8, 's-', color=BLUE, ms=4, label='8 coefficients [3/4]')
    ax[0].set_xlabel(r'$x=\eta^2=4/\beta$'); ax[0].set_ylabel(r'$f(x)=\mathrm{Var}/\eta^2$')
    ax[0].legend(fontsize=8); ax[0].grid(alpha=0.25, ls=':')
    ax[0].set_title('Resummation vs ground truth', fontsize=10)

    ax[1].axhline(0, color='k', lw=0.6)
    for lim in (10, -10):
        ax[1].axhline(lim, color='#2e7d32', ls=':', lw=1.0)
    ax[1].plot(x, (f7 / fe - 1) * 100, 'o--', color=GREY, ms=4, label='7 coefficients')
    ax[1].plot(x, (f8 / fe - 1) * 100, 's-', color=BLUE, ms=4, label='8 coefficients')
    ax[1].axvline(1.00, color=BLUE, ls='--', lw=0.9)
    ax[1].axvline(0.36, color=GREY, ls='--', lw=0.9)
    ax[1].text(1.03, -60, r'$|err|<10\%$ to $x=1.00$', fontsize=8, color=BLUE)
    ax[1].text(0.39, -78, r'$x=0.36$', fontsize=8, color=GREY)
    ax[1].set_xlabel(r'$x=\eta^2$'); ax[1].set_ylabel('relative error (%)')
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.25, ls=':')
    ax[1].set_title('8-coefficient error: systematic, not erratic', fontsize=10)
    fig.tight_layout(); fig.savefig('../swtl_fig_resum.pdf')
    print('wrote swtl_fig_resum.pdf')


if __name__ == '__main__':
    fig_ladder(); fig_borel(); fig_resum()
