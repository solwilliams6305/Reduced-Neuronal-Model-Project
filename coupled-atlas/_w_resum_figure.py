"""Robustness of the median Borel-Pade resummation across Pade orders, and the
validation figure: resummed f(x)=Var/eta^2 vs the fp_cusp ground truth vs the
(divergent) naive partial sums, across x=eta^2 in [0,2.4].  beta=2 is x=2."""
import json, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from _w_resum import V, median_resum, naive_partial

plt.rcParams.update({'font.family':'serif','font.size':11,'mathtext.fontset':'dejavuserif'})

gt = {round(r['eta2'],4): r['f'] for r in json.load(open('_w_groundtruth.json'))}
gx = np.array(sorted(gt)); gy = np.array([gt[x] for x in gx])

# ---- robustness across Pade orders at beta=2 (x=2) and beta=2.67 (x=1.5) ----
print("resummed f at x=2 (beta=2) across Pade [L/M]:")
for (L,M) in [(3,3),(2,4),(4,2),(2,3),(3,2)]:
    try:
        med,_,_ = median_resum(V, 2.0, phi_deg=25.0, L=L, M=M)
        print(f"  [{L}/{M}]  f_resum(2) = {med.real:+.4f}   (ground truth 0.2370)")
    except Exception as e:
        print(f"  [{L}/{M}] failed: {e}")

# ---- curve for the figure (median resum, [3/3]) ----
xs = np.linspace(0.05, 2.4, 90)
fr = np.array([median_resum(V, x, phi_deg=25.0)[0].real for x in xs])
# naive partial sums (last order) -- diverges
fn = np.array([naive_partial(V, x)[-1] for x in xs])

fig, ax = plt.subplots(figsize=(6.2, 3.9))
ax.axvline(1.2, color='0.7', ls=':', lw=1.0)
ax.text(1.19, 0.05, 'pert. radius $x_c\\approx1.2$', rotation=90, va='bottom',
        ha='right', fontsize=8, color='0.4')
ax.axvline(2.0, color='#c0392b', ls='--', lw=1.0, alpha=0.7)
ax.text(2.0, 0.055, r'$\beta=2$', rotation=90, va='bottom', ha='right',
        fontsize=9, color='#c0392b')

ax.plot(xs, fr, '-', color='#1f5fb4', lw=2.0, zorder=3,
        label='median Borel resummation (7 terms)')
ax.plot(gx, gy, 'o', color='k', ms=6, zorder=4, label='ground truth (FP-PDE, MC-free)')
# naive partial sums (clipped, to show divergence)
ax.plot(xs, fn, '--', color='0.55', lw=1.3, zorder=2,
        label='naive 7-term partial sum')

ax.set_xlabel(r'$x=\eta^2=4/\beta$')
ax.set_ylabel(r'$f(x)=\mathrm{Var}(Y^\star)/\eta^2$')
ax.set_title('Resummed trans-series reconstructs $\\mathcal{W}$ beyond the radius', fontsize=10.5)
ax.set_xlim(0, 2.4); ax.set_ylim(-0.35, 0.45)
ax.axhline(0, color='k', lw=0.5)
ax.legend(loc='lower left', fontsize=8.3, framealpha=0.93)
ax.grid(True, ls=':', lw=0.5, alpha=0.6)
fig.subplots_adjust(left=0.12, right=0.97, top=0.91, bottom=0.14)
fig.savefig('../W_fig_resum.pdf'); fig.savefig('../W_fig_resum.png', dpi=150)
print("\nwrote ../W_fig_resum.pdf")

# print the beta=2 comparison line for the paper
med2,_,_ = median_resum(V, 2.0, phi_deg=25.0)
print(f"\nAT beta=2:  resummed {med2.real:.4f}  vs ground truth {gt[2.0]:.4f}  "
      f"(rel err {abs(med2.real-gt[2.0])/gt[2.0]*100:.1f}%);  naive 7-term = {naive_partial(V,2.0)[-1]:.1f}")
