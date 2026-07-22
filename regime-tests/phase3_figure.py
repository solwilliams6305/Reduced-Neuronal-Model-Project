"""phase3_figure.py — render the Phase-3 early-warning figure from the cached ensembles."""
from __future__ import annotations
import os
import numpy as np
import phase3_tipping as p3

HERE = os.path.dirname(os.path.abspath(__file__))


def _nu_true(t):
    return p3.NU_HI + (p3.NU_LO - p3.NU_HI) * (np.asarray(t) / p3.T_TOTAL)


def make_figure():
    if not os.path.exists(p3.CACHE):
        p3.build_sims()
    res, ap, nl, cf = p3.benchmark(verbose=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    C_NU, C_V, C_A = "#16a34a", "#dc2626", "#7c3aed"
    fig, ax = plt.subplots(2, 2, figsize=(13.2, 9.4))

    # ---- A: one approach trial — nu_hat tracks the drifting true nu down to the SNIC -------------
    a = ax[0, 0]
    r = max(ap, key=len)                                   # a well-sampled trial
    t = r[:, 0]; a.plot(t, _nu_true(t), color="#0f172a", lw=2.2, label=r"true $\nu(t)$ (hidden)")
    a.plot(t, r[:, 1], "o-", color=C_NU, ms=4, lw=1.2, label=r"inverted $\hat\nu(t)$ (from ISIs)")
    a.axhline(0, color="#94a3b8", lw=1.0, ls=":")
    a.axhline(p3.NU_ABS, color=C_NU, lw=1.0, ls="--")
    a.axvline(p3.T_CROSS, color="#dc2626", lw=1.2, ls="--")
    a.text(p3.T_CROSS - 30, 1.7, "SNIC\n(firing ceases)", color="#dc2626", fontsize=8, ha="right")
    a.text(t[0] + 20, p3.NU_ABS + 0.06, r"$\hat\nu$ alarm", color=C_NU, fontsize=8)
    a.set_xlabel("time"); a.set_ylabel(r"$\nu$ (distance to bifurcation)")
    a.set_title("(A) The inversion tracks the approach in real time")
    a.legend(fontsize=8.6, framealpha=0.95, edgecolor="#cbd5e1", loc="upper right")

    # ---- B: ensemble signals vs true nu — nu_hat falls, Var rises, AC1 flat ----------------------
    b = ax[0, 1]
    pts = []
    for rr in ap:
        nt = _nu_true(rr[:, 0]); keep = nt > -0.05
        for nu_t, row in zip(nt[keep], rr[keep]):
            pts.append((nu_t, row[1], row[2], row[3]))
    pts = np.array(pts)
    edges = np.linspace(0, 2.0, 11); ctr = 0.5 * (edges[:-1] + edges[1:])
    binmean = lambda col: np.array([np.nanmean(pts[(pts[:, 0] >= edges[i]) & (pts[:, 0] < edges[i+1]), col])
                                    for i in range(len(ctr))])
    nu_b, var_b, ac_b = binmean(1), binmean(2), binmean(3)
    norm = lambda v: (v - np.nanmin(v)) / (np.nanmax(v) - np.nanmin(v))
    b.plot(ctr, nu_b, "o-", color=C_NU, lw=1.8, label=r"$\hat\nu$ (calibrated)")
    b.plot(ctr, [c for c in ctr], ":", color="#94a3b8", lw=1.0, label="ideal $\\hat\\nu=\\nu$")
    b.set_xlabel(r"true $\nu$ (approach $\rightarrow$ 0)"); b.set_ylabel(r"$\hat\nu$", color=C_NU)
    b.tick_params(axis="y", labelcolor=C_NU); b.invert_xaxis()
    b2 = b.twinx()
    b2.plot(ctr, norm(var_b), "s--", color=C_V, ms=5, lw=1.6, label="ISI Var (norm.)")
    b2.plot(ctr, norm(ac_b), "^--", color=C_A, ms=5, lw=1.6, label="lag-1 AC (norm.)")
    b2.set_ylabel("normalised EWS")
    b.set_title("(B) Ensemble: $\\hat\\nu$ tracks $\\nu$; Var rises; AC1 ~flat")
    h1, l1 = b.get_legend_handles_labels(); h2, l2 = b2.get_legend_handles_labels()
    b.legend(h1 + h2, l1 + l2, fontsize=8.0, framealpha=0.95, edgecolor="#cbd5e1", loc="upper left")

    # ---- C: matched-FPR benchmark — detection (approach) & false alarms (null + confound) --------
    c = ax[1, 0]
    keys = ["nu", "var", "ac1", "nu_abs"]; labs = [r"$\hat\nu$ trend", "ISI Var", "lag-1 AC", r"$\hat\nu<$thr"]
    x = np.arange(len(keys)); w = 0.27
    det = [res[k]["detect"] * 100 for k in keys]
    fnull = [res[k]["fpr"] * 100 for k in keys]
    fconf = [res[k]["fpr_cf"] * 100 for k in keys]
    c.bar(x - w, det, w, color="#16a34a", label="detection (approach)")
    c.bar(x, fnull, w, color="#94a3b8", label="false alarm (null)")
    c.bar(x + w, fconf, w, color="#dc2626", label="false alarm (rate-drift confound)")
    c.set_xticks(x); c.set_xticklabels(labs, fontsize=9); c.set_ylabel("% of trials")
    c.set_title("(C) Matched-FPR benchmark")
    c.legend(fontsize=8.0, framealpha=0.95, edgecolor="#cbd5e1", loc="upper center")
    for xi, v in zip(x - w, det):
        c.text(xi, v + 1.5, f"{v:.0f}", ha="center", fontsize=8)
    for xi, v in zip(x + w, fconf):
        c.text(xi, v + 1.5, f"{v:.0f}", ha="center", fontsize=8)

    # ---- D: the confound — rate drift inflates Var (false alarm) but not the shape -> nu_hat flat -
    d = ax[1, 1]
    def ens_by_window(rows_list, col):
        L = min(len(r) for r in rows_list)
        return np.mean([r[:L, col] for r in rows_list], axis=0)
    vv = ens_by_window(cf, 2); nn = ens_by_window(cf, 1)
    wi = np.arange(len(vv))
    d.plot(wi, (vv - vv.min()) / (vv.max() - vv.min()), "s-", color=C_V, ms=4, label="ISI Var (norm.) — rises")
    d.set_xlabel("window (no approach; firing rate drifts)"); d.set_ylabel("normalised ISI Var", color=C_V)
    d.tick_params(axis="y", labelcolor=C_V)
    d2 = d.twinx()
    d2.plot(wi, nn, "o-", color=C_NU, ms=4, label=r"$\hat\nu$ — flat (not fooled)")
    d2.axhline(p3.NU_NULL, color=C_NU, lw=1.0, ls=":")
    d2.set_ylabel(r"$\hat\nu$", color=C_NU); d2.tick_params(axis="y", labelcolor=C_NU)
    d2.set_ylim(0, p3.NU_HI + 0.5)
    d.set_title("(D) Rate-drift confound fools Var, not $\\hat\\nu$")
    h1, l1 = d.get_legend_handles_labels(); h2, l2 = d2.get_legend_handles_labels()
    d.legend(h1 + h2, l1 + l2, fontsize=8.2, framealpha=0.95, edgecolor="#cbd5e1", loc="center right")

    fig.tight_layout()
    out = os.path.abspath(os.path.join(HERE, "..", "figures", "phase3_tipping.png"))
    fig.savefig(out, dpi=140); print("\nsaved", out)


if __name__ == "__main__":
    make_figure()
