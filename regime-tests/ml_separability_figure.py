"""ml_separability_figure.py — the kill-or-confirm verdict figure."""
from __future__ import annotations
import os
import numpy as np
import ml_atlas_separability as m
from bifurcation_classifier import skew


def make_figure():
    if not os.path.exists(m.CACHE):
        m.build()
    atlas, tests = m._load()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    C_I, C_II = "#dc2626", "#2563eb"
    fig, ax = plt.subplots(1, 3, figsize=(16.2, 4.9))

    # ---- A: per-drive shape-only accuracy: good near rheobase, collapses at high drive ----------
    a = ax[0]
    for cls, drives, col in [("I", m.DRIVES_I, C_I), ("II", m.DRIVES_II, C_II)]:
        accs = []
        for di in range(len(drives)):
            sw = tests[(cls, di)]
            accs.append(np.mean([m.classify_shape(s, atlas)[0] == cls for s in sw]) * 100 if sw else np.nan)
        a.plot(range(len(drives)), accs, "o-", color=col, ms=8, lw=1.8, label=f"Type-{cls}")
    a.axhline(50, color="#94a3b8", lw=1.0, ls=":"); a.text(2.7, 52, "chance", fontsize=8)
    a.axvspan(-0.3, 1.3, color="#dcfce7", alpha=0.7, zorder=0)
    a.text(0.0, 8, "near\nrheobase", fontsize=8, color="#15803d")
    a.set_xlabel("drive level (0 = near rheobase $\\to$ high)"); a.set_ylabel("shape-only class accuracy (%)")
    a.set_ylim(0, 105); a.set_title("(A) Class separates near threshold, collapses at high drive")
    a.legend(fontsize=9, framealpha=0.95, edgecolor="#cbd5e1")

    # ---- B: near-rheobase ISI shapes DO differ ---------------------------------------------------
    b = ax[1]; bins = np.linspace(0, 3.0, 38)
    for cls, col in [("I", C_I), ("II", C_II)]:
        pooled = np.concatenate(tests[(cls, 0)]); skv = np.mean([skew(s) for s in tests[(cls, 0)]])
        b.hist(pooled / pooled.mean(), bins=bins, density=True, histtype="step", lw=2.2, color=col,
               label=fr"Type-{cls} near rheobase (skew {skv:.2f})")
    b.set_xlabel(r"ISI / $\langle$ISI$\rangle$"); b.set_ylabel("density")
    b.set_title("(B) Near rheobase the shapes are distinct")
    b.legend(fontsize=8.6, framealpha=0.95, edgecolor="#cbd5e1")

    # ---- C: the verdict --------------------------------------------------------------------------
    c = ax[2]; c.axis("off")
    conf, acc, margins, rows = m.test(verbose=False)
    nt = []; al = []
    for (cls, di), sw in tests.items():
        for s in sw:
            ok = m.classify_shape(s, atlas)[0] == cls; al.append(ok)
            if di <= 1:
                nt.append(ok)
    c.text(0.5, 0.96, "VERDICT: PARTIAL", ha="center", va="top", fontsize=14, fontweight="bold",
           color="#b45309")
    txt = (
        "The Type-I/Type-II distinction SURVIVES in a conductance\n"
        "model's ISIs — but weakly, and only near threshold.\n\n"
        f"  shape-only accuracy (all drive):   {np.mean(al)*100:.0f}%\n"
        f"  shape-only accuracy (near rheobase): {np.mean(nt)*100:.0f}%\n"
        f"  vs idealised normal forms:          89-100%\n\n"
        "  Type-II is the hard class (regular-firing\n"
        "  regime at mid drive reads as Type-I).\n"
        "  f-I onset still separates directionally\n"
        "  (0.32 vs 0.56), gap narrowed from 0.47/0.87.\n\n"
        "Not killed. Realistic spec: biophysically-built atlas,\n"
        "near-rheobase sweeps, shape + f-I onset combined,\n"
        "~80% (population characterisation, not per-cell certainty)."
    )
    c.text(0.04, 0.84, txt, ha="left", va="top", fontsize=9.6, color="#1e293b", family="monospace")

    fig.tight_layout()
    out = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures",
                                       "ml_atlas_separability.png"))
    fig.savefig(out, dpi=140); print("saved", out)


if __name__ == "__main__":
    make_figure()
