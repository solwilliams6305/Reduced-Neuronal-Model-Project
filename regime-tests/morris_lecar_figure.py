"""morris_lecar_figure.py — the conductance-model out-of-distribution figure."""
from __future__ import annotations
import os
import numpy as np
import morris_lecar_check as ml
import bifurcation_classifier as bc
from bifurcation_classifier import skew

HERE = os.path.dirname(os.path.abspath(__file__))


def make_figure():
    if not os.path.exists(ml.CACHE):
        ml.build_sims()
    reps, infam, garbage = ml.check(verbose=True)
    cells = ml.load_sims()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    C_I, C_II, C_AT, C_G = "#dc2626", "#2563eb", "#16a34a", "#64748b"
    fig, ax = plt.subplots(1, 3, figsize=(16.2, 4.9))

    # ---- A: the goodness-of-fit residual puts ML in the reject zone, near garbage ---------------
    a = ax[0]
    groups = [("in-family\n(normal forms)", [infam], C_AT),
              ("ML Type-I", [r["residual"] for l, r in reps if l == "Type-I"], C_I),
              ("ML Type-II", [r["residual"] for l, r in reps if l == "Type-II"], C_II),
              ("garbage\n(Poisson)", [garbage], C_G)]
    for x, (lab, vals, col) in enumerate(groups):
        a.scatter([x] * len(vals), vals, c=col, s=90, edgecolor="white", lw=1.0, zorder=3)
    a.axhline(0.07, color="#0f172a", lw=1.3, ls="--"); a.text(3.0, 0.073, "reject gate", fontsize=9)
    a.fill_between([-0.5, 3.5], 0.07, 0.2, color="#fee2e2", alpha=0.6, zorder=0)
    a.set_xticks(range(4)); a.set_xticklabels([g[0] for g in groups], fontsize=8.6)
    a.set_ylabel(r"goodness-of-fit residual ($W_1$)"); a.set_ylim(0, 0.13)
    a.set_title("(A) Conductance model is out-of-family")

    # ---- B: the ISI-shape mismatch ---------------------------------------------------------------
    b = ax[1]; bins = np.linspace(0, 3.0, 40)
    atI = bc.snic_intervals(0.6, N=3000, seed=3); atII = bc.hopf_intervals(0.4, seed=3)
    b.hist(atI / atI.mean(), bins=bins, density=True, histtype="step", lw=2.0, color=C_AT,
           label="atlas SNIC (Type-I)")
    b.hist(atII / atII.mean(), bins=bins, density=True, histtype="step", lw=2.0, color="#0891b2",
           label="atlas Hopf (Type-II)")
    mlI = np.concatenate([sw[-1] for l, sw, c in cells if l == "Type-I"])      # higher-drive sweep
    b.hist(mlI / mlI.mean(), bins=bins, density=True, histtype="stepfilled", alpha=0.4, color=C_I,
           label="ML Type-I (higher drive)")
    b.set_xlabel(r"ISI / $\langle$ISI$\rangle$"); b.set_ylabel("density")
    b.set_title("(B) ML leaves the first-passage shape")
    b.legend(fontsize=8.4, framealpha=0.95, edgecolor="#cbd5e1")

    # ---- C: the mechanism — skew flips negative as ML enters regular firing ----------------------
    c = ax[2]
    for lab, sw, cur in cells:
        rheo = cur.min(); sks = [skew(s) for s in sw]
        c.plot(np.array(cur) - rheo, sks, "o-", ms=5, lw=1.2,
               color=C_I if lab == "Type-I" else C_II, label=lab)
    c.axhline(0, color="#94a3b8", lw=1.0, ls=":")
    c.text(0.5, 0.15, "first-passage regime (atlas)", fontsize=8, color="#16a34a")
    c.text(0.5, -1.4, "regular firing (out-of-atlas)", fontsize=8, color="#b91c1c")
    c.set_xlabel("I - rheobase"); c.set_ylabel("ISI skewness")
    c.set_title("(C) Why: ML exits the noise-driven regime with drive")
    # dedupe legend
    hands, labs = c.get_legend_handles_labels(); seen = dict(zip(labs, hands))
    c.legend(seen.values(), seen.keys(), fontsize=8.6, framealpha=0.95, edgecolor="#cbd5e1")

    fig.tight_layout()
    out = os.path.abspath(os.path.join(HERE, "..", "figures", "morris_lecar_check.png"))
    fig.savefig(out, dpi=140); print("\nsaved", out)


if __name__ == "__main__":
    make_figure()
