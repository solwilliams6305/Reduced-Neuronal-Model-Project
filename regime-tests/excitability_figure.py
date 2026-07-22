"""excitability_figure.py — render the excitability-characterizer validation figure."""
from __future__ import annotations
import os
import numpy as np
import excitability_characterizer as ec
import bifurcation_classifier as bc

HERE = os.path.dirname(os.path.abspath(__file__))
COL = {"Type-I": "#dc2626", "Type-II": "#2563eb", "out-of-family": "#64748b"}


def make_figure():
    if not os.path.exists(ec.CACHE):
        ec.build_cohort()
    reps, conf, acc, nu_err = ec.validate(verbose=True)
    cells = ec.load_cohort()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    fig, ax = plt.subplots(2, 2, figsize=(13.2, 9.4))

    # ---- A: the decision space — f-I onset vs goodness-of-fit residual --------------------------
    a = ax[0, 0]
    for lab in ["Type-I", "Type-II", "out-of-family"]:
        xs = [r["fI_onset"] for l, r in reps if l == lab and "fI_onset" in r]
        ys = [r["residual"] for l, r in reps if l == lab and "fI_onset" in r]
        a.scatter(xs, ys, c=COL[lab], s=70, edgecolor="white", lw=1.0, label=lab)
    a.axhline(0.07, color="#0f172a", lw=1.2, ls="--"); a.text(0.42, 0.082, "reject (out-of-family)", fontsize=8)
    a.axvline(0.6, color="#0f172a", lw=1.0, ls=":"); a.text(0.61, 0.0, "Hodgkin\nClass 1|2", fontsize=8)
    a.set_xlabel("f-I onset  (min/max sweep rate)"); a.set_ylabel(r"goodness-of-fit residual ($W_1$)")
    a.set_title("(A) Three features separate the classes")
    a.set_yscale("log"); a.legend(fontsize=8.6, framealpha=0.95, edgecolor="#cbd5e1", loc="center right")

    # ---- B: the ISI-shape signatures -------------------------------------------------------------
    b = ax[0, 1]; bins = np.linspace(0, 3.2, 44)
    exemplars = [("Type-I (quartic-FPT, skewed)", bc.snic_intervals(0.6, N=3000, seed=3), COL["Type-I"]),
                 ("Type-II (Hopf jitter, ~symmetric)", bc.hopf_intervals(0.4, seed=3), COL["Type-II"]),
                 ("out-of-family (Poisson)", np.random.default_rng(1).exponential(1.0, 3000), "#16a34a")]
    for nm, s, c in exemplars:
        s = np.asarray(s); s = s[s > 0]
        b.hist(s / s.mean(), bins=bins, density=True, histtype="step", lw=2.2, color=c, label=nm)
    b.set_xlabel(r"ISI / $\langle$ISI$\rangle$"); b.set_ylabel("density")
    b.set_title("(B) Excitability class is written in the ISI shape")
    b.legend(fontsize=8.2, framealpha=0.95, edgecolor="#cbd5e1")

    # ---- C: confusion matrix ---------------------------------------------------------------------
    c = ax[1, 0]; labels = ["Type-I", "Type-II", "out-of-fam"]
    c.imshow(conf, cmap="Blues", vmin=0)
    c.set_xticks(range(3)); c.set_yticks(range(3))
    c.set_xticklabels(labels, fontsize=9); c.set_yticklabels(labels, fontsize=9)
    for i in range(3):
        for j in range(3):
            c.text(j, i, conf[i, j], ha="center", va="center", fontsize=13,
                   color="white" if conf[i, j] > conf.max() / 2 else "#0f172a")
    c.set_xlabel("predicted"); c.set_ylabel("true")
    c.set_title(f"(C) Class recovery from spikes: {acc*100:.0f}% (n={conf.sum()})")

    # ---- D: f-I onset curves (Hodgkin Class 1 vs 2) ----------------------------------------------
    d = ax[1, 1]

    def norm_rate(cell):
        lab, sw, cur = cell
        rates = []
        for s in sw:
            s = np.asarray(s); s = s[s > 0]
            if len(s) >= 8:
                rates.append(1.0 / s.mean())
        rates = np.array(rates); return rates / rates.max()
    tI = next(c for c in cells if c[0] == "Type-I"); tII = next(c for c in cells if c[0] == "Type-II")
    d.plot(np.linspace(0, 1, len(norm_rate(tI))), norm_rate(tI), "o-", color=COL["Type-I"], ms=7,
           label="Type-I: rate $\\to$ 0 at onset (Class 1)")
    d.plot(np.linspace(0, 1, len(norm_rate(tII))), norm_rate(tII), "s-", color=COL["Type-II"], ms=7,
           label="Type-II: bounded onset rate (Class 2)")
    d.set_xlabel("normalised drive (rheobase $\\to$ max)"); d.set_ylabel("normalised firing rate")
    d.set_title(r"(D) f-I onset distinguishes Hodgkin classes")
    d.set_ylim(0, 1.05); d.legend(fontsize=8.6, framealpha=0.95, edgecolor="#cbd5e1", loc="upper left")

    fig.tight_layout()
    out = os.path.abspath(os.path.join(HERE, "..", "figures", "excitability_characterizer.png"))
    fig.savefig(out, dpi=140); print("\nsaved", out)


if __name__ == "__main__":
    make_figure()
