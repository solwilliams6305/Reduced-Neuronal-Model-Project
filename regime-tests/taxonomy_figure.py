"""taxonomy_figure.py — render the bifurcation-taxonomy early-warning map."""
from __future__ import annotations
import os
import numpy as np
import taxonomy_sweep as ts

HERE = os.path.dirname(os.path.abspath(__file__))


def make_figure():
    if not os.path.exists(ts.CACHE):
        ts.build_sims()
    S = ts.analyze(verbose=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    fig, ax = plt.subplots(2, 2, figsize=(13.4, 9.6))
    C_SNIC, C_SH, C_H, C_F = "#dc2626", "# ea580c".replace(" ", ""), "#2563eb", "#7c3aed"

    # ---- A: PHASE edges — period diverges, nu_hat falls (timing precursor) -----------------------
    a = ax[0, 0]
    for name, col in [("SNIC", C_SNIC), ("Saddle-homoclinic", C_SH)]:
        r = S[name]["rows"]; per = r[:, 2] / np.median(r[:5, 2])
        a.plot(r[:, 0], per, "o-", color=col, ms=4, lw=1.6, label=f"{name} period")
    a.axhline(1.0, color="#94a3b8", lw=1.0, ls=":")
    a.set_xlabel("position (0 = start, 1 = end of drift)"); a.set_ylabel("period / initial period")
    a.set_title("(A) PHASE edge: period DIVERGES toward the bifurcation")
    a2 = a.twinx()
    rs = S["SNIC"]["rows"]; a2.plot(rs[:, 0], rs[:, 1], "s--", color="#0f172a", ms=4, lw=1.4,
                                    label=r"SNIC $\hat\nu$")
    a2.axhline(ts.NU_ALARM, color="#16a34a", lw=1.0, ls="--"); a2.set_ylabel(r"SNIC $\hat\nu$")
    a2.set_ylim(-0.5, 3.2)
    h1, l1 = a.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels()
    a.legend(h1 + h2, l1 + l2, fontsize=8.2, framealpha=0.95, edgecolor="#cbd5e1", loc="upper left")

    # ---- B: AMPLITUDE edges — amplitude collapses (Hopf gradual, fold abrupt) ---------------------
    b = ax[0, 1]
    for name, col in [("Supercrit Hopf", C_H), ("Fold-of-cycles", C_F)]:
        s = S[name]; at = s["at"]; av = s["av"]; tt = at[-1]
        ae = np.mean(av[at < 0.35 * s["frac"] * tt])
        b.plot(at / tt, av / ae, color=col, lw=1.8, label=name)
        b.axvline(s["frac"], color=col, lw=1.0, ls=":")
    b.axhline(1.0, color="#94a3b8", lw=1.0, ls=":")
    b.set_xlabel("position (0 = start, 1 = end of drift)"); b.set_ylabel("amplitude / initial amplitude")
    b.set_title("(B) AMPLITUDE edge: Hopf declines (precursor) · fold JUMPS (none)")
    b.legend(fontsize=8.6, framealpha=0.95, edgecolor="#cbd5e1", loc="lower left")

    # ---- C: the early-warning MAP ----------------------------------------------------------------
    c = ax[1, 0]
    order = ["SNIC", "Saddle-homoclinic", "Supercrit Hopf", "Fold-of-cycles"]
    cols = ["period\ndiverges", r"timing $\hat\nu$" + "\nwarns", "amplitude\nprecursor", r"amplitude$\to$0"]
    G = 0.5  # gray = n/a
    M = np.array([
        [1, 1, G, 0],     # SNIC
        [1, 1, G, 0],     # SH
        [0, 0, 1, 1],     # Hopf
        [0, 0, 0, 1],     # fold
    ], float)
    cmap = plt.cm.RdYlGn
    c.imshow(M, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    c.set_xticks(range(4)); c.set_xticklabels(cols, fontsize=8.6)
    c.set_yticks(range(4)); c.set_yticklabels(order, fontsize=9)
    lab = {1: "yes", 0: "no", G: "n/a"}
    for i in range(4):
        for j in range(4):
            c.text(j, i, lab[M[i, j]], ha="center", va="center", fontsize=9.5,
                   color="#0f172a", fontweight="bold")
    c.axhline(1.5, color="#0f172a", lw=2.0)
    c.text(3.6, 0.5, "PHASE\nedge", rotation=90, va="center", fontsize=8.5, color="#b91c1c")
    c.text(3.6, 2.5, "AMPLITUDE\nedge", rotation=90, va="center", fontsize=8.5, color="#1d4ed8")
    c.set_title("(C) Early-warning splits along the two-edge line")

    # ---- D: tie to the two-edge theory -----------------------------------------------------------
    d = ax[1, 1]; d.axis("off")
    d.text(0.5, 0.97, "The early-warning observable IS the two-channel split", ha="center",
           va="top", fontsize=11.5, fontweight="bold", color="#0f172a")
    txt = (
        r"$\bf{PHASE\ edge}$  (Type-I, $\sigma^{2/3}$):  SNIC, saddle-homoclinic" "\n"
        "   the firing PERIOD diverges as the parameter -> threshold;\n"
        r"   the timing inversion $\hat\nu\to0$ gives the early warning" "\n"
        "   (the noisy-folded-cycle phase edge / quartic FPT channel).\n\n"
        r"$\bf{AMPLITUDE\ edge}$  (Type-II, $\sigma^{1/2}$):  Hopf, fold-of-cycles" "\n"
        "   the period stays finite — the timing inversion is BLIND;\n"
        "   the oscillation AMPLITUDE is the early-warning observable\n"
        "   (the canard peel-off / Tracy-Widom level channel).\n\n"
        r"$\bf{Catastrophe\ caveat}$:  the fold-of-cycles holds amplitude" "\n"
        "   then JUMPS — no precursor in either channel: the genuinely\n"
        "   warning-free tipping. Knowing the bifurcation CLASS tells you\n"
        "   which observable (if any) can warn — the value of the inversion."
    )
    d.text(0.02, 0.86, txt, ha="left", va="top", fontsize=9.4, color="#1e293b", family="monospace")

    fig.tight_layout()
    out = os.path.abspath(os.path.join(HERE, "..", "figures", "taxonomy_sweep.png"))
    fig.savefig(out, dpi=140); print("\nsaved", out)


if __name__ == "__main__":
    make_figure()
