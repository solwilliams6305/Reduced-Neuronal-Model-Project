"""neural_mass_figure.py — the collective two-edge figure (theory-track opener)."""
from __future__ import annotations
import os
import numpy as np
import neural_mass_edge as nm
import bifurcation_classifier as bc
from bifurcation_classifier import skew


def make_figure():
    if not os.path.exists(nm.CACHE):
        nm.build_sims()
    d = np.load(nm.CACHE)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    C_R, C_AT = "#7c3aed", "#16a34a"
    fig, ax = plt.subplots(1, 3, figsize=(16.2, 4.9))

    # ---- A: collective bursting — the macroscopic limit cycle (firing rate r) -------------------
    a = ax[0]; dt = float(d["trace_dt"][0]); rtr = d["trace_r"]; t = np.arange(len(rtr)) * dt
    a.plot(t, rtr, color=C_R, lw=0.9)
    a.set_xlabel("time"); a.set_ylabel("collective firing rate $r(t)$")
    a.set_title("(A) Network bursts: macroscopic limit cycle")

    # ---- B: macroscopic period diverges toward the network SNIC/homoclinic ----------------------
    b = ax[1]
    es = d["eta_scan"]; pers = d["period"]; ok = np.isfinite(pers)
    b.plot(es[ok], pers[ok], "o-", color="#2563eb", ms=8, lw=1.8)
    b.axvline(nm.ETA_STAR, color="#dc2626", lw=1.3, ls="--")
    b.text(nm.ETA_STAR - 0.002, pers[ok].min(), r"$\bar\eta^*$ (phase edge)", color="#dc2626",
           fontsize=9, rotation=90, va="bottom", ha="right")
    b.set_xlabel(r"$\bar\eta$ (mean drive)"); b.set_ylabel("macroscopic burst period")
    b.set_title("(B) Period DIVERGES: network-level phase edge")
    b.invert_xaxis()

    # ---- C: collective inter-burst intervals = the single-neuron phase-edge shape ---------------
    c = ax[2]; bins = np.linspace(0, 5, 44)
    ibi = d["ibi_2"]; cls, p, w, per = __import__("harden_inverter").classify3_interp(ibi)
    c.hist(ibi / ibi.mean(), bins=bins, density=True, histtype="stepfilled", alpha=0.45, color=C_R,
           label=fr"collective IBIs (skew {skew(ibi):+.1f}, class {cls})")
    sn = bc.snic_intervals(-0.5, N=4000, seed=3)                       # single-neuron phase-edge (Kramers)
    c.hist(sn / sn.mean(), bins=bins, density=True, histtype="step", lw=2.2, color=C_AT,
           label=fr"single-neuron SNIC law (skew {skew(sn):+.1f})")
    c.set_xlabel(r"interval / $\langle$interval$\rangle$"); c.set_ylabel("density")
    c.set_title("(C) Same phase-edge shape: neuron $\\equiv$ network")
    c.legend(fontsize=8.4, framealpha=0.95, edgecolor="#cbd5e1")

    fig.tight_layout()
    out = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures",
                                       "neural_mass_edge.png"))
    fig.savefig(out, dpi=140); print("saved", out)


if __name__ == "__main__":
    make_figure()
