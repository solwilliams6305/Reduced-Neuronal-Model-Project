"""mpr_finite_size_figure.py — the finite-size Langevin validation figure."""
from __future__ import annotations
import os
import numpy as np
import mpr_finite_size as mfs


def make_figure():
    if not os.path.exists(mfs.SCALE_CACHE):
        mfs.scaling()
    d = np.load(mfs.SCALE_CACHE)
    Ns = d["Ns"]; Var = d["Var"]; slope = float(d["slope"][0])
    diag = mfs.diag()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    fig, ax = plt.subplots(1, 2, figsize=(12.6, 5.0))

    # ---- A: Var(r) ~ 1/N  =>  finite-size noise is O(N^{-1/2}) ----------------------------------
    a = ax[0]
    a.loglog(Ns, Var, "o", color="#2563eb", ms=10, mec="white", mew=1.2, label="network Var$(r)$")
    ref = Var[0] * (Ns[0] / Ns)                                   # 1/N reference through first point
    a.loglog(Ns, ref, "--", color="#dc2626", lw=1.6, label=r"$\propto 1/N$ (Langevin)")
    a.set_xlabel("network size $N$"); a.set_ylabel("macroscopic rate variance Var$(r)$")
    a.set_title(fr"(A) Finite-size noise is $O(N^{{-1/2}})$: slope ${slope:.2f}$")
    a.legend(fontsize=9.5, framealpha=0.95, edgecolor="#cbd5e1")
    a.grid(True, which="both", ls=":", alpha=0.4)

    # ---- B: the macroscopic bifurcation is a FOLD OF CYCLES (amplitude edge), not a phase edge ---
    b = ax[1]
    es = np.array([x[0] for x in diag]); per = np.array([x[1] for x in diag]); amp = np.array([x[2] for x in diag])
    okp = np.isfinite(per)
    b.plot(es[okp], per[okp], "o-", color="#2563eb", ms=7, lw=1.7, label="period (finite $\\sim$17)")
    b.set_xlabel(r"$\bar\eta$ (mean drive)"); b.set_ylabel("macroscopic period", color="#2563eb")
    b.tick_params(axis="y", labelcolor="#2563eb"); b.invert_xaxis(); b.set_ylim(0, 25)
    b2 = b.twinx()
    b2.plot(es, amp, "s--", color="#16a34a", ms=6, lw=1.5, label="amplitude (steady $\\sim$2.1, then 0)")
    b2.set_ylabel("oscillation amplitude", color="#16a34a"); b2.tick_params(axis="y", labelcolor="#16a34a")
    b2.set_ylim(0, 2.6)
    # mark the abrupt disappearance
    if (~okp).any():
        b.axvline(es[~okp][0], color="#dc2626", lw=1.2, ls=":")
        b.text(es[~okp][0], 2, "cycle\nvanishes", color="#dc2626", fontsize=8, ha="left")
    b.set_title("(B) Macroscopic edge = FOLD OF CYCLES (amplitude), not SNIC")
    h1, l1 = b.get_legend_handles_labels(); h2, l2 = b2.get_legend_handles_labels()
    b.legend(h1 + h2, l1 + l2, fontsize=8.4, framealpha=0.95, edgecolor="#cbd5e1", loc="center left")

    fig.tight_layout()
    out = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures",
                                       "mpr_finite_size.png"))
    fig.savefig(out, dpi=140); print("saved", out)


if __name__ == "__main__":
    make_figure()
