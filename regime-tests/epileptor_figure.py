"""epileptor_figure.py — render the Phase-3b transfer figure from the cached Epileptor run."""
from __future__ import annotations
import os
import numpy as np
import epileptor_phase3b as ep

HERE = os.path.dirname(os.path.abspath(__file__))


def make_figure():
    if not os.path.exists(ep.CACHE):
        ep.build_sims()
    rt, rx, rz, pk = ep.load_sims()
    ct = ep.complexes(pk); sz = ep.seizures(ct)
    info = ep.analyze(verbose=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    C_X, C_Z, C_NU, C_V = "#0f172a", "#2563eb", "#16a34a", "#dc2626"
    fig, ax = plt.subplots(2, 2, figsize=(13.2, 9.4))

    # ---- A: Epileptor trace — fast discharges x1 and the slow drift z ----------------------------
    a = ax[0, 0]; sel = rt < 8000
    a.plot(rt[sel], rx[sel], color=C_X, lw=0.5)
    a.set_xlabel("time"); a.set_ylabel(r"$x_1$ (fast discharges)", color=C_X)
    a.set_title("(A) Epileptor: autonomous seizures driven by slow $z$")
    a2 = a.twinx(); a2.plot(rt[sel], rz[sel], color=C_Z, lw=1.8)
    a2.set_ylabel(r"$z$ (slow permittivity)", color=C_Z); a2.tick_params(axis="y", labelcolor=C_Z)

    # ---- B: ISI structure — raw spike train (spike-and-wave, multimodal) vs complex level --------
    b = ax[0, 1]
    raw = np.diff(pk); raw = raw[(raw > 0) & (raw < 9)]
    comp = np.concatenate([np.diff(st) for st in sz])
    b.hist(raw, bins=np.linspace(0, 9, 46), density=True, color="#94a3b8", alpha=0.7,
           label=fr"raw discharge ISIs ({info['frac_short']*100:.0f}% < 1: spike-and-wave)")
    b.hist(comp, bins=np.linspace(0, 9, 46), density=True, histtype="step", lw=2.2, color=C_NU,
           label=fr"complex (dominant-rhythm) ISIs: Type-I, $W_1{{=}}${info['resid']:.3f}")
    b.set_xlabel("inter-discharge interval"); b.set_ylabel("density")
    b.set_title("(B) Discharge train is spike-and-wave; clean only at complex level")
    b.legend(fontsize=8.0, framealpha=0.95, edgecolor="#cbd5e1")

    # ---- C: the negative result — nu_hat stays high & flat, never warns -------------------------
    c = ax[1, 0]; np_ = info["nu_pos"]
    c.scatter(np_[:, 0], np_[:, 1], s=14, color=C_NU, alpha=0.35, edgecolor="none")
    c.plot(info["ctrs"], info["nu_bin"], "o-", color=C_NU, ms=8, lw=2.0, label=r"binned $\hat\nu$")
    c.axhline(0.4, color="#0f172a", lw=1.2, ls="--", label=r"$\hat\nu$ alarm (Phase 3a)")
    c.axhline(0.0, color="#94a3b8", lw=1.0, ls=":")
    c.annotate("Phase 3a: $\\hat\\nu\\rightarrow0$ here\n(adiabatic SNIC)", xy=(0.5, 0.2),
               xytext=(0.32, 1.1), fontsize=8.4, color="#0f172a",
               arrowprops=dict(arrowstyle="->", color="#0f172a", lw=0.9))
    c.set_xlabel("within-seizure position (0 = onset, 1 = offset)"); c.set_ylabel(r"$\hat\nu$")
    c.set_ylim(-0.2, 3.4); c.set_title(r"(C) $\hat\nu$ stays high & flat: no $\nu\rightarrow0$ warning")
    c.legend(fontsize=8.4, framealpha=0.95, edgecolor="#cbd5e1", loc="center right")

    # ---- D: why — the offset accelerates, it does NOT period-diverge (no precursor) --------------
    d = ax[1, 1]
    L = min(len(np.diff(st)) for st in sz)
    aligned = np.array([np.diff(st)[-L:] / np.median(np.diff(st)) for st in sz])  # align to offset, normalise
    xrel = np.arange(-L, 0)
    md = aligned.mean(0); sd = aligned.std(0)
    d.fill_between(xrel, md - sd, md + sd, color=C_V, alpha=0.18)
    d.plot(xrel, md, "o-", color=C_V, ms=4, lw=1.8, label="Epileptor inter-complex period")
    d.axhline(1.0, color="#94a3b8", lw=1.0, ls=":")
    hyp = 1.0 - 0.9 * np.log((np.abs(xrel) + 0.6) / (L + 0.6))      # illustrative homoclinic divergence
    d.plot(xrel, hyp, "--", color="#0f172a", lw=1.6, label="homoclinic would DIVERGE")
    d.set_xlabel("complexes before offset"); d.set_ylabel("period / median")
    d.set_title(fr"(D) Offset accelerates (ratio {info['ratio']:.2f}), no diverging precursor")
    d.legend(fontsize=8.4, framealpha=0.95, edgecolor="#cbd5e1", loc="upper left")
    d.set_ylim(0.4, 2.6)

    fig.tight_layout()
    out = os.path.abspath(os.path.join(HERE, "..", "figures", "epileptor_phase3b.png"))
    fig.savefig(out, dpi=140); print("\nsaved", out)


if __name__ == "__main__":
    make_figure()
