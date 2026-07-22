"""
allen_phase2_inversion.py — MVP Phase 2: run the Phase-1 inverter on per-sweep spike trains.

Two modes share ONE per-sweep routine `invert_sweep` (class via classify3 + confound-corrected nu via
fit_adapt, from bifurcation_classifier):

  dryrun :  Allen-LIKE surrogate cells (varying rheobase/scale/noise/adaptation), each swept across
            currents, with REALISTIC small per-sweep ISI counts (~25) — the stress test that decides
            real-data feasibility.  Runs in-sandbox.
  allen  :  the real thing — pull Long-Square sweeps from the Allen Cell Types DB, invert each.
            RUN LOCALLY (`pip install allensdk`, network).  Same code path as dryrun.

Phase-2 success criteria:
  H1  nu_hat rises MONOTONICALLY with injected current (internal consistency of the inversion);
  H2  the cell classifies Type-I (SNIC/homoclinic, not Hopf) for regular-spiking cells;
  H3  confound-corrected nu_hat agrees across cells at matched (I - rheobase) despite different scale.

    python3 allen_phase2_inversion.py dryrun     # -> figures/allen_phase2_dryrun.png
    python3 allen_phase2_inversion.py allen      # LOCAL: -> figures/allen_phase2_inversion.png
"""
from __future__ import annotations
import os
import sys
import numpy as np
from bifurcation_classifier import classify3, fit_adapt
from offcritical_phase_edge import J_of_nu

HERE = os.path.dirname(os.path.abspath(__file__))


# ---- the shared per-sweep inversion --------------------------------------------------------------
def invert_sweep(isis):
    """ISIs of one sweep -> (class, nu_naive, nu_corrected, b_hat, residual)."""
    isis = np.asarray(isis, float); isis = isis[isis > 0]
    if len(isis) < 8:
        return None
    cls, p, w, per = classify3(isis)
    nu_n, _, _ = fit_adapt(isis, joint=False)
    nu_c, b_c, _ = fit_adapt(isis, joint=True)
    return cls, nu_n, nu_c, b_c, w


def boot_nu(isis, nboot=40, seed=0):
    rng = np.random.default_rng(seed); isis = np.asarray(isis, float)
    out = []
    for _ in range(nboot):
        s = rng.choice(isis, len(isis), replace=True)
        r = fit_adapt(s, joint=True)
        out.append(r[0])
    return np.percentile(out, [16, 84])


# ---- surrogate Long-Square sweeps (vectorised across all cell x current units) --------------------
def sweep_isis_batch(nus, sigmas, bs, tau_w=6.0, n_isi=26, dt=2.5e-3,
                     v_th=14.0, v_reset=-14.0, seed=0):
    M = len(nus); nus = np.asarray(nus, float); sigmas = np.asarray(sigmas, float); bs = np.asarray(bs, float)
    I = nus * sigmas ** (4.0 / 3.0); rng = np.random.default_rng(seed)
    v = np.full(M, v_reset, float); w = np.zeros(M); last = np.full(M, np.nan)
    cnt = np.zeros(M, int); isis = [[] for _ in range(M)]; sdt = np.sqrt(dt)
    means = np.array([sigmas[i] ** (-2 / 3) * max(J_of_nu(nus[i]), 2.0) for i in range(M)])
    Tmax = float((n_isi + 3) * means.max())
    for k in range(int(Tmax / dt)):
        t = k * dt
        v += (v * v + I - w) * dt + sigmas * sdt * rng.standard_normal(M)
        w -= (w / tau_w) * dt
        fired = np.where(v > v_th)[0]
        for i in fired:
            if cnt[i] >= 2 and np.isfinite(last[i]) and len(isis[i]) < n_isi:
                isis[i].append(t - last[i])
            last[i] = t; cnt[i] += 1; w[i] += bs[i]; v[i] = v_reset
    return [np.array(x) for x in isis]


def dryrun():
    # four surrogate Type-I cells: different noise & adaptation (i.e. different "scale")
    cells = [("cell A  σ=.3 b=0",   0.30, 0.00, "#dc2626"),
             ("cell B  σ=.4 b=.1",  0.40, 0.10, "#ea580c"),
             ("cell C  σ=.5 b=.2",  0.50, 0.20, "#2563eb"),
             ("cell D  σ=.6 b=.3",  0.60, 0.30, "#7c3aed")]
    true_nu = np.array([0.3, 0.7, 1.1, 1.6, 2.2])           # increasing injected current
    NU, SG, BS = [], [], []
    for _, sg, b, _c in cells:
        for nu in true_nu:
            NU.append(nu); SG.append(sg); BS.append(b)
    sweeps = sweep_isis_batch(NU, SG, BS, seed=7)

    rows = []; idx = 0
    for ci, (lab, sg, b, col) in enumerate(cells):
        for nu in true_nu:
            isis = sweeps[idx]; idx += 1
            res = invert_sweep(isis)
            if res is None:
                rows.append((ci, nu, np.nan, np.nan, "none", np.nan, np.nan, len(isis))); continue
            cls, nn, nc, bh, w = res
            lo, hi = boot_nu(isis, seed=idx)
            rows.append((ci, nu, nn, nc, cls, lo, hi, len(isis)))
    R = rows
    # H1: monotonicity of corrected nu_hat vs current, per cell
    print("  H1 monotonicity (corrected ν̂ increasing with current):")
    for ci, (lab, sg, b, col) in enumerate(cells):
        nc = [r[3] for r in R if r[0] == ci]
        mono = all(nc[i] <= nc[i + 1] + 1e-9 for i in range(len(nc) - 1))
        print(f"    {lab}: ν̂={[round(x,2) for x in nc]}  monotone={mono}  (median {int(np.median([r[7] for r in R if r[0]==ci]))} ISIs/sweep)")
    classes = [r[4] for r in R]
    typeI = np.mean([c in ("SNIC", "Homoclinic") for c in classes]) * 100
    print(f"  H2 class: {typeI:.0f}% of sweeps called Type-I (SNIC/Homoclinic), "
          f"{np.mean([c=='Hopf' for c in classes])*100:.0f}% Hopf")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    fig, ax = plt.subplots(1, 3, figsize=(15.8, 4.8))

    # Panel A: corrected nu_hat vs true nu (current), per cell -> tracking + cross-cell consistency
    ax[0].plot([0, 2.4], [0, 2.4], "--", color="#94a3b8", lw=1.2, label="ideal")
    for ci, (lab, sg, b, col) in enumerate(cells):
        sub = [r for r in R if r[0] == ci and np.isfinite(r[3])]
        x = [r[1] for r in sub]; y = [r[3] for r in sub]
        yerr = np.abs(np.array([[r[3]-r[5] for r in sub], [r[6]-r[3] for r in sub]]))
        ax[0].errorbar(x, y, yerr=yerr, fmt="o", color=col, ms=6, mec="white", mew=0.8,
                       capsize=2.5, lw=0.9, label=lab)
    ax[0].set_xlabel(r"true $\nu$ (injected current)"); ax[0].set_ylabel(r"corrected $\hat\nu$")
    ax[0].set_title("H1/H3: ν̂ tracks current across cells")
    ax[0].legend(fontsize=7.6, framealpha=0.95, edgecolor="#cbd5e1")

    # Panel B: naive vs corrected (pooled) -> confound correction on realistic small samples
    nn = np.array([r[2] for r in R]); nc = np.array([r[3] for r in R]); tn = np.array([r[1] for r in R])
    ax[1].plot([0, 2.4], [0, 2.4], "--", color="#94a3b8", lw=1.2, label="ideal")
    ax[1].plot(tn, nn, "o", color="#dc2626", ms=6, mec="white", mew=0.8, alpha=0.7, label="naive (b=0)")
    ax[1].plot(tn, nc, "s", color="#16a34a", ms=6, mec="white", mew=0.8, alpha=0.8, label="corrected (ν,b)")
    ax[1].set_xlabel(r"true $\nu$"); ax[1].set_ylabel(r"$\hat\nu$")
    ax[1].set_title("Confound correction (small samples)")
    ax[1].legend(fontsize=8.4, framealpha=0.95, edgecolor="#cbd5e1")

    # Panel C: class calls
    from collections import Counter
    cnt = Counter(classes); order = ["SNIC", "Homoclinic", "Hopf", "none"]
    vals = [cnt.get(k, 0) for k in order]
    ax[2].bar(order, vals, color=["#dc2626", "#16a34a", "#2563eb", "#94a3b8"])
    ax[2].set_ylabel("# sweeps"); ax[2].set_title(f"H2: class calls ({typeI:.0f}% Type-I)")
    for i, v in enumerate(vals):
        ax[2].text(i, v + 0.2, str(v), ha="center", fontsize=10)

    fig.tight_layout()
    out = os.path.abspath(os.path.join(HERE, "..", "figures", "allen_phase2_dryrun.png"))
    fig.savefig(out, dpi=140); print("saved", out)


def allen():
    """REAL Allen Cell Types inversion (run locally). Mirrors allen_qif_validation's pull, then
    inverts each Long-Square sweep with invert_sweep."""
    from allensdk.core.cell_types_cache import CellTypesCache
    ctc = CellTypesCache(manifest_file=os.path.join(HERE, "cell_types", "manifest.json"))
    cells = [c for c in ctc.get_cells(require_morphology=False)
             if c.get("dendrite_type") == "spiny"][:40]
    per_cell = []
    for c in cells:
        cid = c["id"]
        try:
            ds = ctc.get_ephys_data(cid); sw = ctc.get_ephys_sweeps(cid)
        except Exception:
            continue
        ls = [s for s in sw if s.get("stimulus_name") == "Long Square"
              and s.get("stimulus_absolute_amplitude") is not None and (s.get("num_spikes") or 0) >= 6]
        if len(ls) < 3:
            continue
        rheo = min(s["stimulus_absolute_amplitude"] for s in ls)
        pts = []
        for s in sorted(ls, key=lambda z: z["stimulus_absolute_amplitude"]):
            try:
                st = np.asarray(ds.get_spike_times(s["sweep_number"]), float)
            except Exception:
                continue
            isis = np.diff(st); isis = isis[isis > 0]
            r = invert_sweep(isis)
            if r:
                pts.append((s["stimulus_absolute_amplitude"] - rheo, r[2], r[0]))  # (I-rheo, nu_corr, class)
        if len(pts) >= 3:
            per_cell.append((cid, np.array([p[:2] for p in pts], float), [p[2] for p in pts]))
    # H1: fraction of cells with monotone nu_hat vs current; H2: Type-I fraction
    mono = np.mean([np.all(np.diff(p[:, 1]) >= -1e-6) for _, p, _ in per_cell]) if per_cell else np.nan
    allcls = [c for _, _, cs in per_cell for c in cs]
    typeI = np.mean([c in ("SNIC", "Homoclinic") for c in allcls]) * 100 if allcls else np.nan
    print(f"  {len(per_cell)} cells; H1 monotone fraction={mono:.2f}; H2 Type-I={typeI:.0f}%")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11.2, 4.8))
    for cid, p, cs in per_cell:
        ax[0].plot(p[:, 0], p[:, 1], "o-", ms=4, lw=0.8, alpha=0.6)
    ax[0].set_xlabel("I - rheobase (pA)"); ax[0].set_ylabel(r"corrected $\hat\nu$")
    ax[0].set_title(f"H1: ν̂ vs current ({len(per_cell)} Allen cells)")
    from collections import Counter
    cnt = Counter(allcls); order = ["SNIC", "Homoclinic", "Hopf"]
    ax[1].bar(order, [cnt.get(k, 0) for k in order], color=["#dc2626", "#16a34a", "#2563eb"])
    ax[1].set_title(f"H2: {typeI:.0f}% Type-I"); ax[1].set_ylabel("# sweeps")
    fig.tight_layout()
    out = os.path.abspath(os.path.join(HERE, "..", "figures", "allen_phase2_inversion.png"))
    fig.savefig(out, dpi=140); print("saved", out)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) >= 2 else ""
    if cmd == "dryrun":
        dryrun()
    elif cmd == "allen":
        allen()
    else:
        print("usage: allen_phase2_inversion.py dryrun | allen")
