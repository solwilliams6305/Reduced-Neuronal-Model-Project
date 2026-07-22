"""
allen_excitability.py — REAL Allen Cell Types excitability characterization (RUN LOCALLY).

Same code path as the synthetic cohort validation in excitability_characterizer.py: for each cell,
collect its Long-Square sweeps (one ISI array per drive level), call `characterize_cell`, and tabulate
(excitability class, proximity-to-threshold nu_hat, goodness-of-fit residual, accept/reject) against
the cell's metadata. The synthetic cohort (94% three-way accuracy) is the dry run; this is the test on
data.

This does a network pull and so is NOT executed in-sandbox (policy). Run on your machine:

    pip install allensdk --break-system-packages
    python3 allen_excitability.py            # -> figures/allen_excitability.png + prints a table

Notes:
  * Allen has no direct "Type-I/II" label; the honest checks are (a) internal consistency — does nu_hat
    rise monotonically with injected current? (b) does the f-I onset (Hodgkin class) agree with the
    cell's known f-I curve / spiny-vs-aspiny + transgenic line? (c) are non-spiking/garbage sweeps
    rejected? Treat the class call as a HYPOTHESIS to be checked against the biophysics, not a label.
"""
from __future__ import annotations
import os
import sys
import numpy as np
from excitability_characterizer import characterize_cell

HERE = os.path.dirname(os.path.abspath(__file__))
N_CELLS = 60


def run():
    from allensdk.core.cell_types_cache import CellTypesCache
    ctc = CellTypesCache(manifest_file=os.path.join(HERE, "cell_types", "manifest.json"))
    cells = ctc.get_cells(require_morphology=False)[:N_CELLS]
    reports = []
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
        ls = sorted(ls, key=lambda z: z["stimulus_absolute_amplitude"])
        rheo = ls[0]["stimulus_absolute_amplitude"]
        sweeps, currents = [], []
        for s in ls:
            try:
                st = np.asarray(ds.get_spike_times(s["sweep_number"]), float)
            except Exception:
                continue
            isis = np.diff(st); isis = isis[isis > 0]
            if len(isis) >= 8:
                sweeps.append(isis); currents.append(s["stimulus_absolute_amplitude"] - rheo)
        if len(sweeps) < 2:
            continue
        rep = characterize_cell(sweeps, currents)
        rep["cid"] = cid; rep["dendrite_type"] = c.get("dendrite_type")
        rep["line"] = c.get("transgenic_line") or c.get("line_name")
        reports.append(rep)

    usable = [r for r in reports if "klass" in r and r["klass"] != "uncharacterizable"]
    from collections import Counter
    cnt = Counter(r["klass"] for r in usable)
    print(f"  {len(usable)} cells characterized: {dict(cnt)}")
    rej = np.mean([r.get("out_of_family", False) for r in usable]) * 100
    print(f"  out-of-family rejected: {rej:.0f}%")
    # internal-consistency check: nu_hat monotone vs current (excludes rejected)
    mono = []
    for r in usable:
        if r.get("out_of_family"):
            continue
        nu = np.array(r["nu_hat"])
        if len(nu) >= 3:
            mono.append(np.corrcoef(np.arange(len(nu)), nu)[0, 1])
    print(f"  median nu_hat-vs-drive correlation (should be > 0): {np.median(mono):+.2f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11})
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
    COL = {"Type-I (SNIC/homoclinic)": "#dc2626", "Type-II (Hopf)": "#2563eb", "out-of-family": "#64748b"}
    for r in usable:
        ax[0].scatter(r["fI_onset"], r["residual"], c=COL.get(r["klass"], "#999"), s=45,
                      edgecolor="white", lw=0.6)
    ax[0].axhline(0.07, color="#0f172a", ls="--", lw=1.0)
    ax[0].set_xlabel("f-I onset"); ax[0].set_ylabel("residual (W1)"); ax[0].set_yscale("log")
    ax[0].set_title("Decision space (real Allen cells)")
    order = list(COL); ax[1].bar(range(len(order)), [cnt.get(k, 0) for k in order],
                                 color=[COL[k] for k in order])
    ax[1].set_xticks(range(len(order))); ax[1].set_xticklabels([k.split(" ")[0] for k in order], fontsize=9)
    ax[1].set_ylabel("# cells"); ax[1].set_title("Class calls")
    for r in usable:
        if not r.get("out_of_family") and r.get("currents"):
            ax[2].plot(r["currents"], r["nu_hat"], "o-", ms=3, lw=0.7, alpha=0.5,
                       color=COL.get(r["klass"], "#999"))
    ax[2].set_xlabel("I - rheobase (pA)"); ax[2].set_ylabel(r"$\hat\nu$ (proximity)")
    ax[2].set_title("Proximity-to-threshold vs drive")
    fig.tight_layout()
    out = os.path.abspath(os.path.join(HERE, "..", "figures", "allen_excitability.png"))
    fig.savefig(out, dpi=140); print("  saved", out)


if __name__ == "__main__":
    run()
