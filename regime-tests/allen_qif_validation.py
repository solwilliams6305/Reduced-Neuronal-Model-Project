"""
allen_qif_validation.py — REAL-DATA test of the QIF phase-edge predictions on the Allen Cell Types
Database (open patch-clamp recordings).  RUN LOCALLY: needs network access to the Allen API and
`pip install allensdk` (pulls h5py/scipy/pandas).  It is the real-data twin of the synthetic
ground-truth harness `qif_validation.py`, and emits the SAME 3-panel figure
(`figures/allen_qif_validation.png`) so the two can be read side by side.

    pip install allensdk
    python3 allen_qif_validation.py            # pulls ~N_CELLS cells, caches under ./cell_types/

What it tests (see QIF_PHASE_EDGE.md / QIF_DATA_VALIDATION.md):
  P2  CV crossover:   ISI CV vs firing rate, pooled over Type-I cells; near-rheobase CV ~ 0.57?
  P3  ISI shape:      near-rheobase normalized-ISI density vs exponential (Poisson null) and vs the
                      white-QIF QUARTIC-FPT reference (skewed, mode away from 0).
  P1  sigma^{2/3}:    PARTIAL only — the Allen 'Noise' protocol fixes the noise CV (0.2) and scales
                      the mean to 0.75/1/1.5x rheobase, so sigma co-varies with the mean and a clean
                      independent sigma-sweep is NOT available.  A genuine P1 collapse needs
                      dynamic-clamp data with independent injected-noise control; here we only show
                      the f-I noise-smoothing (Noise vs Long-Square onset).  Honest caveat, by design.

Type-I cohort: spiny / regular-spiking excitatory cells are a Type-I (SNIC, class-1 excitability)
proxy.  This is a PROXY — gold-standard Type-I/II separation is by phase-response curve (all-positive
PRC = Type-I), which the database does not tabulate.  Set CRE_LINE / DENDRITE to refine.
"""
from __future__ import annotations
import os
import numpy as np

N_CELLS = 40                      # cap for runtime; raise for a fuller pool
DENDRITE = "spiny"                # Type-I proxy: spiny == excitatory; set None to disable
CRE_LINE = None                   # e.g. "Cux2-CreERT2" to restrict to a line; None = any
LONG_SQUARE = "Long Square"
MIN_SPIKES = 4                    # need >=4 spikes (>=3 ISIs) for a CV estimate
HERE = os.path.dirname(os.path.abspath(__file__))


# ----- estimators (identical definitions to qif_validation.py) -----------------------------------
def skew(x):
    x = np.asarray(x, float); m = x.mean(); sd = x.std()
    return float(np.mean(((x - m) / sd) ** 3)) if sd > 0 and len(x) > 2 else np.nan


def quartic_fpt_reference(n=40000):
    """White-noise QIF threshold ISIs (nu=0) = the quartic FPT law, normalized to unit mean."""
    rng = np.random.default_rng(0)
    N = 4000; dt = 2.5e-3; v = np.full(N, -14.0); last = np.zeros(N); isis = []
    sigma = 0.5
    for k in range(int(160 / dt)):
        t = k * dt
        v += (v * v) * dt + sigma * np.sqrt(dt) * rng.standard_normal(N)
        fired = v > 14.0
        if fired.any():
            idx = np.where(fired)[0]; isis.extend((t - last[idx]).tolist())
            last[idx] = t; v[fired] = -14.0
    isis = np.array(isis); isis = isis[isis > 0]
    return isis / isis.mean()


# ----- Allen pull --------------------------------------------------------------------------------
def main():
    from allensdk.core.cell_types_cache import CellTypesCache
    ctc = CellTypesCache(manifest_file=os.path.join(HERE, "cell_types", "manifest.json"))
    cells = ctc.get_cells(require_morphology=False, require_reconstruction=False)

    def keep(c):
        if DENDRITE and c.get("dendrite_type") != DENDRITE:
            return False
        if CRE_LINE and CRE_LINE not in (c.get("transgenic_line") or ""):
            return False
        return True

    cells = [c for c in cells if keep(c)]
    print(f"{len(cells)} candidate Type-I-proxy cells; using up to {N_CELLS}")

    rate_cv = []          # (rate_Hz, cv, I_minus_rheo_pA) per supra-threshold Long-Square sweep
    near_isis = []        # normalized ISIs from the lowest supra-threshold sweep per cell
    n_used = 0
    for c in cells:
        if n_used >= N_CELLS:
            break
        cid = c["id"]
        try:
            ds = ctc.get_ephys_data(cid)
            sweeps = ctc.get_ephys_sweeps(cid)
        except Exception as e:
            print(f"  skip {cid}: {type(e).__name__}"); continue
        ls = [s for s in sweeps if s.get("stimulus_name") == LONG_SQUARE
              and s.get("stimulus_absolute_amplitude") is not None]
        if not ls:
            continue
        # rheobase = smallest Long-Square amplitude that fires
        fired = [(s["stimulus_absolute_amplitude"], s["sweep_number"]) for s in ls
                 if (s.get("num_spikes") or 0) >= 1]
        if not fired:
            continue
        rheo = min(a for a, _ in fired)
        cell_rows = []
        for amp, swn in sorted(fired):
            try:
                st = ds.get_spike_times(swn)
            except Exception:
                continue
            st = np.asarray(st, float)
            if len(st) < MIN_SPIKES:
                continue
            isi = np.diff(st)
            isi = isi[isi > 0]
            if len(isi) < 3:
                continue
            r = 1.0 / isi.mean(); cv = isi.std() / isi.mean()
            rate_cv.append((r, cv, amp - rheo))
            cell_rows.append((amp, isi / isi.mean()))
        if cell_rows:                       # lowest supra-threshold sweep = near rheobase
            near_isis.append(min(cell_rows, key=lambda t: t[0])[1])
            n_used += 1
            print(f"  [{n_used:2d}] cell {cid}: {len(cell_rows)} f-I points, rheo {rheo:.0f} pA")

    rate_cv = np.array(rate_cv) if rate_cv else np.zeros((0, 3))
    near = np.concatenate(near_isis) if near_isis else np.zeros(0)
    print(f"\n  pooled: {len(rate_cv)} sweeps, {len(near)} near-rheobase ISIs from {n_used} cells")
    if len(near):
        print(f"  near-rheobase ISI:  CV {near.std():.3f}   skew {skew(near):.3f}   "
              f"(predict CV~0.57, skew~1.8; exponential would be CV 1.0 skew 2.0)")

    # ----- figure (same layout as the synthetic harness) -----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    fig, ax = plt.subplots(1, 3, figsize=(15.6, 4.8))

    # P2: CV vs rate (real), with SNIC 0.57 and the binned trend
    if len(rate_cv):
        ax[0].plot(rate_cv[:, 0], rate_cv[:, 1], "o", color="#64748b", ms=4, alpha=0.4,
                   label="per sweep")
        edges = np.linspace(0, np.percentile(rate_cv[:, 0], 95), 9)
        mids = 0.5 * (edges[:-1] + edges[1:]); binned = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            sel = rate_cv[(rate_cv[:, 0] >= lo) & (rate_cv[:, 0] < hi), 1]
            binned.append(sel.mean() if len(sel) else np.nan)
        ax[0].plot(mids, binned, "o-", color="#0f172a", ms=7, mec="white", mew=1.0,
                   label="binned mean")
    ax[0].axhline(0.57, color="#16a34a", lw=1.2, ls="--"); ax[0].text(0.02, 0.60, "SNIC 0.57",
              transform=ax[0].get_yaxis_transform(), fontsize=8.5, color="#16a34a")
    ax[0].set_xlabel("firing rate (Hz)"); ax[0].set_ylabel("ISI coefficient of variation")
    ax[0].set_title("P2: CV vs rate (Type-I cells)")
    ax[0].set_ylim(0, 1.2); ax[0].legend(fontsize=8.6, framealpha=0.95, edgecolor="#cbd5e1")

    # P2b: CV vs (I - rheobase) -- the rescaled-detuning view
    if len(rate_cv):
        ax[1].plot(rate_cv[:, 2], rate_cv[:, 1], "o", color="#64748b", ms=4, alpha=0.4)
        order = np.argsort(rate_cv[:, 2]); xs = rate_cv[order, 2]
        edges = np.linspace(xs.min(), np.percentile(xs, 95), 9)
        mids = 0.5 * (edges[:-1] + edges[1:]); binned = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            sel = rate_cv[(rate_cv[:, 2] >= lo) & (rate_cv[:, 2] < hi), 1]
            binned.append(sel.mean() if len(sel) else np.nan)
        ax[1].plot(mids, binned, "o-", color="#2563eb", ms=7, mec="white", mew=1.0,
                   label="binned mean")
    ax[1].axhline(0.57, color="#16a34a", lw=1.2, ls="--")
    ax[1].set_xlabel(r"$I-I_{\rm rheo}$ (pA)"); ax[1].set_ylabel("ISI coefficient of variation")
    ax[1].set_title("P2: regularity crossover (real)")
    ax[1].set_ylim(0, 1.2); ax[1].legend(fontsize=8.6, framealpha=0.95, edgecolor="#cbd5e1")

    # P3: near-rheobase ISI shape vs exponential vs quartic-FPT reference
    bins = np.linspace(0, 3.2, 46)
    if len(near):
        ax[2].hist(near, bins=bins, density=True, color="#0f172a", alpha=0.5,
                   label=fr"Allen near-rheobase (CV {near.std():.2f}, skew {skew(near):.2f})")
    ref = quartic_fpt_reference()
    ax[2].hist(ref, bins=bins, density=True, histtype="step", lw=2.2, color="#dc2626",
               label="white QIF = quartic FPT (prediction)")
    xx = np.linspace(1e-3, 3.2, 200)
    ax[2].plot(xx, np.exp(-xx), "--", color="#7c3aed", lw=1.8, label="exponential (Poisson null)")
    ax[2].set_xlabel(r"ISI $/\ \langle$ISI$\rangle$"); ax[2].set_ylabel("density")
    ax[2].set_title("P3: near-rheobase ISI law")
    ax[2].set_xlim(0, 3.2); ax[2].legend(fontsize=8.0, framealpha=0.95, edgecolor="#cbd5e1")

    fig.tight_layout()
    out = os.path.abspath(os.path.join(HERE, "..", "figures", "allen_qif_validation.png"))
    fig.savefig(out, dpi=140); print("saved", out)


if __name__ == "__main__":
    main()
