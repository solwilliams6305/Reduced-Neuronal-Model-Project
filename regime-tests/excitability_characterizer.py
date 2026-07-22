"""
excitability_characterizer.py — read EXCITABILITY CLASS and PROXIMITY-TO-THRESHOLD off a cell's
extracellular spike statistics, with a goodness-of-fit that knows when it can't.

The reframe (see EXCITABILITY_CHARACTERIZER.md): instead of "early-warning of a tipping point", point
the edge-theory machinery at a question computational neuroscience actually asks and usually needs
intracellular access to answer — is this neuron Type-I (SNIC) or Type-II (Hopf) excitable, and how far
above its spiking threshold is it? — using only spike trains across a few drive levels (e.g. Allen
Long-Square sweeps).

A cell report combines THREE features (no single one suffices; the smoke test showed f-I onset alone
conflates Type-II with Poisson):
  1. ISI-SHAPE class — near-threshold-weighted vote of the off-grid 3-class classifier (Type-I = the
     skewed quartic-FPT shape, SNIC/homoclinic; Type-II = the near-symmetric Hopf-rotation jitter);
  2. f-I ONSET (Hodgkin Class 1 vs 2) — Type-I fires at arbitrarily low rate near rheobase (min/max
     sweep rate -> small); Type-II jumps to a bounded onset frequency (min/max -> O(1));
  3. GOODNESS-OF-FIT residual — the min W1 to the in-family atlas; large => OUT-OF-FAMILY (reject), the
     self-diagnosis that keeps the tool honest.

Plus the confound-corrected nu_hat across drive (proximity-to-threshold) and an adaptation flag.

    python3 excitability_characterizer.py demo       # characterize a few example cells
    python3 excitability_characterizer.py validate   # synthetic-cohort validation (prints metrics)
    python3 excitability_characterizer.py fig        # -> figures/excitability_characterizer.png
"""
from __future__ import annotations
import os
import sys
import numpy as np

import bifurcation_classifier as bc
import harden_inverter as h
from harden_inverter import classify3_interp, aggregate_class, fit_adapt_interp, TYPE_I

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "_excitability_cohort.npz")


# =================================================================================================
# the cell characterizer
# =================================================================================================
def characterize_cell(sweeps, currents=None, reject_thr=0.07, onset_thr=0.72, n_min=8):
    """sweeps = list of ISI arrays (one per drive level). Returns a cell report dict."""
    rate, res, nu, bvals, valid = [], [], [], [], []
    for s in sweeps:
        s = np.asarray(s, float); s = s[s > 0]
        if len(s) < n_min:
            continue
        cls, p, w, per = classify3_interp(s)
        res.append(min(v[1] for v in per.values()))
        nu_c, b_c, _ = fit_adapt_interp(s, joint=True)
        nu.append(nu_c); bvals.append(b_c); rate.append(1.0 / s.mean()); valid.append(s)
    if len(valid) < 2:
        return dict(klass="uncharacterizable", reason="too few usable sweeps")
    rate = np.array(rate); res = np.array(res); nu = np.array(nu)
    # cell-level ISI-shape class (near-threshold weighted)
    call, P, _ = aggregate_class(valid)
    p_typeI = sum(P.get(c, 0.0) for c in TYPE_I)
    # f-I onset feature (Hodgkin class): bounded min rate => Type-II
    fI_onset = float(rate.min() / rate.max())
    # goodness-of-fit (median residual) -> accept/reject
    resid = float(np.median(res))
    out_of_family = resid > reject_thr
    # fused class call. The f-I ONSET (Hodgkin class) is the adaptation-robust arbiter — adaptation
    # regularises the ISI shape toward Hopf (flipping the shape vote) but does NOT move the onset
    # bifurcation, so a Type-I cell keeps its low-rate Class-1 onset. We therefore require f-I onset
    # AND the shape vote to BOTH indicate Type-II; on conflict we default to Type-I (the bias that
    # neutralises the adaptation confound).
    if out_of_family:
        klass = "out-of-family"
    elif fI_onset > onset_thr and p_typeI < 0.5:
        klass = "Type-II (Hopf)"
    else:
        klass = "Type-I (SNIC/homoclinic)"
    return dict(klass=klass, p_typeI=float(p_typeI), fI_onset=fI_onset, residual=resid,
                out_of_family=bool(out_of_family), nu_hat=nu.tolist(),
                adaptation=float(np.median(bvals)), rate=rate.tolist(),
                currents=(list(currents) if currents is not None else None))


# =================================================================================================
# synthetic cohort with KNOWN ground truth (mimics the Allen validation design)
# =================================================================================================
def typeI_cell(rheo_nu=0.3, span=1.9, n=5, sigma=0.4, b=0.0, seed=0):
    from allen_phase2_inversion import sweep_isis_batch
    nus = np.linspace(rheo_nu, rheo_nu + span, n)
    return sweep_isis_batch(nus, [sigma] * n, [b] * n, n_isi=30, seed=seed), nus


def typeII_cell(mu0=-0.2, span=1.6, n=5, seed=0):
    mus = np.linspace(mu0, mu0 + span, n)
    return [bc.hopf_intervals(mu, ntraj=12, T=120.0, dt=4.0e-3, seed=seed + i)
            for i, mu in enumerate(mus)], mus


def poisson_cell(n=5, seed=0):
    rng = np.random.default_rng(seed)
    return [rng.exponential(1.0, 350) for _ in range(n)], None


def bursting_cell(n=5, seed=0):
    rng = np.random.default_rng(seed)
    return [np.concatenate([rng.exponential(0.12, 180), rng.exponential(2.2, 180)]) for _ in range(n)], None


def build_cohort(n_each=6):
    cells = []
    for k in range(n_each):
        sw, cur = typeI_cell(rheo_nu=0.2 + 0.1 * (k % 3), sigma=0.35 + 0.03 * (k % 4),
                             b=0.0 if k % 2 else 0.2, seed=100 + k)
        cells.append(("Type-I", sw, cur))
    for k in range(n_each):
        sw, cur = typeII_cell(mu0=-0.3 + 0.05 * (k % 3), seed=300 + k)
        cells.append(("Type-II", sw, cur))
    for k in range(n_each // 2):
        sw, cur = poisson_cell(seed=500 + k); cells.append(("out-of-family", sw, cur))
        sw, cur = bursting_cell(seed=700 + k); cells.append(("out-of-family", sw, cur))
    # pack (ragged) to cache
    store = {"labels": np.array([c[0] for c in cells])}
    for i, (_, sw, cur) in enumerate(cells):
        counts = np.array([len(s) for s in sw]); cat = np.concatenate([np.asarray(s, float) for s in sw])
        store[f"c{i}_counts"] = counts; store[f"c{i}_cat"] = cat
        store[f"c{i}_cur"] = np.array(cur, float) if cur is not None else np.array([])
    store["n"] = np.array([len(cells)])
    np.savez(CACHE, **store)
    print(f"  built cohort: {len(cells)} cells "
          f"({sum(c[0]=='Type-I' for c in cells)} Type-I, {sum(c[0]=='Type-II' for c in cells)} Type-II, "
          f"{sum(c[0]=='out-of-family' for c in cells)} out-of-family)")
    print(f"  saved {CACHE}")


def load_cohort():
    d = np.load(CACHE, allow_pickle=True); n = int(d["n"][0]); cells = []
    for i in range(n):
        counts = d[f"c{i}_counts"]; cat = d[f"c{i}_cat"]; sw = []; o = 0
        for c in counts:
            sw.append(cat[o:o + c]); o += c
        cur = d[f"c{i}_cur"]; cur = cur if len(cur) else None
        cells.append((str(d["labels"][i]), sw, cur))
    return cells


def validate(verbose=True):
    if not os.path.exists(CACHE):
        build_cohort()
    cells = load_cohort()
    reps = [(lab, characterize_cell(sw, cur)) for lab, sw, cur in cells]
    # map fused class to coarse label
    def coarse(k):
        if k.startswith("Type-II"):          # check Type-II first ("Type-II" also starts with "Type-I")
            return "Type-II"
        if k.startswith("Type-I"):
            return "Type-I"
        return "out-of-family"
    labels = ["Type-I", "Type-II", "out-of-family"]
    conf = np.zeros((3, 3), int)
    for lab, rep in reps:
        conf[labels.index(lab), labels.index(coarse(rep["klass"]))] += 1
    acc = np.trace(conf) / conf.sum()
    # nu recovery on Type-I cells (truth = current nu)
    nu_err = []
    for (lab, sw, cur), (_, rep) in zip(cells, reps):
        if lab == "Type-I" and cur is not None and not rep["out_of_family"]:
            nu_err.append(np.mean(np.abs(np.array(rep["nu_hat"]) - np.array(cur))))
    if verbose:
        print(f"  cohort: {len(cells)} cells; fused-class accuracy = {acc*100:.0f}%")
        print(f"  confusion (rows=true {labels}, cols=pred):\n{conf}")
        print(f"  Type-I nu recovery: mean |nu_hat - nu_true| = {np.mean(nu_err):.3f}")
        # feature separation
        fI = {l: [] for l in labels}; rr = {l: [] for l in labels}
        for lab, rep in reps:
            if "fI_onset" in rep:
                fI[lab].append(rep["fI_onset"]); rr[lab].append(rep["residual"])
        print("  feature means by true class:")
        for l in labels:
            print(f"    {l:<14} f-I onset={np.mean(fI[l]):.2f}  residual={np.mean(rr[l]):.3f}")
    return reps, conf, acc, np.mean(nu_err)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) >= 2 else ""
    if cmd == "cohort":
        build_cohort()
    elif cmd == "demo":
        for lab, gen in [("Type-I", typeI_cell), ("Type-II", typeII_cell)]:
            sw, cur = gen(seed=1); print(f"  {lab}:", characterize_cell(sw, cur)["klass"])
    elif cmd == "validate":
        validate()
    elif cmd == "fig":
        from excitability_figure import make_figure
        make_figure()
    else:
        print("usage: excitability_characterizer.py demo | validate | fig")
