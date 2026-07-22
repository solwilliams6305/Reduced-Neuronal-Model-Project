"""
ml_atlas_separability.py — the kill-or-confirm test for the excitability direction.

The Morris-Lecar check showed the NORMAL-FORM atlas does not represent conductance-model ISIs. The
deeper question that gates everything downstream: does the Type-I/Type-II distinction SURVIVE in a
conductance model's ISI statistics at all — i.e. if we TRAIN the atlas on Morris-Lecar itself, do
held-out ML cells of known class separate by SHAPE ALONE (scale-free quantile / W1, no rate cheat)?

  * confirm  -> the excitability-from-extracellular-spikes premise is real; rebuild the atlas
                biophysically and go to Allen.
  * kill     -> ISI shape is too weak in biophysics; the class must come from another observable
                (f-I onset / subthreshold resonance), or the direction retires.

Design (honest, leakage-free):
  - ATLAS: for each class, at several drive levels, pool many trials -> a mean-normalised quantile
    function (the scale-free shape), indexed by (class, drive) just like classify3's atlas.
  - TEST: HELD-OUT trials (different seeds), small realistic per-sweep samples; classify each by
    argmin W1 over ALL (class, drive) -> predicted class. Shape only; rate is normalised away.
  - Also report the complementary f-I ONSET feature (Hodgkin class) on the same cells.

    python3 ml_atlas_separability.py build    # integrate + cache ML atlas + held-out test (~30 s)
    python3 ml_atlas_separability.py test      # the separability verdict
    python3 ml_atlas_separability.py fig       # -> figures/ml_atlas_separability.png
"""
from __future__ import annotations
import os
import sys
import numpy as np

from bifurcation_classifier import qfun, W1, skew
from morris_lecar_check import ml_isis, PAR_I, PAR_II

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "_ml_separability_cache.npz")

DRIVES_I = [41.0, 45.0, 50.0, 57.0, 66.0]
DRIVES_II = [89.0, 94.0, 100.0, 110.0, 125.0]


def ml_per_neuron(I, par, M=8, T=2400.0, dt=0.05, sig=1.2, refr=6.0, seed=0):
    """One vectorised integration; returns PER-NEURON ISI arrays (held-out test sweeps)."""
    rng = np.random.default_rng(seed); n = int(T / dt); sdt = np.sqrt(dt)
    V1, V2, V3, V4, gCa, gK, gL, VCa, VK, VL, C, phi = par
    V = np.full(M, -60.0); w = np.full(M, 0.01); last = np.full(M, -1e9)
    spikes = [[] for _ in range(M)]; prev = V.copy()
    for k in range(n):
        prev = V; t = k * dt
        minf = 0.5 * (1 + np.tanh((V - V1) / V2)); winf = 0.5 * (1 + np.tanh((V - V3) / V4))
        tw = 1.0 / np.cosh((V - V3) / (2 * V4))
        V = V + dt * ((I - gL * (V - VL) - gCa * minf * (V - VCa) - gK * w * (V - VK)) / C) \
            + sig * sdt * rng.standard_normal(M)
        w = w + dt * (phi * (winf - w) / tw)
        up = (prev < 0.0) & (V >= 0.0) & ((t - last) > refr)
        if up.any():
            for i in np.where(up)[0]:
                spikes[i].append(t)
            last[up] = t
    return [np.diff(np.array(s)) for s in spikes if len(s) > 8]


def _pool_sweep(I, par, M, seed):
    return ml_per_neuron(I, par, M=M, seed=seed)


def build():
    store = {}
    # atlas: pool ~30 trials per (class, drive) for a clean shape
    for cls, par, drives in [("I", PAR_I, DRIVES_I), ("II", PAR_II, DRIVES_II)]:
        for di, I in enumerate(drives):
            pooled = ml_isis(I, par, M=30, T=1800.0, seed=1000 + di)
            store[f"atlas_{cls}_{di}"] = qfun(pooled)
            # held-out TEST sweeps: 8 independent small samples at the same drive (different seeds)
            tests = _pool_sweep(I, par, M=8, seed=9000 + 17 * di + (0 if cls == "I" else 99))
            cnt = np.array([len(s) for s in tests]); cat = np.concatenate(tests) if tests else np.array([])
            store[f"test_{cls}_{di}_cnt"] = cnt; store[f"test_{cls}_{di}_cat"] = cat
            print(f"  class {cls} drive {I}: atlas pooled, {len(tests)} held-out test sweeps")
    store["drives_I"] = np.array(DRIVES_I); store["drives_II"] = np.array(DRIVES_II)
    np.savez(CACHE, **store); print(f"  saved {CACHE}")


def _load():
    d = np.load(CACHE)
    atlas = {}
    for cls, drives in [("I", DRIVES_I), ("II", DRIVES_II)]:
        atlas[cls] = [d[f"atlas_{cls}_{di}"] for di in range(len(drives))]
    tests = {}
    for cls, drives in [("I", DRIVES_I), ("II", DRIVES_II)]:
        for di in range(len(drives)):
            cnt = d[f"test_{cls}_{di}_cnt"]; cat = d[f"test_{cls}_{di}_cat"]; o = 0; sw = []
            for c in cnt:
                sw.append(cat[o:o + c]); o += c
            tests[(cls, di)] = sw
    return atlas, tests


def classify_shape(isis, atlas):
    """argmin W1 over all (class, drive) -> (pred_class, w_best, w_I, w_II)."""
    qd = qfun(isis)
    wI = min(W1(qd, Q) for Q in atlas["I"]); wII = min(W1(qd, Q) for Q in atlas["II"])
    return ("I" if wI <= wII else "II"), min(wI, wII), wI, wII


def test(verbose=True):
    if not os.path.exists(CACHE):
        build()
    atlas, tests = _load()
    conf = np.zeros((2, 2), int); margins = []; rows = []
    for (cls, di), sweeps in tests.items():
        for s in sweeps:
            pred, wb, wI, wII = classify_shape(s, atlas)
            conf[0 if cls == "I" else 1, 0 if pred == "I" else 1] += 1
            # signed margin: how much closer the CORRECT class is (positive = correct separates)
            m = (wII - wI) if cls == "I" else (wI - wII)
            margins.append(m); rows.append((cls, di, pred, wI, wII, skew(s), len(s)))
    acc = np.trace(conf) / conf.sum(); margins = np.array(margins)
    if verbose:
        print(f"  ML-trained, shape-ONLY (scale-free) classification of held-out ML cells:")
        print(f"  confusion (rows=true I/II, cols=pred):\n{conf}")
        print(f"  per-sweep class accuracy = {acc*100:.0f}%  (n={conf.sum()})")
        print(f"  signed W1 margin (correct-class advantage): median={np.median(margins):+.4f}, "
              f"frac>0 = {np.mean(margins > 0)*100:.0f}%")
        # complementary f-I onset on the same held-out cells (rate, not shape)
        def onset(cls, drives):
            mins = [np.mean([1.0 / np.mean(s) for s in tests[(cls, di)] if len(s) > 0])
                    for di in range(len(drives))]
            mins = np.array(mins); return mins.min() / mins.max()
        print(f"  complementary f-I onset (rate): Type-I={onset('I', DRIVES_I):.2f}  "
              f"Type-II={onset('II', DRIVES_II):.2f}  (idealised normal-form gap was 0.47 vs 0.87)")
        verdict = "CONFIRM (shape separates in biophysics)" if acc >= 0.75 else \
                  ("PARTIAL (shape weak; lean on f-I onset)" if acc >= 0.6 else
                   "KILL (shape does NOT separate — need another observable)")
        print(f"  ===> VERDICT: {verdict}")
    return conf, acc, margins, rows


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) >= 2 else ""
    if cmd == "build":
        build()
    elif cmd == "test":
        test()
    elif cmd == "fig":
        from ml_separability_figure import make_figure
        make_figure()
    else:
        print("usage: ml_atlas_separability.py build | test | fig")
