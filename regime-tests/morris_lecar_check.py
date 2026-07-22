"""
morris_lecar_check.py — the biophysical out-of-distribution test for the excitability characterizer.

The characterizer is trained on the reduced noise-driven FIRST-PASSAGE normal forms (the atlas). Real
neurons are conductance-based. This asks the honest question: does it transfer to the Morris-Lecar
conductance model — the simplest 2D biophysical neuron, with standard Type-I (SNIC) and Type-II (Hopf)
parameter sets (Rinzel-Ermentrout)?

Verdict (run `check`): NO, not as-is. Near rheobase the ML ISIs are first-passage-like (positive skew),
but ML LEAVES that regime at moderate drive and fires regularly with a hard period ceiling (NEGATIVE
skew, low CV) — outside the atlas. The cell-level goodness-of-fit residual is elevated (~0.07-0.15 vs
<=0.03 for in-family normal forms and ~0.30 for true garbage): ML sits in a MARGINAL band, and the
reject mechanism (partly) fires — it refuses the cells rather than confidently mislabelling them, but
where it does call a class it is unreliable (ML Type-II can read as SNIC). So the 100% synthetic
in-family accuracy does NOT imply conductance-model accuracy; external validity is the open problem,
and the fix is a biophysically-BUILT atlas (or direct Allen validation). The self-diagnosing residual
is what makes this failure visible instead of silent.

    python3 morris_lecar_check.py sims    # integrate + cache ML Type-I/II cells (~30 s)
    python3 morris_lecar_check.py check    # the transfer diagnostics
    python3 morris_lecar_check.py fig      # -> figures/morris_lecar_check.png
"""
from __future__ import annotations
import os
import sys
import numpy as np

import bifurcation_classifier as bc
import harden_inverter as h
from excitability_characterizer import characterize_cell

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "_morris_lecar_cache.npz")

# (V1,V2,V3,V4,gCa,gK,gL,VCa,VK,VL,C,phi)
PAR_I = (-1.2, 18, 12, 17.4, 4.0, 8, 2, 120, -84, -60, 20, 1.0 / 15.0)   # Type-I, SNIC onset
PAR_II = (-1.2, 18, 2, 30, 4.4, 8, 2, 120, -84, -60, 20, 0.04)          # Type-II, Hopf onset


def ml_isis(I, par, M=18, T=2200.0, dt=0.05, sig=1.2, refr=6.0, seed=0):
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
    return np.concatenate([np.diff(np.array(s)) for s in spikes if len(s) > 1])


# Type-I cells span just above SNIC rheobase (~40); Type-II just above Hopf rheobase (~88)
CELLS = [("Type-I", PAR_I, [40.5, 43, 46, 50, 56], 10),
         ("Type-I", PAR_I, [41, 44, 48, 53, 60], 60),
         ("Type-II", PAR_II, [88.5, 92, 97, 104, 116], 30),
         ("Type-II", PAR_II, [89.5, 93, 99, 108, 122], 80)]


def build_sims():
    store = {"labels": np.array([c[0] for c in CELLS])}
    for i, (lab, par, cur, sd) in enumerate(CELLS):
        sweeps = [ml_isis(I, par, seed=sd + j) for j, I in enumerate(cur)]
        counts = np.array([len(s) for s in sweeps]); cat = np.concatenate(sweeps)
        store[f"c{i}_counts"] = counts; store[f"c{i}_cat"] = cat; store[f"c{i}_cur"] = np.array(cur, float)
        print(f"  cell {i} ({lab}): {[len(s) for s in sweeps]} ISIs/sweep")
    store["n"] = np.array([len(CELLS)])
    np.savez(CACHE, **store); print(f"  saved {CACHE}")


def load_sims():
    d = np.load(CACHE, allow_pickle=True); n = int(d["n"][0]); out = []
    for i in range(n):
        counts = d[f"c{i}_counts"]; cat = d[f"c{i}_cat"]; sw = []; o = 0
        for c in counts:
            sw.append(cat[o:o + c]); o += c
        out.append((str(d["labels"][i]), sw, d[f"c{i}_cur"]))
    return out


def check(verbose=True):
    if not os.path.exists(CACHE):
        build_sims()
    cells = load_sims()
    reps = []
    for lab, sw, cur in cells:
        r = characterize_cell(sw, cur); reps.append((lab, r))
    # in-family and garbage residual references
    infam = [min(v[1] for v in h.classify3_interp(bc.snic_intervals(nu, N=1500, seed=s))[3].values())
             for nu in (0.4, 1.0, 1.6) for s in (1, 2)]
    rng = np.random.default_rng(0)
    garbage = [min(v[1] for v in h.classify3_interp(rng.exponential(1.0, 600))[3].values()) for _ in range(4)]
    if verbose:
        print("  Morris-Lecar transfer test:")
        for lab, r in reps:
            print(f"    true {lab:<8} -> call={r['klass']:<26} resid={r['residual']:.3f} "
                  f"fI_onset={r['fI_onset']:.2f}")
        print(f"  residual reference: in-family normal forms = {np.mean(infam):.3f} "
              f"(reject gate 0.07); ML = {np.mean([r['residual'] for _, r in reps]):.3f}; "
              f"garbage(Poisson) = {np.mean(garbage):.3f}")
        rejected = np.mean([r['out_of_family'] for _, r in reps]) * 100
        print(f"  => ML cells: {rejected:.0f}% auto-rejected (out-of-family); the normal-form atlas does "
              f"NOT faithfully represent conductance-model ISIs. Self-diagnosis makes the gap visible.")
    return reps, float(np.mean(infam)), float(np.mean(garbage))


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) >= 2 else ""
    if cmd == "sims":
        build_sims()
    elif cmd == "check":
        check()
    elif cmd == "fig":
        from morris_lecar_figure import make_figure
        make_figure()
    else:
        print("usage: morris_lecar_check.py sims | check | fig")
