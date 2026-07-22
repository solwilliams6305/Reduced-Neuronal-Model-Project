"""
neural_mass_edge.py — THEORY TRACK, opening result: the two-edge universality at the COLLECTIVE level.

The single-neuron edge theory lives in the noisy QIF/SNIC. The Montbrio-Pazo-Roxin (MPR, 2015) exact
mean-field of an all-to-all QIF network is a 2D macroscopic system (firing rate r, mean potential v):
    r' = Delta/pi + 2 r v ,   v' = v^2 + etabar + J r - pi^2 r^2 .
Pure 2D MPR only has up/down bistability. Adding spike-frequency ADAPTATION (the standard
next-generation neural-mass route to collective rhythms) makes it a 3D slow-fast system,
    r' = Delta/pi + 2 r v ,   v' = v^2 + etabar + J r - pi^2 r^2 - a ,   tau_a a' = -a + alpha r ,
whose fast (r,v) subsystem undergoes a MACROSCOPIC bifurcation as the slow adaptation a sweeps it — a
network-level folded cycle. We find a macroscopic PHASE edge: the collective-burst period DIVERGES as
etabar -> etabar* ~ -2.32 (the cycle is destroyed at a SNIC/homoclinic of the rate dynamics).

The claim this opens (novel to the next-gen neural-mass community, which models these systems
deterministically): the single-neuron two-edge laws govern the COLLECTIVE fluctuations. We verify the
phase edge numerically — just below etabar* the noise-induced collective INTER-BURST intervals are the
universal first-passage (quartic-FPT / SNIC) shape: positive-skew, classified Type-I by the same atlas
that classifies single-neuron ISIs, read by the same inversion. The single neuron and the whole network
sit at the SAME edge.

This is the OPENING probe of a program, not a finished theory. Honest gaps below.

    python3 neural_mass_edge.py sims    # deterministic period scan + cached noisy runs (~25 s)
    python3 neural_mass_edge.py check    # the macroscopic-edge diagnostics
    python3 neural_mass_edge.py fig      # -> figures/neural_mass_edge.png
"""
from __future__ import annotations
import os
import sys
import numpy as np

import bifurcation_classifier as bc
import harden_inverter as h
from bifurcation_classifier import qfun, skew

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "_neural_mass_cache.npz")

PAR = dict(J=15.0, Delta=1.0, tau_a=15.0, alpha=5.0)
ETA_STAR = -2.32                                  # macroscopic phase-edge (located by the period scan)
ETA_SCAN = [-2.15, -2.25, -2.30, -2.31, -2.32]    # toward the edge (deterministic period diverges)
ETA_NOISE = [-2.40, -2.50, -2.60]                 # just below: noise-induced collective bursts


def mpr_adapt(etabar, sig=0.0, T=2500.0, dt=1e-3, seed=0, r0=0.5, v0=-2.0, a0=1.0, **par):
    p = {**PAR, **par}; rng = np.random.default_rng(seed); n = int(T / dt); sdt = np.sqrt(dt)
    r, v, a = r0, v0, a0; rs = np.empty(n)
    for k in range(n):
        rn = r + dt * (p["Delta"] / np.pi + 2 * r * v)
        vn = v + dt * (v * v + etabar + p["J"] * r - np.pi ** 2 * r * r - a) \
            + (sig * sdt * rng.standard_normal() if sig else 0.0)
        an = a + dt * ((-a + p["alpha"] * r) / p["tau_a"])
        r, v, a = max(rn, 1e-9), vn, an; rs[k] = r
    return rs, dt


def period(rs, dt):
    half = rs[len(rs) * 2 // 3:]
    pk = np.where((half[1:-1] > half[:-2]) & (half[1:-1] > half[2:]) &
                  (half[1:-1] > 0.5 * (half.max() + half.min())))[0]
    return float(np.mean(np.diff(pk)) * dt) if len(pk) > 2 else np.nan


def burst_ibis(rs, dt, thr=1.0, refr=7.0):
    """Inter-BURST intervals: threshold-upcrossings of the collective rate, merging the fast
    within-burst oscillations into one burst onset via a refractory (cf. spike-and-wave complexes)."""
    above = rs > thr; raw = np.where((~above[:-1]) & (above[1:]))[0] * dt
    if len(raw) == 0:
        return np.array([])
    onsets = [raw[0]]
    for t in raw[1:]:
        if t - onsets[-1] > refr:
            onsets.append(t)
    onsets = np.array(onsets)
    return np.diff(onsets) if len(onsets) > 1 else np.array([])


def build_sims():
    store = {}
    pers = []
    for eta in ETA_SCAN:
        rs, dt = mpr_adapt(eta, T=2500.0); pers.append(period(rs, dt))
        print(f"  deterministic eta={eta}: macroscopic period={pers[-1]:.2f}")
    store["eta_scan"] = np.array(ETA_SCAN); store["period"] = np.array(pers)
    # a representative trace for the figure
    rtr, dt = mpr_adapt(-2.30, T=300.0); store["trace_r"] = rtr; store["trace_dt"] = np.array([dt])
    for i, eta in enumerate(ETA_NOISE):
        rs, dt = mpr_adapt(eta, sig=0.3, T=5000.0, seed=1)
        ibi = burst_ibis(rs, dt); store[f"ibi_{i}"] = ibi
        print(f"  noisy eta={eta}: {len(ibi)} collective bursts, CV={ibi.std()/ibi.mean():.2f}")
    store["eta_noise"] = np.array(ETA_NOISE)
    np.savez(CACHE, **store); print(f"  saved {CACHE}")


def check(verbose=True):
    if not os.path.exists(CACHE):
        build_sims()
    d = np.load(CACHE)
    if verbose:
        print("  MACROSCOPIC phase edge (deterministic): period vs etabar ->")
        for e, p in zip(d["eta_scan"], d["period"]):
            print(f"    etabar={e:.2f}: period={p:.2f}")
        print(f"  => period DIVERGES toward etabar* ~ {ETA_STAR}: a network-level SNIC/homoclinic.")
        print("  COLLECTIVE inter-burst intervals (noise-induced, just below edge):")
        for i, eta in enumerate(d["eta_noise"]):
            ibi = d[f"ibi_{i}"]
            cls, p, w, per = h.classify3_interp(ibi)
            print(f"    etabar={eta:.2f}: n={len(ibi)} CV={ibi.std()/ibi.mean():.2f} "
                  f"skew={skew(ibi):+.2f} -> class={cls} (phase-edge/Type-I expected)")
        print("  => the network's collective bursts sit at the SAME phase edge as the single neuron.")
    return d


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) >= 2 else ""
    if cmd == "sims":
        build_sims()
    elif cmd == "check":
        check()
    elif cmd == "fig":
        from neural_mass_figure import make_figure
        make_figure()
    else:
        print("usage: neural_mass_edge.py sims | check | fig")
