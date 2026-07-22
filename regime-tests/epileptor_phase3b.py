"""
epileptor_phase3b.py — Phase 3b: does the edge-law early-warning transfer to a PUBLISHED seizure
model? We port the Phase-3 pipeline to the **Epileptor** (Jirsa, Stacey, Quilichini, Ivanov, Bernard,
*Brain* 2014) — the canonical 5(+1)-D neural-mass model of seizure dynamics — and ask, honestly,
whether streaming the ictal discharge ISIs and inverting nu_hat in sliding windows gives an
early-warning of seizure offset comparable to Phase 3a's ramped SNIC.

Model (standard parameters, x0=-1.6):
  x1' = y1 - f1(x1,x2) - z + I1 ;  y1' = c - d x1^2 - y1
  z'  = (1/tau0)(4(x1 - x0) - z)                                   (slow permittivity -> the drift)
  x2' = -y2 + x2 - x2^3 + I2 + 2 g - 0.3(z-3.5) ;  y2' = (1/tau2)(-y2 + f2(x2))
  g'  = -gamma(g - 0.1 x1)
  f1 = x1^3-3x1^2 (x1<0) | (x2-0.6(z-4)^2)x1 (x1>=0) ;  f2 = 0 (x2<-0.25) | 6(x2+0.25) (else)

The slow z autonomously sweeps the fast subsystem through seizure onset (SNIC) and offset, so the
"drifting parameter" is INTERNAL — a stronger test than 3a's externally-imposed ramp.

What we find (run `analyze`): the discharge train is spike-and-wave (multimodal at the spike level);
at the dominant-rhythm (complex) level the ISIs are a clean Type-I shape with a LOW classifier
residual, so the inversion does not produce garbage. BUT nu_hat stays HIGH and ~FLAT across the
seizure and never approaches 0 at offset — because the Epileptor's termination is NOT a
period-diverging SNIC/homoclinic of the discharge cycle (the inter-complex period is ~constant and
ACCELERATES into an abrupt, multi-variable offset). So the nu_hat->0 warning of Phase 3a does not
transfer. This delineates the method's domain of validity: it needs an ADIABATIC approach — dwelling
near threshold long enough for the universal first-passage shape (CV->0.57, diverging period) to
develop. (Phase 3a's noise-driven ramped SNIC satisfies this; a fast relaxation-oscillator burster
does not.)

    python3 epileptor_phase3b.py sims      # integrate + cache Epileptor run (~10 s)
    python3 epileptor_phase3b.py analyze    # the transfer diagnostics (prints the table)
    python3 epileptor_phase3b.py fig        # -> figures/epileptor_phase3b.png
"""
from __future__ import annotations
import os
import sys
import numpy as np

import bifurcation_classifier as bc
from bifurcation_classifier import qfun
from harden_inverter import fit_interp_1d, classify3_interp

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "_epileptor_cache.npz")

# standard Epileptor parameters (Jirsa et al. 2014)
PAR = dict(x0=-1.6, tau0=2857.0, tau2=10.0, I1=3.1, I2=0.45, c=1.0, d=5.0, gam=0.01, sig=0.0025)
PK_THR = 0.3        # x1 peak threshold for a discharge
MERGE = 2.0         # refractory: merge sub-discharges into one spike-and-wave complex
GAP = 8.0           # inter-complex gap that separates seizures
MIN_DISCHARGES = 25
WIN, STEP = 40, 10  # sliding windows over complex ISIs within a seizure


def integrate_with_spikes(T=16000.0, dt=0.05, rec_stride=10, seed=0, **par):
    """Single run; returns (rec_t, rec_x1, rec_z) downsampled trace + full-res discharge peak times."""
    p = {**PAR, **par}
    rng = np.random.default_rng(seed); n = int(T / dt); sdt = np.sqrt(dt)
    x1 = 0.0; y1 = -5.0; z = 3.0; x2 = 0.0; y2 = 0.0; g = 0.0
    sig = p["sig"]; xm1 = x1; xm2 = x1
    peaks = []; rt, rx, rz = [], [], []
    for k in range(n):
        xm2 = xm1; xm1 = x1
        x1 = x1 + dt * (y1 - (x1 ** 3 - 3 * x1 ** 2 if x1 < 0 else (x2 - 0.6 * (z - 4) ** 2) * x1) - z + p["I1"]) \
            + sig * sdt * rng.standard_normal()
        y1 = y1 + dt * (p["c"] - p["d"] * x1 ** 2 - y1)
        z = z + dt * ((1.0 / p["tau0"]) * (4 * (x1 - p["x0"]) - z))
        x2 = x2 + dt * (-y2 + x2 - x2 ** 3 + p["I2"] + 2 * g - 0.3 * (z - 3.5)) \
            + sig * sdt * rng.standard_normal()
        y2 = y2 + dt * ((1.0 / p["tau2"]) * (-y2 + (0.0 if x2 < -0.25 else 6 * (x2 + 0.25))))
        g = g + dt * (-p["gam"] * (g - 0.1 * x1))
        if xm1 > PK_THR and xm1 > xm2 and xm1 >= x1:           # local max above threshold
            peaks.append((k - 1) * dt)
        if k % rec_stride == 0:
            rt.append(k * dt); rx.append(x1); rz.append(z)
    return np.array(rt), np.array(rx), np.array(rz), np.array(peaks)


def complexes(peaks, merge=MERGE):
    if len(peaks) == 0:
        return peaks
    out = [peaks[0]]
    for t in peaks[1:]:
        if t - out[-1] > merge:
            out.append(t)
    return np.array(out)


def seizures(ct, gap=GAP, min_n=MIN_DISCHARGES):
    if len(ct) < 2:
        return []
    isi = np.diff(ct); brk = np.where(isi > gap)[0]
    segs = np.split(np.arange(len(ct)), brk + 1)
    return [ct[s] for s in segs if len(s) >= min_n]


def build_sims(seed=0):
    rt, rx, rz, pk = integrate_with_spikes(seed=seed)
    np.savez(CACHE, rt=rt, rx=rx, rz=rz, peaks=pk)
    ct = complexes(pk); sz = seizures(ct)
    print(f"  {len(pk)} discharges, {len(ct)} complexes, {len(sz)} seizures (>= {MIN_DISCHARGES})")
    print(f"  saved {CACHE}")


def load_sims():
    d = np.load(CACHE)
    return d["rt"], d["rx"], d["rz"], d["peaks"]


# =================================================================================================
# transfer diagnostics
# =================================================================================================
def _nu_hat(isis):
    d = np.load(bc.ATLAS3)
    return fit_interp_1d(qfun(isis), np.asarray(d["snic_nu"], float), np.asarray(d["snic_Q"], float))[0]


def analyze(verbose=True):
    rt, rx, rz, pk = load_sims()
    ct = complexes(pk); sz = seizures(ct)
    spike_isi = np.diff(pk); spike_isi = spike_isi[spike_isi > 0]
    frac_short = np.mean(spike_isi < 1.0)         # spike-and-wave fast intervals
    # per-seizure: complex-ISI shape, residual, nu_hat trend, offset period ratio
    nu_pos = []         # (norm_position, nu_hat) across all windows of all seizures
    resids, calls, ratios = [], [], []
    for st in sz:
        isis = np.diff(st)
        cls, p, w, per = classify3_interp(isis)
        resids.append(min(v[1] for v in per.values())); calls.append(cls)
        med = np.median(isis); ratios.append(np.mean(isis[-5:]) / med)      # >1 diverge, <1 accelerate
        for a in range(0, len(isis) - WIN + 1, STEP):
            wv = isis[a:a + WIN]
            pos = (a + WIN / 2) / len(isis)
            nu_pos.append((pos, _nu_hat(wv)))
    nu_pos = np.array(nu_pos)
    # nu_hat binned by within-seizure position (does it fall toward offset?)
    edges = np.linspace(0, 1, 5); ctrs = 0.5 * (edges[:-1] + edges[1:])
    nu_bin = [np.mean(nu_pos[(nu_pos[:, 0] >= edges[i]) & (nu_pos[:, 0] < edges[i + 1]), 1])
              for i in range(len(ctrs))]

    if verbose:
        print(f"  {len(sz)} seizures; discharge train:")
        print(f"    spike-level multimodality: {frac_short*100:.0f}% of raw ISIs < 1.0 (spike-and-wave)")
        print(f"    complex-ISI classifier: calls={[c for c in calls[:6]]}")
        print(f"    complex-ISI residual minW1={np.mean(resids):.3f} (genuine SNIC/homoclinic ~0.004-0.011"
              f" -> IN-distribution, not garbage)")
        print(f"    offset period ratio mean(last5)/median = {np.mean(ratios):.2f} "
              f"(>1 = homoclinic slowing; <1 = ACCELERATION into offset)")
        print(f"  nu_hat vs within-seizure position (0=onset, 1=offset):")
        for c, nb in zip(ctrs, nu_bin):
            print(f"    pos {c:.2f}: nu_hat = {nb:.2f}")
        print(f"  => nu_hat stays high (~{np.nanmean(nu_bin):.1f}, far from 0) and does NOT fall toward "
              f"offset: no nu->0 early-warning. The Epileptor offset is non-adiabatic for the discharge "
              f"cycle, outside the edge-law's domain of validity (cf. Phase 3a, where it works).")
    return dict(frac_short=frac_short, resid=float(np.mean(resids)), ratio=float(np.mean(ratios)),
                ctrs=ctrs, nu_bin=nu_bin, nseiz=len(sz), nu_pos=nu_pos)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) >= 2 else ""
    if cmd == "sims":
        build_sims()
    elif cmd == "analyze":
        if not os.path.exists(CACHE):
            build_sims()
        analyze()
    elif cmd == "fig":
        from epileptor_figure import make_figure
        make_figure()
    else:
        print("usage: epileptor_phase3b.py sims | analyze | fig")
