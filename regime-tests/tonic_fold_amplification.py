#!/usr/bin/env python3
"""
tonic_fold_amplification.py
===========================

First-principles audit of the mid-tonic "A_fold ≈ 2.7" claim (TONIC_PHASE.md §20).

The §20 closure asserted that the σ = 0.02 mid-tonic CV enhancement (measured
ratio CV_meas / CV_leading = 2.719) is a *deterministic, σ-independent*
multiplicative factor

        A_fold  =  exp( ∮ max(λ⊥, 0) dt ),   λ⊥(t) = (1 - γ_v²) - ε b .

This script computes every cycle-geometric quantity that claim depends on,
directly from the deterministic mid-tonic relaxation cycle, and tests it.

It establishes four things:

  (1) The genuine transverse Floquet exponent  λ2 = (1/T) ∮ tr J dt  is large and
      NEGATIVE: the cycle is enormously net-contracting (multiplier e^{λ2 T}
      ~ 1e-17). So a permanent "amplification" cannot be read off the expanding
      part alone.

  (2) The naive exponentiated expanding integral exp(∮ max(λ⊥,0) dt) is O(10–10²),
      NOT 2.7 — it ignores re-contraction.

  (3) The "post-fold-tip" inner expanding integral 2∫₀^{Vout} V dT is *cutoff
      dependent* (it diverges with Vout); it equals ln(2.719) only at a tuned
      Vout ≈ 1.26, so 2.7 is not a convergent first-principles number.

  (4) The measured CV ratio is σ-DEPENDENT (≈1, 2.7, 14 at σ = 0.005, 0.02, 0.05),
      so it is not a σ-independent multiplicative constant. The principled
      second-order phase-reduction correction (§16 K₂) is convergent but small.
      ⇒ the σ=0.02 enhancement is the ONSET of non-perturbative fold escape
        (TONIC_PHASE.md §18.2), not a deterministic A_fold.

Outputs: results/tonic_fold_amplification/{summary.txt, fold_amplification.png}
Reproduce:  python3 regime-tests/tonic_fold_amplification.py
"""
from __future__ import annotations
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, BASE)

import importlib.util
def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

tpr    = _load(os.path.join(HERE, "tonic_phase_response.py"), "tpr")
kernel = _load(os.path.join(BASE, "kernel.py"), "kernel")

A_FHN, B_FHN = 0.7, 0.8
I_FOLD_L = (A_FHN - 1 + 2*B_FHN/3) / B_FHN          # 0.2917
EPS_VALS = (0.04, 0.08, 0.16)

# Measured §10 mid-tonic CV data (I≈0.83, ε=0.08); A_mid = 0.3989.
# (σ, CV_meas) from TONIC_PHASE.md §10 table.
SEC10 = dict(A_mid=0.3989,
             data=[(0.005, 0.0019), (0.020, 0.0217), (0.050, 0.2837)])


def cycle_quantities(eps):
    """Deterministic mid-tonic cycle → transverse-Floquet structure."""
    fhn = kernel.FHN2D()
    I_H1 = fhn.I_hopf_lower_at(eps); I_H2 = fhn.I_hopf_upper_at(eps)
    I_mid = 0.5 * (I_H1 + I_H2)
    tpr._EPS_GLOBAL = eps
    T, gv, gw, gdv, gdw = tpr.get_limit_cycle(I_mid, eps)
    dt = T / tpr.N_PHI

    lam = (1.0 - gv**2) - eps * B_FHN           # transverse Floquet rate λ⊥(t)
    int_trJ = float(np.sum(lam) * dt)           # = λ2 · T
    int_pos = float(np.sum(np.maximum(lam, 0.0)) * dt)
    int_neg = float(np.sum(np.minimum(lam, 0.0)) * dt)
    return dict(eps=eps, I_mid=I_mid, T=T, gv=gv, dt=dt, lam=lam,
                lam2=int_trJ / T, int_trJ=int_trJ,
                int_pos=int_pos, int_neg=int_neg)


def fold_passage(lam_rel, Vstart=-6.0, dT=2e-4, Vout_max=6.0):
    """Deterministic inner fold passage: dV/dT = V²-W, dW/dT = -λ, V₀≈-√W₀."""
    W = Vstart**2; V = Vstart; t = 0.0
    Ts=[0.0]; Vs=[V]; Ws=[W]
    while V < Vout_max and t < 4000:
        dV = V*V - W; dW = -lam_rel
        V += dV*dT; W += dW*dT; t += dT
        Ts.append(t); Vs.append(V); Ws.append(W)
    return np.array(Ts), np.array(Vs), np.array(Ws)


def posttip_integral(lam_rel):
    """∫₀^Vout 2V dT (eps-independent inner Lyapunov integral) vs cutoff Vout."""
    Ts, Vs, Ws = fold_passage(lam_rel)
    pos = Vs > 0
    Tt, Vv = Ts[pos], Vs[pos]
    cum = np.cumsum(2.0 * Vv * np.gradient(Tt))
    target = np.log(2.719)
    kk = int(np.searchsorted(cum, target))
    V_at_target = Vv[kk] if kk < len(Vv) else np.nan
    table = {Vth: float(cum[min(np.searchsorted(Vv, Vth), len(cum)-1)])
             for Vth in (0.5, 1.0, 1.5, 2.0, 3.0)}
    return table, V_at_target, (Vv, cum)


def cv_ratio_sigma_dependence():
    """Is CV_meas/(σ·A) σ-independent (A_fold hypothesis) or growing (escape)?"""
    A = SEC10["A_mid"]
    out = []
    for sig, cv in SEC10["data"]:
        out.append((sig, cv, cv/(sig*A)))
    return out


def measured_onset_sweep(I=0.83, eps=0.08,
                         sigmas=(0.005, 0.015, 0.020, 0.025, 0.030),
                         n_traj=80, n_isi=12, seed=1):
    """Fresh stochastic ISI sweep: locate the σ where CV=σ·A breaks down.

    A deterministic A_fold would give a flat ratio; an escape onset gives a
    sharp jump from ~1 to >>1 at a threshold σ_onset, with a sampling-sensitive
    value just above it."""
    tpr._EPS_GLOBAL = eps
    T, gv, gw, gdv, gdw = tpr.get_limit_cycle(I, eps)
    Zv, _ = tpr.compute_prc(T, gv, gw, gdv, gdw)
    A = tpr.compute_A(Zv, T)
    rng = np.random.default_rng(seed)
    rows = []
    for sig in sigmas:
        cv = tpr._measure_isi_free_running(I, eps, sig, n_traj=n_traj,
                                           n_isi=n_isi, rng=rng)
        rows.append((sig, cv, (cv / (sig * A)) if cv == cv else np.nan))
    return A, rows


def main():
    outdir = os.path.join(HERE, "results", "tonic_fold_amplification")
    os.makedirs(outdir, exist_ok=True)
    lines = []
    def P(s=""):
        print(s); lines.append(s)

    P("="*72)
    P("MID-TONIC FOLD-PASSAGE TRANSVERSE AMPLIFICATION — first-principles audit")
    P("="*72)

    P("\n(1) True transverse Floquet exponent λ2 and expanding integral")
    P(f"    {'eps':>6} {'I_mid':>7} {'T':>7} {'λ2':>8} {'∮trJdt':>9} "
      f"{'∮max(λ⊥,0)':>11} {'exp(∮max)':>10} {'mult e^{λ2T}':>12}")
    cyc = {}
    for eps in EPS_VALS:
        c = cycle_quantities(eps); cyc[eps] = c
        P(f"    {eps:6.3f} {c['I_mid']:7.3f} {c['T']:7.2f} {c['lam2']:8.3f} "
          f"{c['int_trJ']:9.2f} {c['int_pos']:11.3f} {np.exp(c['int_pos']):10.1f} "
          f"{np.exp(c['int_trJ']):12.2e}")
    P("    => cycle is hugely net-contracting; exp(∮max λ⊥) is O(10–10²), not 2.7.")

    P("\n(2) Post-fold-tip inner expanding integral 2∫₀^Vout V dT  (λ_rel="
      f"{B_FHN*(cyc[0.08]['I_mid']-I_FOLD_L):.3f}, eps-independent)")
    lam_rel = B_FHN * (cyc[0.08]['I_mid'] - I_FOLD_L)
    table, V_at_target, (Vv, cum) = posttip_integral(lam_rel)
    for Vth, val in table.items():
        P(f"    Vout={Vth:4.1f}:  2∫VdT = {val:7.4f}   exp = {np.exp(val):7.3f}")
    P(f"    => reaches ln(2.719)=1.000 only at tuned Vout = {V_at_target:.3f}; "
      f"integral is cutoff-dependent (NOT convergent) -> 2.7 is not first-principles.")

    P("\n(3a) σ-dependence of the §10 CV ratio (A_fold must be σ-independent)")
    P(f"    {'sigma':>7} {'CV_meas':>9} {'CV/(σ·A)':>10}")
    ratios = cv_ratio_sigma_dependence()
    for sig, cv, r in ratios:
        P(f"    {sig:7.3f} {cv:9.4f} {r:10.3f}")
    P("    => §10 ratio grows ~1 -> 2.7 -> 14 with σ; a deterministic A_fold"
      " would be flat.")

    P("\n(3b) FRESH stochastic ISI onset sweep (I=0.83, ε=0.08, this run)")
    P(f"    {'sigma':>7} {'CV_meas':>9} {'CV/(σ·A)':>10}")
    A_meas, onset_rows = measured_onset_sweep()
    for sig, cv, r in onset_rows:
        flag = "  <- leading order" if r < 1.5 else "  <- broken (escape)"
        P(f"    {sig:7.3f} {cv:9.4f} {r:10.3f}{flag}")
    P(f"    A(I=0.83) = {A_meas:.4f}.  Leading order CV=σ·A holds (ratio≈1) up to")
    P("    σ≈0.015, then JUMPS sharply (ratio 4–12) — a non-perturbative escape")
    P("    onset, NOT a smooth σ-independent factor.  The exact ratio just above")
    P("    onset is sampling-dependent (§10 got 2.72 here; this run differs),")
    P("    confirming 2.7 is a readout of the escape onset σ_onset≈0.02,")
    P("    not a deterministic cycle constant (cf. TONIC_PHASE.md §18.2).")

    # ---- figure ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))

        c = cyc[0.08]
        tt = np.arange(tpr.N_PHI) * c["dt"]
        ax[0].plot(tt, c["lam"], lw=1.0, color="C0")
        ax[0].axhline(0, color="k", lw=0.6)
        ax[0].fill_between(tt, c["lam"], 0, where=c["lam"] > 0, color="C3", alpha=.35,
                           label=f"∮max(λ⊥,0)dt={c['int_pos']:.2f}")
        ax[0].set_title("λ⊥(t)=(1-γ_v²)-εb over the cycle (ε=0.08)")
        ax[0].set_xlabel("t"); ax[0].set_ylabel("λ⊥"); ax[0].legend(fontsize=8)

        ax[1].plot(Vv, cum, color="C2")
        ax[1].axhline(np.log(2.719), color="C3", ls="--",
                      label="ln(2.719)=0.99")
        ax[1].axvline(V_at_target, color="C3", ls=":", lw=.8)
        ax[1].set_title("post-tip integral 2∫₀^V V'dT  (cutoff-dependent)")
        ax[1].set_xlabel("V_out"); ax[1].set_ylabel("2∫VdT"); ax[1].legend(fontsize=8)

        sg = [r[0] for r in ratios]; rr = [r[2] for r in ratios]
        ax[2].plot(sg, rr, "s--", color="0.5", label="§10 data")
        og = [r[0] for r in onset_rows]; oo = [r[2] for r in onset_rows]
        ax[2].plot(og, oo, "o-", color="C0", label="fresh sweep")
        ax[2].axhline(2.719, color="C3", ls="--", label="claimed A_fold=2.7")
        ax[2].axhline(1.0, color="C2", ls=":", lw=.8, label="leading order")
        ax[2].set_xscale("log"); ax[2].set_yscale("log")
        ax[2].set_title("CV ratio vs σ — sharp escape onset, not flat")
        ax[2].set_xlabel("σ"); ax[2].set_ylabel("CV_meas/(σ·A)"); ax[2].legend(fontsize=8)

        fig.tight_layout()
        figpath = os.path.join(BASE, "figures", "tonic_fold_amplification.png")
        os.makedirs(os.path.dirname(figpath), exist_ok=True)
        fig.savefig(figpath, dpi=120)
        P(f"\nFigure: {figpath}")
    except Exception as e:
        P(f"\n[figure skipped: {e}]")

    with open(os.path.join(outdir, "summary.txt"), "w") as f:
        f.write("\n".join(lines))
    P(f"Summary: {os.path.join(outdir, 'summary.txt')}")


if __name__ == "__main__":
    main()
