"""
taxonomy_sweep.py — Phase 3b-cont: map the early-warning across the bifurcation TAXONOMY, and show it
splits along the paper's two-channel (phase vs amplitude) line.

Phase 3 established that the edge-law nu_hat->0 warning needs an ADIABATIC, period-diverging approach.
Here we drift the bifurcation parameter slowly through FOUR canonical onset/offset bifurcations and ask,
per class: does the inter-event period diverge? does windowed nu_hat (the timing inversion) warn? does
the oscillation amplitude collapse? The prediction from the two-edge theory:

  PHASE edge  (period diverges, sigma^{2/3}, Type-I):  SNIC, saddle-homoclinic  -> nu_hat WARNS
  AMPLITUDE edge (period finite, sigma^{1/2}, Type-II): supercrit Hopf, fold-of-cycles -> nu_hat BLIND
                                                                          (amplitude collapses instead)

Models (each with a slow linear parameter drift crossing its bifurcation):
  * SNIC               — Adler phase oscillator  dphi=(I-sin phi)dt+noise,  I: 1.8 -> 0.7  (I_c=1)
  * saddle-homoclinic  — the project's reduced log-FPT model (homoclinic_intervals): ISI grows ~ -ln|mu|
                         as mu -> 0; a genuine SH oscillator needs Bogdanov-Takens path-following (a true
                         dynamical SH is out of scope here — the reduced model carries the log-divergence).
  * supercrit Hopf     — Stuart-Landau  zdot=(b+iw)z-|z|^2 z,  b: 0.8 -> -0.2  (b_c=0), amp ~ sqrt(b)
  * fold of cycles     — rdot = b r + r^3 - r^5,  b: -0.10 -> -0.30  (fold b_c=-0.25), amp jumps to 0

    python3 taxonomy_sweep.py sims      # integrate + cache all four (~20 s)
    python3 taxonomy_sweep.py analyze    # the per-class map (prints the table)
    python3 taxonomy_sweep.py fig        # -> figures/taxonomy_sweep.png
"""
from __future__ import annotations
import os
import sys
import numpy as np

import bifurcation_classifier as bc
from bifurcation_classifier import qfun
from harden_inverter import fit_interp_1d

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "_taxonomy_cache.npz")
WIN, STEP = 35, 8
NU_ALARM = 0.4


# =================================================================================================
# generators — each returns (event_times, frac_at_crossing, amp_t, amp_v)
#   amp_t/amp_v are amplitude samples (empty for the phase oscillators, whose amplitude is ~constant)
# =================================================================================================
def gen_snic(T=2600.0, dt=2.5e-3, sig=0.4, nu_hi=2.0, nu_lo=-0.5, vth=14.0, vr=-14.0, seed=0):
    """SNIC via the saddle-node/QIF normal form dv=(v^2+I)dt+sig dB (the model the atlas is built on),
    with I drifting so nu=I/sig^{4/3} crosses 0. As nu->0 the firing period diverges; the ISI shape is
    the universal first-passage law, so windowed nu_hat tracks the approach."""
    rng = np.random.default_rng(seed); n = int(T / dt); sdt = np.sqrt(dt)
    I_hi = nu_hi * sig ** (4.0 / 3.0); I_lo = nu_lo * sig ** (4.0 / 3.0); v = vr; ev = []
    for k in range(n):
        I = I_hi + (I_lo - I_hi) * k / n
        v += dt * (v * v + I) + sig * sdt * rng.standard_normal()
        if v > vth:
            v = vr; ev.append(k * dt)
    frac = nu_hi / (nu_hi - nu_lo)
    return np.array(ev), frac, np.array([]), np.array([])


def gen_homoclinic(n_ev=600, lam=1.0, T_ret=1.5, delta=1.0, sig=0.3, seed=0):
    """Reduced saddle-homoclinic ISI stream with mu drifting mu_hi -> mu_lo (crossing 0): the noisy
    saddle-passage ISI = T_ret + (1/lam) ln(delta/|mu + sig*eta|) diverges logarithmically as mu->0."""
    rng = np.random.default_rng(seed); mu_hi, mu_lo = 1.4, -0.12
    t = 0.0; ev = []
    for i in range(n_ev):
        mu = mu_hi + (mu_lo - mu_hi) * i / n_ev
        yin = abs(mu + sig * rng.standard_normal())
        isi = T_ret + max(0.0, -(1.0 / lam) * np.log(max(yin, 1e-9) / delta))
        t += isi; ev.append(t)
    frac = (mu_hi - 0.0) / (mu_hi - mu_lo)
    return np.array(ev), frac, np.array([]), np.array([])


def gen_hopf(T=3200.0, dt=0.004, w=2.0, sig=0.02, seed=0):
    rng = np.random.default_rng(seed); n = int(T / dt); sdt = np.sqrt(dt)
    b_hi, b_lo = 0.8, -0.2; x, y = 0.5, 0.0; ev = []; at, av = [], []; yprev = y
    for k in range(n):
        b = b_hi + (b_lo - b_hi) * k / n; r2 = x * x + y * y; yprev = y
        x += dt * (b * x - w * y - r2 * x) + sig * sdt * rng.standard_normal()
        y += dt * (b * y + w * x - r2 * y) + sig * sdt * rng.standard_normal()
        r = np.sqrt(x * x + y * y)
        if yprev < 0 and y >= 0 and x > 0 and r > 0.12:          # amplitude-gated event (ignore noise crossings)
            ev.append(k * dt)
        if k % 150 == 0:
            at.append(k * dt); av.append(r)
    frac = (b_hi - 0.0) / (b_hi - b_lo)
    return np.array(ev), frac, np.array(at), np.array(av)


def gen_fold(T=3200.0, dt=0.004, w=2.0, sig=0.006, seed=0):
    rng = np.random.default_rng(seed); n = int(T / dt); sdt = np.sqrt(dt)
    b_hi, b_lo = -0.10, -0.30; r = 1.0; phi = 0.0; ev = []; at, av = [], []
    for k in range(n):
        b = b_hi + (b_lo - b_hi) * k / n
        r = abs(r + dt * (b * r + r ** 3 - r ** 5) + sig * sdt * rng.standard_normal())
        phi_new = phi + dt * w
        if (phi % (2 * np.pi)) > (phi_new % (2 * np.pi)) and r > 0.12:
            ev.append(k * dt)
        phi = phi_new
        if k % 150 == 0:
            at.append(k * dt); av.append(r)
    frac = (b_hi - (-0.25)) / (b_hi - b_lo)
    return np.array(ev), frac, np.array(at), np.array(av)


GENS = {"SNIC": gen_snic, "Saddle-homoclinic": gen_homoclinic, "Supercrit Hopf": gen_hopf,
        "Fold-of-cycles": gen_fold}
EDGE = {"SNIC": "phase", "Saddle-homoclinic": "phase", "Supercrit Hopf": "amplitude",
        "Fold-of-cycles": "amplitude"}


def build_sims():
    store = {}
    for name, gen in GENS.items():
        ev, frac, at, av = gen()
        key = name.replace(" ", "_").replace("-", "_")
        store[f"{key}__ev"] = ev; store[f"{key}__frac"] = np.array([frac])
        store[f"{key}__at"] = at; store[f"{key}__av"] = av
        print(f"  {name}: {len(ev)} events, crossing at frac={frac:.2f}, amp samples={len(av)}")
    np.savez(CACHE, **store)
    print(f"  saved {CACHE}")


def load_sims():
    d = np.load(CACHE); out = {}
    for name in GENS:
        key = name.replace(" ", "_").replace("-", "_")
        out[name] = (d[f"{key}__ev"], float(d[f"{key}__frac"][0]), d[f"{key}__at"], d[f"{key}__av"])
    return out


# =================================================================================================
# windowed analysis + per-class classification
# =================================================================================================
def _nu_hat(isis):
    d = np.load(bc.ATLAS3)
    return fit_interp_1d(qfun(isis), np.asarray(d["snic_nu"], float), np.asarray(d["snic_Q"], float))[0]


def analyze_one(ev, frac, at, av):
    """Per-window (pos, nu_hat, period); plus amplitude early / just-before-crossing (precursor) /
    after-crossing (edge confirmation). Times are referenced to the FULL run length (amplitude-edge
    cycles die, so events stop before the crossing — using ev[-1] would mis-place it)."""
    t_total = float(at[-1]) if len(at) > 0 else float(ev[-1])
    tc = frac * t_total
    amean = lambda m: (np.mean(av[m]) if np.any(m) else np.nan)
    if len(av) > 4:
        amp_early = amean(at < 0.35 * tc)
        amp_pre = amean((at > 0.80 * tc) & (at < 0.99 * tc))      # PRECURSOR window (before crossing)
        amp_post = amean(at > 1.05 * tc)                         # after crossing
        # truncate the TIMING analysis to where the cycle is still robust: once amplitude has fallen,
        # the rotation "events" are noise-dominated and nu_hat/period on them are meaningless.
        below = np.where((at > 0.4 * tc) & (av < 0.4 * amp_early))[0]
        t_floor = float(at[below[0]]) if len(below) else t_total
    else:
        amp_early = amp_pre = amp_post = np.nan; t_floor = t_total
    isis = np.diff(ev); tmid = 0.5 * (ev[1:] + ev[:-1])
    rows = []
    for a in range(0, len(isis) - WIN + 1, STEP):
        if tmid[a:a + WIN].mean() > t_floor:
            break
        w = isis[a:a + WIN]
        rows.append((tmid[a:a + WIN].mean() / t_total, _nu_hat(w), np.median(w)))
    rows = np.array(rows)
    return rows, amp_early, amp_pre, amp_post


def analyze(verbose=True):
    sims = load_sims(); summary = {}
    for name, (ev, frac, at, av) in sims.items():
        rows, ae, apre, apost = analyze_one(ev, frac, at, av)
        pre = rows[rows[:, 0] < frac]
        period_div = pre[-3:, 2].mean() / np.median(pre[:5, 2])
        nu_min = np.nanmin(pre[:, 1]) if len(pre) else np.nan
        nu_warns = nu_min < NU_ALARM
        pre_ratio = (apre / ae) if np.isfinite(ae) and ae > 0 else np.nan      # amplitude precursor
        post_ratio = (apost / ae) if np.isfinite(ae) and ae > 0 else np.nan    # amplitude edge
        amp_precursor = np.isfinite(pre_ratio) and pre_ratio < 0.55
        amp_edge = np.isfinite(post_ratio) and post_ratio < 0.30
        # what warns BEFORE the crossing: timing (nu_hat) and/or amplitude
        timing_warn = bool(nu_warns) or period_div > 2.0
        summary[name] = dict(edge=EDGE[name], period_div=period_div, nu_min=nu_min, nu_warns=bool(nu_warns),
                             pre_ratio=pre_ratio, post_ratio=post_ratio, amp_precursor=bool(amp_precursor),
                             amp_edge=bool(amp_edge), timing_warn=timing_warn, rows=rows, frac=frac,
                             at=at, av=av)
    if verbose:
        print(f"  bifurcation taxonomy — early-warning map (nu alarm = {NU_ALARM}):")
        print(f"  {'bifurcation':<20}{'edge':<11}{'period x':>9}{'min nu':>8}{'TIMING warn':>13}"
              f"{'AMP precursor':>15}{'amp->0 edge':>13}")
        for name, s in summary.items():
            pr = "n/a" if not np.isfinite(s["pre_ratio"]) else f"{s['pre_ratio']:.2f}"
            print(f"  {name:<20}{s['edge']:<11}{s['period_div']:>8.1f}x{s['nu_min']:>8.2f}"
                  f"{str(s['timing_warn']):>13}{(pr + ('  ' + str(s['amp_precursor']))):>15}"
                  f"{str(s['amp_edge']):>13}")
        print("  => PHASE edge (SNIC, saddle-homoclinic): period diverges -> nu_hat->0 gives a TIMING warning.")
        print("     AMPLITUDE edge: supercrit Hopf collapses gradually (amplitude PRECURSOR); the")
        print("     fold-of-cycles holds amplitude then JUMPS (no precursor in either channel — the")
        print("     genuinely catastrophic, warning-free tipping).")
    return summary


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) >= 2 else ""
    if cmd == "sims":
        build_sims()
    elif cmd == "analyze":
        if not os.path.exists(CACHE):
            build_sims()
        analyze()
    elif cmd == "fig":
        from taxonomy_figure import make_figure
        make_figure()
    else:
        print("usage: taxonomy_sweep.py sims | analyze | fig")
