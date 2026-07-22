"""
phase3_tipping.py — MVP Phase 3: the early-warning demo.

A SNIC neuron (QIF, dv=(v^2+I)dt+sigma dB, reset) whose input current I(t) drifts slowly DOWN toward
the SNIC at I_c=0: the neuron fires regularly, the firing slows as nu=I/sigma^{4/3} -> 0, then ceases
(the limit cycle is destroyed at the saddle-node-on-invariant-circle). The spike train is the only
observable. We stream inter-spike intervals, and in sliding windows compute three early-warning
signals FROM THE SAME ISI DATA:

  * nu_hat   — the edge-theory inversion (off-grid SNIC shape fit, harden_inverter.fit_interp_1d):
               a CALIBRATED, sigma-free distance-to-bifurcation with a meaningful zero;
  * Var      — ISI variance (the classic critical-slowing-down indicator);
  * AC1      — lag-1 autocorrelation of the ISI sequence (the textbook CSD indicator).

Benchmark (the value proposition): at a MATCHED false-alarm rate — d_crit for each method calibrated
on a STATIONARY-NULL ensemble (I fixed, no approach) — which method detects the approach earliest?
Each method alarms on a trailing-window Kendall-tau trend in its own signal (nu_hat decreasing; Var,
AC1 increasing). nu_hat additionally admits an ABSOLUTE-threshold alarm (nu_hat < nu_abs) that the
scale-free Var/AC1 cannot.

    python3 phase3_tipping.py sims     # integrate + cache approach & null spike ensembles (~15 s)
    python3 phase3_tipping.py bench    # the matched-FPR benchmark (prints the table)
    python3 phase3_tipping.py fig      # -> figures/phase3_tipping.png
"""
from __future__ import annotations
import os
import sys
import numpy as np

import bifurcation_classifier as bc
from bifurcation_classifier import qfun
from harden_inverter import fit_interp_1d

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "_phase3_cache.npz")

# ---- experiment constants ------------------------------------------------------------------------
SIGMA = 0.4
NU_HI, NU_LO = 2.0, -0.5            # approach ramp: nu from 2.0 down through 0 to -0.5
NU_NULL = 2.0                       # stationary null sits at the approach's safe starting point
T_TOTAL = 2500.0
DT = 2.5e-3
M_TRIALS = 120
WIN, STEP = 50, 12                  # fixed-COUNT ISI windows (stable estimator as rate collapses)
KTRAIL = 8                          # trailing windows for the Kendall-tau trend detector
NU_ABS = 0.4                        # absolute nu_hat alarm threshold
T_CROSS = T_TOTAL * NU_HI / (NU_HI - NU_LO)     # time when I(t)=0 (the true SNIC crossing)


# =================================================================================================
# nonstationary vectorised integrator
# =================================================================================================
def integrate(mode, M=M_TRIALS, T=T_TOTAL, dt=DT, sigma=SIGMA, v_th=14.0, v_reset=-14.0, seed=0):
    """M parallel SNIC-QIF trials under a shared schedule. mode in {'approach','null'}.
    Returns list of M sorted spike-time arrays."""
    rng = np.random.default_rng(seed); sdt = np.sqrt(dt); nsteps = int(T / dt)
    I_hi = NU_HI * sigma ** (4.0 / 3.0); I_lo = NU_LO * sigma ** (4.0 / 3.0)
    I_null = NU_NULL * sigma ** (4.0 / 3.0)
    v = np.full(M, float(v_reset)); spikes = [[] for _ in range(M)]
    for k in range(nsteps):
        t = k * dt
        I = (I_hi + (I_lo - I_hi) * (t / T)) if mode == "approach" else I_null
        v += (v * v + I) * dt + sigma * sdt * rng.standard_normal(M)
        fired = v > v_th
        if fired.any():
            for i in np.where(fired)[0]:
                spikes[i].append(t)
            v[fired] = v_reset
    return [np.asarray(s, float) for s in spikes]


def build_sims():
    ap = integrate("approach", seed=1); nu = integrate("null", seed=2)
    # cache ragged spike trains as concatenation + per-trial counts, for both ensembles
    def pack(trains):
        counts = np.array([len(s) for s in trains]); cat = np.concatenate(trains) if trains else np.array([])
        return counts, cat
    ac, aa = pack(ap); nc, na = pack(nu)
    np.savez(CACHE, ap_counts=ac, ap_cat=aa, nl_counts=nc, nl_cat=na)
    print(f"  approach: median {int(np.median(ac))} spikes/trial; null: median {int(np.median(nc))}")
    print(f"  saved {CACHE}")


def _unpack(counts, cat):
    out = []; o = 0
    for c in counts:
        out.append(cat[o:o + c]); o += c
    return out


def load_sims():
    d = np.load(CACHE)
    return _unpack(d["ap_counts"], d["ap_cat"]), _unpack(d["nl_counts"], d["nl_cat"])


# =================================================================================================
# windowed early-warning signals
# =================================================================================================
_SNIC = None


def _snic_atlas():
    global _SNIC
    if _SNIC is None:
        d = np.load(bc.ATLAS3); _SNIC = (np.asarray(d["snic_nu"], float), np.asarray(d["snic_Q"], float))
    return _SNIC


def nu_hat_snic(isis):
    params, Q = _snic_atlas()
    return fit_interp_1d(qfun(isis), params, Q)[0]


def ac1(x):
    x = np.asarray(x, float)
    if len(x) < 4 or x.std() == 0:
        return np.nan
    return float(np.corrcoef(x[:-1], x[1:])[0, 1])


def windowed_signals(spikes, win=WIN, step=STEP):
    """Slide fixed-count ISI windows -> per-window (center_time, nu_hat, var, ac1)."""
    spikes = np.sort(spikes)
    if len(spikes) < win + 2:
        return None
    isis = np.diff(spikes); tmid = 0.5 * (spikes[1:] + spikes[:-1])     # time of each ISI
    rows = []
    for a in range(0, len(isis) - win + 1, step):
        w = isis[a:a + win]
        rows.append((tmid[a:a + win].mean(), nu_hat_snic(w), np.var(w), ac1(w)))
    return np.array(rows)                                              # (nwin, 4)


# =================================================================================================
# trend detector + matched-FPR calibration
# =================================================================================================
def kendall_tau(y):
    """Kendall's tau of y vs its (monotone) index in [-1, 1] (the EWS-standard trend statistic)."""
    y = np.asarray(y, float); n = len(y)
    if n < 3:
        return 0.0
    s = 0
    for i in range(n):
        s += np.sign(y[i + 1:] - y[i]).sum()
    return float(s / (n * (n - 1) / 2))


def _pearson_trend(seg):
    """Continuous trend = Pearson correlation of seg vs index. (Kendall tau gives the same ranking but
    is discrete over k=8 windows, which prevents exact false-alarm-rate matching across methods.)"""
    seg = np.asarray(seg, float); x = np.arange(len(seg))
    if seg.std() == 0:
        return 0.0
    return float(np.corrcoef(x, seg)[0, 1])


def directed_trend(series, k, direction):
    """Trailing continuous trend over the last k windows, signed so the EXPECTED approach direction is
    POSITIVE (direction=-1 for nu_hat which falls, +1 for Var/AC1 which rise)."""
    s = np.asarray(series, float); out = np.full(len(s), -np.inf)
    for m in range(k, len(s) + 1):
        seg = s[m - k:m]
        if np.all(np.isfinite(seg)):
            out[m - 1] = direction * _pearson_trend(seg)
    return out


def first_alarm_time(times, dseries, d_crit):
    idx = np.where(dseries >= d_crit)[0]
    return float(times[idx[0]]) if len(idx) else np.inf


SIGNAL_COL = {"nu": 1, "var": 2, "ac1": 3}
SIGNAL_DIR = {"nu": -1, "var": +1, "ac1": +1}


def _trend_series(rows, name):
    return rows[:, 0], directed_trend(rows[:, SIGNAL_COL[name]], KTRAIL, SIGNAL_DIR[name])


def calibrate_dcrit(null_rows, name, target_fpr=0.15):
    """d_crit = the (1-target) quantile of the per-trial MAX trailing-trend on the null ensemble,
    so that exactly ~target_fpr of null trials raise an alarm."""
    peaks = []
    for rows in null_rows:
        _, d = _trend_series(rows, name)
        peaks.append(np.max(d[np.isfinite(d)]) if np.any(np.isfinite(d)) else -np.inf)
    return float(np.quantile(peaks, 1.0 - target_fpr))


# =================================================================================================
# benchmark
# =================================================================================================
def _all_rows(trains):
    out = []
    for s in trains:
        r = windowed_signals(s)
        if r is not None and len(r) > KTRAIL:
            out.append(r)
    return out


def confound_trains(null_trains, gmax=2.4):
    """Apply a slow MULTIPLICATIVE time-rescaling (firing rate drifts ~1 -> 1/gmax over the record):
    an electrode-drift / slow-gain confound. It inflates ISI variance (a spurious 'critical slowing'
    trend) while leaving the sigma-free ISI SHAPE — hence nu_hat — invariant."""
    out = []
    for s in null_trains:
        if len(s) < 2:
            out.append(s); continue
        isis = np.diff(s); n = len(isis)
        g = 1.0 + (gmax - 1.0) * np.arange(n) / n
        out.append(np.concatenate([[s[0]], s[0] + np.cumsum(isis * g)]))
    return out


def benchmark(target_fpr=0.15, verbose=True):
    ap_tr, nl_tr = load_sims()
    ap = _all_rows(ap_tr); nl = _all_rows(nl_tr)
    cf = _all_rows(confound_trains(nl_tr))                 # rate-drift confound null (no approach)
    res = {}
    for name in ("nu", "var", "ac1"):
        dcrit = calibrate_dcrit(nl, name, target_fpr)
        fp = np.mean([first_alarm_time(*_trend_series(r, name), dcrit) < np.inf for r in nl])
        fp_cf = np.mean([first_alarm_time(*_trend_series(r, name), dcrit) < np.inf for r in cf])
        leads, det = [], 0
        for r in ap:
            ta = first_alarm_time(*_trend_series(r, name), dcrit)
            if ta < T_CROSS:
                det += 1; leads.append(T_CROSS - ta)
        res[name] = dict(dcrit=dcrit, fpr=fp, fpr_cf=fp_cf, detect=det / len(ap),
                         lead=float(np.median(leads)) if leads else np.nan)
    # bonus: nu_hat absolute-threshold rule (unavailable to scale-free Var/AC1)
    leads_abs, det_abs = [], 0
    for r in ap:
        idx = np.where(r[:, 1] < NU_ABS)[0]
        ta = r[idx[0], 0] if len(idx) else np.inf
        if ta < T_CROSS:
            det_abs += 1; leads_abs.append(T_CROSS - ta)
    fp_abs = np.mean([np.any(r[:, 1] < NU_ABS) for r in nl])
    fp_abs_cf = np.mean([np.any(r[:, 1] < NU_ABS) for r in cf])
    res["nu_abs"] = dict(dcrit=NU_ABS, fpr=fp_abs, fpr_cf=fp_abs_cf, detect=det_abs / len(ap),
                         lead=float(np.median(leads_abs)) if leads_abs else np.nan)

    if verbose:
        print(f"  matched-FPR benchmark (target FPR={target_fpr}, true crossing t={T_CROSS:.0f}):")
        print(f"  {'method':<11}{'null FPR':>9}{'detect':>9}{'med lead':>10}{'  | confound FPR':>17}")
        names = {"nu": "nu_hat", "var": "ISI Var", "ac1": "lag-1 AC", "nu_abs": "nu_hat<thr"}
        for k in ("nu", "var", "ac1", "nu_abs"):
            v = res[k]
            print(f"  {names[k]:<11}{v['fpr']*100:>7.0f}% {v['detect']*100:>7.0f}% {v['lead']:>9.0f}"
                  f"   {v['fpr_cf']*100:>13.0f}%")
    return res, ap, nl, cf


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) >= 2 else ""
    if cmd == "sims":
        build_sims()
    elif cmd == "bench":
        if not os.path.exists(CACHE):
            build_sims()
        benchmark()
    elif cmd == "fig":
        from phase3_figure import make_figure
        make_figure()
    else:
        print("usage: phase3_tipping.py sims | bench | fig")
