"""
t2_3_coherence.py — T2.3: noise near the canard in a coupled excitable pair.
============================================================================

Two electrically coupled FitzHugh–Nagumo units (the kernel.py model) in the
EXCITABLE regime just below the canard/fold onset:
    v_i' = (v_i - v_i^3/3 - w_i + I + g(v_j - v_i)) + σ ξ_i ,
    w_i' = ε (v_i + a - b w_i) ,        a=0.7, b=0.8, I < I_fold (excitable).

Noise drives spikes; each large excursion passes near the fold (a canard-like
slow passage). Coherence resonance (Pikovsky–Kurths 1997): the inter-spike-interval
coefficient of variation CV(σ) has a MINIMUM at an optimal noise — firing is most
regular there. The coupling question (array-enhanced CR; Neiman et al. 1999): does
coupling lower that minimum / synchronize the firing?

We sweep σ for the uncoupled (g=0) and coupled (g>0) pair, pooling ISIs over an
ensemble of independent-noise realizations (online mean/var). Euler–Maruyama.

Output: figures/t2_3_coherence.png + printed summary.
"""
from __future__ import annotations
import os, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

A, B, EPS, I0 = 0.7, 0.8, 0.08, 0.20      # excitable (I0 < I_fold ≈ 0.29)
VTHR, VARM = 1.0, -0.5


def run_cr(sig_vals, g_vals, M=24, T=1500.0, dt=0.015, seed=0):
    ns, ng = len(sig_vals), len(g_vals)
    shape = (ns, ng, M, 2)
    sig = np.asarray(sig_vals)[:, None, None, None]
    g = np.asarray(g_vals)[None, :, None, None]
    v = np.full(shape, -1.2); w = np.full(shape, (-1.2 + A) / B)
    last_t = np.full(shape, np.nan); cnt = np.zeros(shape); ssum = np.zeros(shape); ssq = np.zeros(shape)
    armed = np.ones(shape, bool)
    rng = np.random.default_rng(seed); sdt = np.sqrt(dt)
    n = int(T / dt)
    for step in range(n):
        t = step * dt
        coup = g * (v[..., ::-1] - v)
        v = v + (v - v**3 / 3.0 - w + I0 + coup) * dt + sig * sdt * rng.standard_normal(shape)
        w = w + EPS * (v + A - B * w) * dt
        spike = armed & (v > VTHR)
        good = spike & np.isfinite(last_t)
        if good.any():
            isi = t - last_t
            cnt[good] += 1; ssum[good] += isi[good]; ssq[good] += isi[good] ** 2
        last_t[spike] = t
        armed[spike] = False
        armed[v < VARM] = True
    C = cnt.sum(axis=(2, 3)); S = ssum.sum(axis=(2, 3)); Q = ssq.sum(axis=(2, 3))
    C = np.where(C < 5, np.nan, C)
    mean = S / C; var = Q / C - mean**2
    CV = np.sqrt(np.clip(var, 0, None)) / mean
    rate = C / (M * 2 * T)
    return CV, rate                         # (ns, ng)


def trace(sigma, g, T=700.0, dt=0.015, seed=1):
    rng = np.random.default_rng(seed); sdt = np.sqrt(dt)
    v = np.array([-1.2, -1.2]); w = np.array([(-1.2 + A) / B] * 2)
    n = int(T / dt); ts = np.empty(n); vv = np.empty((n, 2))
    for k in range(n):
        coup = g * (v[::-1] - v)
        v = v + (v - v**3 / 3.0 - w + I0 + coup) * dt + sigma * sdt * rng.standard_normal(2)
        w = w + EPS * (v + A - B * w) * dt
        ts[k] = k * dt; vv[k] = v
    return ts, vv


def main():
    t0 = time.time()
    sig_vals = np.logspace(np.log10(0.03), np.log10(1.3), 11)
    g_vals = np.array([0.0, 0.15])
    CV, rate = run_cr(sig_vals, g_vals)

    print("=" * 70)
    print(f"T2.3 — coherence resonance in a coupled excitable FHN pair (I={I0})")
    print("=" * 70)
    print(f"  {'σ':>7} | {'CV g=0':>9} {'CV g=0.15':>10} | {'rate g=0':>9} {'rate g=.15':>10}")
    for i, s in enumerate(sig_vals):
        print(f"  {s:7.3f} | {CV[i,0]:9.3f} {CV[i,1]:10.3f} | {rate[i,0]:9.4f} {rate[i,1]:10.4f}")

    i0 = np.nanargmin(CV[:, 0]); i1 = np.nanargmin(CV[:, 1])
    print(f"\n  coherence resonance (CV minimum):")
    print(f"    g=0   : min CV {CV[i0,0]:.3f} at σ={sig_vals[i0]:.3f}")
    print(f"    g=0.15: min CV {CV[i1,1]:.3f} at σ={sig_vals[i1]:.3f}")
    valid = np.where(np.isfinite(CV[:, 0]))[0]
    imin = valid[np.argmin(CV[valid, 0])]
    cr = imin not in (valid[0], valid[-1])      # interior minimum among valid σ
    enh = np.nanmin(CV[:, 1]) < np.nanmin(CV[:, 0]) - 0.01
    print(f"  CV has an interior minimum (coherence resonance): {'PASS' if cr else 'CHECK'}")
    print(f"  coupling deepens it (array-enhanced CR): {'YES' if enh else 'no (≈ or weaker)'}")

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    ax[0].semilogx(sig_vals, CV[:, 0], "o-", color="#888780", label="g=0 (uncoupled)")
    ax[0].semilogx(sig_vals, CV[:, 1], "s-", color="#b3402b", label="g=0.15 (coupled)")
    ax[0].plot(sig_vals[i0], CV[i0, 0], "v", color="#444", ms=11)
    ax[0].plot(sig_vals[i1], CV[i1, 1], "v", color="#b3402b", ms=11)
    ax[0].set_xlabel("noise amplitude σ"); ax[0].set_ylabel("ISI coefficient of variation")
    ax[0].set_title("(A) coherence resonance: CV(σ) minimum")
    ax[0].legend(fontsize=9, frameon=False); ax[0].grid(alpha=0.3, which="both")

    sig_lo, sig_opt, sig_hi = 0.03, sig_vals[i0], 0.5
    for sg, lab, off, col in [(sig_lo, "low σ (irregular)", 8, "#2c7fb8"),
                              (sig_opt, f"optimal σ≈{sig_opt:.2f} (regular)", 4, "#2c7d59"),
                              (sig_hi, "high σ (irregular)", 0, "#d95f0e")]:
        ts, vv = trace(sg, 0.0)
        win = ts < 700
        ax[1].plot(ts[win], vv[win, 0] + off, lw=0.7, color=col, label=lab)
    ax[1].set_xlabel("t"); ax[1].set_yticks([])
    ax[1].set_title("(B) single-unit firing regularity vs σ (offset)")
    ax[1].legend(fontsize=8, frameon=False, loc="upper right")

    fig.suptitle("T2.3 — noise near the canard: coherence resonance in a coupled excitable pair",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fp = os.path.join(FIG, "t2_3_coherence.png")
    fig.savefig(fp, dpi=140, bbox_inches="tight")
    print(f"\n  Figure: {fp}")
    print(f"  TOTAL {time.time()-t0:.1f}s — DONE")


if __name__ == "__main__":
    main()
