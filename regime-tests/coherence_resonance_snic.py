"""
coherence_resonance_snic.py — Direction A (neuro): the noise that best regularizes a SNIC / Type-I
(QIF) neuron scales as  sigma* ~ |I_c - I|^{3/4}.

Coherence resonance = an intermediate noise level optimises spike-train regularity.  Edge prediction:
the ISI CV of the noisy saddle-node is a function of the single rescaled current
        CV(I, sigma) = F(nu) ,   nu = I / sigma^{4/3}   (I_c = 0).
So ANY feature of the regularity-vs-noise curve sits at a fixed nu, hence at
        sigma_c(I) = (|I| / |nu_c|)^{3/4}  ~  |I|^{3/4}.
This is a falsifiable 3/4 scaling law for the coherence-optimising noise of a Type-I neuron.

Test: fix several sub-threshold currents I<0, sweep sigma, measure CV; check (1) the CV(sigma)
curves COLLAPSE onto one F(nu); (2) the iso-regularity noise sigma_c(I) ~ |I|^{3/4} (log-log slope).

    python3 coherence_resonance_snic.py compute 0   # ... 4  (one current each)
    python3 coherence_resonance_snic.py plot
"""
from __future__ import annotations
import os
import sys
import numpy as np
from offcritical_phase_edge import J_of_nu

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "_cr_cache")
os.makedirs(CACHE, exist_ok=True)

IS = -np.array([0.12, 0.22, 0.40, 0.72, 1.30])     # sub-threshold currents (I_c = 0)
NSIG = 9                                            # sigma points per current
CV_ISO = 0.70                                       # iso-regularity level for sigma_c(I)


def simulate_cv(I, sigma, N=1500, dt=2.5e-3, fpt_factor=5.0, Tcap=110.0,
                v_th=40.0, v_reset=-14.0, seed=0):
    """ISI CV = CV of the first-passage time from reset, over the whole population (unbiased:
    one sample per neuron, no fast-neuron selection bias). Censored (non-firing) neurons excluded;
    point flagged unreliable if firing fraction < 0.6."""
    rng = np.random.default_rng(seed)
    v = np.full(N, float(v_reset)); tau = np.full(N, np.nan); alive = np.ones(N, bool)
    nu = I / sigma ** (4.0 / 3.0)
    Tmax = min(fpt_factor * sigma ** (-2.0 / 3.0) * J_of_nu(nu), Tcap)
    nsteps = int(Tmax / dt); sdt = np.sqrt(dt)
    for k in range(nsteps):
        t = k * dt; na = int(alive.sum())
        if na == 0:
            break
        v[alive] += (v[alive] ** 2 + I) * dt + sigma * sdt * rng.standard_normal(na)
        cross = alive & (v > v_th)
        tau[cross] = t; alive[cross] = False
    tau = tau[np.isfinite(tau)]; frac = len(tau) / N
    if len(tau) < 30 or frac < 0.6:
        return np.nan, frac, len(tau)
    return tau.std() / tau.mean(), frac, len(tau)


def compute(i):
    I = float(IS[i]); scale = abs(I) ** 0.75
    sig_grid = np.logspace(np.log10(0.80 * scale), np.log10(6.0 * scale), NSIG)
    rows = []
    for j, sg in enumerate(sig_grid):
        cv, frac, n = simulate_cv(I, sg, seed=1000 * i + j)
        nu = I / sg ** (4.0 / 3.0)
        rows.append((I, sg, nu, cv, frac, n))
        print(f"  I={I:6.3f} sigma={sg:7.4f} nu={nu:7.3f} CV={cv:6.3f} frac={frac:5.2f} n={int(n)}")
    np.savez(os.path.join(CACHE, f"cr_{i}.npz"), rows=np.array(rows, float))
    print(f"  saved cr_{i}.npz")


def iso_sigma(sig, cv, level):
    """sigma where CV crosses `level` (descending), by log-linear interpolation."""
    m = np.isfinite(cv)
    sig, cv = sig[m], cv[m]
    for k in range(len(cv) - 1):
        if (cv[k] - level) * (cv[k + 1] - level) <= 0 and cv[k] != cv[k + 1]:
            f = (level - cv[k]) / (cv[k + 1] - cv[k])
            return np.exp(np.log(sig[k]) + f * (np.log(sig[k + 1]) - np.log(sig[k])))
    return np.nan


def plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    data = []
    for i in range(len(IS)):
        f = os.path.join(CACHE, f"cr_{i}.npz")
        if os.path.exists(f):
            data.append(np.load(f)["rows"])
    cols = ["#dc2626", "#ea580c", "#16a34a", "#2563eb", "#7c3aed"]

    def clean(d):                       # unbiased branch: well-sampled (frac>=0.95) & finite CV
        m = (d[:, 4] >= 0.95) & np.isfinite(d[:, 3])
        return d[m]

    fig, ax = plt.subplots(1, 3, figsize=(15.6, 4.8))

    # Panel A: CV vs sigma (raw) -- the regularity curves, shifted by I
    for d, c in zip(data, cols):
        I = d[0, 0]; dc = clean(d)
        ax[0].plot(dc[:, 1], dc[:, 3], "o-", color=c, ms=6, mec="white", mew=0.8,
                   label=fr"$I={I:.2f}$")
    ax[0].axhline(0.57, color="#94a3b8", lw=1.0, ls="--"); ax[0].text(ax[0].get_xlim()[1]*0.6,
              0.59, "SNIC 0.57", fontsize=8.4, color="#475569")
    ax[0].set_xscale("log"); ax[0].set_xlabel(r"noise $\sigma$"); ax[0].set_ylabel("ISI CV")
    ax[0].set_title("Regularity vs noise (per current)")
    ax[0].set_ylim(0.4, 1.05); ax[0].legend(fontsize=8.2, framealpha=0.95, edgecolor="#cbd5e1")

    # Panel B: COLLAPSE -- CV vs nu = I/sigma^{4/3}
    for d, c in zip(data, cols):
        dc = clean(d)
        ax[1].plot(dc[:, 2], dc[:, 3], "o", color=c, ms=7, mec="white", mew=0.8,
                   label=fr"$I={d[0,0]:.2f}$")
    ax[1].axhline(0.57, color="#94a3b8", lw=1.0, ls="--")
    ax[1].axhline(CV_ISO, color="#0f172a", lw=1.0, ls=":")
    ax[1].text(-2.3, CV_ISO + 0.012, fr"iso-CV {CV_ISO}", fontsize=8.2, color="#0f172a")
    ax[1].set_xlabel(r"$\nu = I/\sigma^{4/3}$"); ax[1].set_ylabel("ISI CV")
    ax[1].set_title(r"Collapse $\Rightarrow$ $\sigma$ scales as $|I|^{3/4}$")
    ax[1].set_ylim(0.4, 1.05); ax[1].legend(fontsize=8.2, framealpha=0.95, edgecolor="#cbd5e1")

    # Panel C: the 3/4 law -- sigma_c(I) at iso-CV vs |I|, log-log
    Iabs, sigc = [], []
    for d in data:
        dc = clean(d)
        sc = iso_sigma(dc[:, 1], dc[:, 3], CV_ISO)
        if np.isfinite(sc):
            Iabs.append(abs(d[0, 0])); sigc.append(sc)
    Iabs = np.array(Iabs); sigc = np.array(sigc)
    slope = b = np.nan
    if len(Iabs) >= 2:
        slope, b = np.polyfit(np.log(Iabs), np.log(sigc), 1)
        xx = np.linspace(Iabs.min() * 0.8, Iabs.max() * 1.2, 50)
        ax[2].plot(xx, np.exp(b) * xx ** slope, "-", color="#0f172a", lw=2.0,
                   label=fr"fit slope ${slope:.2f}$")
        ax[2].plot(xx, sigc[0] * (xx / Iabs[0]) ** 0.75, "--", color="#16a34a", lw=1.8,
                   label=r"$|I|^{3/4}$ (predicted)")
    ax[2].plot(Iabs, sigc, "o", color="#dc2626", ms=9, mec="white", mew=1.0, zorder=5)
    ax[2].set_xscale("log"); ax[2].set_yscale("log")
    ax[2].set_xlabel(r"$|I-I_c|$"); ax[2].set_ylabel(r"coherence-optimising $\sigma_c$")
    ax[2].set_title(r"The $3/4$ scaling law")
    ax[2].legend(fontsize=9, framealpha=0.95, edgecolor="#cbd5e1")

    fig.tight_layout()
    out = os.path.abspath(os.path.join(HERE, "..", "figures", "coherence_resonance_snic.png"))
    fig.savefig(out, dpi=140)
    print(f"saved {out}\n  iso-CV={CV_ISO}: sigma_c(|I|) log-log slope = {slope:.3f}  (predict 0.75)")


if __name__ == "__main__":
    if len(sys.argv) >= 3 and sys.argv[1] == "compute":
        compute(int(sys.argv[2]))
    elif len(sys.argv) >= 2 and sys.argv[1] == "plot":
        plot()
    else:
        print("usage: coherence_resonance_snic.py compute {0..4} | plot")
