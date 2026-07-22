"""
bifurcation_classifier.py — MVP rung 2 (Phase 1): full-density, MULTI-CLASS inversion.

Upgrades rung 1 (`bifurcation_inverter.py`: CV+skew, SNIC only) to:
  * full normalized-density matching via the 1-Wasserstein distance on quantile functions
    (the distance at the best-fit parameter IS the class-check / goodness-of-fit);
  * a multi-class atlas {SNIC (Type-I), Hopf (Type-II)} indexed by each class's distance parameter;
  * classify(intervals) -> (class, param_hat, W1) = argmin over ALL (class, param).

Headline test: at MATCHED CV, SNIC and Hopf are confusable by CV alone but cleanly SEPARATED by the
full density — SNIC intervals are the skewed quartic FPT (mode away from 0, right tail); Hopf
intervals are the near-symmetric jitter of the rotation period. Density matching sees the difference.

    python3 bifurcation_classifier.py atlas      # build + cache {SNIC, Hopf} atlas (one-time)
    python3 bifurcation_classifier.py phase1     # -> figures/bifurcation_classifier.png
"""
from __future__ import annotations
import os
import sys
import numpy as np
from offcritical_phase_edge import fpt_sim

HERE = os.path.dirname(os.path.abspath(__file__))
ATLAS = os.path.join(HERE, "_classifier_atlas.npz")
NQ = 100
PGRID = np.linspace(0.5 / NQ, 1 - 0.5 / NQ, NQ)          # quantile probabilities


def skew(x):
    x = np.asarray(x, float); m = x.mean(); s = x.std()
    return float(np.mean(((x - m) / s) ** 3)) if s > 0 and len(x) > 2 else np.nan


def qfun(samples):
    """Mean-normalized quantile function (the sigma-free / scale-free shape)."""
    s = np.asarray(samples, float); s = s[s > 0]
    return np.quantile(s / s.mean(), PGRID)


# ---- interval generators -------------------------------------------------------------------------
def snic_intervals(nu, N=5000, seed=0):
    Te, m, sd, sk = fpt_sim(float(nu), N=N, seed=seed)
    return Te


def hopf_intervals(mu, sigma=0.5, omega=2.5, ntraj=40, T=240.0, dt=2.0e-3, seed=0):
    """Rotation-period intervals of a noisy Stuart-Landau oscillator: times of +x-axis up-crossings
    (y: -> +, x>0), one per cycle, with a refractory of 0.4 period."""
    rng = np.random.default_rng(seed)
    x = 0.1 * rng.standard_normal(ntraj); y = 0.1 * rng.standard_normal(ntraj)
    sdt = np.sqrt(dt); refr = 0.4 * (2 * np.pi / omega)
    last_t = np.full(ntraj, -1e9); prev_t = np.full(ntraj, np.nan); isis = []
    yprev = y.copy()
    for k in range(int(T / dt)):
        t = k * dt; r2 = x * x + y * y
        x += (mu * x - omega * y - r2 * x) * dt + sigma * sdt * rng.standard_normal(ntraj)
        y += (mu * y + omega * x - r2 * y) * dt + sigma * sdt * rng.standard_normal(ntraj)
        cross = (yprev < 0) & (y >= 0) & (x > 0) & (t - last_t > refr)
        if cross.any():
            idx = np.where(cross)[0]
            for i in idx:
                if np.isfinite(prev_t[i]):
                    isis.append(t - prev_t[i])
                prev_t[i] = t; last_t[i] = t
        yprev = y.copy()
    return np.array(isis)


# ---- atlas ---------------------------------------------------------------------------------------
def build_atlas():
    snic_nu = np.linspace(-0.6, 3.0, 14)
    hopf_mu = np.linspace(-0.6, 1.4, 12)
    Qs = np.array([qfun(snic_intervals(nu, seed=11)) for nu in snic_nu])
    print("  SNIC atlas built")
    Qh = np.array([qfun(hopf_intervals(mu, seed=11)) for mu in hopf_mu])
    print("  Hopf atlas built")
    np.savez(ATLAS, snic_nu=snic_nu, snic_Q=Qs, hopf_mu=hopf_mu, hopf_Q=Qh)
    print(f"  saved {ATLAS}")


def _atlas():
    d = np.load(ATLAS)
    return {"SNIC": (d["snic_nu"], d["snic_Q"]), "Hopf": (d["hopf_mu"], d["hopf_Q"])}


def W1(qa, qb):
    return float(np.mean(np.abs(qa - qb)))


def classify(intervals):
    """-> (best_class, best_param, best_W1, {class: (param, W1)})."""
    qd = qfun(intervals); atl = _atlas(); per = {}
    for cls, (params, Q) in atl.items():
        d = np.array([W1(qd, Q[j]) for j in range(len(params))])
        j = int(np.argmin(d)); per[cls] = (float(params[j]), float(d[j]))
    best = min(per.items(), key=lambda kv: kv[1][1])
    return best[0], best[1][0], best[1][1], per


# ---- Phase-1 validation --------------------------------------------------------------------------
def phase1():
    if not os.path.exists(ATLAS):
        build_atlas()

    # (1) find a SNIC nu and a Hopf mu with MATCHED CV (~0.25) for the confusion panel
    def cv_of(samp):
        s = samp[samp > 0]; return s.std() / s.mean()
    snic_demo = snic_intervals(1.45, seed=5); hopf_demo = hopf_intervals(0.2, seed=5)
    print(f"  matched-CV demo:  SNIC CV={cv_of(snic_demo):.3f} skew={skew(snic_demo):.2f} | "
          f"Hopf CV={cv_of(hopf_demo):.3f} skew={skew(hopf_demo):.2f}")

    # (2) classification test: many series from each true class, varied params + seeds
    tests = []   # (true_class, true_param, pred_class, w1_snic, w1_hopf, param_hat)
    for nu in [0.4, 0.9, 1.45, 2.0]:
        for sd in range(4):
            pred, ph, w, per = classify(snic_intervals(nu, N=2500, seed=100 + sd))
            tests.append(("SNIC", nu, pred, per["SNIC"][1], per["Hopf"][1], ph))
    for mu in [-0.3, 0.1, 0.5, 1.0]:
        for sd in range(4):
            pred, ph, w, per = classify(hopf_intervals(mu, ntraj=30, seed=100 + sd))
            tests.append(("Hopf", mu, pred, per["SNIC"][1], per["Hopf"][1], ph))
    acc = np.mean([t[0] == t[2] for t in tests])
    print(f"  classification accuracy = {acc*100:.0f}%  ({len(tests)} test series)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    fig, ax = plt.subplots(1, 3, figsize=(15.8, 4.8))
    CS, CH = "#dc2626", "#2563eb"

    # Panel A: matched-CV densities -- same CV, different shape
    bins = np.linspace(0, 3.0, 42)
    ax[0].hist(snic_demo / snic_demo.mean(), bins=bins, density=True, color=CS, alpha=0.5,
               label=fr"SNIC (CV {cv_of(snic_demo):.2f}, skew {skew(snic_demo):.1f})")
    ax[0].hist(hopf_demo / hopf_demo.mean(), bins=bins, density=True, histtype="step", lw=2.4,
               color=CH, label=fr"Hopf (CV {cv_of(hopf_demo):.2f}, skew {skew(hopf_demo):.1f})")
    ax[0].set_xlabel(r"interval $/\ \langle$interval$\rangle$"); ax[0].set_ylabel("density")
    ax[0].set_title("Matched CV, different density")
    ax[0].legend(fontsize=8.2, framealpha=0.95, edgecolor="#cbd5e1")

    # Panel B: classification map -- W1 to each atlas, colored by true class
    T = tests
    for cls, col, mk in [("SNIC", CS, "o"), ("Hopf", CH, "D")]:
        xs = [t[3] for t in T if t[0] == cls]; ys = [t[4] for t in T if t[0] == cls]
        ax[1].scatter(xs, ys, c=col, marker=mk, s=55, edgecolor="white", lw=0.8, label=f"true {cls}")
    lim = [0, max(max(t[3], t[4]) for t in T) * 1.1]
    ax[1].plot(lim, lim, "--", color="#94a3b8", lw=1.2, label="W1$_{SNIC}$=W1$_{Hopf}$")
    ax[1].set_xlabel(r"W1 to SNIC atlas"); ax[1].set_ylabel(r"W1 to Hopf atlas")
    ax[1].set_title(fr"Density classifier: {acc*100:.0f}% correct")
    ax[1].set_xlim(lim); ax[1].set_ylim(lim)
    ax[1].legend(fontsize=8.2, framealpha=0.95, edgecolor="#cbd5e1", loc="upper right")

    # Panel C: within-class recovery (SNIC nu_hat vs true)
    snic_T = [t for t in T if t[0] == "SNIC" and t[2] == "SNIC"]
    tn = np.array([t[1] for t in snic_T]); ph = np.array([t[5] for t in snic_T])
    ax[2].plot([0.2, 2.2], [0.2, 2.2], "--", color="#94a3b8", lw=1.2, label="ideal")
    ax[2].plot(tn + np.random.default_rng(0).normal(0, 0.01, len(tn)), ph, "o", color=CS, ms=7,
               mec="white", mew=1.0)
    ax[2].set_xlabel(r"true $\nu$ (SNIC)"); ax[2].set_ylabel(r"recovered $\hat\nu$")
    ax[2].set_title("Within-class recovery")
    ax[2].legend(fontsize=9, framealpha=0.95, edgecolor="#cbd5e1")

    fig.tight_layout()
    out = os.path.abspath(os.path.join(HERE, "..", "figures", "bifurcation_classifier.png"))
    fig.savefig(out, dpi=140); print("saved", out)


# =================================================================================================
# Phase 1+ : a third class (saddle-homoclinic) and a confound (adaptation) axis
# =================================================================================================
ATLAS3 = os.path.join(HERE, "_classifier_atlas3.npz")
ADAPT = os.path.join(HERE, "_adapt_atlas.npz")


def homoclinic_intervals(mu, sigma=0.3, lam=1.0, T_ret=1.5, delta=1.0, N=6000, seed=0):
    """Saddle-homoclinic onset (one of Izhikevich's four canonical spiking bifurcations): the ISI is
    the noisy saddle-passage time, ISI = T_ret + (1/lam) ln(delta/|mu + sigma*eta|).  Signature: a
    HARD lower edge at T_ret plus a LOGARITHMIC/exponential right tail (period diverges as mu->0)."""
    rng = np.random.default_rng(seed)
    yin = np.abs(mu + sigma * rng.standard_normal(N))
    Ts = np.clip(-(1.0 / lam) * np.log(np.clip(yin, 1e-12, None) / delta), 0.0, None)
    return T_ret + Ts


def build_atlas3():
    d = np.load(ATLAS)
    homo_mu = np.linspace(0.04, 1.3, 12)
    Qc = np.array([qfun(homoclinic_intervals(mu, seed=11)) for mu in homo_mu])
    np.savez(ATLAS3, snic_nu=d["snic_nu"], snic_Q=d["snic_Q"], hopf_mu=d["hopf_mu"],
             hopf_Q=d["hopf_Q"], homo_mu=homo_mu, homo_Q=Qc)
    print(f"  saved {ATLAS3} (SNIC + Hopf + Homoclinic)")


def _atlas3():
    d = np.load(ATLAS3)
    return {"SNIC": (d["snic_nu"], d["snic_Q"]), "Hopf": (d["hopf_mu"], d["hopf_Q"]),
            "Homoclinic": (d["homo_mu"], d["homo_Q"])}


def classify3(intervals):
    qd = qfun(intervals); per = {}
    for cls, (params, Q) in _atlas3().items():
        dd = np.array([W1(qd, Q[j]) for j in range(len(params))])
        j = int(np.argmin(dd)); per[cls] = (float(params[j]), float(dd[j]))
    best = min(per.items(), key=lambda kv: kv[1][1])
    return best[0], best[1][0], best[1][1], per


def adapt_snic_intervals(nu, b, sigma=0.4, tau_w=6.0, v_th=14.0, v_reset=-14.0,
                         N=400, dt=2.5e-3, n_isi=9, seed=0):
    """Steady-state ISIs of a SNIC neuron WITH spike-frequency adaptation (increment b, decay tau_w);
    b=0 recovers the plain SNIC.  First two ISIs per neuron (the adaptation transient) are dropped."""
    from offcritical_phase_edge import J_of_nu
    rng = np.random.default_rng(seed); I = nu * sigma ** (4.0 / 3.0)
    v = np.full(N, float(v_reset)); w = np.zeros(N); last = np.full(N, np.nan)
    cnt = np.zeros(N, int); isis = []; sdt = np.sqrt(dt)
    Tmax = (n_isi + 3) * sigma ** (-2.0 / 3.0) * max(J_of_nu(nu), 2.0)
    for k in range(int(Tmax / dt)):
        t = k * dt
        v += (v * v + I - w) * dt + sigma * sdt * rng.standard_normal(N)
        w -= (w / tau_w) * dt
        fired = v > v_th
        if fired.any():
            ok = fired & (cnt >= 2) & np.isfinite(last)
            if ok.any():
                isis.extend((t - last[ok]).tolist())
            last[fired] = t; cnt[fired] += 1; w[fired] += b; v[fired] = v_reset
    return np.array(isis)


def build_adapt_atlas():
    nu_g = np.array([0.0, 0.4, 0.8, 1.2, 1.7, 2.3]); b_g = np.array([0.0, 0.12, 0.25, 0.45])
    Q = np.empty((len(nu_g), len(b_g), NQ))
    for i, nu in enumerate(nu_g):
        for jb, b in enumerate(b_g):
            Q[i, jb] = qfun(adapt_snic_intervals(nu, b, seed=20 + i * 7 + jb))
        print(f"  adapt atlas nu={nu:+.1f} done")
    np.savez(ADAPT, nu_g=nu_g, b_g=b_g, Q=Q); print(f"  saved {ADAPT} (SNIC x adaptation)")


def fit_adapt(intervals, joint=True):
    d = np.load(ADAPT); qd = qfun(intervals); nu_g, b_g, Q = d["nu_g"], d["b_g"], d["Q"]
    if joint:
        best = (1e9, np.nan, np.nan)
        for i in range(len(nu_g)):
            for jb in range(len(b_g)):
                w = W1(qd, Q[i, jb])
                if w < best[0]:
                    best = (w, float(nu_g[i]), float(b_g[jb]))
        return best[1], best[2], best[0]
    dd = np.array([W1(qd, Q[i, 0]) for i in range(len(nu_g))])     # naive: b=0 row only
    i = int(np.argmin(dd)); return float(nu_g[i]), 0.0, float(dd[i])


def phase1plus():
    if not os.path.exists(ATLAS3):
        build_atlas3()
    if not os.path.exists(ADAPT):
        build_adapt_atlas()
    classes = ["SNIC", "Hopf", "Homoclinic"]
    gen = {"SNIC": lambda p, s: snic_intervals(p, N=2500, seed=s),
           "Hopf": lambda p, s: hopf_intervals(p, ntraj=30, seed=s),
           "Homoclinic": lambda p, s: homoclinic_intervals(p, seed=s)}
    tp = {"SNIC": [0.4, 1.0, 1.8], "Hopf": [-0.2, 0.4, 1.0], "Homoclinic": [0.1, 0.4, 0.9]}

    conf = np.zeros((3, 3), int)
    for ti, tc in enumerate(classes):
        for p in tp[tc]:
            for s in range(4):
                pred, ph, w, per = classify3(gen[tc](p, 200 + s))
                conf[ti, classes.index(pred)] += 1
    acc = np.trace(conf) / conf.sum()
    print(f"  3-class confusion (rows=true SNIC/Hopf/Homo, cols=pred):\n{conf}\n  accuracy={acc*100:.0f}%")

    true_nu = np.array([0.4, 0.8, 1.2, 1.7]); b_true = 0.25
    naive, joint, bhat = [], [], []
    for nu in true_nu:
        data = adapt_snic_intervals(nu, b_true, seed=int(50 + 100 * nu))
        nn, _, _ = fit_adapt(data, joint=False); nj, bj, _ = fit_adapt(data, joint=True)
        naive.append(nn); joint.append(nj); bhat.append(bj)
        print(f"  adapt ν={nu:+.2f} (b={b_true}): naive ν̂={nn:+.2f} | joint ν̂={nj:+.2f}, b̂={bj:.2f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    fig, ax = plt.subplots(1, 3, figsize=(15.8, 4.8))
    bins = np.linspace(0, 3.2, 44)
    demos = [("SNIC", snic_intervals(1.0, seed=5), "#dc2626"),
             ("Hopf", hopf_intervals(0.4, seed=5), "#2563eb"),
             ("Homoclinic", homoclinic_intervals(0.2, seed=5), "#16a34a")]
    for nm, dat, col in demos:
        ax[0].hist(dat / dat.mean(), bins=bins, density=True, histtype="step", lw=2.3, color=col,
                   label=fr"{nm} (skew {skew(dat):.1f})")
    ax[0].set_xlabel(r"interval $/\ \langle$interval$\rangle$"); ax[0].set_ylabel("density")
    ax[0].set_title("Three canonical onset shapes")
    ax[0].legend(fontsize=8.6, framealpha=0.95, edgecolor="#cbd5e1")

    im = ax[1].imshow(conf, cmap="Blues", vmin=0)
    ax[1].set_xticks(range(3)); ax[1].set_yticks(range(3))
    ax[1].set_xticklabels(classes, fontsize=9); ax[1].set_yticklabels(classes, fontsize=9)
    for i in range(3):
        for j in range(3):
            ax[1].text(j, i, conf[i, j], ha="center", va="center",
                       color="white" if conf[i, j] > conf.max() / 2 else "#0f172a", fontsize=12)
    ax[1].set_xlabel("predicted"); ax[1].set_ylabel("true")
    ax[1].set_title(fr"3-class confusion: {acc*100:.0f}% correct")

    ax[2].plot([0.2, 2.0], [0.2, 2.0], "--", color="#94a3b8", lw=1.2, label="ideal")
    ax[2].plot(true_nu, naive, "o", color="#dc2626", ms=9, mec="white", mew=1.0,
               label="naive (b=0 atlas) — biased")
    ax[2].plot(true_nu, joint, "s", color="#16a34a", ms=9, mec="white", mew=1.0,
               label="joint (ν,b) fit — corrected")
    ax[2].set_xlabel(r"true $\nu$ (adapting data, b=0.25)"); ax[2].set_ylabel(r"recovered $\hat\nu$")
    ax[2].set_title("Confound correction via joint fit")
    ax[2].legend(fontsize=8.4, framealpha=0.95, edgecolor="#cbd5e1")

    fig.tight_layout()
    out = os.path.abspath(os.path.join(HERE, "..", "figures", "bifurcation_classifier_phase1plus.png"))
    fig.savefig(out, dpi=140); print("saved", out)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) >= 2 else ""
    if cmd == "atlas":
        build_atlas()
    elif cmd == "atlas3":
        build_atlas3()
    elif cmd == "adaptatlas":
        build_adapt_atlas()
    elif cmd == "phase1":
        phase1()
    elif cmd == "phase1plus":
        phase1plus()
    else:
        print("usage: bifurcation_classifier.py atlas | atlas3 | adaptatlas | phase1 | phase1plus")
