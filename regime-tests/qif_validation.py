"""
qif_validation.py — synthetic GROUND-TRUTH test of the QIF phase-edge predictions, with the two
confounds that will bite in real data: COLORED (OU) noise and SPIKE-FREQUENCY ADAPTATION (sAHP).

Purpose: before touching real recordings, show (a) the estimators recover the universal sigma^{2/3}
collapse and the CV crossover for an IDEAL white-noise QIF, and (b) exactly how colored noise and
adaptation distort each signature — so a real f-I / ISI dataset can be read correctly.

Three predictions under test (from QIF_PHASE_EDGE.md):
  P1  f-I collapse:   r * sigma^{-2/3}  vs  nu = I/sigma^{4/3}   ->  Phi(nu)=1/J(nu)   (universal)
  P2  CV crossover:   ISI CV ~1 (sub) -> ~0.57 (SNIC) -> ->0 (supra)
  P3  ISI shape:      at threshold the ISI law is the QUARTIC FPT law (skewed; NOT exponential,
                      NOT Gaussian) — the white-QIF threshold ISIs ARE that law by construction.

Model (adaptive QIF with colored input noise):
    dv = (v^2 + I - w) dt + xi_noise ,   reset v->v_reset, w->w+b  when v>v_th
    dw = -(w/tau_w) dt                                            (sAHP adaptation)
  white :  xi_noise = sigma dW
  ou    :  dxi = -(xi/tau_s) dt + s sqrt(2/tau_s) dW ,  s = sigma/sqrt(2 tau_s)   (matches white as tau_s->0)

Run:
    python3 qif_validation.py compute white
    python3 qif_validation.py compute ou
    python3 qif_validation.py compute adapt
    python3 qif_validation.py plot
"""
from __future__ import annotations
import os
import sys
import numpy as np
from offcritical_phase_edge import J_of_nu

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "_qif_validation_cache")
os.makedirs(CACHE, exist_ok=True)

SIGMAS = [0.5, 0.3]
NUS = np.array([-0.4, 0.0, 0.5, 1.0, 1.6, 2.4])
# confound parameters (physical units; ISIs are O(sigma^{-2/3} J) ~ O(10))
TAU_S = 0.6          # OU correlation time (colored noise)
TAU_W = 6.0          # sAHP decay time (slow)
B_ADAPT = 0.18       # spike-triggered adaptation increment


def Phi(nu):
    nu = np.asarray(nu, float)
    if nu.ndim == 0:
        return 1.0 / J_of_nu(float(nu))
    return np.array([1.0 / J_of_nu(float(x)) for x in nu.ravel()]).reshape(nu.shape)


def simulate(I, sigma, cond, N=1300, dt=2.5e-3, n_isi=13, Tcap=200.0,
             v_th=14.0, v_reset=-14.0, seed=0):
    """Return ISI samples for one (I, sigma, cond). cond in {white, ou, adapt}."""
    rng = np.random.default_rng(seed)
    v = np.full(N, float(v_reset)); last = np.zeros(N); isis = []
    w = np.zeros(N); xi = np.zeros(N)
    nu = I / sigma ** (4.0 / 3.0)
    Tmax = min(n_isi * sigma ** (-2.0 / 3.0) * J_of_nu(nu), Tcap)
    nsteps = int(Tmax / dt); sdt = np.sqrt(dt)
    adapt = (cond == "adapt"); ou = (cond == "ou")
    s_ou = sigma / np.sqrt(2.0 * TAU_S); ou_a = s_ou * np.sqrt(2.0 / TAU_S)
    for k in range(nsteps):
        t = k * dt
        drift = v * v + I
        if adapt:
            drift = drift - w
            w += -(w / TAU_W) * dt
        if ou:
            xi += -(xi / TAU_S) * dt + ou_a * sdt * rng.standard_normal(N)
            v += drift * dt + xi * dt
        else:
            v += drift * dt + sigma * sdt * rng.standard_normal(N)
        fired = v > v_th
        if fired.any():
            idx = np.where(fired)[0]
            isis.extend((t - last[idx]).tolist())
            last[idx] = t; v[fired] = v_reset
            if adapt:
                w[fired] += B_ADAPT
    isis = np.array(isis); isis = isis[isis > 0]
    return isis


def skew(x):
    x = np.asarray(x, float); m = x.mean(); sd = x.std()
    return float(np.mean(((x - m) / sd) ** 3)) if sd > 0 and len(x) > 2 else np.nan


def compute(cond):
    rows = []                                   # sigma, I, nu, r, cv, skew, n
    thr_isis = {}                               # sigma -> normalized ISIs at nu~0 (shape panel)
    for sg in SIGMAS:
        for nu in NUS:
            I = nu * sg ** (4.0 / 3.0)
            isis = simulate(I, sg, cond, seed=int(1000 * sg) + int(100 * nu) + hash(cond) % 97)
            if len(isis) < 5:
                rows.append((sg, I, nu, np.nan, np.nan, np.nan, len(isis))); continue
            r = 1.0 / isis.mean(); cv = isis.std() / isis.mean()
            rows.append((sg, I, nu, r, cv, skew(isis), len(isis)))
            if abs(nu) < 1e-9:
                thr_isis[sg] = isis / isis.mean()
    rows = np.array(rows, float)
    np.savez(os.path.join(CACHE, f"{cond}.npz"), rows=rows,
             thr05=thr_isis.get(0.5, np.array([])), thr03=thr_isis.get(0.3, np.array([])))
    print(f"[{cond}]  sigma     I      nu     r      CV     skew     n")
    for sg, I, nu, r, cv, sk, n in rows:
        print(f"        {sg:5.2f} {I:7.3f} {nu:6.2f} {r:6.3f} {cv:6.3f} {sk:6.2f}  {int(n):6d}")
    print(f"  saved {cond}.npz")


def plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    data = {}
    for cond in ("white", "ou", "adapt"):
        f = os.path.join(CACHE, f"{cond}.npz")
        if os.path.exists(f):
            data[cond] = np.load(f)
    cstyle = {"white": ("#0f172a", "ideal QIF (white)"),
              "ou":    ("#2563eb", fr"colored noise (OU $\tau_s={TAU_S}$)"),
              "adapt": ("#dc2626", fr"adaptation (sAHP $\tau_w={TAU_W}$)")}

    fig, ax = plt.subplots(1, 3, figsize=(15.6, 4.8))
    nu_grid = np.linspace(-0.6, 2.6, 120); Phi_grid = Phi(nu_grid)

    # Panel A: P1 -- the sigma^{2/3} collapse, ideal vs confounds
    ax[0].plot(nu_grid, Phi_grid, color="#16a34a", lw=2.6, zorder=1,
               label=r"$\Phi(\nu)=1/J(\nu)$ (prediction)")
    for cond, d in data.items():
        col, lab = cstyle[cond]; rows = d["rows"]
        ok = np.isfinite(rows[:, 3])
        ax[0].plot(rows[ok, 2], rows[ok, 3] * rows[ok, 0] ** (-2 / 3), "o", color=col, ms=7,
                   mec="white", mew=1.0, zorder=3, label=lab)
    ax[0].set_xlabel(r"$\nu = I/\sigma^{4/3}$"); ax[0].set_ylabel(r"$r\,\sigma^{-2/3}$")
    ax[0].set_title(r"P1: $\sigma^{2/3}$ f–I collapse")
    ax[0].set_ylim(0, 0.9); ax[0].legend(fontsize=8.4, framealpha=0.95, edgecolor="#cbd5e1")

    # Panel B: P2 -- the CV crossover, ideal vs confounds (average over the two sigmas)
    for cond, d in data.items():
        col, lab = cstyle[cond]; rows = d["rows"]
        cvm = []
        for nu in NUS:
            sel = rows[(np.abs(rows[:, 2] - nu) < 1e-6) & np.isfinite(rows[:, 4])]
            cvm.append(sel[:, 4].mean() if len(sel) else np.nan)
        ax[1].plot(NUS, cvm, "o-", color=col, ms=7, mec="white", mew=1.0, label=lab)
    ax[1].axhline(1.0, color="#94a3b8", lw=1.0, ls=":")
    ax[1].axhline(0.57, color="#16a34a", lw=1.2, ls="--")
    ax[1].text(1.7, 0.60, "SNIC 0.57", fontsize=8.5, color="#16a34a")
    ax[1].set_xlabel(r"$\nu = I/\sigma^{4/3}$"); ax[1].set_ylabel("ISI coefficient of variation")
    ax[1].set_title("P2: spike-train regularity crossover")
    ax[1].set_ylim(0, 1.1); ax[1].legend(fontsize=8.4, framealpha=0.95, edgecolor="#cbd5e1")

    # Panel C: P3 -- threshold ISI shape (white QIF = quartic FPT) vs exponential vs adapting
    bins = np.linspace(0, 3.2, 46)
    if "white" in data and len(data["white"]["thr05"]):
        xw = data["white"]["thr05"]
        ax[2].hist(xw, bins=bins, density=True, color="#0f172a", alpha=0.5,
                   label=fr"white QIF = quartic FPT (CV {xw.std():.2f}, skew {skew(xw):.2f})")
    xx = np.linspace(1e-3, 3.2, 200)
    ax[2].plot(xx, np.exp(-xx), "--", color="#7c3aed", lw=1.8, label="exponential (Poisson null)")
    if "adapt" in data and len(data["adapt"]["thr05"]):
        xa = data["adapt"]["thr05"]
        ax[2].hist(xa, bins=bins, density=True, histtype="step", lw=2.2, color="#dc2626",
                   label=fr"adapting QIF (CV {xa.std():.2f}, skew {skew(xa):.2f})")
    ax[2].set_xlabel(r"ISI $/\ \langle$ISI$\rangle$"); ax[2].set_ylabel("density")
    ax[2].set_title("P3: threshold ISI law (shape)")
    ax[2].set_xlim(0, 3.2); ax[2].legend(fontsize=8.0, framealpha=0.95, edgecolor="#cbd5e1")

    fig.tight_layout()
    out = os.path.abspath(os.path.join(HERE, "..", "figures", "qif_validation.png"))
    fig.savefig(out, dpi=140); print("saved", out)


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "compute":
        compute(sys.argv[2] if len(sys.argv) > 2 else "white")
    elif len(sys.argv) >= 2 and sys.argv[1] == "plot":
        plot()
    else:
        print("usage: qif_validation.py compute {white|ou|adapt} | plot")
