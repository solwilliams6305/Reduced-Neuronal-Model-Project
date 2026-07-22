"""
bifurcation_inverter.py — MVP rung 1 of the moonshot: invert the edge theory to read the rescaled
distance-to-threshold nu off a spike train, with NO knowledge of the noise level or time-scale.

Key identifiability fact: the MEAN inter-event interval carries sigma and the (unknown) time-scale,
but the DIMENSIONLESS SHAPE of the interval distribution — CV and skew — is a sigma-free, monotone
function of nu = (distance past threshold)/sigma^{4/3}.  So shape statistics pin nu directly; the
early-warning signal is "nu -> 0" (approaching the SNIC/Type-I threshold).

  atlas   : nu -> (CV, skew)   from the canonical inner FPT  drho=(nu+rho^2)ds+dW   (offcritical_phase_edge.fpt_sim)
  invert  : ISIs -> nu_hat     by matching (CV*, skew*) to the atlas, + bootstrap CI + fit residual (class check)
  validate: physical QIF spike trains at known nu -> recover nu_hat
            (white at two sigmas: recovery should be sigma-invariant; adaptation: a known bias)

    python3 bifurcation_inverter.py atlas      # build + cache the SNIC atlas (one-time)
    python3 bifurcation_inverter.py validate   # -> figures/bifurcation_inverter.png
"""
from __future__ import annotations
import os
import sys
import numpy as np
from offcritical_phase_edge import fpt_sim, J_of_nu
from qif_validation import simulate as qif_isis        # physical QIF ISIs (white/ou/adapt)

HERE = os.path.dirname(os.path.abspath(__file__))
ATLAS = os.path.join(HERE, "_inverter_atlas.npz")


def skew(x):
    x = np.asarray(x, float); m = x.mean(); s = x.std()
    return float(np.mean(((x - m) / s) ** 3)) if s > 0 and len(x) > 2 else np.nan


# ---- atlas: nu -> (CV, skew) from the universal inner FPT -----------------------------------------
def build_atlas(nu_min=-0.9, nu_max=3.0, ngrid=16, N=4500):
    nus = np.linspace(nu_min, nu_max, ngrid); CV = []; SK = []
    for nu in nus:
        Te, m, sd, sk = fpt_sim(float(nu), N=N, seed=7)
        CV.append(sd / m); SK.append(sk)
        print(f"  atlas nu={nu:+.3f}  CV={sd/m:.3f}  skew={sk:.3f}")
    np.savez(ATLAS, nus=nus, CV=np.array(CV), SK=np.array(SK))
    print(f"  saved {ATLAS}")


def _load_atlas():
    d = np.load(ATLAS)
    nf = np.linspace(d["nus"][0], d["nus"][-1], 600)
    return nf, np.interp(nf, d["nus"], d["CV"]), np.interp(nf, d["nus"], d["SK"])


def invert(isis, w_skew=0.5):
    """ISIs -> (nu_hat, residual). residual is the min normalised shape-mismatch (a SNIC class check)."""
    isis = np.asarray(isis, float); isis = isis[isis > 0]
    if len(isis) < 20:
        return np.nan, np.nan
    cv = isis.std() / isis.mean(); sk = skew(isis)
    nf, CVf, SKf = _load_atlas()
    sCV = CVf.std(); sSK = SKf.std()
    cost = ((CVf - cv) / sCV) ** 2 + w_skew * ((SKf - sk) / sSK) ** 2
    j = int(np.argmin(cost))
    return float(nf[j]), float(np.sqrt(cost[j]))


def invert_ci(isis, nboot=120, seed=0):
    rng = np.random.default_rng(seed); isis = np.asarray(isis, float)
    base, res = invert(isis)
    bs = []
    for _ in range(nboot):
        s = rng.choice(isis, size=len(isis), replace=True)
        nu, _ = invert(s)
        if np.isfinite(nu):
            bs.append(nu)
    lo, hi = np.percentile(bs, [16, 84]) if bs else (np.nan, np.nan)
    return base, lo, hi, res


def validate():
    if not os.path.exists(ATLAS):
        build_atlas()
    true_nus = np.array([-0.3, 0.0, 0.5, 1.0, 1.8])
    runs = [("white", 0.4, "#2563eb", "QIF white  σ=0.4"),
            ("white", 0.7, "#16a34a", "QIF white  σ=0.7"),
            ("adapt", 0.4, "#dc2626", "QIF adapt  σ=0.4")]
    results = {}
    for cond, sg, col, lab in runs:
        out = []
        for nu in true_nus:
            I = nu * sg ** (4.0 / 3.0)
            isis = qif_isis(I, sg, cond, seed=int(500 * sg) + int(100 * nu))
            nu_hat, lo, hi, res = invert_ci(isis)
            out.append((nu, nu_hat, lo, hi, res))
            print(f"  [{cond:5s} σ={sg}] true ν={nu:+.2f}  ->  ν̂={nu_hat:+.2f} "
                  f"[{lo:+.2f},{hi:+.2f}]  resid={res:.2f}")
        results[(cond, sg)] = (np.array(out), col, lab)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.9))
    d = np.load(ATLAS)

    # Panel A: the atlas -- the invertible (monotone) shape-vs-nu map
    ax[0].plot(d["nus"], d["CV"], "o-", color="#0f172a", ms=5, label="CV(ν)")
    ax[0].plot(d["nus"], d["SK"], "s--", color="#7c3aed", ms=5, label="skew(ν)")
    ax[0].axvline(0, color="#94a3b8", lw=1.0, ls=":"); ax[0].text(0.03, 0.05,
              "threshold", rotation=90, fontsize=8, color="#475569", transform=ax[0].get_xaxis_transform())
    ax[0].set_xlabel(r"$\nu$ (rescaled distance to threshold)"); ax[0].set_ylabel("dimensionless shape")
    ax[0].set_title("The SNIC atlas (σ-free, monotone)")
    ax[0].legend(fontsize=9, framealpha=0.95, edgecolor="#cbd5e1")

    # Panel B: recovered vs true nu
    lim = (-0.7, 2.2)
    ax[1].plot(lim, lim, "--", color="#94a3b8", lw=1.2, label="ideal")
    for (cond, sg), (out, col, lab) in results.items():
        yerr = np.abs(out[:, 1] - out[:, 2:4].T)
        mk = "o" if cond == "white" else "D"
        ax[1].errorbar(out[:, 0], out[:, 1], yerr=yerr, fmt=mk, color=col, ms=7, mec="white",
                       mew=1.0, capsize=3, lw=1.0, label=lab)
    ax[1].set_xlim(lim); ax[1].set_ylim(lim)
    ax[1].set_xlabel(r"true $\nu$"); ax[1].set_ylabel(r"recovered $\hat\nu$ (from ISI shape only)")
    ax[1].set_title("Inversion: ν from shape, no σ needed")
    ax[1].legend(fontsize=8.4, framealpha=0.95, edgecolor="#cbd5e1")

    fig.tight_layout()
    out = os.path.abspath(os.path.join(HERE, "..", "figures", "bifurcation_inverter.png"))
    fig.savefig(out, dpi=140); print("saved", out)


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "atlas":
        build_atlas()
    elif len(sys.argv) >= 2 and sys.argv[1] == "validate":
        validate()
    else:
        print("usage: bifurcation_inverter.py atlas | validate")
