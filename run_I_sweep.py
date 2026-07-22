"""
run_I_sweep.py
--------------
Sweeps the external current I across the excitable → tonic spiking transition,
running a full (sigma, eps) regime map at each I value.

For each I:
  - I < I_SNIC  : mode='first_passage', diagnostic = log(T_drift / MFPT)
  - I > I_SNIC  : mode='isi_sequence',  diagnostic = CV + log(T_cycle / ISI_mean)

Produces:
  data/I_sweep/   — .npz files per I value
  figures/I_sweep/ — regime map per I value + summary figure

Usage
-----
    python run_I_sweep.py                  # default I values, fine grid
    python run_I_sweep.py --quick          # coarse grid, fast check
    python run_I_sweep.py --n-traj 1000   # override trajectories
    python run_I_sweep.py --I-vals -0.1 0.3 0.5 0.8   # custom I values
"""

from __future__ import annotations
import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

import multiprocessing
multiprocessing.set_start_method("fork", force=True)

sys.path.insert(0, os.path.dirname(__file__))
from kernel import (FHN2D, sweep_grid, failure_boundary,
                    reduced_prediction, log_ratio)

os.makedirs("figures/I_sweep", exist_ok=True)
os.makedirs("data/I_sweep",    exist_ok=True)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--quick",   action="store_true")
parser.add_argument("--n-traj",  type=int,   default=None)
parser.add_argument("--n-isi",   type=int,   default=40)
parser.add_argument("--seed",    type=int,   default=42)
parser.add_argument("--delta-ratio", type=float, default=2.0)
parser.add_argument("--cv-threshold", type=float, default=0.15,
                    help="CV threshold for tonic failure boundary (default 0.15)")
parser.add_argument("--I-vals",  type=float, nargs="+", default=None)
args = parser.parse_args()

# ---------------------------------------------------------------------------
# I values to sweep
# ---------------------------------------------------------------------------
ref    = FHN2D(I=-0.1, a=0.7, b=0.8)
I_SNIC = ref.I_snic
print(f"\nReference model: {ref}")
print(f"I_SNIC ≈ {I_SNIC:.4f}  (fixed point annihilates at left fold)")

if args.I_vals is not None:
    I_VALUES = args.I_vals
else:
    I_VALUES = [
        -0.10,            # deep excitable (baseline)
         0.15,            # excitable, moderate
         0.22,            # excitable, approaching SNIC
         0.26,            # excitable, close
         0.28,            # excitable, very close
         0.30,            # just below SNIC
         I_SNIC + 0.03,   # just above SNIC
         I_SNIC + 0.06,   # tonic, near
         0.44,            # tonic, mid
         0.69,            # tonic, robust
    ]

# ---------------------------------------------------------------------------
# Grid parameters
# ---------------------------------------------------------------------------
if args.quick:
    sigma_vals = np.linspace(0.02, 0.45, 10)
    eps_vals   = np.linspace(0.01, 0.20, 8)
    T_max      = 200.0
    n_traj     = args.n_traj or 400
    tag        = "quick"
else:
    sigma_vals = np.linspace(0.02, 0.45, 20)
    eps_vals   = np.linspace(0.01, 0.20, 15)
    T_max      = 600.0
    n_traj     = args.n_traj or 1000
    tag        = "fine"

print(f"\nGrid: {len(sigma_vals)}×{len(eps_vals)},  "
      f"sigma=[{sigma_vals[0]:.2f},{sigma_vals[-1]:.2f}],  "
      f"eps=[{eps_vals[0]:.3f},{eps_vals[-1]:.3f}]")
print(f"N_traj={n_traj},  N_isi={args.n_isi},  T_max={T_max},  tag={tag}")
print(f"CV threshold (tonic failure): {args.cv_threshold}")
print(f"\nI values: {[f'{v:.4f}' for v in I_VALUES]}\n")

# ---------------------------------------------------------------------------
# Helper: classify I
# ---------------------------------------------------------------------------

def classify(I: float) -> tuple[str, str]:
    return ("first_passage", "excitable") if I < I_SNIC else ("isi_sequence", "tonic")


# ---------------------------------------------------------------------------
# Helper: CV-based failure boundary for tonic mode
# ---------------------------------------------------------------------------

def cv_failure_boundary(
    result: dict,
    cv_threshold: float = 0.15,
) -> tuple[np.ndarray, np.ndarray]:
    """
    For tonic mode: find smallest sigma where CV > cv_threshold for each eps.
    Returns (eps_boundary, sigma_boundary).
    """
    sigma_vals = result["sigma_vals"]
    eps_vals   = result["eps_vals"]
    cv_arr     = result.get("cv", np.full((len(sigma_vals), len(eps_vals)), np.nan))

    eps_bd, sig_bd = [], []
    for j, eps in enumerate(eps_vals):
        col   = cv_arr[:, j]
        above = np.where(~np.isnan(col) & (col > cv_threshold))[0]
        if len(above) > 0:
            idx = above[0]
            if idx > 0 and not np.isnan(col[idx - 1]):
                frac = (cv_threshold - col[idx-1]) / (col[idx] - col[idx-1])
                sig  = sigma_vals[idx-1] + frac*(sigma_vals[idx]-sigma_vals[idx-1])
            else:
                sig = sigma_vals[idx]
            eps_bd.append(eps)
            sig_bd.append(sig)

    return np.array(eps_bd), np.array(sig_bd)


# ---------------------------------------------------------------------------
# Helper: plot single regime map
# ---------------------------------------------------------------------------

def plot_regime_map_I(
    result: dict,
    model:  FHN2D,
    mode:   str,
    eps_bd_lr: np.ndarray,   # log-ratio boundary
    sig_bd_lr: np.ndarray,
    eps_bd_cv: np.ndarray,   # CV boundary (tonic only)
    sig_bd_cv: np.ndarray,
    save_path: str,
) -> None:
    sv  = result["sigma_vals"]
    ev  = result["eps_vals"]
    lr  = result["log_ratio"]
    cv  = result.get("cv", np.full_like(lr, np.nan))

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # --- Left: log-ratio heatmap ---
    ax = axes[0]
    vmax = np.nanpercentile(np.abs(lr), 95)
    vmax = max(vmax, 0.5)
    norm = TwoSlopeNorm(vmin=-vmax/2, vcenter=0, vmax=vmax)
    im   = ax.pcolormesh(ev, sv, lr, cmap="RdBu_r", norm=norm, shading="auto")
    fig.colorbar(im, ax=ax, pad=0.02).set_label(
        r"$\log(T_{\rm ref}\,/\,{\rm stat}_{\rm full})$", fontsize=9)

    try:
        ax.contour(ev, sv, lr, levels=[0],
                   colors=["k"], linewidths=0.8, linestyles="--", alpha=0.5)
    except Exception:
        pass

    # Log-ratio boundary
    if len(eps_bd_lr) > 0:
        ax.plot(eps_bd_lr, sig_bd_lr, "o-", color="#1a9641",
                ms=5, lw=1.8, zorder=5,
                label=f"LR boundary (×{args.delta_ratio:.0f})")

    # CV boundary on left panel too (tonic only)
    if mode == "isi_sequence" and len(eps_bd_cv) > 0:
        ax.plot(eps_bd_cv, sig_bd_cv, "s--", color="#984ea3",
                ms=5, lw=1.8, zorder=5,
                label=f"CV boundary (>{args.cv_threshold})")

    # √ε reference fit
    eps_ref = np.linspace(ev.min(), ev.max(), 200)
    all_bd_sig = np.concatenate([sig_bd_lr, sig_bd_cv]) if len(sig_bd_cv) > 0 else sig_bd_lr
    all_bd_eps = np.concatenate([eps_bd_lr, eps_bd_cv]) if len(eps_bd_cv) > 0 else eps_bd_lr
    if len(all_bd_eps) >= 2:
        C_fit = float(np.median(all_bd_sig / np.sqrt(all_bd_eps)))
        sig_ref = C_fit * np.sqrt(eps_ref)
        mask = (sig_ref >= sv.min()) & (sig_ref <= sv.max())
        ax.plot(eps_ref[mask], sig_ref[mask], "--",
                color="#fdae61", lw=1.5, alpha=0.7,
                label=f"$\\sigma={C_fit:.2f}\\sqrt{{\\varepsilon}}$  (ref)")

    ax.legend(fontsize=7, loc="upper left")
    ax.set_xlabel(r"$\varepsilon$", fontsize=11)
    ax.set_ylabel(r"$\sigma$", fontsize=11)
    ax.set_title(
        f"$I={model.I:.3f}$  ({model.regime})  —  log-ratio\n"
        f"$\\Delta I_{{SNIC}} = {model.I-I_SNIC:+.3f}$,  "
        f"$\\Delta U={model.delta_U:.5f}$",
        fontsize=10)

    # --- Right: CV map ---
    ax2 = axes[1]
    if np.all(np.isnan(cv)):
        ax2.text(0.5, 0.5, "CV not available\n(excitable mode)",
                 ha="center", va="center", transform=ax2.transAxes, fontsize=12)
        ax2.set_title("CV — not applicable (excitable)", fontsize=10)
    else:
        vmax_cv = min(np.nanpercentile(cv, 98), 1.0)
        vmax_cv = max(vmax_cv, 0.1)
        im2 = ax2.pcolormesh(ev, sv, cv, cmap="viridis_r",
                             vmin=0, vmax=vmax_cv, shading="auto")
        fig.colorbar(im2, ax=ax2, pad=0.02).set_label("CV", fontsize=9)

        # CV threshold contour
        try:
            ax2.contour(ev, sv, cv, levels=[args.cv_threshold],
                        colors=["white"], linewidths=1.5)
        except Exception:
            pass

        if len(eps_bd_cv) > 0:
            ax2.plot(eps_bd_cv, sig_bd_cv, "s--", color="#984ea3",
                     ms=5, lw=1.8, zorder=5,
                     label=f"CV={args.cv_threshold} boundary")
            ax2.legend(fontsize=7)

        ax2.set_title(
            f"CV map  ($I={model.I:.3f}$)\n"
            f"CV=0: perfect clock,  CV=1: Poisson,  "
            f"threshold={args.cv_threshold}",
            fontsize=9)

    ax2.set_xlabel(r"$\varepsilon$", fontsize=11)
    ax2.set_ylabel(r"$\sigma$", fontsize=11)

    fig.suptitle(
        f"Regime map  |  $I={model.I:.3f}$  ({model.regime})"
        f"  |  $\\Delta I_{{SNIC}}={model.I-I_SNIC:+.3f}$",
        fontsize=12, y=1.01)
    fig.tight_layout()
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {save_path}")


# ---------------------------------------------------------------------------
# Helper: summary figure
# ---------------------------------------------------------------------------

def plot_summary(all_results: list[dict], save_path: str) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    cmap   = plt.cm.RdYlBu_r
    I_list = [r["params"]["I"] for r in all_results]
    vmin_I = min(I_list)
    vmax_I = max(I_list)

    eps_ref   = np.linspace(eps_vals.min(), eps_vals.max(), 200)
    sigma_ref = 0.8 * np.sqrt(eps_ref)

    for ax_idx, (ax, bd_key) in enumerate(zip(axes, ["lr", "cv"])):
        ax.plot(eps_ref, sigma_ref, "k:", lw=1.2, alpha=0.4,
                label=r"$\sigma\sim\sqrt{\varepsilon}$")

        for res, bsum in zip(all_results, boundary_summary):
            I_val  = res["params"]["I"]
            mode   = res["params"]["mode"]
            col    = cmap((I_val - vmin_I) / max(vmax_I - vmin_I, 1e-6))
            label  = f"$I={I_val:.2f}$ ({I_val-I_SNIC:+.2f})"
            ls     = "-" if mode == "first_passage" else "--"

            if bd_key == "lr":
                eps_bd, sig_bd = bsum["eps_bd_lr"], bsum["sig_bd_lr"]
            else:
                eps_bd, sig_bd = bsum["eps_bd_cv"], bsum["sig_bd_cv"]

            if len(eps_bd) >= 2:
                ax.plot(eps_bd, sig_bd, ls, color=col, lw=2.0, label=label)
            elif len(eps_bd) == 1:
                ax.plot(eps_bd, sig_bd, "o", color=col, ms=6)

        ax.axvline(0, color="gray", lw=0.5)
        ax.set_xlabel(r"$\varepsilon$", fontsize=12)
        ax.set_ylabel(r"$\sigma^*$", fontsize=12)
        title = ("Log-ratio failure boundary" if bd_key == "lr"
                 else f"CV>{args.cv_threshold} failure boundary")
        ax.set_title(
            f"{title}\nSolid: excitable,  Dashed: tonic  |  "
            f"$I_{{SNIC}}\\approx{I_SNIC:.3f}$",
            fontsize=10)
        ax.legend(fontsize=7, loc="upper left", ncol=2)
        ax.set_xlim(eps_vals.min(), eps_vals.max())
        ax.set_ylim(sigma_vals.min(), sigma_vals.max())

    fig.tight_layout()
    fig.savefig(save_path, dpi=150)
    plt.close(fig)
    print(f"\nSummary saved: {save_path}")


# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------

all_results    = []
boundary_summary = []

print("=" * 64)
for I_val in I_VALUES:
    model = FHN2D(I=I_val, a=0.7, b=0.8)
    mode, regime = classify(I_val)

    print(f"\nI = {I_val:.4f}  |  regime={regime}  |  mode={mode}")
    print(f"  {model}")

    result = sweep_grid(
        model=model,
        sigma_vals=sigma_vals,
        eps_vals=eps_vals,
        mode=mode,
        n_trajectories=n_traj,
        n_isi=args.n_isi,
        T=T_max,
        seed=args.seed,
        verbose=True,
    )

    # Save
    I_tag    = f"I_{I_val:.4f}".replace("-", "m").replace(".", "p")
    npz_path = f"data/I_sweep/{I_tag}_{tag}.npz"
    np.savez(npz_path,
             **{k: v for k, v in result.items() if isinstance(v, np.ndarray)})
    print(f"  Saved: {npz_path}")

    # Failure boundaries
    eps_bd_lr, sig_bd_lr = failure_boundary(result, delta_ratio=args.delta_ratio)
    eps_bd_cv, sig_bd_cv = (cv_failure_boundary(result, args.cv_threshold)
                             if mode == "isi_sequence"
                             else (np.array([]), np.array([])))

    print(f"  LR boundary pts: {len(eps_bd_lr)},  "
          f"CV boundary pts: {len(eps_bd_cv)}")

    # Plot
    fig_path = f"figures/I_sweep/regime_map_{I_tag}_{tag}.png"
    plot_regime_map_I(result, model, mode,
                      eps_bd_lr, sig_bd_lr,
                      eps_bd_cv, sig_bd_cv,
                      fig_path)

    all_results.append(result)
    boundary_summary.append({
        "I": I_val, "regime": regime,
        "eps_bd_lr": eps_bd_lr, "sig_bd_lr": sig_bd_lr,
        "eps_bd_cv": eps_bd_cv, "sig_bd_cv": sig_bd_cv,
    })

    # Print fit
    for label, ebd, sbd in [("LR", eps_bd_lr, sig_bd_lr),
                              ("CV", eps_bd_cv, sig_bd_cv)]:
        if len(ebd) >= 3:
            alpha, log_C = np.polyfit(np.log(ebd), np.log(sbd), 1)
            print(f"  {label} fit: σ = {np.exp(log_C):.3f} · ε^{alpha:.3f}")

print("\n" + "=" * 64)

# Summary figure
plot_summary(all_results, "figures/I_sweep/boundary_summary_fine.png")

# Summary table
print(f"\n{'I':>7}  {'regime':>10}  {'ΔI':>7}  "
      f"{'α_LR':>7}  {'α_CV':>7}  {'n_LR':>5}  {'n_CV':>5}")
print("-" * 60)
for b in boundary_summary:
    def fit(e, s):
        if len(e) >= 3:
            a, lc = np.polyfit(np.log(e), np.log(s), 1)
            return f"{a:.3f}"
        return "  nan"
    print(f"  {b['I']:5.3f}  {b['regime']:>10}  "
          f"{b['I']-I_SNIC:+7.3f}  "
          f"{fit(b['eps_bd_lr'], b['sig_bd_lr']):>7}  "
          f"{fit(b['eps_bd_cv'], b['sig_bd_cv']):>7}  "
          f"{len(b['eps_bd_lr']):5d}  "
          f"{len(b['eps_bd_cv']):5d}")

print("\nDone.")
