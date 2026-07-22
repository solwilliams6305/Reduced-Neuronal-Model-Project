"""
investigate_snic.py
-------------------
Dedicated fine investigation of the SNIC transition region.

Two parts:

  Part 1 — Fine I sweep around I_SNIC
    I values clustered tightly around the bifurcation, σ grid extended
    down to 0.005 to resolve the low boundary in the tonic regime.
    Fits α_CV for each I and tracks how it evolves through the transition.

  Part 2 — Scaling analysis at I = I_SNIC + 0.008 (just past SNIC)
    Fixed I just past bifurcation, fine (σ, ε) grid, fits α_CV precisely.
    Tests whether α_CV → 0.5 is a genuine result or a coarse-grid artefact.

Usage
-----
    python investigate_snic.py          # full investigation
    python investigate_snic.py --quick  # faster, fewer trajectories
"""

from __future__ import annotations
import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import linregress

import multiprocessing
multiprocessing.set_start_method("fork", force=True)

sys.path.insert(0, os.path.dirname(__file__))
from kernel import FHN2D, sweep_grid, failure_boundary

os.makedirs("figures/snic", exist_ok=True)
os.makedirs("data/snic",    exist_ok=True)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--quick",  action="store_true")
parser.add_argument("--n-traj", type=int, default=None)
parser.add_argument("--seed",   type=int, default=42)
parser.add_argument("--cv-threshold", type=float, default=0.15)
args = parser.parse_args()

ref    = FHN2D(I=-0.1, a=0.7, b=0.8)
I_SNIC = ref.I_snic
print(f"\nI_SNIC = {I_SNIC:.5f}")
print(f"CV threshold = {args.cv_threshold}")

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
if args.quick:
    n_traj  = args.n_traj or 600
    n_sigma = 12
    n_eps   = 8
    n_isi   = 40
    T_max   = 300.0
    tag     = "quick"
else:
    n_traj  = args.n_traj or 1500
    n_sigma = 20
    n_eps   = 15
    n_isi   = 60
    T_max   = 800.0
    tag     = "fine"

# σ grid extended down to 0.005 to capture low tonic boundary
sigma_vals = np.concatenate([
    np.linspace(0.005, 0.05,  8),    # fine low end
    np.linspace(0.06,  0.45, n_sigma - 8),  # coarser high end
])
sigma_vals = np.unique(np.round(sigma_vals, 4))

eps_vals = np.linspace(0.01, 0.20, n_eps)

print(f"\nGrid: {len(sigma_vals)}×{len(eps_vals)}")
print(f"σ: [{sigma_vals[0]:.4f}, {sigma_vals[-1]:.3f}]  ({len(sigma_vals)} pts)")
print(f"ε: [{eps_vals[0]:.3f}, {eps_vals[-1]:.3f}]  ({len(eps_vals)} pts)")
print(f"N_traj={n_traj}, N_isi={n_isi}, T_max={T_max}\n")

# ---------------------------------------------------------------------------
# I values — tight cluster around SNIC
# ---------------------------------------------------------------------------
PART1_I_VALS = [
    I_SNIC - 0.050,   # excitable, approaching
    I_SNIC - 0.020,   # excitable, very close
    I_SNIC - 0.008,   # excitable, extremely close
    I_SNIC + 0.004,   # tonic, just past
    I_SNIC + 0.008,   # tonic, just past (key point)
    I_SNIC + 0.015,   # tonic, near
    I_SNIC + 0.030,   # tonic, moderate
    I_SNIC + 0.060,   # tonic, further
    I_SNIC + 0.150,   # tonic, well past
]

# ---------------------------------------------------------------------------
# CV-based failure boundary
# ---------------------------------------------------------------------------

def cv_failure_boundary(
    result: dict,
    cv_threshold: float,
) -> tuple[np.ndarray, np.ndarray]:
    sv  = result["sigma_vals"]
    ev  = result["eps_vals"]
    cv  = result.get("cv", np.full((len(sv), len(ev)), np.nan))
    eps_bd, sig_bd = [], []
    for j in range(len(ev)):
        col   = cv[:, j]
        above = np.where(~np.isnan(col) & (col > cv_threshold))[0]
        if len(above) > 0:
            idx = above[0]
            if idx > 0 and not np.isnan(col[idx-1]):
                frac = (cv_threshold - col[idx-1]) / (col[idx] - col[idx-1])
                sig  = sv[idx-1] + frac * (sv[idx] - sv[idx-1])
            else:
                sig = sv[idx]
            eps_bd.append(ev[j])
            sig_bd.append(sig)
    return np.array(eps_bd), np.array(sig_bd)


def fit_powerlaw(eps: np.ndarray, sig: np.ndarray) -> tuple[float, float, float]:
    """Returns (alpha, C, r²)."""
    if len(eps) < 3:
        return np.nan, np.nan, np.nan
    mask = (eps > 0) & (sig > 0)
    if mask.sum() < 3:
        return np.nan, np.nan, np.nan
    res  = linregress(np.log(eps[mask]), np.log(sig[mask]))
    return float(res.slope), float(np.exp(res.intercept)), float(res.rvalue**2)


# ---------------------------------------------------------------------------
# Part 1 — Fine I sweep around SNIC
# ---------------------------------------------------------------------------
print("=" * 64)
print("PART 1 — Fine I sweep around I_SNIC")
print("=" * 64)

part1_results = []

for I_val in PART1_I_VALS:
    model = FHN2D(I=I_val, a=0.7, b=0.8)
    # Switch to isi_sequence when ΔU ≈ 0 even if technically excitable —
    # near the SNIC the barrier vanishes and first_passage trajectories
    # wander indefinitely, making the sweep extremely slow.
    if I_val >= I_SNIC or model.delta_U < 1e-4:
        mode = "isi_sequence"
    else:
        mode = "first_passage"
    print(f"\n  I={I_val:.5f}  ΔI={I_val-I_SNIC:+.5f}  "
          f"ΔU={model.delta_U:.6f}  mode={mode}")

    result = sweep_grid(
        model=model,
        sigma_vals=sigma_vals,
        eps_vals=eps_vals,
        mode=mode,
        n_trajectories=n_traj,
        n_isi=n_isi,
        T=T_max,
        seed=args.seed,
        verbose=False,
    )

    I_tag = f"I_{I_val:.5f}".replace("-","m").replace(".","p")
    np.savez(f"data/snic/{I_tag}_{tag}.npz",
             **{k: v for k, v in result.items() if isinstance(v, np.ndarray)})

    # Boundaries
    eps_lr, sig_lr = failure_boundary(result, delta_ratio=2.0)
    eps_cv, sig_cv = (cv_failure_boundary(result, args.cv_threshold)
                      if mode == "isi_sequence"
                      else (np.array([]), np.array([])))

    alpha_lr, C_lr, r2_lr = fit_powerlaw(eps_lr, sig_lr)
    alpha_cv, C_cv, r2_cv = fit_powerlaw(eps_cv, sig_cv)

    part1_results.append(dict(
        I=I_val, dI=I_val-I_SNIC, mode=mode,
        delta_U=model.delta_U,
        eps_lr=eps_lr, sig_lr=sig_lr,
        eps_cv=eps_cv, sig_cv=sig_cv,
        alpha_lr=alpha_lr, C_lr=C_lr, r2_lr=r2_lr,
        alpha_cv=alpha_cv, C_cv=C_cv, r2_cv=r2_cv,
        n_lr=len(eps_lr), n_cv=len(eps_cv),
    ))

    print(f"    LR: α={alpha_lr:.3f}, C={C_lr:.3f}, r²={r2_lr:.3f}  (n={len(eps_lr)})")
    print(f"    CV: α={alpha_cv:.3f}, C={C_cv:.3f}, r²={r2_cv:.3f}  (n={len(eps_cv)})")

# ---------------------------------------------------------------------------
# Part 1 figure — α_CV and α_LR vs ΔI
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

dI_vals    = [r["dI"]      for r in part1_results]
alpha_lr   = [r["alpha_lr"] for r in part1_results]
alpha_cv   = [r["alpha_cv"] for r in part1_results]
dU_vals    = [r["delta_U"]  for r in part1_results]

# Panel 1: exponent vs ΔI
ax = axes[0]
exc = [r for r in part1_results if r["mode"] == "first_passage" and not np.isnan(r["alpha_lr"])]
ton = [r for r in part1_results if r["mode"] == "isi_sequence"  and not np.isnan(r["alpha_cv"])]

if exc:
    ax.plot([r["dI"] for r in exc], [r["alpha_lr"] for r in exc],
            "o-", color="#2166ac", lw=2, ms=7, label="α_LR (excitable)")
if ton:
    ax.plot([r["dI"] for r in ton], [r["alpha_cv"] for r in ton],
            "s--", color="#d6604d", lw=2, ms=7, label="α_CV (tonic)")

ax.axvline(0, color="k", lw=1.0, ls=":", alpha=0.5, label="I_SNIC")
ax.axhline(0.5, color="gray", lw=1.2, ls="--", alpha=0.7, label="α=0.5")
ax.set_xlabel(r"$\Delta I = I - I_{SNIC}$", fontsize=11)
ax.set_ylabel(r"Power-law exponent $\alpha$", fontsize=11)
ax.set_title("Exponent evolution through SNIC", fontsize=11)
ax.legend(fontsize=8)
ax.set_ylim(0, 1.8)

# Panel 2: prefactor C vs ΔI
ax2 = axes[1]
if exc:
    ax2.semilogy([r["dI"] for r in exc], [r["C_lr"] for r in exc],
                 "o-", color="#2166ac", lw=2, ms=7, label="C_LR (excitable)")
if ton:
    ax2.semilogy([r["dI"] for r in ton], [r["C_cv"] for r in ton],
                 "s--", color="#d6604d", lw=2, ms=7, label="C_CV (tonic)")
ax2.axvline(0, color="k", lw=1.0, ls=":", alpha=0.5)
ax2.set_xlabel(r"$\Delta I$", fontsize=11)
ax2.set_ylabel("Prefactor $C$  (log scale)", fontsize=11)
ax2.set_title("Prefactor evolution through SNIC", fontsize=11)
ax2.legend(fontsize=8)

# Panel 3: boundary curves overlaid
ax3 = axes[2]
cmap  = plt.cm.RdYlBu_r
dI_arr = np.array([r["dI"] for r in part1_results])
vmin, vmax = dI_arr.min(), dI_arr.max()
eps_ref = np.linspace(eps_vals.min(), eps_vals.max(), 200)

for r in part1_results:
    col = cmap((r["dI"] - vmin) / max(vmax - vmin, 1e-6))
    if r["mode"] == "first_passage" and len(r["eps_lr"]) >= 2:
        ax3.plot(r["eps_lr"], r["sig_lr"], "-", color=col, lw=1.8,
                 label=f"$\\Delta I={r['dI']:+.3f}$ (exc)")
    elif len(r["eps_cv"]) >= 2:
        ax3.plot(r["eps_cv"], r["sig_cv"], "--", color=col, lw=1.8,
                 label=f"$\\Delta I={r['dI']:+.3f}$ (ton)")

ax3.plot(eps_ref, 0.5*np.sqrt(eps_ref), "k:", lw=1.2, alpha=0.5,
         label=r"$0.5\sqrt{\varepsilon}$")
ax3.set_xlabel(r"$\varepsilon$", fontsize=11)
ax3.set_ylabel(r"$\sigma^*$", fontsize=11)
ax3.set_title("Failure boundaries near SNIC", fontsize=11)
ax3.legend(fontsize=6, ncol=2)
ax3.set_ylim(0, sigma_vals.max())

fig.suptitle(f"SNIC transition  ($I_{{SNIC}}={I_SNIC:.4f}$)", fontsize=13)
fig.tight_layout()
fig.savefig("figures/snic/part1_snic_transition.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"\nPart 1 saved: figures/snic/part1_snic_transition.png")

# ---------------------------------------------------------------------------
# Part 2 — Scaling analysis at fixed I just past SNIC
# ---------------------------------------------------------------------------
print("\n" + "=" * 64)
print("PART 2 — Precise scaling at I = I_SNIC + 0.008")
print("=" * 64)

I_focus = I_SNIC + 0.008
model2  = FHN2D(I=I_focus, a=0.7, b=0.8)
print(f"\n  I={I_focus:.5f}, ΔU={model2.delta_U:.6f}")

# Finer grids for precise exponent measurement
sigma_fine = np.concatenate([
    np.linspace(0.005, 0.05,  12),
    np.linspace(0.06,  0.30,  10),
])
sigma_fine = np.unique(np.round(sigma_fine, 4))
eps_fine   = np.logspace(np.log10(0.008), np.log10(0.20), 18)  # log-spaced

result2 = sweep_grid(
    model=model2,
    sigma_vals=sigma_fine,
    eps_vals=eps_fine,
    mode="isi_sequence",
    n_trajectories=n_traj,
    n_isi=n_isi,
    T=T_max,
    seed=args.seed + 99,
    verbose=True,
)

np.savez(f"data/snic/focus_I{I_focus:.5f}_{tag}.npz",
         **{k: v for k, v in result2.items() if isinstance(v, np.ndarray)})

eps_cv2, sig_cv2 = cv_failure_boundary(result2, args.cv_threshold)
alpha2, C2, r2_2 = fit_powerlaw(eps_cv2, sig_cv2)

print(f"\n  CV boundary: {len(eps_cv2)} points")
print(f"  Fit: σ = {C2:.4f} · ε^{alpha2:.4f}  (r²={r2_2:.4f})")
print(f"  Predicted: α → 0.5 at SNIC")
print(f"  Gap from 0.5: {abs(alpha2 - 0.5):.4f}")

# Sliding window local exponent
if len(eps_cv2) >= 6:
    window = 4
    eps_c, alphas_local = [], []
    for k in range(len(eps_cv2) - window + 1):
        ew = eps_cv2[k:k+window]
        sw = sig_cv2[k:k+window]
        a, _, _ = fit_powerlaw(ew, sw)
        eps_c.append(np.sqrt(ew[0]*ew[-1]))
        alphas_local.append(a)
    eps_c      = np.array(eps_c)
    alphas_local = np.array(alphas_local)
    print(f"\n  Local exponents (window={window}):")
    for ec, al in zip(eps_c, alphas_local):
        print(f"    ε={ec:.4f}  α={al:.4f}")

# Part 2 figure
fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5))

# Left: CV map
ax = axes2[0]
cv2 = result2.get("cv", np.full((len(sigma_fine), len(eps_fine)), np.nan))
vmax_cv = min(np.nanpercentile(cv2, 98), 1.0)
im = ax.pcolormesh(eps_fine, sigma_fine, cv2, cmap="viridis_r",
                   vmin=0, vmax=max(vmax_cv, 0.1), shading="auto")
fig2.colorbar(im, ax=ax).set_label("CV", fontsize=9)
try:
    ax.contour(eps_fine, sigma_fine, cv2, levels=[args.cv_threshold],
               colors=["white"], linewidths=1.5)
except Exception:
    pass
if len(eps_cv2) >= 2:
    ax.plot(eps_cv2, sig_cv2, "s--", color="#984ea3", ms=5, lw=1.8,
            label=f"CV={args.cv_threshold} boundary")
    ax.legend(fontsize=8)
ax.set_xlabel(r"$\varepsilon$  (log-spaced)", fontsize=11)
ax.set_ylabel(r"$\sigma$", fontsize=11)
ax.set_title(f"CV map at $I={I_focus:.4f}$  ($\\Delta I_{{SNIC}}=+{I_focus-I_SNIC:.4f}$)",
             fontsize=10)

# Right: log-log boundary fit
ax2b = axes2[1]
if len(eps_cv2) >= 3:
    ax2b.loglog(eps_cv2, sig_cv2, "s", color="#984ea3", ms=8, zorder=5,
                label=f"CV boundary  (n={len(eps_cv2)})")
    eps_plot = np.logspace(np.log10(eps_cv2.min()*0.8),
                           np.log10(eps_cv2.max()*1.2), 200)
    if not np.isnan(alpha2):
        ax2b.loglog(eps_plot, C2 * eps_plot**alpha2, "-",
                    color="#984ea3", lw=2,
                    label=f"Fit: $\\sigma={C2:.3f}\\varepsilon^{{{alpha2:.3f}}}$  "
                          f"($r^2={r2_2:.3f}$)")
    ax2b.loglog(eps_plot, 0.5*eps_plot**0.5, "k--", lw=1.5, alpha=0.6,
                label=r"$0.5\,\varepsilon^{0.5}$  (predicted)")

    if len(eps_c) >= 2:
        ax2b_twin = ax2b.twinx()
        ax2b_twin.semilogx(eps_c, alphas_local, "o-",
                           color="#fdae61", lw=1.5, ms=5, alpha=0.8)
        ax2b_twin.axhline(0.5, color="gray", lw=1.0, ls="--", alpha=0.5)
        ax2b_twin.set_ylabel("Local exponent α", fontsize=9, color="#fdae61")
        ax2b_twin.tick_params(axis="y", labelcolor="#fdae61")
        ax2b_twin.set_ylim(0, 1.5)

ax2b.set_xlabel(r"$\varepsilon$  (log scale)", fontsize=11)
ax2b.set_ylabel(r"$\sigma^*$  (log scale)", fontsize=11)
ax2b.set_title(
    f"Log-log scaling at $I_{{SNIC}}+{I_focus-I_SNIC:.4f}$\n"
    f"Global fit: $\\alpha={alpha2:.4f}$,  predicted: $\\alpha=0.5$",
    fontsize=10)
ax2b.legend(fontsize=8)

fig2.suptitle(
    f"Part 2: Precise CV scaling just past SNIC  ($\\Delta I=+{I_focus-I_SNIC:.4f}$)",
    fontsize=12)
fig2.tight_layout()
fig2.savefig("figures/snic/part2_scaling_focus.png", dpi=150, bbox_inches="tight")
plt.close(fig2)
print(f"\nPart 2 saved: figures/snic/part2_scaling_focus.png")

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------
print(f"\n{'='*72}")
print(f"  SUMMARY")
print(f"  I_SNIC = {I_SNIC:.5f}")
print(f"\n  {'I':>8}  {'ΔI':>8}  {'ΔU':>10}  {'mode':>12}  "
      f"{'α_LR':>7}  {'α_CV':>7}  {'n_LR':>5}  {'n_CV':>5}")
print(f"  {'-'*72}")
for r in part1_results:
    print(f"  {r['I']:8.5f}  {r['dI']:+8.5f}  {r['delta_U']:10.6f}  "
          f"{r['mode']:>12}  "
          f"{r['alpha_lr']:7.3f}  {r['alpha_cv']:7.3f}  "
          f"{r['n_lr']:5d}  {r['n_cv']:5d}")
print(f"\n  Focus point (Part 2):")
print(f"  I={I_focus:.5f},  α_CV={alpha2:.4f},  C={C2:.4f},  r²={r2_2:.4f}")
print(f"{'='*72}\n")
