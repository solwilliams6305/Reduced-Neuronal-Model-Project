"""
run_diagnostics.py
------------------
Three-part investigation of the 0.43 vs 0.5 exponent gap and
Kramers/fold-layer crossover:

  Step 1 — Overlay sigma_cross = sqrt(DeltaU / log(1/eps)) on regime map
  Step 2 — Fine sweep: log-spaced eps down to 0.001, sigma down to 0.01
  Step 3 — w_escape heatmap: empirical crossover between mechanisms

Usage
-----
    python run_diagnostics.py           # standard (20x15)
    python run_diagnostics.py --quick   # coarse   (8x6)
    python run_diagnostics.py --fine    # fine, log-spaced eps (30x20)
"""

import argparse, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# macOS uses 'spawn' by default which causes Pool to re-import __main__
# and crash. 'fork' avoids this and is safe here (no GUI before fork).
import multiprocessing
multiprocessing.set_start_method("fork", force=True)

from simulate import simulate_fhn, potential_barrier, w_fold_left, W_FP, W_FP
from sweep   import (sweep_grid, failure_boundary, fit_power_law,
                     local_exponent, kramers_crossover_curve)
from plot    import (plot_regime_map, plot_boundary_scaling,
                     plot_local_exponent, plot_w_spike_map)

os.makedirs("figures", exist_ok=True)
os.makedirs("data",    exist_ok=True)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--quick", action="store_true")
parser.add_argument("--fine",  action="store_true")
parser.add_argument("--n-traj", type=int,   default=None)
parser.add_argument("--seed",   type=int,   default=42)
parser.add_argument("--delta-ratio", type=float, default=2.0)
args = parser.parse_args()

if args.quick:
    sigma_vals = np.linspace(0.05, 0.50, 8)
    eps_vals   = np.linspace(0.01, 0.25, 6)
    n_traj     = args.n_traj or 400
    T          = 200.0
    tag        = "quick"
elif args.fine:
    # Log-spaced eps: probes fold-layer regime at very small eps
    # Small sigma: probes below Kramers crossover
    sigma_vals = np.concatenate([
        np.linspace(0.01, 0.08, 10),   # below crossover — fold-layer regime
        np.linspace(0.10, 0.50, 20),   # above crossover — Kramers regime
    ])
    eps_vals   = np.logspace(-3, np.log10(0.25), 20)   # 0.001 to 0.25
    n_traj     = args.n_traj or 800
    T          = 800.0    # longer T needed at small eps/sigma
    tag        = "fine"
else:
    sigma_vals = np.linspace(0.05, 0.50, 20)
    eps_vals   = np.linspace(0.01, 0.25, 15)
    n_traj     = args.n_traj or 600
    T          = 300.0
    tag        = "standard"

# ---------------------------------------------------------------------------
# Compute Delta_U analytically at the fixed point w value
# ---------------------------------------------------------------------------
I, a, b    = -0.1, 0.7, 0.8
w_fold     = w_fold_left(I)
delta_U    = potential_barrier(W_FP, I=I, a=a, b=b)

print(f"\n{'='*62}")
print(f"  FHN Diagnostics — Kramers/Fold-layer crossover  [{tag}]")
print(f"  sigma: {sigma_vals[0]:.3f} – {sigma_vals[-1]:.3f}  ({len(sigma_vals)} pts)")
print(f"  eps:   {eps_vals[0]:.4f} – {eps_vals[-1]:.4f}  ({len(eps_vals)} pts)")
print(f"  N_traj={n_traj},  T={T},  seed={args.seed}")
print(f"  W_FP={W_FP:.4f},  w_fold={w_fold:.4f}")
print(f"  Delta_U = {delta_U:.5f}  (barrier height at W_FP)")
print(f"{'='*62}\n")

# Step 1 preview: print crossover sigma at each eps
print("Step 1 — Predicted crossover curve sigma_cross(eps):")
for eps in eps_vals[::max(1, len(eps_vals)//6)]:
    sc = np.sqrt(delta_U / np.log(1/eps)) if eps < 1 else np.nan
    print(f"  eps={eps:.4f}  sigma_cross={sc:.4f}")
print()

# ---------------------------------------------------------------------------
# Run sweep
# ---------------------------------------------------------------------------
print("Running sweep...")
res = sweep_grid(
    sigma_vals, eps_vals,
    n_trajectories=n_traj, T=T, seed=args.seed, verbose=True,
)

np.savez(f"data/diagnostics_{tag}.npz",
         **{k: v for k, v in res.items() if isinstance(v, np.ndarray)})

# ---------------------------------------------------------------------------
# Failure boundary + power law
# ---------------------------------------------------------------------------
eps_bd, sig_bd = failure_boundary(res, delta_ratio=args.delta_ratio)
print(f"\nBoundary points: {len(eps_bd)}")

alpha_global, C_global = fit_power_law(eps_bd, sig_bd)
print(f"Global fit:  sigma = {C_global:.3f} * eps^{alpha_global:.3f}")
print(f"Fold-layer:  sigma ~ eps^0.5")
print(f"BG:          sigma ~ eps^1.5")

# Compare boundary to crossover curve
if len(eps_bd) > 0:
    sigma_cross_at_bd = kramers_crossover_curve(eps_bd, delta_U)
    print(f"\nStep 1 — Boundary vs crossover curve:")
    print(f"  {'eps':>8}  {'sigma_bd':>10}  {'sigma_cross':>12}  {'ratio':>8}")
    for ep, sb, sc in zip(eps_bd, sig_bd, sigma_cross_at_bd):
        print(f"  {ep:8.4f}  {sb:10.4f}  {sc:12.4f}  {sb/sc:8.3f}")

# ---------------------------------------------------------------------------
# Step 1: Regime map with crossover curve overlay
# ---------------------------------------------------------------------------
fig = plot_regime_map(
    res,
    boundary_eps=eps_bd, boundary_sigma=sig_bd,
    bg_C=1.14,
    crossover_delta_U=delta_U,
    save_path=f"figures/regime_map_{tag}.png",
)
plt.close(fig)
print(f"\nStep 1 saved: figures/regime_map_{tag}.png")

# Step 1: Boundary scaling with crossover overlay
if len(eps_bd) >= 2:
    fig = plot_boundary_scaling(
        eps_bd, sig_bd, alpha_global, C_global,
        crossover_delta_U=delta_U,
        save_path=f"figures/boundary_scaling_{tag}.png",
    )
    plt.close(fig)
    print(f"Step 1 saved: figures/boundary_scaling_{tag}.png")

# ---------------------------------------------------------------------------
# Step 2: Local exponent (Diagnostic 1)
# ---------------------------------------------------------------------------
for win in [3, 4]:
    eps_c, alphas = local_exponent(eps_bd, sig_bd, window=win)
    if len(alphas) >= 2:
        fig = plot_local_exponent(
            eps_c, alphas, alpha_global,
            save_path=f"figures/local_exponent_w{win}_{tag}.png",
        )
        plt.close(fig)
        print(f"Step 2 saved: figures/local_exponent_w{win}_{tag}.png")
        print(f"  Local exponents (window={win}):")
        for ec, al in zip(eps_c, alphas):
            sc = np.sqrt(delta_U / np.log(1/ec)) if ec < 1 else np.nan
            region = "Kramers" if sig_bd[np.argmin(np.abs(eps_bd - ec))] > sc else "fold-layer"
            print(f"    eps={ec:.4f}  alpha={al:.3f}  [{region} regime]")

# ---------------------------------------------------------------------------
# Step 3: w_escape heatmap — empirical mechanism discriminator
# ---------------------------------------------------------------------------
fig = plot_w_spike_map(
    res,
    boundary_eps=eps_bd, boundary_sigma=sig_bd,
    use_escape=True,
    save_path=f"figures/w_escape_map_{tag}.png",
)
plt.close(fig)
print(f"\nStep 3 saved: figures/w_escape_map_{tag}.png")

# Step 3: print w_escape summary, labelling regime
print(f"\nStep 3 — w_escape summary:")
print(f"  W_FP={W_FP:.4f} (Kramers ref),  w_fold={w_fold:.4f} (fold-layer ref)")
print(f"  {'eps':>6}  {'sigma':>7}  {'w_esc':>8}  {'dist_FP':>9}  {'dist_fold':>10}  regime")
for j, eps in enumerate(eps_vals[::max(1, len(eps_vals)//8)]):
    j_full = np.argmin(np.abs(res["eps_vals"] - eps))
    sc     = np.sqrt(delta_U / np.log(1/eps)) if eps < 1 else np.nan
    # Find boundary sigma for this eps
    if len(eps_bd) > 0:
        bd_j  = np.argmin(np.abs(eps_bd - eps))
        i_sig = np.argmin(np.abs(sigma_vals - sig_bd[bd_j]))
    else:
        i_sig = len(sigma_vals) // 2
    we = res["w_escape_mean"][i_sig, j_full]
    if not np.isnan(we):
        regime = "Kramers" if (we - W_FP) > -0.01 else "fold-layer"
        print(f"  {eps:6.4f}  {sigma_vals[i_sig]:7.4f}  {we:8.4f}  "
              f"{we-W_FP:+9.4f}  {we-w_fold:+10.4f}  {regime}")

print(f"\nAll outputs saved to figures/ and data/")
