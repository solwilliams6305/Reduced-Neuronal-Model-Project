#!/usr/bin/env python3
"""
plot_exposure_scaling_boundary.py

Tests the "slow-time exposure" theory for the stochastic FHN failure boundary.

Motivation
----------
The extracted log_ratio=0 contour has fitted exponent around 0.48, close to 1/2.
That suggests

    sigma_fail^2 ~ epsilon,

which is the natural scaling if failure is controlled by accumulated fast-noise
exposure over an O(1/epsilon) slow-time window:

    accumulated exposure ~ sigma^2 * T_drift(eps) ~ sigma^2 / epsilon.

A fixed failure probability therefore gives

    sigma_fail ~ C sqrt(epsilon).

This script overlays several simple exposure-scaling laws on the existing regime
map and prints quantitative fits:

    Model A: pure square-root exposure law
        sigma = C sqrt(epsilon)

    Model B: affine variance law
        sigma^2 = alpha + beta epsilon

    Model C: free power law
        sigma = C epsilon^p

    Model D: log-corrected exposure law
        sigma^2 = beta epsilon / log(K/epsilon)
        or, if K is not fitted, sigma = C sqrt(epsilon)

The point is not to "prove" the theory, but to test whether the grey zero contour
is better described by sqrt(epsilon)-type exposure scaling than by frozen Kramers.

Usage
-----
    python plot_exposure_scaling_boundary.py --tau 1.1615

    python plot_exposure_scaling_boundary.py --fit-range 0.02 0.16

    python plot_exposure_scaling_boundary.py --rerun --n-traj 600
"""

from __future__ import annotations

import argparse
import os
import sys
import multiprocessing

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from scipy.optimize import curve_fit


try:
    multiprocessing.set_start_method("fork", force=True)
except RuntimeError:
    pass

sys.path.insert(0, os.path.dirname(__file__))

from simulate import (
    potential_barrier,
    w_fold_left,
    W_FP,
    V_FP,
    reduced_drift_time,
)

from sweep import (
    sweep_grid,
    failure_boundary,
    kramers_crossover_curve,
)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description="Test sqrt(epsilon) exposure-scaling theory against FHN regime map."
)

parser.add_argument("--tau", type=float, default=None,
                    help="Override tau_v for OU spread line.")
parser.add_argument("--rerun", action="store_true",
                    help="Rerun coarse sweep instead of loading data/sweep_result.npz.")
parser.add_argument("--n-traj", type=int, default=600,
                    help="Trajectories per grid point if --rerun is used.")
parser.add_argument("--seed", type=int, default=42)
parser.add_argument("--sweep-path", type=str, default="data/sweep_result.npz")
parser.add_argument("--out", type=str,
                    default="figures/regime_map_exposure_scaling.png")
parser.add_argument("--fit-range", type=float, nargs=2, default=None,
                    metavar=("EPS_MIN", "EPS_MAX"),
                    help="Fit only zero-contour points with EPS_MIN <= eps <= EPS_MAX.")
parser.add_argument("--bootstrap", type=int, default=2000,
                    help="Bootstrap samples for exponent uncertainty. Set 0 to disable.")
parser.add_argument("--show-empirical-x2", action="store_true",
                    help="Also plot the x2 MFPT-error empirical boundary in green.")
parser.add_argument("--no-frozen", action="store_true",
                    help="Hide frozen/basic Kramers comparison curves.")

args = parser.parse_args()


# ---------------------------------------------------------------------------
# Parameters and data
# ---------------------------------------------------------------------------

I, a, b = -0.1, 0.7, 0.8

os.makedirs("figures", exist_ok=True)
os.makedirs("data", exist_ok=True)

w_fold = w_fold_left(I)
delta_U0 = potential_barrier(W_FP, I=I, a=a, b=b)
fold_dist = abs(W_FP - w_fold)

if args.tau is not None:
    TAU_V = float(args.tau)
else:
    TAU_V = 1.0 / abs(1.0 - V_FP**2)

print("\nParameters")
print(f"  I={I}, a={a}, b={b}")
print(f"  V_FP={V_FP:.6f}, W_FP={W_FP:.6f}, w_fold={w_fold:.6f}")
print(f"  Delta U(W_FP)={delta_U0:.6f}")
print(f"  tau_v={TAU_V:.6f}")


if args.rerun or not os.path.exists(args.sweep_path):
    print("\nRunning coarse sweep...")
    sigma_vals_run = np.linspace(0.05, 0.45, 20)
    eps_vals_run = np.linspace(0.01, 0.25, 15)

    result = sweep_grid(
        sigma_vals_run,
        eps_vals_run,
        n_trajectories=args.n_traj,
        T=400.0,
        seed=args.seed,
        verbose=True,
    )

    np.savez(
        args.sweep_path,
        sigma_vals=result["sigma_vals"],
        eps_vals=result["eps_vals"],
        mfpt_full=result["mfpt_full"],
        mfpt_reduced=result["mfpt_reduced"],
        log_ratio=result["log_ratio"],
        rel_error=result["rel_error"],
        drift_times=result["drift_times"],
    )
else:
    print(f"\nLoading sweep data from {args.sweep_path}")
    raw = np.load(args.sweep_path)
    result = {k: raw[k] for k in raw.files}
    result["params"] = dict(I=I, a=a, b=b)

sigma_vals = np.asarray(result["sigma_vals"], dtype=float)
eps_vals = np.asarray(result["eps_vals"], dtype=float)
log_rat = np.asarray(result["log_ratio"], dtype=float)


# ---------------------------------------------------------------------------
# Zero-contour extraction
# ---------------------------------------------------------------------------

def extract_zero_contour(eps_vals: np.ndarray, sigma_vals: np.ndarray, log_rat: np.ndarray):
    eps_zero = []
    sigma_zero = []

    for j, eps in enumerate(eps_vals):
        col = log_rat[:, j]

        for i in range(len(sigma_vals) - 1):
            y0, y1 = col[i], col[i + 1]

            if not (np.isfinite(y0) and np.isfinite(y1)):
                continue

            if (y0 <= 0 < y1) or (y0 >= 0 > y1):
                if abs(y1 - y0) < 1e-14:
                    continue

                frac = -y0 / (y1 - y0)
                s = sigma_vals[i] + frac * (sigma_vals[i + 1] - sigma_vals[i])
                eps_zero.append(float(eps))
                sigma_zero.append(float(s))
                break

    return np.asarray(eps_zero), np.asarray(sigma_zero)


eps_zero_all, sigma_zero_all = extract_zero_contour(eps_vals, sigma_vals, log_rat)

if len(eps_zero_all) < 3:
    raise RuntimeError("Too few zero-contour points to fit exposure scaling.")

fit_mask = np.isfinite(eps_zero_all) & np.isfinite(sigma_zero_all) & (eps_zero_all > 0) & (sigma_zero_all > 0)

if args.fit_range is not None:
    lo, hi = args.fit_range
    fit_mask &= (eps_zero_all >= lo) & (eps_zero_all <= hi)

eps_fit = eps_zero_all[fit_mask]
sigma_fit = sigma_zero_all[fit_mask]

if len(eps_fit) < 3:
    raise RuntimeError("Too few zero-contour points inside requested fit range.")

print(f"\nUsing {len(eps_fit)} zero-contour points for fits.")
print(f"  eps range: {eps_fit.min():.4f} to {eps_fit.max():.4f}")


# ---------------------------------------------------------------------------
# Fit models
# ---------------------------------------------------------------------------

def rmse(y, yhat):
    return float(np.sqrt(np.mean((np.asarray(y) - np.asarray(yhat))**2)))


def aic_gaussian(y, yhat, k):
    y = np.asarray(y)
    yhat = np.asarray(yhat)
    n = len(y)
    rss = float(np.sum((y - yhat)**2))
    rss = max(rss, 1e-16)
    return float(n * np.log(rss / n) + 2 * k)


# Model A: sigma = C sqrt(eps), least squares in sigma.
C_sqrt = float(np.sum(sigma_fit * np.sqrt(eps_fit)) / np.sum(eps_fit))
pred_sqrt_fit = C_sqrt * np.sqrt(eps_fit)

# Model B: sigma^2 = alpha + beta eps.
X = np.column_stack([np.ones_like(eps_fit), eps_fit])
alpha_var, beta_var = np.linalg.lstsq(X, sigma_fit**2, rcond=None)[0]
alpha_var = float(alpha_var)
beta_var = float(beta_var)

def pred_affine_var(eps):
    return np.sqrt(np.maximum(alpha_var + beta_var * eps, 0.0))

# Model C: free power law sigma = C eps^p.
log_eps = np.log(eps_fit)
log_sig = np.log(sigma_fit)
p_power, logC_power = np.polyfit(log_eps, log_sig, 1)
C_power = float(np.exp(logC_power))
p_power = float(p_power)

def pred_power(eps):
    return C_power * eps**p_power

# Model D: log-corrected exposure law.
# sigma^2 = beta * eps / log(K/eps), K > max(eps).
# Let K be fitted nonlinearly, but keep it modestly constrained.
def logcorr_model(eps, beta, K):
    denom = np.log(np.maximum(K / eps, 1.000001))
    return np.sqrt(np.maximum(beta * eps / denom, 0.0))

try:
    popt, pcov = curve_fit(
        logcorr_model,
        eps_fit,
        sigma_fit,
        p0=[C_sqrt**2 * np.log(1.0 / np.mean(eps_fit)), 1.0],
        bounds=([1e-8, float(np.max(eps_fit) * 1.01)], [100.0, 100.0]),
        maxfev=20000,
    )
    beta_logcorr, K_logcorr = map(float, popt)
except Exception as exc:
    print(f"  Warning: log-corrected fit failed: {exc}")
    beta_logcorr, K_logcorr = np.nan, np.nan

def pred_logcorr(eps):
    if not np.isfinite(beta_logcorr):
        return np.full_like(np.asarray(eps), np.nan, dtype=float)
    return logcorr_model(np.asarray(eps), beta_logcorr, K_logcorr)


# Bootstrap exponent uncertainty for free power law.
p_boot = []
if args.bootstrap and args.bootstrap > 0:
    rng = np.random.default_rng(args.seed + 12345)
    n = len(eps_fit)
    for _ in range(args.bootstrap):
        idx = rng.integers(0, n, size=n)
        if len(np.unique(idx)) < 2:
            continue
        try:
            pp, _ = np.polyfit(np.log(eps_fit[idx]), np.log(sigma_fit[idx]), 1)
            if np.isfinite(pp):
                p_boot.append(float(pp))
        except Exception:
            pass

if p_boot:
    p_ci = np.percentile(p_boot, [2.5, 50, 97.5])
else:
    p_ci = [np.nan, np.nan, np.nan]


# ---------------------------------------------------------------------------
# Print diagnostics
# ---------------------------------------------------------------------------

models = []

models.append((
    "sqrt exposure",
    pred_sqrt_fit,
    1,
    f"sigma = {C_sqrt:.4f} sqrt(eps)"
))

pred_aff_fit = pred_affine_var(eps_fit)
models.append((
    "affine variance",
    pred_aff_fit,
    2,
    f"sigma^2 = {alpha_var:.5f} + {beta_var:.5f} eps"
))

pred_pow_fit = pred_power(eps_fit)
models.append((
    "free power",
    pred_pow_fit,
    2,
    f"sigma = {C_power:.4f} eps^{p_power:.4f}"
))

pred_log_fit = pred_logcorr(eps_fit)
models.append((
    "log-corrected exposure",
    pred_log_fit,
    2,
    f"sigma^2 = {beta_logcorr:.5f} eps / log({K_logcorr:.4f}/eps)"
))

print("\nFit summary on extracted zero contour")
print(f"{'model':>24}  {'RMSE':>10}  {'AIC':>10}  formula")
print("-" * 92)

for name, pred, k, formula in models:
    print(f"{name:>24}  {rmse(sigma_fit, pred):10.5f}  {aic_gaussian(sigma_fit, pred, k):10.3f}  {formula}")

print("\nExponent diagnostic")
print(f"  free power exponent p = {p_power:.4f}")
if np.isfinite(p_ci[0]):
    print(f"  bootstrap 95% CI     = [{p_ci[0]:.4f}, {p_ci[2]:.4f}]")
print(f"  distance from 1/2    = {p_power - 0.5:+.4f}")

print("\nExposure interpretation")
print("  If failure is controlled by accumulated exposure sigma^2*T_drift")
print("  and T_drift ~ const/eps, then fixed failure probability gives")
print("      sigma^2/eps ~ const, i.e. sigma ~ sqrt(eps).")


# ---------------------------------------------------------------------------
# Other comparison curves
# ---------------------------------------------------------------------------

eps_plot = np.linspace(float(eps_vals.min()), float(eps_vals.max()), 600)

sigma_sqrt = C_sqrt * np.sqrt(eps_plot)
sigma_affine_var = pred_affine_var(eps_plot)
sigma_power = pred_power(eps_plot)
sigma_logcorr = pred_logcorr(eps_plot)

sigma_spread = fold_dist * np.sqrt(2.0 * b / (eps_plot * TAU_V))
sigma_cross = kramers_crossover_curve(eps_plot, delta_U0)

# Frozen exact Kramers match.
A_WFP = 1.0
try:
    # local curvature prefactor, copied lightly to avoid extra dependencies
    from scipy.optimize import brentq
    def _pref_local(w):
        def h(v): return v - v**3/3 - w + I
        vm = brentq(h, -3.0, -1.001)
        vs = brentq(h, -0.999, 0.999)
        return (1.0/(2.0*np.pi))*np.sqrt(abs(vm**2 - 1.0)*abs(vs**2 - 1.0))
    A_WFP = float(_pref_local(W_FP))
except Exception:
    A_WFP = 1.0

T_plot = np.array([reduced_drift_time(eps=float(e), I=I, a=a, b=b) for e in eps_plot])
log_arg = T_plot * A_WFP
sigma_star = np.full_like(eps_plot, np.nan)
valid = log_arg > 1.0
sigma_star[valid] = np.sqrt(2.0 * delta_U0 / np.log(log_arg[valid]))

eps_bd, sigma_bd = failure_boundary(result, delta_ratio=2.0)


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(8.8, 6.5))

abs_vals = np.abs(log_rat[np.isfinite(log_rat)])
vmax = float(np.nanpercentile(abs_vals, 95)) if len(abs_vals) else 1.0
vmax = max(vmax, 1.0)
norm = TwoSlopeNorm(vmin=-vmax/2.0, vcenter=0.0, vmax=vmax)

im = ax.pcolormesh(
    eps_vals,
    sigma_vals,
    log_rat,
    cmap="RdBu_r",
    norm=norm,
    shading="auto",
)

cbar = fig.colorbar(im, ax=ax, pad=0.02)
cbar.set_label(
    r"$\log(T_{\rm drift}/{\rm MFPT}_{\rm full})$"
    "\n(+ve = reduction fails)",
    fontsize=9,
)

# Grey zero contour from heatmap.
try:
    ax.contour(
        eps_vals,
        sigma_vals,
        log_rat,
        levels=[0],
        colors=["0.35"],
        linewidths=2.0,
        linestyles="--",
        alpha=0.85,
    )
except Exception as exc:
    print(f"Warning: could not draw zero contour: {exc}")

# Zero contour interpolation points.
ax.plot(
    eps_zero_all,
    sigma_zero_all,
    "x",
    color="0.15",
    ms=5,
    mew=1.2,
    alpha=0.8,
    label="Extracted zero-contour points",
)

# Optional x2 boundary.
if args.show_empirical_x2 and len(eps_bd) > 0:
    ax.plot(
        eps_bd,
        sigma_bd,
        "o-",
        color="#1a9641",
        ms=5,
        lw=1.8,
        zorder=6,
        label=r"Empirical boundary ($\times 2$ MFPT error)",
    )

# Main exposure theory curves.
mask = np.isfinite(sigma_sqrt) & (sigma_sqrt >= sigma_vals.min()) & (sigma_sqrt <= sigma_vals.max())
ax.plot(
    eps_plot[mask],
    sigma_sqrt[mask],
    color="#984ea3",
    lw=3.2,
    ls="-",
    zorder=9,
    label=(
        r"Exposure law: $\sigma=C\sqrt{\varepsilon}$"
        f"\n  $C={C_sqrt:.3f}$"
    ),
)

mask = np.isfinite(sigma_affine_var) & (sigma_affine_var >= sigma_vals.min()) & (sigma_affine_var <= sigma_vals.max())
ax.plot(
    eps_plot[mask],
    sigma_affine_var[mask],
    color="#984ea3",
    lw=2.1,
    ls="--",
    zorder=8,
    label=(
        r"Affine variance: $\sigma^2=\alpha+\beta\varepsilon$"
        f"\n  $\\alpha={alpha_var:.3g}$, $\\beta={beta_var:.3g}$"
    ),
)

mask = np.isfinite(sigma_power) & (sigma_power >= sigma_vals.min()) & (sigma_power <= sigma_vals.max())
ax.plot(
    eps_plot[mask],
    sigma_power[mask],
    color="#000000",
    lw=1.7,
    ls="-",
    alpha=0.7,
    zorder=7,
    label=(
        r"Free power fit"
        f"\n  $p={p_power:.3f}$"
    ),
)

if np.isfinite(beta_logcorr):
    mask = np.isfinite(sigma_logcorr) & (sigma_logcorr >= sigma_vals.min()) & (sigma_logcorr <= sigma_vals.max())
    ax.plot(
        eps_plot[mask],
        sigma_logcorr[mask],
        color="#7b3294",
        lw=1.8,
        ls=":",
        zorder=7,
        label=r"Log-corrected exposure",
    )

# Older comparison curves.
if not args.no_frozen:
    mask = np.isfinite(sigma_cross) & (sigma_cross >= sigma_vals.min()) & (sigma_cross <= sigma_vals.max())
    ax.plot(
        eps_plot[mask],
        sigma_cross[mask],
        color="#fdae61",
        lw=2.0,
        ls="-.",
        label=r"Basic Kramers",
    )

    mask = np.isfinite(sigma_star) & (sigma_star >= sigma_vals.min()) & (sigma_star <= sigma_vals.max())
    ax.plot(
        eps_plot[mask],
        sigma_star[mask],
        color="#d6604d",
        lw=1.8,
        ls=":",
        label=r"Frozen Kramers match",
    )

ax.set_xlabel(r"$\varepsilon$  (timescale separation)", fontsize=12)
ax.set_ylabel(r"$\sigma$  (noise amplitude)", fontsize=12)

title = (
    "Regime map: testing slow-time exposure scaling\n"
    rf"Extracted contour exponent $p={p_power:.3f}$; square-root theory predicts $p=1/2$"
)
ax.set_title(title, fontsize=12)

ax.set_xlim(float(eps_vals.min()), float(eps_vals.max()))
ax.set_ylim(float(sigma_vals.min()), float(sigma_vals.max()))
ax.grid(alpha=0.18, lw=0.6)

leg = ax.legend(fontsize=7.2, loc="upper left", frameon=True, framealpha=0.92)
leg.get_frame().set_linewidth(0.5)

fig.tight_layout()

out_dir = os.path.dirname(args.out)
if out_dir:
    os.makedirs(out_dir, exist_ok=True)
fig.savefig(args.out, dpi=180, bbox_inches="tight")
plt.close(fig)

print(f"\nSaved: {args.out}")


# ---------------------------------------------------------------------------
# Secondary diagnostic plot: sigma^2 vs epsilon
# ---------------------------------------------------------------------------

fig2, ax2 = plt.subplots(figsize=(6.4, 4.4))

ax2.plot(eps_zero_all, sigma_zero_all**2, "x", color="0.15", label="zero contour data")
ax2.plot(eps_plot, sigma_sqrt**2, lw=2.5, label=r"$C^2\varepsilon$")
ax2.plot(eps_plot, sigma_affine_var**2, lw=2.0, ls="--",
         label=r"$\alpha+\beta\varepsilon$")

ax2.set_xlabel(r"$\varepsilon$")
ax2.set_ylabel(r"$\sigma_{\rm fail}^2$")
ax2.set_title(r"Variance scaling test: does $\sigma_{\rm fail}^2$ grow linearly with $\varepsilon$?")
ax2.grid(alpha=0.25)
ax2.legend(fontsize=8)

fig2.tight_layout()
secondary = os.path.join(os.path.dirname(args.out) or ".", "exposure_variance_scaling.png")
fig2.savefig(secondary, dpi=180, bbox_inches="tight")
plt.close(fig2)

print(f"Saved: {secondary}")
