#!/usr/bin/env python3
"""
plot_predicted_boundaries_random_barrier.py

Regime-map plotting script for the stochastic FHN reduction project.

This version implements the more mechanistic correction:

    fast noise in v
        -> induced instantaneous movement in w
        -> stochastic barrier process Delta U(W_t)
        -> random Kramers intensity
             kappa_t = A(W_t) exp(-2 Delta U(W_t)/sigma^2)
        -> survival functional
             S = E[ exp(- int_0^T kappa_t dt) ]
        -> effective hazard
             H_eff = -log S

The corrected boundary is

    H_eff(sigma, eps; T_drift(eps)) = H_c.

Important modelling choice
--------------------------
We do NOT add direct noise to dw. Instead we simulate the linearised
subthreshold dynamics around the fixed point:

    d v_tilde = ((1 - V_FP^2) v_tilde - w_tilde) dt + sigma dB_t
    d w_tilde = eps (v_tilde - b w_tilde) dt.

Thus w is random because it integrates noisy v. This is the mechanism
you described: fast noise zigs v around, w instantaneously responds, and
the potential barrier Delta U(W_t) is shifted at the same time escape is
being attempted.

Usage
-----
    python plot_predicted_boundaries_random_barrier.py --tau 1.1615

    python plot_predicted_boundaries_random_barrier.py --tau 1.1615 --fit-hazard-threshold

    python plot_predicted_boundaries_random_barrier.py --tau 1.1615 --n-paths 512 --dt-path 0.05

If you want to overwrite your old plotting script:

    cp plot_predicted_boundaries_random_barrier.py plot_predicted_boundaries.py
"""

from __future__ import annotations

import argparse
import os
import sys
import multiprocessing
from functools import lru_cache
from typing import Optional

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from scipy.optimize import brentq
from scipy.interpolate import interp1d


# ---------------------------------------------------------------------------
# Multiprocessing / imports
# ---------------------------------------------------------------------------

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
    description="Plot stochastic FHN regime map with random-barrier survival boundary."
)

parser.add_argument("--tau", type=float, default=None,
                    help="Override tau_v for OU spread line only.")
parser.add_argument("--rerun", action="store_true",
                    help="Rerun coarse sweep instead of loading data/sweep_result.npz.")
parser.add_argument("--n-traj", type=int, default=600,
                    help="Trajectories per grid point if --rerun is used.")
parser.add_argument("--seed", type=int, default=42,
                    help="Random seed.")
parser.add_argument("--sweep-path", type=str, default="data/sweep_result.npz",
                    help="Path to saved sweep_result.npz.")
parser.add_argument("--out", type=str,
                    default="figures/regime_map_predicted_boundaries_random_barrier.png",
                    help="Output figure path.")

# Random-barrier solver controls.
parser.add_argument("--hazard-threshold", type=float, default=1.0,
                    help="Boundary threshold H_c for H_eff = H_c.")
parser.add_argument("--fit-hazard-threshold", action="store_true",
                    help="Fit H_c from extracted zero contour using median H_eff.")
parser.add_argument("--n-paths", type=int, default=384,
                    help="Number of linearised subthreshold paths per (sigma, eps).")
parser.add_argument("--dt-path", type=float, default=0.05,
                    help="Time step for random barrier paths.")
parser.add_argument("--eps-curve-points", type=int, default=90,
                    help="Number of eps points used for purple random-barrier curve.")
parser.add_argument("--sigma-scan-points", type=int, default=50,
                    help="Number of sigma scan points used to bracket root.")
parser.add_argument("--max-path-steps", type=int, default=5000,
                    help="Cap on path steps. If T/dt is larger, dt is enlarged.")
parser.add_argument("--clip-w-low", type=float, default=0.015,
                    help="How far below fold to allow W before treating barrier as zero.")
parser.add_argument("--clip-w-high", type=float, default=0.12,
                    help="How far above W_FP to precompute barriers.")
parser.add_argument("--reuse-noise", action="store_true",
                    help="Use same standard Brownian increments across sigma for each eps.")
parser.add_argument("--no-random-barrier", action="store_true",
                    help="Skip the random-barrier purple line.")

args = parser.parse_args()


# ---------------------------------------------------------------------------
# Parameters
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
print(f"  V_FP={V_FP:.4f}, W_FP={W_FP:.4f}, w_fold={w_fold:.4f}")
print(f"  |W_FP - w_fold|={fold_dist:.4f}")
print(f"  Delta U(W_FP)={delta_U0:.5f}")
print(f"  tau_v={TAU_V:.4f}")


# ---------------------------------------------------------------------------
# Load or rerun sweep
# ---------------------------------------------------------------------------

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
    print(f"Saved sweep to {args.sweep_path}")
else:
    print(f"\nLoading sweep data from {args.sweep_path}")
    raw = np.load(args.sweep_path)
    result = {k: raw[k] for k in raw.files}
    result["params"] = dict(I=I, a=a, b=b)

sigma_vals = np.asarray(result["sigma_vals"], dtype=float)
eps_vals = np.asarray(result["eps_vals"], dtype=float)
log_rat = np.asarray(result["log_ratio"], dtype=float)


# ---------------------------------------------------------------------------
# Frozen-w potential machinery
# ---------------------------------------------------------------------------

def v_saddle(w: float) -> Optional[float]:
    """Middle unstable root of v - v^3/3 - w + I = 0."""
    def f(v: float) -> float:
        return v - v**3 / 3.0 - w + I

    try:
        return brentq(f, -0.999, 0.999)
    except ValueError:
        return None


def v_min_left(w: float) -> Optional[float]:
    """Left stable root of v - v^3/3 - w + I = 0."""
    def f(v: float) -> float:
        return v - v**3 / 3.0 - w + I

    try:
        return brentq(f, -3.0, -1.001)
    except ValueError:
        return None


def U(v: float, w: float) -> float:
    """Fast potential U(v,w) with dV/dt = -dU/dv plus noise."""
    return -v**2 / 2.0 + v**4 / 12.0 + (w - I) * v


def delta_U(w: float) -> float:
    """Barrier height Delta U(w) = U(v_sad,w) - U(v_min,w)."""
    vm = v_min_left(w)
    vs = v_saddle(w)

    if vm is None or vs is None:
        # Past a fold the Kramers barrier has collapsed.
        return 0.0

    return max(U(vs, w) - U(vm, w), 0.0)


def kramers_prefactor(w: float) -> Optional[float]:
    """
    Frozen-w Kramers prefactor

        A(w) = (1/(2pi)) sqrt(|U''(v_min)| |U''(v_sad)|),

    where U''(v)=v^2-1.
    """
    vm = v_min_left(w)
    vs = v_saddle(w)

    if vm is None or vs is None:
        return None

    curv_min = abs(vm**2 - 1.0)
    curv_sad = abs(vs**2 - 1.0)

    return (1.0 / (2.0 * np.pi)) * np.sqrt(curv_min * curv_sad)


A_WFP = kramers_prefactor(W_FP)
A_WFP = float(A_WFP) if A_WFP is not None and np.isfinite(A_WFP) else 1.0
print(f"  A(W_FP)={A_WFP:.5f}")


# ---------------------------------------------------------------------------
# Precompute Delta U(w), A(w) interpolation for speed
# ---------------------------------------------------------------------------

w_grid_min = w_fold - args.clip_w_low
w_grid_max = W_FP + args.clip_w_high
w_pre = np.linspace(w_grid_min, w_grid_max, 5000)

dU_pre = np.array([delta_U(float(w)) for w in w_pre], dtype=float)
A_pre_raw = np.array(
    [
        np.nan if kramers_prefactor(float(w)) is None else float(kramers_prefactor(float(w)))
        for w in w_pre
    ],
    dtype=float,
)

# Fill invalid A values. When the barrier is gone, Delta U=0. Use a conservative
# finite prefactor rather than exploding or producing NaNs.
valid_A = np.isfinite(A_pre_raw) & (A_pre_raw > 0)
if np.any(valid_A):
    A_fill = A_pre_raw.copy()

    # Forward fill.
    for i in range(1, len(A_fill)):
        if not (np.isfinite(A_fill[i]) and A_fill[i] > 0):
            A_fill[i] = A_fill[i - 1]

    # Backward fill.
    for i in range(len(A_fill) - 2, -1, -1):
        if not (np.isfinite(A_fill[i]) and A_fill[i] > 0):
            A_fill[i] = A_fill[i + 1]

    A_pre = A_fill
else:
    A_pre = np.full_like(w_pre, A_WFP)

# Interpolators clamp outside the precomputed range.
dU_interp = interp1d(
    w_pre, dU_pre,
    kind="linear",
    bounds_error=False,
    fill_value=(float(dU_pre[0]), float(dU_pre[-1])),
    assume_sorted=True,
)

A_interp = interp1d(
    w_pre, A_pre,
    kind="linear",
    bounds_error=False,
    fill_value=(float(A_pre[0]), float(A_pre[-1])),
    assume_sorted=True,
)


# ---------------------------------------------------------------------------
# Empirical / diagnostic boundaries
# ---------------------------------------------------------------------------

eps_bd, sigma_bd = failure_boundary(result, delta_ratio=2.0)
eps_bd = np.asarray(eps_bd, dtype=float)
sigma_bd = np.asarray(sigma_bd, dtype=float)

eps_zero_list: list[float] = []
sigma_zero_list: list[float] = []

for j, eps in enumerate(eps_vals):
    col = log_rat[:, j]

    for i in range(len(sigma_vals) - 1):
        y0 = col[i]
        y1 = col[i + 1]

        if not (np.isfinite(y0) and np.isfinite(y1)):
            continue

        if (y0 <= 0 < y1) or (y0 >= 0 > y1):
            if abs(y1 - y0) < 1e-14:
                continue
            frac = -y0 / (y1 - y0)
            s = sigma_vals[i] + frac * (sigma_vals[i + 1] - sigma_vals[i])
            eps_zero_list.append(float(eps))
            sigma_zero_list.append(float(s))
            break

eps_zero = np.asarray(eps_zero_list, dtype=float)
sigma_zero = np.asarray(sigma_zero_list, dtype=float)

if len(eps_zero) >= 3:
    alpha_zero, log_C_zero = np.polyfit(np.log(eps_zero), np.log(sigma_zero), 1)
    C_zero = float(np.exp(log_C_zero))
    print(f"\nExtracted zero contour fit: sigma = {C_zero:.3f} eps^{alpha_zero:.3f}")
else:
    alpha_zero, C_zero = np.nan, np.nan


# ---------------------------------------------------------------------------
# Standard analytical curves
# ---------------------------------------------------------------------------

eps_plot = np.linspace(float(eps_vals.min()), float(eps_vals.max()), 400)

sigma_spread = fold_dist * np.sqrt(2.0 * b / (eps_plot * TAU_V))
sigma_cross = kramers_crossover_curve(eps_plot, delta_U0)

T_drift_plot = np.array(
    [reduced_drift_time(eps=float(e), I=I, a=a, b=b) for e in eps_plot],
    dtype=float,
)

log_arg = T_drift_plot * A_WFP
valid_star = log_arg > 1.0

sigma_star = np.full_like(eps_plot, np.nan, dtype=float)
sigma_star[valid_star] = np.sqrt(2.0 * delta_U0 / np.log(log_arg[valid_star]))


# ---------------------------------------------------------------------------
# Random-barrier survival model
# ---------------------------------------------------------------------------

def make_noise(seed: int, n_paths: int, n_steps: int) -> np.ndarray:
    """Standard normal increments, shape (n_steps, n_paths)."""
    rng = np.random.default_rng(seed)
    return rng.standard_normal((n_steps, n_paths)).astype(np.float64)


def simulate_linearised_W_paths(
    sigma: float,
    eps: float,
    T: float,
    n_paths: int,
    dt_requested: float,
    seed: int,
    max_steps: int,
    base_noise: Optional[np.ndarray] = None,
) -> tuple[np.ndarray, float]:
    """
    Simulate W_t paths from the linearised subthreshold dynamics around FP:

        d vtilde = ((1 - V_FP^2) vtilde - wtilde) dt + sigma dB_t
        d wtilde = eps (vtilde - b wtilde) dt.

    Returns
    -------
    W_paths : array, shape (n_steps+1, n_paths)
        Actual W = W_FP + wtilde.
    dt : float
        Possibly enlarged dt so that n_steps <= max_steps.
    """
    if T <= 0:
        return np.full((1, n_paths), W_FP, dtype=float), T

    n_steps = int(np.ceil(T / dt_requested))
    if n_steps > max_steps:
        n_steps = int(max_steps)

    n_steps = max(n_steps, 2)
    dt = float(T / n_steps)
    sqrt_dt = np.sqrt(dt)

    if base_noise is None or base_noise.shape != (n_steps, n_paths):
        Z = make_noise(seed, n_paths, n_steps)
    else:
        Z = base_noise

    vtilde = np.zeros(n_paths, dtype=np.float64)
    wtilde = np.zeros(n_paths, dtype=np.float64)

    W_paths = np.empty((n_steps + 1, n_paths), dtype=np.float64)
    W_paths[0, :] = W_FP

    alpha = 1.0 - V_FP**2

    # Euler-Maruyama. Use old v,w for both updates.
    for k in range(n_steps):
        dv = (alpha * vtilde - wtilde) * dt + sigma * sqrt_dt * Z[k]
        dw = eps * (vtilde - b * wtilde) * dt

        vtilde = vtilde + dv
        wtilde = wtilde + dw

        W_paths[k + 1, :] = W_FP + wtilde

    return W_paths, dt


@lru_cache(maxsize=4096)
def _random_barrier_hazard_cached(
    sigma_round: float,
    eps_round: float,
    H_seed: int,
) -> float:
    """
    Cached effective hazard. Rounded arguments are used to avoid pathological
    cache misses from brentq floating point calls.
    """
    sigma = float(sigma_round)
    eps = float(eps_round)

    T = float(reduced_drift_time(eps=eps, I=I, a=a, b=b))

    n_steps = int(np.ceil(T / args.dt_path))
    if n_steps > args.max_path_steps:
        n_steps = int(args.max_path_steps)
    n_steps = max(n_steps, 2)
    dt = float(T / n_steps)

    base_noise = None
    if args.reuse_noise:
        # Same underlying standard Brownian path for each eps. sigma still scales it.
        base_noise = make_noise(H_seed, args.n_paths, n_steps)

    W_paths, dt_used = simulate_linearised_W_paths(
        sigma=sigma,
        eps=eps,
        T=T,
        n_paths=args.n_paths,
        dt_requested=dt,
        seed=H_seed,
        max_steps=args.max_path_steps,
        base_noise=base_noise,
    )

    # Evaluate random instantaneous barrier and prefactor along the paths.
    W_eval = np.clip(W_paths, w_grid_min, w_grid_max)

    dU_vals = dU_interp(W_eval)
    A_vals = A_interp(W_eval)

    # Kramers intensity kappa_t = A(W_t) exp(-2 Delta U(W_t)/sigma^2).
    expo = np.exp(np.clip(-2.0 * dU_vals / max(sigma**2, 1e-14), -700, 0))
    kappa = A_vals * expo
    kappa = np.where(np.isfinite(kappa), kappa, 0.0)

    # Pathwise integrated hazard H_i = int kappa_t dt.
    # Trapezoid over time axis.
    H_paths = dt_used * (
        0.5 * kappa[0, :]
        + np.sum(kappa[1:-1, :], axis=0)
        + 0.5 * kappa[-1, :]
    )

    H_paths = np.where(np.isfinite(H_paths), H_paths, np.inf)

    # Effective survival S = E exp(-H_i), H_eff = -log S.
    survival = float(np.mean(np.exp(-np.clip(H_paths, 0, 745))))

    if survival <= 0 or not np.isfinite(survival):
        return float(np.nanmax(np.clip(H_paths, 0, 700)))

    H_eff = -np.log(survival)

    return float(H_eff)


def random_barrier_hazard(sigma: float, eps: float) -> float:
    """
    Wrapper around cached hazard.

    Rounding keeps brent root finding from triggering a completely new Monte
    Carlo simulation for numerically indistinguishable sigma values.
    """
    sigma_round = round(float(sigma), 6)
    eps_round = round(float(eps), 6)

    # A deterministic seed per eps keeps the curve smooth-ish while still
    # changing samples across eps.
    H_seed = int(args.seed + 100000 * eps_round) % (2**32 - 1)

    return _random_barrier_hazard_cached(sigma_round, eps_round, H_seed)


def solve_random_barrier_boundary(
    eps: float,
    H_c: float,
    sigma_min: Optional[float] = None,
    sigma_max: Optional[float] = None,
) -> float:
    """Solve H_eff(sigma, eps) = H_c."""
    if sigma_min is None:
        sigma_min = max(0.005, float(sigma_vals.min()) * 0.2)
    if sigma_max is None:
        sigma_max = max(float(sigma_vals.max()) * 1.8, 0.8)

    def obj(sig: float) -> float:
        H = random_barrier_hazard(float(sig), float(eps))
        if not np.isfinite(H):
            return np.nan
        return H - H_c

    scan = np.geomspace(sigma_min, sigma_max, args.sigma_scan_points)
    vals = np.array([obj(float(s)) for s in scan], dtype=float)

    # If already above threshold at smallest sigma, return smallest.
    if np.isfinite(vals[0]) and vals[0] >= 0:
        return float(scan[0])

    for i in range(len(scan) - 1):
        f0 = vals[i]
        f1 = vals[i + 1]

        if not (np.isfinite(f0) and np.isfinite(f1)):
            continue

        if f0 == 0:
            return float(scan[i])

        if f0 * f1 < 0:
            try:
                return float(brentq(obj, float(scan[i]), float(scan[i + 1]), maxiter=60))
            except ValueError:
                continue

    return np.nan


def fit_random_barrier_threshold() -> float:
    """Fit H_c from the extracted zero-contour points."""
    if len(eps_zero) < 2:
        return float(args.hazard_threshold)

    H_vals = []
    for e, s in zip(eps_zero, sigma_zero):
        H = random_barrier_hazard(float(s), float(e))
        if np.isfinite(H) and H > 0:
            H_vals.append(float(H))

    if not H_vals:
        return float(args.hazard_threshold)

    return float(np.median(H_vals))


if args.no_random_barrier:
    H_c = float(args.hazard_threshold)
    eps_rb = np.array([])
    sigma_rb = np.array([])
else:
    if args.fit_hazard_threshold:
        print("\nFitting H_c from extracted zero contour...")
        H_c = fit_random_barrier_threshold()
        print(f"  Fitted H_c = {H_c:.4g}")
    else:
        H_c = float(args.hazard_threshold)
        print(f"\nUsing H_c = {H_c:.4g}")

    eps_rb = np.linspace(float(eps_vals.min()), float(eps_vals.max()), args.eps_curve_points)

    print("\nSolving random-barrier boundary...")
    sigma_rb_list = []
    for idx, e in enumerate(eps_rb):
        s = solve_random_barrier_boundary(float(e), H_c=H_c)
        sigma_rb_list.append(s)
        if idx % max(1, args.eps_curve_points // 10) == 0:
            print(f"  eps={e:.4f}, sigma_rb={s:.4f}")

    sigma_rb = np.asarray(sigma_rb_list, dtype=float)


# ---------------------------------------------------------------------------
# Diagnostics table
# ---------------------------------------------------------------------------

print("\nBoundary comparison")
print(f"{'eps':>8}  {'OU spread':>10}  {'basic K':>10}  {'frozen K':>10}  {'zero':>10}  {'RB':>10}")
print("-" * 72)

for eps in eps_vals[::max(1, len(eps_vals) // 8)]:
    sp = fold_dist * np.sqrt(2.0 * b / (eps * TAU_V))
    sc = np.sqrt(delta_U0 / np.log(1.0 / eps)) if eps < 1 else np.nan

    td = reduced_drift_time(eps=float(eps), I=I, a=a, b=b)
    ss = np.sqrt(2.0 * delta_U0 / np.log(td * A_WFP)) if td * A_WFP > 1 else np.nan

    if len(eps_zero) > 0:
        k = np.argmin(np.abs(eps_zero - eps))
        sz = sigma_zero[k] if abs(eps_zero[k] - eps) < 0.02 else np.nan
    else:
        sz = np.nan

    if len(eps_rb) > 0:
        krb = np.argmin(np.abs(eps_rb - eps))
        srb = sigma_rb[krb] if abs(eps_rb[krb] - eps) < 0.02 else np.nan
    else:
        srb = np.nan

    print(f"{eps:8.4f}  {sp:10.4f}  {sc:10.4f}  {ss:10.4f}  {sz:10.4f}  {srb:10.4f}")


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(8.6, 6.4))

abs_vals = np.abs(log_rat[np.isfinite(log_rat)])
vmax = float(np.nanpercentile(abs_vals, 95)) if len(abs_vals) else 1.0
vmax = max(vmax, 1.0)

norm = TwoSlopeNorm(vmin=-vmax / 2.0, vcenter=0.0, vmax=vmax)

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

# Grey dashed zero contour.
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

# Empirical x2 MFPT-error boundary.
if len(eps_bd) > 0:
    ax.plot(
        eps_bd,
        sigma_bd,
        "o-",
        color="#1a9641",
        ms=6,
        lw=2.0,
        zorder=6,
        label=r"Empirical boundary ($\times 2$ MFPT error)",
    )

# OU spread.
mask = np.isfinite(sigma_spread) & (sigma_spread >= sigma_vals.min()) & (sigma_spread <= sigma_vals.max())
ax.plot(
    eps_plot[mask],
    sigma_spread[mask],
    color="#d6604d",
    lw=2.2,
    ls="-",
    label=(
        r"OU spread: $\sigma_w=|W_{FP}-w_{\rm fold}|$"
        f"\n  $\\sigma=\\Delta w\\sqrt{{2b/(\\varepsilon\\tau_v)}}$, $\\tau_v={TAU_V:.3f}$"
    ),
)

# Basic Kramers.
mask = np.isfinite(sigma_cross) & (sigma_cross >= sigma_vals.min()) & (sigma_cross <= sigma_vals.max())
ax.plot(
    eps_plot[mask],
    sigma_cross[mask],
    color="#fdae61",
    lw=2.2,
    ls="-.",
    label=(
        r"Basic Kramers: $\sigma=\sqrt{\Delta U/\log(1/\varepsilon)}$"
        f"\n  $\\Delta U(W_{{FP}})={delta_U0:.4f}$"
    ),
)

# Frozen exact Kramers.
mask = np.isfinite(sigma_star) & (sigma_star >= sigma_vals.min()) & (sigma_star <= sigma_vals.max())
ax.plot(
    eps_plot[mask],
    sigma_star[mask],
    color="#984ea3",
    lw=2.0,
    ls=":",
    label=(
        r"Frozen Kramers match: "
        r"$\sigma=\sqrt{2\Delta U/\log(T_{\rm drift}A)}$"
    ),
)

# Random-barrier survival boundary.
if len(eps_rb) > 0:
    mask = np.isfinite(sigma_rb) & (sigma_rb >= sigma_vals.min()) & (sigma_rb <= sigma_vals.max())
    ax.plot(
        eps_rb[mask],
        sigma_rb[mask],
        color="#984ea3",
        lw=3.2,
        ls="-",
        zorder=8,
        label=(
            r"Random-barrier survival: "
            r"$-\log\,\mathbb{E}e^{-\int\kappa_tdt}=H_c$"
            f"\n  $H_c={H_c:.3g}$, paths={args.n_paths}, dt≈{args.dt_path:g}"
        ),
    )

# Extracted zero contour points.
if len(eps_zero) > 0:
    ax.plot(
        eps_zero,
        sigma_zero,
        "x",
        color="0.15",
        ms=5,
        mew=1.2,
        alpha=0.75,
        label="Extracted zero-contour points",
    )

ax.set_xlabel(r"$\varepsilon$  (timescale separation)", fontsize=12)
ax.set_ylabel(r"$\sigma$  (noise amplitude)", fontsize=12)
ax.set_title(
    "Regime map with predicted failure boundaries\n"
    "Purple solid: random-barrier survival from fast-noise-induced barrier fluctuations",
    fontsize=12,
)

ax.set_xlim(float(eps_vals.min()), float(eps_vals.max()))
ax.set_ylim(float(sigma_vals.min()), float(sigma_vals.max()))
ax.grid(alpha=0.18, lw=0.6)

leg = ax.legend(
    fontsize=7.4,
    loc="upper left",
    frameon=True,
    framealpha=0.92,
)
leg.get_frame().set_linewidth(0.5)

fig.tight_layout()

out_dir = os.path.dirname(args.out)
if out_dir:
    os.makedirs(out_dir, exist_ok=True)

fig.savefig(args.out, dpi=180, bbox_inches="tight")
plt.close(fig)

print(f"\nSaved: {args.out}")
