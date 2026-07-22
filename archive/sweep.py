"""
sweep.py
--------
Parameter sweep over (sigma, eps) grid.

Discrepancy metric
------------------
We use the log-ratio:

    D(sigma, eps) = log( T_drift(eps) / MFPT_full(sigma, eps) )

  D >> 0 : noise dominates, MFPT_full << T_drift  →  REDUCTION FAILS
           (sigma >> eps^(3/2), above BG curve)
  D ~ 0  : good reduction  (sigma << eps^(3/2), below BG curve)
  D < 0  : rare; full model slower than drift (very small sigma + finite-T bias)

The FAILURE BOUNDARY is defined by D = log(delta_ratio), e.g. delta_ratio = 2
(full MFPT is less than half the drift time).
"""

from __future__ import annotations

import numpy as np
import warnings
from typing import Callable
from multiprocessing import Pool, cpu_count

from simulate import simulate_fhn, simulate_reduced, reduced_drift_time


# ---------------------------------------------------------------------------
# Worker function for parallel sweep (must be top-level for pickling)
# ---------------------------------------------------------------------------

def _sweep_worker(args: tuple) -> tuple:
    """
    Single grid point evaluation. Returns (i, j, mfpt, w_spike, w_escape).
    Each worker gets its own RNG seeded from (global_seed, i, j) for
    reproducibility and independence.
    """
    (i, j, sigma, eps, td, I, a, b, n_trajectories,
     dt, T_adaptive, threshold, global_seed) = args

    # Independent RNG per grid point — reproducible and non-overlapping
    rng = np.random.default_rng([global_seed, i, j])

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = simulate_fhn(
            sigma=sigma, eps=eps, I=I, a=a, b=b,
            n_trajectories=n_trajectories, dt=dt, T=T_adaptive,
            threshold=threshold, rng=rng,
        )

    return (
        i, j,
        res["mfpt"],
        res.get("w_spike_mean",  np.nan),
        res.get("w_escape_mean", np.nan),
    )


# ---------------------------------------------------------------------------
# Discrepancy metric
# ---------------------------------------------------------------------------

def log_ratio(mfpt_full: float, mfpt_reduced: float) -> float:
    """
    log(T_drift / MFPT_full).
    Positive = noise shortens MFPT (reduction overestimates) = FAILURE regime.
    Returns nan if either is nan or non-positive.
    """
    if np.isnan(mfpt_full) or np.isnan(mfpt_reduced):
        return np.nan
    if mfpt_full <= 0 or mfpt_reduced <= 0:
        return np.nan
    return float(np.log(mfpt_reduced / mfpt_full))


def relative_error(mfpt_full: float, mfpt_reduced: float) -> float:
    """
    Symmetric relative error: |MFPT_full - T_drift| / T_drift.
    Kept for backward compatibility.
    """
    if np.isnan(mfpt_full) or np.isnan(mfpt_reduced):
        return np.nan
    if mfpt_reduced == 0:
        return np.nan
    return abs(mfpt_full - mfpt_reduced) / mfpt_reduced


# ---------------------------------------------------------------------------
# Grid sweep
# ---------------------------------------------------------------------------

def sweep_grid(
    sigma_vals: np.ndarray,
    eps_vals:   np.ndarray,
    I: float = -0.1,
    a: float = 0.7,
    b: float = 0.8,
    n_trajectories: int = 1000,
    dt: float = 5e-3,
    T: float = 300.0,
    threshold: float = 1.0,
    seed: int = 42,
    verbose: bool = True,
    n_jobs: int = -1,
    callback: Callable | None = None,
) -> dict:
    """
    Sweep over (sigma, eps) grid.

    n_jobs : number of parallel workers.
             -1 = use all available cores (default)
              1 = serial (useful for debugging)
              N = use N cores

    Returns
    -------
    dict with sigma_vals, eps_vals, mfpt_full, mfpt_reduced,
    log_ratio, rel_error, w_spike_mean, w_escape_mean, drift_times, params.
    """
    n_s  = len(sigma_vals)
    n_e  = len(eps_vals)

    mfpt_full      = np.full((n_s, n_e), np.nan)
    mfpt_reduced   = np.full((n_s, n_e), np.nan)
    log_rat        = np.full((n_s, n_e), np.nan)
    rel_err        = np.full((n_s, n_e), np.nan)
    w_spike_mean   = np.full((n_s, n_e), np.nan)
    w_escape_mean  = np.full((n_s, n_e), np.nan)

    # Pre-compute drift times (sigma-independent)
    drift_times = np.array([reduced_drift_time(eps=eps, I=I, a=a, b=b)
                             for eps in eps_vals])

    # Build flat list of work items
    work = []
    for i, sigma in enumerate(sigma_vals):
        for j, eps in enumerate(eps_vals):
            T_adaptive = float(np.clip(8.0 * drift_times[j], 50.0, T))
            work.append((i, j, sigma, eps, drift_times[j],
                         I, a, b, n_trajectories, dt, T_adaptive,
                         threshold, seed))

    total   = len(work)
    n_cores = cpu_count() if n_jobs == -1 else max(1, n_jobs)
    n_cores = min(n_cores, total)

    if verbose:
        print(f"  Parallelising over {n_cores} cores  ({total} grid points)")

    done = 0

    if n_cores == 1:
        # Serial path — easier to debug
        for item in work:
            i, j, mfpt, ws, we = _sweep_worker(item)
            td = drift_times[j]
            mfpt_full[i, j]     = mfpt
            mfpt_reduced[i, j]  = td
            log_rat[i, j]       = log_ratio(mfpt, td)
            rel_err[i, j]       = relative_error(mfpt, td)
            w_spike_mean[i, j]  = ws
            w_escape_mean[i, j] = we
            done += 1
            if verbose:
                print(f"  [{done}/{total}]  sigma={item[2]:.3f}  eps={item[3]:.4f}",
                      end="\r")
    else:
        with Pool(processes=n_cores) as pool:
            for result in pool.imap_unordered(_sweep_worker, work):
                i, j, mfpt, ws, we = result
                td = drift_times[j]
                mfpt_full[i, j]     = mfpt
                mfpt_reduced[i, j]  = td
                log_rat[i, j]       = log_ratio(mfpt, td)
                rel_err[i, j]       = relative_error(mfpt, td)
                w_spike_mean[i, j]  = ws
                w_escape_mean[i, j] = we
                done += 1
                if verbose:
                    print(f"  [{done}/{total}] done", end="\r")

    if verbose:
        print(f"\nDone. Grid: {n_s}×{n_e} = {total} points.")

    return {
        "sigma_vals":    sigma_vals,
        "eps_vals":      eps_vals,
        "mfpt_full":     mfpt_full,
        "mfpt_reduced":  mfpt_reduced,
        "log_ratio":     log_rat,
        "rel_error":     rel_err,
        "w_spike_mean":  w_spike_mean,
        "w_escape_mean": w_escape_mean,
        "drift_times":   drift_times,
        "params": dict(I=I, a=a, b=b, n_trajectories=n_trajectories,
                       dt=dt, T=T, threshold=threshold, seed=seed),
    }


# ---------------------------------------------------------------------------
# Failure boundary: smallest sigma where log_ratio > log(delta_ratio)
# ---------------------------------------------------------------------------

def failure_boundary(
    sweep_result: dict,
    delta_ratio: float = 2.0,
) -> tuple[np.ndarray, np.ndarray]:
    """
    For each eps, find the smallest sigma where log_ratio > log(delta_ratio)
    (i.e., MFPT_full < T_drift / delta_ratio — noise shortens MFPT by more than 2x).

    Returns (eps_boundary, sigma_boundary).
    """
    sigma_vals = sweep_result["sigma_vals"]
    eps_vals   = sweep_result["eps_vals"]
    log_rat    = sweep_result["log_ratio"]
    threshold  = np.log(delta_ratio)

    eps_bd   = []
    sigma_bd = []

    for j, eps in enumerate(eps_vals):
        col   = log_rat[:, j]   # varies over sigma
        above = np.where(~np.isnan(col) & (col > threshold))[0]
        if len(above) > 0:
            idx = above[0]
            if idx > 0 and not np.isnan(col[idx - 1]):
                frac = (threshold - col[idx - 1]) / (col[idx] - col[idx - 1])
                sig  = sigma_vals[idx - 1] + frac * (sigma_vals[idx] - sigma_vals[idx - 1])
            else:
                sig = sigma_vals[idx]
            eps_bd.append(eps)
            sigma_bd.append(sig)

    return np.array(eps_bd), np.array(sigma_bd)


# ---------------------------------------------------------------------------
# BG reference curve and power-law fit
# ---------------------------------------------------------------------------

def berglund_gentz_curve(eps_vals: np.ndarray, C: float = 1.0) -> np.ndarray:
    """Reference scaling: sigma = C * eps^(3/2)."""
    return C * eps_vals ** 1.5


def kramers_crossover_curve(
    eps_vals: np.ndarray,
    delta_U: float,
) -> np.ndarray:
    """
    Crossover curve: sigma_cross(eps) = sqrt(Delta_U / log(1/eps)).

    Above this: Kramers escape dominates (noise faster than drift).
    Below this: drift dominates, fold-layer mechanism relevant.

    Derived by setting Kramers escape time ~ exp(Delta_U / sigma^2)
    equal to deterministic drift time ~ 1/eps:
        Delta_U / sigma^2 = log(1/eps)  =>  sigma = sqrt(Delta_U / log(1/eps))
    """
    log_term = np.log(1.0 / np.where(eps_vals > 0, eps_vals, np.nan))
    return np.sqrt(delta_U / log_term)


def fit_power_law(
    eps_bd: np.ndarray, sigma_bd: np.ndarray
) -> tuple[float, float]:
    """
    Fit log(sigma) = log(C) + alpha*log(eps) via OLS.
    Returns (alpha, C).
    """
    if len(eps_bd) < 2:
        return np.nan, np.nan
    log_eps = np.log(eps_bd)
    log_sig = np.log(sigma_bd)
    alpha, log_C = np.polyfit(log_eps, log_sig, 1)
    return float(alpha), float(np.exp(log_C))


# ---------------------------------------------------------------------------
# Diagnostic 1: local exponent analysis
# ---------------------------------------------------------------------------

def local_exponent(
    eps_bd: np.ndarray,
    sigma_bd: np.ndarray,
    window: int = 4,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Fit power law in a sliding window of `window` consecutive boundary points.
    Returns (eps_centres, local_alphas).

    If the true exponent is 0.5 (fold-layer), local_alpha -> 0.5 as eps -> 0.
    If it stays flat at ~0.43, it is structural rather than a finite-eps artefact.

    Requires at least `window` boundary points; returns empty arrays otherwise.
    """
    n = len(eps_bd)
    if n < window:
        return np.array([]), np.array([])

    eps_centres   = np.empty(n - window + 1)
    local_alphas  = np.empty(n - window + 1)

    for k in range(n - window + 1):
        eps_w   = eps_bd[k : k + window]
        sigma_w = sigma_bd[k : k + window]
        alpha, _ = fit_power_law(eps_w, sigma_w)
        eps_centres[k]  = np.sqrt(eps_w[0] * eps_w[-1])   # geometric centre
        local_alphas[k] = alpha

    return eps_centres, local_alphas
