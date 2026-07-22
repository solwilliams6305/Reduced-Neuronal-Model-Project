"""
simulate.py
-----------
Euler–Maruyama simulation of the full stochastic FitzHugh–Nagumo (FHN) system,
and the reduced model MFPT predictions for regime-map comparison.

Full 2D SDE:
    dv = (v - v^3/3 - w + I) dt + sigma dW
    dw = eps*(v + a - b*w) dt

Default parameters: I=-0.1, a=0.7, b=0.8
  Fixed point: (v, w) ≈ (-1.2563, -0.6954) on left stable branch (v < -1).
  Potential barrier DeltaU ≈ 0.0255 (at v_saddle ≈ -0.720 on middle branch).
  Left fold: v=-1, w_fold = I + 2/3 ≈ -0.767.

Reduction framework
-------------------
The reduced MFPT is the DETERMINISTIC slow-drift time to the fold:

    T_drift(eps) = integral_{w_FP}^{w_fold} dw / F_w(w),  F_w = eps*(v_s + a - bw)

This is the sigma→0 reduced prediction. For fixed eps, as sigma increases:
  - sigma << eps^(3/2): noise negligible, MFPT_full ≈ T_drift  (good reduction)
  - sigma >> eps^(3/2): noise drives v over the barrier far earlier than the
    deterministic drift reaches the fold — MFPT_full << T_drift  (bad reduction)

The regime map shows log(T_drift / MFPT_full) across (sigma, eps) space.
Positive = noise shortens MFPT (reduction overestimates) — the FAILURE regime.
Near zero = good reduction.
"""

import numpy as np
from scipy.optimize import brentq

# ---------------------------------------------------------------------------
# Default parameters
# ---------------------------------------------------------------------------
DEFAULTS = dict(I=-0.1, a=0.7, b=0.8)
V_FP = -1.2563
W_FP = -0.6954
V_SADDLE = -0.7196   # middle-branch saddle at w=W_FP (barrier top)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _v_left_newton(w_arr: np.ndarray, I: float, n_iter: int = 14) -> np.ndarray:
    """Vectorised Newton for left-branch v_s(w) < -1."""
    v = np.full_like(w_arr, -1.5)
    for _ in range(n_iter):
        f  = v - v**3 / 3.0 - w_arr + I
        df = 1.0 - v**2
        df = np.where(np.abs(df) < 1e-12, 1e-12, df)
        v  = v - f / df
    return v


def w_fold_left(I: float) -> float:
    """w-value of the left fold (v = -1)."""
    return -1.0 - (-1.0)**3 / 3.0 + I   # = I + 2/3


def potential_barrier(w_val: float, I: float, a: float = 0.7, b: float = 0.8) -> float:
    """
    DeltaU(w) = U(v_saddle, w) - U(v_min, w) at the given w.
    U(v, w) = -(v^2/2 - v^4/12) + (w - I)*v
    Returns 0 if no saddle exists (past fold).
    """
    def U(v): return -(v**2/2 - v**4/12) + (w_val - I)*v
    def vnull(v): return v - v**3/3 - w_val + I

    # v_min on left branch
    try:
        v_min = brentq(vnull, -3.0, -1.001)
    except ValueError:
        return 0.0

    # v_saddle on middle branch (exists only when cubic has 3 real roots)
    try:
        v_sad = brentq(vnull, -0.999, 0.999)
    except ValueError:
        return 0.0

    return max(U(v_sad) - U(v_min), 0.0)


# ---------------------------------------------------------------------------
# Full 2D SDE  (Euler–Maruyama)
# ---------------------------------------------------------------------------

def simulate_fhn(
    sigma: float,
    eps: float,
    I: float = -0.1,
    a: float = 0.7,
    b: float = 0.8,
    v0: float = V_FP,
    w0: float = W_FP,
    dt: float = 5e-3,
    T: float = 500.0,
    n_trajectories: int = 1000,
    threshold: float = 1.0,
    rng: np.random.Generator | None = None,
) -> dict:
    """
    Euler–Maruyama simulation of the full stochastic FHN.
    Firing: v >= threshold (right fold, v=+1).
    """
    if rng is None:
        rng = np.random.default_rng()

    n_steps   = int(T / dt)
    sqrt_dt   = np.sqrt(dt)
    n_save    = min(10, n_trajectories)
    save_every = max(1, n_steps // 5000)

    fpts       = np.full(n_trajectories, np.nan)
    w_spikes   = np.full(n_trajectories, np.nan)   # w at v >= threshold (firing)
    w_escapes  = np.full(n_trajectories, np.nan)   # w when v first crosses saddle (escape initiation)
    crossed    = np.zeros(n_trajectories, dtype=bool)
    escaped    = np.zeros(n_trajectories, dtype=bool)   # has crossed saddle

    v = np.full(n_trajectories, v0, dtype=float)
    w = np.full(n_trajectories, w0, dtype=float)

    v_save_list = [v[:n_save].copy()]
    w_save_list = [w[:n_save].copy()]
    t_save_list = [0.0]
    last_step   = 0

    early_exit_frac = 0.90   # stop once 90% of trajectories have fired

    for step in range(n_steps):
        noise = rng.standard_normal(n_trajectories)
        v = v + (v - v**3 / 3.0 - w + I) * dt + sigma * sqrt_dt * noise
        w = w + eps * (v + a - b * w) * dt

        # Track w at saddle crossing (escape initiation — v first crosses V_SADDLE upward)
        new_escape = (~escaped) & (v >= V_SADDLE)
        w_escapes[new_escape] = w[new_escape]
        escaped |= new_escape

        newly = (~crossed) & (v >= threshold)
        fpts[newly]     = (step + 1) * dt
        w_spikes[newly] = w[newly]        # record w at moment of firing
        crossed |= newly
        last_step = step + 1

        if step % save_every == 0:
            v_save_list.append(v[:n_save].copy())
            w_save_list.append(w[:n_save].copy())
            t_save_list.append((step + 1) * dt)

        if crossed.mean() >= early_exit_frac:
            break

    valid          = fpts[~np.isnan(fpts)]
    valid_w_spike  = w_spikes[~np.isnan(w_spikes)]
    valid_w_escape = w_escapes[~np.isnan(w_escapes)]
    return {
        "mfpt":             float(np.mean(valid))          if len(valid) > 0          else np.nan,
        "mfpt_std":         float(np.std(valid))           if len(valid) > 1          else np.nan,
        "fpts":             fpts,
        "fraction_fired":   float(np.mean(crossed)),
        "w_spike_mean":     float(np.mean(valid_w_spike))  if len(valid_w_spike) > 0  else np.nan,
        "w_spike_std":      float(np.std(valid_w_spike))   if len(valid_w_spike) > 1  else np.nan,
        "w_spikes":         w_spikes,
        "w_escape_mean":    float(np.mean(valid_w_escape)) if len(valid_w_escape) > 0 else np.nan,
        "w_escape_std":     float(np.std(valid_w_escape))  if len(valid_w_escape) > 1 else np.nan,
        "w_escapes":        w_escapes,
        "v_sample":         np.array(v_save_list),
        "w_sample":         np.array(w_save_list),
        "t":                np.array(t_save_list),
        "params":           dict(sigma=sigma, eps=eps, I=I, a=a, b=b,
                                 dt=dt, T=T, n_trajectories=n_trajectories,
                                 threshold=threshold),
    }


# ---------------------------------------------------------------------------
# Reduced model: deterministic drift time (sigma -> 0 limit)
# ---------------------------------------------------------------------------

def reduced_drift_time(
    eps: float,
    I: float  = -0.1,
    a: float  = 0.7,
    b: float  = 0.8,
    w0: float = W_FP,
    n_pts: int = 2000,
) -> float:
    """
    Deterministic reduced MFPT: T_drift = integral_{w0}^{w_fold} dw / F_w(w).

    F_w(w) = eps*(v_s(w) + a - bw) < 0 on this interval (drifts toward fold).
    Returns a positive number (integrates the absolute time).
    Scales exactly as 1/eps.
    """
    wf    = w_fold_left(I)
    w_arr = np.linspace(w0 - 1e-4, wf + 2e-3, n_pts)
    v_s   = _v_left_newton(w_arr, I)
    F_w   = eps * (v_s + a - b * w_arr)
    F_safe = np.where(np.abs(F_w) < 1e-12, -1e-12, F_w)
    integrand = 1.0 / F_safe
    # np.trapz works on NumPy < 2.0; np.trapezoid on >= 2.0
    try:
        t_drift = np.trapezoid(integrand, w_arr)
    except AttributeError:
        t_drift = np.trapz(integrand, w_arr)
    return float(abs(t_drift))


def simulate_reduced(
    sigma: float,
    eps:   float,
    I: float = -0.1,
    a: float = 0.7,
    b: float = 0.8,
    **kwargs,
) -> dict:
    """Wrapper: reduced MFPT = deterministic drift time (sigma-independent)."""
    t_drift = reduced_drift_time(eps=eps, I=I, a=a, b=b)
    return {
        "mfpt_reduced":   t_drift,
        "mfpt_std":       np.nan,
        "fraction_fired": 1.0,
        "w_fold":         w_fold_left(I),
        "params": dict(sigma=sigma, eps=eps, I=I, a=a, b=b),
    }
