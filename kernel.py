"""
kernel.py
---------
Core machinery for the FHN stochastic reduction project.

Three layers:

  1. MODEL       — FHN2D dataclass: parameters, geometry, bifurcation values
  2. SIMULATOR   — simulate_kernel(): EM integrator, two modes:
                     'first_passage'  → MFPT (excitable regime)
                     'isi_sequence'   → ISI array + CV (tonic spiking regime)
  3. ANALYSIS    — reduced_prediction(), sweep_grid(), failure_boundary(),
                     compute_cv(), log_ratio()

Old simulate.py / sweep.py remain as thin wrappers importing from here.
New code should import from kernel directly.

Regimes supported
-----------------
  Excitable     : I < I_snic   — single stable fixed point, noise-induced escape
  Tonic spiking : I_snic < I < I_hopf — stable limit cycle, noise perturbs ISI
  Resonator     : I > I_hopf   — stable spiral on right branch (future work)
"""

from __future__ import annotations

import warnings
import numpy as np
from dataclasses import dataclass, field
from typing import Literal
from multiprocessing import Pool, cpu_count
from scipy.optimize import brentq

# ---------------------------------------------------------------------------
# Helpers (geometry)
# ---------------------------------------------------------------------------

def _v_left_newton(w_arr: np.ndarray, I: float, n_iter: int = 14) -> np.ndarray:
    """Vectorised Newton solver for left-branch v_s(w) < -1."""
    v = np.full_like(w_arr, -1.5, dtype=float)
    for _ in range(n_iter):
        f  = v - v**3 / 3.0 - w_arr + I
        df = 1.0 - v**2
        df = np.where(np.abs(df) < 1e-12, 1e-12, df)
        v  = v - f / df
    return v


def _find_fixed_point(I: float, a: float, b: float) -> tuple[float, float]:
    """
    Numerically locate the left-branch fixed point (V_FP, W_FP).
    Returns (v_fp, w_fp).
    """
    def fp_eq(v):
        return v - v**3 / 3.0 - (v + a) / b + I
    try:
        v_fp = brentq(fp_eq, -3.0, -1.001)
    except ValueError:
        # Try middle or right branch
        v_fp = brentq(fp_eq, -3.0, 3.0)
    w_fp = (v_fp + a) / b
    return float(v_fp), float(w_fp)


def _find_saddle(w: float, I: float) -> float | None:
    """Middle-branch saddle v_sad(w). Returns None past the fold."""
    def f(v): return v - v**3 / 3.0 - w + I
    try:
        return float(brentq(f, -0.999, 0.999))
    except ValueError:
        return None


def _potential_barrier(w: float, I: float) -> float:
    """ΔU(w) = U(v_sad) - U(v_min).  Returns 0 past fold."""
    def U(v): return -v**2 / 2.0 + v**4 / 12.0 + (w - I) * v
    def vnull(v): return v - v**3 / 3.0 - w + I
    try:
        v_min = brentq(vnull, -3.0, -1.001)
        v_sad = brentq(vnull, -0.999,  0.999)
    except ValueError:
        return 0.0
    return float(max(U(v_sad) - U(v_min), 0.0))


# ---------------------------------------------------------------------------
# Layer 1 — Model
# ---------------------------------------------------------------------------

@dataclass
class FHN2D:
    """
    FitzHugh-Nagumo model parameters and derived geometry.

    Parameters
    ----------
    I : float   External current.  Controls regime:
                  I < I_snic  → excitable
                  I_snic < I < I_hopf → tonic spiking
                  I > I_hopf  → resonator (stable spiral, right branch)
    a, b : float  Recovery parameters (defaults: 0.7, 0.8)

    Derived quantities are computed lazily on first access via properties.
    """
    I: float = -0.1
    a: float =  0.7
    b: float =  0.8

    # ---- geometry (computed once, cached) ----------------------------------

    @property
    def fixed_point(self) -> tuple[float, float]:
        """(V_FP, W_FP) — stable fixed point location."""
        return _find_fixed_point(self.I, self.a, self.b)

    @property
    def V_FP(self) -> float:
        return self.fixed_point[0]

    @property
    def W_FP(self) -> float:
        return self.fixed_point[1]

    @property
    def w_fold_left(self) -> float:
        """w-value of left fold (v = -1): w = I - 2/3."""
        return self.I - 2.0 / 3.0

    @property
    def w_fold_right(self) -> float:
        """w-value of right fold (v = +1): w = I + 2/3."""
        return self.I + 2.0 / 3.0

    @property
    def V_SADDLE(self) -> float:
        """Middle-branch saddle at W_FP."""
        vs = _find_saddle(self.W_FP, self.I)
        return vs if vs is not None else float("nan")

    @property
    def delta_U(self) -> float:
        """Potential barrier ΔU at the fixed point W_FP."""
        return _potential_barrier(self.W_FP, self.I)

    def barrier_at(self, w: float) -> float:
        """ΔU(w) — barrier height as a function of slow variable."""
        return _potential_barrier(w, self.I)

    # ---- bifurcation values ------------------------------------------------
    #
    # Bifurcation structure (settled analytically — see hopf_analysis.py).
    # The Jacobian of the deterministic flow is
    #     J = [[1 - v², -1], [ε, -ε b]]
    #     trace(J) = 1 - v² - ε b ,   det(J) = ε(1 - b + b v²) > 0  (always).
    # A Hopf bifurcation requires trace = 0, i.e.
    #     v*_Hopf = ±√(1 - ε b).
    # As the fixed point climbs the middle branch with increasing I it crosses
    # BOTH roots, giving TWO Hopf points that bracket the tonic-spiking window:
    #
    #     I < I_fold_left            v* < -1                 excitable (left branch)
    #     I_fold_left < I < I_H1     -1 < v* < -√(1-εb)      stable spiral  (lower resonator)
    #     I_H1 < I < I_H2            -√(1-εb) < v* < √(1-εb) UNSTABLE → tonic limit cycle
    #     I_H2 < I < I_fold_right    √(1-εb) < v* < 1        stable spiral  (upper resonator)
    #     I > I_fold_right           v* > 1                  excitable (right branch)
    #
    # I_H1 (tonic birth) uses the NEGATIVE root; I_H2 (re-stabilisation) the positive.
    # det > 0 at both, and the Hopf frequency is ω = √det = √(ε(1 - ε b²)),
    # so the limit cycle is born with a FINITE period 2π/ω — confirming Hopf
    # (not SNIC, which would give a diverging period).

    def _I_at_v(self, v: float) -> float:
        """Current I that places the fixed point exactly at fast-variable value v."""
        w = (v + self.a) / self.b
        return float(w - v + v**3 / 3.0)

    def _I_at_fold_left(self) -> float:
        """I value that places fixed point exactly at left fold v* = -1."""
        return self._I_at_v(-1.0)

    def _I_at_fold_right(self) -> float:
        """I value that places fixed point exactly at right fold v* = +1."""
        return self._I_at_v(1.0)

    def I_hopf_at(self, eps: float) -> float:
        """
        Lower (tonic-birth) Hopf current at finite ε.
        trace(J)=0 ⟹ v*_Hopf = -√(1 - bε).  Tonic spiking begins here.
        """
        v_hopf = -np.sqrt(max(1.0 - self.b * eps, 0.0))
        return self._I_at_v(v_hopf)

    # explicit, self-documenting aliases
    def I_hopf_lower_at(self, eps: float) -> float:
        """Lower Hopf (tonic spiking is BORN here), v* = -√(1-bε)."""
        return self.I_hopf_at(eps)

    def I_hopf_upper_at(self, eps: float) -> float:
        """Upper Hopf (limit cycle DIES, spiral re-stabilises), v* = +√(1-bε)."""
        v_hopf = +np.sqrt(max(1.0 - self.b * eps, 0.0))
        return self._I_at_v(v_hopf)

    def hopf_frequency_at(self, eps: float) -> float:
        """Hopf angular frequency ω = √det = √(ε(1 - ε b²)) at trace=0."""
        return float(np.sqrt(max(eps * (1.0 - eps * self.b**2), 0.0)))

    def hopf_at(self, eps: float) -> dict:
        """
        Full description of the tonic Hopf transition at this ε.
        Returns I_lower, I_upper (the two Hopf currents bracketing tonic),
        ω, and the nascent limit-cycle period 2π/ω.
        """
        w = self.hopf_frequency_at(eps)
        return {
            "I_hopf_lower":  self.I_hopf_lower_at(eps),   # tonic spiking begins
            "I_hopf_upper":  self.I_hopf_upper_at(eps),   # tonic spiking ends
            "omega":         w,
            "nascent_period": float(2.0 * np.pi / w) if w > 0 else float("inf"),
            "v_hopf":        float(np.sqrt(max(1.0 - self.b * eps, 0.0))),
        }

    @property
    def I_hopf(self) -> float:
        """ε→0 limit of the (lower) Hopf current: fixed point at left fold v*=-1."""
        return self._I_at_fold_left()

    @property
    def I_snic(self) -> float:
        """
        DEPRECATED label.  This is the LEFT-FOLD current where the left-branch
        fixed point annihilates (excitable → lower-resonator boundary).
        It is NOT a SNIC: the tonic limit cycle is born later, at a Hopf
        (I_hopf_lower_at(eps)), with finite period.  Kept for backward compat.
        """
        return self._I_at_fold_left()

    @property
    def regime(self) -> str:
        """
        ε→0 regime label (uses Hopf points evaluated at a small reference ε).
        For finite ε use regime_at(eps), which is exact.
        """
        return self.regime_at(1e-3)

    def regime_at(self, eps: float) -> str:
        """
        Exact regime classification at timescale separation ε, based on where
        the operating fixed point sits relative to the folds and the two Hopf
        currents.  See the table above.
        """
        I_fold_l = self._I_at_fold_left()
        I_fold_r = self._I_at_fold_right()
        I_lo, I_hi = sorted((I_fold_l, I_fold_r))   # fold_left < fold_right normally
        I_H1 = self.I_hopf_lower_at(eps)
        I_H2 = self.I_hopf_upper_at(eps)

        if self.I < I_lo:
            return "excitable"            # left branch
        if self.I > I_hi:
            return "excitable"            # right branch
        # fixed point on the middle branch
        if self.I < I_H1:
            return "resonator"            # lower stable-spiral window
        if self.I < I_H2:
            return "tonic"                # limit cycle (fixed point unstable)
        return "resonator"               # upper stable-spiral window

    # ---- vector field ------------------------------------------------------

    def dv(self, v: np.ndarray, w: np.ndarray,
           sigma: float, noise: np.ndarray, dt: float) -> np.ndarray:
        """Fast variable increment (Euler-Maruyama)."""
        return (v - v**3 / 3.0 - w + self.I) * dt + sigma * np.sqrt(dt) * noise

    def dw(self, v: np.ndarray, w: np.ndarray,
           eps: float, dt: float) -> np.ndarray:
        """Slow variable increment (deterministic)."""
        return eps * (v + self.a - self.b * w) * dt

    def __repr__(self) -> str:
        v, w = self.fixed_point
        return (f"FHN2D(I={self.I}, a={self.a}, b={self.b} | "
                f"regime={self.regime}, "
                f"V_FP={v:.4f}, W_FP={w:.4f}, "
                f"ΔU={self.delta_U:.5f})")


# ---------------------------------------------------------------------------
# Layer 2 — Simulator
# ---------------------------------------------------------------------------

def simulate_kernel(
    model:          FHN2D,
    sigma:          float,
    eps:            float,
    mode:           Literal["first_passage", "isi_sequence"] = "first_passage",
    n_trajectories: int   = 1000,
    dt:             float = 5e-3,
    T:              float = 500.0,
    n_isi:          int   = 20,
    threshold:      float = 1.0,
    reset_v:        float | None = None,
    reset_w:        float | None = None,
    rng:            np.random.Generator | None = None,
    early_exit_frac: float = 0.90,
) -> dict:
    """
    Euler-Maruyama simulation of the full stochastic FHN.

    Parameters
    ----------
    model          : FHN2D instance (parameters + geometry)
    sigma          : noise amplitude
    eps            : timescale separation
    mode           : 'first_passage' — record first spike time per trajectory
                     'isi_sequence'  — record n_isi consecutive ISIs per trajectory
    n_trajectories : number of independent realisations
    dt             : Euler-Maruyama timestep
    T              : max simulation time
    n_isi          : number of ISIs to collect per trajectory (isi_sequence mode)
    threshold      : spike detection threshold (v >= threshold)
    reset_v/reset_w: reset point after spike in isi_sequence mode
                     (defaults to model fixed point)
    rng            : numpy random Generator (created if None)
    early_exit_frac: stop when this fraction of trajectories have fired

    Returns
    -------
    dict with keys depending on mode:

    Both modes:
        params, fraction_fired, v_sample, w_sample, t

    first_passage mode:
        mfpt, mfpt_std, fpts, w_spike_mean, w_spike_std, w_spikes,
        w_escape_mean, w_escape_std, w_escapes

    isi_sequence mode:
        isis          — (n_trajectories, n_isi) array, nan for uncollected
        isi_mean      — mean ISI across all valid ISIs
        isi_std       — std of ISI
        cv            — coefficient of variation = isi_std / isi_mean
        n_spikes      — total spikes per trajectory
    """
    if rng is None:
        rng = np.random.default_rng()

    # Reset point for isi_sequence mode
    v_reset = reset_v if reset_v is not None else model.V_FP
    w_reset = reset_w if reset_w is not None else model.W_FP

    n_steps    = int(T / dt)
    sqrt_dt    = np.sqrt(dt)
    n_save     = min(10, n_trajectories)
    save_every = max(1, n_steps // 5000)

    v = np.full(n_trajectories, v_reset, dtype=float)
    w = np.full(n_trajectories, w_reset, dtype=float)

    v_save_list = [v[:n_save].copy()]
    w_save_list = [w[:n_save].copy()]
    t_save_list = [0.0]

    # ---- first_passage bookkeeping ----
    fpts      = np.full(n_trajectories, np.nan)
    w_spikes  = np.full(n_trajectories, np.nan)
    w_escapes = np.full(n_trajectories, np.nan)
    crossed   = np.zeros(n_trajectories, dtype=bool)
    escaped   = np.zeros(n_trajectories, dtype=bool)

    # ---- isi_sequence bookkeeping ----
    isis         = np.full((n_trajectories, n_isi), np.nan)
    n_spikes_arr = np.zeros(n_trajectories, dtype=int)
    last_spike_t = np.full(n_trajectories, np.nan)
    in_refractory = np.zeros(n_trajectories, dtype=bool)  # brief post-spike pause

    v_saddle = model.V_SADDLE

    for step in range(n_steps):
        t = (step + 1) * dt
        noise = rng.standard_normal(n_trajectories)
        v = v + model.dv(v, w, sigma, noise, dt)
        w = w + model.dw(v, w, eps, dt)

        if mode == "first_passage":
            # Saddle crossing (escape initiation)
            new_escape = (~escaped) & (v >= v_saddle)
            w_escapes[new_escape] = w[new_escape]
            escaped |= new_escape

            # Threshold crossing (spike)
            newly = (~crossed) & (v >= threshold)
            fpts[newly]    = t
            w_spikes[newly] = w[newly]
            crossed |= newly

            if crossed.mean() >= early_exit_frac:
                break

        else:  # isi_sequence
            # Detect threshold crossings (not in refractory)
            newly = (~in_refractory) & (v >= threshold)

            if newly.any():
                idx = np.where(newly)[0]
                for k in idx:
                    if not np.isnan(last_spike_t[k]):
                        isi_idx = n_spikes_arr[k] - 1
                        if isi_idx < n_isi:
                            isis[k, isi_idx] = t - last_spike_t[k]
                    last_spike_t[k] = t
                    n_spikes_arr[k] += 1

                # Reset to fixed point
                v[newly] = v_reset
                w[newly] = w_reset
                in_refractory[newly] = True

            # Clear refractory after one step
            in_refractory[:] = False

            # Stop when all trajectories have collected enough ISIs
            if np.all(n_spikes_arr >= n_isi + 1):
                break

        if step % save_every == 0:
            v_save_list.append(v[:n_save].copy())
            w_save_list.append(w[:n_save].copy())
            t_save_list.append(t)

    # ---- assemble output ----
    out = {
        "params": dict(sigma=sigma, eps=eps, I=model.I, a=model.a, b=model.b,
                       dt=dt, T=T, mode=mode, n_trajectories=n_trajectories),
        "fraction_fired": float(np.mean(~np.isnan(fpts))) if mode == "first_passage"
                          else float(np.mean(n_spikes_arr >= 1)),
        "v_sample": np.array(v_save_list),
        "w_sample": np.array(w_save_list),
        "t":        np.array(t_save_list),
    }

    if mode == "first_passage":
        valid         = fpts[~np.isnan(fpts)]
        valid_w_spike  = w_spikes[~np.isnan(w_spikes)]
        valid_w_escape = w_escapes[~np.isnan(w_escapes)]
        out.update({
            "mfpt":           float(np.mean(valid))           if len(valid) > 0 else np.nan,
            "mfpt_std":       float(np.std(valid))            if len(valid) > 1 else np.nan,
            "fpts":           fpts,
            "w_spike_mean":   float(np.mean(valid_w_spike))   if len(valid_w_spike) > 0 else np.nan,
            "w_spike_std":    float(np.std(valid_w_spike))    if len(valid_w_spike) > 1 else np.nan,
            "w_spikes":       w_spikes,
            "w_escape_mean":  float(np.mean(valid_w_escape))  if len(valid_w_escape) > 0 else np.nan,
            "w_escape_std":   float(np.std(valid_w_escape))   if len(valid_w_escape) > 1 else np.nan,
            "w_escapes":      w_escapes,
        })

    else:  # isi_sequence
        valid_isis = isis[~np.isnan(isis)]
        cv = (float(np.std(valid_isis) / np.mean(valid_isis))
              if len(valid_isis) > 1 and np.mean(valid_isis) > 0 else np.nan)
        out.update({
            "isis":      isis,
            "isi_mean":  float(np.mean(valid_isis)) if len(valid_isis) > 0 else np.nan,
            "isi_std":   float(np.std(valid_isis))  if len(valid_isis) > 1 else np.nan,
            "cv":        cv,
            "n_spikes":  n_spikes_arr,
        })

    return out


# ---------------------------------------------------------------------------
# Layer 3 — Analysis
# ---------------------------------------------------------------------------

def reduced_prediction(
    model: FHN2D,
    eps:   float,
    mode:  Literal["first_passage", "isi_sequence"] = "first_passage",
    n_pts: int = 2000,
) -> float:
    """
    Reduced (sigma→0) prediction for the timescale diagnostic.

    first_passage : T_drift = ∫ dw / F_w(w)  from W_FP to w_fold
    isi_sequence  : T_cycle = deterministic limit cycle period (numerical ODE)
    """
    if mode == "first_passage":
        w0  = model.W_FP
        wf  = model.w_fold_left
        w_arr = np.linspace(w0 - 1e-4, wf + 2e-3, n_pts)
        v_s   = _v_left_newton(w_arr, model.I)
        F_w   = eps * (v_s + model.a - model.b * w_arr)
        F_safe = np.where(np.abs(F_w) < 1e-12, -1e-12, F_w)
        try:
            t_drift = float(abs(np.trapezoid(1.0 / F_safe, w_arr)))
        except AttributeError:
            t_drift = float(abs(np.trapz(1.0 / F_safe, w_arr)))
        return t_drift

    else:  # isi_sequence — integrate deterministic FHN limit cycle
        dt    = 1e-3
        T_max = 5000.0
        v, w  = 0.0, 0.0   # start near limit cycle
        threshold = 1.0
        spike_times = []
        t = 0.0
        last_above = False
        while t < T_max and len(spike_times) < 6:
            v = v + (v - v**3 / 3.0 - w + model.I) * dt
            w = w + eps * (v + model.a - model.b * w) * dt
            t += dt
            above = v >= threshold
            if above and not last_above:
                spike_times.append(t)
            last_above = above
        if len(spike_times) >= 3:
            isis = np.diff(spike_times)
            return float(np.mean(isis))
        return np.nan


def compute_cv(isis: np.ndarray) -> float:
    """Coefficient of variation of ISI distribution."""
    valid = isis[~np.isnan(isis)]
    if len(valid) < 2 or np.mean(valid) == 0:
        return np.nan
    return float(np.std(valid) / np.mean(valid))


def log_ratio(mfpt_full: float, mfpt_reduced: float) -> float:
    """
    log(T_reduced / MFPT_full).
    +ve = noise shortens MFPT → reduction FAILS.
    """
    if np.isnan(mfpt_full) or np.isnan(mfpt_reduced):
        return np.nan
    if mfpt_full <= 0 or mfpt_reduced <= 0:
        return np.nan
    return float(np.log(mfpt_reduced / mfpt_full))


# ---------------------------------------------------------------------------
# Parallel sweep worker (top-level for pickling)
# ---------------------------------------------------------------------------

def _kernel_worker(args: tuple) -> tuple:
    """Single grid point. Returns (i, j, stat_full, stat_reduced, w_spike, w_escape, cv)."""
    (i, j, sigma, eps, model_params, mode, n_trajectories,
     n_isi, dt, T_adaptive, threshold, global_seed) = args

    model = FHN2D(**model_params)
    rng   = np.random.default_rng([global_seed, i, j])

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = simulate_kernel(
            model=model, sigma=sigma, eps=eps, mode=mode,
            n_trajectories=n_trajectories, n_isi=n_isi,
            dt=dt, T=T_adaptive, threshold=threshold, rng=rng,
        )

    t_reduced = reduced_prediction(model, eps, mode=mode)

    if mode == "first_passage":
        stat_full = res["mfpt"]
        w_spike   = res.get("w_spike_mean",  np.nan)
        w_escape  = res.get("w_escape_mean", np.nan)
        cv        = np.nan
    else:
        stat_full = res["isi_mean"]
        w_spike   = np.nan
        w_escape  = np.nan
        cv        = res["cv"]

    return (i, j, stat_full, t_reduced, w_spike, w_escape, cv)


# ---------------------------------------------------------------------------
# Grid sweep
# ---------------------------------------------------------------------------

def sweep_grid(
    model:          FHN2D,
    sigma_vals:     np.ndarray,
    eps_vals:       np.ndarray,
    mode:           Literal["first_passage", "isi_sequence"] = "first_passage",
    n_trajectories: int   = 1000,
    n_isi:          int   = 20,
    dt:             float = 5e-3,
    T:              float = 300.0,
    threshold:      float = 1.0,
    seed:           int   = 42,
    verbose:        bool  = True,
    n_jobs:         int   = -1,
) -> dict:
    """
    Sweep (sigma, eps) grid for a given model and mode.

    Returns dict with sigma_vals, eps_vals, stat_full, stat_reduced,
    log_ratio, w_spike_mean, w_escape_mean, cv, drift_times, params.
    """
    n_s, n_e = len(sigma_vals), len(eps_vals)

    stat_full     = np.full((n_s, n_e), np.nan)
    stat_reduced  = np.full((n_s, n_e), np.nan)
    log_rat       = np.full((n_s, n_e), np.nan)
    w_spike_mean  = np.full((n_s, n_e), np.nan)
    w_escape_mean = np.full((n_s, n_e), np.nan)
    cv_arr        = np.full((n_s, n_e), np.nan)

    # Pre-compute reduced predictions
    drift_times = np.array([reduced_prediction(model, eps, mode=mode)
                             for eps in eps_vals])

    model_params = dict(I=model.I, a=model.a, b=model.b)

    work = []
    for i, sigma in enumerate(sigma_vals):
        for j, eps in enumerate(eps_vals):
            T_adaptive = float(np.clip(8.0 * drift_times[j], 50.0, T)) \
                         if not np.isnan(drift_times[j]) else T
            work.append((i, j, sigma, eps, model_params, mode,
                         n_trajectories, n_isi, dt, T_adaptive,
                         threshold, seed))

    total   = len(work)
    n_cores = cpu_count() if n_jobs == -1 else max(1, n_jobs)
    n_cores = min(n_cores, total)

    if verbose:
        print(f"  sweep_grid: {n_s}×{n_e} grid, mode={mode}, "
              f"I={model.I}, {n_cores} cores")

    done = 0

    if n_cores == 1:
        for item in work:
            i, j, sf, sr, ws, we, cv = _kernel_worker(item)
            stat_full[i, j]     = sf
            stat_reduced[i, j]  = sr
            log_rat[i, j]       = log_ratio(sf, sr)
            w_spike_mean[i, j]  = ws
            w_escape_mean[i, j] = we
            cv_arr[i, j]        = cv
            done += 1
            if verbose:
                print(f"  [{done}/{total}]  sigma={item[2]:.3f}  eps={item[3]:.4f}",
                      end="\r")
    else:
        with Pool(processes=n_cores) as pool:
            for result in pool.imap_unordered(_kernel_worker, work):
                i, j, sf, sr, ws, we, cv = result
                stat_full[i, j]     = sf
                stat_reduced[i, j]  = sr
                log_rat[i, j]       = log_ratio(sf, sr)
                w_spike_mean[i, j]  = ws
                w_escape_mean[i, j] = we
                cv_arr[i, j]        = cv
                done += 1
                if verbose:
                    print(f"  [{done}/{total}] done", end="\r")

    if verbose:
        print(f"\n  Done. {n_s}×{n_e} = {total} points.")

    return {
        "sigma_vals":    sigma_vals,
        "eps_vals":      eps_vals,
        "stat_full":     stat_full,
        "stat_reduced":  stat_reduced,
        "log_ratio":     log_rat,
        # legacy keys for backward compat with plot.py
        "mfpt_full":     stat_full,
        "mfpt_reduced":  stat_reduced,
        "w_spike_mean":  w_spike_mean,
        "w_escape_mean": w_escape_mean,
        "cv":            cv_arr,
        "drift_times":   drift_times,
        "params": dict(I=model.I, a=model.a, b=model.b,
                       mode=mode, n_trajectories=n_trajectories,
                       dt=dt, T=T, threshold=threshold, seed=seed),
    }


# ---------------------------------------------------------------------------
# Failure boundary extraction  (unchanged from sweep.py)
# ---------------------------------------------------------------------------

def failure_boundary(
    sweep_result: dict,
    delta_ratio:  float = 2.0,
) -> tuple[np.ndarray, np.ndarray]:
    """
    For each eps, smallest sigma where log_ratio > log(delta_ratio).
    Returns (eps_boundary, sigma_boundary).
    """
    sigma_vals = sweep_result["sigma_vals"]
    eps_vals   = sweep_result["eps_vals"]
    log_rat    = sweep_result["log_ratio"]
    threshold  = np.log(delta_ratio)

    eps_bd, sigma_bd = [], []
    for j, eps in enumerate(eps_vals):
        col   = log_rat[:, j]
        above = np.where(~np.isnan(col) & (col > threshold))[0]
        if len(above) > 0:
            idx = above[0]
            if idx > 0 and not np.isnan(col[idx - 1]):
                frac = (threshold - col[idx - 1]) / (col[idx] - col[idx - 1])
                sig  = sigma_vals[idx-1] + frac * (sigma_vals[idx] - sigma_vals[idx-1])
            else:
                sig = sigma_vals[idx]
            eps_bd.append(eps)
            sigma_bd.append(sig)

    return np.array(eps_bd), np.array(sigma_bd)
