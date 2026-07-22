"""
validate_commitment.py
----------------------
Validates the non-homogeneous Poisson escape model against empirical
w_escape distributions from the full stochastic FHN simulation.

Framework
---------
Treating w as frozen at each instant, the Kramers escape rate out of the
left well is:

    k(w) = A(w) · exp(-2·ΔU(w) / σ²)

where:
    ΔU(w)  = U(v_sad(w)) - U(v_min(w))   [barrier height at this w]
    A(w)   = (1/2π)·√(|U''(v_min)| · |U''(v_sad)|)  [Kramers prefactor]
    U(v)   = -v²/2 + v⁴/12 + (w - I)·v

Changing variable from time t to slow variable w via ẇ = ε·F_w(w):

    λ(w)   = k(w) / (ε · |F_w(w)|)       [escape rate per unit w-drift]

Survival probability (prob. of not having escaped by the time slow var = w):

    S(w)   = exp( -∫_{W_FP}^{w} λ(w') dw' )

Commitment density (prob. density of escape occurring at w):

    p(w)   = λ(w) · S(w)

If the model is correct, p(w) should match the empirical histogram of
w_escape values recorded in simulate_fhn().

Usage
-----
    python validate_commitment.py                   # default test pairs
    python validate_commitment.py --n-traj 2000     # more trajectories
    python validate_commitment.py --quick           # fewer pairs, fast check
"""

from __future__ import annotations
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.integrate import cumulative_trapezoid

# Make sure simulate.py is importable from same directory
sys.path.insert(0, os.path.dirname(__file__))
from simulate import (
    simulate_fhn, potential_barrier,
    w_fold_left, W_FP, V_FP, V_SADDLE,
    _v_left_newton,
)

os.makedirs("figures", exist_ok=True)

# ---------------------------------------------------------------------------
# Default parameters
# ---------------------------------------------------------------------------
I, a, b = -0.1, 0.7, 0.8

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--n-traj", type=int, default=3000,
                    help="Trajectories per (sigma, eps) pair (default 3000)")
parser.add_argument("--quick", action="store_true",
                    help="Fewer test pairs for fast iteration")
parser.add_argument("--seed", type=int, default=42)
args = parser.parse_args()

# ---------------------------------------------------------------------------
# (sigma, eps) pairs to validate
# All sit in the Kramers-dominated region (above sigma_cross, below fold layer)
# ---------------------------------------------------------------------------
if args.quick:
    TEST_PAIRS = [
        (0.20, 0.05),
        (0.30, 0.05),
        (0.25, 0.08),
    ]
else:
    TEST_PAIRS = [
        (0.15, 0.03),   # near boundary, small eps
        (0.20, 0.05),   # core Kramers region
        (0.28, 0.05),   # higher sigma, same eps
        (0.25, 0.08),   # larger eps
        (0.35, 0.08),   # high sigma, larger eps
        (0.30, 0.12),   # well inside failure region
    ]

# ---------------------------------------------------------------------------
# Analytical machinery
# ---------------------------------------------------------------------------

def _v_saddle(w: float) -> float | None:
    """Middle-branch saddle v_sad(w). Returns None past the fold."""
    def f(v): return v - v**3 / 3.0 - w + I
    try:
        return brentq(f, -0.999, 0.999)
    except ValueError:
        return None


def _v_min(w: float) -> float | None:
    """Left-branch minimum v_min(w). Returns None past the fold."""
    def f(v): return v - v**3 / 3.0 - w + I
    try:
        return brentq(f, -3.0, -1.001)
    except ValueError:
        return None


def U(v: float, w: float) -> float:
    """Potential U(v, w) = -v²/2 + v⁴/12 + (w - I)·v."""
    return -v**2 / 2.0 + v**4 / 12.0 + (w - I) * v


def U_curvature(v: float) -> float:
    """U''(v) = -1 + v²."""
    return v**2 - 1.0


def kramers_prefactor(w: float) -> float | None:
    """
    A(w) = (1/2π) · √(|U''(v_min)| · |U''(v_sad)|)

    Returns None if no saddle exists (past fold).
    """
    vm = _v_min(w)
    vs = _v_saddle(w)
    if vm is None or vs is None:
        return None
    curv_min = abs(U_curvature(vm))   # > 0 since |v_min| > 1
    curv_sad = abs(U_curvature(vs))   # > 0 since |v_sad| < 1
    return (1.0 / (2.0 * np.pi)) * np.sqrt(curv_min * curv_sad)


def delta_U(w: float) -> float:
    """ΔU(w) = U(v_sad) - U(v_min). Returns 0 past fold."""
    vm = _v_min(w)
    vs = _v_saddle(w)
    if vm is None or vs is None:
        return 0.0
    return max(U(vs, w) - U(vm, w), 0.0)


def F_w(w: float, v_s: float) -> float:
    """Slow flow: ẇ = ε · F_w,  F_w = v_s + a - b·w."""
    return v_s + a - b * w


# ---------------------------------------------------------------------------
# Compute λ(w), S(w), p(w) on a grid
# ---------------------------------------------------------------------------

def commitment_model(
    sigma: float,
    eps: float,
    n_pts: int = 2000,
) -> dict:
    """
    Compute the non-homogeneous Poisson commitment density p(w).

    The slow variable w drifts DOWNWARD from W_FP toward w_fold (w_fold < W_FP).
    F_w = v_s + a - b*w < 0 throughout (negative drift), and F_w = 0 exactly
    at the fixed point W_FP.

    Integration convention
    ----------------------
    We parameterise by arc-length along the trajectory: s = W_FP - w ≥ 0.
    All integrals are against ds = -dw > 0, so w_grid is DECREASING and
    step sizes dw_abs = -diff(w_grid) > 0.  The cumulative hazard is:

        H(k) = Σ_{j<k} λ(w_j) * dw_abs_j   (trapezoid rule, all positive)

    Singularity at W_FP
    -------------------
    Near W_FP, |F_w| ≈ |α|*(W_FP - w) where α = d(F_w)/dw|_{W_FP} ≈ -2.524.
    This gives λ ~ 1/|w - W_FP|, a logarithmic singularity in H.
    We start the grid at W_FP - 1e-3 (not exactly W_FP) and clip |F_w| from
    below at 1e-4 to regularise; the contribution missed is negligible for
    the σ values we care about.

    Returns dict with:
        w_grid   : DECREASING w values from W_FP-1e-3 to w_fold-1e-3
        lambda_w : λ(w) = k(w)/(ε|F_w|)  [rate per unit |dw|, always ≥ 0]
        S_w      : survival probability S(w) ∈ [0,1]
        p_w      : normalised commitment density (per unit |dw|)
        dU_w     : ΔU(w) barrier profile
        A_w      : Kramers prefactor A(w)
        norm     : ∫ p_raw |dw| (should be ≈ 1 in Kramers regime)
        w_star   : mode of p(w)  (most probable escape w)
        mean_w   : mean of p(w)
    """
    w_fold = w_fold_left(I)

    # DECREASING grid: starts just below W_FP (skip singularity), ends before fold
    w_start = W_FP  - 1e-3
    w_end   = w_fold - 1e-3
    w_grid  = np.linspace(w_start, w_end, n_pts)   # w_start > w_end, decreasing

    # Step sizes along trajectory (all positive)
    dw_abs = -np.diff(w_grid)                        # length n_pts-1, > 0

    # Vectorised left-branch v_s(w)
    v_s_arr = _v_left_newton(w_grid, I)

    # ΔU and Kramers prefactor — explicit float64
    dU_arr = np.array([delta_U(w) for w in w_grid], dtype=np.float64)

    A_arr = np.empty(n_pts, dtype=np.float64)
    for k, w in enumerate(w_grid):
        val = kramers_prefactor(w)
        A_arr[k] = (float(val) if (val is not None and not np.isnan(float(val)))
                    else (A_arr[k - 1] if k > 0 else 1.0))

    # Slow-flow magnitude |F_w(w)| — clip away from zero to regularise W_FP singularity
    Fw_abs = np.abs(v_s_arr + a - b * w_grid)
    Fw_abs = np.maximum(Fw_abs, 1e-4)

    # Kramers rate k(w) = A(w)·exp(-2ΔU(w)/σ²)
    k_arr = A_arr * np.exp(np.clip(-2.0 * dU_arr / sigma**2, -700, 0))

    # λ(w) = k(w) / (ε·|F_w|)  — positive, rate per unit |dw|
    lambda_arr = k_arr / (eps * Fw_abs)

    # Cumulative hazard H(k) = ∫_0^{s_k} λ ds (trapezoid, all positive)
    lam_mid = 0.5 * (lambda_arr[:-1] + lambda_arr[1:])
    H = np.concatenate([[0.0], np.cumsum(lam_mid * dw_abs)])

    # Survival and commitment density
    S_arr = np.exp(-H)                    # guaranteed ∈ (0, 1]
    p_arr = lambda_arr * S_arr            # guaranteed ≥ 0

    # Normalise: ∫ p |dw| using trapezoidal rule over decreasing grid
    p_mid = 0.5 * (p_arr[:-1] + p_arr[1:])
    norm  = float(np.sum(p_mid * dw_abs))
    p_norm = p_arr / max(norm, 1e-12)

    # Summary statistics (integrals against |dw|)
    w_star = w_grid[np.argmax(p_norm)]
    # p_norm already integrates to 1, so no further division needed
    wpm    = 0.5 * (w_grid[:-1] * p_norm[:-1] + w_grid[1:] * p_norm[1:])
    mean_w = float(np.sum(wpm * dw_abs))

    return {
        "w_grid":   w_grid,
        "lambda_w": lambda_arr,
        "S_w":      S_arr,
        "p_w":      p_norm,
        "p_raw":    p_arr,
        "norm":     norm,
        "dU_w":     dU_arr,
        "A_w":      A_arr,
        "w_star":   w_star,
        "mean_w":   mean_w,
    }


# ---------------------------------------------------------------------------
# Run simulations and validate
# ---------------------------------------------------------------------------

def run_validation(sigma: float, eps: float, n_traj: int, seed: int) -> dict:
    """Run simulate_fhn and return w_escape array for this (sigma, eps) pair."""
    rng = np.random.default_rng(seed)
    T   = max(200.0, 8.0 / eps)   # long enough to capture slow escapes
    res = simulate_fhn(
        sigma=sigma, eps=eps, I=I, a=a, b=b,
        v0=V_FP, w0=W_FP,
        n_trajectories=n_traj,
        T=T, rng=rng,
    )
    w_esc = res["w_escapes"]
    w_esc_valid = w_esc[~np.isnan(w_esc)]
    return {
        "w_escapes":      w_esc_valid,
        "fraction_fired": res["fraction_fired"],
        "mfpt":           res["mfpt"],
        "n_valid":        len(w_esc_valid),
    }


# ---------------------------------------------------------------------------
# Figure: one panel per (sigma, eps) pair
# ---------------------------------------------------------------------------

def plot_single_validation(
    sigma: float,
    eps: float,
    model: dict,
    sim: dict,
    ax: plt.Axes,
    tau_v: float | None = None,
) -> dict:
    """
    Overlay non-homogeneous Poisson p(w) and Gaussian OU prediction
    against empirical w_escape histogram.
    """
    from scipy.stats import norm as scipy_norm

    w_fold = w_fold_left(I)
    w_esc  = sim["w_escapes"]

    g_mean, g_std = gaussian_w_escape(sigma, eps, tau_v=tau_v)
    x_lo = min(W_FP - 4*g_std, w_fold - 0.005)
    x_hi = max(W_FP + 4*g_std, W_FP + 0.01)
    w_plot = np.linspace(x_lo, x_hi, 500)

    # Empirical histogram
    if len(w_esc) >= 10:
        ax.hist(w_esc, bins=40, density=True, alpha=0.40,
                color="#2166ac", label=f"Empirical  (N={sim['n_valid']})")

    # Non-homogeneous Poisson p(w)
    ax.plot(model["w_grid"], model["p_w"], color="#d6604d", lw=2.0,
            label=f"Poisson  [norm={model['norm']:.2f}]")

    # Gaussian OU prediction
    gauss = gaussian_pdf(w_plot, g_mean, g_std)
    ax.plot(w_plot, gauss, color="#1a9641", lw=2.0, ls="--",
            label=f"Gaussian  $\\mu=W_{{FP}}$, $\\sigma_w$={g_std:.4f}")

    # Reference lines
    ax.axvline(W_FP,  color="gray",  lw=1.0, ls="--", alpha=0.6,
               label=f"$W_{{FP}}$={W_FP:.3f}")
    if len(w_esc) >= 5:
        ax.axvline(np.mean(w_esc), color="#2166ac", lw=1.0, ls=":",
                   label=f"$\\langle w_{{sim}}\\rangle$={np.mean(w_esc):.3f}")

    ax.set_xlabel("$w$ at escape", fontsize=10)
    ax.set_ylabel("Density", fontsize=10)
    ax.set_title(
        f"$\\sigma={sigma}$,  $\\varepsilon={eps}$\n"
        f"fired={sim['fraction_fired']:.0%},  MFPT={sim['mfpt']:.1f}",
        fontsize=9)
    ax.legend(fontsize=7, loc="upper left")
    ax.set_xlim(x_lo, x_hi)

    # KS for Poisson model
    ks_poisson = np.nan
    if len(w_esc) >= 20:
        w_g   = model["w_grid"]
        w_inc = w_g[::-1]
        p_inc = model["p_w"][::-1]
        dw_i  = np.diff(w_inc)
        cdf   = np.concatenate([[0.0],
                                np.cumsum(0.5*(p_inc[:-1]+p_inc[1:])*dw_i)])
        cdf  /= max(cdf[-1], 1e-12)
        s_obs = np.sort(w_esc)
        ks_poisson = float(np.max(np.abs(
            np.arange(1, len(s_obs)+1)/len(s_obs) -
            np.interp(s_obs, w_inc, cdf, left=0., right=1.))))

    # KS for Gaussian model
    ks_gaussian = np.nan
    if len(w_esc) >= 20:
        s_obs = np.sort(w_esc)
        pred_cdf = scipy_norm.cdf(s_obs, loc=g_mean, scale=g_std)
        ks_gaussian = float(np.max(np.abs(
            np.arange(1, len(s_obs)+1)/len(s_obs) - pred_cdf)))

    return {"ks_poisson": ks_poisson, "ks_gaussian": ks_gaussian,
            "g_mean": g_mean, "g_std": g_std}



# ---------------------------------------------------------------------------
# Gaussian OU prediction for w_escape
# ---------------------------------------------------------------------------

def estimate_tau_v(sigma: float, eps: float, n_traj: int = 200,
                   T: float = 50.0, dt: float = 5e-3) -> float:
    """
    Estimate the autocorrelation time τ_v of v fluctuations near the fixed
    point by running short trajectories at low noise (so the system stays
    near the fixed point) and computing the v autocorrelation function.

    τ_v = ∫_0^∞ C(t)/C(0) dt  ≈  Σ_k C(k·dt)/C(0) · dt

    This is the timescale over which individual v fluctuations persist before
    reversing — the correct clock for how much w gets jostled before escape.
    """
    rng    = np.random.default_rng(999)
    n_steps = int(T / dt)
    sqrt_dt = np.sqrt(dt)

    # Use a small sigma so the system stays near FP (don't want escape events)
    sigma_probe = min(sigma * 0.3, 0.05)

    v = np.full(n_traj, V_FP, dtype=float)
    w = np.full(n_traj, W_FP, dtype=float)

    v_trace = np.empty((n_steps, n_traj))
    for step in range(n_steps):
        noise = rng.standard_normal(n_traj)
        v = v + (v - v**3/3 - w + I) * dt + sigma_probe * sqrt_dt * noise
        w = w + eps * (v + a - b*w) * dt
        # Reset any trajectories that escaped (keep near FP)
        escaped = v >= -0.5
        v[escaped] = V_FP
        w[escaped] = W_FP
        v_trace[step] = v

    # Mean-subtract
    v_fluct = v_trace - V_FP

    # Autocorrelation via FFT (averaged over trajectories)
    max_lag = min(n_steps // 4, int(20.0 / dt))
    acf = np.zeros(max_lag)
    for k in range(n_traj):
        x = v_fluct[:, k]
        c0 = float(np.mean(x**2))
        if c0 < 1e-12:
            continue
        for lag in range(max_lag):
            acf[lag] += float(np.mean(x[:n_steps-lag] * x[lag:n_steps])) / c0
    acf /= n_traj

    # Integrate until ACF drops below 1/e or goes negative
    tau = 0.0
    for lag in range(max_lag):
        if acf[lag] < 0:
            break
        tau += acf[lag] * dt

    return max(tau, dt)


def gaussian_w_escape(sigma: float, eps: float,
                      tau_v: float | None = None) -> tuple[float, float]:
    """
    Predicted w_escape distribution from the OU approximation.

    The w equation near W_FP is an OU process driven by v-noise:
        dw̃ = -bε·w̃·dt + εσ·dW_t

    Stationary variance (applies when T_esc >> 1/(bε)):
        Var(w̃) = εσ²τ_v / (2b)

    Short-time variance (applies when T_esc << 1/(bε)):
        Var(w̃) ≈ ε²σ²τ_v·T_esc

    Full formula interpolating between both regimes:
        Var(w̃) = εσ²τ_v/(2b) · (1 - exp(-2bε·T_esc))

    where T_esc = 1/k(W_FP) is the mean Kramers escape time.
    """
    if tau_v is None:
        tau_v = estimate_tau_v(sigma, eps)

    # Kramers escape rate at W_FP
    k0    = kramers_prefactor(W_FP)
    k0    = float(k0) if k0 is not None else 1.0
    k_WFP = k0 * np.exp(-2.0 * delta_U(W_FP) / sigma**2)
    T_esc = 1.0 / max(k_WFP, 1e-12)

    # Full OU variance
    stationary_var = eps * sigma**2 * tau_v / (2.0 * b)
    decay          = 1.0 - np.exp(-2.0 * b * eps * T_esc)
    var_w          = stationary_var * decay

    return float(W_FP), float(np.sqrt(var_w))


def gaussian_pdf(w: np.ndarray, mean: float, std: float) -> np.ndarray:
    """Evaluate Gaussian density at w."""
    return np.exp(-0.5 * ((w - mean) / std)**2) / (std * np.sqrt(2 * np.pi))




def mean_w_correction(eps: float, v_sad_offset: float = 0.02) -> float:
    """
    Second-order correction to mean w_escape from the escape path integral.

    δw = ε · ∫_{V_FP+δ}^{V_SADDLE-δ} (v + a - bW_FP) / (v - v³/3 - W_FP + I) dv

    The integrand has a removable singularity at V_FP (0/0, fixed-point condition)
    and a simple pole at V_SADDLE (velocity → 0 on nullcline).  We avoid both
    endpoints by a small offset v_sad_offset; the excluded region contributes
    negligibly since noise dominates near the saddle.
    """
    from scipy.integrate import quad

    def integrand(v):
        num = v + a - b * W_FP
        den = v - v**3 / 3.0 - W_FP + I
        return num / den if abs(den) > 1e-10 else 0.0

    v_lo = V_FP     + 0.001
    v_hi = V_SADDLE - v_sad_offset
    gamma, _ = quad(integrand, v_lo, v_hi, limit=200)
    return eps * V_SADDLE * gamma



def plot_barrier_decay(ax: plt.Axes) -> None:
    w_fold = w_fold_left(I)
    w_grid = np.linspace(W_FP, w_fold - 1e-3, 500)
    dU_arr = np.array([delta_U(w) for w in w_grid])
    A_arr  = np.array([kramers_prefactor(w) or 0.0 for w in w_grid])

    ax.plot(w_grid, dU_arr, color="#1a9641", lw=2.2, label=r"$\Delta U(w)$")
    ax2 = ax.twinx()
    ax2.plot(w_grid, A_arr, color="#fdae61", lw=1.8, ls="--", label="$A(w)$ prefactor")
    ax2.set_ylabel("$A(w)$  (Kramers prefactor)", fontsize=9, color="#fdae61")
    ax2.tick_params(axis="y", labelcolor="#fdae61")

    ax.axvline(W_FP,   color="gray",  lw=1.0, ls="--", alpha=0.6)
    ax.axvline(w_fold, color="black", lw=1.0, ls=":",  alpha=0.6)
    ax.set_xlabel("$w$", fontsize=10)
    ax.set_ylabel(r"$\Delta U(w)$  (barrier height)", fontsize=10)
    ax.set_title(
        r"Barrier decay: $\Delta U(w)$ from $W_{FP}$ to $w_{fold}$"
        f"\n$\\Delta U(W_{{FP}}) = {delta_U(W_FP):.5f}$,"
        f"  $w_{{fold}} = {w_fold:.3f}$",
        fontsize=9)
    ax.legend(fontsize=8, loc="upper right")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("Estimating v autocorrelation time τ_v near fixed point...")
# Use the middle pair's (sigma, eps) as representative — tau_v is not very
# sensitive to the exact parameters as long as the system stays near W_FP
_sig_ref, _eps_ref = TEST_PAIRS[len(TEST_PAIRS)//2]
TAU_V = estimate_tau_v(_sig_ref, _eps_ref)
print(f"  τ_v = {TAU_V:.4f}  (eps={_eps_ref}, sigma_probe={min(_sig_ref*0.3,0.05):.3f})\n")


print(f"  I={I}, a={a}, b={b}")
print(f"  W_FP={W_FP:.4f},  w_fold={w_fold_left(I):.4f},  ΔU(W_FP)={delta_U(W_FP):.5f}")
print(f"  N_traj={args.n_traj}  per pair,  seed={args.seed}")
print(f"  Pairs: {TEST_PAIRS}")
print(f"{'='*64}\n")

# --- Figure 1: barrier decay (standalone diagnostic) ---
fig_bd, ax_bd = plt.subplots(figsize=(6, 4))
plot_barrier_decay(ax_bd)
fig_bd.tight_layout()
fig_bd.savefig("figures/barrier_decay.png", dpi=150)
plt.close(fig_bd)
print("Saved: figures/barrier_decay.png")

# --- Figure 2: validation grid ---
n_pairs = len(TEST_PAIRS)
n_cols  = 3
n_rows  = int(np.ceil(n_pairs / n_cols))
fig, axes = plt.subplots(n_rows, n_cols, figsize=(5.5 * n_cols, 4.5 * n_rows))
axes_flat = axes.flatten() if n_pairs > 1 else [axes]

# Summary table header
print(f"{'sigma':>7}  {'eps':>6}  {'norm':>7}  {'w*_pred':>9}  "
      f"{'<w>_pred':>9}  {'<w>_sim':>9}  {'KS_stat':>8}  fired")
print("-" * 72)

ks_results = []   # list of dicts

# Summary table header
print(f"{'sigma':>7}  {'eps':>6}  {'std_sim':>8}  {'std_OU':>8}  "
      f"{'KS_Pois':>9}  {'KS_Gauss':>9}  {'KS_2nd':>8}  fired")
print("-" * 80)

for idx, (sigma, eps) in enumerate(TEST_PAIRS):
    print(f"  Running sigma={sigma}, eps={eps} ...", flush=True)

    model = commitment_model(sigma, eps)
    sim   = run_validation(sigma, eps, n_traj=args.n_traj, seed=args.seed + idx)

    mean_w_sim = float(np.mean(sim["w_escapes"])) if len(sim["w_escapes"]) > 0 else np.nan
    std_w_sim  = float(np.std(sim["w_escapes"]))  if len(sim["w_escapes"]) > 1 else np.nan
    _, g_std   = gaussian_w_escape(sigma, eps, tau_v=TAU_V)

# Hide unused axes
for idx in range(n_pairs, len(axes_flat)):
    axes_flat[idx].set_visible(False)

fig.suptitle(
    "Non-homogeneous Poisson commitment model vs. empirical $w_{\\rm escape}$\n"
    "Red curve: predicted $p(w)$ — Blue histogram: simulated escapes",
    fontsize=12, y=1.01)
fig.tight_layout()
fig.savefig("figures/commitment_validation.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"\nSaved: figures/commitment_validation.png")

# --- Figure 3: ΔU(w) contribution — how much of barrier remains at w* ---
fig3, ax3 = plt.subplots(figsize=(6, 4))
w_fold = w_fold_left(I)
w_grid = np.linspace(W_FP, w_fold - 1e-3, 500)
dU_arr = np.array([delta_U(w) for w in w_grid])

ax3.plot(w_grid, dU_arr, color="#1a9641", lw=2.0,
         label=r"$\Delta U(w)$")
ax3.axhline(delta_U(W_FP), color="gray", lw=1.0, ls="--",
            label=f"$\\Delta U(W_{{FP}}) = {delta_U(W_FP):.5f}$  (naive Kramers)")

colors = plt.cm.RdBu(np.linspace(0.1, 0.9, len(TEST_PAIRS)))
for (sigma, eps), ks, col in zip(TEST_PAIRS, ks_results, colors):
    model = commitment_model(sigma, eps)
    dU_star = delta_U(model["w_star"])
    ax3.axvline(model["w_star"], color=col, lw=1.5, ls=":",
                label=f"$w^*(\\sigma={sigma},\\varepsilon={eps})$  "
                      f"$\\Delta U={dU_star:.4f}$")

ax3.axvline(W_FP,   color="gray",  lw=1.0, ls="--", alpha=0.5)
ax3.axvline(w_fold, color="black", lw=1.0, ls=":",  alpha=0.5)
ax3.set_xlabel("$w$", fontsize=11)
ax3.set_ylabel(r"$\Delta U(w)$", fontsize=11)
ax3.set_title(
    "Effective barrier at predicted commitment point $w^*$\n"
    "Naive Kramers uses $\\Delta U(W_{FP})$ — model uses $\\Delta U(w^*)$",
    fontsize=10)
ax3.legend(fontsize=7, loc="upper right")
fig3.tight_layout()
fig3.savefig("figures/effective_barrier.png", dpi=150)
plt.close(fig3)
print("Saved: figures/effective_barrier.png\n")

# --- Summary ---
print(f"\n{'='*78}")
print("  Model comparison: KS statistics  (lower = better fit to w_escape)")
print(f"  {'sigma':>7}  {'eps':>6}  {'KS_Poisson':>11}  {'KS_Gauss':>10}  {'KS_2nd':>8}  winner")
print(f"  {'-'*66}")
for (sigma, eps), ks in zip(TEST_PAIRS, ks_results):
    kp  = ks["ks_poisson"]
    kg  = ks["ks_gaussian"]
    k2  = ks.get("ks_2nd", np.nan)
    best = min((v for v in [kp, kg, k2] if not np.isnan(v)), default=np.nan)
    winner = ("Poisson"  if best == kp else
              "Gauss"    if best == kg else
              "Gauss+2nd" if best == k2 else "—")
    print(f"  {sigma:7.3f}  {eps:6.3f}  {kp:11.4f}  {kg:10.4f}  {k2:8.4f}  {winner}")
print(f"{'='*78}\n")
