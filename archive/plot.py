"""
plot.py
-------
Matplotlib visualisations for the FHN project.

plot_phase_portrait    — nullclines, fixed point, sample trajectories
plot_regime_map        — heatmap of log_ratio + failure boundary
plot_mfpt_comparison   — full MFPT vs drift time across sigma at fixed eps
plot_boundary_scaling  — log-log fit vs BG prediction
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.optimize import brentq

from simulate import _v_left_newton, w_fold_left

COLORS = {
    "full":        "#2166ac",
    "reduced":     "#d6604d",
    "boundary":    "#1a9641",
    "bg_curve":    "#fdae61",
    "nullcline_v": "#4d4d4d",
    "nullcline_w": "#969696",
}


def _fig(w=7, h=5, **kw):
    return plt.subplots(figsize=(w, h), **kw)


# ---------------------------------------------------------------------------
# Phase portrait
# ---------------------------------------------------------------------------

def plot_phase_portrait(
    result: dict,
    ax=None,
    show_trajectories: int = 5,
    title: str = "",
) -> plt.Figure:
    p    = result["params"]
    I, a, b = p["I"], p["a"], p["b"]

    if ax is None:
        fig, ax = _fig(6, 5)
    else:
        fig = ax.get_figure()

    v_range  = np.linspace(-2.5, 2.5, 500)
    w_v_null = v_range - v_range**3 / 3 + I
    w_range  = np.linspace(-1.5, 1.5, 200)
    v_w_null = b * w_range - a

    ax.plot(v_range, w_v_null, color=COLORS["nullcline_v"], lw=1.8, label="$v$-nullcline")
    ax.plot(v_w_null, w_range, color=COLORS["nullcline_w"], lw=1.8, ls="--", label="$w$-nullcline")

    def fp_eq(v):
        return v - v**3/3 - (v + a)/b + I

    try:
        v_fp = brentq(fp_eq, -2.5, -1.0)
        w_fp = (v_fp + a) / b
        ax.plot(v_fp, w_fp, "ko", ms=7, zorder=5, label=f"Fixed point ({v_fp:.2f}, {w_fp:.2f})")
    except Exception:
        pass

    for v_fold in [-1.0, 1.0]:
        wf = v_fold - v_fold**3/3 + I
        ax.axvline(v_fold, color="gray", lw=0.8, ls=":", alpha=0.5)
        ax.plot(v_fold, wf, "r^", ms=8, zorder=5)

    n = min(show_trajectories, result["v_sample"].shape[1])
    for k in range(n):
        ax.plot(result["v_sample"][:, k], result["w_sample"][:, k],
                color=COLORS["full"], alpha=0.4, lw=0.8)

    ax.set_xlabel("$v$ (voltage)", fontsize=12)
    ax.set_ylabel("$w$ (recovery)", fontsize=12)
    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim(-1.5, 1.5)
    ax.legend(fontsize=9, loc="upper left")
    ax.set_title(title or
        f"Phase portrait  ($\\sigma$={p['sigma']}, $\\varepsilon$={p['eps']})", fontsize=11)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Regime map  (log-ratio heatmap)
# ---------------------------------------------------------------------------

def plot_regime_map(
    sweep_result: dict,
    delta_ratio:    float = 2.0,
    boundary_eps:   np.ndarray | None = None,
    boundary_sigma: np.ndarray | None = None,
    bg_C:           float | None = None,
    crossover_delta_U: float | None = None,
    save_path:      str | None = None,
) -> plt.Figure:
    """
    Heatmap of log(T_drift / MFPT_full) in (eps, sigma) space.
    Positive (red) = noise shortens MFPT = reduction FAILS.
    Near zero (white) = good reduction.
    """
    sigma_vals = sweep_result["sigma_vals"]
    eps_vals   = sweep_result["eps_vals"]
    log_rat    = sweep_result["log_ratio"]   # (n_sigma, n_eps)

    fig, ax = _fig(7.5, 5.5)

    vmax = np.nanpercentile(np.abs(log_rat), 95)
    vmax = max(vmax, 1.0)
    norm = TwoSlopeNorm(vmin=-vmax/2, vcenter=0, vmax=vmax)

    im = ax.pcolormesh(eps_vals, sigma_vals, log_rat,
                       cmap="RdBu_r", norm=norm, shading="auto")
    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label(
        r"$\log\!\left(T_{\rm drift}\,/\,\mathrm{MFPT}_{\rm full}\right)$"
        "\n(+ve = noise shortens MFPT, reduction fails)",
        fontsize=9)

    # Zero contour
    try:
        ax.contour(eps_vals, sigma_vals, log_rat, levels=[0],
                   colors=["k"], linewidths=1.0, linestyles="--")
    except Exception:
        pass

    # Failure boundary
    log_thresh = np.log(delta_ratio)
    try:
        ax.contour(eps_vals, sigma_vals, log_rat, levels=[log_thresh],
                   colors=[COLORS["boundary"]], linewidths=2.5)
    except Exception:
        pass

    if boundary_eps is not None and len(boundary_eps) > 0:
        ax.plot(boundary_eps, boundary_sigma, "o",
                color=COLORS["boundary"], ms=6, zorder=5,
                label=f"Failure boundary ($\\times${delta_ratio:.0f} MFPT error)")

    if bg_C is not None:
        eps_plot  = np.linspace(eps_vals.min(), eps_vals.max(), 300)
        sigma_bg  = bg_C * eps_plot**1.5
        mask      = sigma_bg <= sigma_vals.max()
        ax.plot(eps_plot[mask], sigma_bg[mask], "--",
                color=COLORS["bg_curve"], lw=2.5,
                label=f"$\\sigma = {bg_C:.2f}\\,\\varepsilon^{{3/2}}$  (BG)")

    if crossover_delta_U is not None:
        from sweep import kramers_crossover_curve
        eps_plot   = np.linspace(eps_vals.min(), eps_vals.max(), 300)
        sigma_cross = kramers_crossover_curve(eps_plot, crossover_delta_U)
        mask = (sigma_cross >= sigma_vals.min()) & (sigma_cross <= sigma_vals.max())
        ax.plot(eps_plot[mask], sigma_cross[mask], "-.",
                color="yellow", lw=2.0, alpha=0.85,
                label=f"$\\sigma_{{cross}} = \\sqrt{{\\Delta U / \\log(1/\\varepsilon)}}$"
                      f"  ($\\Delta U={crossover_delta_U:.4f}$)")

    ax.set_xlabel("$\\varepsilon$  (timescale separation)", fontsize=12)
    ax.set_ylabel("$\\sigma$  (noise amplitude)", fontsize=12)
    ax.set_title("Regime map: slow-manifold reduction reliability", fontsize=12)
    if boundary_eps is not None or bg_C is not None:
        ax.legend(fontsize=9, loc="upper left")

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150)
    return fig


# ---------------------------------------------------------------------------
# MFPT comparison at fixed eps
# ---------------------------------------------------------------------------

def plot_mfpt_comparison(
    sweep_result: dict,
    eps_fixed: float | None = None,
    save_path:  str | None = None,
) -> plt.Figure:
    sigma_vals = sweep_result["sigma_vals"]
    eps_vals   = sweep_result["eps_vals"]
    drift_times = sweep_result.get("drift_times",
                  sweep_result["mfpt_reduced"][0, :])

    if eps_fixed is None:
        eps_fixed = eps_vals[len(eps_vals) // 2]

    j      = np.argmin(np.abs(eps_vals - eps_fixed))
    mfpt_f = sweep_result["mfpt_full"][:, j]
    td     = drift_times[j]

    fig, ax = _fig(6, 4)
    ax.semilogy(sigma_vals, mfpt_f, "o-", color=COLORS["full"],
                label="MFPT$_{\\rm full}$ (simulated)", ms=5)
    ax.axhline(td, color=COLORS["reduced"], ls="--", lw=2,
               label=f"$T_{{\\rm drift}}$ (reduced,  $\\varepsilon={eps_fixed:.3f}$) = {td:.1f}")

    ax.set_xlabel("$\\sigma$  (noise amplitude)", fontsize=12)
    ax.set_ylabel("MFPT  (log scale)", fontsize=12)
    ax.set_title(f"MFPT comparison at $\\varepsilon = {eps_fixed:.3f}$", fontsize=11)
    ax.legend(fontsize=10)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150)
    return fig


# ---------------------------------------------------------------------------
# Boundary power-law scaling
# ---------------------------------------------------------------------------

def plot_boundary_scaling(
    boundary_eps: np.ndarray,
    boundary_sigma: np.ndarray,
    alpha_fit: float,
    C_fit: float,
    crossover_delta_U: float | None = None,
    save_path: str | None = None,
) -> plt.Figure:
    fig, ax = _fig(6, 4)

    ax.loglog(boundary_eps, boundary_sigma, "o",
              color=COLORS["boundary"], ms=8, zorder=5, label="Empirical boundary")

    eps_plot = np.linspace(boundary_eps.min() * 0.7,
                           boundary_eps.max() * 1.5, 300)

    if not np.isnan(alpha_fit):
        ax.loglog(eps_plot, C_fit * eps_plot**alpha_fit, "-",
                  color=COLORS["full"], lw=2,
                  label=f"Fit: $\\sigma = {C_fit:.2f}\\,\\varepsilon^{{{alpha_fit:.2f}}}$")

    ax.loglog(eps_plot, eps_plot**1.5, "--",
              color=COLORS["bg_curve"], lw=2,
              label="BG: $\\sigma \\sim \\varepsilon^{3/2}$")

    if crossover_delta_U is not None:
        from sweep import kramers_crossover_curve
        sigma_cross = kramers_crossover_curve(eps_plot, crossover_delta_U)
        ax.loglog(eps_plot, sigma_cross, "-.",
                  color=COLORS["reduced"], lw=2,
                  label=f"Crossover: $\\sqrt{{\\Delta U / \\log(1/\\varepsilon)}}$"
                        f"  ($\\Delta U={crossover_delta_U:.4f}$)")

    ax.set_xlabel("$\\varepsilon$", fontsize=12)
    ax.set_ylabel("$\\sigma^*$ (failure boundary)", fontsize=12)
    ax.set_title("Failure boundary scaling", fontsize=11)
    ax.legend(fontsize=10)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150)
    return fig

# ---------------------------------------------------------------------------
# Diagnostic 1: local exponent vs eps
# ---------------------------------------------------------------------------

def plot_local_exponent(
    eps_centres:  np.ndarray,
    local_alphas: np.ndarray,
    alpha_global: float,
    save_path: str | None = None,
) -> plt.Figure:
    """
    Plot the sliding-window power-law exponent as a function of eps.

    If 0.43 is a finite-eps correction, alpha should drift toward 0.5 as
    eps -> 0 (left side of plot).  If it stays flat, the gap is structural.

    Reference lines at alpha = 0.5 (fold-layer) and alpha = 1.5 (BG).
    """
    fig, ax = _fig(6, 4)

    ax.semilogx(eps_centres, local_alphas, "o-",
                color=COLORS["boundary"], ms=7, lw=2,
                label="Local exponent (sliding window)")

    ax.axhline(0.5,  color=COLORS["full"],     lw=1.5, ls="--",
               label=r"$\alpha = 0.5$  (fold-layer $\sqrt{\varepsilon}$)")
    ax.axhline(1.5,  color=COLORS["bg_curve"], lw=1.5, ls="--",
               label=r"$\alpha = 1.5$  (Berglund–Gentz $\varepsilon^{3/2}$)")
    if not np.isnan(alpha_global):
        ax.axhline(alpha_global, color=COLORS["reduced"], lw=1.5, ls=":",
                   label=f"Global fit $\\alpha = {alpha_global:.2f}$")

    ax.set_xlabel(r"$\varepsilon$  (window centre, log scale)", fontsize=12)
    ax.set_ylabel(r"Local power-law exponent $\alpha$", fontsize=12)
    ax.set_title(r"Local exponent: does $\alpha \to 0.5$ as $\varepsilon \to 0$?",
                 fontsize=11)
    ax.legend(fontsize=9)
    ax.set_ylim(0, 2.0)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150)
    return fig


# ---------------------------------------------------------------------------
# Diagnostic 2: w_spike heatmap
# ---------------------------------------------------------------------------

def plot_w_spike_map(
    sweep_result: dict,
    boundary_eps:   np.ndarray | None = None,
    boundary_sigma: np.ndarray | None = None,
    use_escape: bool = True,
    save_path: str | None = None,
) -> plt.Figure:
    """
    Heatmap of mean w at escape initiation across (eps, sigma) space.

    use_escape=True  (default): uses w_escape_mean — w when v first crosses
                                V_SADDLE. This is the true escape initiation
                                point before right-branch drift contaminates w.
    use_escape=False           : uses w_spike_mean — w at v=1.0 threshold.

    Key reference values:
      W_FP    ≈ -0.695  (fixed point — fast noise escape at constant w)
      w_fold  ≈ -0.767  (left fold  — deterministic drift-to-fold mechanism)

    If w_escape ≈ W_FP:    pre-fold escape dominates (hazard / Kramers mechanism)
    If w_escape ≈ w_fold:  fold-layer crossing is the failure mode
    """
    from simulate import w_fold_left, W_FP

    sigma_vals = sweep_result["sigma_vals"]
    eps_vals   = sweep_result["eps_vals"]

    key = "w_escape_mean" if use_escape else "w_spike_mean"
    w_arr = sweep_result.get(key)
    label = "escape (saddle crossing)" if use_escape else "spike (v=1 threshold)"

    if w_arr is None:
        raise ValueError(f"sweep_result missing '{key}' — re-run sweep.")

    I      = sweep_result["params"]["I"]
    w_fold = w_fold_left(I)

    # Shift so 0 = fixed-point w; fold is at w_fold - W_FP (negative)
    # Colour by distance from W_FP: near 0 = fast escape, near (w_fold - W_FP) = fold mechanism
    w_rel = w_arr - W_FP   # 0 = escaped at fixed-point w, negative = escaped after slow drift

    fig, ax = _fig(7.5, 5.5)

    vext = max(np.nanpercentile(np.abs(w_rel), 95), 0.05)
    norm = TwoSlopeNorm(vmin=-vext, vcenter=0, vmax=vext)

    im = ax.pcolormesh(eps_vals, sigma_vals, w_rel,
                       cmap="RdBu_r", norm=norm, shading="auto")
    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label(
        rf"$\langle w_{{\rm {('escape' if use_escape else 'spike')}}}\rangle - W_{{FP}}$"
        f"\n(0 = fast escape, {w_fold - W_FP:.3f} = fold mechanism)",
        fontsize=9)

    # Mark the fold-mechanism level
    w_fold_line = w_fold - W_FP
    try:
        ax.contour(eps_vals, sigma_vals, w_rel, levels=[w_fold_line],
                   colors=["k"], linewidths=1.0, linestyles=":",
                   alpha=0.7)
    except Exception:
        pass

    if boundary_eps is not None and len(boundary_eps) > 0:
        ax.plot(boundary_eps, boundary_sigma, "o-",
                color=COLORS["boundary"], ms=5, lw=1.5, zorder=5,
                label="Failure boundary")
        ax.legend(fontsize=9, loc="upper left")

    ax.set_xlabel(r"$\varepsilon$  (timescale separation)", fontsize=12)
    ax.set_ylabel(r"$\sigma$  (noise amplitude)", fontsize=12)
    ax.set_title(
        f"Mean $w$ at {label}\n"
        rf"$W_{{FP}}={W_FP:.3f}$,  $w_{{\rm fold}}={w_fold:.3f}$  "
        r"(near 0 $\Rightarrow$ pre-fold escape)",
        fontsize=10)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150)
    return fig
