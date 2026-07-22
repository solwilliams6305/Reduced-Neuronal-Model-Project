"""
folded_cycle_normal_form_map.py
-------------------------------
Channel-A (amplitude-escape) hit-location regime map for the STOCHASTIC
folded LIMIT CYCLE, on the Jelbart-Kuehn-Kuntz (2024) rescaling-chart (K2)
normal form.

Deterministic scaffold (JKK 2024, arXiv:2208.01361, eq. (14) -> K2 chart,
frozen-phase case alpha = 2). With degenerate noise sigma dW on the fast
RADIAL variable r, pushed through the JKK blow-up
    (r, theta, y, eps) = (rho2 r2, theta2, rho2^2 y2, rho2),   rho2 = eps,
    eps2 = eps^3   =>   eta := sigma / sqrt(eps2),
the inner-chart SDE is

    dr2 = ( b(th) r2^2 - a(th) y2 ) dt2 + eta dB_t2,     (fast radius)
    dy2 = - c(th) dt2,                                   (parameter drift)

with theta frozen (alpha = 2), a, b, c > 0 the JKK normal-form coefficients
at the escape phase th. This is the Krupa-Szmolyan / canard Riccati fold
with theta-dependent coefficients.

PREDICTED inner scales (this work, FW accumulated-variance argument):
    y2_*   = (c^2 / (a b))^{1/3}        (drift window; analog of W_* = lam^{2/3})
    r2_*   = (a c / b^2)^{1/3}          (escape amplitude; analog of V_* = lam^{1/3})
    T_win  = (a b c)^{-1/3}             (passage time in t2)
    eta_*  = C_q * sqrt(a c / b)        (Channel-A escape scale; analog of C_q lam^{1/2})

=> physical critical noise   sigma_*^A(th) = C_q * sqrt(eps2) * sqrt( a(th) c(th) / b(th) ).

Scaled (collapse) variables, exactly as canard_normal_form_map.py:
    Theta = eta / sqrt(a c / b)        (eta_* = sqrt(a c / b),  C_q absorbed)
    R_hit = y2_hit / y2_*

TEST.  If sqrt(a c / b) is the correct geometric noise scale, the distribution
of R_hit at fixed Theta must be INDEPENDENT of the geometry (a, b, c) -- the
collapse criterion.  Under the exact rescaling to canonical form
    dR = (R^2 - Y) dT + Theta dB,  dY = -dT,  cross R = 1,
every (a, b, c) maps onto the SAME canard normal form, so the median R_hit
must cross zero at the SAME Theta_crit ~ 2.8 found for the fold/canard
(C_q inheritance, cross-chapter consistency check).

FALSIFICATION CONTROL.  We rerun with the WRONG normalization eta = Theta *
sqrt(a / b) (the c-factor dropped).  If c genuinely belongs in eta_*, the
collapse must break across geometries with different c.
"""
from __future__ import annotations

import os
import time

import numpy as np


# ---------------------------------------------------------------------------
# Inner-chart integrator (vectorised Euler-Maruyama), frozen-phase alpha = 2
# ---------------------------------------------------------------------------

def simulate_inner(
    a: float,
    b: float,
    c: float,
    eta: float,
    n_traj: int = 500,
    y0_over_star: float = 5.0,
    T_canon_max: float = 40.0,
    dT_canon: float = 1e-3,
    Rcross_canon: float = 1.0,
    rng: np.random.Generator | None = None,
) -> tuple[np.ndarray, float]:
    """
    Integrate  dr2 = (b r2^2 - a y2) dt2 + eta dB,  dy2 = -c dt2
    from the attracting branch r2 = -sqrt((a/b) y2) at y2 = y0_over_star * y2_*.

    Resolution is set in CANONICAL time T2 = (abc)^{1/3} t2 so every geometry
    is integrated at identical canonical step/length:
        dt2 = dT_canon * tau,  T2_max -> t2_max = T_canon_max * tau,
        tau = (abc)^{-1/3}.
    Crossing threshold r2_cross = Rcross_canon * r2_*  (geometry-scaled "spike").

    Returns (y2_hit, y2_star); y2_hit is NaN where no crossing occurred.
    """
    if rng is None:
        rng = np.random.default_rng()

    tau = (a * b * c) ** (-1.0 / 3.0)
    dt2 = dT_canon * tau
    t2_max = T_canon_max * tau

    y2_star = (c * c / (a * b)) ** (1.0 / 3.0)
    r2_star = (a * c / (b * b)) ** (1.0 / 3.0)
    r2_cross = Rcross_canon * r2_star

    y2_0 = y0_over_star * y2_star
    r2_0 = -np.sqrt((a / b) * y2_0)

    r2 = np.full(n_traj, r2_0, dtype=float)
    y2 = np.full(n_traj, y2_0, dtype=float)
    hit = np.zeros(n_traj, dtype=bool)
    y2_hit = np.full(n_traj, np.nan)

    sqrt_dt = np.sqrt(dt2)
    n_steps = int(t2_max / dt2)

    for _ in range(n_steps):
        active = ~hit
        if not active.any():
            break
        noise = rng.standard_normal(n_traj)
        r2_new = r2.copy()
        y2_new = y2.copy()
        r2_new[active] = (
            r2[active]
            + (b * r2[active] ** 2 - a * y2[active]) * dt2
            + eta * sqrt_dt * noise[active]
        )
        y2_new[active] = y2[active] - c * dt2

        new_hit = active & (r2 < r2_cross) & (r2_new >= r2_cross)
        if new_hit.any():
            frac = (r2_cross - r2[new_hit]) / (r2_new[new_hit] - r2[new_hit])
            y2_hit[new_hit] = y2[new_hit] + frac * (y2_new[new_hit] - y2[new_hit])
            hit |= new_hit

        r2, y2 = r2_new, y2_new

    return y2_hit, y2_star


# ---------------------------------------------------------------------------
# Sweep over (geometry, Theta), correct vs wrong normalization
# ---------------------------------------------------------------------------

GEOMS = [
    (1.0, 1.0, 1.0),
    (2.0, 1.0, 1.0),
    (1.0, 2.0, 1.0),
    (1.0, 1.0, 2.0),
    (0.5, 1.5, 1.0),
    (1.7, 0.8, 1.3),
]


def run_sweep(
    geoms=GEOMS,
    theta_vals=np.logspace(-0.5, 1.0, 10),
    n_traj=500,
    seed=42,
    normalization="correct",
    verbose=True,
):
    """normalization='correct' -> eta = Theta sqrt(ac/b);
       'wrong' -> eta = Theta sqrt(a/b)  (c dropped, falsification control)."""
    rng_master = np.random.default_rng(seed)
    grid = np.full((len(geoms), len(theta_vals), n_traj), np.nan)
    y2_star = np.zeros(len(geoms))

    t0 = time.time()
    for j, (a, b, c) in enumerate(geoms):
        if normalization == "correct":
            eta_star = np.sqrt(a * c / b)
        elif normalization == "wrong":
            eta_star = np.sqrt(a / b)
        else:
            raise ValueError(normalization)
        for k, theta in enumerate(theta_vals):
            eta = theta * eta_star
            rng = np.random.default_rng(rng_master.integers(0, 2**63 - 1))
            y2_hit, ys = simulate_inner(a, b, c, eta, n_traj=n_traj, rng=rng)
            grid[j, k, :] = y2_hit
            y2_star[j] = ys
        if verbose:
            fired = np.mean(~np.isnan(grid[j]))
            print(f"  geom (a,b,c)=({a:.2f},{b:.2f},{c:.2f})  "
                  f"eta_*={eta_star:.3f}  y2_*={y2_star[j]:.3f}  "
                  f"fired={fired:.2f}  [{normalization}]")
    if verbose:
        print(f"  wall time: {time.time() - t0:.1f}s")

    return dict(
        geoms=np.asarray(geoms),
        theta_vals=np.asarray(theta_vals),
        y2_hit=grid,
        y2_star=y2_star,
        normalization=normalization,
    )


# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------

def median_R_table(data):
    """median R_hit at each (geometry, Theta). Returns array [n_geom, n_theta]."""
    geoms = data["geoms"]
    theta_vals = data["theta_vals"]
    y2_star = data["y2_star"]
    Y = data["y2_hit"]
    tab = np.full((len(geoms), len(theta_vals)), np.nan)
    for j in range(len(geoms)):
        for k in range(len(theta_vals)):
            R = Y[j, k, :] / y2_star[j]
            R = R[~np.isnan(R)]
            if R.size:
                tab[j, k] = np.median(R)
    return tab


def pooled_summary(data):
    geoms = data["geoms"]
    theta_vals = data["theta_vals"]
    y2_star = data["y2_star"]
    Y = data["y2_hit"]

    bins_edges = np.array([-np.inf, 0.1, 0.5, 2.0, 4.0, np.inf])
    bin_labels = ["R<0.1", "0.1-0.5", "0.5-2", "2-4", "R>4"]

    per = {}
    for k, theta in enumerate(theta_vals):
        R_all = []
        for j in range(len(geoms)):
            R = Y[j, k, :] / y2_star[j]
            R_all.append(R[~np.isnan(R)])
        R_all = np.concatenate(R_all) if R_all else np.array([])
        if R_all.size:
            counts, _ = np.histogram(R_all, bins=bins_edges)
            per[float(theta)] = dict(
                n_fired=int(R_all.size),
                median=float(np.median(R_all)),
                q25=float(np.quantile(R_all, 0.25)),
                q75=float(np.quantile(R_all, 0.75)),
                bin_probs=counts / counts.sum(),
            )
        else:
            per[float(theta)] = dict(n_fired=0, median=np.nan, q25=np.nan,
                                     q75=np.nan, bin_probs=np.zeros(len(bin_labels)))
    return dict(per_theta=per, bin_labels=bin_labels)


def theta_crit(theta_vals, median_pooled):
    """Linear-interpolate (in log Theta) the Theta at which median R crosses 0."""
    lt = np.log10(theta_vals)
    m = median_pooled
    for i in range(len(m) - 1):
        if np.isfinite(m[i]) and np.isfinite(m[i + 1]) and m[i] < 0 <= m[i + 1]:
            f = (0.0 - m[i]) / (m[i + 1] - m[i])
            return 10 ** (lt[i] + f * (lt[i + 1] - lt[i]))
    return np.nan


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def make_figure(data, summary, tab, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    theta_vals = data["theta_vals"]
    geoms = data["geoms"]
    per = summary["per_theta"]
    bin_labels = summary["bin_labels"]

    medians = np.array([per[float(t)]["median"] for t in theta_vals])
    q25s = np.array([per[float(t)]["q25"] for t in theta_vals])
    q75s = np.array([per[float(t)]["q75"] for t in theta_vals])
    tc = theta_crit(theta_vals, medians)

    fig, axes = plt.subplots(2, 1, figsize=(7.5, 8.4))

    ax = axes[0]
    ax.fill_between(theta_vals, q25s, q75s, alpha=0.20, color="C0",
                    label="pooled IQR")
    for j, (a, b, c) in enumerate(geoms):
        ax.plot(theta_vals, tab[j], "-", lw=1.0, alpha=0.7,
                label=f"(a,b,c)=({a:g},{b:g},{c:g})")
    ax.plot(theta_vals, medians, "ko-", lw=2.2, ms=5, label="pooled median",
            zorder=5)
    ax.axhline(0.0, color="k", lw=0.7, ls="--", alpha=0.5)
    if np.isfinite(tc):
        ax.axvline(tc, color="r", lw=1.2, ls=":",
                   label=fr"$\Theta_{{\rm crit}}\approx{tc:.2f}$")
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=0.5)
    ax.set_xlabel(r"$\Theta = \eta\,/\,\sqrt{a c / b}$")
    ax.set_ylabel(r"$R_{\rm hit}=y_{2,\rm hit}/y_{2,*}$,  $y_{2,*}=(c^2/ab)^{1/3}$")
    ax.set_title("Folded limit cycle, Channel-A regime map (JKK K2 normal form)\n"
                 "collapse across geometry $(a,b,c)$ at fixed $\\Theta$")
    ax.legend(loc="upper right", frameon=False, fontsize=7.5, ncol=2)
    ax.grid(True, which="both", alpha=0.25)

    ax = axes[1]
    bottom = np.zeros_like(theta_vals, dtype=float)
    colors = ["#2c7fb8", "#7fcdbb", "#cccccc", "#fdae6b", "#d94701"]
    for bnum, label in enumerate(bin_labels):
        probs = np.array([per[float(t)]["bin_probs"][bnum] for t in theta_vals])
        ax.bar(theta_vals, probs, bottom=bottom, width=theta_vals * 0.16,
               color=colors[bnum], label=label, edgecolor="white", lw=0.5)
        bottom += probs
    ax.set_xscale("log")
    ax.set_xlabel(r"$\Theta$")
    ax.set_ylabel(r"$P(R_{\rm hit}\ \mathrm{in\ bin})$")
    ax.set_title("Stacked $R_{\\rm hit}$-bin probability vs $\\Theta$ (pooled)")
    ax.set_ylim(0, 1)
    ax.legend(ncol=5, loc="lower center", bbox_to_anchor=(0.5, -0.30),
              frameon=False, fontsize=8)
    ax.grid(True, axis="y", alpha=0.25)

    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return tc


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)          # project root (regime-tests/ -> root)
    out_data = os.path.join(root, "data")
    out_fig = os.path.join(root, "figures")
    os.makedirs(out_data, exist_ok=True)
    os.makedirs(out_fig, exist_ok=True)

    theta_vals = np.logspace(-0.5, 1.0, 10)

    print("\n=== Folded limit cycle: Channel-A normal-form regime map ===\n")
    print("Tests  eta_* = C_q sqrt(a c / b)  on the JKK K2 inner SDE.\n")

    print("CORRECT normalization eta = Theta sqrt(ac/b):")
    data = run_sweep(theta_vals=theta_vals, normalization="correct", seed=42)
    print("\nWRONG normalization eta = Theta sqrt(a/b)  (c dropped; control):")
    data_w = run_sweep(theta_vals=theta_vals, normalization="wrong", seed=43)

    tab = median_R_table(data)
    tab_w = median_R_table(data_w)
    summ = pooled_summary(data)

    np.savez(
        os.path.join(out_data, "folded_cycle_normal_form_hit_location.npz"),
        geoms=data["geoms"], theta_vals=theta_vals,
        y2_hit=data["y2_hit"], y2_star=data["y2_star"],
        y2_hit_wrong=data_w["y2_hit"], y2_star_wrong=data_w["y2_star"],
    )
    print("\n  data saved -> data/folded_cycle_normal_form_hit_location.npz")

    fig_path = os.path.join(out_fig, "folded_cycle_normal_form_regime_map.png")
    tc = make_figure(data, summ, tab, fig_path)
    print(f"  figure saved -> figures/folded_cycle_normal_form_regime_map.png")

    per = summ["per_theta"]
    print("\nPooled Theta sweep (correct normalization):")
    print(f"  {'Theta':>8s}  {'n_fired':>8s}  {'median R':>10s}  {'IQR':>20s}")
    for t in theta_vals:
        d = per[float(t)]
        print(f"  {t:8.3f}  {d['n_fired']:8d}  {d['median']:+10.3f}  "
              f"[{d['q25']:+6.3f}, {d['q75']:+6.3f}]")

    print(f"\n  Theta_crit (median R crosses 0) = {tc:.2f}   "
          f"(canard/fold C_q ~ 2.8 expected)")

    # Collapse quality: cross-geometry spread of median R at each Theta
    print("\nCollapse quality  (std of median R across the 6 geometries):")
    print(f"  {'Theta':>8s}  {'correct sqrt(ac/b)':>20s}  {'wrong sqrt(a/b)':>18s}")
    sc_corr, sc_wrong = [], []
    for k, t in enumerate(theta_vals):
        s_c = np.nanstd(tab[:, k])
        s_w = np.nanstd(tab_w[:, k])
        sc_corr.append(s_c)
        sc_wrong.append(s_w)
        print(f"  {t:8.3f}  {s_c:20.3f}  {s_w:18.3f}")
    print(f"\n  mean cross-geometry std:  correct={np.nanmean(sc_corr):.3f}   "
          f"wrong={np.nanmean(sc_wrong):.3f}")
    print("  (correct << wrong  =>  the c-factor in eta_* = sqrt(ac/b) is required)\n")

    print("Collapse table  median R at each (geometry, Theta), correct norm:")
    hdr = "  ".join([f"({a:g},{b:g},{c:g})" for (a, b, c) in data["geoms"]])
    print(f"  {'Theta':>8s}  " + hdr)
    for k, t in enumerate(theta_vals):
        row = "  ".join([f"{tab[j, k]:+7.3f}" for j in range(len(data['geoms']))])
        print(f"  {t:8.3f}  " + row)


if __name__ == "__main__":
    main()
