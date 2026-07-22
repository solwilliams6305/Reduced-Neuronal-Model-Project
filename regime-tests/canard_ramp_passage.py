"""
canard_ramp_passage.py
----------------------
Ramped canard passage experiment (Roadmap §5–7).

Dynamic complement to the autonomous canard chapter. Measures where
noise-driven trajectories fire (I_spike) when I(t) = I_0 + rho * t
ramps through the canard-explosion window near I_H1.

Observable:
    I_spike = I(t_spike),  t_spike = first time v >= +1 from below.
    Delay   = E[I_spike] - I_H1      (> 0 even at sigma = 0: BG delay)

Sweep:
    eps  in {0.02, 0.04, 0.08, 0.16}
    rho  in {eps^1.5, eps^1.0, eps^0.5}   (3 rates per eps)
    sigma in logspace(-3, -0.5, 12)
    N = 400 trajectories per (eps, rho, sigma)

Outputs:
    data/canard_ramp_passage.npz
    figures/canard_ramp_passage.png  (2x2 layout)
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _ROOT)
sys.path.insert(0, _HERE)

import _shim   # noqa: F401
from kernel import FHN2D

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
A_FHN    = 0.7
B_FHN    = 0.8
I_FOLD_L = (A_FHN - 1.0 + 2.0 * B_FHN / 3.0) / B_FHN   # ≈ 0.2917

V_SPIKE  = +1.0   # threshold that defines "a spike"


# ---------------------------------------------------------------------------
# Integrator
# ---------------------------------------------------------------------------

def simulate_ramp(
    eps:    float,
    I_0:    float,
    rho:    float,
    sigma:  float,
    I_H1:   float,
    v0:     float,
    w0:     float,
    n_traj: int = 400,
    dt:     float = 5e-3,
    rng:    np.random.Generator | None = None,
) -> np.ndarray:
    """
    Vectorised Euler-Maruyama for the ramped FHN SDE.

    I(t) = I_0 + rho * t

    Records I_spike = I(t_spike) at first upward crossing of v = V_SPIKE.
    Returns NaN for trajectories that reach I > I_H1 + 8*eps without spiking.

    Returns
    -------
    I_spike : shape (n_traj,)
    """
    if rng is None:
        rng = np.random.default_rng()

    I_stop = I_H1 + 8.0 * eps

    v = np.full(n_traj, v0, dtype=float)
    w = np.full(n_traj, w0, dtype=float)
    t = 0.0

    hit    = np.zeros(n_traj, dtype=bool)
    I_spike = np.full(n_traj, np.nan)

    sqrt_dt = np.sqrt(dt)
    # max steps: run until slowest trajectory would exceed I_stop
    t_max   = (I_stop - I_0) / rho if rho > 0 else 1e6
    n_steps = int(t_max / dt) + 1

    for _ in range(n_steps):
        if hit.all():
            break

        I_t = I_0 + rho * t

        # bail-out: any active trajectory whose I(t) exceeds I_stop gets NaN
        bail = (~hit) & (I_t >= I_stop)
        if bail.any():
            # mark as done but keep NaN — already initialised
            hit[bail] = True
            if hit.all():
                break

        active = ~hit
        if not active.any():
            break

        noise = rng.standard_normal(n_traj)

        v_new = np.copy(v)
        w_new = np.copy(w)
        v_new[active] = (v[active]
                         + (v[active] - v[active]**3 / 3.0 - w[active] + I_t) * dt
                         + sigma * sqrt_dt * noise[active])
        w_new[active] = (w[active]
                         + eps * (v[active] + A_FHN - B_FHN * w[active]) * dt)

        # spike detector: v crosses V_SPIKE from below
        spike = active & (v < V_SPIKE) & (v_new >= V_SPIKE)
        if spike.any():
            # linear interpolate t_cross, then I_spike = I_0 + rho*t_cross
            dv    = v_new[spike] - v[spike]
            frac  = np.where(np.abs(dv) > 1e-15,
                             (V_SPIKE - v[spike]) / dv,
                             0.5)
            frac  = np.clip(frac, 0.0, 1.0)
            t_cross = t + frac * dt
            I_spike[spike] = I_0 + rho * t_cross
            hit[spike] = True

        v, w = v_new, w_new
        t   += dt

    return I_spike


# ---------------------------------------------------------------------------
# Sweep
# ---------------------------------------------------------------------------

def run_sweep(
    eps_vals   = (0.02, 0.04, 0.08, 0.16),
    n_rho      = 3,          # rho exponents: 1.5, 1.0, 0.5
    sigma_vals = np.logspace(-3, -0.5, 12),
    n_traj     = 400,
    seed       = 42,
    verbose    = True,
):
    fhn_template = FHN2D()
    rng_master   = np.random.default_rng(seed)

    n_eps   = len(eps_vals)
    n_sigma = len(sigma_vals)

    # rho for each eps: shape (n_eps, n_rho)
    rho_exponents = [1.5, 1.0, 0.5]

    # storage
    I_spike_grid = np.full((n_eps, n_rho, n_sigma, n_traj), np.nan)
    rho_grid     = np.full((n_eps, n_rho), np.nan)
    I_H1_grid    = np.full(n_eps, np.nan)
    I_0_grid     = np.full(n_eps, np.nan)

    t0    = time.time()
    total = n_eps * n_rho * n_sigma
    done  = 0

    for i, eps in enumerate(eps_vals):
        dt = 2e-3 if eps <= 0.02 else 5e-3

        I_H1 = fhn_template.I_hopf_lower_at(eps)
        I_0  = max(I_FOLD_L, I_H1 - 4.0 * eps)
        I_H1_grid[i] = I_H1
        I_0_grid[i]  = I_0

        # initial condition: deterministic FP at I_0
        fhn_I0 = FHN2D(I=I_0)
        try:
            v0, w0 = fhn_I0.V_FP, fhn_I0.W_FP
        except Exception:
            v0, w0 = -1.5, (I_0 - 2.0 / 3.0)   # fallback: near left fold

        if verbose:
            print(f"\n  eps={eps:.2f}  I_H1={I_H1:.4f}  I_0={I_0:.4f}  "
                  f"v0={v0:.4f}  w0={w0:.4f}  dt={dt}")

        for jr, exp in enumerate(rho_exponents):
            rho = eps ** exp
            rho_grid[i, jr] = rho

            for k, sigma in enumerate(sigma_vals):
                rng = np.random.default_rng(rng_master.integers(0, 2**63 - 1))
                Is  = simulate_ramp(
                    eps=eps, I_0=I_0, rho=rho, sigma=sigma,
                    I_H1=I_H1, v0=v0, w0=w0,
                    n_traj=n_traj, dt=dt, rng=rng,
                )
                I_spike_grid[i, jr, k, :] = Is
                done += 1
                if verbose:
                    fired = np.mean(~np.isnan(Is))
                    delay = np.nanmean(Is) - I_H1 if fired > 0 else np.nan
                    print(
                        f"    [{done:4d}/{total}] rho=eps^{exp:.1f}={rho:.5f} "
                        f"sigma={sigma:.4f}  fired={fired:.2f}  "
                        f"mean_delay={delay:+.4f}"
                    )

    if verbose:
        print(f"\n  total wall time: {time.time() - t0:.1f}s")

    return dict(
        eps_vals   = np.asarray(eps_vals),
        rho_vals   = rho_grid,          # (n_eps, n_rho)
        rho_exponents = np.asarray(rho_exponents),
        sigma_vals = np.asarray(sigma_vals),
        I_spike    = I_spike_grid,      # (n_eps, n_rho, n_sigma, n_traj)
        I_H1_vals  = I_H1_grid,
        I_0_vals   = I_0_grid,
    )


# ---------------------------------------------------------------------------
# Figure  (2 × 2)
#
#   [0,0]  D(sigma) curves, one line per (eps, rho), log-log
#   [0,1]  D(sigma=0) vs eps for each rho exponent — power-law check
#   [1,0]  sigma_crit(eps) from D crossing threshold — exponent fit
#   [1,1]  firing rate (fraction non-NaN) vs sigma, per eps
# ---------------------------------------------------------------------------

def make_figure(data: dict, out_path: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    eps_vals      = data["eps_vals"]
    rho_vals      = data["rho_vals"]       # (n_eps, n_rho)
    rho_exponents = data["rho_exponents"]
    sigma_vals    = data["sigma_vals"]
    I_spike       = data["I_spike"]        # (n_eps, n_rho, n_sigma, n_traj)
    I_H1_vals     = data["I_H1_vals"]

    n_eps, n_rho = rho_vals.shape
    colors  = plt.cm.viridis(np.linspace(0.1, 0.9, n_eps))
    lstyles = ["-", "--", ":"]
    exp_labels = [f"ρ=ε^{e:.1f}" for e in rho_exponents]

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    # ---- panel [0,0]: mean delay D(sigma) ----
    ax = axes[0, 0]
    for i, eps in enumerate(eps_vals):
        I_H1 = I_H1_vals[i]
        for jr in range(n_rho):
            D = np.array([
                np.nanmean(I_spike[i, jr, k, :]) - I_H1
                for k in range(len(sigma_vals))
            ])
            mask = np.isfinite(D) & (D > 0)
            if mask.any():
                ax.plot(sigma_vals[mask], D[mask],
                        color=colors[i], ls=lstyles[jr], lw=1.5,
                        label=f"ε={eps:.2f} {exp_labels[jr]}" if i == 0 else None)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\sigma$")
    ax.set_ylabel(r"$D(\sigma) = \langle I_{\rm spike}\rangle - I_{H1}$")
    ax.set_title("Mean delay vs noise")
    ax.grid(True, which="both", alpha=0.2)
    # custom legend: eps by color, rho by linestyle
    color_handles = [Line2D([0],[0], color=colors[i], lw=2,
                             label=f"ε={eps_vals[i]:.2f}")
                     for i in range(n_eps)]
    ls_handles    = [Line2D([0],[0], color="k", ls=lstyles[jr], lw=1.5,
                             label=exp_labels[jr])
                     for jr in range(n_rho)]
    ax.legend(handles=color_handles + ls_handles, fontsize=7,
              frameon=False, ncol=2)

    # ---- panel [0,1]: deterministic delay D(sigma→0) vs eps ----
    ax = axes[0, 1]
    for jr in range(n_rho):
        D0 = []
        for i, eps in enumerate(eps_vals):
            I_H1 = I_H1_vals[i]
            # use smallest sigma point that has >50% firing
            for k in range(len(sigma_vals)):
                fired = np.mean(~np.isnan(I_spike[i, jr, k, :]))
                if fired >= 0.5:
                    D0.append(np.nanmean(I_spike[i, jr, k, :]) - I_H1)
                    break
            else:
                D0.append(np.nan)
        D0 = np.array(D0)
        mask = np.isfinite(D0) & (D0 > 0)
        if mask.sum() >= 2:
            ax.plot(eps_vals[mask], D0[mask], "o-", color=f"C{jr}",
                    lw=1.8, label=exp_labels[jr])
            # power-law fit
            log_e = np.log(eps_vals[mask])
            log_D = np.log(D0[mask])
            slope, intercept = np.polyfit(log_e, log_D, 1)
            e_fit = eps_vals[mask]
            ax.plot(e_fit, np.exp(intercept) * e_fit**slope, "--",
                    color=f"C{jr}", lw=1, alpha=0.6,
                    label=f"  slope={slope:.2f}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$\varepsilon$")
    ax.set_ylabel(r"$D_0 = D(\sigma \to 0)$")
    ax.set_title("Deterministic delay vs ε")
    ax.legend(fontsize=7, frameon=False)
    ax.grid(True, which="both", alpha=0.2)

    # ---- panel [1,0]: sigma_crit(eps) — where D drops by half ----
    ax = axes[1, 0]
    for jr in range(n_rho):
        sig_crit = []
        for i, eps in enumerate(eps_vals):
            I_H1 = I_H1_vals[i]
            D = np.array([
                np.nanmean(I_spike[i, jr, k, :]) - I_H1
                for k in range(len(sigma_vals))
            ])
            # D0: value at smallest sigma with >50% firing
            D0 = np.nan
            for k in range(len(sigma_vals)):
                fired = np.mean(~np.isnan(I_spike[i, jr, k, :]))
                if fired >= 0.5:
                    D0 = float(np.nanmean(I_spike[i, jr, k, :]) - I_H1)
                    break
            if not np.isfinite(D0) or D0 <= 0:
                sig_crit.append(np.nan)
                continue
            # find sigma where D first drops below D0/2
            sc = np.nan
            for k in range(len(sigma_vals)):
                if np.isfinite(D[k]) and D[k] < D0 / 2.0:
                    sc = sigma_vals[k]
                    break
            sig_crit.append(sc)
        sig_crit = np.array(sig_crit)
        mask = np.isfinite(sig_crit) & (sig_crit > 0)
        if mask.sum() >= 2:
            ax.plot(eps_vals[mask], sig_crit[mask], "s-", color=f"C{jr}",
                    lw=1.8, label=exp_labels[jr])
            log_e  = np.log(eps_vals[mask])
            log_sc = np.log(sig_crit[mask])
            slope, intercept = np.polyfit(log_e, log_sc, 1)
            e_fit = eps_vals[mask]
            ax.plot(e_fit, np.exp(intercept) * e_fit**slope, "--",
                    color=f"C{jr}", lw=1, alpha=0.6,
                    label=f"  slope={slope:.2f}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$\varepsilon$")
    ax.set_ylabel(r"$\sigma_{\rm crit}(\varepsilon)$")
    ax.set_title(r"$\sigma_{\rm crit}$ where $D$ halves")
    ax.legend(fontsize=7, frameon=False)
    ax.grid(True, which="both", alpha=0.2)

    # ---- panel [1,1]: firing rate vs sigma, per eps (rho=eps^1.0) ----
    ax = axes[1, 1]
    jr_mid = 1   # rho = eps^1.0
    for i, eps in enumerate(eps_vals):
        fired_rate = np.array([
            np.mean(~np.isnan(I_spike[i, jr_mid, k, :]))
            for k in range(len(sigma_vals))
        ])
        ax.plot(sigma_vals, fired_rate, "o-", color=colors[i],
                lw=1.5, label=f"ε={eps:.2f}")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\sigma$")
    ax.set_ylabel("Firing fraction")
    ax.set_title(r"Firing rate vs $\sigma$  ($\rho = \varepsilon^{1.0}$)")
    ax.set_ylim(-0.05, 1.05)
    ax.axhline(0.5, color="k", lw=0.7, ls="--", alpha=0.4)
    ax.legend(fontsize=8, frameon=False)
    ax.grid(True, which="both", alpha=0.2)

    fig.suptitle("Ramped canard passage: FHN with I(t) = I₀ + ρt", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure saved -> {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    out_data = os.path.join(_ROOT, "data")
    out_fig  = os.path.join(_ROOT, "figures")
    os.makedirs(out_data, exist_ok=True)
    os.makedirs(out_fig, exist_ok=True)

    print("\n=== Ramped canard passage experiment ===\n")
    print(f"  I_FOLD_L = {I_FOLD_L:.6f}")

    # quick sanity: print initial conditions for each eps
    fhn = FHN2D()
    print("\n  Initial conditions:")
    for eps in (0.02, 0.04, 0.08, 0.16):
        I_H1 = fhn.I_hopf_lower_at(eps)
        I_0  = max(I_FOLD_L, I_H1 - 4.0 * eps)
        fhn0 = FHN2D(I=I_0)
        v0, w0 = fhn0.V_FP, fhn0.W_FP
        print(f"    eps={eps:.2f}  I_H1={I_H1:.4f}  I_0={I_0:.4f}  "
              f"v0={v0:.4f}  w0={w0:.4f}")

    print()
    data = run_sweep(
        eps_vals   = (0.02, 0.04, 0.08, 0.16),
        sigma_vals = np.logspace(-3, -0.5, 12),
        n_traj     = 400,
        seed       = 42,
    )

    npz_path = os.path.join(out_data, "canard_ramp_passage.npz")
    np.savez(
        npz_path,
        eps_vals      = data["eps_vals"],
        rho_vals      = data["rho_vals"],
        rho_exponents = data["rho_exponents"],
        sigma_vals    = data["sigma_vals"],
        I_spike       = data["I_spike"],
        I_H1_vals     = data["I_H1_vals"],
        I_0_vals      = data["I_0_vals"],
    )
    print(f"  data saved -> {npz_path}")

    fig_path = os.path.join(out_fig, "canard_ramp_passage.png")
    make_figure(data, fig_path)

    # ---- print summary table ----
    eps_vals      = data["eps_vals"]
    rho_exponents = data["rho_exponents"]
    rho_vals      = data["rho_vals"]
    sigma_vals    = data["sigma_vals"]
    I_spike       = data["I_spike"]
    I_H1_vals     = data["I_H1_vals"]

    print("\nDelay D = mean(I_spike) - I_H1  at each (eps, rho, sigma):")
    print(f"  {'eps':>5s}  {'rho_exp':>7s}  {'rho':>9s}  "
          + "  ".join(f"σ={s:.3f}" for s in sigma_vals[[0, 3, 6, 9, 11]]))
    for i, eps in enumerate(eps_vals):
        I_H1 = I_H1_vals[i]
        for jr, exp in enumerate(rho_exponents):
            rho = rho_vals[i, jr]
            row = []
            for k in [0, 3, 6, 9, 11]:
                D = np.nanmean(I_spike[i, jr, k, :]) - I_H1
                fired = np.mean(~np.isnan(I_spike[i, jr, k, :]))
                row.append(f"{D:+.4f}({fired:.0%})")
            print(f"  {eps:5.2f}  {exp:7.1f}  {rho:9.5f}  "
                  + "  ".join(row))


if __name__ == "__main__":
    main()
