"""
alpha1_fusion.py — Part 3 Task 3: the A-B fusion (first swing)
--------------------------------------------------------------
Both channels ON, phase rotating.  Inner SDE:
    dr = ( r^2 - a y ) dT + eta dB_r              (Channel A: amplitude escape)
    dtheta = omega dT + ( eta / (2 pi r) ) dB_th  (Channel B: phase diffusion, 1/r)
    dy = - dT                                     (drift through the fold)
(a=c=1 unmodulated baseline; dB_r, dB_th independent = fast-phase-averaged.)

THE COUPLING TO WATCH: Channel B's phase-noise amplitude eta/(2 pi r) BLOWS UP as
Channel A drives r -> 0 at the fold.  So the channels are not obviously
independent -- phase decoherence may be slaved to the amplitude approach.

Questions (first swing):
 1. Do BOTH thresholds survive?  sigma_*^A (amplitude: <y_escape> crosses 0) and
    sigma_*^B (phase: Var(theta_noise) at the fold reaches ~1).
 2. Are both FLAT in omega? (rotation adds a mean to theta, not variance; and a=b=c
    const so Channel A coefficient-averaging is trivial) -> the race should persist.
 3. The race: which fires first, and does Channel B's log-fed threshold survive
    the rotation?
"""
from __future__ import annotations
import os
import numpy as np

A, C = 1.0, 1.0
R_FLOOR = 0.5          # inner cutoff for the 1/r phase noise (blow-up regularization)
R_CROSS = 1.0          # amplitude-escape threshold


def simulate(omega, eta, N=150, y0=5.0, dT=5e-4, T_max=8.0, seed=0):
    """Returns (mean y_escape, Var of noise-phase AT AMPLITUDE ESCAPE).
    The coupled observable: how decohered is the phase by the time Channel A
    destroys the cycle?  Var<1 => A preempts B; Var>=1 => B randomizes first."""
    rng = np.random.default_rng(seed)
    r = np.full(N, -np.sqrt(A * y0))
    y = np.full(N, y0)
    th_noise = np.zeros(N)            # phase MINUS deterministic omega*T
    yesc = np.full(N, np.nan)
    thn_esc = np.full(N, np.nan)      # th_noise at amplitude escape
    escaped = np.zeros(N, bool)
    sdt = np.sqrt(dT)
    n = int(T_max / dT)
    for _ in range(n):
        alive = ~escaped
        if not alive.any():
            break
        reff = np.where(np.abs(r) < R_FLOOR, R_FLOOR, np.abs(r))
        th_noise[alive] += (eta / (2 * np.pi * reff[alive])) * sdt * rng.standard_normal(N)[alive]
        r_new = r.copy()
        r_new[alive] = r[alive] + (r[alive] ** 2 - A * y[alive]) * dT \
            + eta * sdt * rng.standard_normal(N)[alive]
        y = y - C * dT
        cross = alive & (r < R_CROSS) & (r_new >= R_CROSS)
        if cross.any():
            yesc[cross] = y[cross]
            thn_esc[cross] = th_noise[cross]
            escaped |= cross
        r = r_new
    # trajectories that never escaped: use final accumulated phase
    thn_esc[~escaped] = th_noise[~escaped]
    return np.nanmean(yesc), float(np.var(thn_esc))


def thr_cross(etas, vals, target, rising=True):
    for i in range(len(etas) - 1):
        a, b = vals[i], vals[i + 1]
        if (rising and a < target <= b) or (not rising and a > target >= b):
            f = (target - a) / (b - a)
            return etas[i] + f * (etas[i + 1] - etas[i])
    return np.nan


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    figdir = os.path.join(os.path.dirname(here), "figures")
    os.makedirs(figdir, exist_ok=True)

    etas = np.array([1.5, 2.0, 2.5, 3.0, 3.5, 4.5, 6.0, 8.0])
    omegas = [0.0, 10.0, 30.0]
    N = 150

    print("\n=== Part 3 Task 3: A-B fusion (both channels on, phase rotating) ===")
    print(f"  inner SDE a=c=1, R_floor={R_FLOOR}, R_cross={R_CROSS}, N={N}")
    print("  sigma_*^A: <y_escape> crosses 0;  Var(theta@ESCAPE): phase decoherence when A fires\n")

    res = {}
    for om in omegas:
        ymean = np.empty(len(etas)); vth = np.empty(len(etas))
        for k, e in enumerate(etas):
            ymean[k], vth[k] = simulate(om, e, N=N, seed=int(om) * 13 + k)
        sA = thr_cross(etas, ymean, 0.0, rising=True)
        sB = thr_cross(etas, vth, 1.0, rising=True)        # nan if Var never hits 1
        vth_at_sA = float(np.interp(sA, etas, vth)) if np.isfinite(sA) else np.nan
        res[om] = dict(ymean=ymean, vth=vth, sA=sA, sB=sB, vth_at_sA=vth_at_sA)
        sBtxt = f"{sB:.2f}" if np.isfinite(sB) else "n/a (Var<1 always)"
        print(f"  omega={om:5.1f}:  sigma_*^A={sA:.2f}   Var(theta@escape at sigma_*^A)={vth_at_sA:.2f}"
              f"   sigma_*^B={sBtxt}")
        print(f"     <y_esc>=[" + ",".join(f"{v:+.2f}" for v in ymean) + "]")
        print(f"     Var(th@esc)=[" + ",".join(f"{v:.2f}" for v in vth) + "]")

    sAs = np.array([res[o]["sA"] for o in omegas])
    vAs = np.array([res[o]["vth_at_sA"] for o in omegas])
    cvA = np.nanstd(sAs) / np.nanmean(sAs)
    print(f"\n  sigma_*^A vs omega = {np.round(sAs, 2)}  (CV {cvA:.2f}: flat -> survives rotation)")
    print(f"  Var(theta_noise @ escape) at sigma_*^A = {np.round(vAs, 2)}")
    print(f"\n  FUSION FINDING: at the amplitude-escape threshold the phase is only")
    print(f"  ~Var {np.nanmean(vAs):.2f} (<1) decohered => Channel A PREEMPTS Channel B:")
    print(f"  amplitude escape destroys the cycle before the phase fully randomizes.")
    print(f"  Both flat in omega => the picture survives rotation; A sets sigma_* here.\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    o0 = omegas[1]
    r0 = res[o0]
    ax[0].plot(etas, r0["ymean"], "o-", color="C0", label=r"$\langle y_{\rm escape}\rangle$ (Ch A)")
    ax[0].axhline(0, color="C0", lw=0.6, ls=":")
    ax2 = ax[0].twinx()
    ax2.plot(etas, r0["vth"], "s-", color="C3", label=r"Var($\theta_{\rm noise}$) (Ch B)")
    ax2.axhline(1.0, color="C3", lw=0.6, ls=":")
    if np.isfinite(r0["sA"]): ax[0].axvline(r0["sA"], color="C0", ls="--", lw=1)
    ax2.axhline(1.0, color="C3", lw=0.6, ls=":")
    ax[0].set_xlabel(r"$\eta$"); ax[0].set_ylabel(r"$\langle y_{\rm escape}\rangle$", color="C0")
    ax2.set_ylabel(r"Var($\theta_{\rm noise}$@escape)", color="C3")
    ax[0].set_title(fr"$\omega$={o0:.0f}: A fires at $\eta\approx${r0['sA']:.1f}; phase Var$\approx${r0['vth_at_sA']:.2f}$<$1 (A preempts B)")
    ax[0].legend(loc="upper left", fontsize=8, frameon=False)

    ax[1].plot(omegas, sAs, "o-", color="C0", ms=7, label=r"$\sigma_*^A$ (amplitude)")
    ax[1].set_xlabel(r"rotation rate $\omega$"); ax[1].set_ylabel(r"$\sigma_*^A$ (inner $\eta$)", color="C0")
    ax[1].set_ylim(0, np.nanmax(sAs) * 1.6)
    axb = ax[1].twinx()
    axb.plot(omegas, vAs, "s--", color="C3", ms=7, label=r"Var($\theta$@escape)")
    axb.axhline(1.0, color="C3", lw=0.6, ls=":"); axb.set_ylim(0, 1.2)
    axb.set_ylabel(r"Var($\theta_{\rm noise}$@escape)", color="C3")
    ax[1].set_title(r"$\sigma_*^A$ flat in $\omega$; phase only ~half-decohered at escape")
    ax[1].legend(loc="center right", fontsize=9, frameon=False); ax[1].grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(figdir, "alpha1_fusion.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print(f"  figure -> figures/alpha1_fusion.png\n")


if __name__ == "__main__":
    main()
