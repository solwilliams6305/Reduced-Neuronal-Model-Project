"""
alpha1_averaged_law.py — Part 3 Task 1-2: which average does alpha=1 select?
----------------------------------------------------------------------------
ALPHA1_METHODOLOGY.md refuted the MIN-barrier route; fast rotation AVERAGES.
This pins WHICH average, by derivation + a discriminating test.

DERIVATION (Khasminskii stochastic averaging).  Rotating-phase inner SDE
    dr = ( b(theta) r^2 - a y ) dT + eta dB,   theta = omega T,   dy = -c dT.
When theta rotates fast, (r,y) see the time-AVERAGED drift, so the DRIFT
averages term-by-term:  dr = ( <b> r^2 - a y ) dT + eta dB.  The effective
autonomous canard therefore has the ARITHMETIC-MEAN coefficient <b>, giving
    eta_*(alpha=1) = C_q sqrt( a c / <b> )            (COEFFICIENT-averaging),
NOT the average of the barrier G^2 = ac/b  (that would give ac*<1/b>), and NOT
min.  We modulate b (it enters G as 1/b) so these candidates separate:

    b(theta) = b0 (1 + delta cos 2pi theta),  <b> = b0,  <1/b> = 1/(b0 sqrt(1-delta^2)).

Predicted ratio  eta_*(fast, modulated b) / eta_*(frozen, b=b0):
    coefficient-avg (THIS):  sqrt(b0/<b>)        = 1.000
    barrier-avg (RMS of G):  sqrt(b0*<1/b>)      = 1/(1-delta^2)^{1/4}
    min  (REFUTED):          sqrt(b0/b_max)      = 1/sqrt(1+delta)

TEST: measure eta_* at fast omega and check the ratio -> 1.0 (coeff-avg), flat in
omega.  Observable: <y_escape> (y at first r-crossing of r_cross) crosses 0 at
eta_* (the canard convention; W_* = 1 here since a=c=1).
"""
from __future__ import annotations
import os
import numpy as np

A, C, B0, DELTA = 1.0, 1.0, 1.0, 0.8


def b_of(theta):
    return B0 * (1.0 + DELTA * np.cos(2 * np.pi * theta))


def simulate(omega, delta, eta, N=120, y0=5.0, dT=8e-4, T_max=7.0,
             r_cross=1.0, seed=0):
    """dr=(b r^2 - a y)dT + eta dB, dtheta=omega dT, dy=-c dT.
    Returns mean y at first up-crossing r>=r_cross (NaN-safe mean over escaped)."""
    rng = np.random.default_rng(seed)
    bmean = B0                                   # <b> = b0 (cos has zero mean)
    r = np.full(N, -np.sqrt(A * y0 / bmean))
    y = np.full(N, y0)
    th = rng.random(N)                            # random initial phase
    yesc = np.full(N, np.nan)
    done = np.zeros(N, bool)
    sdt = np.sqrt(dT)
    n = int(T_max / dT)
    for _ in range(n):
        act = ~done
        if not act.any():
            break
        b = B0 * (1.0 + delta * np.cos(2 * np.pi * th))
        r_new = r.copy()
        r_new[act] = r[act] + (b[act] * r[act] ** 2 - A * y[act]) * dT \
            + eta * sdt * rng.standard_normal(N)[act]
        y[act] = y[act] - C * dT
        th[act] = th[act] + omega * dT
        cross = act & (r < r_cross) & (r_new >= r_cross)
        if cross.any():
            yesc[cross] = y[cross]
            done |= cross
        done |= (y < -2.0)
        r = r_new
    return np.nanmean(yesc)


def eta_star(etas, ymean, target=0.0):
    """eta at which <y_escape> crosses target (rising)."""
    for i in range(len(etas) - 1):
        if ymean[i] < target <= ymean[i + 1]:
            f = (target - ymean[i]) / (ymean[i + 1] - ymean[i])
            return etas[i] + f * (etas[i + 1] - etas[i])
    return np.nan


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    figdir = os.path.join(os.path.dirname(here), "figures")
    os.makedirs(figdir, exist_ok=True)

    inv_b_mean = 1.0 / (B0 * np.sqrt(1 - DELTA ** 2))     # <1/b>
    r_coeff = 1.0
    r_barrier = (1 - DELTA ** 2) ** (-0.25)
    r_min = 1.0 / np.sqrt(1 + DELTA)
    print("\n=== Part 3 Task 1-2: which average does alpha=1 select? ===")
    print(f"  b(theta)=b0(1+{DELTA}cos),  <b>={B0}, <1/b>={inv_b_mean:.3f}, b_max={B0*(1+DELTA)}")
    print(f"  predicted eta_*(fast)/eta_*(frozen b0):")
    print(f"     coefficient-avg (THIS derivation): {r_coeff:.3f}")
    print(f"     barrier-avg (RMS of G):            {r_barrier:.3f}")
    print(f"     min (REFUTED):                     {r_min:.3f}\n")

    etas = np.array([1.5, 2.0, 2.5, 3.0, 3.5, 4.5, 6.0])
    N = 120

    # reference: frozen, unmodulated b=b0
    yref = np.array([simulate(0.0, 0.0, e, N=N, seed=k) for k, e in enumerate(etas)])
    es_ref = eta_star(etas, yref)
    print(f"  REFERENCE (frozen b=b0):  eta_* = {es_ref:.3f}")
    print(f"     <y_esc>(eta) = [" + ", ".join(f"{v:+.3f}" for v in yref) + "]\n")

    # modulated b, sweep omega
    print(f"  MODULATED b (delta={DELTA}), eta_*(omega) and ratio to reference:")
    omegas = [0.0, 5.0, 20.0, 50.0]
    es_list = []
    for om in omegas:
        ym = np.array([simulate(om, DELTA, e, N=N, seed=100 + int(om) + k)
                       for k, e in enumerate(etas)])
        es = eta_star(etas, ym)
        es_list.append(es)
        print(f"     omega={om:5.1f}:  eta_*={es:.3f}   ratio={es/es_ref:.3f}")

    es_fast = es_list[-1]
    ratio = es_fast / es_ref
    print(f"\n  FAST (omega=50) ratio = {ratio:.3f}")
    cands = {"coefficient-avg (THIS)": r_coeff, "barrier-avg (RMS-G)": r_barrier,
             "min (refuted)": r_min}
    best = min(cands, key=lambda k: abs(cands[k] - ratio))
    print(f"  closest candidate: {best}  (predicted {cands[best]:.3f})")
    print(f"  => alpha=1 selects {'COEFFICIENT' if best.startswith('coeff') else best} averaging:"
          f"  eta_* = C_q sqrt(a c / <b>)\n")

    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    om_arr = np.array(omegas); es_arr = np.array(es_list)
    ax[0].plot(om_arr, es_arr / es_ref, "o-", color="C0", ms=6, label="measured ratio")
    ax[0].axhline(r_coeff, color="C2", ls="-", lw=1.2, label=f"coeff-avg {r_coeff:.2f} (THIS)")
    ax[0].axhline(r_barrier, color="C3", ls="--", lw=1.2, label=f"barrier-avg (RMS-G) {r_barrier:.2f}")
    ax[0].axhline(r_min, color="C1", ls=":", lw=1.2, label=f"min {r_min:.2f} (refuted)")
    ax[0].set_xlabel(r"rotation rate $\omega$")
    ax[0].set_ylabel(r"$\eta_*(\omega)\,/\,\eta_*^{\rm frozen}$")
    ax[0].set_title(r"$\alpha=1$: which average? (modulate $b$)")
    ax[0].legend(fontsize=8, frameon=False); ax[0].grid(alpha=0.3)

    ax[1].plot(etas, yref, "ks-", lw=2, ms=5, label="frozen $b_0$ (ref)")
    for om, es in zip(omegas, es_list):
        ym = np.array([simulate(om, DELTA, e, N=N, seed=400 + int(om) + k)
                       for k, e in enumerate(etas)])
        ax[1].plot(etas, ym, "o-", ms=3, alpha=0.8, label=fr"$\omega$={om:.0f}")
    ax[1].axhline(0.0, color="k", lw=0.6)
    ax[1].set_xlabel(r"$\eta$"); ax[1].set_ylabel(r"$\langle y_{\rm escape}\rangle$")
    ax[1].set_title(r"Fast-$\omega$ curve matches frozen-$b_0$ (coeff-avg)")
    ax[1].legend(fontsize=8, frameon=False); ax[1].grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(figdir, "alpha1_averaged_law.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print(f"  figure -> figures/alpha1_averaged_law.png\n")


if __name__ == "__main__":
    main()
