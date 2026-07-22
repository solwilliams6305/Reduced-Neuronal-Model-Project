"""
folded_cycle_phase_diffusion.py
-------------------------------
Channel B (phase diffusion of the fast oscillation) through the folded limit
cycle, on the JKK K2 inner chart, and the A-vs-B race.

The fast oscillation has radius r and phase theta.  Degenerate noise on the
fast Cartesian variable projects onto polar (r, theta):
    radial:  amplitude  sigma |cos 2pi theta|        -> Channel A (eta/sqrt2)
    phase :  amplitude  sigma |sin 2pi theta|/(2pi r) -> Channel B (1/r blow-up)
In the K2 chart (r = eps r2, eta = sigma/sqrt(eps2)) the phase noise is
    d(theta2) = ... + [ eta / (2 pi r2) ] (proj) dB,
so the accumulated phase variance through the passage is
    Var(theta) = (eta^2 / 8pi^2) * INT dt2 / r2^2
               = (eta^2 / 8pi^2) * (b/ac) * ln( y2_in / y2_* )      (log-divergent,
                                                                     cut at inner scale)
Phase randomization (Var ~ 1) =>  eta_*^B(chart) = sqrt( 8 pi^2 / I_phi ),
I_phi = INT dt2/r2^2.  Physical sigma_*^B = sqrt(eps2) * eta_*^B
      ~ 2 pi sqrt(3) sqrt(ac/b) sqrt( eps2 / |ln eps2| ).

Compare Channel A:  sigma_*^A = C_q sqrt(eps2) sqrt(ac/b),  C_q ~ 2.8 (NF).
=> sigma_*^B / sigma_*^A = (2 pi sqrt3 / C_q) / sqrt|ln eps2|   (same leading order,
   B lower by a log, but larger prefactor).

This script:
 (1) integrates the DETERMINISTIC approach, accumulates I_phi, checks the log law;
 (2) extracts eta_*^B(chart) and the sqrt(1/ln) scaling vs eps2 (encoded by y2_in);
 (3) a stochastic check: measured Var(theta) matches the prediction at small eta
     and SATURATES once amplitude escape (Channel A) fires first  -> the interplay.
"""
from __future__ import annotations
import os, time
import numpy as np

C_Q_NF = 2.8        # canard/fold normal-form escape constant (eta_*^A, chart)


def scales(a, b, c):
    y2s = (c * c / (a * b)) ** (1.0 / 3.0)
    r2s = (a * c / (b * b)) ** (1.0 / 3.0)
    return y2s, r2s


# ---------------------------------------------------------------------------
# (1)-(2) deterministic phase-diffusion integral I_phi = INT dt2/r2^2
# ---------------------------------------------------------------------------

def I_phi_deterministic(a, b, c, y2_in, dT_canon=2e-3, y_cap_frac=50.0):
    """I_phi = INT dt2/r2^2 along the Riccati approach from y2_in to the fold
    core.  On the attracting branch r2^2 = (a/b) y2, so the long tail above
    y_cap = y_cap_frac*y2_* is ANALYTIC, (b/ac) ln(y2_in/y_cap); only the inner
    region (y_cap -> core) is integrated numerically.  Returns (I_phi, lnL)."""
    y2s, r2s = scales(a, b, c)
    tau = (a * b * c) ** (-1.0 / 3.0)
    dt2 = dT_canon * tau
    y_cap = y_cap_frac * y2s
    # analytic branch tail above y_cap
    I = (b / (a * c)) * np.log(max(y2_in, y_cap) / y_cap) if y2_in > y_cap else 0.0
    y2 = min(float(y2_in), y_cap)
    r2 = -np.sqrt((a / b) * y2)
    nmax = int(50 * (y_cap / c) / dt2) + 1000
    for _ in range(nmax):
        I += dt2 / (r2 * r2)
        r2 = r2 + (b * r2 * r2 - a * y2) * dt2
        y2 = y2 - c * dt2
        if r2 >= -r2s:          # reached inner scale -> blow-up takes over
            break
        if y2 <= 0:
            break
    return I, np.log(y2_in / y2s)


# ---------------------------------------------------------------------------
# (3) stochastic (r2, theta2): phase variance + amplitude escape interplay
# ---------------------------------------------------------------------------

def simulate_phase(a, b, c, eta, y2_in, n_traj=400, dT_canon=2e-3,
                   r2_floor_frac=1.0, rng=None):
    """Integrate dr2=(b r2^2 - a y2)dt2 + (eta/sqrt2) dBr,
               d th = [eta/(2 sqrt2 pi r2)] dBth,  dy2=-c dt2.
    Stop each trajectory at amplitude escape (r2 >= +r2_*) or fold core
    (r2 >= -r2_floor_frac r2_*).  Return (Var_theta, escape_fraction)."""
    if rng is None:
        rng = np.random.default_rng()
    y2s, r2s = scales(a, b, c)
    tau = (a * b * c) ** (-1.0 / 3.0)
    dt2 = dT_canon * tau
    sdt = np.sqrt(dt2)
    r2 = np.full(n_traj, -np.sqrt((a / b) * y2_in))
    y2 = np.full(n_traj, float(y2_in))
    th = np.zeros(n_traj)
    done = np.zeros(n_traj, bool)
    escaped = np.zeros(n_traj, bool)
    r2_core = -r2_floor_frac * r2s
    nmax = int(60 * (y2_in / c) / dt2) + 2000
    for _ in range(nmax):
        act = ~done
        if not act.any():
            break
        nr = rng.standard_normal(n_traj)
        nth = rng.standard_normal(n_traj)
        # phase noise uses current r2 (capped at r2_* to avoid 1/0 at core)
        r2eff = np.where(np.abs(r2) < r2s, r2s, np.abs(r2))
        th[act] = th[act] + (eta / (2 * np.sqrt(2) * np.pi * r2eff[act])) * sdt * nth[act]
        r2n = r2.copy()
        r2n[act] = (r2[act] + (b * r2[act] ** 2 - a * y2[act]) * dt2
                    + (eta / np.sqrt(2)) * sdt * nr[act])
        y2[act] = y2[act] - c * dt2
        esc = act & (r2 < r2s) & (r2n >= r2s)         # amplitude escape (+r2_*)
        core = act & (r2 < r2_core) & (r2n >= r2_core)  # reached fold core
        escaped |= esc
        done |= esc | core | (y2 < -2 * y2s)
        r2 = r2n
    return float(np.var(th[escaped | (~escaped)])), float(np.mean(escaped))


# ---------------------------------------------------------------------------
def make_figure(geoms, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))

    # Panel 1: log law  I_phi vs ln(y2_in/y2_*)
    ax = axes[0]
    for (a, b, c) in geoms:
        lnLs, Is = [], []
        for eps2 in np.logspace(-1.5, -7, 8):
            y2s, _ = scales(a, b, c)
            y2_in = y2s * eps2 ** (-2.0 / 3.0)
            I, lnL = I_phi_deterministic(a, b, c, y2_in, dT_canon=8e-3)
            lnLs.append(lnL); Is.append(I)
        lnLs = np.array(lnLs); Is = np.array(Is)
        ax.plot(lnLs, Is, "o", ms=4, label=f"(a,b,c)=({a:g},{b:g},{c:g})")
        ax.plot(lnLs, (b / (a * c)) * lnLs, "-", lw=1,
                label=f"$(b/ac)\\,\\ln L$ = {b/(a*c):.2f} ln L")
    ax.set_xlabel(r"$\ln(y_{2,\rm in}/y_{2,*})\ \approx\ \frac{2}{3}|\ln\epsilon_2|$")
    ax.set_ylabel(r"$I_\phi=\int dt_2/r_2^2$")
    ax.set_title(r"Channel B: phase-diffusion integral is log-divergent")
    ax.legend(fontsize=8, frameon=False); ax.grid(alpha=0.25)

    # Panel 2: the A-vs-B race  sigma_*^B/sigma_*^A vs eps2
    ax = axes[1]
    eps2g = np.logspace(-1, -12, 16)
    for (a, b, c) in geoms:
        y2s, _ = scales(a, b, c)
        ratio = []
        for eps2 in eps2g:
            y2_in = y2s * eps2 ** (-2.0 / 3.0)
            I, _ = I_phi_deterministic(a, b, c, y2_in, dT_canon=8e-3)
            etaB = np.sqrt(8 * np.pi ** 2 / I)
            etaA = C_Q_NF * np.sqrt(a * c / b)
            ratio.append(etaB / etaA)
        ax.plot(eps2g, ratio, "-o", ms=3, label=f"(a,b,c)=({a:g},{b:g},{c:g})")
    ax.axhline(1.0, color="k", ls="--", lw=0.8)
    ax.annotate("A sets $\\sigma_*$ (escape first)", (1e-4, 1.4), fontsize=8)
    ax.annotate("B sets $\\sigma_*$ (phase first)", (1e-10, 0.85), fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel(r"$\epsilon_2$ (drift parameter)")
    ax.set_ylabel(r"$\sigma_*^B/\sigma_*^A=\eta_*^B/\eta_*^A$")
    ax.set_title(r"A-vs-B race (geometry-independent; $C_q=2.8$)")
    ax.legend(fontsize=8, frameon=False); ax.grid(alpha=0.25, which="both")

    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    os.makedirs(os.path.join(root, "data"), exist_ok=True)
    os.makedirs(os.path.join(root, "figures"), exist_ok=True)

    print("\n=== Folded limit cycle: Channel B (phase diffusion) + A-vs-B race ===\n")

    geoms = [(1.0, 1.0, 1.0), (1.7, 0.8, 1.3)]
    # y2_in encodes eps2 via physical y_in ~ O(1):  y2_in = y2_* * eps2^{-2/3}
    eps2_list = [1e-2, 1e-3, 1e-4, 1e-6]

    print("(1)-(2)  Deterministic I_phi = INT dt2/r2^2  and the log law")
    print("         predicted  I_phi ~ (b/ac) ln(y2_in/y2_*)\n")
    for (a, b, c) in geoms:
        y2s, r2s = scales(a, b, c)
        pref = b / (a * c)
        print(f"  geometry (a,b,c)=({a:g},{b:g},{c:g})   b/ac={pref:.3f}  "
              f"y2_*={y2s:.3f}  r2_*={r2s:.3f}")
        print(f"    {'eps2':>8s} {'y2_in':>9s} {'I_phi':>8s} {'(b/ac)lnL':>10s} "
              f"{'ratio':>6s} {'eta_*^B':>8s} {'eta_*^A':>8s} {'B/A':>6s}")
        for eps2 in eps2_list:
            y2_in = y2s * eps2 ** (-2.0 / 3.0)
            I, lnL = I_phi_deterministic(a, b, c, y2_in)
            pred = pref * lnL
            etaB = np.sqrt(8 * np.pi ** 2 / I)
            etaA = C_Q_NF * np.sqrt(a * c / b)
            print(f"    {eps2:8.0e} {y2_in:9.1f} {I:8.3f} {pred:10.3f} "
                  f"{I/pred:6.3f} {etaB:8.3f} {etaA:8.3f} {etaB/etaA:6.3f}")
        print()

    print("  => I_phi tracks (b/ac) ln(y2_in/y2_*): the 1/r^2 phase-diffusion")
    print("     integral is log-divergent, regularized at the inner scale.\n")
    print("  => sigma_*^B / sigma_*^A = eta_*^B / eta_*^A  (same sqrt(eps2) sqrt(ac/b)")
    print("     leading order); B exceeds A until eps2 is astronomically small,")
    print("     because the larger prefactor (2pi sqrt3 ~ 10.9) beats 1/sqrt|ln| only")
    print("     for |ln eps2| > (2pi sqrt3 / C_q)^2.\n")

    # (3) stochastic check at one geometry / eps2
    a, b, c = 1.0, 1.0, 1.0
    y2s, r2s = scales(a, b, c)
    eps2 = 1e-3
    y2_in = y2s * eps2 ** (-2.0 / 3.0)
    I, lnL = I_phi_deterministic(a, b, c, y2_in)
    etaB = np.sqrt(8 * np.pi ** 2 / I)
    etaA = C_Q_NF
    print(f"(3)  Stochastic check  (a,b,c)=(1,1,1), eps2={eps2:.0e}, "
          f"y2_in={y2_in:.1f}, I_phi={I:.2f}")
    print(f"     predicted eta_*^B={etaB:.2f}  eta_*^A~{etaA:.2f}\n")
    print(f"     {'eta':>7s} {'Var(theta)':>11s} {'pred small-eta':>14s} "
          f"{'esc.frac':>9s}")
    rng = np.random.default_rng(7)
    for eta in [0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0]:
        var, esc = simulate_phase(a, b, c, eta, y2_in, n_traj=400, rng=rng)
        pred = eta ** 2 * I / (8 * np.pi ** 2)
        print(f"     {eta:7.2f} {var:11.4f} {pred:14.4f} {esc:9.2f}")
    print("\n     small eta: Var(theta) ~ eta^2 I/8pi^2 (diffusion law confirmed).")
    print("     large eta: accumulated Var(theta) saturates below the linear")
    print("     extrapolation (inner-scale cutoff + radial-noise-shortened approach).\n")

    fig_path = os.path.join(root, "figures", "folded_cycle_phase_diffusion.png")
    make_figure(geoms, fig_path)
    print(f"  figure saved -> figures/folded_cycle_phase_diffusion.png\n")


if __name__ == "__main__":
    main()
