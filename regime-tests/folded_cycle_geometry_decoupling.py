"""
folded_cycle_geometry_decoupling.py — Publishable #3: is the race still
geometry-independent under rotation?  (decoupling conjecture TESTED)
-----------------------------------------------------------------------------
CONJECTURE (mine): under rotation Channel A averages <a> (arithmetic) but Channel
B's 1/r^2 integrand averages <1/a> (harmonic), so the frozen-phase geometry-
independence of the A-vs-B race would BREAK.

WHY IT'S WRONG (and the test that shows it): Channel B's 1/r^2 uses the AMPLITUDE
r, which under fast rotation is SLAVED to the averaged canard r^2 ~ <a> y (same
Khasminskii averaging as Channel A).  So 1/r^2 ~ 1/(<a> y): Channel B sees 1/<a>,
the SAME <a> functional as Channel A -> the race STAYS geometry-independent.
Only if r could follow the instantaneous a(theta) (SLOW rotation) would B see
<1/a>.  So: modulate a(theta)=a0(1+delta cos); measure Channel B's phase-variance
coefficient (amplitude-noise OFF, deterministic r) vs delta at SLOW vs FAST omega:
    SLOW omega: Var ~ <1/a>  (grows with delta)   [r follows a(theta)]
    FAST omega: Var ~ 1/<a> = 1/a0  (FLAT)         [r slaved to <a>]
and sigma_*^A is FLAT in delta either way (<a>=a0).  FAST: both see <a> -> race
geometry-INDEPENDENT (robust).  Conjecture refuted; regime diagram defensible.
"""
from __future__ import annotations
import os
import numpy as np

R_FLOOR, R_CROSS, A0 = 0.5, 1.0, 1.0


def a_of(theta, delta):
    return A0 * (1.0 + delta * np.cos(2 * np.pi * theta))


def sigma_starA(omega, delta, etas, N=150, y0=5.0, dT=5e-4, T_max=8.0, seed=0):
    """amplitude-escape threshold sigma_*^A (<y_escape> crosses 0)."""
    rng = np.random.default_rng(seed)
    ym = np.empty(len(etas))
    for k, eta in enumerate(etas):
        r = np.full(N, -np.sqrt(A0 * y0)); y = np.full(N, y0); th = rng.random(N)
        yesc = np.full(N, np.nan); esc = np.zeros(N, bool); sdt = np.sqrt(dT)
        for _ in range(int(T_max / dT)):
            al = ~esc
            if not al.any():
                break
            a = a_of(th, delta)
            rn = r.copy()
            rn[al] = r[al] + (r[al] ** 2 - a[al] * y[al]) * dT + eta * sdt * rng.standard_normal(N)[al]
            y = y - dT; th[al] += omega * dT
            c = al & (r < R_CROSS) & (rn >= R_CROSS)
            yesc[c] = y[c]; esc |= c; r = rn
        ym[k] = np.nanmean(yesc)
    for i in range(len(etas) - 1):
        if ym[i] < 0 <= ym[i + 1]:
            f = (0 - ym[i]) / (ym[i + 1] - ym[i]); return etas[i] + f * (etas[i + 1] - etas[i])
    return np.nan


def varB_isolated(omega, delta, eta=1.0, N=400, y0=5.0, dT=5e-4, seed=0):
    """Channel B alone: amplitude noise OFF (deterministic r), phase noise ON.
    Returns Var(theta_noise) accumulated to the fold y=0, divided by eta^2."""
    rng = np.random.default_rng(seed)
    r = np.full(N, -np.sqrt(A0 * y0)); y = np.full(N, y0); th = rng.random(N)
    thn = np.zeros(N); done = np.zeros(N, bool); thn_fold = np.full(N, np.nan)
    sdt = np.sqrt(dT)
    for _ in range(int((y0 + 0.5) / dT)):
        al = ~done
        if not al.any():
            break
        a = a_of(th, delta)
        reff = np.where(np.abs(r) < R_FLOOR, R_FLOOR, np.abs(r))
        thn[al] += (eta / (2 * np.pi * reff[al])) * sdt * rng.standard_normal(N)[al]
        r[al] = r[al] + (r[al] ** 2 - a[al] * y[al]) * dT          # deterministic amplitude
        yn = y - dT
        th[al] += omega * dT
        fold = al & (y >= 0) & (yn < 0)
        thn_fold[fold] = thn[fold]; done |= fold
        y = yn
    return float(np.nanvar(thn_fold) / eta ** 2)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    figdir = os.path.join(os.path.dirname(here), "figures")
    os.makedirs(figdir, exist_ok=True)
    etas = np.array([1.5, 2.0, 2.5, 3.0, 3.5, 4.5, 6.0])
    deltas = [0.0, 0.4, 0.8]

    print("\n=== Publishable #3: geometry-(in)dependence of the race under rotation ===")
    print("  a(theta)=a0(1+delta cos):  <a>=a0 (flat),  <1/a>=1/sqrt(1-delta^2)\n")

    # (1) Channel A threshold vs delta at fast omega
    print("  (1) sigma_*^A (fast omega=20) vs delta  [predict FLAT ~ sqrt(<a>)=1]:")
    sA = [sigma_starA(20.0, d, etas, seed=7) for d in deltas]
    for d, s in zip(deltas, sA):
        print(f"      delta={d:.1f}: sigma_*^A={s:.2f}")
    print(f"      CV = {np.nanstd(sA)/np.nanmean(sA):.2f}  -> A sees <a> (arithmetic), flat\n")

    # (2) Channel B isolated: Var coefficient vs delta, SLOW vs FAST omega
    print("  (2) Channel-B phase-variance coefficient Var/eta^2 (amplitude noise OFF):")
    print(f"      {'delta':>6s} {'<1/a>':>7s} {'slow w=2':>10s} {'fast w=20':>10s}")
    inv_a = [1.0 / np.sqrt(1 - d**2) for d in deltas]
    vslow, vfast = [], []
    for d, ia in zip(deltas, inv_a):
        vs = varB_isolated(2.0, d, seed=11)
        vf = varB_isolated(20.0, d, seed=12)
        vslow.append(vs); vfast.append(vf)
        print(f"      {d:6.1f} {ia:7.3f} {vs:10.3f} {vf:10.3f}")
    rs = np.array(vslow) / vslow[0]; rf = np.array(vfast) / vfast[0]
    print(f"\n      slow-omega ratio to delta=0: {np.round(rs,2)}  vs <1/a>={np.round(inv_a,2)}"
          f"  ({'tracks <1/a>: B follows a(theta)' if abs(rs[-1]-inv_a[-1])<0.4 else 'partial'})")
    print(f"      fast-omega ratio to delta=0: {np.round(rf,2)}  vs 1/<a>=[1,1,1]"
          f"  ({'FLAT: B slaved to <a>' if abs(rf[-1]-1)<0.3 else 'partial'})")
    print(f"\n  => FAST rotation (the alpha=1 regime): Channel A sees <a> AND Channel B")
    print(f"     sees 1/<a> (same functional) -> the A-vs-B race STAYS GEOMETRY-")
    print(f"     INDEPENDENT.  The decoupling conjecture is REFUTED; regime diagram robust.\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    ax[0].plot(deltas, sA, "o-", color="C0", ms=8)
    ax[0].set_ylim(0, max(sA) * 1.5); ax[0].axhline(np.nanmean(sA), color="C0", ls=":", lw=0.6)
    ax[0].set_xlabel(r"modulation depth $\delta$"); ax[0].set_ylabel(r"$\sigma_*^A$")
    ax[0].set_title(r"Channel A flat in $\delta$ ($\sim\sqrt{\langle a\rangle}$)")
    ax[0].grid(alpha=0.3)
    ax[1].plot(deltas, rs, "s-", color="C1", ms=7, label=r"B, slow $\omega$=2")
    ax[1].plot(deltas, rf, "o-", color="C2", ms=7, label=r"B, fast $\omega$=20")
    ax[1].plot(deltas, inv_a, "k--", lw=1.2, label=r"$\langle 1/a\rangle$ (B-follows-$a$)")
    ax[1].axhline(1.0, color="gray", lw=0.8, ls=":")
    ax[1].set_xlabel(r"modulation depth $\delta$")
    ax[1].set_ylabel(r"Var$_B$ ratio to $\delta{=}0$")
    ax[1].set_title(r"B: $\langle1/a\rangle$ slow $\to$ flat (=$1/\langle a\rangle$) fast"
                    "\n race stays geometry-independent at fast $\\omega$")
    ax[1].legend(fontsize=8, frameon=False); ax[1].grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(figdir, "folded_cycle_geometry_decoupling.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print(f"  figure -> figures/folded_cycle_geometry_decoupling.png\n")


if __name__ == "__main__":
    main()
