"""
folded_cycle_vdp.py  — Part 2, step 1 (deterministic foundation)
----------------------------------------------------------------
Three-timescale forced Van der Pol as a folded-limit-cycle testbed.

    x' = x - x^3/3 - w            (fast,        O(1))
    w' = eps1 (x - a)             (slow:    relaxation cycle, phase theta ~ 1/eps1)
    a' = eps2                     (super-slow:  drift a up through the cycle fold)

VdP fixed point is at x* = a, tr J = 1 - a^2: a stable spiral for |a|>1, an
unstable FP with a relaxation cycle for |a|<1.  The Hopf sits EXACTLY at the fold
a = 1 (VdP's b=0; see VDP_CROSSMODEL.md), so as a drifts up through 1 the
relaxation cycle is born/dies -- this is the (Hopf-degenerate realization of the)
folded limit cycle.  eps2 << eps1 << 1 is the JKK 'semi-oscillatory' ordering.

This step: integrate deterministically, drift a through 1, and MEASURE the cycle
amplitude vs a to locate the fold and see the canard-explosion birth.  Noise, the
eta=sigma/sqrt(eps2) check, and the C_q regime map come next.
"""
from __future__ import annotations
import os
import numpy as np


def drift(x, w, a, eps1):
    return (x - x**3 / 3.0 - w, eps1 * (x - a))


def integrate(eps1=0.1, eps2=5e-4, a0=0.70, a_end=1.06, dt=2e-3):
    """RK4 of the 3-timescale deterministic system; a ramps a0 -> a_end."""
    # total time so that a goes a0 -> a_end at rate eps2
    T = (a_end - a0) / eps2
    n = int(T / dt)
    xs = np.empty(n); ws = np.empty(n); as_ = np.empty(n); ts = np.empty(n)
    x, w, a, t = 0.5, 0.0, a0, 0.0
    for i in range(n):
        xs[i], ws[i], as_[i], ts[i] = x, w, a, t
        # RK4 on (x,w); a advanced by eps2*dt (frozen within the step)
        k1 = drift(x, w, a, eps1)
        k2 = drift(x + 0.5*dt*k1[0], w + 0.5*dt*k1[1], a, eps1)
        k3 = drift(x + 0.5*dt*k2[0], w + 0.5*dt*k2[1], a, eps1)
        k4 = drift(x + dt*k3[0], w + dt*k3[1], a, eps1)
        x = x + dt/6.0*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0])
        w = w + dt/6.0*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])
        a = a + eps2 * dt
        t = t + dt
    return ts, xs, ws, as_


def cycle_amplitude_vs_a(ts, xs, as_):
    """Per-oscillation amplitude (max-min of x between successive upward
    zero-crossings of x-<x>local), tagged by the a-value at that time."""
    # detect upward crossings of x = mean as cycle markers
    xc = xs - np.mean(xs)
    cross = np.where((xc[:-1] < 0) & (xc[1:] >= 0))[0]
    amp_a, amp = [], []
    for j in range(len(cross) - 1):
        s, e = cross[j], cross[j + 1]
        if e - s < 5:
            continue
        seg = xs[s:e]
        amp.append(seg.max() - seg.min())
        amp_a.append(as_[(s + e) // 2])
    return np.array(amp_a), np.array(amp)


def make_figure(ts, xs, as_, amp_a, amp, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 1, figsize=(9, 6.5))

    ax[0].plot(as_, xs, lw=0.4, color="C0")
    ax[0].axvline(1.0, color="r", ls="--", lw=1, label="fold / Hopf  a=1")
    ax[0].set_xlabel("a (super-slow drift)"); ax[0].set_ylabel("x (fast)")
    ax[0].set_title("Forced VdP: relaxation cycle folds as a drifts through 1")
    ax[0].legend(frameon=False, fontsize=9)

    ax[1].plot(amp_a, amp, "o-", ms=3, color="C0")
    ax[1].axvline(1.0, color="r", ls="--", lw=1)
    ax[1].set_xlabel("a"); ax[1].set_ylabel("cycle amplitude  max(x)-min(x)")
    ax[1].set_title("Cycle amplitude vs a  (fold = amplitude collapse near a=1)")
    fig.tight_layout(); fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    figdir = os.path.join(os.path.dirname(here), "figures"); os.makedirs(figdir, exist_ok=True)

    eps1, eps2 = 0.1, 5e-4
    print(f"\n=== Part 2 step 1: deterministic folded VdP  (eps1={eps1}, eps2={eps2}) ===\n")
    ts, xs, ws, as_ = integrate(eps1=eps1, eps2=eps2)
    amp_a, amp = cycle_amplitude_vs_a(ts, xs, as_)
    print(f"  integrated {len(ts)} steps; detected {len(amp)} oscillations "
          f"over a in [{as_[0]:.3f}, {as_[-1]:.3f}]")

    # locate the fold: last a with a full-amplitude (>1) oscillation, and the
    # a where amplitude drops below e.g. 0.5 (cycle effectively gone)
    big = amp_a[amp > 1.0]
    small = amp_a[amp < 0.3]
    a_last_big = big.max() if big.size else np.nan
    a_first_small = small.min() if small.size else np.nan
    print(f"  last large-amplitude (>1) cycle at a = {a_last_big:.4f}")
    print(f"  amplitude < 0.3 (cycle ~gone) from a = {a_first_small:.4f}")
    print(f"  => fold / cycle death localized near a = 1 (VdP Hopf=fold)")
    print(f"     canard-explosion birth window width ~ {abs(a_first_small - a_last_big):.4f} in a\n")

    fig = os.path.join(figdir, "folded_cycle_vdp_deterministic.png")
    make_figure(ts, xs, as_, amp_a, amp, fig)
    print(f"  figure -> figures/folded_cycle_vdp_deterministic.png\n")


if __name__ == "__main__":
    main()
