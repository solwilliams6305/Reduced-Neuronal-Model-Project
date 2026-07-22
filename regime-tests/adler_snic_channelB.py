"""
adler_snic_channelB.py — Channel B at a SNIC (finite-amplitude degeneration)
----------------------------------------------------------------------------
The autonomous VdP attempt failed because its cycle dies at a HOPF (amplitude->0
=> resonator).  A SNIC (saddle-node on invariant circle) instead kills the cycle
at FINITE amplitude (the PERIOD diverges, the trajectory lingers in a bottleneck),
so there is no resonator.  Cleanest controlled SNIC oscillator = the Adler eqn
(phase on a circle, fixed amplitude):

    dphi = ( omega - A sin(phi) ) dt + sigma dW,   A' = eps2   (drift through SNIC)

SNIC exactly at A = omega: for A<omega phi rotates (period 2pi/sqrt(omega^2-A^2));
as A->omega the bottleneck (near phi=pi/2) makes the passage diverge and become
maximally noise-sensitive; for A>omega two fixed points (excitable, noise-induced
rotations continue).  Amplitude is fixed -> NO resonator, and rotation-counting does
NOT saturate (noise-induced rotations persist past the SNIC).

TEST: drift A:0.5->1.3 through omega=1; count rotations (phi through 2pi) per
trajectory; Var(rotation count) = phase decoherence.  sigma_*^B (Var~0.25) should
be cleanly measurable and scale ~ sqrt(eps2) (the Channel-B drift scaling) --
grounding Channel B's existence + eps2-scaling on a finite-amplitude autonomous
folded cycle, isolating the Hopf as the VdP obstruction.
"""
from __future__ import annotations
import os
import numpy as np

OMEGA = 1.0


def simulate(eps2, sigma, A0=0.5, A_end=1.3, N=200, dt=1e-2, seed=0):
    """drift A0->A_end through the SNIC; return Var(rotation count)."""
    rng = np.random.default_rng(seed)
    phi = np.zeros(N)
    A = A0
    sdt = np.sqrt(dt)
    n = int((A_end - A0) / eps2 / dt)
    for _ in range(n):
        phi = phi + (OMEGA - A * np.sin(phi)) * dt + sigma * sdt * rng.standard_normal(N)
        A = A + eps2 * dt
    return float(np.var(phi / (2 * np.pi)))      # variance of (unwrapped) rotation count


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    figdir = os.path.join(os.path.dirname(here), "figures")
    os.makedirs(figdir, exist_ok=True)

    eps2_list = [2e-3, 4e-3, 8e-3]
    sigmas = np.array([0.02, 0.04, 0.07, 0.12, 0.20, 0.32, 0.50])
    N = 200
    thr = 0.25

    print("\n=== Channel B at a SNIC (Adler eqn, finite amplitude, no resonator) ===")
    print(f"  omega={OMEGA}; drift A:0.5->1.3 through SNIC at A=1; Var(rotation count).\n")
    print(f"  {'eps2':>8s}   Var(rot) vs sigma ->")
    sB = []
    for e2 in eps2_list:
        row = np.array([simulate(e2, s, N=N, seed=int(1e4 * e2) + k) for k, s in enumerate(sigmas)])
        s_star = np.nan
        for i in range(len(sigmas) - 1):
            if row[i] < thr <= row[i + 1]:
                f = (thr - row[i]) / (row[i + 1] - row[i])
                s_star = sigmas[i] * (sigmas[i + 1] / sigmas[i]) ** f; break
        sB.append(s_star)
        print(f"  {e2:8.0e}   [" + ",".join(f"{v:.2f}" for v in row) + f"]   sigma_*^B={s_star:.3f}")

    e2a = np.array(eps2_list); sBa = np.array(sB); good = ~np.isnan(sBa)
    slope = np.polyfit(np.log(e2a[good]), np.log(sBa[good]), 1)[0] if good.sum() >= 2 else np.nan
    print(f"\n  sigma_*^B = {np.round(sBa,3)}")
    print(f"  sigma_*^B/sqrt(eps2) = {np.round(sBa/np.sqrt(e2a),2)}  (flat => ~sqrt(eps2))")
    print(f"  fitted slope d ln sigma_*^B / d ln eps2 = {slope:.3f}   (predict ~0.5)")
    if good.all() and abs(slope - 0.5) < 0.2:
        print(f"\n  => CLEAN: sigma_*^B measurable and ~sqrt(eps2). Channel B is grounded on a")
        print(f"     finite-amplitude autonomous SNIC cycle (no resonator, no saturation).")
        print(f"     The VdP obstruction was the Hopf (amplitude->0), now isolated.\n")
    else:
        print(f"\n  => partial; see slope/values above.\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    for e2 in eps2_list:
        row = np.array([simulate(e2, s, N=N, seed=int(2e4 * e2) + k) for k, s in enumerate(sigmas)])
        ax[0].plot(sigmas, row, "o-", ms=4, label=fr"$\epsilon_2$={e2:.0e}")
    ax[0].axhline(thr, color="k", ls=":", lw=0.8)
    ax[0].set_xlabel(r"$\sigma$"); ax[0].set_ylabel(r"Var(rotation count) = phase decoherence")
    ax[0].set_title("Channel B at a SNIC: clean phase decoherence (no saturation)")
    ax[0].legend(fontsize=8, frameon=False); ax[0].grid(alpha=0.3)
    if good.sum() >= 2:
        ax[1].loglog(e2a[good], sBa[good], "o-", color="C0", ms=7, label=f"measured (slope {slope:.2f})")
        ax[1].loglog(e2a[good], sBa[good][0] * (e2a[good] / e2a[good][0]) ** 0.5, "k--", lw=1, label="slope 1/2")
        ax[1].legend(fontsize=9, frameon=False)
    ax[1].set_xlabel(r"$\epsilon_2$ (drift through SNIC)"); ax[1].set_ylabel(r"$\sigma_*^B$")
    ax[1].set_title(r"$\sigma_*^B \sim \sqrt{\epsilon_2}$ — Channel B grounded (SNIC)")
    ax[1].grid(alpha=0.3, which="both")
    fig.tight_layout()
    out = os.path.join(figdir, "adler_snic_channelB.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print(f"  figure -> figures/adler_snic_channelB.png\n")


if __name__ == "__main__":
    main()
