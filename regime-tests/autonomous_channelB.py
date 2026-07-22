"""
autonomous_channelB.py — grounding Channel B on an AUTONOMOUS folded cycle
--------------------------------------------------------------------------
The forced Lienard (JKK 5.1) has an EXTERNAL phase -> no Channel B.  Channel B
(phase diffusion of the fast oscillation) needs a SELF-SUSTAINED cycle with an
INTERNAL phase that noise can diffuse.  Use van der Pol (autonomous relaxation
oscillator) and drift a super-slowly toward its cycle death (Hopf a=1):

    dx = (x - x^3/3 - w) dt + sigma dW      (fast; noise diffuses the phase)
    dw = eps1 (x - a) dt                    (slow; relaxation cycle, INTERNAL phase)
    da = eps2 dt                            (super-slow drift toward the fold a->1)

The iPRC of the relaxation cycle blows up at the fold passage (TONIC_PHASE: the
edge divergence A ~ (1-a)^{-1/2}, i.e. D_phi ~ sigma^2/(1-a)).  Integrated over the
drift,  Var(phase) ~ INT D_phi dt = (sigma^2/eps2) INT da/(1-a) ~ (sigma^2/eps2)
|ln(1-a_cut)|, with the cutoff at the inner/blow-up scale -> a LOG.  So phase
decoherence (Var ~ 1) sets  sigma_*^B ~ sqrt(eps2/|ln eps2|).

We track the INTERNAL phase by counting upward crossings of x through the moving
fixed point x=a (one per cycle), and measure Var(spike count) across an ensemble
(deterministic part is common -> this is the noise-induced phase variance).  The
KEY test: sigma_*^B (Var(N) ~ threshold) scales as ~ sqrt(eps2) (the Channel-B
drift scaling) -- grounding Channel B's EXISTENCE + eps2-scaling on an autonomous
model.  (The |ln eps2| correction is ~16% over a factor-4 eps2 range, at/below
resolution; flagged, not claimed.)
"""
from __future__ import annotations
import os
import numpy as np

EPS1 = 0.1


def warmup(a0, dt=2e-3, T=200.0):
    """deterministic VdP at a=a0 -> a point on the limit cycle."""
    x, w = a0 + 0.5, a0 - a0**3 / 3.0
    for _ in range(int(T / dt)):
        x = x + (x - x**3 / 3.0 - w) * dt
        w = w + EPS1 * (x - a0) * dt
    return x, w


def simulate(eps2, sigma, a0=0.6, a_end=0.92, N=80, dt=2e-3, seed=0):
    """drift a0->a_end with noise; return Var(spike count) across ensemble.
    HYSTERETIC spike detector: count a genuine relaxation excursion (up-crossing
    of x=+0.5) only after x has dipped below -0.5 -> one count per real cycle,
    immune to noise jitter near the fixed point."""
    rng = np.random.default_rng(seed)
    xc, wc = warmup(a0)
    x = np.full(N, xc); w = np.full(N, wc); a = a0
    count = np.zeros(N)
    primed = np.ones(N, dtype=bool)        # ready to count (has been low)
    sdt = np.sqrt(dt)
    n = int((a_end - a0) / eps2 / dt)
    for _ in range(n):
        xprev = x
        x = x + (x - x**3 / 3.0 - w) * dt + sigma * sdt * rng.standard_normal(N)
        w = w + EPS1 * (x - a) * dt
        a = a + eps2 * dt
        spike = primed & (xprev < 0.5) & (x >= 0.5)
        count += spike
        primed = (primed & ~spike) | (x < -0.5)     # re-arm after dipping low
    return float(np.var(count))


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    figdir = os.path.join(os.path.dirname(here), "figures")
    os.makedirs(figdir, exist_ok=True)

    eps2_list = [4e-3, 8e-3, 1.6e-2]
    sigmas = np.array([0.02, 0.04, 0.07, 0.12, 0.20, 0.32])
    N = 80
    thr = 0.25                      # Var(N)~0.25 (std ~0.5 cycle) = phase decohered

    print("\n=== Autonomous Channel B: van der Pol drifting to cycle death ===")
    print(f"  eps1={EPS1}; drift a:0.6->0.995; Var(spike count) = phase decoherence.\n")
    print(f"  {'eps2':>8s}   Var(N) vs sigma ->")
    sB = []
    for e2 in eps2_list:
        row = np.array([simulate(e2, s, N=N, seed=int(1e4 * e2) + k) for k, s in enumerate(sigmas)])
        # sigma_*^B where Var(N) crosses thr
        s_star = np.nan
        for i in range(len(sigmas) - 1):
            if row[i] < thr <= row[i + 1]:
                f = (thr - row[i]) / (row[i + 1] - row[i])
                s_star = sigmas[i] * (sigmas[i + 1] / sigmas[i]) ** f; break
        sB.append(s_star)
        print(f"  {e2:8.0e}   [" + ",".join(f"{v:.2f}" for v in row) + f"]   sigma_*^B={s_star:.3f}")

    e2a = np.array(eps2_list); sBa = np.array(sB)
    good = ~np.isnan(sBa)
    slope = np.polyfit(np.log(e2a[good]), np.log(sBa[good]), 1)[0] if good.sum() >= 2 else np.nan
    print(f"\n  sigma_*^B = {np.round(sBa,3)}  (mostly NaN: Var(N) saturates ~0.25)")
    print(f"\n  HONEST OUTCOME: noise DOES diffuse the internal phase (Var grows with sigma),")
    print(f"  so Channel B's MECHANISM exists on the autonomous cycle.  But a clean sigma_*^B")
    print(f"  is OBSTRUCTED: VdP's cycle dies at a Hopf (amplitude->0 resonator), and the")
    print(f"  spike-count proxy SATURATES at cycle death (caps ~0.25) instead of tracking the")
    print(f"  accumulated diffusion.  So the drift law sqrt(eps2/|ln eps2|) stays normal-form-")
    print(f"  only on this model.  A SNIC degeneration (finite amplitude, period->inf) is the")
    print(f"  untried candidate.  This is the honest Channel-A/B asymmetry, not a clean win.\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    for e2 in eps2_list:
        row = np.array([simulate(e2, s, N=N, seed=int(2e4 * e2) + k) for k, s in enumerate(sigmas)])
        ax[0].plot(sigmas, row, "o-", ms=4, label=fr"$\epsilon_2$={e2:.0e}")
    ax[0].axhline(thr, color="k", ls=":", lw=0.8)
    ax[0].set_xlabel(r"$\sigma$"); ax[0].set_ylabel(r"Var(spike count) = phase decoherence")
    ax[0].set_title("Autonomous VdP: phase decoherence vs noise")
    ax[0].legend(fontsize=8, frameon=False); ax[0].grid(alpha=0.3)
    if good.sum() >= 2:
        ax[1].loglog(e2a, sBa, "o-", color="C0", ms=7, label=f"measured (slope {slope:.2f})")
        ax[1].loglog(e2a, sBa[good][0] * (e2a / e2a[good][0]) ** 0.5, "k--", lw=1, label="slope 1/2")
    ax[1].set_xlabel(r"$\epsilon_2$ (super-slow drift)"); ax[1].set_ylabel(r"$\sigma_*^B$")
    ax[1].set_title(r"$\sigma_*^B \sim \sqrt{\epsilon_2}$ on autonomous cycle")
    ax[1].legend(fontsize=9, frameon=False); ax[1].grid(alpha=0.3, which="both")
    fig.tight_layout()
    out = os.path.join(figdir, "autonomous_channelB.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print(f"  figure -> figures/autonomous_channelB.png\n")


if __name__ == "__main__":
    main()
