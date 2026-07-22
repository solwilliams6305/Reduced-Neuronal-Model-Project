"""
coupled_peeloff_fhr.py — T1.1 physical cross-check in the REAL coupled FHR.
===========================================================================

The inner-model result (coupled_peeloff.py) showed that electrical coupling
correlates the two folded cycles' peel-off levels (corr 0 → 0.86). Here we check
the SAME phenomenon in the full noisy coupled FitzHugh–Rinzel model — no blow-up
chart, real spikes, real gap-junction current.

Setup: both units in the MMO band (c=-0.85), additive noise σ on v, electrical
coupling g(v_j - v_i). From a common post-spike SAO-phase initial condition, each
realization runs until each unit's FIRST large spike (v > v_thr) — that escape is
the physical peel-off. We record the peel-off TIME and the slow variable y at peel-off
(the physical analogue of the inner Y_node), and measure their cross-unit correlation
vs coupling g across an ensemble of independent-noise realizations.

Prediction (from the inner model): corr ≈ 0 at g=0 (independent noise), rising toward
1 as g increases.  Euler–Maruyama, vectorized over (g, realization).

Output: figures/coupled_peeloff_fhr.png + printed summary.
"""
from __future__ import annotations
import os, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

A, B, EPS, DELTA, I0 = 0.7, 0.8, 0.08, 0.2, 0.30
C = -0.85
VTHR = 1.0                       # large-spike (peel-off) threshold


def _f(v, w, y):
    return (v - v**3 / 3.0 - w + y + I0, EPS * (v + A - B * w), EPS * DELTA * (C - v))


def find_ic(dt=0.02):
    """Deterministic warmup → a common post-spike SAO-phase initial condition."""
    v, w, y = 0.5, 0.0, 0.0
    def rk4(v, w, y):
        k1 = _f(v, w, y); k2 = _f(v+.5*dt*k1[0], w+.5*dt*k1[1], y+.5*dt*k1[2])
        k3 = _f(v+.5*dt*k2[0], w+.5*dt*k2[1], y+.5*dt*k2[2]); k4 = _f(v+dt*k3[0], w+dt*k3[1], y+dt*k3[2])
        return (v+dt*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6,
                w+dt*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6,
                y+dt*(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6)
    for _ in range(int(300/dt)):
        v, w, y = rk4(v, w, y)
    spiked = False
    for _ in range(int(300/dt)):
        v, w, y = rk4(v, w, y)
        if v > VTHR:
            spiked = True
        if spiked and v < -0.3:
            return (v, w, y)
    return (v, w, y)


def ensemble_first_peeloff(g_arr, sigma, ic, N=2000, T=600.0, dt=0.025, seed=0):
    G = len(g_arr)
    g = np.asarray(g_arr, float)[:, None, None]
    v = np.full((G, N, 2), ic[0]); w = np.full((G, N, 2), ic[1]); y = np.full((G, N, 2), ic[2])
    done = np.zeros((G, N, 2), bool)
    tpk = np.full((G, N, 2), np.nan); ypk = np.full((G, N, 2), np.nan)
    rng = np.random.default_rng(seed); sdt = np.sqrt(dt)
    n = int(T / dt)
    for step in range(n):
        t = step * dt
        coup = g * (v[:, :, ::-1] - v)
        v = v + (v - v**3 / 3.0 - w + y + I0 + coup) * dt + sigma * sdt * rng.standard_normal((G, N, 2))
        w = w + EPS * (v + A - B * w) * dt
        y = y + EPS * DELTA * (C - v) * dt
        new = (~done) & (v > VTHR)
        if new.any():
            tpk[new] = t; ypk[new] = y[new]; done |= new
        if done.all():
            break
    return tpk, ypk


def main():
    t0 = time.time()
    ic = find_ic()
    print("=" * 72)
    print(f"T1.1 physical cross-check — coupled FHR (c={C}), peel-off = first large spike")
    print("=" * 72)
    print(f"  common IC (post-spike SAO phase): v={ic[0]:.3f} w={ic[1]:.3f} y={ic[2]:.3f}")

    sigma = 0.03
    g_arr = np.array([0.0, 0.05, 0.1, 0.2, 0.4])
    tpk, ypk = ensemble_first_peeloff(g_arr, sigma, ic, seed=5)

    print(f"  σ={sigma}, ensemble N={tpk.shape[1]}")
    print(f"  {'g':>6} | {'corr(t1,t2)':>12} {'corr(y1,y2)':>12} {'esc.frac':>9} {'std(t)':>8}")
    ct, cy = [], []
    for gi, gg in enumerate(g_arr):
        t1, t2 = tpk[gi, :, 0], tpk[gi, :, 1]
        y1, y2 = ypk[gi, :, 0], ypk[gi, :, 1]
        ok = np.isfinite(t1) & np.isfinite(t2)
        rt = float(np.corrcoef(t1[ok], t2[ok])[0, 1])
        ry = float(np.corrcoef(y1[ok], y2[ok])[0, 1])
        ct.append(rt); cy.append(ry)
        print(f"  {gg:6.3f} | {rt:12.3f} {ry:12.3f} {ok.mean():9.2f} {np.std(t1[ok]):8.2f}")
    ct = np.array(ct); cy = np.array(cy)

    rises = ct[-1] > ct[0] + 0.2 and cy[-1] > cy[0] + 0.2
    indep0 = abs(ct[0]) < 0.15
    print(f"\n  g=0 independent: corr(t1,t2)={ct[0]:+.3f} {'PASS' if indep0 else 'CHECK'}")
    print(f"  coupling correlates physical peel-offs: corr_t {ct[0]:+.2f}→{ct[-1]:+.2f}, "
          f"corr_y {cy[0]:+.2f}→{cy[-1]:+.2f}  {'PASS' if rises else 'CHECK'}")
    print(f"  ⇒ inner-model prediction confirmed in the full coupled FHR." if rises else
          "  ⇒ revisit σ / threshold.")

    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.5))
    ax[0].plot(g_arr, ct, "o-", color="#b3402b", label="corr(t¹,t²)  peel-off time")
    ax[0].plot(g_arr, cy, "s--", color="#1f3b73", label="corr(y¹,y²)  peel-off level")
    ax[0].set_ylim(-0.1, 1.02); ax[0].axhline(0, color="grey", lw=0.5, ls=":")
    ax[0].set_xlabel("electrical coupling g"); ax[0].set_ylabel("cross-unit correlation")
    ax[0].set_title("(A) coupling correlates physical peel-offs (coupled FHR)")
    ax[0].legend(fontsize=9, frameon=False)

    for gi, col, lab in [(0, "#888780", f"g=0 (corr {ct[0]:+.2f})"),
                         (len(g_arr) - 1, "#1f3b73", f"g={g_arr[-1]} (corr {ct[-1]:+.2f})")]:
        t1, t2 = tpk[gi, :, 0], tpk[gi, :, 1]
        ok = np.isfinite(t1) & np.isfinite(t2)
        ax[1].plot(t1[ok], t2[ok], ".", ms=2.5, alpha=0.3, color=col, label=lab)
    ax[1].set_xlabel("unit 1 peel-off time t¹"); ax[1].set_ylabel("unit 2 peel-off time t²")
    ax[1].set_title("(B) joint peel-off timing: g=0 vs strong g")
    ax[1].legend(fontsize=9, frameon=False)

    fig.suptitle("T1.1 cross-check — coupling correlates peel-offs in the real coupled FHR",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fp = os.path.join(FIG, "coupled_peeloff_fhr.png")
    fig.savefig(fp, dpi=140, bbox_inches="tight")
    print(f"\n  Figure: {fp}")
    print(f"  TOTAL {time.time()-t0:.1f}s — DONE")


if __name__ == "__main__":
    main()
