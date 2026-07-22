"""
kp_cusp_explore.py — locate the MMO/cusp regime of the Kristiansen–Pedersen coupled FHN.
========================================================================================

Exact KP model (arXiv:2202.12027), gap-junction coupling on the fast variable:
    v_i' = -v_i^3 + 3 v_i - w_i + g (v_j - v_i)
    w_i' = ε (v_i - c)
Cubic fast nullcline w = -v^3 + 3v (folds at v = ±1, w = ∓2). Single unit (g=0) is a
relaxation oscillator for |c| < 1. KP prove that symmetric (repulsive) coupling makes the
pair do MMOs organised by a CUSP of the critical manifold.

Here we just SCAN (c, g) at small ε and detect mixed-mode oscillations (small-amplitude
SAOs + large relaxation spikes) in v_1, to find the cusp/canard regime where a noise-induced
peel-off law can be measured. Vectorised RK4 over the grid.

Output: figures/kp_cusp_explore.png + the located MMO cell printed.
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
EPS = 0.015


def scan(c_vals, g_vals, T=2200.0, dt=0.02, warmup=600.0, stride=3):
    nc, ng = len(c_vals), len(g_vals)
    C = np.repeat(c_vals, ng)[:, None]            # (R,1)
    G = np.tile(g_vals, nc)[:, None]
    R = nc * ng
    v = np.tile([0.2, -0.2], (R, 1)).astype(float)
    w = np.zeros((R, 2))

    def drift(v, w):
        coup = G * (v[:, ::-1] - v)
        return (-v**3 + 3*v - w + coup, EPS * (v - C))

    def rk4(v, w):
        k1 = drift(v, w); k2 = drift(v+.5*dt*k1[0], w+.5*dt*k1[1])
        k3 = drift(v+.5*dt*k2[0], w+.5*dt*k2[1]); k4 = drift(v+dt*k3[0], w+dt*k3[1])
        return (v+dt*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6, w+dt*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6)

    for _ in range(int(warmup/dt)):
        v, w = rk4(v, w)
    n = int(T/dt); nrec = n//stride
    V = np.empty((nrec, R)); ri = 0
    for k in range(n):
        v, w = rk4(v, w)
        if k % stride == 0 and ri < nrec:
            V[ri] = v[:, 0]; ri += 1
    return V[:ri], C.ravel(), G.ravel(), (nc, ng)


def count_modes(v):
    dv = np.diff(v); pk = np.where((dv[:-1] > 0) & (dv[1:] <= 0))[0] + 1
    p = v[pk]
    nL = int(np.sum(p > 1.5))                      # large relaxation spikes
    nS = int(np.sum((p > 0.4) & (p <= 1.5)))       # small-amplitude oscillations
    return nL, nS


def main():
    t0 = time.time()
    c_vals = np.linspace(0.80, 1.02, 12)
    g_vals = np.linspace(-0.7, 0.3, 13)
    V, C, G, (nc, ng) = scan(c_vals, g_vals)
    R = V.shape[1]
    nL = np.zeros(R, int); nS = np.zeros(R, int)
    for r in range(R):
        nL[r], nS[r] = count_modes(V[:, r])
    mmo = (nL >= 2) & (nS >= 1)
    ratio = np.where(nL > 0, nS / np.maximum(nL, 1), 0.0)

    print("=" * 70)
    print(f"KP coupled FHN — MMO scan (ε={EPS})")
    print("=" * 70)
    print(f"  MMO cells (nL≥2 & nS≥1): {mmo.sum()} / {R}")
    if mmo.any():
        # pick the MMO cell with the most SAOs (deepest canard)
        idx = np.where(mmo)[0]
        best = idx[np.argmax(nS[idx])]
        print(f"  richest MMO: c={C[best]:.3f}, g={G[best]:.3f}  (nL={nL[best]}, nS={nS[best]})")
    else:
        print("  no MMO cells found — adjust ε / ranges")
    print(f"  SAO-ratio range: {ratio.min():.2f} … {ratio.max():.2f}")

    REG = ratio.reshape(nc, ng)
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    ext = [g_vals[0], g_vals[-1], c_vals[0], c_vals[-1]]
    im = ax[0].imshow(REG, origin="lower", aspect="auto", extent=ext, cmap="magma")
    fig.colorbar(im, ax=ax[0], label="SAO per large spike (nS/nL)")
    ax[0].set_xlabel("coupling g"); ax[0].set_ylabel("control c")
    ax[0].set_title(f"(A) MMO map of KP coupled FHN (ε={EPS})")
    if mmo.any():
        ax[0].plot(G[best], C[best], "c*", ms=16, markeredgecolor="k")

    if mmo.any():
        # show the richest-MMO trajectory
        tg = np.linspace(0, 1, V.shape[0])
        ax[1].plot(tg, V[:, best], lw=0.7, color="#1f3b73")
        ax[1].axhline(1.5, color="grey", lw=0.5, ls=":")
        ax[1].set_title(f"(B) MMO trace at c={C[best]:.3f}, g={G[best]:.3f}")
        ax[1].set_xlabel("time (norm.)"); ax[1].set_ylabel("v₁")
    fig.suptitle("Kristiansen–Pedersen coupled FHN — locating the cusp/MMO regime", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fp = os.path.join(FIG, "kp_cusp_explore.png")
    fig.savefig(fp, dpi=130)
    print(f"\n  Figure: {fp}")
    print(f"  TOTAL {time.time()-t0:.1f}s — DONE")


if __name__ == "__main__":
    main()
