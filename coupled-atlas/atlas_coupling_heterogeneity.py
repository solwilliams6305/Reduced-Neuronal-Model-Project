"""
atlas_coupling_heterogeneity.py — first deterministic regime atlas.
===================================================================

Maps the (coupling strength g) x (frequency heterogeneity ΔI) plane for two
electrically (linear gap-junction) coupled FitzHugh–Rinzel units in the MMO
band (c = -0.82). Units detuned via the drive I:
    I_1 = I0 - ΔI/2 ,  I_2 = I0 + ΔI/2 .

Classification uses two robust, unbiased diagnostics:
  • ρ  = zero-lag Pearson correlation of v1,v2  (+1 in-phase, -1 anti-phase, 0 incoherent)
  • Δf = |f1 - f2| / mean(f)  spike-rate mismatch  (→0 means 1:1 frequency-locked)
plus the Golomb–Hansel χ and oscillation amplitude.

Regimes: death / in-phase / anti-phase / locked-other / unlocked(drifting).
Panels: (1) regime map  (2) synchrony χ  (3) rate mismatch Δf (the Arnold tongue).

All grid cells stepped together in one vectorized RK4 pass.
Outputs: figures/atlas_coupling_heterogeneity.png + atlas_coupling_heterogeneity.npz
"""
from __future__ import annotations
import os, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

from coupled_fhr import sync_chi, large_spike_times

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

A, B, EPS, DELTA = 0.7, 0.8, 0.08, 0.2
C0, I0 = -0.82, 0.30


def integrate_grid(C, Ivec, g, T=2000.0, dt=0.04, warmup=900.0, stride=4):
    R = C.shape[0]
    V = np.tile([0.6, -0.6], (R, 1)).astype(float)
    W = np.zeros((R, 2)); Y = np.zeros((R, 2))

    def drift(V, W, Y):
        coup = g * (V[:, ::-1] - V)
        return (V - V**3 / 3.0 - W + Y + Ivec + coup,
                EPS * (V + A - B * W),
                EPS * DELTA * (C - V))

    def rk4(V, W, Y):
        k1 = drift(V, W, Y)
        k2 = drift(V + .5*dt*k1[0], W + .5*dt*k1[1], Y + .5*dt*k1[2])
        k3 = drift(V + .5*dt*k2[0], W + .5*dt*k2[1], Y + .5*dt*k2[2])
        k4 = drift(V + dt*k3[0], W + dt*k3[1], Y + dt*k3[2])
        return (V + dt*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6,
                W + dt*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6,
                Y + dt*(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6)

    for _ in range(int(warmup / dt)):
        V, W, Y = rk4(V, W, Y)
    n = int(T / dt); nrec = n // stride
    ts = np.empty(nrec); Vr = np.empty((nrec, R, 2)); ri = 0
    for k in range(n):
        V, W, Y = rk4(V, W, Y)
        if k % stride == 0 and ri < nrec:
            ts[ri] = k * dt; Vr[ri] = V; ri += 1
    return ts[:ri], Vr[:ri]


def _rate(t, v):
    tt = large_spike_times(t, v)
    return tt.size / (t[-1] - t[0]) if tt.size >= 2 else 0.0


def classify(ts, Vr, death_amp=0.8, lock_tol=0.06):
    R = Vr.shape[1]
    amp = Vr.max(0) - Vr.min(0)
    reg = np.zeros(R, int); chi = np.empty(R); rho = np.empty(R); df = np.empty(R)
    for r in range(R):
        v = Vr[:, r, :]
        chi[r] = sync_chi(v)
        v0 = v - v.mean(0)
        s1, s2 = v0[:, 0].std(), v0[:, 1].std()
        rho[r] = float(np.mean(v0[:, 0] * v0[:, 1]) / (s1 * s2)) if s1 > 0 and s2 > 0 else np.nan
        f1, f2 = _rate(ts, v[:, 0]), _rate(ts, v[:, 1])
        fm = 0.5 * (f1 + f2)
        df[r] = abs(f1 - f2) / fm if fm > 0 else np.nan
        locked = (not np.isnan(df[r])) and (df[r] < lock_tol)
        if amp[r].max() < death_amp:
            reg[r] = 0                               # death / quiescent
        elif locked and rho[r] > 0.3:
            reg[r] = 1                               # in-phase
        elif locked and rho[r] < -0.3:
            reg[r] = 2                               # anti-phase
        elif locked:
            reg[r] = 4                               # locked, intermediate phase
        else:
            reg[r] = 3                               # unlocked / drifting
    return reg, chi, rho, df


def main():
    t0 = time.time()
    ng, nh = 22, 22
    g_vals = np.linspace(0.0, 0.25, ng)
    dc_vals = np.linspace(0.0, 0.14, nh)
    GG, HH = np.meshgrid(g_vals, dc_vals, indexing="xy")
    g_flat = GG.ravel(); dc_flat = HH.ravel(); R = g_flat.size

    C = np.empty((R, 2))
    C[:, 0] = C0 - dc_flat / 2.0           # heterogeneity via c — real period detuning
    C[:, 1] = C0 + dc_flat / 2.0
    Ivec = np.full((R, 2), I0)

    ts, Vr = integrate_grid(C, Ivec, g_flat[:, None])
    reg, chi, rho, df = classify(ts, Vr)
    REG = reg.reshape(nh, ng); CHI = chi.reshape(nh, ng); DF = df.reshape(nh, ng)

    np.savez(os.path.join(HERE, "atlas_coupling_heterogeneity.npz"),
             g=g_vals, dc=dc_vals, reg=REG, chi=CHI, rho=rho.reshape(nh, ng), df=DF)

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    ext = [g_vals[0], g_vals[-1], dc_vals[0], dc_vals[-1]]

    labels = ["death", "in-phase", "anti-phase", "unlocked", "locked-other"]
    colors = ["#cfcfcf", "#2c7fb8", "#d95f0e", "#f7e463", "#7fbf7b"]
    order = [0, 1, 2, 3, 4]
    cmap = ListedColormap([colors[i] for i in order])
    norm = BoundaryNorm([-.5, .5, 1.5, 2.5, 3.5, 4.5], cmap.N)
    im0 = ax[0].imshow(REG, origin="lower", aspect="auto", extent=ext, cmap=cmap, norm=norm)
    cb = fig.colorbar(im0, ax=ax[0], ticks=order)
    cb.ax.set_yticklabels([labels[i] for i in order])
    ax[0].set_title("(1) phase-relationship regime")

    im1 = ax[1].imshow(CHI, origin="lower", aspect="auto", extent=ext,
                       cmap="viridis", vmin=0, vmax=1)
    fig.colorbar(im1, ax=ax[1]); ax[1].set_title("(2) synchrony χ")

    im2 = ax[2].imshow(DF, origin="lower", aspect="auto", extent=ext,
                       cmap="magma_r", vmin=0, vmax=0.5)
    fig.colorbar(im2, ax=ax[2])
    ax[2].set_title("(3) rate mismatch Δf  (dark = 1:1 locked → Arnold tongue)")

    for a_ in ax:
        a_.set_xlabel("coupling strength  g")
        a_.set_ylabel("heterogeneity  Δc")

    fig.suptitle(f"Coupled FHR atlas — c₀={C0}, I₀={I0}:  coupling × period heterogeneity (Δc) "
                 f"(linear gap junction)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fp = os.path.join(FIG, "atlas_coupling_heterogeneity.png")
    fig.savefig(fp, dpi=130)

    counts = {labels[k]: int((reg == k).sum()) for k in range(5)}
    print("regime cell counts:", counts, flush=True)
    print(f"χ range [{np.nanmin(chi):.2f}, {np.nanmax(chi):.2f}];  "
          f"ρ range [{np.nanmin(rho):.2f}, {np.nanmax(rho):.2f}]", flush=True)
    print(f"Figure: {fp}", flush=True)
    print(f"TOTAL {time.time()-t0:.1f}s — DONE", flush=True)


if __name__ == "__main__":
    main()
