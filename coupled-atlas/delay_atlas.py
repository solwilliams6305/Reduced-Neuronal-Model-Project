"""
delay_atlas.py — T2.4: delay in pure electrical (gap-junction) coupling.
========================================================================

Two IDENTICAL FitzHugh–Rinzel units, electrical coupling WITH delay:
    v_i' += g ( v_j(t - τ) - v_i(t) ) .

We map the (coupling g, delay τ) plane and classify the asymptotic state, looking
for the textbook delay phenomena that are well established for limit-cycle / neural
oscillators but comparatively unmapped for slow-fast pairs:
  • in-phase ↔ anti-phase switching as τ grows (phase-flip; Crook et al. 1997,
    Prasad et al. 2008),
  • delay-induced oscillation/amplitude DEATH in islands of (g, τ)
    (Reddy–Sen–Johnston 1998).

Classification (deterministic): amplitude (death if small), and the zero-lag voltage
correlation ρ (+1 in-phase, −1 anti-phase). All (g,τ) cells stepped together; delay
handled by a shared circular history buffer with per-replica lag (method of steps,
delayed value frozen across the RK4 substeps).

Output: figures/delay_atlas.png + printed summary.
"""
from __future__ import annotations
import os, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

A, B, EPS, DELTA, I0 = 0.7, 0.8, 0.08, 0.2, 0.30
C = -0.72                                   # clean relaxation oscillator (period ~54)


def delay_grid(g_flat, tau_flat, T=500.0, warmup=350.0, dt=0.04, stride=5):
    R = len(g_flat)
    lag = np.maximum(1, np.round(tau_flat / dt).astype(int))   # (R,)
    buflen = int(lag.max()) + 1
    g = g_flat[:, None]
    V = np.tile([0.6, -0.6], (R, 1)).astype(float)
    W = np.zeros((R, 2)); Y = np.zeros((R, 2))
    Vbuf = np.repeat(V[None, :, :], buflen, axis=0)            # (buflen,R,2), filled with IC
    head = 0
    ridx = np.arange(R)

    def partner_delayed():
        idx = (head - lag) % buflen
        Vd = Vbuf[idx, ridx, :]                                # (R,2) delayed self
        return Vd[:, ::-1]                                     # delayed partner

    def drift(V, W, Y, Vpd):
        coup = g * (Vpd - V)
        return (V - V**3 / 3.0 - W + Y + I0 + coup,
                EPS * (V + A - B * W), EPS * DELTA * (C - V))

    def rk4(V, W, Y):
        nonlocal head
        Vpd = partner_delayed()
        k1 = drift(V, W, Y, Vpd)
        k2 = drift(V + .5*dt*k1[0], W + .5*dt*k1[1], Y + .5*dt*k1[2], Vpd)
        k3 = drift(V + .5*dt*k2[0], W + .5*dt*k2[1], Y + .5*dt*k2[2], Vpd)
        k4 = drift(V + dt*k3[0], W + dt*k3[1], Y + dt*k3[2], Vpd)
        V = V + dt*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6
        W = W + dt*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6
        Y = Y + dt*(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6
        head = (head + 1) % buflen
        Vbuf[head] = V
        return V, W, Y

    for _ in range(int(warmup / dt)):
        V, W, Y = rk4(V, W, Y)
    n = int(T / dt); nrec = n // stride
    ts = np.empty(nrec); Vr = np.empty((nrec, R, 2)); ri = 0
    for k in range(n):
        V, W, Y = rk4(V, W, Y)
        if k % stride == 0 and ri < nrec:
            ts[ri] = k*dt; Vr[ri] = V; ri += 1
    return ts[:ri], Vr[:ri]


def classify(Vr, death_amp=0.8):
    R = Vr.shape[1]
    amp = Vr.max(0) - Vr.min(0)
    reg = np.zeros(R, int); rho = np.empty(R); chi = np.empty(R)
    for r in range(R):
        v = Vr[:, r, :]; v0 = v - v.mean(0)
        s1, s2 = v0[:, 0].std(), v0[:, 1].std()
        rho[r] = float(np.mean(v0[:, 0]*v0[:, 1])/(s1*s2)) if s1 > 0 and s2 > 0 else np.nan
        vbar = v.mean(1); chi[r] = np.sqrt(np.var(vbar)/np.mean(np.var(v, 0))) if np.mean(np.var(v,0))>0 else np.nan
        if amp[r].max() < death_amp:
            reg[r] = 0
        elif rho[r] > 0.3:
            reg[r] = 1
        elif rho[r] < -0.3:
            reg[r] = 2
        else:
            reg[r] = 3
    return reg, rho, chi


def main():
    t0 = time.time()
    ng, ntau = 16, 18
    g_vals = np.linspace(0.0, 0.30, ng)
    tau_vals = np.linspace(0.0, 50.0, ntau)
    GG, TT = np.meshgrid(g_vals, tau_vals, indexing="xy")
    g_flat = GG.ravel(); tau_flat = TT.ravel()

    ts, Vr = delay_grid(g_flat, tau_flat)
    reg, rho, chi = classify(Vr)
    REG = reg.reshape(ntau, ng); RHO = rho.reshape(ntau, ng)

    counts = {k: int((reg == k).sum()) for k in range(4)}
    names = {0: "death", 1: "in-phase", 2: "anti-phase", 3: "other"}
    print("=" * 70)
    print(f"T2.4 — delay atlas, two identical FHR (c={C}), electrical coupling")
    print("=" * 70)
    print("  regime cell counts:", {names[k]: counts[k] for k in range(4)})
    print(f"  ρ range [{np.nanmin(rho):+.2f}, {np.nanmax(rho):+.2f}]")
    flip = counts[2] > 5
    death = counts[0] > 2
    print(f"  delay-induced anti-phase present: {'YES' if flip else 'no'};  "
          f"delay-induced death present: {'YES' if death else 'no'}")

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    ext = [g_vals[0], g_vals[-1], tau_vals[0], tau_vals[-1]]
    cmap = ListedColormap(["#cfcfcf", "#2c7fb8", "#d95f0e", "#f7e463"])
    norm = BoundaryNorm([-.5, .5, 1.5, 2.5, 3.5], cmap.N)
    im0 = ax[0].imshow(REG, origin="lower", aspect="auto", extent=ext, cmap=cmap, norm=norm)
    cb = fig.colorbar(im0, ax=ax[0], ticks=[0, 1, 2, 3])
    cb.ax.set_yticklabels(["death", "in-phase", "anti-phase", "other"])
    ax[0].set_title("(1) regime in (g, τ)")

    im1 = ax[1].imshow(RHO, origin="lower", aspect="auto", extent=ext, cmap="RdBu_r",
                       vmin=-1, vmax=1)
    fig.colorbar(im1, ax=ax[1]); ax[1].set_title("(2) phase: ρ>0 in-phase, ρ<0 anti-phase")

    # time-series at fixed g showing the phase-flip
    gi = np.argmin(np.abs(g_vals - 0.15))
    col = REG[:, gi]
    ti_in = np.where(col == 1)[0]; ti_anti = np.where(col == 2)[0]
    win = ts < (ts[0] + 250)
    if ti_in.size and ti_anti.size:
        r_in = ti_in[0]*ng + gi; r_anti = ti_anti[len(ti_anti)//2]*ng + gi
        ax[2].plot(ts[win], Vr[win, r_in, 0], color="#2c7fb8", lw=0.9,
                   label=f"τ={tau_vals[ti_in[0]]:.0f}: in-phase")
        ax[2].plot(ts[win], Vr[win, r_in, 1], color="#2c7fb8", lw=0.9, ls=":")
        ax[2].plot(ts[win], Vr[win, r_anti, 0] + 4.5, color="#d95f0e", lw=0.9,
                   label=f"τ={tau_vals[ti_anti[len(ti_anti)//2]]:.0f}: anti-phase")
        ax[2].plot(ts[win], Vr[win, r_anti, 1] + 4.5, color="#d95f0e", lw=0.9, ls=":")
        ax[2].legend(fontsize=8.5, frameon=False)
    ax[2].set_title(f"(3) phase-flip at g={g_vals[gi]:.2f} (offset for clarity)")
    ax[2].set_xlabel("t"); ax[2].set_yticks([])

    for a_ in ax[:2]:
        a_.set_xlabel("coupling g"); a_.set_ylabel("delay τ")
    fig.suptitle(f"T2.4 — delay in electrical coupling (two identical FHR, c={C}, period≈54)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fp = os.path.join(FIG, "delay_atlas.png")
    fig.savefig(fp, dpi=130)
    print(f"\n  Figure: {fp}")
    print(f"  TOTAL {time.time()-t0:.1f}s — DONE")


if __name__ == "__main__":
    main()
