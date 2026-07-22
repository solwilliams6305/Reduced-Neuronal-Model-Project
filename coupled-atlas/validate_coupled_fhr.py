"""
validate_coupled_fhr.py — validation gate for the coupled FHR engine.
=====================================================================

All conditions are stepped TOGETHER as replicas in one vectorized RK4 loop
(state shape (R,2)), so the whole gate runs in a single pass.

Checks:
  [1] MMO staircase — a single unit shows more small oscillations per large
      spike as c -> -1 (folded-node staircase; cf. mmo_fhr_alpha_fine.py).
  [2] g=0 decoupling — unit 1 of a pair at g=0 equals the isolated unit.
  [3] Known coupled result — two detuned units phase-lock as electrical
      coupling g grows: Golomb–Hansel χ rises from ~0 toward ~1.

Outputs: figures/validation_coupled_fhr.png + printed summary.
"""
from __future__ import annotations
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from coupled_fhr import mmo_counts, sync_chi, large_spike_times, phase_difference

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

A, B, EPS, DELTA, I0 = 0.7, 0.8, 0.08, 0.2, 0.30   # repo FHR tuning


def run_all(replicas, T=2200.0, dt=0.04, warmup=900.0, stride=4):
    """Step all (c1,c2,g) replicas together. Returns ts, Vr with Vr (nrec,R,2)."""
    R = len(replicas)
    C = np.array([[r[0], r[1]] for r in replicas], float)   # (R,2)
    g = np.array([r[2] for r in replicas], float)[:, None]  # (R,1)
    V = np.tile([0.6, -0.6], (R, 1)).astype(float)
    W = np.zeros((R, 2)); Y = np.zeros((R, 2))

    def drift(V, W, Y):
        coup = g * (V[:, ::-1] - V)                         # linear gap junction
        return (V - V**3 / 3.0 - W + Y + I0 + coup,
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


def period_of(t, v):
    ts = large_spike_times(t, v)
    return float(np.median(np.diff(ts))) if ts.size >= 3 else np.nan


def main():
    import time
    t0 = time.time()
    print("=" * 70, flush=True)
    print("VALIDATION — coupled FitzHugh–Rinzel (single vectorized pass)", flush=True)
    print("=" * 70, flush=True)

    mmo_cs = [-0.78, -0.86, -0.92]
    g_arr = np.linspace(0.0, 0.30, 13)
    c1, c2 = -0.82, -0.72
    replicas = ([(c, c, 0.0) for c in mmo_cs] +
                [(c1, c2, g) for g in g_arr] +
                [(c1, c1, 0.0)])                 # last = isolated reference
    iM = {c: i for i, c in enumerate(mmo_cs)}
    iS0 = len(mmo_cs)                            # first sync replica (g=0)
    iSync = list(range(iS0, iS0 + len(g_arr)))
    iRef = len(replicas) - 1

    ts, Vr = run_all(replicas)
    print(f"(integration {time.time()-t0:.1f}s, {Vr.shape[1]} replicas)\n", flush=True)

    # [1] MMO staircase
    print("[1] Single-unit MMO staircase (more S per L as c -> -1):", flush=True)
    print(f"    {'c':>7}{'nL':>5}{'nS':>5}{'nS/nL':>8}{'period':>9}", flush=True)
    prev_ratio = -1; mono = True
    for c in mmo_cs:
        v = Vr[:, iM[c], 0]
        mc = mmo_counts(v)
        r = mc["ratio_S_per_L"]
        if not np.isnan(r):
            mono &= (r >= prev_ratio - 0.2); prev_ratio = r
        print(f"    {c:7.2f}{mc['nL']:5d}{mc['nS']:5d}{r:8.2f}{period_of(ts, v):9.1f}",
              flush=True)
    print(f"    staircase monotone (S/L grows toward c=-1): {'PASS' if mono else 'CHECK'}\n",
          flush=True)

    # [2] g=0 decoupling
    dev = np.max(np.abs(Vr[:, iS0, 0] - Vr[:, iRef, 0]))   # both c1, g=0
    print("[2] g=0 decoupling correctness:", flush=True)
    print(f"    max|v1_pair(g=0) - v1_isolated| = {dev:.2e}  "
          f"({'PASS' if dev < 1e-9 else 'CHECK'})\n", flush=True)

    # [3] synchronization vs coupling
    p1 = period_of(ts, Vr[:, iS0, 0]); p2 = period_of(ts, Vr[:, iS0, 1])
    chi = np.array([sync_chi(Vr[:, i, :]) for i in iSync])
    phi = np.array([phase_difference(ts, Vr[:, i, :]) for i in iSync])
    print("[3] Synchronization of two detuned units (χ should rise with g):", flush=True)
    print(f"    detuned via c: c1={c1} (solo period {p1:.1f}), "
          f"c2={c2} (solo period {p2:.1f})", flush=True)
    for gg, ch, ph in zip(g_arr, chi, phi):
        sph = "nan" if np.isnan(ph) else f"{ph:.3f}"
        print(f"      g={gg:5.3f}   χ={ch:5.3f}   phasediff={sph}", flush=True)
    print(f"    χ(g=0)={chi[0]:.3f} -> χ(g={g_arr[-1]:.2f})={chi[-1]:.3f}  "
          f"({'PASS' if chi[-1] > chi[0] + 0.2 else 'CHECK'})\n", flush=True)

    # ---- figure ----
    fig, ax = plt.subplots(2, 2, figsize=(13, 8))
    vA = Vr[:, iM[-0.86], 0]; wA = ts < (ts[0] + 1200)
    ax[0, 0].plot(ts[wA], vA[wA], lw=0.8, color="#1f3b73")
    ax[0, 0].axhline(0, color="grey", lw=0.5, ls=":")
    ax[0, 0].set_title("(A) single FHR unit, c=-0.86 — MMO staircase")
    ax[0, 0].set_xlabel("t"); ax[0, 0].set_ylabel("v")

    win = ts < (ts[0] + 1400)
    ax[0, 1].plot(ts[win], Vr[win, iSync[0], 0], lw=0.8, label="unit 1")
    ax[0, 1].plot(ts[win], Vr[win, iSync[0], 1], lw=0.8, label="unit 2", alpha=0.8)
    ax[0, 1].set_title(f"(B) detuned pair, g=0 — independent (χ={chi[0]:.2f})")
    ax[0, 1].set_xlabel("t"); ax[0, 1].set_ylabel("v"); ax[0, 1].legend(fontsize=8)

    ax[1, 0].plot(ts[win], Vr[win, iSync[-1], 0], lw=0.8, label="unit 1")
    ax[1, 0].plot(ts[win], Vr[win, iSync[-1], 1], lw=0.8, label="unit 2", alpha=0.8)
    ax[1, 0].set_title(f"(C) same pair, g={g_arr[-1]:.2f} — locked (χ={chi[-1]:.2f})")
    ax[1, 0].set_xlabel("t"); ax[1, 0].set_ylabel("v"); ax[1, 0].legend(fontsize=8)

    ax[1, 1].plot(g_arr, chi, "o-", color="#b3402b")
    ax[1, 1].set_ylim(0, 1.02)
    ax[1, 1].set_title("(D) synchrony χ vs coupling g (detuned units)")
    ax[1, 1].set_xlabel("coupling strength g"); ax[1, 1].set_ylabel("Golomb–Hansel χ")
    ax[1, 1].axhline(1.0, color="grey", lw=0.5, ls=":")

    fig.suptitle("Coupled FitzHugh–Rinzel — validation gate", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fp = os.path.join(FIG, "validation_coupled_fhr.png")
    fig.savefig(fp, dpi=130)
    print(f"Figure: {fp}", flush=True)
    print(f"TOTAL {time.time()-t0:.1f}s — DONE", flush=True)


if __name__ == "__main__":
    main()
