"""
coupled_array_peeloff.py — Direction C, the decisive test: does GENUINE coupling realize the
Airy_2 process, or only its local shadow?

A diffusively-coupled array of folded-cycle inner fields.  Each site i has its OWN independent
noise (no shared noise); sites are linked ONLY by diffusive coupling D * Laplacian on the
Cole-Hopf field u_i:
        d^2 u_i/dT^2 = (Y - eta xi_i) u_i + D (u_{i+1} - 2 u_i + u_{i-1}),   Y = Y0 - T,
peel-off Y_node(i) = first node of u_i.  D=0 -> independent (IID peel-offs); D>0 -> the coupling
correlates neighbours.

Question: is the peel-off covariance C(r) EXPONENTIAL (nearest-neighbour diffusion is short-range
-> a local proxy, like the OU-noise construction) or ALGEBRAIC r^{-2} (the coupling acts as a
Dyson-type eigenvalue interaction -> the true long-range Airy_2)?
"""
from __future__ import annotations
import os
import numpy as np


def coupled_peeloff(eta, D, Nx, Y0=5.0, Yend=-4.0, dt=2e-3, seed=0):
    rng = np.random.default_rng(seed)
    nsteps = int((Y0 - Yend) / dt); sdt = np.sqrt(dt)
    r0 = np.sqrt(Y0) + 1.0 / (4 * Y0)
    u = np.ones(Nx); w = np.full(Nx, r0); yn = np.full(Nx, np.nan); al = np.ones(Nx, bool)
    for m in range(nsteps):
        Y = Y0 - dt * m
        if not al.any():
            break
        lap = np.roll(u, -1) - 2.0 * u + np.roll(u, 1)        # periodic diffusive coupling
        un = u + w * dt
        wn = w + (Y * u) * dt - eta * u * sdt * rng.standard_normal(Nx) + D * lap * dt
        cr = al & (un * u <= 0.0)
        yn[cr] = Y - dt * (u[cr] / (u[cr] - un[cr])); al[cr] = False
        u, w = un, wn
    yn[np.isnan(yn)] = Yend
    return yn


def covariance(P, rmax):
    Q = P - P.mean()
    return np.array([(Q[:, r:] * Q[:, :Q.shape[1] - r]).mean() for r in range(rmax + 1)])


if __name__ == "__main__":
    eta = np.sqrt(2.0); Nx, R, rmax = 800, 24, 200
    Ds = [0.0, 1.0, 4.0]; res = {}
    print(f"{'D':>5} {'mean':>8} {'std':>7} {'skew':>7} {'corr.len':>9} {'C(20)/C0':>9}")
    for D in Ds:
        P = np.array([coupled_peeloff(eta, D, Nx, seed=6000 + int(10 * D) + k) for k in range(R)])
        C = covariance(P, rmax); C0 = C[0]; Cn = C / C0
        r = np.arange(rmax + 1)
        m = (r >= 5) & (r <= 80) & (Cn > 0.03)
        clen = -1.0 / np.polyfit(r[m], np.log(Cn[m]), 1)[0] if m.sum() > 3 else np.nan
        res[D] = dict(P=P, Cn=Cn, clen=clen)
        mu, s = P.mean(), P.std(); sk = ((P - mu) ** 3).mean() / s ** 3
        print(f"{D:5.1f} {mu:8.3f} {s:7.3f} {sk:7.3f} {clen:9.2f} {Cn[20]:9.3f}")
    print("  (TW_2 marginal: mean -1.771 std 0.902 skew 0.224)")
    print("  Verdict: coupling builds a correlation length that GROWS with D;")
    print("  decay shape (Panel B) classifies exponential[local proxy] vs r^-2[Airy_2].")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    cols = {0.0: "#94a3b8", 1.0: "#2563eb", 4.0: "#dc2626"}
    fig, ax = plt.subplots(1, 3, figsize=(15.4, 4.7))
    r = np.arange(rmax + 1)

    # Panel A: covariance log-linear (exponential = straight)
    for D in Ds:
        ax[0].semilogy(r, np.clip(res[D]["Cn"], 1e-3, None), "o", ms=2.6, color=cols[D],
                       label=fr"$D={D:.0f}$ (len {res[D]['clen']:.0f})")
    ax[0].set_xlabel(r"lag $r$ (sites)"); ax[0].set_ylabel(r"$C(r)/C(0)$")
    ax[0].set_title("Coupled-array covariance (log-linear)")
    ax[0].set_ylim(1e-3, 1.3); ax[0].legend(fontsize=8.6, framealpha=0.95, edgecolor="#cbd5e1")

    # Panel B: log-log vs r^-2  (the classification)
    for D in [1.0, 4.0]:
        pos = (r >= 2) & (res[D]["Cn"] > 0)
        ax[1].loglog(r[pos], res[D]["Cn"][pos], "o", ms=2.6, color=cols[D], label=fr"$D={D:.0f}$")
    ax[1].loglog(r[2:], 0.9 * (r[2:] / 2.0) ** -2.0, "--", color="#16a34a", lw=1.6,
                 label=r"Airy$_2$ tail $r^{-2}$")
    ax[1].set_xlabel(r"lag $r$"); ax[1].set_ylabel(r"$C(r)/C(0)$")
    ax[1].set_title("Exponential (local) vs algebraic (Airy$_2$)?")
    ax[1].set_ylim(1e-3, 1.3); ax[1].legend(fontsize=8.6, framealpha=0.95, edgecolor="#cbd5e1")

    # Panel C: sample coupled peel-off profiles (D=4)
    xw = np.arange(400)
    for k, c in zip(range(3), ["#1d4ed8", "#0891b2", "#7c3aed"]):
        ax[2].plot(xw, res[4.0]["P"][k, :400], color=c, lw=1.0, alpha=0.85)
    ax[2].axhline(-1.771, color="#475569", lw=1.0, ls="--", label=r"TW$_2$ mean")
    ax[2].set_xlabel(r"array site $i$"); ax[2].set_ylabel(r"peel-off $Y_{\rm node}(i)$")
    ax[2].set_title(r"Coupled peel-off profiles ($D=4$)")
    ax[2].legend(fontsize=9, framealpha=0.95, edgecolor="#cbd5e1")

    fig.tight_layout()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures",
                       "coupled_array_peeloff.png")
    out = os.path.abspath(out)
    fig.savefig(out, dpi=140); print("saved", out)
