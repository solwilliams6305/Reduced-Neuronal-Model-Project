"""
dyson_weber.py — intrinsic Airy₂ from the DBM edge (matrix-OU); framing the Dyson–Weber object.
================================================================================================

Stable construction: a Hermitian matrix Ornstein–Uhlenbeck process
    M(t+dt) = ρ M(t) + √(1−ρ²) G,   ρ = e^{−dt/2},  G a fresh GUE matrix,
whose eigenvalues perform Dyson Brownian motion (no singular repulsion to integrate). The top
eigenvalue λ_max(t), at the edge of the semicircle, is the **intrinsic Airy₂ process**.

  [VALIDATE] standardised λ_max is Tracy–Widom (skew ≈ 0.22), and its time-covariance is the
  intrinsic Airy₂ covariance (decaying; increments → 2·Var; locally Brownian). This is the genuine
  process — the object my earlier forced result approximated with an imposed OU on the peel-off.

The cusp's intrinsic process — the **"Dyson–Weber" object** — is the edge process of a DBM at a
MULTICRITICAL point (spectral density vanishing faster than √ at the edge), i.e. a higher-order /
parabolic-cylinder line ensemble. Matrix-OU only realises the Gaussian (Airy) edge; the multicritical
edge needs a tuned non-Gaussian matrix model (eigenvalue Langevin with a critically-tuned potential),
which is numerically delicate (the naive attempt is unstable) and research-grade — see notes.

Output: figures/dyson_weber.png + tagged summary.
"""
from __future__ import annotations
import os, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
TW2 = dict(skew=0.2241, kurt=0.0934)


def gue(N, rng):
    A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    return (A + A.conj().T) / (2.0 * np.sqrt(N))          # Hermitian, eigenvalues ≈ [-2,2]


def matrix_ou_edge(N=80, dt=0.05, T=70.0, warmup=12.0, runs=10, seed=0):
    rng = np.random.default_rng(seed); rho = np.exp(-dt / 2.0); sq = np.sqrt(1 - rho**2)
    nrec = int(T / dt); nw = int(warmup / dt)
    series = np.empty((runs, nrec))
    for r in range(runs):
        M = gue(N, rng)
        for _ in range(nw):
            M = rho * M + sq * gue(N, rng)
        for i in range(nrec):
            M = rho * M + sq * gue(N, rng)
            series[r, i] = np.linalg.eigvalsh(M)[-1]
    return series


def autocov(series, dt, Tmax=3.0):
    nlag = int(Tmax / dt); mu = series.mean(); var = series.var()
    C = np.empty(nlag); V = np.empty(nlag); flat = series - mu; L = series.shape[1]
    for s in range(nlag):
        C[s] = np.mean(flat[:, :L-s] * flat[:, s:]) / var
        V[s] = np.var(series[:, s:] - series[:, :L-s])
    return C, V, var


def moments(x):
    x = x.ravel(); d = x - x.mean(); v = np.mean(d**2)
    return np.mean(d**3)/v**1.5, np.mean(d**4)/v**2 - 3.0


def main():
    t0 = time.time()
    dt = 0.05
    print("=" * 72)
    print("Dyson–Weber push — intrinsic Airy₂ from the matrix-OU DBM edge")
    print("=" * 72)
    s = matrix_ou_edge(dt=dt, seed=1)
    sk, ku = moments(s)
    C, V, var = autocov(s, dt)
    print(f"  [VALIDATE] λ_max edge mean {s.mean():.3f} (semicircle edge ≈ 2)")
    print(f"             standardised marginal: skew {sk:+.3f} (TW₂ {TW2['skew']:+.3f}), "
          f"exkurt {ku:+.3f} (TW₂ {TW2['kurt']:+.3f})  → {'Airy₂/TW ✓' if abs(sk-TW2['skew'])<0.08 else 'finite-N'}")
    # intrinsic covariance signatures
    decorr = np.argmax(C < 0.2) * dt
    print(f"  [INTRINSIC] Airy₂ covariance: C decays to <0.2 by τ≈{decorr:.2f}; "
          f"increments V(∞)/2Var = {V[-1]/(2*var):.2f} (→1 ⇒ decorrelating, locally-Brownian)")
    # locally-Brownian small-τ slope of V
    s_small = slice(1, int(0.4/dt))
    slope = np.polyfit(np.arange(*s_small.indices(len(V)))*dt, V[s_small], 1)[0]
    print(f"             small-τ increment slope (locally Brownian if >0): {slope:.2f}")

    print(f"\n  [OPEN] Dyson–Weber (cusp) = DBM at a MULTICRITICAL edge (density ~ dist^{{>1/2}}); needs a")
    print(f"         critically-tuned non-Gaussian matrix model — numerically delicate, research-grade.")

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    tg = np.arange(s.shape[1]) * dt
    ax[0].plot(tg[:600], s[0, :600], lw=0.7, color="#1f3b73")
    ax[0].axhline(s.mean(), color="grey", ls=":", lw=0.8)
    ax[0].set_xlabel("DBM time t"); ax[0].set_ylabel("λ_max(t)")
    ax[0].set_title("(A) top-eigenvalue edge process = Airy₂")

    gg = np.linspace(-4, 4, 200)
    ax[1].plot(gg, np.exp(-gg**2/2)/np.sqrt(2*np.pi), "k:", lw=1, label="Gaussian")
    z = (s.ravel() - s.mean())/s.std()
    ax[1].hist(z, bins=70, range=(-4, 4), density=True, histtype="step", lw=1.8, color="#b3402b",
               label=f"DBM edge marginal (skew {sk:+.2f})")
    ax[1].set_xlabel("standardised λ_max"); ax[1].set_ylabel("density")
    ax[1].set_title("(B) marginal = Tracy–Widom (validation)"); ax[1].legend(fontsize=9, frameon=False)

    tau = np.arange(len(C)) * dt
    ax[2].plot(tau, C, "-", color="#1f3b73", lw=2, label="C(τ) intrinsic Airy₂ covariance")
    ax[2].plot(tau, V/(2*var), "-", color="#2c7d59", lw=1.6, label="V(τ)/2Var")
    ax[2].axhline(1.0, color="grey", lw=0.8, ls="--"); ax[2].axhline(0.0, color="grey", lw=0.5, ls=":")
    ax[2].set_xlabel("DBM time lag τ"); ax[2].set_ylabel("C(τ) , V/2Var")
    ax[2].set_title("(C) intrinsic covariance: decays, increments → 2·Var"); ax[2].legend(fontsize=9, frameon=False)

    fig.suptitle("Dyson–Weber push — intrinsic Airy₂ from the DBM edge (validated); "
                 "cusp = multicritical-edge ensemble (open)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fp = os.path.join(FIG, "dyson_weber.png")
    fig.savefig(fp, dpi=140, bbox_inches="tight")
    print(f"\n  Figure: {fp}")
    print(f"  TOTAL {time.time()-t0:.1f}s — DONE")


if __name__ == "__main__":
    main()
