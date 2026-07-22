"""
airy_covariance.py — Direction C, next step: is the eigenvalue field the FULL Airy_2, or only a
local proxy?

The single-point marginal is TW_2 and the local roughness -> Hölder-1/2 (Airy_2) as eta->0.
But the LONG-RANGE structure is the real test.  Linear response: the SAO ground state
Lambda0[xi(.,x)] is a local functional of the noise column, so if the noise is OU (exponential
correlation in x) the eigenvalue covariance is EXPONENTIAL, C(r) ~ e^{-r/ell} -- whereas the
true Airy_2 covariance has an ALGEBRAIC tail C(r) ~ r^{-2}.  We test this directly:

  (A) SAO-OU eigenvalue Lambda0(x) (Sturm) and node Y_node(x): covariance C(r) -- exponential?
  (B) Dyson-OU GUE edge lambda_max(t): the CANONICAL Airy_2 -- algebraic tail, for contrast.

Verdict: correlated-noise coupling realizes Airy_2's marginal + local roughness, but its
long-range correlations are OU (exponential), not Airy_2 (algebraic); the true long-range Airy_2
needs Dyson-type eigenvalue dynamics (genuine spectral interaction), not just correlated noise.
"""
from __future__ import annotations
import os
import numpy as np
from eigenvalue_vs_node import ou_field
from airy_process_peeloffs import peeloff_profile


def lam0_sturm(eta, Nx, ell, M=600, L=7.0, seed=0):
    rng = np.random.default_rng(seed); h = L / M; Yj = np.arange(1, M + 1) * h
    Z = ou_field(M, Nx, ell, rng); diag = 2 / h ** 2 + Yj[:, None] + eta * Z / np.sqrt(h)
    e2 = (1 / h ** 2) ** 2

    def sturm(lam):
        q = diag[0] - lam; cnt = (q < 0).astype(np.int32)
        for i in range(1, M):
            q = (diag[i] - lam) - e2 / np.where(np.abs(q) < 1e-300, -1e-300, q); cnt += (q < 0)
        return cnt
    lo = np.full(Nx, -5.0); hi = np.full(Nx, 12.0)
    for _ in range(44):
        mid = 0.5 * (lo + hi); has = sturm(mid) >= 1; hi = np.where(has, mid, hi); lo = np.where(has, lo, mid)
    return 0.5 * (lo + hi)


def covariance(profiles, rmax):
    P = profiles - profiles.mean()
    C = np.array([(P[:, r:] * P[:, :P.shape[1] - r]).mean() for r in range(rmax + 1)])
    return C


def gue(rng, N):
    X = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    return (X + X.conj().T) / 2.0 / np.sqrt(N)              # Hermitian, semicircle edge ~ 2


def dyson_lmax(N, T, dt=0.04, tau=1.0, seed=0):
    rng = np.random.default_rng(seed); a = np.exp(-dt / tau); b = np.sqrt(1 - a * a)
    H = gue(rng, N); lm = np.empty(T)
    for t in range(T):
        lm[t] = np.linalg.eigvalsh(H)[-1]
        H = a * H + b * gue(rng, N)
    return lm


if __name__ == "__main__":
    eta = np.sqrt(2.0); Nx, ell, R = 700, 80.0, 24
    L0 = np.array([lam0_sturm(eta, Nx, ell, seed=3000 + k) for k in range(R)])
    YN = np.array([peeloff_profile(eta, Nx, ell, seed=4000 + k) for k in range(R)])
    rmax = 360
    Ce = covariance(L0, rmax); Ce /= Ce[0]
    Cn = covariance(YN, rmax); Cn /= Cn[0]
    r = np.arange(rmax + 1)
    # exponential-fit decay length of the eigenvalue covariance
    m = (r >= 10) & (r <= 200) & (Ce > 0.02)
    decay = -1.0 / np.polyfit(r[m], np.log(Ce[m]), 1)[0]
    print(f"[SAO-OU] eigenvalue C(r): exp-fit decay length = {decay:.1f}  (cf OU ell={ell:.0f})")
    print(f"         C(50)={Ce[50]:.3f}, C(150)={Ce[150]:.3f}, C(300)={Ce[300]:.4f}  "
          f"(exp e^-r/ell: {np.exp(-50/ell):.3f}, {np.exp(-150/ell):.3f}, {np.exp(-300/ell):.4f})")

    # Dyson-OU edge: canonical Airy_2 reference
    N = 60; lm = np.concatenate([dyson_lmax(N, 9000, seed=10 + k) for k in range(3)])
    lm0 = lm - lm.mean(); Tcov = 400
    Cd = np.array([(lm0[:len(lm0) - t] * lm0[t:]).mean() for t in range(Tcov + 1)])
    Cd /= Cd[0]
    print(f"[Dyson]  lambda_max C(t): C(20)={Cd[20]:.3f}, C(80)={Cd[80]:.3f}, C(200)={Cd[200]:.4f}"
          f"  (slower-than-exponential => algebraic Airy_2 tail)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    fig, ax = plt.subplots(1, 3, figsize=(15.4, 4.7))

    # Panel A: SAO-OU covariances on log-LINEAR -> exponential = straight line
    ax[0].semilogy(r, Ce, "o", color="#2563eb", ms=3, label=r"eigenvalue $\Lambda_0$")
    ax[0].semilogy(r, np.clip(Cn, 1e-3, None), "s", color="#dc2626", ms=3, label=r"node $Y_{\rm node}$")
    ax[0].semilogy(r, np.exp(-r / ell), "--", color="#0f172a", lw=1.6, label=fr"$e^{{-r/\ell}}$, $\ell={ell:.0f}$")
    ax[0].set_xlabel(r"lag $r$"); ax[0].set_ylabel(r"$C(r)/C(0)$")
    ax[0].set_title("SAO-OU covariance is EXPONENTIAL (log-linear)")
    ax[0].set_ylim(1e-3, 1.3); ax[0].legend(fontsize=8.6, framealpha=0.95, edgecolor="#cbd5e1")

    # Panel B: same on log-LOG with r^-2 reference -> not algebraic
    pos = (r >= 2) & (Ce > 0)
    ax[1].loglog(r[pos], Ce[pos], "o", color="#2563eb", ms=3, label=r"SAO-OU $\Lambda_0$ (exponential)")
    ax[1].loglog(r[pos], 0.9 * (r[pos] / 2.0) ** -2.0, "--", color="#16a34a", lw=1.6,
                 label=r"Airy$_2$ tail $r^{-2}$")
    ax[1].set_xlabel(r"lag $r$"); ax[1].set_ylabel(r"$C(r)/C(0)$")
    ax[1].set_title(r"SAO-OU is NOT the algebraic Airy$_2$ tail")
    ax[1].set_ylim(1e-3, 1.3); ax[1].legend(fontsize=8.6, framealpha=0.95, edgecolor="#cbd5e1")

    # Panel C: Dyson-OU edge (canonical Airy_2) -> slower/algebraic tail
    td = np.arange(Tcov + 1); pd = (td >= 2) & (Cd > 0)
    ax[2].loglog(td[pd], Cd[pd], "o", color="#7c3aed", ms=3, label=r"Dyson edge $\lambda_{\max}$")
    ax[2].loglog(td[pd], 0.8 * (td[pd] / 2.0) ** -2.0, "--", color="#16a34a", lw=1.6, label=r"$t^{-2}$")
    ax[2].set_xlabel(r"lag $t$"); ax[2].set_ylabel(r"$C(t)/C(0)$")
    ax[2].set_title(r"Dyson (true Airy$_2$): heavier, algebraic tail")
    ax[2].set_ylim(1e-3, 1.3); ax[2].legend(fontsize=8.6, framealpha=0.95, edgecolor="#cbd5e1")

    fig.tight_layout()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures",
                       "airy_covariance.png")
    out = os.path.abspath(out)
    fig.savefig(out, dpi=140); print("saved", out)
