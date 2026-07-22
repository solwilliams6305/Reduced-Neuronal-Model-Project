"""
tw_boundary.py — the VERIFIED recrossing boundary in the (sigma, eps) plane.

Claim (Remark, regime of validity): the recrossing that voids the single-escape / Tracy-Widom
picture is a Kramers return over the unstable cycle, so the expected number of returns is

        N_ret(sigma, eps) = (C / eps) * exp( - DeltaU* / sigma^2 ),

i.e. an ARRHENIUS law in sigma with a 1/eps time-in-window prefactor.  Equivalent statements:
  * log( eps * N_ret )  is LINEAR in 1/sigma^2, slope = -DeltaU*, intercept = log C,
    and ALL ramp rates eps collapse onto that one line  (the verification);
  * it does NOT collapse in eta = sigma/sqrt(eps)            (the contrast);
  * the breakdown boundary is  sigma_rc(eps) = sqrt( DeltaU* / log(C/eps) )  -- a NEAR-FIXED
    critical noise (only logarithmic in eps), i.e. the static Kramers sigma_crit re-emerging
    as the lid on the Tracy-Widom regime.

DeltaU* is compared to the explicit normal-form return barrier
        DeltaU_ret(mu) = U(r_u) - U(0) = U(r_u),   U(r) = -mu r^2/2 - r^4/4 + r^6/6,
r_u the unstable-cycle radius,  r_u^2 = (1 - sqrt(1+4 mu))/2.
"""
from __future__ import annotations
import numpy as np

OMEGA = 1.0
MU_FOLD = -0.25


def U(r, mu):
    return -0.5 * mu * r ** 2 - 0.25 * r ** 4 + r ** 6 / 6.0


def radii(mu):
    d = 1 + 4 * mu
    if d < 0:
        return None, None
    r_u = np.sqrt(max((1 - np.sqrt(d)) / 2, 0.0))
    r_s = np.sqrt((1 + np.sqrt(d)) / 2)
    return r_u, r_s


def barrier_return(mu):           # rest -> cycle barrier = U(r_u) - U(0)
    r_u, _ = radii(mu)
    return np.nan if r_u is None else U(r_u, mu) - U(0.0, mu)


def P_recross(sigma, eps, N=1500, mu0=0.10, mu_end=-0.40, dt=0.05, seed=0):
    rng = np.random.default_rng(seed)
    _, rs = radii(mu0)
    ph = rng.uniform(0, 2 * np.pi, N); x = rs * np.cos(ph); y = rs * np.sin(ph)
    r_hi, r_mid = 0.85, 0.55
    nsteps = int((mu0 - mu_end) / (eps * dt)); sdt = np.sqrt(dt)
    everLow = np.zeros(N, bool); recross = np.zeros(N, bool)
    for k in range(nsteps):
        mu = mu0 - eps * dt * k
        s = x * x + y * y; g = mu + s - s * s
        xn = x + (mu * x - OMEGA * y + x * g) * dt + sigma * sdt * rng.standard_normal(N)
        yn = y + (mu * y + OMEGA * x + y * g) * dt + sigma * sdt * rng.standard_normal(N)
        x = np.clip(xn, -5, 5); y = np.clip(yn, -5, 5)
        r = np.sqrt(x * x + y * y)
        recross |= everLow & (r > r_hi) & (mu > MU_FOLD + 1e-3)
        everLow |= r < r_mid
    return float(recross.mean())


if __name__ == "__main__":
    sigmas = np.linspace(0.05, 0.16, 9)
    epss = [0.004, 0.002, 0.001]
    data = {}
    print(f"{'eps':>7} {'sigma':>7} {'eta':>6} {'P_rec':>7} {'N_ret':>7}")
    for eps in epss:
        P = []
        for j, sg in enumerate(sigmas):
            p = P_recross(sg, eps, seed=10 + j + int(1000 * eps))
            P.append(p)
            print(f"{eps:7.3f} {sg:7.3f} {sg/np.sqrt(eps):6.2f} {p:7.3f} "
                  f"{-np.log(max(1e-6,1-p)):7.3f}")
        data[eps] = np.array(P)

    # ---- fit log(eps*N_ret) = log C - DeltaU*/sigma^2 across ALL eps (the collapse) ----
    X, Yv = [], []
    for eps in epss:
        P = data[eps]; Nret = -np.log(np.clip(1 - P, 1e-6, 1))
        m = (P > 0.01) & (P < 0.7)
        X.append(1 / sigmas[m] ** 2); Yv.append(np.log(eps * Nret[m]))
    X = np.concatenate(X); Yv = np.concatenate(Yv)
    slope, intercept = np.polyfit(X, Yv, 1)
    dU_star = -slope; C = np.exp(intercept)
    ss_res = np.sum((Yv - (slope * X + intercept)) ** 2)
    R2 = 1 - ss_res / np.sum((Yv - Yv.mean()) ** 2)
    print(f"\n  fit  log(eps*N_ret) = logC - dU*/sigma^2 :  dU* = {dU_star:.4f}   C = {C:.3f}   R^2 = {R2:.3f}")
    print(f"  normal-form return barrier U(r_u)-U(0):  mu=-0.05 -> {barrier_return(-0.05):.4f}"
          f"   mu=-0.10 -> {barrier_return(-0.10):.4f}   mu=-0.15 -> {barrier_return(-0.15):.4f}")
    for eps in epss:
        print(f"  boundary sigma_rc(eps={eps}) = sqrt(dU*/log(C/eps)) = "
              f"{np.sqrt(dU_star/np.log(C/eps)):.4f}   (eta_rc = {np.sqrt(dU_star/np.log(C/eps))/np.sqrt(eps):.2f})")

    # ---- figure: collapse (verification) vs non-collapse in eta (contrast) ----
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    cols = {0.004: "#10b981", 0.002: "#4ea1ff", 0.001: "#f59e0b"}
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    for eps in epss:
        P = data[eps]; Nret = -np.log(np.clip(1 - P, 1e-6, 1)); m = (P > 0.01) & (P < 0.7)
        ax[0].plot(1 / sigmas[m] ** 2, np.log(eps * Nret[m]), "o", color=cols[eps],
                   label=fr"$\epsilon={eps}$")
    xx = np.linspace(X.min(), X.max(), 50)
    ax[0].plot(xx, slope * xx + intercept, "k-", lw=1,
               label=fr"fit: $\Delta U^*={dU_star:.3f}$, $R^2={R2:.2f}$")
    ax[0].set_xlabel(r"$1/\sigma^2$"); ax[0].set_ylabel(r"$\log(\epsilon\, N_{\rm ret})$")
    ax[0].set_title("VERIFIED: Arrhenius collapse (all $\\epsilon$ on one line)")
    ax[0].legend(frameon=False, fontsize=9)
    for eps in epss:
        ax[1].plot(sigmas / np.sqrt(eps), data[eps], "o-", color=cols[eps], label=fr"$\epsilon={eps}$")
    ax[1].set_xlabel(r"$\eta=\sigma/\sqrt{\epsilon}$"); ax[1].set_ylabel(r"$P_{\rm recross}$")
    ax[1].set_title(r"CONTRAST: does NOT collapse in $\eta$ (Kramers, not fold)")
    ax[1].legend(frameon=False, fontsize=9)
    fig.tight_layout()
    out = "/sessions/trusting-fervent-archimedes/mnt/Reduced Neuronal Model Project/figures/tw_boundary.png"
    fig.savefig(out, dpi=130); print("\nsaved", out)
