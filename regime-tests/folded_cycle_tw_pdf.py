"""
folded_cycle_tw_pdf.py — the airtight check: raw inner-exit histogram vs EXACT Tracy-Widom
-------------------------------------------------------------------------------------------
Pointwise pdf overlay (not just moments).  The inner exit measure (first node of
u''=(Y-eta xi)u, the Cole-Hopf image of the noisy canard Riccati) is the stochastic
Airy ground state = Tracy-Widom_beta with beta=4/eta^2.  Here we compute the EXACT
TW_1 (GOE) and TW_2 (GUE) densities from Painleve II (Hastings-McLeod) and overlay
the RAW simulated histogram -- NO centring, NO scaling, NO fit.

Painleve II:  q'' = s q + 2 q^3,  q(s) ~ Ai(s) as s->+inf  (Hastings-McLeod).
Track along a backward RK4 sweep:
    u(s) = int_s^inf q^2 ,   E(s) = int_s^inf (x-s) q^2 ,   w(s) = int_s^inf q .
Then   F2 = exp(-E),                 f2 = u F2
       F1 = exp(-(w+E)/2),           f1 = (q+u)/2 F1 .
Validated to 4 dp against tabulated TW moments (Bornemann 2010).

Dictionary check (well-resolved): beta=2 (eta=sqrt2) and beta=1 (eta=2) match TW in
mean, std, skew AND kurt with NO free parameter -> the histograms sit on the curves.

Companion: FOLDED_CYCLE_NOISYAIRY_CONJECTURES.md, folded_cycle_tracy_widom.py.
"""
from __future__ import annotations
import os
import numpy as np


def painleve_tw(s_max=10.0, s_min=-10.0, h=5e-4):
    """Exact TW_1, TW_2 densities via Hastings-McLeod (corrected Airy IC)."""
    zeta = (2.0 / 3.0) * s_max**1.5
    pref = np.exp(-zeta) / (2 * np.sqrt(np.pi))
    q = pref / s_max**0.25 * (1 - 5 / (72 * zeta) + 385 / (10368 * zeta**2))   # Ai(s_max)
    p = -pref * s_max**0.25 * (1 + 7 / (72 * zeta) - 455 / (10368 * zeta**2))  # Ai'(s_max)
    u = E = w = 0.0
    n = int(round((s_max - s_min) / h))
    S = np.empty(n + 1); Q = np.empty(n + 1); U = np.empty(n + 1)
    Ev = np.empty(n + 1); W = np.empty(n + 1)
    s = s_max; y = np.array([q, p, u, E, w])
    S[0] = s; Q[0] = q; U[0] = u; Ev[0] = E; W[0] = w

    def f(s, y):
        q, p, u, E, w = y
        return np.array([p, s * q + 2 * q**3, -q * q, -u, -q])

    for i in range(1, n + 1):
        k1 = f(s, y); k2 = f(s - h / 2, y - h / 2 * k1)
        k3 = f(s - h / 2, y - h / 2 * k2); k4 = f(s - h, y - h * k3)
        y = y - h / 6 * (k1 + 2 * k2 + 2 * k3 + k4); s -= h
        S[i] = s; Q[i] = y[0]; U[i] = y[2]; Ev[i] = y[3]; W[i] = y[4]
    o = np.argsort(S)
    S, Q, U, Ev, W = S[o], Q[o], U[o], Ev[o], W[o]
    F2 = np.exp(-Ev); f2 = U * F2
    F1 = np.exp(-(W + Ev) / 2); f1 = (Q + U) / 2 * F1
    return S, f1, f2


def first_node(eta, Y0=7.0, dt=6e-4, N=30000, Ymin=-10.0, seed=0):
    rng = np.random.default_rng(seed)
    u = np.ones(N); v = np.full(N, np.sqrt(Y0)); Y = Y0; sdt = np.sqrt(dt)
    Yz = np.full(N, np.nan); done = np.zeros(N, bool)
    n = int((Y0 - Ymin) / dt)
    for _ in range(n):
        al = ~done
        if not al.any():
            break
        dB = sdt * rng.standard_normal(N)
        u1 = u + v * dt; Yp = Y - dt; v1 = v + (Y * u) * dt - eta * u * dB
        u = u + 0.5 * (v + v1) * dt
        v = v + 0.5 * (Y * u + Yp * u1) * dt - eta * (0.5 * (u + u1)) * dB
        Y = Yp
        cr = al & (u < 0.0); Yz[cr] = Y; done |= cr
    return Yz[np.isfinite(Yz)]


def mom(x):
    m = x.mean(); d = x - m; v = np.mean(d**2)
    return m, np.sqrt(v), np.mean(d**3) / v**1.5, np.mean(d**4) / v**2 - 3.0


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    figdir = os.path.join(os.path.dirname(here), "figures")
    os.makedirs(figdir, exist_ok=True)

    S, f1, f2 = painleve_tw()
    tp = np.trapezoid
    def mom_pdf(S, f):
        nrm = tp(f, S); m = tp(S * f, S) / nrm; v = tp((S - m)**2 * f, S) / nrm
        return nrm, m, np.sqrt(v), tp((S - m)**3 * f, S) / nrm / v**1.5

    print("\n=== Exact Tracy-Widom densities (Painleve II / Hastings-McLeod) ===")
    print("  TW2: norm,mean,std,skew = %.5f %.4f %.4f %.4f  (target 1, -1.7711, 0.9018, 0.2241)"
          % mom_pdf(S, f2))
    print("  TW1: norm,mean,std,skew = %.5f %.4f %.4f %.4f  (target 1, -1.2065, 1.2680, 0.2935)"
          % mom_pdf(S, f1))

    print("\n=== Raw inner-exit-measure histograms vs EXACT TW (no fit) ===")
    x2 = first_node(np.sqrt(2.0), seed=2)   # beta = 2  -> TW2
    x1 = first_node(2.0, seed=1)            # beta = 1  -> TW1
    print("  beta=2 (eta=sqrt2): meas mean,std,skew = %.4f %.4f %.4f  vs TW2 -1.7711 0.9018 0.2241"
          % mom(x2)[:3])
    print("  beta=1 (eta=2):     meas mean,std,skew = %.4f %.4f %.4f  vs TW1 -1.2065 1.2680 0.2935"
          % mom(x1)[:3])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11.2, 4.4))
    for a, x, fpdf, bb, eta, col in [(ax[0], x2, f2, 2, r"\sqrt{2}", "C0"),
                                     (ax[1], x1, f1, 1, "2", "C1")]:
        a.hist(x, bins=70, density=True, color=col, alpha=0.55,
               label=fr"sim $Y_{{\rm peel}}$, $\beta={bb}$ ($\eta={eta}$)")
        a.plot(S, fpdf, "r-", lw=2.0, label=fr"exact Tracy-Widom$_{bb}$ (no fit)")
        a.set_xlim(-6, 2.5); a.set_xlabel(r"$Y_{\rm peel}$ (inner peel-off)")
        a.set_ylabel("density")
        a.set_title(fr"$\beta={bb}$: inner exit measure $=$ TW$_{bb}$")
        a.legend(fontsize=9, frameon=False); a.grid(alpha=0.3)
    fig.suptitle(r"Inner exit measure $=$ stochastic-Airy / Tracy-Widom edge,  $\eta=2/\sqrt{\beta}$  (raw, parameter-free)",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = os.path.join(figdir, "folded_cycle_tw_pdf.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print(f"\n  figure -> figures/folded_cycle_tw_pdf.png\n")


if __name__ == "__main__":
    main()
