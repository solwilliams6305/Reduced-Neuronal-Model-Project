#!/usr/bin/env python3
"""
normal_form_barrier.py

Hopf normal-form derivation of the escape barrier B(delta) for the stochastic FHN
resonator, tested against the multi-seed measured B(delta).

Theory
------
Near the (lower-branch) Hopf bifurcation the resting focus loses stability at I_H,
trace J = 0. Write the deficit delta = I_H - I. The linearization has eigenvalues
lambda = alpha +- i*omega with alpha = 1/2 trace J. We reduce the planar field to
the Hopf normal form in the complex amplitude z:

    z' = (alpha + i omega) z + c1 |z|^2 z ,    c1 = a3 + i b3.

The radial amplitude r = |z| then obeys a gradient flow r' = -U'(r) with the
amplitude potential
    U(r) = -1/2 alpha r^2 - 1/4 a3 r^4 .

For a SUBCRITICAL Hopf (a3 > 0) and a STABLE focus (alpha < 0) this is a single-well
+ barrier: a well at r = 0 and a barrier maximum at the unstable limit cycle
    r_u^2 = -alpha / a3 = |alpha| / a3 ,
with barrier height
    dU = U(r_u) - U(0) = alpha^2 / (4 a3) .

Degenerate (voltage-only) noise sigma enters the amplitude equation through the
projection of the v-direction onto the normal-form coordinate. With noise only in v
the phase-averaged radial diffusion is D_r = 1/2 * sigma^2 * |q_v|^2 * P, where q is
the (normalized) Hopf eigenvector and P=1/2 is the cos^2 phase average. Kramers'
law k ~ exp(-dU / D_r) then gives the measured-style barrier in exp(-B/sigma^2):

    B_nf = dU / (1/2 |q_v|^2 P) = alpha^2 / (4 a3) * 4 / (|q_v|^2)
         = alpha^2 / (a3 |q_v|^2) .

Since alpha(delta) ~ -k*delta is linear near the Hopf, the LEADING normal-form
prediction is  B_nf ~ delta^2.  We test that against the data.

a3 (first-Lyapunov-type cubic coefficient) is computed by the standard planar Hopf
formula (Kuznetsov) at the Hopf point, using that the only nonlinearity is -v^3/3.
"""
from __future__ import annotations
import csv, math
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

A_PAR, B_PAR, EPS = 0.7, 0.8, 0.08


def fixed_point(I):
    """Lower-branch resting fixed point: w=(v+a)/b, v - v^3/3 - w + I = 0."""
    # v - v^3/3 - (v+a)/b + I = 0  ->  -(1/3)v^3 + (1 - 1/b) v + (I - a/b) = 0
    c = [-1.0/3.0, 0.0, (1.0 - 1.0/B_PAR), (I - A_PAR/B_PAR)]
    roots = np.roots(c)
    real = sorted(r.real for r in roots if abs(r.imag) < 1e-8)
    v = real[0]                      # most-negative (resting) branch
    w = (v + A_PAR) / B_PAR
    return v, w


def jacobian(v):
    return np.array([[1.0 - v*v, -1.0],
                     [EPS,       -EPS*B_PAR]])


def alpha_omega(I):
    v, w = fixed_point(I)
    J = jacobian(v)
    ev = np.linalg.eigvals(J)
    # complex pair expected
    return float(ev[0].real), float(abs(ev[0].imag)), v, w


def hopf_current():
    """trace J = (1-v^2) - eps*b = 0 -> v_H = -sqrt(1-eps*b); solve I_H."""
    vH = -math.sqrt(1.0 - EPS*B_PAR)
    wH = (vH + A_PAR) / B_PAR
    I_H = -(vH - vH**3/3.0 - wH)      # from f=0
    return I_H, vH, wH


def lyapunov_a3(vH, omega0):
    """
    Planar Hopf first-Lyapunov cubic coefficient a3 in r' = alpha r + a3 r^3,
    using Kuznetsov's formula. Only nonlinearity: f1 = ... - v^3/3, so
      B(x,y)  has 1st comp  f_vv * x1 y1 = (-2 vH) x1 y1
      C(x,y,z)has 1st comp  f_vvv x1 y1 z1 = (-2) x1 y1 z1
    """
    A = jacobian(vH)
    w0 = omega0
    # eigenvector q: A q = i w0 q ; adjoint p: A^T p = -i w0 p ; <p,q>=1
    evals, evecs = np.linalg.eig(A)
    k = int(np.argmax(evals.imag))            # +i w0
    q = evecs[:, k]
    evalsT, evecsT = np.linalg.eig(A.T)
    kt = int(np.argmin(evalsT.imag))          # -i w0
    p = evecsT[:, kt]
    p = p / np.conj(np.vdot(p, q))            # <p,q> = conj(p).q = 1
    # forms
    def Bf(x, y):
        return np.array([(-2.0*vH) * x[0]*y[0], 0.0], dtype=complex)
    def Cf(x, y, z):
        return np.array([(-2.0) * x[0]*y[0]*z[0], 0.0], dtype=complex)
    I2 = np.eye(2, dtype=complex)
    qb = np.conj(q)
    h11 = np.linalg.solve(A, -Bf(q, qb))                  # -A^{-1} B(q,qbar)
    h20 = np.linalg.solve(2j*w0*I2 - A, Bf(q, q))         # (2i w0 - A)^{-1} B(q,q)
    g21 = np.vdot(p, Cf(q, q, qb)) + 2.0*np.vdot(p, Bf(q, h11)) + np.vdot(p, Bf(qb, h20))
    c1 = g21 / 2.0
    l1 = c1.real / w0          # first Lyapunov coefficient (sign = criticality)
    # amplitude cubic coefficient a3 = Re(c1)  (coefficient of r^3 in r')
    a3 = c1.real
    return a3, l1, q


def main():
    out = Path(__file__).resolve().parent / "results" / "normal_form_barrier"
    out.mkdir(parents=True, exist_ok=True)
    csv_in = Path(__file__).resolve().parent / "results" / "multiseed_delta" / "barriers.csv"

    I_H, vH, wH = hopf_current()
    aH, oH, _, _ = alpha_omega(I_H)
    a3, l1, qH = lyapunov_a3(vH, oH if oH > 0 else math.sqrt(abs(np.linalg.det(jacobian(vH)))))
    qv2 = abs(qH[0])**2 / (abs(qH[0])**2 + abs(qH[1])**2)   # normalized v-weight of eigenvector
    crit = "SUBcritical (a3>0): unstable limit cycle = threshold, barrier exists" if a3 > 0 \
           else "SUPERcritical (a3<0): no unstable cycle from normal form"
    print(f"Hopf: I_H={I_H:.6f}  v_H={vH:.5f}  omega0={oH:.5f}")
    print(f"first Lyapunov l1={l1:.4e}  cubic a3={a3:.4e}  -> {crit}")
    print(f"eigenvector v-weight |q_v|^2 (normalized) = {qv2:.4f}\n")

    # aggregate measured B per delta
    agg = {}
    with open(csv_in) as f:
        for r in csv.DictReader(f):
            d = float(r["delta"]); agg.setdefault(d, []).append(float(r["B"]))
    deltas = sorted(agg)
    Bm = np.array([np.mean(agg[d]) for d in deltas])
    Bs = np.array([np.std(agg[d], ddof=1) for d in deltas])
    deltas = np.array(deltas)

    # alpha(delta) numerically
    alphas = np.array([alpha_omega(I_H - d)[0] for d in deltas])

    # normal-form barrier prediction (shape ~ alpha^2). Use measured-style B:
    #   B_nf = alpha^2 / (a3 * |q_v|^2)
    if a3 > 0:
        B_nf = alphas**2 / (a3 * qv2)
    else:
        B_nf = np.full_like(alphas, np.nan)

    # log-log slopes
    def slope(x, y):
        return float(np.polyfit(np.log(x), np.log(y), 1)[0])
    s_delta = slope(deltas, Bm)
    s_alpha = slope(np.abs(alphas), Bm)
    s_alpha_delta = slope(deltas, np.abs(alphas))
    print("MEASURED scaling")
    print(f"  d ln B / d ln delta   = {s_delta:.3f}   (normal form predicts 2)")
    print(f"  d ln B / d ln|alpha|  = {s_alpha:.3f}   (normal form predicts 2)")
    print(f"  d ln|alpha|/d ln delta= {s_alpha_delta:.3f}   (near-Hopf predicts 1)\n")

    # is B/alpha^2 constant?  (normal form) -> implied a3*qv2
    implied = alphas**2 / Bm          # = a3*qv2 if normal form held
    print("  delta    alpha       B_meas      B/alpha^2 (=a3*qv2 if NF)   B_nf")
    for i, d in enumerate(deltas):
        bnf = B_nf[i] if a3 > 0 else float('nan')
        print(f"  {d:.3f}  {alphas[i]:+.5f}  {Bm[i]:.3e}   {implied[i]:.4f}            {bnf:.3e}")
    print(f"\n  normal-form constant a3*|q_v|^2 = {a3*qv2:.4f}")
    print(f"  B/alpha^2 ranges {implied.min():.3f} -> {implied.max():.3f}"
          f"  (factor {implied.max()/implied.min():.1f})")

    # plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.3))
    ax1.errorbar(deltas, Bm, yerr=Bs, fmt="o", capsize=3, label="measured B (multi-seed)")
    if a3 > 0:
        ax1.plot(deltas, B_nf, "r--", label=r"normal form $B=\alpha^2/(a_3|q_v|^2)$")
        # fit constant scale to data for shape comparison
        scale = np.sum(Bm*B_nf)/np.sum(B_nf*B_nf)
        ax1.plot(deltas, scale*B_nf, "g:", label=f"NF shape x{scale:.2f}")
    ax1.set_xlabel("Hopf deficit delta"); ax1.set_ylabel("barrier B")
    ax1.set_title("B(delta): measured vs normal form"); ax1.legend(fontsize=8)

    ax2.loglog(deltas, Bm, "o", label=f"measured (slope {s_delta:.2f})")
    ax2.loglog(deltas, Bm[0]*(deltas/deltas[0])**2, "r--", label="slope 2 (normal form)")
    ax2.loglog(deltas, Bm[0]*(deltas/deltas[0])**0.5, "b:", label="slope 0.5")
    ax2.set_xlabel("delta"); ax2.set_ylabel("B")
    ax2.set_title("log-log: B vs delta"); ax2.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out/"normal_form_barrier.png", dpi=160); plt.close(fig)

    with open(out/"normal_form_barrier.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["delta", "alpha", "B_meas", "B_std", "B_over_alpha2", "B_nf"])
        for i, d in enumerate(deltas):
            wr.writerow([d, alphas[i], Bm[i], Bs[i], implied[i],
                         B_nf[i] if a3 > 0 else ""])
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
