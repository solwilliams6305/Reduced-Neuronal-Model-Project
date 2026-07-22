#!/usr/bin/env python3
"""
instanton_exit.py

Most-probable escape path (instanton) for the resonator regime, and the phase
at which it crosses a mid radius -- to be compared with the *measured* ~-0.55 rad
commitment lead from test_commitment_gate.py / resonance_scaling.py.

Why this test
-------------
test_commitment_gate.py found a genuine, sigma-dependent ~0.5 rad lead of the
escape-commitment phase over the deterministic steepest-outward ("geometric
danger") direction, and resonance_scaling.py showed that lead does NOT scale with
the linear damping/frequency ratio kappa = |alpha|/omega. That rules out a simple
linear-resonance explanation and points at the *nonlinear* large-deviation
geometry: in the small-noise limit the system escapes along the minimiser of the
Freidlin-Wentzell action (the instanton). If the measured lead is set by that
geometry, the instanton should cross the mid shell at roughly the same phase the
ensemble does -- i.e. its predicted lead over the geometric-danger angle should be
of the same sign and rough magnitude as the measured -0.55 rad.

Model and action
----------------
2D FHN, noise ONLY in v (degenerate diffusion):
    v' = f(v,w) = v - v^3/3 - w + I
    w' = g(v,w) = eps (v + a - b w)         (no noise)
Freidlin-Wentzell rate functional for additive noise with diffusion B=diag(1,0):
    S[X] = (1/2) integral ( v' - f )^2 dt      (the w-equation must hold exactly)
Because the w-channel is noiseless, any finite-action path must satisfy
w' = g(v,w) identically; deviations there are infinitely costly. We enforce this
with a large anisotropic penalty lambda on the w-residual, which (a) makes the
constraint effectively hard and (b) gives a clean analytic gradient.

Discrete action (fixed endpoints, free interior nodes X_1..X_{M-1}):
    S = sum_m dt * (1/2) [ (vdot - f)^2 + lam * (wdot - g)^2 ]
with vdot,wdot the centred/forward finite differences. Minimised by vectorised
gradient descent; we scan a few path durations T_path and keep the lowest action.

Endpoints: start = stable spiral fixed point; end = a point on the right
(spiking) branch past the separatrix (v_end ~ 0.9). The interior is initialised
as a straight line and relaxed.

Output: instanton_exit.png (phase plane with nullclines + path), and a printed
comparison of the instanton's mid-shell crossing phase, the geometric-danger
angle, and their difference vs the measured lead.

Leading order only: this is the sigma -> 0 most-probable path. The measured lead
is at finite, near-threshold sigma, so we expect agreement in sign and rough
magnitude, not to many digits.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_resonator_phase_difference import (  # noqa: E402
    FHNParams, choose_stable_spiral_fp, make_phase_basis,
    phase_and_radius, fhn_rhs,
)
from test_commitment_gate import geometric_danger_angle  # noqa: E402


def f_v(v, w, p: FHNParams):
    return v - v ** 3 / 3.0 - w + p.I


def g_w(v, w, p: FHNParams):
    return p.eps * (v + p.a - p.b * w)


def action_and_grad(X, dt, p: FHNParams, lam: float):
    """
    X: (M+1, 2) path, rows are (v,w). Endpoints X[0], X[-1] fixed.
    Returns (S, grad) with grad shape (M+1,2); grad on endpoints is zeroed.

    Discrete action with forward differences:
       Xdot[m] = (X[m+1]-X[m])/dt,   residual r[m] = Xdot[m] - F(X[m])
       S = sum_m dt * 0.5 * ( r_v[m]^2 + lam r_w[m]^2 )
    F(X)=(f,g). Gradient via chain rule; metric Wt = diag(1,lam).
    """
    v = X[:, 0]; w = X[:, 1]
    F = np.column_stack([f_v(v, w, p), g_w(v, w, p)])      # (M+1,2)
    Xdot = (X[1:] - X[:-1]) / dt                            # (M,2)
    r = Xdot - F[:-1]                                       # residual at left node of each segment
    Wt = np.array([1.0, lam])
    S = float(np.sum(0.5 * dt * (r ** 2 * Wt).sum(axis=1)))

    M = X.shape[0] - 1
    grad = np.zeros_like(X)
    Wr = r * Wt                                             # (M,2), weighted residual per segment

    # d/dX[m] of segment-m term: residual r[m] depends on X[m] (via -F and -1/dt)
    #   dr/dX[m] = -(1/dt) I - JF(X[m])
    # and on X[m+1] (via +1/dt):  dr/dX[m+1] = +(1/dt) I
    # contribution to grad from term m (= dt * Wr[m] . dr/d.):
    for m in range(M):
        vm, wm = X[m, 0], X[m, 1]
        JF = np.array([[1.0 - vm ** 2, -1.0],
                       [p.eps,        -p.eps * p.b]])
        # wrt X[m]:
        grad[m] += dt * (Wr[m] @ (-(1.0 / dt) * np.eye(2) - JF))
        # wrt X[m+1]:
        grad[m + 1] += dt * (Wr[m] @ ((1.0 / dt) * np.eye(2)))

    grad[0] = 0.0
    grad[-1] = 0.0
    return S, grad


def minimise(X0, dt, p, lam, iters, lr):
    X = X0.copy()
    S_hist = []
    step = lr
    S_prev, _ = action_and_grad(X, dt, p, lam)
    for it in range(iters):
        S, g = action_and_grad(X, dt, p, lam)
        S_hist.append(S)
        gn = np.sqrt((g ** 2).sum())
        if gn < 1e-12:
            break
        Xtry = X - step * g
        Stry, _ = action_and_grad(Xtry, dt, p, lam)
        if Stry < S:                 # accept, gently grow step
            X = Xtry
            step *= 1.05
        else:                        # backtrack
            step *= 0.5
            if step < 1e-12:
                break
    Sf, _ = action_and_grad(X, dt, p, lam)
    S_hist.append(Sf)
    return X, S_hist


def drift(X, p: FHNParams):
    v = X[:, 0]; w = X[:, 1]
    return np.column_stack([f_v(v, w, p), g_w(v, w, p)])


def geo_action_and_grad(X, p: FHNParams, lam: float):
    """
    Parametrisation-free geometric (Freidlin-Wentzell) action with metric
    M = diag(1, lam) enforcing the noiseless w-channel:
        S = sum_i ( ||dX_i||_M * ||b_i||_M  -  dX_i^T M b_i )
    where dX_i = X[i+1]-X[i], b_i = drift(X[i]). Each term >= 0, zero where the
    segment is parallel to the drift (cost-free relaxation). Analytic gradient.
    """
    M = np.array([1.0, lam])
    b = drift(X, p)                       # (N+1,2)
    dX = X[1:] - X[:-1]                    # (N,2)
    bi = b[:-1]                            # drift at left node of each segment
    # M-norms
    np_ = np.sqrt((dX * dX * M).sum(axis=1))      # ||dX||_M
    nq_ = np.sqrt((bi * bi * M).sum(axis=1))      # ||b||_M
    dot = (dX * bi * M).sum(axis=1)               # dX^T M b
    eps_n = 1e-12
    terms = np_ * nq_ - dot
    S = float(terms.sum())

    # d term_i / d(dX_i) = (nq/np) M dX - M b
    # d term_i / d(b_i)  = (np/nq) M b  - M dX
    dterm_dp = (nq_ / (np_ + eps_n))[:, None] * (dX * M) - (bi * M)
    dterm_dq = (np_ / (nq_ + eps_n))[:, None] * (bi * M) - (dX * M)

    grad = np.zeros_like(X)
    N = X.shape[0] - 1
    # dX_i contributes -> X[i] (with -) and X[i+1] (with +)
    grad[:-1] += -dterm_dp
    grad[1:]  +=  dterm_dp
    # b_i = drift(X[i]); chain through Jacobian J(X[i])
    for i in range(N):
        vi = X[i, 0]
        J = np.array([[1.0 - vi ** 2, -1.0],
                      [p.eps,         -p.eps * p.b]])
        grad[i] += J.T @ dterm_dq[i]
    grad[0] = 0.0
    grad[-1] = 0.0
    return S, grad


def reparametrize(X, lam):
    """Redistribute interior nodes to uniform M-metric arc length."""
    M = np.array([1.0, lam])
    seg = np.sqrt(((X[1:] - X[:-1]) ** 2 * M).sum(axis=1))
    s = np.concatenate([[0.0], np.cumsum(seg)])
    if s[-1] <= 0:
        return X
    s /= s[-1]
    snew = np.linspace(0.0, 1.0, X.shape[0])
    Xv = np.interp(snew, s, X[:, 0])
    Xw = np.interp(snew, s, X[:, 1])
    out = np.column_stack([Xv, Xw])
    out[0] = X[0]; out[-1] = X[-1]
    return out


def gmam_minimise(X0, p, lam, iters, lr, reparam_every=10):
    X = reparametrize(X0.copy(), lam)
    S_hist = []
    step = lr
    for it in range(iters):
        S, g = geo_action_and_grad(X, p, lam)
        S_hist.append(S)
        gn = np.sqrt((g ** 2).sum())
        if gn < 1e-12:
            break
        Xtry = X - step * g
        Stry, _ = geo_action_and_grad(Xtry, p, lam)
        if Stry < S:
            X = Xtry
            step *= 1.05
        else:
            step *= 0.5
            if step < 1e-14:
                break
        if (it + 1) % reparam_every == 0:
            X = reparametrize(X, lam)
    X = reparametrize(X, lam)
    S_hist.append(geo_action_and_grad(X, p, lam)[0])
    return X, S_hist


def crossing_phase_at_radius(X, basis, R_m, which="last"):
    """Phase where the path crosses local radius R_m outbound.
    which='first' -> first time it reaches >=R_m; 'last' -> the final outbound
    crossing before the end (the committed leg, robust to winding)."""
    ph, rad = phase_and_radius(X[:, 0], X[:, 1], basis)
    below = np.flatnonzero(rad < R_m)
    above = np.flatnonzero(rad >= R_m)
    if above.size == 0:
        return float("nan"), float("nan")
    if which == "first" or below.size == 0:
        k = int(above[0])
    else:
        # last index that is still below R_m; the next step is the committed crossing
        k = int(below[-1])
        k = min(k + 1, len(rad) - 1)
    return float(ph[k]), float(rad[k])


def count_windings(X, basis):
    """Total unwrapped phase swept (in revolutions) along the path."""
    ph, _ = phase_and_radius(X[:, 0], X[:, 1], basis)
    return float(np.sum(np.abs(np.diff(np.unwrap(ph)))) / (2.0 * math.pi))


def wrapped(d):
    return (d + math.pi) % (2.0 * math.pi) - math.pi


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--I", type=float, default=0.30)
    ap.add_argument("--a", type=float, default=0.7)
    ap.add_argument("--b", type=float, default=0.8)
    ap.add_argument("--eps", type=float, default=0.08)
    ap.add_argument("--v-end", type=float, default=0.9,
                    help="v of target endpoint on the spiking branch.")
    ap.add_argument("--M", type=int, default=400, help="path segments")
    ap.add_argument("--T-paths", type=float, nargs="*",
                    default=[40.0, 70.0, 110.0, 160.0])
    ap.add_argument("--lam", type=float, default=2000.0)
    ap.add_argument("--iters", type=int, default=8000)
    ap.add_argument("--lr", type=float, default=2e-4)
    ap.add_argument("--mid-frac", type=float, default=0.55,
                    help="mid shell as a fraction of |fp->target| local radius.")
    ap.add_argument("--method", choices=["mam", "gmam"], default="mam")
    ap.add_argument("--spiral-init", action="store_true",
                    help="initialise the path as an outward spiral (lets it wind).")
    ap.add_argument("--n-turns", type=float, default=2.0,
                    help="number of turns in the spiral initialisation.")
    ap.add_argument("--measured-lead", type=float, default=-0.55,
                    help="ensemble commitment lead (rad) for reference line.")
    ap.add_argument("--outdir", type=str, default="results/instanton_exit")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    p = FHNParams(I=args.I, a=args.a, b=args.b, eps=args.eps)

    v_fp, w_fp, cls = choose_stable_spiral_fp(p)
    if not str(cls["kind"]).startswith("stable"):
        print(f"Fixed point is {cls['kind']}, not a stable spiral. Abort.")
        return
    basis = make_phase_basis(v_fp, w_fp, p)
    print(f"stable spiral fp=({v_fp:.4f},{w_fp:.4f})  "
          f"alpha={basis.alpha:.4f} omega={basis.omega:.4f} "
          f"kappa={abs(basis.alpha)/basis.omega:.4f}")

    start = np.array([v_fp, w_fp])
    # target on the spiking branch: put it on the v-nullcline (w = v - v^3/3 + I)
    # so it is a natural attracting point of the fast field past the fold.
    w_end = args.v_end - args.v_end ** 3 / 3.0 + p.I
    target = np.array([args.v_end, w_end])
    # local radius of the target (for choosing the mid shell)
    _, R_target = phase_and_radius(np.array([target[0]]), np.array([target[1]]), basis)
    R_mid = args.mid_frac * float(R_target[0])
    print(f"target=({target[0]:.4f},{target[1]:.4f})  "
          f"R_target={float(R_target[0]):.4f}  R_mid={R_mid:.4f}")

    # initial winding phase of the target (its angular position from the focus)
    ph_t, _ = phase_and_radius(np.array([target[0]]), np.array([target[1]]), basis)
    ph_target = float(ph_t[0])

    if args.method == "gmam":
        s = np.linspace(0.0, 1.0, args.M + 1)
        results = {}
        for tag, n_turns in (("straight", 0.0), ("spiral2", 2.0), ("spiral4", 4.0)):
            rr = s * float(R_target[0])
            ang = ph_target + 2.0 * math.pi * n_turns * (s - 1.0)
            local = np.column_stack([rr * np.cos(ang), rr * np.sin(ang)])
            X0 = basis.fp[None, :] + local @ basis.B.T
            X0[0] = start; X0[-1] = target
            X, S_hist = gmam_minimise(X0, p, args.lam, args.iters, args.lr)
            nw = count_windings(X, basis)
            ph_last, _ = crossing_phase_at_radius(X, basis, R_mid, which="last")
            gd = geometric_danger_angle(basis, p, R_mid, 0.01)
            results[tag] = (S_hist[-1], nw, ph_last, wrapped(ph_last - gd))
        print("\n  gMAM (parametrisation-free) from different initialisations:")
        print("   init      S_geo      windings   last_phase   lead(inst-geom)")
        for tag, (Sg, nw, phl, ld) in results.items():
            print(f"   {tag:8s} {Sg:.4e}   {nw:6.2f}    {phl:+7.3f}     {ld:+7.3f}")
        leads = [results[t][3] for t in results]
        spread = max(leads) - min(leads)
        gd = geometric_danger_angle(basis, p, R_mid, 0.01)
        print("\n================ gMAM VERDICT ================")
        print(f"geometric danger angle at mid shell = {gd:+.3f} rad")
        print(f"lead spread across initialisations  = {spread:.3f} rad")
        med_lead = float(np.median(leads))
        print(f"median predicted lead               = {med_lead:+.3f} rad")
        print(f"measured ensemble lead              = {args.measured_lead:+.3f} rad")
        if spread < 0.25:
            print("=> initialisation-independent: this IS the geometric instanton lead.")
        else:
            print("=> still init-dependent: not yet converged to a single minimiser.")
        if abs(wrapped(med_lead - args.measured_lead)) < 0.3:
            print("   and it matches the measured commitment lead.")
        print("=============================================")
        # plot the spiral2 result
        Xp = None
        rr = s * float(R_target[0]); ang = ph_target + 2.0 * math.pi * 2.0 * (s - 1.0)
        local = np.column_stack([rr * np.cos(ang), rr * np.sin(ang)])
        X0 = basis.fp[None, :] + local @ basis.B.T; X0[0] = start; X0[-1] = target
        Xp, Sh = gmam_minimise(X0, p, args.lam, args.iters, args.lr)
        fig, (axp, axc) = plt.subplots(1, 2, figsize=(12, 5.2))
        vv = np.linspace(-2.2, 2.2, 400)
        axp.plot(vv, vv - vv ** 3 / 3.0 + p.I, "g-", lw=1, label="v-nullcline")
        axp.plot(vv, (vv + p.a) / p.b, "m-", lw=1, label="w-nullcline")
        axp.plot(Xp[:, 0], Xp[:, 1], "b-", lw=2, label="gMAM instanton")
        axp.plot(*start, "ko", ms=7); axp.plot(*target, "rs", ms=7)
        axp.set_xlim(-2.2, 2.2); axp.set_ylim(-1.0, 1.6)
        axp.set_xlabel("v"); axp.set_ylabel("w")
        axp.set_title(f"gMAM instanton (I={p.I}, eps={p.eps})"); axp.legend(fontsize=8)
        axc.semilogy(Sh); axc.set_xlabel("iteration"); axc.set_ylabel("geometric action")
        axc.set_title("gMAM convergence")
        fig.tight_layout(); fig.savefig(outdir / "instanton_gmam.png", dpi=170); plt.close(fig)
        print(f"\nOutputs -> {outdir}")
        return

    best = None
    print("\n  per-duration: windings, first vs LAST outbound crossing of mid shell")
    print("   T_path   dt     S_final     windings   first_phase   last_phase")
    for Tp in args.T_paths:
        dt = Tp / args.M
        s = np.linspace(0.0, 1.0, args.M + 1)[:, None]
        if args.spiral_init:
            # spiral outward: radius grows 0->R_target, phase sweeps n_turns then
            # lands on the target's angular position.
            rr = s[:, 0] * float(R_target[0])
            ang = ph_target + 2.0 * math.pi * args.n_turns * (s[:, 0] - 1.0)
            local = np.column_stack([rr * np.cos(ang), rr * np.sin(ang)])
            X0 = basis.fp[None, :] + local @ basis.B.T
            X0[0] = start; X0[-1] = target
        else:
            X0 = (1 - s) * start[None, :] + s * target[None, :]
            X0[1:-1, 1] += 0.02 * np.sin(np.pi * s[1:-1, 0])
        X, S_hist = minimise(X0, dt, p, args.lam, args.iters, args.lr)
        S_final = S_hist[-1]
        nw = count_windings(X, basis)
        ph_first, _ = crossing_phase_at_radius(X, basis, R_mid, which="first")
        ph_last, r_last = crossing_phase_at_radius(X, basis, R_mid, which="last")
        print(f"   {Tp:6.1f} {dt:.4f}  {S_final:.4e}   {nw:6.2f}    {ph_first:+7.3f}     {ph_last:+7.3f}")
        if best is None or S_final < best["S"]:
            best = {"S": S_final, "X": X, "S_hist": S_hist, "Tp": Tp, "dt": dt,
                    "ph_cross": ph_last, "r_cross": r_last}

    X = best["X"]
    # Compare instanton vs geometric danger at several radii spanning the
    # ensemble's measurement band (ensemble used ~0.15-0.35 of threshold radius).
    print("\n  radius-resolved crossing (instanton vs geometric danger):")
    print("   R_frac   R       inst_phase   geom_phase   lead(inst-geom)")
    for rf in (0.12, 0.18, 0.25, 0.35, 0.55):
        Rr = rf * float(R_target[0])
        ph_i, _ = crossing_phase_at_radius(X, basis, Rr)
        gd = geometric_danger_angle(basis, p, Rr, best["dt"])
        print(f"   {rf:4.2f}   {Rr:5.3f}   {ph_i:+7.3f}     {gd:+7.3f}     {wrapped(ph_i-gd):+7.3f}")

    geom = geometric_danger_angle(basis, p, R_mid, best["dt"])
    pred_lead = wrapped(best["ph_cross"] - geom)
    print("\n================ INSTANTON EXIT VERDICT ================")
    print(f"best T_path = {best['Tp']:.1f}   action S = {best['S']:.4e}")
    print(f"instanton crosses mid shell (R={R_mid:.4f}) at phase = {best['ph_cross']:+.3f} rad")
    print(f"geometric danger angle at that shell             = {geom:+.3f} rad")
    print(f"=> predicted instanton lead over geometry        = {pred_lead:+.3f} rad")
    print(f"   measured ensemble commitment lead             = {args.measured_lead:+.3f} rad")
    same_sign = (pred_lead < 0) == (args.measured_lead < 0)
    close = abs(wrapped(pred_lead - args.measured_lead)) < 0.4
    if same_sign and close:
        print("   MATCH (sign + magnitude): the lead is set by the large-deviation")
        print("   instanton geometry, not by linear resonance.")
    elif same_sign:
        print("   SAME SIGN but magnitude differs: instanton geometry is the right")
        print("   mechanism; finite-sigma / leading-order corrections account for the gap.")
    else:
        print("   sign disagrees: the leading-order instanton does not explain the")
        print("   measured lead at this operating point (revisit endpoints / sigma).")
    print("========================================================")

    # ---- plot ----
    fig, (axp, axc) = plt.subplots(1, 2, figsize=(12, 5.2))

    vv = np.linspace(-2.2, 2.2, 400)
    axp.plot(vv, vv - vv ** 3 / 3.0 + p.I, "g-", lw=1, label="v-nullcline")
    axp.plot(vv, (vv + p.a) / p.b, "m-", lw=1, label="w-nullcline")
    axp.plot(X[:, 0], X[:, 1], "b-", lw=2, label="instanton")
    axp.plot(*start, "ko", ms=7, label="spiral fp")
    axp.plot(*target, "rs", ms=7, label="target")
    # mid shell circle in (v,w)
    phis = np.linspace(-math.pi, math.pi, 200)
    ring = np.array([basis.fp + basis.B @ np.array([R_mid * math.cos(t), R_mid * math.sin(t)])
                     for t in phis])
    axp.plot(ring[:, 0], ring[:, 1], "k:", lw=0.8, label="mid shell")
    axp.plot(*basis.fp + basis.B @ np.array([R_mid * math.cos(best["ph_cross"]),
                                             R_mid * math.sin(best["ph_cross"])]),
             "b^", ms=9)
    axp.plot(*basis.fp + basis.B @ np.array([R_mid * math.cos(geom),
                                             R_mid * math.sin(geom)]),
             "rv", ms=9, label="geom danger")
    axp.set_xlim(-2.2, 2.2); axp.set_ylim(-1.0, 1.6)
    axp.set_xlabel("v"); axp.set_ylabel("w")
    axp.set_title(f"Instanton exit path (I={p.I}, eps={p.eps})")
    axp.legend(fontsize=8, loc="upper left")

    axc.semilogy(best["S_hist"])
    axc.set_xlabel("descent iteration"); axc.set_ylabel("action S")
    axc.set_title(f"action convergence (best T_path={best['Tp']:.0f})")
    fig.tight_layout()
    fig.savefig(outdir / "instanton_exit.png", dpi=170)
    plt.close(fig)
    print(f"\nOutputs -> {outdir}")


if __name__ == "__main__":
    main()
