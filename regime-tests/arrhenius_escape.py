#!/usr/bin/env python3
"""
arrhenius_escape.py

Measure the escape-rate law k(sigma) for the resonator focus and extract the
effective barrier B in  k(sigma) ~ A * exp( -B / sigma^2 )  (Kramers/Arrhenius).
This is the load-bearing probabilistic quantity: B = quasipotential barrier ΔV.
We MEASURE it from survival curves (not from the instanton solver, whose absolute
action did not converge), then compare to the analytic frozen-w (eps->0) barrier.

Method
------
For each sigma we run an ensemble, record first-passage (spike) times, and build
the survival curve S(t) = fraction not yet escaped. After the initial transient
the cloud reaches a quasi-stationary state and escapes as a Poisson process, so
ln S(t) is linear with slope -k. We fit k in the window where S in [S_hi, S_lo].
Then Arrhenius: ln k = ln A - B / sigma^2  ->  slope of ln k vs 1/sigma^2 is -B.

Frozen-w barrier
----------------
With w held at the fixed-point value w*, the v-dynamics is gradient,
dv = -dPhi/dv dt + sigma dW, Phi(v) = v^4/12 - v^2/2 - (I - w*) v. The 1D barrier
is Phi(v_saddle) - Phi(v_well). The FW prefactor for a 1D double well gives
k ~ exp(-2 dPhi / sigma^2), i.e. B_pred = 2 * (Phi_saddle - Phi_well).

Outputs: arrhenius_escape.csv, arrhenius_escape.png (survival + Arrhenius fit).
"""
from __future__ import annotations
import argparse, csv, math
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_resonator_phase_difference import (
    FHNParams, choose_stable_spiral_fp, make_phase_basis, simulate_many,
)


def survival_rate(spike_time, spiked, t, S_hi, S_lo):
    """Fit Poisson rate k from ln S(t) over the window S in [S_lo, S_hi]."""
    N = spiked.size
    st = spike_time[spiked]
    if st.size < 20:
        return None
    # survival on the time grid
    esc_by_t = np.searchsorted(np.sort(st), t, side="right")  # cumulative escapes
    S = 1.0 - esc_by_t / N
    mask = (S <= S_hi) & (S >= S_lo) & (t > 0)
    if mask.sum() < 8:
        return None
    tt = t[mask]; lnS = np.log(S[mask])
    A = np.vstack([tt, np.ones_like(tt)]).T
    slope, intercept = np.linalg.lstsq(A, lnS, rcond=None)[0]
    k = -float(slope)
    # R^2 of the fit
    pred = A @ np.array([slope, intercept])
    ss_res = float(((lnS - pred) ** 2).sum())
    ss_tot = float(((lnS - lnS.mean()) ** 2).sum()) + 1e-30
    r2 = 1.0 - ss_res / ss_tot
    return {"k": k, "r2": r2, "n_esc": int(st.size), "S_grid": S}


def frozen_w_barrier(p: FHNParams, w_star: float):
    """1D double-well barrier at fixed w = w_star. Returns (B_pred, roots)."""
    # f(v)=0: v - v^3/3 - w + I = 0  ->  v^3 - 3 v + 3(w - I) = 0
    coeffs = [1.0, 0.0, -3.0, 3.0 * (w_star - p.I)]
    roots = np.roots(coeffs)
    real = sorted(r.real for r in roots if abs(r.imag) < 1e-7)
    if len(real) < 3:
        return None, real
    v_well = real[0]              # left well (resting fixed point branch)
    v_saddle = real[1]            # middle (unstable) root
    Phi = lambda v: v ** 4 / 12.0 - v ** 2 / 2.0 - (p.I - w_star) * v
    dPhi = Phi(v_saddle) - Phi(v_well)
    return 2.0 * dPhi, (v_well, v_saddle, real[2])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--I", type=float, default=0.30)
    ap.add_argument("--a", type=float, default=0.7)
    ap.add_argument("--b", type=float, default=0.8)
    ap.add_argument("--eps", type=float, default=0.08)
    ap.add_argument("--sigmas", type=float, nargs="*",
                    default=[0.035, 0.04, 0.045, 0.05, 0.06, 0.07])
    ap.add_argument("--N", type=int, default=400)
    ap.add_argument("--T", type=float, default=250.0)
    ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--r0", type=float, default=0.06)
    ap.add_argument("--S-hi", type=float, default=0.85)
    ap.add_argument("--S-lo", type=float, default=0.10)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--outdir", type=str, default="results/arrhenius_escape")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    p = FHNParams(I=args.I, a=args.a, b=args.b, eps=args.eps)
    v_fp, w_fp, cls = choose_stable_spiral_fp(p)
    basis = make_phase_basis(v_fp, w_fp, p)
    B_pred, roots = frozen_w_barrier(p, w_fp)
    print(f"canonical point I={p.I} eps={p.eps}  fp=({v_fp:.4f},{w_fp:.4f})")
    print(f"alpha={basis.alpha:.4f} omega={basis.omega:.4f}")
    if B_pred is not None:
        print(f"frozen-w roots v=({roots[0]:.3f},{roots[1]:.3f},{roots[2]:.3f})")
        print(f"frozen-w predicted barrier B_pred = 2*dPhi = {B_pred:.4e}\n")

    rows = []
    surv = {}
    print("  sigma   n_esc   k(rate)     R^2")
    for s in sorted(args.sigmas):
        sim = simulate_many(p=p, basis=basis, sigma=s, N=args.N, T=args.T,
                            dt=args.dt, r0=args.r0, phi0=0.0,
                            spike_threshold=1.0, seed=args.seed)
        res = survival_rate(sim["spike_time"], sim["spiked"], sim["t"],
                            args.S_hi, args.S_lo)
        if res is None:
            print(f"  {s:.3f}    --     insufficient escapes in window, skip")
            continue
        rows.append({"sigma": s, "inv_sig2": 1.0 / s ** 2,
                     "k": res["k"], "ln_k": math.log(res["k"]),
                     "r2": res["r2"], "n_esc": res["n_esc"]})
        surv[s] = (sim["t"], res["S_grid"])
        print(f"  {s:.3f}   {res['n_esc']:4d}   {res['k']:.4e}   {res['r2']:.3f}")

    if len(rows) < 3:
        print("\nToo few rate points for an Arrhenius fit."); return
    with open(outdir / "arrhenius_escape.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wr.writeheader(); wr.writerows(rows)

    x = np.array([r["inv_sig2"] for r in rows])
    y = np.array([r["ln_k"] for r in rows])
    A = np.vstack([x, np.ones_like(x)]).T
    slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
    B_meas = -float(slope)
    A_pref = math.exp(intercept)
    pred = A @ np.array([slope, intercept])
    r2 = 1.0 - ((y - pred) ** 2).sum() / (((y - y.mean()) ** 2).sum() + 1e-30)

    print("\n================ ARRHENIUS VERDICT ================")
    print(f"fit  ln k = ln A - B / sigma^2")
    print(f"  measured barrier B = {B_meas:.4e}")
    print(f"  prefactor       A = {A_pref:.4e}")
    print(f"  fit R^2           = {r2:.4f}")
    if B_pred is not None:
        print(f"  frozen-w (eps->0) prediction B_pred = {B_pred:.4e}")
        print(f"  ratio B_meas / B_pred = {B_meas / B_pred:.3f}")
    if r2 > 0.95:
        print("=> escape rate follows a clean Arrhenius law in 1/sigma^2:")
        print("   k(sigma) ~ A exp(-B/sigma^2). The barrier B is well defined and")
        print("   the sigma-dependence of escape is fully captured by it.")
    else:
        print("=> Arrhenius fit is not clean; rate law may not be simple exp(-B/sigma^2)")
        print("   over this range (near-Hopf shallow barrier / pre-asymptotic).")
    print("===================================================")

    fig, (axs, axa) = plt.subplots(1, 2, figsize=(12, 5.2))
    for s, (t, S) in surv.items():
        axs.semilogy(t, np.clip(S, 1e-3, 1), label=f"sigma={s}")
    axs.set_xlabel("time"); axs.set_ylabel("survival S(t)")
    axs.set_title("First-passage survival curves"); axs.legend(fontsize=8)

    axa.plot(x, y, "o", label="measured ln k")
    xs = np.linspace(x.min(), x.max(), 50)
    axa.plot(xs, slope * xs + intercept, "r-",
             label=f"fit: B={B_meas:.3f}, R^2={r2:.3f}")
    axa.set_xlabel(r"$1/\sigma^2$"); axa.set_ylabel("ln k")
    axa.set_title("Arrhenius plot: ln k vs 1/sigma^2"); axa.legend(fontsize=9)
    fig.tight_layout(); fig.savefig(outdir / "arrhenius_escape.png", dpi=170)
    plt.close(fig)
    print(f"\nOutputs -> {outdir}")


if __name__ == "__main__":
    main()
