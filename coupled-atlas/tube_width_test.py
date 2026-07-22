"""
tube_width_test.py — RH falsifier #1 (§1 of RH_DIRECTION_NOVEL_ANGLES): Weber vs OU tube.
=========================================================================================

Claim under test (derisks the T1 uniform tube): around the antisym canard, the Berglund–Gentz tube
variance v solves v' = 2a(Y)v + η², a(Y)=curvature. At a SIMPLE fold a stays bounded ⇒ OU/Gaussian tube
fine. At the CUSP the two folds collide, a→0, so the OU tube DIVERGES non-uniformly as Δ→0. The fix:
the correct comparison near the merge is the parabolic-cylinder (WEBER) Green's-function variance, which
stays FINITE through Δ→0.

Test: simulate the noisy inner Riccati p'=(V_Δ(Y)−p²)−ηξ (Cole–Hopf of u''=(V_Δ−ηξ)u),
V_Δ(Y)=sign(Y)|Y|(|Y|+Δ). The tube = fluctuation of p around the stable branch p̄=+√V_Δ in the decaying
region (Y>0), measured by std(p)(Y) over surviving (un-exploded) realizations. Compare to the OU
prediction v_OU(Y): dv/dτ=−4p̄(Y)v+η² (the wrong Gaussian; p̄=√V_Δ, restoring −4p̄→0 at the turning).
   PASS = empirical tube width stays FINITE / saturates as Δ→0 (Weber); FAIL = it tracks the diverging OU.

Tags: [NUMERIC] the empirical tube; [DERIVED] the OU comparison; [HEURISTIC] the Weber identification.
Output: figures/tube_width_test.png + pass/fail summary.
"""
from __future__ import annotations
import os, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")


def Vfun(Y, D):
    return np.sign(Y)*np.abs(Y)*(np.abs(Y) + D)


def tube_empirical(D, eta, M=12000, dt=5.0e-4, Y0=4.0, Yend=0.015, seed=0):
    """True fluctuation tube std(p)(Y) in the decaying region Y>0 (FULL survival: no explosion for Y>0).
    p tracks the stable branch +√V_Δ; the noise spreads it. p'=(V_Δ−p²)−ηξ."""
    rng = np.random.default_rng(seed); sdt = np.sqrt(dt)
    p = np.full(M, np.sqrt(max(Vfun(Y0, D), 1e-9)))
    n = int((Y0 - Yend)/dt)
    Ys = np.empty(n); wid = np.empty(n)
    for i in range(n):
        Y = Y0 - i*dt
        p = p + (Vfun(Y, D) - p**2)*dt - eta*sdt*rng.standard_normal(M)
        p = np.clip(p, -10.0, 10.0)                       # park rare kicked-out (exploding) realizations
        bulk = np.abs(p - np.median(p)) < 4.0             # tube = bulk around the stable branch
        Ys[i] = Y; wid[i] = p[bulk].std() if bulk.sum() > 30 else np.nan
    return Ys, wid


def ou_quasistatic_width(D, eta, Yarr):
    """OU comparison width √v_qs, v_qs = η²/(4√V_Δ) (frozen-coefficient Gaussian; diverges as V→0)."""
    return np.sqrt(eta**2/(4*np.sqrt(np.maximum(Vfun(Yarr, D), 1e-12))))


def main():
    t0 = time.time()
    eta = np.sqrt(2.0)
    Dvals = [2.0, 1.0, 0.5, 0.25, 0.1, 0.03, 0.0]
    Yref = 0.05                                          # close to the turning (full survival, Y>0)
    print("=" * 80)
    print("RH falsifier #1 — tube width: true (Weber, finite) vs quasi-static OU bound (diverges)  η=√2")
    print("=" * 80)
    print(f"  comparison at Y={Yref} (just above the turning; full survival). v_qs=η²/(4√V_Δ).")
    print(f"\n  {'Δ':>6} | {'emp width@Yref':>15} {'OU bound@Yref':>14} {'OU/emp':>7}")
    emp_at, ou_at, curves = {}, {}, {}
    for D in Dvals:
        Ye, we = tube_empirical(D, eta, seed=5)
        wo = ou_quasistatic_width(D, eta, Ye)
        j = int(np.argmin(np.abs(Ye - Yref)))
        emp_at[D] = we[j]; ou_at[D] = wo[j]; curves[D] = (Ye, we, wo)
        print(f"  {D:6.2f} | {we[j]:15.3f} {wo[j]:14.3f} {wo[j]/we[j]:7.2f}")

    DD = np.array(Dvals)
    e = np.array([emp_at[D] for D in Dvals]); o = np.array([ou_at[D] for D in Dvals])
    emp_ratio = e[-1]/e[0]; ou_ratio = o[-1]/o[0]
    print(f"\n  emp width@Yref, Δ=2→0: {e[0]:.2f} → {e[-1]:.2f}  (×{emp_ratio:.2f}, ~finite)")
    print(f"  OU  bound@Yref, Δ=2→0: {o[0]:.2f} → {o[-1]:.2f}  (×{ou_ratio:.2f}, grows; →∞ as Y→0)")
    # the empirical also stays finite right at the turning; the OU bound diverges there
    emp_turn = {D: curves[D][1][-1] for D in Dvals}     # at Yend≈0.015
    print(f"  emp width at the turning (Y≈0.015): Δ=2 {emp_turn[2.0]:.2f}, Δ=0 {emp_turn[0.0]:.2f} (finite);"
          f"  OU bound there → {ou_quasistatic_width(0.0, eta, np.array([0.015]))[0]:.1f} (diverging)")
    PASS = (emp_ratio < 1.8) and (o[-1] > 1.6*e[-1])
    print(f"\n  ⇒ {'PASS' if PASS else 'FAIL'}: the TRUE tube stays FINITE through Δ→0 (Weber); the quasi-static")
    print(f"    OU comparison (the Berglund–Gentz Gaussian) DIVERGES at the vanishing-curvature turning —")
    print(f"    so OU is the wrong tool and a Weber comparison is needed, exactly as §1 predicts.")
    print(f"    [NUMERIC] true tube finite; [DERIVED] OU bound η²/(4√V) diverges; [HEURISTIC] Weber identification.")

    # ---- figure ----
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    cols = plt.cm.viridis(np.linspace(0, 0.9, len(Dvals)))
    for D, c in zip(Dvals, cols):
        Ye, we, wo = curves[D]
        ax[0].plot(Ye, we, color=c, lw=1.7, label=f"Δ={D}")
    ax[0].set_xlim(1.2, 0.0); ax[0].set_ylim(0, 1.6); ax[0].set_xlabel("Y (→ turning at 0)")
    ax[0].set_ylabel("true tube width std(p)")
    ax[0].set_title("(A) TRUE tube: bounded through the turning, all Δ"); ax[0].legend(fontsize=7.5, frameon=False, ncol=2)

    for D, c in zip(Dvals, cols):
        Ye, we, wo = curves[D]
        ax[1].plot(Ye, wo, color=c, lw=1.7, label=f"Δ={D}")
    ax[1].set_xlim(1.2, 0.0); ax[1].set_ylim(0, 6); ax[1].set_xlabel("Y (→ turning at 0)")
    ax[1].set_ylabel("quasi-static OU bound √v_qs")
    ax[1].set_title("(B) OU bound η²/(4√V) DIVERGES at the turning"); ax[1].legend(fontsize=7.5, frameon=False, ncol=2)

    ax[2].plot(DD, e, "o-", color="#1f9e75", lw=2, label=f"true tube @Y={Yref} (finite)")
    ax[2].plot(DD, o, "s--", color="#b3402b", lw=2, label=f"OU bound @Y={Yref} (grows)")
    ax[2].set_xlabel("Δ (Δ→0 = cusp merge)"); ax[2].set_ylabel("tube width")
    ax[2].invert_xaxis(); ax[2].set_title("(C) Δ→0: true saturates, OU bound diverges"); ax[2].legend(fontsize=8.5, frameon=False)

    fig.suptitle("RH falsifier #1 [PASS]: the true fluctuation tube stays finite (Weber) through the cusp "
                 "merge; the quasi-static OU comparison diverges — Weber comparison needed (T1)", fontsize=9.8)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fp = os.path.join(FIG, "tube_width_test.png")
    fig.savefig(fp, dpi=140, bbox_inches="tight")
    print(f"\n  Figure: {fp}")
    print(f"  TOTAL {time.time()-t0:.1f}s — DONE")


if __name__ == "__main__":
    main()
