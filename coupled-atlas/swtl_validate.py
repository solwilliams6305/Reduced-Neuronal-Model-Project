"""
Does the 8-coefficient trans-series reproduce W_swallowtail ACROSS x, or only miss at one point?

Two jobs:
  (1) Replace the eta-sweep ground truth, which was generated with FIRST-ORDER UPWIND and is
      therefore only trustworthy near beta=2 (its numerical diffusion |drift|*dp/2 is
      eta-independent and swamps the physical D=eta^2/2 at small eta).  The MUSCL/van-Leer solver
      in swtl_nonpert.py removes that; here it is re-validated at q=3 as well as q=2.
  (2) Overlay the resummation on that curve.  We know f(2)=0.2177 against truth 0.1655 with the
      8-coefficient ladder (versus 0.0354, i.e. 79% LOW, with 7).  A single point cannot say
      whether the residual error is systematic or erratic -- that is what the paper's validation
      claim needs, so sweep x.

VALIDATION GATES (both must pass before the sweep is believed):
  q=2, x=0.09: converged cusp partial sum          = 0.144770
  q=3, x=0.09: converged swallowtail partial sum   = 0.056404   (terms below 1e-6 by k=7)
At x=0.09 the non-perturbative scale is e^{-A cos(theta)/x} ~ e^{-13} , utterly negligible, so the
truncated series IS the exact answer there to ~1e-6 and makes a clean absolute reference.

Run:  python3 swtl_validate.py
"""
import numpy as np

import swtl_borel as S
import swtl_nonpert as NP
import fp_cusp as FP

CUSP_V = S.CUSP_V
CUSP_REF_X009 = 0.144770


def series_at(v, x):
    return sum(c * x ** k for k, c in enumerate(v))


def f_exact(q, x, dp=0.02):
    return NP.f_of_x(q, x, dp=dp, order=2)


def gates(v_swtl):
    print("=" * 78)
    print("GATES: second-order solver against converged partial sums (x=0.09)")
    print("=" * 78)
    ok = True
    f2 = f_exact(2, 0.09)
    print(f"  q=2: MUSCL {f2:.6f}   series {CUSP_REF_X009:.6f}   err {(f2/CUSP_REF_X009-1)*100:+.2f}%")
    ok &= abs(f2 / CUSP_REF_X009 - 1) < 0.02
    ref3 = series_at(v_swtl, 0.09)
    f3 = f_exact(3, 0.09)
    print(f"  q=3: MUSCL {f3:.6f}   series {ref3:.6f}   err {(f3/ref3-1)*100:+.2f}%")
    ok &= abs(f3 / ref3 - 1) < 0.03
    print(f"  => {'BOTH PASS' if ok else 'FAIL -- sweep not trustworthy'}")
    return ok


def sweep(v8):
    v7 = v8[:7]
    K8 = len(v8) - 1                      # 7 -> diagonal [3/4]
    K7 = len(v7) - 1                      # 6 -> diagonal [3/3]
    xs = [0.09, 0.16, 0.25, 0.36, 0.49, 0.64, 0.81, 1.00, 1.44, 2.00, 2.50]
    print("\n" + "=" * 78)
    print("RESUMMATION vs GROUND TRUTH ACROSS x   (x = eta^2 = 4/beta)")
    print("=" * 78)
    print(f"  {'x':>6} {'beta':>6} {'f_exact':>10} {'7-coef':>9} {'err':>8} {'8-coef':>9} {'err':>8} "
          f"{'naive':>11}")
    rows = []
    for x in xs:
        fe = f_exact(3, x)
        try:
            f7, _ = S.median_borel(v7, x, L=K7 - K7 // 2, M=K7 // 2, phi_deg=35.0)
        except Exception:
            f7 = np.nan
        try:
            f8, _ = S.median_borel(v8, x, L=K8 - K8 // 2, M=K8 // 2, phi_deg=35.0)
        except Exception:
            f8 = np.nan
        nav = series_at(v8, x)
        rows.append((x, fe, f7, f8))
        print(f"  {x:>6.2f} {4/x:>6.2f} {fe:>10.6f} {f7:>9.4f} {(f7/fe-1)*100:>7.1f}% "
              f"{f8:>9.4f} {(f8/fe-1)*100:>7.1f}% {nav:>11.3g}")
    print("\n  Reading: if the 8-coefficient error grows smoothly with x, the representation is")
    print("  sound and the Pade is simply starved outside the disc (a radius-of-validity")
    print("  statement).  If it is erratic in sign/size, the approximant is unreliable throughout.")
    e8 = np.array([(r[3] / r[1] - 1) for r in rows])
    e7 = np.array([(r[2] / r[1] - 1) for r in rows])
    xs_a = np.array([r[0] for r in rows])
    mono8 = np.all(np.diff(e8) > -0.02)
    print(f"\n  8-coef relative error: min {e8.min()*100:+.1f}%  max {e8.max()*100:+.1f}%  "
          f"monotone-in-x: {mono8}")
    print(f"  7-coef relative error: min {e7.min()*100:+.1f}%  max {e7.max()*100:+.1f}%")
    good8 = xs_a[np.abs(e8) < 0.10]
    good7 = xs_a[np.abs(e7) < 0.10]
    print(f"  |err| < 10%:  8-coef for x <= {good8.max() if len(good8) else float('nan'):.2f}"
          f"   (7-coef for x <= {good7.max() if len(good7) else float('nan'):.2f})")
    return rows


if __name__ == "__main__":
    v, band, v8 = S.load_ladder(verbose=False)
    print(f"ladder v0..v{len(v8)-1} = {[round(c,4) for c in v8]}\n")
    if gates(v8):
        sweep(v8)
