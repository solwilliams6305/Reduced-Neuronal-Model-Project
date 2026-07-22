"""
hopf_analysis.py
----------------
Settles the tonic Hopf transition of the FitzHugh-Nagumo model "for good".

Strategy (analytic + numerical, no anchoring to a guessed I_c):

  1. For each I, find ALL fixed points (roots of the cubic) and classify each
     via the Jacobian  J = [[1-v^2, -1],[eps, -eps*b]]:
         trace = 1 - v^2 - eps*b ,   det = eps*(1 - b + b*v^2) > 0 (always).
  2. A Hopf needs trace = 0  ⟺  v* = ±sqrt(1 - eps*b).  The fixed point crosses
     BOTH roots as I increases, giving two Hopf currents I_H1 < I_H2 that
     bracket the tonic-spiking window.
  3. Hopf frequency  ω = sqrt(det) = sqrt(eps*(1 - eps*b^2))  ⇒  nascent period
     2π/ω is FINITE  ⇒  Hopf, not SNIC.
  4. Confirm by deterministic integration: sweep I and record limit-cycle
     amplitude & period; the cycle should appear at I_H1 with finite period and
     vanish at I_H2.

Run:  python hopf_analysis.py
"""
from __future__ import annotations
import numpy as np

try:
    from scipy.integrate import solve_ivp
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False

a, b = 0.7, 0.8
eps = 0.08

# ---------------------------------------------------------------------------
def fixed_points(I):
    # v - v^3/3 - (v+a)/b + I = 0  →  -1/3 v^3 + (1 - 1/b) v + (I - a/b) = 0
    coeffs = [-(1.0 / 3.0), 0.0, (1.0 - 1.0 / b), (I - a / b)]
    roots = np.roots(coeffs)
    return sorted(r.real for r in roots if abs(r.imag) < 1e-9)

def classify(v):
    tr = 1 - v * v - eps * b
    det = eps * (1 - b + b * v * v)
    disc = tr * tr - 4 * det
    if det < 0:
        kind = "saddle"
    else:
        osc = disc < 0
        if tr < 0:
            kind = "stable spiral" if osc else "stable node"
        elif tr > 0:
            kind = "unstable spiral" if osc else "unstable node"
        else:
            kind = "center"
    return tr, det, kind

def I_at_v(v):
    return (v + a) / b - v + v ** 3 / 3.0

# ---------------------------------------------------------------------------
def main():
    print(f"\n=== FHN Hopf analysis  (a={a}, b={b}, eps={eps}) ===\n")

    I_fold_l = I_at_v(-1.0)
    I_fold_r = I_at_v(+1.0)
    vH = np.sqrt(1 - eps * b)
    I_H1 = I_at_v(-vH)
    I_H2 = I_at_v(+vH)
    omega = np.sqrt(eps * (1 - eps * b ** 2))
    print("Analytic landmarks:")
    print(f"  left fold  (v=-1):  I = {I_fold_l:.4f}")
    print(f"  Hopf lower (v=-{vH:.4f}):  I_H1 = {I_H1:.4f}   <- tonic BORN")
    print(f"  Hopf upper (v=+{vH:.4f}):  I_H2 = {I_H2:.4f}   <- tonic DIES")
    print(f"  right fold (v=+1):  I = {I_fold_r:.4f}")
    print(f"  Hopf frequency  ω = {omega:.4f}   nascent period 2π/ω = {2*np.pi/omega:.2f}\n")

    print("Stability scan:")
    print("  I      | fixed points (v) and stability")
    for I in np.arange(0.20, 1.55, 0.05):
        parts = []
        for v in fixed_points(I):
            tr, det, kind = classify(v)
            parts.append(f"v={v:+.3f}[{kind}]")
        print(f"  {I:5.3f} | " + "   ".join(parts))

    if not HAVE_SCIPY:
        print("\n(scipy not available — skipping deterministic cycle confirmation)")
        return

    print("\nDeterministic limit-cycle confirmation (amplitude & period):")
    def rhs(t, y, I):
        v, w = y
        return [v - v ** 3 / 3 - w + I, eps * (v + a - b * w)]

    for I in [0.30, 0.33, 0.35, 0.50, 1.0, 1.40, 1.42, 1.45]:
        sol = solve_ivp(rhs, [0, 4000], [-1.0, -0.5], args=(I,),
                        max_step=0.05, rtol=1e-8, atol=1e-10, dense_output=True)
        t = np.linspace(2000, 4000, 200000)   # discard transient
        v = sol.sol(t)[0]
        amp = v.max() - v.min()
        # period from upward crossings of mean
        thr = 0.5 * (v.max() + v.min())
        cr = np.where((v[:-1] < thr) & (v[1:] >= thr))[0]
        per = np.mean(np.diff(t[cr])) if len(cr) > 2 else np.nan
        tag = "TONIC" if amp > 0.5 else "rest/spiral"
        print(f"  I={I:4.2f}  amplitude={amp:5.2f}  period={per if not np.isnan(per) else float('nan'):7.2f}  -> {tag}")

if __name__ == "__main__":
    main()
