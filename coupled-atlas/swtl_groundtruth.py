"""
Independent ground truth Var(Y*) for the SWALLOWTAIL (q=3) from the MC-free FP-PDE.

Why: SWALLOWTAIL_TRANSSERIES_NOTES.md asserts the beta=2 target Var=0.328, but that number has no
traceable computation anywhere in the repo (the only occurrence is the assertion itself).  Since it
is the target the median-Borel resummation is graded against, it needs to be reproduced before any
"resummation fails / succeeds" verdict means anything.

Method: exactly fp_cusp.solve_fp, with the two q-dependent points generalized --
    V = sign(Y)|Y|^q            (was sign(Y) Y^2)
    p0 = sqrt(V(Y0)) = Y0^{q/2} (recessive incoming Riccati value)
Regression gate: at q=2 this must reproduce the cusp ground truth Var(beta=2) ~ 0.237*2 = 0.474
(f = Var/eta^2 = 0.237, skew ~ +0.607, exk ~ -0.237) from _w_groundtruth.json.

Run:  python3 swtl_groundtruth.py [q] [eta ...]
"""
import sys
import numpy as np
import fp_cusp as FP


def solve_fp_q(q, eta=np.sqrt(2.0), Y0=3.0, pmin=-9.0, pmax=9.0, dp=0.02, dt=5.0e-5, tau_max=7.5):
    D = eta * eta / 2.0
    pc = np.arange(pmin + dp / 2, pmax, dp)
    pf = np.arange(pmin, pmax + dp / 2, dp); pf2 = pf * pf
    p0 = np.sqrt(max(np.sign(Y0) * abs(Y0) ** q, 1e-9))
    rho = np.exp(-0.5 * ((pc - p0) / 0.30) ** 2); rho /= rho.sum() * dp
    nst = int(round(tau_max / dt))
    taus = np.empty(nst + 1); S = np.empty(nst + 1)
    taus[0] = 0.0; S[0] = rho.sum() * dp
    for k in range(nst):
        tau = k * dt; Y = Y0 - tau
        V = np.sign(Y) * abs(Y) ** q
        drift = V - pf2
        rho_left = np.concatenate(([0.0], rho))
        rho_right = np.concatenate((rho, [0.0]))
        adv = np.where(drift > 0, rho_left, rho_right) * drift
        diff = -D * (rho_right - rho_left) / dp
        J = adv + diff
        J[-1] = 0.0
        rho = rho - dt * (J[1:] - J[:-1]) / dp
        np.maximum(rho, 0.0, out=rho)
        taus[k + 1] = tau + dt; S[k + 1] = rho.sum() * dp
    y = Y0 - taus
    return y[::-1], S[::-1]


if __name__ == "__main__":
    q = float(sys.argv[1]) if len(sys.argv) > 1 else 3.0
    etas = [float(a) for a in sys.argv[2:]] or [np.sqrt(2.0)]
    print(f"=== FP-PDE ground truth, q={q:g}  (V=sign(Y)|Y|^{q:g}) ===")
    if abs(q - 2.0) < 1e-9:
        print("  [q=2 REGRESSION GATE: expect f=Var/eta^2 ~ 0.237, skew ~ +0.607, exk ~ -0.237 at beta=2]")
    for eta in etas:
        for dp, dt in ((0.02, 5e-5), (0.01, 2.5e-5)):
            y, F = solve_fp_q(q, eta=eta, dp=dp, dt=dt)
            m1, sd, sk, exk, c5, c6 = FP.cumulants_from_F(y, F)
            var = sd * sd; e2 = eta * eta
            print(f"  eta={eta:.5f} (beta={4/e2:.3f})  dp={dp} dt={dt:g}:  mean={m1:+.5f}  "
                  f"Var={var:.5f}  f=Var/eta^2={var/e2:.5f}  skew={sk:+.4f}  exk={exk:+.4f}")
