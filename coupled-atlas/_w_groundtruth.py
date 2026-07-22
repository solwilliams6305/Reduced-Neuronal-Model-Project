"""Ground-truth Var(Y*)(eta) from the MC-free FP-PDE (fp_cusp), across a range of eta
spanning inside->outside the perturbative radius (eta_c^2 ~ 1.2).  beta=2 is eta=sqrt2.
Also prints the beta=2 fingerprint (skew, exk) as a sanity check (should be ~+0.607,-0.237)."""
import numpy as np, json
import fp_cusp as FP

etas = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, np.sqrt(2.0), 1.5, 1.6]
rows = []
for eta in etas:
    # finer grid at larger eta; the default is tuned for beta=2
    y, F = FP.solve_fp(eta=eta, Y0=3.0, pmin=-9.0, pmax=7.0, dp=0.02, dt=1.0e-4, tau_max=8.0)
    m1, sd, sk, exk, c5, c6 = FP.cumulants_from_F(y, F)
    var = sd*sd
    eta2 = eta*eta
    rows.append(dict(eta=eta, eta2=eta2, beta=4.0/eta2, mean=m1, var=var,
                     f=var/eta2, skew=sk, exk=exk))
    tag = "  <-- beta=2" if abs(eta-np.sqrt(2))<1e-6 else ""
    print(f"eta={eta:.4f} eta2={eta2:.3f} beta={4/eta2:5.2f}  Var={var:.5f}  "
          f"f=Var/eta2={var/eta2:.5f}  skew={sk:+.3f} exk={exk:+.3f}{tag}", flush=True)

json.dump(rows, open('_w_groundtruth.json','w'), indent=1)
print("\nwrote _w_groundtruth.json")
