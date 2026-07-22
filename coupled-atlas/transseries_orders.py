"""
Program 2, Route 2b (resurgent trans-series) — extract the WEAK-NOISE PERTURBATIVE SECTORS of W cleanly via the
MC-free FP-PDE (fp_cusp), at a ladder of small eta. Fit the eta-expansions of mean/var/skew to several orders and
inspect the coefficient structure (growth = asymptotic/resurgent signature; Borel singularity at instanton S=s^5/10).
"""
import numpy as np
from fp_cusp import solve_fp, cumulants_from_F
etas=np.sqrt(np.array([0.05,0.08,0.11,0.15,0.20,0.26,0.33,0.42]))
rows=[]
for eta in etas:
    y,F=solve_fp(eta=eta, dp=0.015, dt=6e-5)
    m1,sd,sk,ek,k5,k6=cumulants_from_F(y,F)
    rows.append((eta**2, m1, sd*sd, sk, ek)); 
    print(f"  eta^2={eta**2:.3f} (beta={4/eta**2:5.1f}): mean={m1:+.4f} var={sd*sd:.5f} skew={sk:+.4f} exk={ek:+.4f}")
R=np.array(rows); e2=R[:,0]
# mean = Y0 + m1 e2 + m2 e2^2 + m3 e2^3
cm=np.polyfit(e2,R[:,1],3)[::-1]
# var/e2 = v0 + v1 e2 + v2 e2^2
cv=np.polyfit(e2,R[:,2]/e2,2)[::-1]
# skew/eta = s0 + s1 e2 + s2 e2^2   (skew ~ eta = sqrt(e2))
cs=np.polyfit(e2,R[:,3]/np.sqrt(e2),2)[::-1]
print("\n  perturbative trans-series (eta^2 = 4/beta):")
print(f"    <Y*> = {cm[0]:+.3f} {cm[1]:+.3f} e2 {cm[2]:+.3f} e2^2 {cm[3]:+.3f} e2^3")
print(f"    Var  = e2*( {cv[0]:.3f} {cv[1]:+.3f} e2 {cv[2]:+.3f} e2^2 )")
print(f"    skew = eta*( {cs[0]:+.3f} {cs[1]:+.3f} e2 {cs[2]:+.3f} e2^2 )")
print(f"\n  coefficient ratios (growth => asymptotic/resurgent, Borel-sing at instanton S):")
print(f"    mean |m2/m1|={abs(cm[2]/cm[1]):.2f}, |m3/m2|={abs(cm[3]/cm[2]):.2f}")
print(f"    var  |v1/v0|={abs(cv[1]/cv[0]):.2f}, |v2/v1|={abs(cv[2]/cv[1]):.2f}")
