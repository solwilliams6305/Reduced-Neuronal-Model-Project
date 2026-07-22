"""Median Borel-Pade resummation of the MEAN series
      <Y*>(x) = m_0 + sum_{j>=1} m_j x^j,   x=eta^2,  m_0=-2.188,
using the engine mean ladder m_1..m_6, compared to fp_cusp ground-truth mean, esp beta=2.
A SECOND observable test of the resurgent representation (variance was the first)."""
import json, numpy as np
from _w_resum import median_resum, naive_partial

M0 = -2.188
# engine mean ladder (n=24); m_6 has residual grid drift -> also test extrapolated value
LADDER_N24 = [0.21348, 0.08334, -0.00542, -0.10560, -0.19708, 0.18958]
M6_EXTRAP  = 0.17
def coeffs(m6):
    return np.array([0.0] + LADDER_N24[:5] + [m6])   # a_0=0, a_1..a_6 = m_1..m_6

gt = {round(r['eta2'],4): r['mean'] for r in json.load(open('_w_groundtruth.json'))}

print("MEAN resummation:  <Y*> = -2.188 + median-Borel-resum( sum m_j x^j )")
print(" x=eta^2 beta   naive        resum(m6=.19)  resum(m6=.17)  ground-truth")
for x in [0.5,1.0,1.5,2.0,2.25,2.56]:
    naive = M0 + naive_partial(coeffs(0.18958), x)[-1]
    r1,_,_ = median_resum(coeffs(0.18958), x, phi_deg=25.0); r1 = M0 + r1.real
    r2,_,_ = median_resum(coeffs(M6_EXTRAP), x, phi_deg=25.0); r2 = M0 + r2.real
    g = gt.get(round(x,4)); gs = f"{g:+.4f}" if g is not None else "   -   "
    tag = " <-- beta=2" if abs(x-2.0)<1e-6 else ""
    print(f" {x:5.2f} {4/x:4.2f}  {naive:+10.3f}   {r1:+.4f}      {r2:+.4f}     {gs}{tag}", flush=True)

# Pade-order robustness at beta=2
print("\nresummed <Y*> at beta=2 across Pade [L/M] (m6=0.18958):")
c = coeffs(0.18958)
for (L,Mq) in [(3,3),(2,3),(3,2),(2,4)]:
    try:
        r,_,_ = median_resum(c, 2.0, phi_deg=25.0, L=L, M=Mq)
        print(f"  [{L}/{Mq}]  <Y*>(2) = {M0+r.real:+.4f}   (ground truth -1.6072)")
    except Exception as e:
        print(f"  [{L}/{Mq}] failed: {e}")
