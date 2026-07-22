"""Skewness resummation: skew(eta) = eta * T(x)/V(x)^{3/2},  x=eta^2,
  V(x)=sum v_k x^k  (variance ladder),  T(x)=sum t_j x^j  (kappa_3 ladder /eta^4).
Median-Borel resum T and V separately, form the ratio, compare to fp_cusp ground-truth
skew (+0.601 at beta=2).  Second CENTERED observable test.
"""
import json, sys, numpy as np
from _w_resum import median_resum, naive_partial

V = np.array([0.134, 0.111, 0.104, -0.030, -0.451, -1.19, -1.90])

def load_T(nfile):
    d = json.load(open(nfile)); return np.array(d['t'])

def skew_resum(Tc, x, phi=25.0, LV=3, MV=3, LT=None, MT=None):
    nT = len(Tc)
    if LT is None:  # diagonal-ish for T (nT coeffs)
        LT = (nT-1)//2; MT = (nT-1)-LT
    Vr,_,_ = median_resum(V, x, phi_deg=phi, L=LV, M=MV)
    Tr,_,_ = median_resum(Tc, x, phi_deg=phi, L=LT, M=MT)
    return np.sqrt(x) * Tr.real / (Vr.real**1.5)

if __name__ == "__main__":
    gt = {round(r['eta2'],4): r['skew'] for r in json.load(open('_w_groundtruth.json'))}
    nfile = sys.argv[1] if len(sys.argv)>1 else '_w_kappa3_n24.json'
    Tc = load_T(nfile)
    print(f"T ladder ({nfile}): {np.round(Tc,4)}")
    print("\n x=eta^2 beta   naive-skew    resum-skew   ground-truth")
    for x in [0.5,1.0,1.5,2.0,2.25,2.56]:
        Vn = naive_partial(V,x)[-1]; Tn = naive_partial(Tc,x)[-1]
        naive = np.sqrt(x)*Tn/Vn**1.5 if Vn>0 else float('nan')
        sk = skew_resum(Tc, x)
        g = gt.get(round(x,4)); gs = f"{g:+.4f}" if g is not None else "   -   "
        tag = " <-- beta=2" if abs(x-2.0)<1e-6 else ""
        print(f" {x:5.2f} {4/x:4.2f}   {naive:+9.3f}    {sk:+.4f}     {gs}{tag}", flush=True)

    print("\nresummed skew at beta=2 across V-Pade [L/M] (T diagonal):")
    for (L,M) in [(3,3),(2,3),(3,2),(2,4)]:
        try:
            Vr,_,_ = median_resum(V, 2.0, phi_deg=25.0, L=L, M=M)
            nT=len(Tc); LT=(nT-1)//2; MT=(nT-1)-LT
            Tr,_,_ = median_resum(Tc, 2.0, phi_deg=25.0, L=LT, M=MT)
            sk = np.sqrt(2)*Tr.real/Vr.real**1.5
            print(f"  V[{L}/{M}] T[{LT}/{MT}]  skew(2) = {sk:+.4f}   (ground truth +0.601)")
        except Exception as e:
            print(f"  [{L}/{M}] failed: {e}")
