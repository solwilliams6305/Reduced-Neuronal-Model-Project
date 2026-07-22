"""
Program 2, Route 2b / Step 0 — pin the leading weak-noise coefficients by a CONVERGED perturbative method
off the EXACT deterministic Weber backbone (no FP boundary layer, no rare-escape MC bias).

Weak-noise expansion of the escape location:  Y* = Y*_0 + eta Y1 + eta^2 Y2 + ...
with Y1 (linear in noise, Gaussian) and Y2 (quadratic) FUNCTIONALS of the noise, eta-independent.
Then:   Var(Y*) = C_V eta^2 + O(eta^4),  C_V = <Y1^2>   (should reproduce the closed-form 0.1339)
        skew    = s0 eta   + O(eta^3),  s0 = 3 <Y1^2 Y2>_c / C_V^{3/2}
(<.>_c = connected;  <Y1^3>=0, <Y1 Y2>=0 by W->-W symmetry.)

Field perturbation (noise = random potential, u'' = (V + xi) u,  xi white in depth):
   u1'' - V u1 = xi u0 ,   u2'' - V u2 = xi u1 ,  recessive BC (u1=u1'=0 at Y0).
Node-shift (Y-derivatives, u0(Y*_0)=0):
   Y1 = -u1(Y*_0)/u0'(Y*_0)
   Y2 = -( 0.5 u0''(Y*_0) Y1^2 + u1'(Y*_0) Y1 + u2(Y*_0) ) / u0'(Y*_0)
This is a "perturbative MC": sample xi, integrate the LINEAR perturbation ODEs (fast, converged), average.
Decisively resolves the FP-vs-MC small-eta discrepancy (FP: C_V~0.20, s0~0.66;  MC: C_V~0.134, s0~1.18).
numpy only. Vectorized over realizations.
"""
import numpy as np

def V(Y): return np.sign(Y)*Y*Y

def run(N=200000, Y0=4.0, Yend=-5.0, h=1.0e-3, seed=0):
    rng=np.random.default_rng(seed)
    nY=int(round((Y0-Yend)/h))
    Yg=Y0-np.arange(nY+1)*h
    # deterministic recessive backbone u0 (u0(Y0)=1, u0'(Y0)=-Y0), Y-march
    u0=np.empty(nY+1); u0p=np.empty(nY+1); u0[0]=1.0; u0p[0]=-Y0
    dY=-h
    for i in range(nY):
        Y=Yg[i]
        u0[i+1]=u0[i]+u0p[i]*dY
        u0p[i+1]=u0p[i]+V(Y)*u0[i]*dY
    # first node Y*_0 (Y<0)
    istar=None
    for i in range(1,nY+1):
        if Yg[i]<0 and u0[i-1]*u0[i]<0: istar=i; break
    f=u0[istar-1]/(u0[istar-1]-u0[istar]); Ystar=Yg[istar-1]+f*(Yg[istar]-Yg[istar-1])
    u0p_star=u0p[istar-1]+f*(u0p[istar]-u0p[istar-1])
    u0_star=0.0
    u0pp_star=V(Ystar)*u0_star  # =0 (since u0(Y*)=0) -> the 0.5 u0'' Y1^2 term vanishes!
    # perturbation SDEs, vectorized over realizations, marched in Y
    u1=np.zeros(N); u1p=np.zeros(N); u2=np.zeros(N); u2p=np.zeros(N)
    sh=np.sqrt(h)
    # store values at the two grid points bracketing Ystar to interpolate
    for i in range(nY):
        Y=Yg[i]; Vi=V(Y)
        dB=sh*rng.standard_normal(N)
        # update derivatives (use current u1,u2), then positions (use current u1p,u2p)
        new_u1p=u1p+Vi*u1*dY+u0[i]*dB
        new_u2p=u2p+Vi*u2*dY+u1*dB
        u1=u1+u1p*dY; u2=u2+u2p*dY
        u1p=new_u1p; u2p=new_u2p
        if i==istar-1:
            u1_lo,u1p_lo,u2_lo=u1.copy(),u1p.copy(),u2.copy()
        if i==istar:
            u1_hi,u1p_hi,u2_hi=u1.copy(),u1p.copy(),u2.copy()
    # interpolate perturbation fields to Ystar
    u1s=u1_lo+f*(u1_hi-u1_lo); u1ps=u1p_lo+f*(u1p_hi-u1p_lo); u2s=u2_lo+f*(u2_hi-u2_lo)
    Y1=-u1s/u0p_star
    Y2=-(0.5*u0pp_star*Y1**2 + u1ps*Y1 + u2s)/u0p_star
    C_V=np.mean(Y1**2)
    Y1sq=Y1**2
    conn=np.mean(Y1sq*Y2)-np.mean(Y1sq)*np.mean(Y2)
    s0=3*conn/C_V**1.5
    return Ystar,u0p_star,C_V,s0,Y1,Y2

if __name__=="__main__":
    print("=== weak-noise perturbative sectors off the exact Weber backbone ===")
    print("  (converged: no FP boundary layer, no rare-escape MC bias)")
    for h in [2e-3,1e-3,5e-4]:
        Ystar,u0ps,C_V,s0,Y1,Y2=run(N=150000,h=h,seed=1)
        print(f"  h={h:.0e}: Y*_0={Ystar:.4f} u0'(Y*)={u0ps:.3f}  C_V=<Y1^2>={C_V:.4f}  s0(skew/eta)={s0:.3f}")
    print("\n  closed-form / MC targets:  C_V=0.1339 (closed-form, MC-confirmed);  s0(skew/eta) MC ~1.18")
    print("  FP small-eta (contaminated): C_V~0.20, s0~0.66")
    # seed robustness at finest h
    print("\n  seed robustness (h=5e-4, N=150k):")
    for sd in [1,2,3]:
        _,_,C_V,s0,_,_=run(N=150000,h=5e-4,seed=sd)
        print(f"    seed={sd}: C_V={C_V:.4f}  s0={s0:.3f}")
