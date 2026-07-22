"""
Program 2, Route 2b — the NEXT variance coefficient v1 by the CONVERGED direct-functional method (no eta-ladder,
no truncation).  Var = eta^2 v0 + eta^4 v1 + ...,  with
    v0 = <Y1^2>                         (=0.1339 closed-form)
    v1 = Var(Y2) + 2 <Y1 Y3>
and mean coeff  m1 = <Y2>  (<Y*> = Y*_0 + m1 eta^2 + ...).
Node-shift functionals (derived; using u0(Y*_0)=0 => u0''(Y*_0)=0, u1''(Y*_0)=V* u1(Y*_0)):
    Y1 = -u1/u0'
    Y2 = -(u1' Y1 + u2)/u0'
    Y3 = -( (1/6)V* u0' Y1^3 + u1' Y2 + (1/2)V* u1 Y1^2 + u2' Y1 + u3 )/u0'
all fields at Y*_0.  This tests whether the weak-noise series is asymptotic/resurgent with CONVERGED (not FP) data.
numpy only, vectorized over realizations.
"""
import numpy as np

def V(Y): return np.sign(Y)*Y*Y

def run(N=300000, Y0=4.0, Yend=-5.0, h=1.0e-3, seed=0):
    rng=np.random.default_rng(seed)
    nY=int(round((Y0-Yend)/h)); Yg=Y0-np.arange(nY+1)*h; dY=-h; sh=np.sqrt(h)
    u0=np.empty(nY+1); u0p=np.empty(nY+1); u0[0]=1.0; u0p[0]=-Y0
    for i in range(nY):
        u0[i+1]=u0[i]+u0p[i]*dY; u0p[i+1]=u0p[i]+V(Yg[i])*u0[i]*dY
    istar=next(i for i in range(1,nY+1) if Yg[i]<0 and u0[i-1]*u0[i]<0)
    f=u0[istar-1]/(u0[istar-1]-u0[istar]); Ystar=Yg[istar-1]+f*(Yg[istar]-Yg[istar-1])
    u0p_star=u0p[istar-1]+f*(u0p[istar]-u0p[istar-1]); Vstar=V(Ystar)
    # perturbation fields u1,u2,u3 and derivatives
    u1=np.zeros(N);u1p=np.zeros(N);u2=np.zeros(N);u2p=np.zeros(N);u3=np.zeros(N);u3p=np.zeros(N)
    snap={}
    for i in range(nY):
        Vi=V(Yg[i]); dB=sh*rng.standard_normal(N)
        n1=u1p+Vi*u1*dY+u0[i]*dB
        n2=u2p+Vi*u2*dY+u1*dB
        n3=u3p+Vi*u3*dY+u2*dB
        u1=u1+u1p*dY; u2=u2+u2p*dY; u3=u3+u3p*dY
        u1p,u2p,u3p=n1,n2,n3
        if i in (istar-1,istar):
            snap[i]=(u1.copy(),u1p.copy(),u2.copy(),u2p.copy(),u3.copy())
    def interp(a,b): return a+f*(b-a)
    lo=snap[istar-1]; hi=snap[istar]
    u1s=interp(lo[0],hi[0]); u1ps=interp(lo[1],hi[1])
    u2s=interp(lo[2],hi[2]); u2ps=interp(lo[3],hi[3]); u3s=interp(lo[4],hi[4])
    Y1=-u1s/u0p_star
    Y2=-(u1ps*Y1+u2s)/u0p_star
    Y3=-((1/6)*Vstar*u0p_star*Y1**3 + u1ps*Y2 + 0.5*Vstar*u1s*Y1**2 + u2ps*Y1 + u3s)/u0p_star
    v0=np.mean(Y1**2)
    m1=np.mean(Y2)
    varY2=np.var(Y2)
    Y1Y3=np.mean(Y1*Y3)
    v1=varY2+2*Y1Y3
    s0=3*(np.mean(Y1**2*Y2)-np.mean(Y1**2)*np.mean(Y2))/v0**1.5
    return dict(Ystar=Ystar,v0=v0,m1=m1,varY2=varY2,Y1Y3=Y1Y3,v1=v1,s0=s0)

if __name__=="__main__":
    print("=== converged direct-functional weak-noise coefficients (Y1,Y2,Y3) ===")
    for h in [1.5e-3,1e-3]:
        agg=[run(N=300000,h=h,seed=s) for s in (1,2,3)]
        def mean(k): return np.mean([a[k] for a in agg]);
        def std(k): return np.std([a[k] for a in agg])
        print(f"\n h={h:.1e}, N=300k x3 seeds:  Y*_0={agg[0]['Ystar']:.4f}")
        print(f"   v0=<Y1^2>       = {mean('v0'):.4f} +- {std('v0'):.4f}   (closed-form 0.1339)")
        print(f"   m1=<Y2>         = {mean('m1'):+.4f} +- {std('m1'):.4f}   (FP/notes ~0.212)")
        print(f"   Var(Y2)         = {mean('varY2'):.4f} +- {std('varY2'):.4f}")
        print(f"   <Y1 Y3>         = {mean('Y1Y3'):+.4f} +- {std('Y1Y3'):.4f}")
        print(f"   v1=Var(Y2)+2<Y1Y3> = {mean('v1'):+.4f} +- {std('v1'):.4f}")
        print(f"   s0=skew/eta     = {mean('s0'):.4f} +- {std('s0'):.4f}   (MC ~1.18)")
        print(f"   => Var = eta^2( {mean('v0'):.3f} {mean('v1'):+.3f} eta^2 + ...);  |v1/v0|={abs(mean('v1')/mean('v0')):.2f}")
    print("\n [contaminated FP gave v0=0.200, v1<0 with |v1/v0|=1.23]")
