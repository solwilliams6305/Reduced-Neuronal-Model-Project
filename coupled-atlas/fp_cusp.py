"""
Smooth W_2 via the escape Fokker-Planck PDE (no Monte-Carlo noise) + sigma-form test.

Riccati density rho(p,tau) for  dp=(V-p^2)dtau+eta dW, V=sign(Y)Y^2, Y=Y0-tau:
    d_tau rho = -d_p[(V(Y)-p^2) rho] + (eta^2/2) d_pp rho,
conservative finite-volume (upwind advection + central diffusion), absorbing at p=pmin
(escape). Survival S(tau)=int rho dp gives the escape-location CDF  F(y)=S(Y0-y).
Then test whether sigma=(log F)' obeys a Painleve sigma-form (now with SMOOTH data).
numpy only.
"""
import numpy as np
from itertools import combinations

def solve_fp(eta=np.sqrt(2.0), Y0=3.0, pmin=-8.0, pmax=7.0, dp=0.02, dt=1.0e-4, tau_max=7.5, p0=None):
    D=eta*eta/2.0
    pc=np.arange(pmin+dp/2, pmax, dp); N=len(pc)
    pf=np.arange(pmin, pmax+dp/2, dp); pf2=pf*pf      # faces, size N+1
    if p0 is None: p0=np.sqrt(max(np.sign(Y0)*Y0*Y0,1e-9))
    rho=np.exp(-0.5*((pc-p0)/0.30)**2); rho/=rho.sum()*dp
    nst=int(round(tau_max/dt))
    taus=np.empty(nst+1); S=np.empty(nst+1); taus[0]=0.0; S[0]=rho.sum()*dp
    for k in range(nst):
        tau=k*dt; Y=Y0-tau; V=np.sign(Y)*Y*Y
        drift=V-pf2                                   # at faces, size N+1
        rho_left=np.concatenate(([0.0],rho))          # rho[j-1]
        rho_right=np.concatenate((rho,[0.0]))         # rho[j]
        adv=np.where(drift>0, rho_left, rho_right)*drift
        diff=-D*(rho_right-rho_left)/dp
        J=adv+diff
        J[-1]=0.0                                     # REFLECTING right boundary (no escape rightward)
        rho=rho - dt*(J[1:]-J[:-1])/dp                # left face J[0] absorbing (escape)
        np.maximum(rho,0.0,out=rho)                   # positivity
        taus[k+1]=tau+dt; S[k+1]=rho.sum()*dp
    y=Y0-taus                                         # escape location
    return y[::-1], S[::-1]                            # ascending y, F(y)=S

def cumulants_from_F(y,F):
    F=np.clip(F,0,1); f=np.gradient(F,y); f=np.clip(f,0,None); f/=np.trapz(f,y)
    m1=np.trapz(y*f,y); c=y-m1; v=np.trapz(c**2*f,y); sd=np.sqrt(v); z=c/sd
    m3=np.trapz(z**3*f,y); m4=np.trapz(z**4*f,y); m5=np.trapz(z**5*f,y); m6=np.trapz(z**6*f,y)
    return m1,sd,m3,m4-3, m5-10*m3, m6-15*m4-10*m3**2+30

def derivs(zg,logF,zev,deg=10):
    c=np.polyfit(zg,logF,deg); p=np.poly1d(c)
    return p.deriv(1)(zev),p.deriv(2)(zev),p.deriv(3)(zev)

def basis(z,s,sp):
    return np.column_stack([sp**3,sp**2,sp,np.ones_like(z),z*sp**2,z*sp,z,s*sp,s,s**2,z**2*sp**2,z*s*sp])

def sparse_cv(z,sig,sp,spp,k=5):
    y=spp**2; X=basis(z,sig,sp); sc=np.maximum(np.abs(X).max(0),1e-12); Xn=X/sc
    n=len(z); h=n//2; sst=np.sum((y-y.mean())**2); best=(-1e9,None)
    for cols in combinations(range(Xn.shape[1]),k):
        A=Xn[:,cols]; cf,*_=np.linalg.lstsq(A,y,rcond=None); R2=1-np.sum((y-A@cf)**2)/sst
        if R2>best[0]: best=(R2,cols)
    cols=best[1]; A=Xn[:,cols]; cf,*_=np.linalg.lstsq(A[:h],y[:h],rcond=None)
    cv=1-np.sum((y[h:]-A[h:]@cf)**2)/np.sum((y[h:]-y[h:].mean())**2)
    return cv,best[0],cols

def first_order_R2(z, sig, sp):
    # is sigma' determined by a 1st-order polynomial Riccati in (z, sigma)?  (Gaussian: sp=-sig^2-z sig => R2=1)
    X=np.column_stack([np.ones_like(z), z, z*z, sig, sig*sig, z*sig, z*z*sig, z*sig*sig])
    sc=np.maximum(np.abs(X).max(0),1e-12); Xn=X/sc
    cf,*_=np.linalg.lstsq(Xn, sp, rcond=None); pred=Xn@cf
    return 1-np.sum((sp-pred)**2)/np.sum((sp-sp.mean())**2)

def second_order_gain(z, sig, sp, spp):
    # does adding sigma'' (i.e. a genuine 2nd-order sigma-form) explain (sp residual)? compare 1st vs 2nd order R2 on sp
    r1=first_order_R2(z,sig,sp)
    # 2nd order: sp determined including spp?  fit sp on (1,z,z^2,sig,sig^2,z sig, spp, z spp, sig spp)
    X=np.column_stack([np.ones_like(z),z,z*z,sig,sig*sig,z*sig,spp,z*spp,sig*spp])
    sc=np.maximum(np.abs(X).max(0),1e-12); Xn=X/sc; cf,*_=np.linalg.lstsq(Xn,sp,rcond=None)
    r2=1-np.sum((sp-Xn@cf)**2)/np.sum((sp-sp.mean())**2)
    return r1, r2

NAMES=['sp^3','sp^2','sp','1','z sp^2','z sp','z','sig sp','sig','sig^2','z^2 sp^2','z sig sp']

if __name__=="__main__":
    y,F=solve_fp()
    print(f"FP solved: y in [{y[0]:.2f},{y[-1]:.2f}], F in [{F.min():.4f},{F.max():.4f}]")
    m1,sd,sk,ek,k5,k6=cumulants_from_F(y,F)
    print(f"  FP cumulants: mean={m1:+.3f} std={sd:.3f} skew={sk:+.3f} exk={ek:+.3f} k5={k5:+.2f} k6={k6:+.2f}")
    print(f"  MC target:    mean=-1.628 std=0.684 skew=+0.61  exk=-0.24  k5=-2.2 k6=-2.7")
    # sigma-form test on the SMOOTH F (standardized)
    z=(y-m1)/sd; m=(F>0.02)&(F<0.985)
    zg=z[m]; logF=np.log(F[m])
    zev=np.linspace(zg.min()+0.15, zg.max()-0.15, 60)
    sig,sp,spp=derivs(zg,logF,zev)
    cv,full,cols=sparse_cv(zev,sig,sp,spp)
    print(f"\n  sigma-form (SMOOTH FP):  CV-R^2 = {cv:+.3f}   (in-sample {full:.4f})")
    print(f"    best size-5 monomials: {[NAMES[c] for c in cols]}")
    print(f"  [calibration TW2-smooth gave CV-R^2 ~ +0.99; Gaussian-null negative]")
    np.save("fp_cusp_F.npy", np.vstack([y,F]))

    # ---- DECISIVE: 1st-order (trivial/Gaussian-like) vs genuinely 2nd-order (Painleve) ----
    print("\n=== 1st-order vs 2nd-order discrimination (smooth data) ===")
    print("  (Gaussian sigma obeys a 1st-order Riccati => R2_1st~1, trivial; Painleve needs 2nd order)")
    # cusp (FP smooth)
    r1,r2=second_order_gain(zev,sig,sp,spp)
    print(f"  cusp W_2 (FP) : R2_1st={r1:.4f}   R2_2nd={r2:.4f}   2nd-order gain={r2-r1:+.4f}")
    # TW2 smooth control
    d=np.load("hoTW_F_k1.npy"); sv,Ft=d
    ft=np.gradient(Ft,sv); ft=np.clip(ft,0,None); ft/=np.trapz(ft,sv)
    mu=np.trapz(sv*ft,sv); sd2=np.sqrt(np.trapz((sv-mu)**2*ft,sv)); zt=(sv-mu)/sd2; mt=(Ft>0.02)&(Ft<0.985)
    zte=np.linspace(zt[mt].min()+0.15, zt[mt].max()-0.15,60)
    st,spt,sppt=derivs(zt[mt],np.log(Ft[mt]),zte)
    r1t,r2t=second_order_gain(zte,st,spt,sppt)
    print(f"  TW2 (smooth)  : R2_1st={r1t:.4f}   R2_2nd={r2t:.4f}   2nd-order gain={r2t-r1t:+.4f}")
    # smooth Gaussian control (analytic density integrated)
    zg2=np.arange(-5,5,0.01); fg=np.exp(-zg2**2/2)/np.sqrt(2*np.pi); Fg=np.cumsum(fg)*0.01; Fg/=Fg[-1]
    mg=(Fg>0.02)&(Fg<0.985); zge=np.linspace(zg2[mg].min()+0.15,zg2[mg].max()-0.15,60)
    sg,spg,sppg=derivs(zg2[mg],np.log(Fg[mg]),zge)
    r1g,r2g=second_order_gain(zge,sg,spg,sppg)
    print(f"  Gaussian      : R2_1st={r1g:.4f}   R2_2nd={r2g:.4f}   2nd-order gain={r2g-r1g:+.4f}")
    print("\n  => if cusp R2_1st << 1 (like TW2) and Gaussian R2_1st ~ 1, then W_2 is GENUINELY 2nd-order (Painleve-type).")
