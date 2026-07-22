"""
Tier B empirical rate probe: how fast does the finite-coalescence escape law W_{beta,Delta} converge to the
pure-cusp W_beta as Delta -> 0?   V_Delta(Y) = sign(Y)|Y|(|Y|+Delta)   (COUPLED_CUSP_RESULTS sec 2.1);
Delta=0 is the pure cusp.  Delta ~ eps^{1/4} near the cusp, so a rate |cum(Delta)-cum(0)| ~ Delta^r gives the
outer-correction rate in eps as p = r/4.  This addresses the deep-research OPEN QUESTION (deterministic ceiling
on p), at least empirically.  beta=2 fixed (eta=sqrt2) so FP is reliable (NOT the small-eta boundary-layer regime).
Adapts fp_cusp.solve_fp with a Delta-dependent potential.
"""
import numpy as np

def solve_fp_delta(Delta, eta=np.sqrt(2.0), Y0=3.0, pmin=-8.0, pmax=7.0, dp=0.02, dt=1.0e-4, tau_max=7.5):
    D=eta*eta/2.0
    pc=np.arange(pmin+dp/2, pmax, dp); N=len(pc)
    pf=np.arange(pmin, pmax+dp/2, dp); pf2=pf*pf
    def Vd(Y):
        a=np.abs(Y); return np.sign(Y)*a*(a+Delta)
    p0=np.sqrt(max(Vd(Y0),1e-9))
    rho=np.exp(-0.5*((pc-p0)/0.30)**2); rho/=rho.sum()*dp
    nst=int(round(tau_max/dt))
    for k in range(nst):
        Y=Y0-k*dt; V=Vd(Y)
        drift=V-pf2
        rl=np.concatenate(([0.0],rho)); rr=np.concatenate((rho,[0.0]))
        adv=np.where(drift>0, rl, rr)*drift
        diff=-D*(rr-rl)/dp
        J=adv+diff; J[-1]=0.0
        rho=rho-dt*(J[1:]-J[:-1])/dp
        np.maximum(rho,0.0,out=rho)
    S=rho.sum()*dp
    # reconstruct CDF over escape location
    taus=np.arange(nst+1)*dt
    return None

def cumulants_delta(Delta, eta=np.sqrt(2.0), **kw):
    # replicate fp_cusp.solve_fp but with V_Delta; return (mean,std,skew,exk)
    D=eta*eta/2.0
    Y0=kw.get('Y0',3.0); pmin=kw.get('pmin',-8.0); pmax=kw.get('pmax',7.0)
    dp=kw.get('dp',0.02); dt=kw.get('dt',1.0e-4); tau_max=kw.get('tau_max',7.5)
    pc=np.arange(pmin+dp/2,pmax,dp); pf=np.arange(pmin,pmax+dp/2,dp); pf2=pf*pf
    def Vd(Y):
        a=abs(Y); return np.sign(Y)*a*(a+Delta)
    p0=np.sqrt(max(Vd(Y0),1e-9))
    rho=np.exp(-0.5*((pc-p0)/0.30)**2); rho/=rho.sum()*dp
    nst=int(round(tau_max/dt)); taus=np.empty(nst+1); S=np.empty(nst+1); taus[0]=0; S[0]=rho.sum()*dp
    for k in range(nst):
        Y=Y0-k*dt; V=Vd(Y); drift=V-pf2
        rl=np.concatenate(([0.0],rho)); rr=np.concatenate((rho,[0.0]))
        J=np.where(drift>0,rl,rr)*drift - D*(rr-rl)/dp; J[-1]=0.0
        rho=rho-dt*(J[1:]-J[:-1])/dp; np.maximum(rho,0.0,out=rho)
        taus[k+1]=(k+1)*dt; S[k+1]=rho.sum()*dp
    y=(Y0-taus)[::-1]; F=S[::-1]
    F=np.clip(F,0,1); f=np.gradient(F,y); f=np.clip(f,0,None); f/=np.trapz(f,y)
    m1=np.trapz(y*f,y); c=y-m1; v=np.trapz(c**2*f,y); sd=np.sqrt(v); z=c/sd
    sk=np.trapz(z**3*f,y); ek=np.trapz(z**4*f,y)-3
    return m1,sd,sk,ek

if __name__=="__main__":
    Deltas=[0.0,0.1,0.2,0.35,0.5,0.8,1.2]
    print("=== W_{beta=2,Delta} cumulants vs Delta (V_Delta=sign(Y)|Y|(|Y|+Delta)); Delta=0 is pure cusp ===")
    print(f"  {'Delta':>6} {'mean':>9} {'std':>8} {'skew':>8} {'exk':>8}")
    rows=[]
    for D in Deltas:
        m,s,sk,ek=cumulants_delta(D); rows.append((D,m,s,sk,ek))
        print(f"  {D:6.2f} {m:9.4f} {s:8.4f} {sk:8.4f} {ek:8.4f}")
    R=np.array(rows)
    D0=R[0]  # Delta=0 reference
    print(f"\n  reference (pure cusp, Delta=0): skew={D0[3]:.4f} exk={D0[4]:.4f}  (fp_cusp target +0.601/-0.244)")
    # fit |cum(Delta)-cum(0)| ~ Delta^r  for skew and exk, small-Delta
    m=R[:,0]>0
    for j,name in [(3,'skew'),(4,'exk'),(1,'mean')]:
        d=np.abs(R[m,j]-D0[j]); Dv=R[m,0]
        good=d>1e-4
        if good.sum()>=3:
            r=np.polyfit(np.log(Dv[good]),np.log(d[good]),1)[0]
            print(f"  |{name}(Delta)-{name}(0)| ~ Delta^{r:.2f}   => eps-rate p={r/4:.2f} (Delta~eps^1/4)")
