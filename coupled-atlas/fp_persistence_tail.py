"""
Frontier item 2, MC-free: extract the deep left-tail (persistence) law of W_beta from the survival PDE
(fp_cusp), which reaches Theta = 15-30 where MC has zero samples. Question: does a clean asymptotic tail
exponent exist there, or is the left tail a genuine crossover with no single universal exponent?

fp_cusp: rho(p,tau) for dp=(V-p^2)dtau+eta dW, absorbing at pmin (escape). S(tau)=int rho dp = survival =
P(escape location Y* < Y0 - tau).  Deep left tail: -log S(tau) vs Theta=(Y0-tau)^2/2  (Y0-tau<0).
Sweep beta=4/eta^2.  Report local power a(Theta)=d log(-logS)/d log Theta deep in the asymptotic region.
numpy only. [NUMERIC].
"""
import numpy as np

def solve_fp(eta, Y0=3.0, pmin=-12.0, pmax=8.0, dp=0.02, dt=6e-5, tau_max=11.0):
    D=eta*eta/2.0
    pc=np.arange(pmin+dp/2,pmax,dp); pf=np.arange(pmin,pmax+dp/2,dp); pf2=pf*pf
    p0=np.sqrt(max(np.sign(Y0)*Y0*Y0,1e-9)); rho=np.exp(-0.5*((pc-p0)/0.30)**2); rho/=rho.sum()*dp
    nst=int(round(tau_max/dt)); taus=np.empty(nst+1); S=np.empty(nst+1); taus[0]=0; S[0]=rho.sum()*dp
    for k in range(nst):
        tau=k*dt; Y=Y0-tau; V=np.sign(Y)*Y*Y; drift=V-pf2
        rl=np.concatenate(([0.0],rho)); rr=np.concatenate((rho,[0.0]))
        J=np.where(drift>0,rl,rr)*drift - D*(rr-rl)/dp; J[-1]=0.0
        rho=rho-dt*(J[1:]-J[:-1])/dp; np.maximum(rho,0.0,out=rho)
        taus[k+1]=tau+dt; S[k+1]=rho.sum()*dp
    return Y0-taus, S    # escape location y (descending), survival S=P(Y*<y)

print("=== MC-free deep left tail from survival PDE (fp_cusp), per beta ===")
print("  local power a(Theta)=d log(-logS)/d log Theta in the DEEP region (Theta=8..25, MC-inaccessible)")
res={}
for beta in [2.0,4.0,8.0,16.0]:
    y,S=solve_fp(2.0/np.sqrt(beta))
    neg=y<-0.2; yy=y[neg]; SS=np.clip(S[neg],1e-300,1)
    Th=yy*yy/2.0; nlS=-np.log(SS)
    # keep the numerically reliable deep window (S below the bulk, above the floor)
    m=(SS<0.2)&(SS>1e-120)&(Th>3)
    Th=Th[m]; nlS=nlS[m]
    order=np.argsort(Th); Th=Th[order]; nlS=nlS[order]
    res[beta]=(Th,nlS)
    # local power in two depth bands
    def locpow(lo,hi):
        mm=(Th>=lo)&(Th<=hi)
        return np.polyfit(np.log(Th[mm]),np.log(nlS[mm]),1)[0] if mm.sum()>5 else np.nan
    a_mid=locpow(4,9); a_deep=locpow(10,20); a_vdeep=locpow(20,40)
    # asymptotic Kramers/action check: is -logS ~ Theta (linear/exponential) or Theta^{3/2} or Theta^2 deep?
    print(f"  beta={beta:>4.0f}: a[4-9]={a_mid:.2f}  a[10-20]={a_deep:.2f}  a[20-40]={a_vdeep:.2f}  "
          f"(deepest Theta={Th.max():.0f}, -logS up to {nlS.max():.0f})")

# Cross-beta: at fixed DEEP Theta, how does -logS scale with beta?  (the true 'persistence exponent' content)
print("\n=== beta-scaling of -logS at fixed deep Theta (MC-free) ===")
for Th0 in [8.0,12.0,18.0]:
    vals=[]; bs=[]
    for beta in [2.0,4.0,8.0,16.0]:
        Th,nlS=res[beta]
        if Th.max()>=Th0:
            vals.append(np.interp(Th0,Th,nlS)); bs.append(beta)
    if len(vals)>=3:
        q=np.polyfit(np.log(bs),np.log(vals),1)[0]
        print(f"  Theta={Th0:>4.0f}: -logS = {[f'{v:.1f}' for v in vals]} at beta={bs}  =>  ~beta^{q:.2f}")
np.savez("fp_persistence_data.npz",**{f"Th_{int(b)}":res[b][0] for b in res},**{f"nlS_{int(b)}":res[b][1] for b in res})
print("\nsaved fp_persistence_data.npz")
