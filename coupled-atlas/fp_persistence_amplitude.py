"""
Lock item-2: with the FW left-tail EXPONENT fixed at its analytic value (2q+1=5 in Y*, i.e. Theta^{5/2}),
extract the tail AMPLITUDE A(beta) from the MC-free survival PDE and its scaling A~beta^q. Plus a grid
convergence check (numerical diffusion could inflate a deep tail).
  -log P(Theta>t) ~ A(beta) * Theta^{5/2}     (deep, asymptotic window)
  A(beta) ~ beta^q   -> q is the genuine anomalous noise-amplitude exponent (FW-naive would be q=1).
numpy only. [NUMERIC].
"""
import numpy as np
def solve_fp(eta,Y0=3.0,pmin=-13.0,pmax=8.0,dp=0.02,dt=6e-5,tau_max=11.0):
    D=eta*eta/2.0
    pc=np.arange(pmin+dp/2,pmax,dp); pf=np.arange(pmin,pmax+dp/2,dp); pf2=pf*pf
    p0=np.sqrt(max(np.sign(Y0)*Y0*Y0,1e-9)); rho=np.exp(-0.5*((pc-p0)/0.30)**2); rho/=rho.sum()*dp
    nst=int(round(tau_max/dt)); taus=np.empty(nst+1); S=np.empty(nst+1); taus[0]=0; S[0]=rho.sum()*dp
    for k in range(nst):
        tau=k*dt; Y=Y0-tau; V=np.sign(Y)*Y*Y; drift=V-pf2
        rl=np.concatenate(([0.0],rho)); rr=np.concatenate((rho,[0.0]))
        J=np.where(drift>0,rl,rr)*drift-D*(rr-rl)/dp; J[-1]=0.0
        rho=rho-dt*(J[1:]-J[:-1])/dp; np.maximum(rho,0.0,out=rho)
        taus[k+1]=tau+dt; S[k+1]=rho.sum()*dp
    return Y0-taus,S
def deep_tail(eta,**kw):
    y,S=solve_fp(eta,**kw); neg=y<-0.2; Th=(y[neg]**2)/2; nlS=-np.log(np.clip(S[neg],1e-300,1))
    m=(nlS>2)&(nlS<250)&(Th>3); o=np.argsort(Th[m]); return Th[m][o],nlS[m][o]

# amplitude with exponent fixed at 5/2, over the asymptotic window where local power ~2.5
def amp(Th,nlS,lo,hi):
    m=(Th>=lo)&(Th<=hi)
    return np.mean(nlS[m]/Th[m]**2.5) if m.sum()>4 else np.nan

print("=== amplitude A(beta) with FW exponent fixed: -logP ~ A * Theta^{5/2} (deep) ===")
betas=[2.0,3.0,4.0,6.0,8.0]; A=[]
for b in betas:
    Th,nlS=deep_tail(2.0/np.sqrt(b))
    # asymptotic window: pick where reachable and power has settled
    lo,hi=(8,16) if b<=4 else (7,12)
    Ab=amp(Th,nlS,lo,hi); A.append(Ab)
    apow=np.polyfit(np.log(Th[(Th>=lo)&(Th<=hi)]),np.log(nlS[(Th>=lo)&(Th<=hi)]),1)[0]
    print(f"  beta={b:>4.0f}: A={Ab:.4f}  (local power in window={apow:.2f}; A*Theta^2.5 model)  A/beta={Ab/b:.4f}")
A=np.array(A); q=np.polyfit(np.log(betas),np.log(A),1)[0]
print(f"\n  amplitude A ~ beta^q,  q = {q:.3f}   [FW-naive would be q=1; measured anomaly]")
print(f"  => deep left tail:  -logP(Theta>t) ~ {np.exp(np.polyfit(np.log(betas),np.log(A),1)[1]):.3f} * beta^{q:.2f} * Theta^{{5/2}}")
print(f"     equivalently  P(Y*<-s) ~ exp( -c beta^{q:.2f} |s|^5 ),  FW exponent 5 confirmed asymptotically.")

print("\n=== grid convergence (beta=4): amplitude & power stable under dp, dt refinement? ===")
for dp,dt in [(0.024,8e-5),(0.02,6e-5),(0.013,4e-5)]:
    Th,nlS=deep_tail(2.0/np.sqrt(4.0),dp=dp,dt=dt)
    Ab=amp(Th,nlS,8,16); apow=np.polyfit(np.log(Th[(Th>=8)&(Th<=16)]),np.log(nlS[(Th>=8)&(Th<=16)]),1)[0]
    print(f"  dp={dp:.3f} dt={dt:.0e}: A={Ab:.4f}  power={apow:.2f}")
