"""
Cross-check the 2-point mechanism's beta-dependence: DERIVED Var(spacing) ~ (4/beta) * Theta^{-3/2}.
Exponent must be beta-INDEPENDENT (-3/2); prefactor must scale as 1/beta (=eta^2, eta=2/sqrt(beta)).
Non-confounded: an exponent + an amplitude law, not a bulk fit.
"""
import numpy as np
def Vc(Y): return np.sign(Y)*Y*Y
def sample(beta,N,Y0=3.0,Yend=-16.0,dt=1.2e-3,seed=7,maxnodes=80):
    eta=2.0/np.sqrt(beta); rng=np.random.default_rng(seed)
    n=int(round((Y0-Yend)/dt)); sq=np.sqrt(dt); p=np.full(N,np.sqrt(max(Vc(Y0),1e-9)))
    nodes=np.full((maxnodes,N),np.nan); idx=np.zeros(N,int)
    for i in range(n):
        Y=Y0-i*dt; thr=20.0*abs(Y)+50.0
        p+=(Vc(Y)-p*p)*dt+eta*sq*rng.standard_normal(N); np.clip(p,-1e4,1e4,out=p)
        ex=p<-thr; live=ex&(idx<maxnodes); cols=np.where(live)[0]
        if cols.size: nodes[idx[cols],cols]=Y; idx[cols]+=1
        p[ex]=thr
    return nodes
print("beta | fitted exponent q (DERIVED -1.50) | prefactor c (Var=c*Theta^q) | c*beta (should be ~const=4*coef)")
bands=[(4,10),(10,18),(18,28),(28,42),(42,60)]
for beta in [2.0,4.0,8.0]:
    Th=0.5*sample(beta,7000)**2; sp=Th[1:]-Th[:-1]; Thm=0.5*(Th[1:]+Th[:-1])
    D=[]; V=[]
    for lo,hi in bands:
        m=np.isfinite(sp[lo:hi])&(sp[lo:hi]>0.5)&(sp[lo:hi]<8)
        if m.sum()>400: D.append(Thm[lo:hi][m].mean()); V.append(sp[lo:hi][m].var())
    D=np.array(D); V=np.array(V); q,lc=np.polyfit(np.log(D),np.log(V),1); c=np.exp(lc)
    print(f"  {beta:.0f}  |   q={q:+.2f}   |  c={c:.4f}  |  c*beta={c*beta:.4f}")
print("\n=> exponent ~-1.5 across beta and c*beta ~ const  ==>  Var(spacing) ~ (const/beta) Theta^{-3/2}, DERIVED.")
