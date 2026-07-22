"""
Cross-check Piece-1 GO via the ADDITIVE-noise Riccati (no multiplicative-noise-vanishes-at-zero artifact):
dp=(V-p^2)dtau + eta dW; each explosion (p<-thr, reset p=+thr) = a node. Count nodes vs phase Theta=Y^2/2.
Var(N) sub-linear/bounded => RIGID (confirms GO). Also phase-spacing CV between consecutive nodes.
Compares to 1-D random-Schrodinger spectral rigidity (Sturm count). numpy only.
"""
import numpy as np
def Vc(Y): return np.sign(Y)*Y*Y
def run(eta,N,Y0=3.0,Yend=-9.0,dt=2e-3,seed=11,thr=25.0):
    rng=np.random.default_rng(seed); n=int(round((Y0-Yend)/dt)); sq=np.sqrt(dt)
    p=np.full(N,np.sqrt(Vc(Y0))); cnt=np.zeros(N)
    cps=np.arange(-1.0,-9.0,-1.0); ci=0; Ncp=np.zeros((len(cps),N))
    n1=np.full(N,np.nan); n2=np.full(N,np.nan); n3=np.full(N,np.nan)
    for i in range(n):
        Y=Y0-i*dt
        p=p+(Vc(Y)-p*p)*dt+eta*sq*rng.standard_normal(N); np.clip(p,-1e3,1e3,out=p)
        ex=p<-thr
        # record first 3 node depths for spacing
        a=ex&np.isnan(n1); n1[a]=Y
        b=ex&(~np.isnan(n1))&np.isnan(n2)&(cnt>=1); n2[b]=Y
        c=ex&(~np.isnan(n2))&np.isnan(n3)&(cnt>=2); n3[c]=Y
        cnt=cnt+ex; p[ex]=thr
        if ci<len(cps) and Y<=cps[ci]: Ncp[ci]=cnt.copy(); ci+=1
    return cps*cps/2.0, Ncp, n1,n2,n3
Th,Ncp,n1,n2,n3=run(np.sqrt(2.0),6000)
Nbar=Ncp.mean(axis=1); Nvar=Ncp.var(axis=1)
print("Theta | mean N | Var(N) | Var/mean (Poisson=1, rigid<<1)")
for t,m,vv in zip(Th,Nbar,Nvar):
    if m>0.5: print(f"  {t:6.2f} | {m:5.2f} | {vv:6.3f} | {vv/m:5.2f}")
ok=Nbar>1.0; alpha=np.polyfit(np.log(Th[ok]),np.log(Nvar[ok]+1e-9),1)[0]
g=~np.isnan(n2); sp=(n1[g]**2-n2[g]**2)/2.0  # phase gap node1->node2
g2=~np.isnan(n3); sp2=(n2[g2]**2-n3[g2]**2)/2.0
cv=sp.std()/sp.mean();
# spacing correlation (consecutive): corr(sp1, sp2)
gg=(~np.isnan(n3))
s1=(n1[gg]**2-n2[gg]**2)/2; s2=(n2[gg]**2-n3[gg]**2)/2
rho=np.corrcoef(s1,s2)[0,1]
print(f"\nVar(N) exponent alpha={alpha:.2f} (alpha<0.75 => RIGID/GO)")
print(f"phase-spacing mean={sp.mean():.2f} (clockwork pi={np.pi:.2f}), CV={cv:.2f}")
print(f"consecutive-spacing correlation rho={rho:+.2f}  (rho<0 => ANTI-correlated => hyperuniform rigidity)")
print("VERDICT:", "RIGID (GO) confirmed" if alpha<0.75 else "Poisson-like (NO-GO)")
