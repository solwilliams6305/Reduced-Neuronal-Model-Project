"""
Characterize the rigid node process: (i) number-variance class (bounded/hyperuniform vs log/GUE) over a wider
depth, at beta=2,4,8 (beta-dependence); (ii) phase-spacing distribution (shape vs Wigner/Gaussian/Poisson).
Riccati node-counting (additive noise). numpy only.
"""
import numpy as np
def Vc(Y): return np.sign(Y)*Y*Y
def run(eta,N,Y0=3.0,Yend=-12.0,dt=2e-3,seed=11,thr=25.0,nrec=8):
    rng=np.random.default_rng(seed); n=int(round((Y0-Yend)/dt)); sq=np.sqrt(dt)
    p=np.full(N,np.sqrt(Vc(Y0))); cnt=np.zeros(N)
    cps=np.array([-1.5,-2.5,-4,-6,-8,-10,-11.8]); ci=0; Ncp=np.zeros((len(cps),N))
    nd=np.full((nrec,N),np.nan); idx=np.zeros(N,int)
    for i in range(n):
        Y=Y0-i*dt
        p=p+(Vc(Y)-p*p)*dt+eta*sq*rng.standard_normal(N); np.clip(p,-1e3,1e3,out=p)
        ex=p<-thr
        st=ex&(idx<nrec); cols=np.where(st)[0]
        if cols.size: nd[idx[cols],cols]=Y; idx[cols]+=1
        cnt=cnt+ex; p[ex]=thr
        if ci<len(cps) and Y<=cps[ci]: Ncp[ci]=cnt.copy(); ci+=1
    return cps*cps/2.0, Ncp, nd
print("(i) number variance vs depth, per beta  (bounded=>hyperuniform; log=>GUE-class)")
for b in [2.0,4.0,8.0]:
    Th,Ncp,nd=run(2.0/np.sqrt(b),8000)
    Nv=Ncp.var(axis=1); Nm=Ncp.mean(axis=1); ok=Nm>1.0
    blog=np.polyfit(np.log(Th[ok]),Nv[ok],1)[0]   # slope of Var vs log(Theta): ~0 bounded, >0 log
    print(f"  beta={b:.0f}: Var(N)=[{', '.join(f'{x:.2f}' for x in Nv)}] ; d Var/d log(Theta) = {blog:+.3f} (max Var/mean={np.max(Nv[ok]/Nm[ok]):.3f})")
    if b==2.0:
        sp=(nd[1:]**2-nd[:-1]**2)/2.0; sp=sp[np.isfinite(sp)]  # phase gaps (deeper node has larger Y^2)
        sp=sp[(sp>0)&(sp<8)]; m=sp.mean(); cv=sp.std()/m
        # shape: P(s->0)? Wigner ~ s e^{-s^2}; rigid Gaussian; Poisson e^{-s}. Check small-gap density.
        h,edges=np.histogram(sp/m,bins=24,range=(0,3),density=True)
        small=h[:3].mean()  # density near s=0 (repulsion=>~0)
        print(f"  [beta=2] phase-spacing mean={m:.2f}(pi={np.pi:.2f}) CV={cv:.2f}; density near s=0 ={small:.3f} (0=>repulsion; ~1 Poisson)")
