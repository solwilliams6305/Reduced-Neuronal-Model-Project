"""
Piece-1 GO/NO-GO: is the zero process of the stochastic Weber field RIGID (intrinsic Dyson/Airy2-type process,
number variance sub-linear) or POISSON-like (phase random-walks, number variance linear)?
Simulate u''=(sign(Y)Y^2 - eta xi)u (recessive at +inf), count sign changes (zeros) vs phase Theta=Y^2/2.
Var(N(Theta)) ~ Theta^alpha:  alpha<1 => RIGID (GO);  alpha~1 => Poisson-like (NO-GO).
Also spacing CV (rigid: <1; Poisson/exponential: =1). numpy only.
"""
import numpy as np
def Vc(Y): return np.sign(Y)*Y*Y
def run(eta,N,Y0=3.0,Yend=-9.0,dt=2e-3,seed=7):
    rng=np.random.default_rng(seed); n=int(round((Y0-Yend)/dt)); sq=np.sqrt(dt)
    u=np.ones(N); v=-3.0*np.ones(N); prev=np.sign(u); cnt=np.zeros(N)
    cps=np.arange(-1.0,-9.0,-1.0); ci=0; Ncp=np.zeros((len(cps),N)); firstY=np.full(N,np.nan)
    secondY=np.full(N,np.nan); pcnt=np.zeros(N)
    for i in range(n):
        Y=Y0-i*dt
        v=v+Vc(Y)*u*dt - eta*u*sq*rng.standard_normal(N)
        u=u+v*dt
        s=np.sign(u); flip=(s*prev<0)
        # record first two zeros (for spacing)
        nz=flip&np.isnan(firstY); firstY[nz]=Y
        nz2=flip&(~np.isnan(firstY))&np.isnan(secondY)&(pcnt>=1); secondY[nz2]=Y
        cnt=cnt+flip; pcnt=pcnt+flip
        prev=np.where(u!=0,s,prev)
        if ci<len(cps) and Y<=cps[ci]: Ncp[ci]=cnt.copy(); ci+=1
        m=np.abs(u)>1e6
        if m.any(): u[m]/=1e6; v[m]/=1e6
    Th=cps*cps/2.0
    return Th, Ncp, firstY, secondY
eta=np.sqrt(2.0)  # beta=2
Th,Ncp,firstY,secondY = run(eta,6000)
Nbar=Ncp.mean(axis=1); Nvar=Ncp.var(axis=1)
print("phase Theta | mean N | Var(N) | Var/mean (=1 Poisson, <1 rigid)")
for t,m,vv in zip(Th,Nbar,Nvar):
    if m>0.5: print(f"  {t:6.2f} | {m:5.2f} | {vv:6.3f} | {vv/m:5.2f}")
ok=Nbar>1.0
alpha=np.polyfit(np.log(Th[ok]),np.log(Nvar[ok]+1e-9),1)[0]
# spacing CV (phase spacing between 1st and 2nd zero)
good=~np.isnan(secondY)
sp=(firstY[good]**2-secondY[good]**2)/2.0  # phase gap = Th2-Th1
cv=sp.std()/sp.mean() if len(sp)>100 else np.nan
print(f"\nNumber-variance exponent alpha = {alpha:.2f}  (alpha~1 => Poisson/NO-GO; alpha<1 => RIGID/GO)")
print(f"phase-spacing mean={sp.mean():.2f} (clockwork pi={np.pi:.2f}), CV={cv:.2f}  (CV~1 Poisson; CV<<1 rigid clockwork)")
print("VERDICT:", "RIGID (GO) — intrinsic process exists" if alpha<0.75 else "POISSON-LIKE (NO-GO) — no intrinsic rigidity, zeros are a noisy renewal process")
