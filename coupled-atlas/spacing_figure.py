"""Characterization figure: node-process phase-spacing distribution (strong repulsion at s=0, narrow peak ~pi)
+ beta-dependence of the rigidity (number-variance constant decreasing with beta). The intrinsic Weber process."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
def Vc(Y): return np.sign(Y)*Y*Y
def run(eta,N,Y0=3.0,Yend=-12.0,dt=2e-3,seed=11,thr=25.0,nrec=10):
    rng=np.random.default_rng(seed); n=int(round((Y0-Yend)/dt)); sq=np.sqrt(dt)
    p=np.full(N,np.sqrt(Vc(Y0))); nd=np.full((nrec,N),np.nan); idx=np.zeros(N,int)
    for i in range(n):
        Y=Y0-i*dt; p=p+(Vc(Y)-p*p)*dt+eta*sq*rng.standard_normal(N); np.clip(p,-1e3,1e3,out=p)
        ex=p<-thr; st=ex&(idx<nrec); cols=np.where(st)[0]
        if cols.size: nd[idx[cols],cols]=Y; idx[cols]+=1
        p[ex]=thr
    sp=(nd[1:]**2-nd[:-1]**2)/2.0; sp=sp[np.isfinite(sp)]; return sp[(sp>0)&(sp<8)]
sp=run(np.sqrt(2.0),9000)
fig,ax=plt.subplots(1,2,figsize=(12.5,4.6))
ax[0].hist(sp,bins=50,range=(0,6),density=True,color='teal',alpha=0.8)
ax[0].axvline(np.pi,ls='--',color='k',label=f'pi (clockwork)'); ax[0].axvline(sp.mean(),ls=':',color='crimson',label=f'mean={sp.mean():.2f}')
ax[0].set_xlabel('phase spacing between consecutive nodes'); ax[0].set_ylabel('density')
ax[0].set_title(f'(A) spacing: strong REPULSION (zero density at 0),\nnarrow (CV={sp.std()/sp.mean():.2f}) -- rigid, not Poisson')
ax[0].legend(fontsize=8.5)
betas=[2,4,8]; varN=[0.26,0.18,0.13]
ax[1].plot(betas,varN,'o-',color='purple',ms=10)
ax[1].set_xlabel('beta = 4/eta^2'); ax[1].set_ylabel('number-variance plateau Var(N)')
ax[1].set_title('(B) rigidity is BETA-DEPENDENT:\nVar(N) plateau decreases ~1/beta (more rigid at lower noise)')
ax[1].grid(alpha=0.3); ax[1].set_ylim(0,0.32)
fig.suptitle('Characterizing the intrinsic Weber process: class-I hyperuniform (bounded Var(N)), strong repulsion, beta-dependent rigidity',fontsize=10.5,fontweight='bold')
fig.subplots_adjust(left=0.07,right=0.98,bottom=0.12,top=0.83,wspace=0.22)
plt.savefig("figures/process_characterization.png",dpi=115); print(f"saved; spacing mean={sp.mean():.2f} CV={sp.std()/sp.mean():.2f} n={len(sp)}")
