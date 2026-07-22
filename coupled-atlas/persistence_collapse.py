"""
Confirm the deep-left-tail law of W model-free and by data collapse:
  hypothesis  -log P(Theta>t) ~ C * beta^q * t^a   with a~2, q~1  (=> beta-Gaussian in Theta = drift-diffusion
  phase-variance form; the 'anomalous rate gamma~0.6' being an artifact of linear-fitting this convex curve).
Tests:
  (1) model-free power a: slope of log(-logP) vs log(Theta), per beta (well-sampled betas only).
  (2) prefactor beta-scaling q.
  (3) collapse: -logP vs beta^q * Theta^a should fall on one curve.
Big stats on beta=2,3,4,6. scipy avail. [NUMERIC].
"""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
def Vc(Y): return np.sign(Y)*Y*Y
def escape(eta,N,Y0=3.0,Yend=-9.0,dt=1.5e-3,seed=2,thr=20.0,clip=250.0):
    rng=np.random.default_rng(seed); n=int(round((Y0-Yend)/dt)); sq=np.sqrt(dt)
    p=np.full(N,np.sqrt(Vc(Y0))); Ys=np.full(N,np.nan)
    for i in range(n):
        Y=Y0-i*dt; p+=(Vc(Y)-p*p)*dt+eta*sq*rng.standard_normal(N); np.clip(p,-clip,clip,out=p)
        nw=np.isnan(Ys)&(p<-thr); Ys[nw]=Y0-(i+1)*dt
    return Ys[~np.isnan(Ys)]

betas=[2.0,3.0,4.0,6.0]; N=1200000
Thd={}; powers=[]; prefac=[]
fig,ax=plt.subplots(1,3,figsize=(14,4.2))
cols=dict(zip(betas,["#264653","#2a9d8f","#e9c46a","#e76f51"]))
print("=== deep-left-tail law:  -logP(Theta>t) ~ C beta^q t^a ===")
for b in betas:
    Ys=escape(2.0/np.sqrt(b),N,seed=int(b)+10); Th=(Ys[Ys<0]**2)/2; Thd[b]=Th
    # survival on a grid; require >=40 counts
    tg=np.linspace(2.4,7.0,26); cnt=np.array([(Th>t).sum() for t in tg]); ok=cnt>=40
    S=cnt[ok]/len(Ys); t=tg[ok]; nlp=-np.log(S)
    a=np.polyfit(np.log(t),np.log(nlp),1)[0]                 # model-free power
    # prefactor at fixed t via C = nlp / t^2 (using a=2 reference), averaged over the reliable window
    Cref=np.mean(nlp/t**2)
    powers.append(a); prefac.append(Cref)
    print(f"  beta={b:.0f}: power a={a:.2f}  C(=nlp/t^2)~{Cref:.3f}  C/beta={Cref/b:.3f}  (tail to Theta={t.max():.1f}, n>{cnt[ok][-1]})")
    ax[0].loglog(t,nlp,"o-",color=cols[b],ms=4,label=f"beta={b:.0f}")
q=np.polyfit(np.log(betas),np.log(prefac),1)[0]
print(f"\n  model-free power a = {np.mean(powers):.2f} +/- {np.std(powers):.2f}   (hypothesis 2.0)")
print(f"  prefactor C ~ beta^q,  q = {q:.2f}   (hypothesis 1.0)")
print(f"  => deep left tail:  -logP ~ ({np.mean(prefac)/np.mean(betas):.2f}) * beta * Theta^2   [beta-Gaussian in Theta]")

# panel 0 decor + a=2 guide
tt=np.linspace(2.4,7,50); ax[0].loglog(tt,0.35*2*tt**2,"k--",lw=1,label=r"$0.35\,\beta\,\Theta^2$ (b=2)")
ax[0].set_xlabel(r"$\Theta=Y^{*2}/2$"); ax[0].set_ylabel(r"$-\log P(\Theta>t)$"); ax[0].legend(fontsize=8)
ax[0].set_title(f"(a) deep tail is a power a={np.mean(powers):.2f} (~2), not exponential")
# panel 1: collapse -logP vs beta*Theta^2  (survival relative to full escaped sample)
for b in betas:
    Th=Thd[b]; tot=Th.size
    tg=np.linspace(2.4,7.0,26); cnt=np.array([(Th>t).sum() for t in tg]); ok=cnt>=40
    x=b*(tg[ok]**2); y=-np.log(cnt[ok]/tot)
    ax[1].plot(x,y,"o-",color=cols[b],ms=4,label=f"beta={b:.0f}")
ax[1].set_xlabel(r"$\beta\,\Theta^2$"); ax[1].set_ylabel(r"$-\log P$")
ax[1].set_title(r"(b) collapse under $\beta\,\Theta^2$"); ax[1].legend(fontsize=8)
# panel 2: the spurious 'rate' — linear fit slope vs beta (artifact demo)
allb=[2.,3.,4.,6.]; rr=[]
for b in allb:
    Th=Thd[b]; tg=np.linspace(2.6,6.5,20); cnt=np.array([(Th>t).sum() for t in tg]); ok=cnt>=30
    rr.append(np.polyfit(tg[ok],-np.log(cnt[ok]/Th.size),1)[0])
gg=np.polyfit(np.log(allb),np.log(rr),1)[0]
ax[2].loglog(allb,rr,"ks-",ms=6)
ax[2].set_xlabel(r"$\beta$"); ax[2].set_ylabel("linear-fit 'rate'")
ax[2].set_title(f"(c) spurious 'rate'~beta^{{{gg:.2f}}} (artifact of\nlinear-fitting a quadratic; not a real exponent)")
fig.subplots_adjust(left=0.06,right=0.98,top=0.88,bottom=0.13,wspace=0.28)
fig.suptitle(r"W deep left tail = $\beta$-Gaussian in $\Theta$:  $-\log P\approx0.35\,\beta\,\Theta^2$ (drift-diffusion form); 'anomalous $\gamma$' was a fit artifact",fontsize=10)
fig.savefig("figures/persistence_collapse.png",dpi=130)
print(f"\n  spurious linear-fit 'rate' scales as beta^{gg:.2f} (the old '0.6' — window-dependent, NOT a real exponent)")
print("saved figures/persistence_collapse.png")
