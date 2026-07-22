import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
d=np.load("node_kernel_data.npz")
k=d["kgrid"]; sbar=float(d["sbar"])
fig,ax=plt.subplots(2,2,figsize=(11,8))

# panel 1: S(k) with controls
a=ax[0,0]
a.semilogy(k,np.clip(d["S_pois"],1e-4,None),color="#888",lw=1.2,label="Poisson (control)")
a.semilogy(k,np.clip(d["S_j6"],1e-4,None),color="#f4a261",lw=1.0,label=r"jitter $\sigma=0.6$")
a.semilogy(k,np.clip(d["S_j3"],1e-4,None),color="#e9c46a",lw=1.0,label=r"jitter $\sigma=0.3$")
a.semilogy(k,np.clip(d["S_j1"],1e-4,None),color="#2a9d8f",lw=1.0,label=r"jitter $\sigma=0.1$")
a.semilogy(k,np.clip(d["S_data"],1e-4,None),color="#264653",lw=2.0,label="Weber nodes (DATA)")
a.axvline(2*np.pi/sbar,color="r",ls=":",lw=1,alpha=0.6)
a.set_xlabel("k"); a.set_ylabel("S(k)"); a.set_title("(a) structure factor: DATA at floor between Bragg\n=> lattice-type (class I), class-II/GUE excluded")
a.legend(fontsize=7,loc="lower right"); a.set_ylim(1e-4,60); a.set_xlim(0,3)

# panel 2: g2(r)
a=ax[0,1]; c=d["cent"]; g=d["g2"]
a.plot(c,g,color="#264653",lw=1.2)
for n in range(1,6): a.axvline(n*sbar,color="r",ls=":",lw=0.8,alpha=0.4)
a.axhline(1,color="#888",ls="--",lw=0.8)
a.set_xlabel(r"r  (phase-depth $\Theta$)"); a.set_ylabel(r"$g_2(r)$")
a.set_title(r"(b) pair correlation: hard core + $\pi$-lattice comb"); a.set_xlim(0,6*np.pi)

# panel 3: Var(spacing) ~ Theta^{-3/2}
a=ax[1,0]; D=d["depths"]; V=d["vars"]
a.loglog(D,V,"o",color="#264653",ms=7,label="measured")
xx=np.linspace(D.min(),D.max(),50)
a.loglog(xx,V[0]*(xx/D[0])**(-1.5),"r--",lw=1.2,label=r"$\Theta^{-3/2}$ (derived)")
a.set_xlabel(r"depth $\Theta$"); a.set_ylabel("Var(spacing)")
a.set_title(f"(c) jitter mechanism: Var(s)~$\\Theta^{{{float(d['q']):.2f}}}$  (derived $-3/2$)"); a.legend(fontsize=8)

# panel 4: beta collapse  (recompute cheaply)
def Vc(Y): return np.sign(Y)*Y*Y
def sample(beta,N,Y0=3.0,Yend=-15.0,dt=1.2e-3,seed=7,maxnodes=75):
    eta=2.0/np.sqrt(beta); rng=np.random.default_rng(seed)
    n=int(round((Y0-Yend)/dt)); sq=np.sqrt(dt); p=np.full(N,np.sqrt(max(Vc(Y0),1e-9)))
    nd=np.full((maxnodes,N),np.nan); idx=np.zeros(N,int)
    for i in range(n):
        Y=Y0-i*dt; thr=20*abs(Y)+50
        p+=(Vc(Y)-p*p)*dt+eta*sq*rng.standard_normal(N); np.clip(p,-1e4,1e4,out=p)
        ex=p<-thr; live=ex&(idx<maxnodes); cols=np.where(live)[0]
        if cols.size: nd[idx[cols],cols]=Y; idx[cols]+=1
        p[ex]=thr
    return nd
a=ax[1,1]; cols={2.0:"#264653",4.0:"#2a9d8f",8.0:"#e76f51"}
bands=[(4,10),(10,18),(18,28),(28,42),(42,58)]
for beta in [2.0,4.0,8.0]:
    Th=0.5*sample(beta,5000)**2; sp=Th[1:]-Th[:-1]; Thm=0.5*(Th[1:]+Th[:-1]); DD=[];VV=[]
    for lo,hi in bands:
        m=np.isfinite(sp[lo:hi])&(sp[lo:hi]>0.5)&(sp[lo:hi]<8)
        if m.sum()>300: DD.append(Thm[lo:hi][m].mean()); VV.append(sp[lo:hi][m].var())
    DD=np.array(DD); VV=np.array(VV)
    a.loglog(DD,VV*beta,"o-",color=cols[beta],ms=5,lw=1,label=f"beta={beta:.0f}")
xx=np.linspace(20,150,50); a.loglog(xx, 2.0*xx**(-1.5), "r--",lw=1.0,label=r"$2\,\Theta^{-3/2}$")
a.set_xlabel(r"depth $\Theta$"); a.set_ylabel(r"$\beta\cdot$Var(spacing)")
a.set_title(r"(d) $\beta$-collapse: Var(s)$\approx(2/\beta)\Theta^{-3/2}$"); a.legend(fontsize=8)

fig.suptitle(r"Intrinsic Weber node process = perturbed $\pi$-lattice (hyperuniform, class I); jitter $\sim(2/\beta)\Theta^{-3/2}$",fontsize=11)
fig.subplots_adjust(left=0.07,right=0.97,top=0.90,bottom=0.08,hspace=0.35,wspace=0.25)
fig.savefig("figures/node_kernel.png",dpi=130)
print("saved figures/node_kernel.png")
