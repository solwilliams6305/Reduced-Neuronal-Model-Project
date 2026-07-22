import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
fp=np.load("fp_persistence_data.npz"); mc=np.load("persistence_data.npz")
betas=[2,4,8,16]; cols={2:"#264653",4:"#2a9d8f",8:"#e9c46a",16:"#e76f51"}
fig,ax=plt.subplots(2,2,figsize=(11,8.4))

# (a) FP deep tail: -logS vs Theta, log-log, with FW Theta^{5/2} guide
a=ax[0,0]
for b in betas:
    Th=fp[f"Th_{b}"]; nlS=fp[f"nlS_{b}"]; m=(nlS>2)&(nlS<250)
    a.loglog(Th[m],nlS[m],"-",color=cols[b],lw=1.6,label=f"beta={b}")
tt=np.linspace(4,20,50); a.loglog(tt,0.19*tt**2.5,"k--",lw=1.1,label=r"$\propto\Theta^{5/2}$ (FW, $|Y^*|^5$)")
a.set_xlabel(r"$\Theta=Y^{*2}/2$"); a.set_ylabel(r"$-\log P(\Theta>t)$")
a.set_title("(a) MC-free deep tail (survival PDE): power law, NOT exponential"); a.legend(fontsize=8)

# (b) local power a(Theta) converging to ~5/2 (beta-independent)
a=ax[0,1]
for b in betas:
    Th=fp[f"Th_{b}"]; nlS=fp[f"nlS_{b}"]; m=(nlS>2)&(nlS<250)&(Th>3); Th=Th[m]; nlS=nlS[m]
    o=np.argsort(Th); Th=Th[o]; nlS=nlS[o]
    ce=[]; av=[]
    for lo in np.arange(4,Th.max()-3,2):
        mm=(Th>=lo)&(Th<=lo+3)
        if mm.sum()>4: ce.append(lo+1.5); av.append(np.polyfit(np.log(Th[mm]),np.log(nlS[mm]),1)[0])
    a.plot(ce,av,"o-",color=cols[b],ms=4,label=f"beta={b}")
a.axhline(2.5,color="k",ls="--",lw=1,label="FW 5/2"); a.set_ylim(2.0,3.4)
a.set_xlabel(r"depth $\Theta$"); a.set_ylabel(r"local power $d\log(-\log P)/d\log\Theta$")
a.set_title(r"(b) deep power $\to\sim5/2$ ($|Y^*|^5$), $\beta$-independent"); a.legend(fontsize=8)

# (c) amplitude A(beta) ~ beta^0.72 (anomalous; FW-naive beta^1)
a=ax[1,0]
bb=np.array([2.,3.,4.,6.,8.]); AA=np.array([0.1898,0.2596,0.3214,0.4208,0.5145])
a.loglog(bb,AA,"ks",ms=7,label=r"measured $A(\beta)$")
xx=np.linspace(2,8,30); a.loglog(xx,0.117*xx**0.715,"r--",lw=1.2,label=r"$\propto\beta^{0.72}$ (fit)")
a.loglog(xx,0.0949*xx**1.0,"b:",lw=1.1,label=r"$\propto\beta^{1}$ (FW-naive)")
a.set_xlabel(r"$\beta$"); a.set_ylabel(r"tail amplitude $A$  ($-\log P\sim A\,\Theta^{5/2}$)")
a.set_title(r"(c) amplitude anomaly: $A\sim\beta^{0.72}\neq\beta^{1}$"); a.legend(fontsize=8)

# (d) artifact demo: MC -logP is a power; linear 'rate' fit gives spurious window/beta-dependent slope
a=ax[1,1]
for b in [2,4,8]:
    Th=mc[f"Th_{b}"]; tg=np.linspace(2.6,6.5,20); cnt=np.array([(Th>t).sum() for t in tg]); ok=cnt>=30
    a.plot(tg[ok],-np.log(cnt[ok]/Th.size),"o",color=cols[b],ms=4,label=f"beta={b} (MC)")
    # spurious linear fit
    r=np.polyfit(tg[ok],-np.log(cnt[ok]/Th.size),1)
    a.plot(tg[ok],np.polyval(r,tg[ok]),"-",color=cols[b],lw=0.8,alpha=0.6)
a.set_xlabel(r"$\Theta$"); a.set_ylabel(r"$-\log P(\Theta>t)$")
a.set_title("(d) 'persistence rate' = linear fit to a convex power\n=> spurious window-dependent gamma (0.47-0.6): an artifact")
a.legend(fontsize=8)

fig.suptitle(r"W left tail (item 2):  $-\log P(Y^*<-s)\sim c\,\beta^{0.72}|s|^{5}$  — FW exponent 5 confirmed (MC-free); 'anomalous rate $\gamma$' retracted as fit artifact",fontsize=10)
fig.subplots_adjust(left=0.07,right=0.97,top=0.90,bottom=0.08,hspace=0.32,wspace=0.24)
fig.savefig("figures/persistence_item2.png",dpi=130)
print("saved figures/persistence_item2.png")
