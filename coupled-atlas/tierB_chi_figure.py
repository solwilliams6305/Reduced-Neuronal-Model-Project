"""Tier B: the susceptibility chi=dY*/dDelta — pathwise coupling + heavy tail (alpha~2) figure."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
d=np.load("tierB_chi.npz"); chi=d["chi1"]; dY1=d["dY1"]; dY2=d["dY2"]
c=chi-chi.mean(); a=np.abs(c)
fig,ax=plt.subplots(1,2,figsize=(11,4.3))
# (a) pathwise linearity: dY*(d2) vs dY*(d1) ~ slope 2
A=ax[0]; s=np.random.default_rng(0).choice(len(dY1),4000,replace=False)
A.plot(dY1[s],dY2[s],".",ms=2,alpha=0.3,color="#2a9d8f")
xx=np.array([dY1.min(),dY1.max()]); A.plot(xx,2*xx,"r-",lw=1.2,label="slope 2 (= d2/d1)")
A.set_xlabel(r"$\delta Y^*(\Delta{=}0.05)$"); A.set_ylabel(r"$\delta Y^*(\Delta{=}0.10)$")
A.set_title(r"(a) coupling is pathwise-linear: $\delta Y^*=\Delta\chi$"); A.legend(fontsize=8)
# (b) tail of |chi|: log-log survival, compare Gaussian, show alpha~2
A=ax[1]; xs=np.sort(a); surv=1-np.arange(len(xs))/len(xs)
A.loglog(xs,surv,"-",color="#264653",lw=1.5,label=r"$P(|\chi-E\chi|>x)$")
xr=np.logspace(np.log10(1),np.log10(30),20); A.loglog(xr,0.3*xr**(-2.0),"r--",lw=1.2,label=r"$x^{-2}$ (Hill $\alpha\approx2$)")
sd=a.std(); A.loglog(xr,np.exp(-xr**2/(2*sd**2)),":",color="gray",lw=1,label="Gaussian tail")
A.set_xlabel(r"$x$ (units of $\chi$)"); A.set_ylabel("survival"); A.set_ylim(1e-6,1.5)
A.set_title(r"(b) $\chi$ heavy tail $\alpha\approx2$: $E|\chi|<\infty$ (clean $W_1$), $E\chi^4=\infty$ (slow kurtosis)")
A.legend(fontsize=8)
fig.suptitle(r"Tier B coupling: $W_1(W_{\beta,\Delta},W_\beta)\leq\Delta\,E|\chi|=O(\varepsilon^{1/4})$; kurtosis limited by $\chi$'s $\alpha\approx2$ tail",fontsize=10)
fig.subplots_adjust(left=0.08,right=0.98,top=0.86,bottom=0.13,wspace=0.25)
fig.savefig("figures/tierB_chi.png",dpi=130); print("saved figures/tierB_chi.png")
