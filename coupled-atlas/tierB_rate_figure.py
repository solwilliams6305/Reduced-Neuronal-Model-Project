import sys; sys.path.insert(0,"/Users/solomonwilliams/Reduced Neuronal Model Project/coupled-atlas")
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from tierB_rate_probe import cumulants_delta
Deltas=np.array([0.0,0.1,0.15,0.25,0.4,0.6,0.9])
C=np.array([cumulants_delta(D) for D in Deltas])  # mean,std,skew,exk
mean,std,skew,exk=C[:,0],C[:,1],C[:,2],C[:,3]
fig,ax=plt.subplots(1,2,figsize=(11,4.3))
a=ax[0]
a.plot(Deltas,skew,"o-",color="#2a9d8f",label="skew")
a.plot(Deltas,exk,"s-",color="#e76f51",label="excess kurtosis")
a.axhline(skew[0],color="#2a9d8f",ls=":",lw=1); a.axhline(exk[0],color="#e76f51",ls=":",lw=1)
a.set_xlabel(r"$\Delta$ (incomplete-merge, $\sim\varepsilon^{1/4}$)"); a.set_ylabel("cumulant")
a.set_title("(a) finite-merge $W_{2,\\Delta}\\to W_2$ as $\\Delta\\to0$\n(dotted = pure-cusp limit)"); a.legend(fontsize=8)
a=ax[1]
d_sk=np.abs(skew[1:]-skew[0]); d_ek=np.abs(exk[1:]-exk[0]); d_mn=np.abs(mean[1:]-mean[0]); Dv=Deltas[1:]
for d,c,lab in [(d_ek,"#e76f51","exk: $\\Delta^{0.8}$ (non-analytic!)"),(d_sk,"#2a9d8f","skew: $\\Delta^{\\gtrsim1}$"),(d_mn,"#264653","mean: $\\Delta^{1.15}$")]:
    a.loglog(Dv,d,"o-",color=c,ms=5,label=lab)
a.loglog(Dv,0.15*Dv**0.8,"--",color="#e76f51",lw=1); a.loglog(Dv,0.1*Dv**1.15,"--",color="#264653",lw=1)
a.set_xlabel(r"$\Delta$"); a.set_ylabel(r"$|{\rm cum}(\Delta)-{\rm cum}(0)|$")
a.set_title("(b) kurtosis converges SLOWEST & non-analytically\n=> rate-limiting p~0.2 (=0.8/4), explains exk overshoot"); a.legend(fontsize=8)
fig.suptitle("Tier B empirical rate: outer correction is kurtosis-limited, non-analytic, $p\\approx0.2$",fontsize=10.5)
fig.subplots_adjust(left=0.08,right=0.98,top=0.85,bottom=0.13,wspace=0.25)
fig.savefig("figures/tierB_rate.png",dpi=130); print("saved figures/tierB_rate.png")
