import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.special import loggamma
def V(Y,lam): return np.sign(Y)*Y*Y - lam
def scan(lam,Y0=6.5,Yend=-8.0):
    def rhs(Y,y): return [y[1],V(Y,lam)*y[0]]
    sol=solve_ivp(rhs,[Y0,Yend],[1.0,np.sqrt(max(V(Y0,lam),1e-9))],max_step=1.5e-3,rtol=1e-11,atol=1e-13,dense_output=True)
    Yg=np.linspace(Y0,Yend,24000); u,_=sol.sol(Yg); npos=0; neg=[]
    for i in range(1,len(Yg)):
        if u[i-1]*u[i]<0:
            Yc=Yg[i-1]+(u[i-1]/(u[i-1]-u[i]))*(Yg[i]-Yg[i-1])
            (neg.append(Yc) if Yc<0 else None); npos+= (Yc>0)
    return npos,(neg[0]**2/2 if neg else np.nan)
lams=np.arange(0.0,8.01,0.06); Th1=[]; npos=[]
for l in lams:
    a,b=scan(l); npos.append(a); Th1.append(b)
Th1=np.array(Th1); npos=np.array(npos)
fig,ax=plt.subplots(1,2,figsize=(11,4.2))
ax[0].plot(lams,Th1,".-",ms=3,color="#264653")
for p in [3,7]: ax[0].axvline(p,color="r",ls="--",lw=1)
ax[0].set_xlabel(r"$\lambda$"); ax[0].set_ylabel(r"$\Theta_1=Y^{*2}_{\rm det}/2$ (first node depth)")
ax[0].set_title(r"first node hits turning ($\Theta_1\!\to\!0$) at $\lambda=3,7$"+"\n= poles of $\\Gamma(3/4-\\lambda/4)$ (red)")
# |1/Gamma(3/4-lam/4)| = U(a,0) shape (zeros at 3,7,11)
lg=np.linspace(0,8,400); invG=np.real(np.exp(-loggamma(0.75-0.25*lg)))
ax[1].plot(lg,invG,color="#2a9d8f",lw=1.8,label=r"$1/\Gamma(3/4-\lambda/4)\propto U(a,0)$")
for p in [3,7]: ax[1].axvline(p,color="r",ls="--",lw=1)
ax[1].axhline(0,color="grey",lw=0.5)
ax[1].set_xlabel(r"$\lambda$"); ax[1].set_ylabel(r"$U(a,0)$ (recessive value at turning)")
ax[1].set_title(r"$U(a,0)=0$ (node at $Y=0$) at $\lambda=3,7,11$"); ax[1].legend(fontsize=9)
fig.suptitle(r"Confining-side connection coefficient of W's skeleton = c=1/Weber $\Gamma(3/4-\lambda/4)$ (pole test, non-confounded)",fontsize=10)
fig.subplots_adjust(left=0.08,right=0.97,top=0.83,bottom=0.13,wspace=0.28)
fig.savefig("figures/connection_poles.png",dpi=130); print("saved figures/connection_poles.png; transitions at lambda=3,7 = Gamma(3/4-lam/4) poles")
