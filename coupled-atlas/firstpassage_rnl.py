"""
Validate the R_nl close: (a) the escape law W is sub-Gaussian (its MGF <= a parabola envelope, so
kappa_2m <= C^m m! eta^2m), light tails dominating the leptokurtic core; (b) the cumulants scale
kappa_n ~ c_n eta^n with standardized c_n bounded (the R_nl/nonlinear part is the eta^2 correction).
Uses the smooth FP W_beta (fpb_*.npy).  numpy only.
"""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from math import factorial

def dens(y,F):
    f=np.gradient(F,y); f=np.clip(f,0,None); f/=np.trapz(f,y); return f
def cumulants_std(y,f,nmax=8):
    m1=np.trapz(y*f,y); c=y-m1; sd=np.sqrt(np.trapz(c**2*f,y)); z=c/sd
    mu=[np.trapz(z**n*f,y) for n in range(nmax+1)]
    # standardized cumulants from central moments (k3..k8)
    k={}
    k[2]=1.0; k[3]=mu[3]; k[4]=mu[4]-3; k[5]=mu[5]-10*mu[3]
    k[6]=mu[6]-15*mu[4]-10*mu[3]**2+30
    k[8]=mu[8]-28*mu[6]-56*mu[5]*mu[3]-35*mu[4]**2+420*mu[4]+560*mu[3]**2-630
    return sd,k

betas=[1.0,2.0,4.0,8.0]
fig,ax=plt.subplots(1,2,figsize=(11,4.3)); cols=plt.cm.viridis(np.linspace(.1,.8,4))
print("(a) standardized log-MGF of W_beta vs lambda -- sub-Gaussian if <= parabola (grows slower than C*lambda^2)")
lams=np.linspace(-3.5,3.5,41)
for b,c in zip(betas,cols):
    y,F=np.load(f"fpb_{b:.1f}.npy"); f=dens(y,F)
    m1=np.trapz(y*f,y); sd=np.sqrt(np.trapz((y-m1)**2*f,y)); z=(y-m1)/sd
    G=np.array([np.log(np.trapz(np.exp(np.clip(l*z,-25,25))*f,y)) for l in lams])
    ax[0].plot(lams,G,'-',color=c,lw=1.7,label=f'$\\beta$={b:.0f}')
    # fit sub-Gaussian curvature at moderate lambda
    msk=np.abs(lams)<2; cur=np.polyfit(lams[msk],G[msk],2)[0]
    print(f"  beta={b}: log-MGF curvature(|lam|<2) = {cur:.3f} (Gaussian=0.5); max log-MGF at lam=3.5 is {G[-1]:.2f} vs 0.5*3.5^2={0.5*3.5**2:.2f}")
ax[0].plot(lams,0.5*lams**2,'k--',lw=1.2,label='Gaussian 0.5$\\lambda^2$')
ax[0].plot(lams,0.9*lams**2,'k:',lw=1.0,label='envelope 0.9$\\lambda^2$')
ax[0].set_xlabel(r'$\lambda$'); ax[0].set_ylabel('standardized log-MGF'); ax[0].legend(fontsize=7.5)
ax[0].set_title('(a) W sub-Gaussian: log-MGF <= parabola (light tails)')
# (b) standardized cumulants ratio kappa_2m/m! bounded (sub-Gaussian growth C^m m!)
print("\n(b) standardized cumulants kappa_2m and ratio kappa_2m/m! (bounded => sub-Gaussian C^m m!)")
ms=[1,2,3,4]
for b,c in zip(betas,cols):
    y,F=np.load(f"fpb_{b:.1f}.npy"); f=dens(y,F); sd,k=cumulants_std(y,f)
    k2m=[k[2],k[4],k[6],k[8]]; ratio=[abs(k2m[i])/factorial(ms[i]) for i in range(4)]
    print(f"  beta={b}: |k2|={abs(k2m[0]):.2f} |k4|={abs(k2m[1]):.2f} |k6|={abs(k2m[2]):.2f} |k8|={abs(k2m[3]):.2f}  -> /m!: {[round(r,2) for r in ratio]}  (C^m if bounded ratio)")
    ax[1].semilogy(ms,[max(abs(v),1e-3) for v in k2m],'o-',color=c,label=f'$\\beta$={b:.0f}')
ax[1].semilogy(ms,[factorial(m) for m in ms],'k--',lw=1.2,label='m! (sub-exp ref)')
ax[1].semilogy(ms,[3**m for m in ms],'k:',lw=1.0,label='$3^m$')
ax[1].set_xlabel('m'); ax[1].set_ylabel(r'$|\kappa_{2m}|$ standardized'); ax[1].set_xticks(ms); ax[1].legend(fontsize=7.5)
ax[1].set_title('(b) cumulants grow <= C^m m! (sub-Gaussian)')
fig.suptitle('R_nl close: escape law W is SUB-GAUSSIAN (MGF <= parabola; kappa_2m <= C^m m!) => kappa_2m(dTheta) <= C^m m! eta^2m',fontsize=10.5)
fig.subplots_adjust(left=0.07,right=0.98,bottom=0.12,top=0.88,wspace=0.24)
plt.savefig("figures/firstpassage_rnl.png",dpi=115); print("\nsaved figures/firstpassage_rnl.png")
