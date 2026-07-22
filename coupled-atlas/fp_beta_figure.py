"""Figure: the beta-family of smooth W_beta, the eta^1 std scaling, and the sigma-form confound."""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from fp_beta import cumulants, joint_R2, skewnormal_F

betas=[1.0,2.0,4.0,8.0]; curves=[]
for b in betas:
    y,F=np.load(f"fpb_{b:.1f}.npy"); eta=2.0/np.sqrt(b); m1,sd,sk,ek=cumulants(y,F)
    curves.append((eta,y,F,m1,sd,sk))

fig,ax=plt.subplots(1,3,figsize=(15,4.3)); cols=plt.cm.viridis(np.linspace(0.1,0.85,4))
for (eta,y,F,m1,sd,sk),b,cc in zip(curves,betas,cols):
    f=np.gradient(F,y); z=(y-m1)/sd
    ax[0].plot(z, f*sd, color=cc, lw=1.8, label=f'$\\beta$={b:.0f} (skew {sk:+.2f})')
ax[0].plot(np.linspace(-4,3,200),np.exp(-np.linspace(-4,3,200)**2/2)/np.sqrt(2*np.pi),'k:',lw=0.8,label='normal')
ax[0].set_xlim(-4,3); ax[0].set_xlabel('standardized s'); ax[0].set_ylabel('density')
ax[0].set_title(r'(A) $\mathcal{W}_\beta$ family (shape varies with $\beta$)'); ax[0].legend(fontsize=8)
# (B) std ~ eta^1
et=np.array([c[0] for c in curves]); sds=np.array([c[4] for c in curves])
p=np.polyfit(np.log(et),np.log(sds),1)[0]
ax[1].loglog(et,sds,'o',color='crimson',ms=9,label=f'data, slope {p:.2f}')
ax[1].loglog(et,sds[2]/et[2]*et,'k--',lw=1.2,label=r'$\propto\eta^1$')
ax[1].set_xlabel(r'noise $\eta$'); ax[1].set_ylabel(r'std$(\mathcal{W}_\beta)$')
ax[1].set_title(r'(B) escape-width $\sim\eta^{1}$ (not $\eta^{2/5}$)'); ax[1].legend(fontsize=8)
# (C) the confound
cv=[(eta,y,F,m1) for (eta,y,F,m1,sd,sk) in curves]
ps=[0.30,0.40,0.50,0.667]
rc=[joint_R2(cv,-2.09,pp) for pp in ps]
xg=np.linspace(-8,8,4000)
null=[(c[0],xg.copy(),skewnormal_F(al,xg),0.0) for c,al in zip(cv,[0.0,2.0,4.0,8.0])]
rn=[joint_R2(null,None,pp) for pp in ps]
ax[2].plot(ps,rc,'o-',color='teal',ms=8,label='cusp $\\mathcal{W}_\\beta$')
ax[2].plot(ps,rn,'s--',color='orange',ms=8,label='skew-normal NULL')
ax[2].set_ylim(0.99,1.001); ax[2].set_xlabel('scale exponent p ($\\eta^p$)'); ax[2].set_ylabel('joint shared-$\\sigma$-form $R^2$')
ax[2].set_title('(C) CONFOUND: null fits as well $\\Rightarrow$ no power'); ax[2].legend(fontsize=8)
fig.suptitle(r'$\beta$-family of $\mathcal{W}_\beta$: genuine varying family, but the $\sigma$-form regression is confounded (cusp $\approx$ null $\approx$1) $\Rightarrow$ PIV undecided by regression',fontsize=10.5)
fig.subplots_adjust(left=0.06,right=0.98,bottom=0.12,top=0.87,wspace=0.27)
plt.savefig("figures/fp_beta_family.png",dpi=115); print("saved figures/fp_beta_family.png")
print(f"cusp joint R^2 over p: {[round(r,4) for r in rc]}")
print(f"null joint R^2 over p: {[round(r,4) for r in rn]}")
