"""Consolidated rung-E figure: (A) cusp linear unfolding -> fold; (B) symmetry-breaking susceptibility;
(C) beta-response (fold/cusp inversion, reused). T2-consistent."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
d=np.load("unfold_v2.npy",allow_pickle=True).item()
fold,cusp=np.load("beta_axis_res.npy"); betas=[1,2,4]; TW=[0.293,0.224,0.165]
fig,ax=plt.subplots(1,3,figsize=(15,4.4))
# (A)
ax[0].plot(d["c1s"],d["skA"],'o-',color='crimson',ms=7)
ax[0].axhline(0.60,ls=':',color='purple',label='cusp $\\mathcal{W}_2$ (0.60)')
ax[0].axhline(0.20,ls=':',color='navy',label='fold TW$_2$ (0.20)')
ax[0].set_xlabel('$c_1$ (linear unfolding term)'); ax[0].set_ylabel('escape-loc skewness')
ax[0].set_title('(A) relevant unfolding: cusp $\\to$ fold'); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
# (B)
ax[1].plot(d["as_"],d["skB"],'s-',color='teal',ms=7,label='measured')
xs=np.linspace(-0.4,0.4,2); ax[1].plot(xs,0.60+d["suscep"]*xs,'--',color='gray',
    label=f'fit: d(skew)/da={d["suscep"]:+.2f}')
ax[1].axhspan(0.5,0.75,color='teal',alpha=0.08)
ax[1].set_xlabel('$a$ (symmetry-breaking: steepness asymmetry)'); ax[1].set_ylabel('skewness')
ax[1].set_title('(B) within-class symmetry-breaking (stays cusp)'); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
# (C)
ax[2].plot(betas,fold,'o-',color='navy',ms=7,label='fold TW$_\\beta$ (skew $\\downarrow$)')
ax[2].plot(betas,TW,'^--',color='navy',ms=6,alpha=0.5,label='TW$_\\beta$ ref')
ax[2].plot(betas,cusp,'s-',color='crimson',ms=7,label='cusp $\\mathcal{W}_\\beta$ (skew $\\uparrow$)')
ax[2].set_xlabel('$\\beta=4/\\eta^2$ (noise-symmetry)'); ax[2].set_ylabel('skewness'); ax[2].set_xticks(betas)
ax[2].set_title('(C) noise-symmetry $\\beta$: fold/cusp inversion'); ax[2].legend(fontsize=8); ax[2].grid(alpha=0.3)
fig.suptitle('Rung E — perturbation/unfolding responses of the cusp edge law (T2-consistent)',fontsize=12.5,fontweight='bold')
fig.subplots_adjust(left=0.06,right=0.985,bottom=0.12,top=0.88,wspace=0.27)
plt.savefig("figures/unfolding_response_v2.png",dpi=115); print("saved figures/unfolding_response_v2.png")
