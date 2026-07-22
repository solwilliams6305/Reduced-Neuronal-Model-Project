"""Falsification figure: the cusp gap determinant det(1-K^theta) is theta-DEPENDENT (two kernel forms),
so it is NOT the right complex-scaled object and does NOT reproduce FP W (skew +0.60). Pipeline verified on
Airy->TW2 (skew +0.20, flat). Diagnosis: W is a FIRST-PASSAGE law, not an eigenvalue-gap determinant."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
th=[0.45,0.55,0.65]
cd_sk=[1.272,1.522,1.103]; cd_ek=[0.381,1.987,0.771]
tj_sk=[0.701,0.000,1.284]; tj_ek=[-1.077,-1.973,-0.332]
fig,ax=plt.subplots(1,2,figsize=(12,4.5))
ax[0].plot(th,cd_sk,'o-',color='crimson',ms=8,label='CD kernel (psi,psi'')')
ax[0].plot(th,tj_sk,'s-',color='teal',ms=8,label='two-Jost kernel (u+,u-)')
ax[0].axhline(0.60,ls='--',color='k',lw=2,label='FP W target +0.60')
ax[0].axhline(0.224,ls=':',color='gray',label='TW2 (+0.22, pipeline OK)')
ax[0].set_xlabel('complex-scaling angle theta'); ax[0].set_ylabel('skew of det(1-K^theta)')
ax[0].set_title('skew is THETA-DEPENDENT => not the invariant object'); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
ax[1].plot(th,cd_ek,'o-',color='crimson',ms=8,label='CD kernel')
ax[1].plot(th,tj_ek,'s-',color='teal',ms=8,label='two-Jost kernel')
ax[1].axhline(-0.24,ls='--',color='k',lw=2,label='FP W target -0.24')
ax[1].set_xlabel('theta'); ax[1].set_ylabel('excess kurtosis'); ax[1].set_title('exk also theta-dependent, wrong')
ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
fig.suptitle('Stage-1 falsification: W != det(1-K^theta). The cusp edge is FIRST-PASSAGE (unbounded-below),\n'
             'not an eigenvalue-gap of a determinantal process (cf. Airy->TW2 which IS).',fontsize=11,fontweight='bold')
fig.subplots_adjust(left=0.07,right=0.98,bottom=0.12,top=0.82,wspace=0.22)
plt.savefig("figures/kernel_falsification.png",dpi=115); print("saved figures/kernel_falsification.png")
