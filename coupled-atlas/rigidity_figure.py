"""Figure: the stochastic Weber node process is RIGID (GO). Var(N) stays bounded (~0.3-0.5) while the mean
count grows to ~12 -- far below the Poisson line Var=mean. Two independent methods (field, Riccati) agree."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
mF=[0.84,1.27,2.08,3.19,4.66,6.39,8.46,10.86]; vF=[0.347,0.471,0.468,0.479,0.448,0.487,0.483,0.450]
mR=[0.63,1.49,2.71,4.34,6.36,8.83,11.75]; vR=[0.263,0.265,0.269,0.254,0.259,0.271,0.277]
fig,ax=plt.subplots(1,2,figsize=(12.5,4.7))
ax[0].plot(mF,vF,'s-',color='navy',ms=7,label='field (multiplicative noise)')
ax[0].plot(mR,vR,'o-',color='crimson',ms=7,label='Riccati (additive noise)')
mm=np.linspace(0,12,50); ax[0].plot(mm,mm,'k--',lw=1.5,label='Poisson: Var = mean')
ax[0].set_xlabel('mean node count  N(Theta)'); ax[0].set_ylabel('Var( N )')
ax[0].set_title('(A) Var(N) stays BOUNDED (~0.3-0.5) vs Poisson (grows)\n=> RIGID node process (GO)')
ax[0].legend(fontsize=8.5); ax[0].grid(alpha=0.3)
ax[1].plot(mF,np.array(vF)/np.array(mF),'s-',color='navy',ms=7,label='field')
ax[1].plot(mR,np.array(vR)/np.array(mR),'o-',color='crimson',ms=7,label='Riccati')
ax[1].axhline(1,ls='--',color='k',lw=1.5,label='Poisson (=1)')
ax[1].set_xlabel('mean node count N(Theta)'); ax[1].set_ylabel('Var(N) / mean(N)')
ax[1].set_title('(B) Var/mean -> 0.02  (hyperuniform-class rigidity,\nstronger than GUE log; far from Poisson)')
ax[1].legend(fontsize=8.5); ax[1].grid(alpha=0.3); ax[1].set_ylim(0,1.1)
fig.suptitle('Piece-1 GO: the stochastic Weber node process is RIGID (1-D spectral rigidity of the swept operator).\n'
             '𝒲 = the EDGE (first node) of this rigid process; the intrinsic "Weber process" EXISTS.',fontsize=10.5,fontweight='bold')
fig.subplots_adjust(left=0.07,right=0.98,bottom=0.12,top=0.80,wspace=0.22)
plt.savefig("figures/node_rigidity.png",dpi=115); print("saved figures/node_rigidity.png")
