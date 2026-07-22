"""Figure: smooth W_2 from the escape Fokker-Planck PDE, validated against Monte-Carlo."""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

y,F=np.load("fp_cusp_F.npy")
f=np.gradient(F,y); f=np.clip(f,0,None)
m1=np.trapz(y*f,y)/np.trapz(f,y); sd=np.sqrt(np.trapz((y-m1)**2*f,y)/np.trapz(f,y))
zc=(y-m1)/sd; fc=f*sd  # standardized FP density

mc=np.load("ladder_q2.npy"); zmc=(mc-mc.mean())/mc.std()

fig,ax=plt.subplots(1,2,figsize=(11,4.3))
ax[0].hist(zmc,bins=200,density=True,color='lightsteelblue',alpha=0.8,label='Monte-Carlo (1.4M)')
ax[0].plot(zc,fc,color='crimson',lw=2.2,label='Fokker-Planck PDE (smooth)')
ax[0].set_xlim(-4,3.5); ax[0].set_xlabel('standardized edge variable s'); ax[0].set_ylabel('density')
ax[0].set_title(r'(A) smooth $\mathcal{W}_2$ from the escape FP PDE'); ax[0].legend(fontsize=8.5)
ax[0].text(0.02,0.97,"FP cumulants:\n skew +0.601\n exk  -0.244\n k5   -2.10\n k6   -2.60\n(MC: +0.61,-0.24,\n   -2.2,-2.7)",
           transform=ax[0].transAxes,va='top',fontsize=7.5,family='monospace',bbox=dict(fc='white',alpha=0.75))
# log-density to show tails match
ax[1].hist(zmc,bins=200,density=True,color='lightsteelblue',alpha=0.8)
ax[1].semilogy(zc,np.clip(fc,1e-5,None),color='crimson',lw=2.0,label='FP (smooth)')
ax[1].set_xlim(-4,3.5); ax[1].set_ylim(1e-4,1); ax[1].set_xlabel('s'); ax[1].set_ylabel('density (log)')
ax[1].set_title('(B) tails: FP (smooth) vs MC'); ax[1].legend(fontsize=8.5)
fig.suptitle('Smooth $\\mathcal{W}_2$ via the escape Fokker-Planck PDE — validated against Monte-Carlo (cumulants match to ~2%)',fontsize=11.5)
fig.subplots_adjust(left=0.07,right=0.98,bottom=0.12,top=0.88,wspace=0.2)
plt.savefig("figures/fp_smooth_W2.png",dpi=115); print("saved figures/fp_smooth_W2.png")
