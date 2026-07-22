"""Figure: k=1 higher-Airy Fredholm det = TW2 (validated); naive k=2 != cusp W_2 (reframe not confirmed)."""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

d1=np.load("hoTW_F_k1.npy"); d2=np.load("hoTW_F_k2.npy")
s1,F1=d1; s2,F2=d2
def dens(s,F):
    F=np.clip(F,0,1); f=np.gradient(F,s); f=np.clip(f,0,None); f/=np.trapz(f,s); return f
f1=dens(s1,F1); f2=dens(s2,F2)
# standardize each
def std_curve(s,f):
    m=np.trapz(s*f,s); v=np.trapz((s-m)**2*f,s); sd=np.sqrt(v)
    return (s-m)/sd, f*sd
z1,g1=std_curve(s1,f1); z2,g2=std_curve(s2,f2)
cs=np.load("cusp_s.npy")  # W_2 samples (already standardized)

fig,ax=plt.subplots(1,2,figsize=(11,4.3))
xs=np.linspace(-5,4,400)
ax[0].plot(z1,g1,color='navy',lw=2,label='k=1 higher-Airy Fredholm det')
ax[0].plot(xs,np.exp(-(xs+0)**2/2)/np.sqrt(2*np.pi),'k:',lw=0.8)
ax[0].set_xlim(-5,4); ax[0].set_title('(A) VALIDATION: k=1 = TW$_2$'); ax[0].set_xlabel('standardized')
ax[0].text(0.02,0.97,"skew +0.225 (TW2 +0.224)\nexk +0.090 (+0.093)\ntail exp 2.97 (=3)",
           transform=ax[0].transAxes,va='top',fontsize=8,family='monospace',bbox=dict(fc='white',alpha=0.7))
ax[0].legend(fontsize=8)
h,e=np.histogram(cs,bins=200,density=True); ctr=0.5*(e[1:]+e[:-1])
ax[1].plot(ctr,h,color='crimson',lw=1.6,label='cusp $\\mathcal{W}_2$ (escape law)')
ax[1].plot(z2,g2,color='teal',lw=2,label='naive k=2 higher-Airy')
ax[1].plot(xs,np.exp(-xs**2/2)/np.sqrt(2*np.pi),'k:',lw=0.8,label='normal')
ax[1].set_xlim(-5,4); ax[1].set_title('(B) k=2 naive kernel $\\neq$ cusp $\\mathcal{W}_2$'); ax[1].set_xlabel('standardized')
ax[1].text(0.02,0.97,"cusp:  skew +0.61, exk -0.24\nnaive k2: skew -0.00, exk -0.05\n         tail exp 2.4 (not 5)\n=> wrong kernel (not multicritical)",
           transform=ax[1].transAxes,va='top',fontsize=7.5,family='monospace',bbox=dict(fc='white',alpha=0.7))
ax[1].legend(fontsize=8)
fig.suptitle('Higher-order Airy test: machinery validated (k=1=TW2); naive k=2 is NOT the cusp law (reframe unresolved)',fontsize=11)
fig.subplots_adjust(left=0.07,right=0.98,bottom=0.12,top=0.88,wspace=0.18)
plt.savefig("figures/hoTW_compare.png",dpi=115); print("saved figures/hoTW_compare.png")
