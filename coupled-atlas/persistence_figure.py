"""Figure: the first-passage decomposition of W. (A) standardized tails at beta=2 -- right heavier (Kramers),
left lighter (persistence). (B) tail exponents vs beta -- beta-DEPENDENT (the W_beta family; no single universal
number). Plus the Gamow-bulk falsification (Rayleigh skew -0.63 vs W +0.60)."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
def Vc(Y): return np.sign(Y)*Y*Y
def escape(eta,N,Y0=3.0,Yend=-8.0,dt=6e-3,seed=2):
    rng=np.random.default_rng(seed); n=int(round((Y0-Yend)/dt)); sq=np.sqrt(dt)
    p=np.full(N,np.sqrt(Vc(Y0))); Ys=np.full(N,np.nan)
    for i in range(n):
        Y=Y0-i*dt; p+=(Vc(Y)-p*p)*dt+eta*sq*rng.standard_normal(N); np.clip(p,-60,60,out=p)
        nw=np.isnan(Ys)&(p<-12.0); Ys[nw]=Y0-(i+1)*dt
    return Ys[~np.isnan(Ys)]
Ys=escape(1.4,300000); s=(Ys-Ys.mean())/Ys.std()
S=np.linspace(0.8,3.0,18); PL=np.array([(s<-x).mean() for x in S]); PR=np.array([(s>x).mean() for x in S])
oL=PL*len(s)>=30; oR=PR*len(s)>=30
betas=[2.0,4.0,8.2]; Lexp=[1.99,1.85,1.68]; Rexp=[1.20,0.87,0.80]
fig,ax=plt.subplots(1,2,figsize=(12.5,4.7))
ax[0].plot(S[oL],-np.log(PL[oL]),'s-',color='navy',ms=5,label='LEFT (late) = persistence, lighter')
ax[0].plot(S[oR],-np.log(PR[oR]),'o-',color='crimson',ms=5,label='RIGHT (early) = Kramers, heavier')
ax[0].set_xlabel('standardized |s|'); ax[0].set_ylabel('-log P'); ax[0].legend(fontsize=8.5)
ax[0].set_title('(A) two mechanisms at beta=2 (skew +0.61):\nright Kramers (heavy) vs left persistence (light)')
ax[0].grid(alpha=0.3)
ax[1].plot(betas,Lexp,'s-',color='navy',ms=8,label='LEFT exp (persistence)')
ax[1].plot(betas,Rexp,'o-',color='crimson',ms=8,label='RIGHT exp (Kramers)')
ax[1].axhline(2,ls=':',color='gray',lw=1); ax[1].axhline(1,ls=':',color='gray',lw=1)
ax[1].set_xlabel('beta = 4/eta^2'); ax[1].set_ylabel('tail exponent'); ax[1].legend(fontsize=8.5)
ax[1].set_title('(B) exponents are BETA-DEPENDENT => W_beta family,\nno single universal number (Gamow rate 1.64 absent)')
ax[1].grid(alpha=0.3)
fig.suptitle('W = first-passage law: right tail Kramers + left tail PERSISTENCE (named), bulk = beta-family.\n'
             'Gamow survival FALSIFIED for the bulk (dominant-resonance Rayleigh skew -0.63 vs W +0.60).',fontsize=10.5,fontweight='bold')
fig.subplots_adjust(left=0.07,right=0.98,bottom=0.12,top=0.80,wspace=0.22)
plt.savefig("figures/persistence_decomposition.png",dpi=115); print("saved figures/persistence_decomposition.png")
