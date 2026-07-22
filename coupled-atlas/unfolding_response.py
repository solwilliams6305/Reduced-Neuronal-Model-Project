"""
E — perturbation/unfolding response: how the edge law flows under the catastrophe unfolding.

q=3 (swallowtail) inner potential, unfolded by the lower-degree terms (codim 2):
    V(Y) = sign(Y)( |Y|^3 + c2 |Y|^2 + c1 |Y| ).
Near the turning Y=0 the LOWEST-degree active term dominates the escape, so:
  c1>0 (linear/fold term active) -> flow toward q=1 (fold/TW, skew~0.20)
  c2>0 (quadratic/cusp term)     -> flow toward q=2 (cusp,   skew~0.61)
  c1=c2=0                        -> pure q=3 (swallowtail,   skew~0.96)
This is Thom genericity / RG-relevance: higher catastrophes need lower terms tuned away.
numpy only.
"""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

def Vunf(Y, c1, c2):
    aY=np.abs(Y); return np.sign(Y)*(aY**3 + c2*aY**2 + c1*aY)

def run_one(c1,c2,N,eta,seed,Y0=3.0,Yend=-4.5,dt=4.5e-3):
    rng=np.random.default_rng(seed); nsteps=int(round((Y0-Yend)/dt)); sq=np.sqrt(dt)
    p=np.full(N, max(Vunf(Y0,c1,c2),0.0)**0.5); Ys=np.full(N,np.nan)
    for i in range(nsteps):
        Y=Y0-i*dt
        p+=(Vunf(Y,c1,c2)-p*p)*dt+eta*sq*rng.standard_normal(N)
        np.clip(p,-60,60,out=p); nw=np.isnan(Ys)&(p<-12.0); Ys[nw]=Y0-(i+1)*dt
    Ys[np.isnan(Ys)]=Yend
    x=(Ys-Ys.mean())/Ys.std(); m3=np.mean(x**3); m4=np.mean(x**4)
    return m3, m4-3.0

if __name__=="__main__":
    eta=np.sqrt(2.0); N=120000
    # c1 flow (c2=0): q3 -> fold
    c1s=[0.0,0.4,1.0,2.5,6.0]; sk1=[]; ek1=[]
    for j,c in enumerate(c1s):
        s,e=run_one(c,0.0,N,eta,seed=400+j); sk1.append(s); ek1.append(e)
        print(f"c1={c:>4.1f} c2=0.0 : skew {s:+.3f}  exk {e:+.3f}")
    # c2 flow (c1=0): q3 -> cusp
    c2s=[0.0,0.6,1.5,4.0,9.0]; sk2=[]; ek2=[]
    for j,c in enumerate(c2s):
        s,e=run_one(0.0,c,N,eta,seed=500+j); sk2.append(s); ek2.append(e)
        print(f"c1=0.0 c2={c:>4.1f} : skew {s:+.3f}  exk {e:+.3f}")
    fig,ax=plt.subplots(1,2,figsize=(11,4.3))
    ax[0].plot(c1s,sk1,'o-',color='crimson',ms=7,label='skew(c1), c2=0')
    ax[0].axhline(0.96,ls=':',color='green',label='q=3 (swallowtail) 0.96')
    ax[0].axhline(0.20,ls='--',color='navy',label='q=1 (fold/TW) 0.20')
    ax[0].set_xlabel('c1 (linear/fold unfolding term)'); ax[0].set_ylabel('skewness')
    ax[0].set_title('(A) c1 on: q=3 $\\to$ fold (q=1)'); ax[0].legend(fontsize=8); ax[0].set_xscale('symlog',linthresh=0.3)
    ax[1].plot(c2s,sk2,'s-',color='purple',ms=7,label='skew(c2), c1=0')
    ax[1].axhline(0.96,ls=':',color='green',label='q=3 0.96')
    ax[1].axhline(0.61,ls='--',color='teal',label='cusp $\\beta$=2 ref 0.61')
    ax[1].set_xlabel('c2 (quadratic/cusp unfolding term)'); ax[1].set_ylabel('skewness')
    ax[1].set_title('(B) c2 on: cusp CLASS, but $\\beta_{\\rm eff}$ rises with c2'); ax[1].legend(fontsize=8); ax[1].set_xscale('symlog',linthresh=0.3)
    fig.suptitle('Unfolding response: the LOWEST-degree term sets the class (Thom genericity); its coefficient also rescales $\\beta_{\\rm eff}$',fontsize=11.5)
    fig.subplots_adjust(left=0.07,right=0.98,bottom=0.13,top=0.88,wspace=0.22)
    plt.savefig("figures/unfolding_response.png",dpi=110); print("saved figures/unfolding_response.png")
