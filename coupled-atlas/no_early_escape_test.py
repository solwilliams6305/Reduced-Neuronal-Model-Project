"""
No-early-escape one-sided tail bound — numerical validation (theory-aligned).

Inner Riccati SDE:  dp = (V_D(Y) - p^2) dtau + eta dW,  Y = Y0 - tau,
V_D(Y)=sign(Y)|Y|(|Y|+D).  Canard pbar tracks +sqrt(V) (attracting branch);
repelling branch at p=-sqrt(V), i.e. fluctuation barrier dp_=p-pbar = -2pbar.
Fluctuation: d(dp_) = -2 pbar dp_ dtau - dp_^2 dtau + eta dW  (bistable cubic).

PHYSICAL early escape = p crosses below the repelling branch (-sqrt(V)) inside the
inner window Y in [Ycut,Yin] (i.e. peels off BEFORE the deterministic turning).
The barrier in sd units is h* = min_{inner} 2*pbar/sigma, sigma^2 = OU var
(dv/dtau=-4pbar v+eta^2). Sweeping eta sweeps h* (sigma ~ eta).

Scale-function/Kramers prediction (cubic potential Phi=pbar dp_^2 + dp_^3/3,
barrier height 4/3 pbar^3):
    -ln P_esc  ~  2*Phi_barrier/eta^2 = (8/3) pbar^3/eta^2 = h*^2 / 6.
TEST:  -ln P linear in h*^2 (sub-Gaussian), slope c ~ 1/6 and UNIFORM in rho=D/ell.
numpy only.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def Vfun(Y, D):
    aY = np.abs(Y); return np.sign(Y)*aY*(aY+D)

def canard(Yg, D):
    n=len(Yg); pb=np.empty(n); pb[0]=np.sqrt(max(Vfun(Yg[0],D),1e-12))
    for i in range(n-1):
        Y=Yg[i]; h=Yg[i+1]-Yg[i]; f=lambda Yv,pv: pv*pv-Vfun(Yv,D)
        k1=f(Y,pb[i]); k2=f(Y+0.5*h,pb[i]+0.5*h*k1)
        k3=f(Y+0.5*h,pb[i]+0.5*h*k2); k4=f(Y+h,pb[i]+h*k3)
        nx=pb[i]+(h/6)*(k1+2*k2+2*k3+k4)
        pb[i+1]=nx if (np.isfinite(nx) and 0<nx<1e6) else pb[i]
    return pb

def variance(Yg, pb, eta):
    n=len(Yg); v=np.empty(n); v[0]=eta**2/(4*max(pb[0],1e-6))
    for i in range(n-1):
        dtau=-(Yg[i+1]-Yg[i]); a=4*max(pb[i],0.0)
        if a>1e-8:
            e=np.exp(-a*dtau); v[i+1]=v[i]*e+(eta**2/a)*(1-e)
        else: v[i+1]=v[i]+eta**2*dtau
    return v

def sim_escape(Yg, sqrtV, D, eta, inner, N, seed):
    """return P(min over inner steps of (p+sqrtV) < 0) = P(crossed repelling branch)."""
    rng=np.random.default_rng(seed); n=len(Yg)
    dtau=-(Yg[1]-Yg[0]); sq=np.sqrt(dtau)
    pb0=np.sqrt(max(Vfun(Yg[0],D),1e-12)); p=pb0*np.ones(N)
    margin=np.full(N, 1e9)
    for i in range(n-1):
        V=Vfun(Yg[i],D)
        p=p+(V-p*p)*dtau+eta*sq*rng.standard_normal(N)
        np.clip(p,-50.0,50.0,out=p)
        if inner[i+1]:
            np.minimum(margin, p+sqrtV[i+1], out=margin)
    return (margin<0.0).mean()

def main():
    Y0,Yend,n = 3.0,0.02,700
    Ycut,Yin = 0.25,2.0
    Yg=np.linspace(Y0,Yend,n)
    inner=(Yg>=Ycut)&(Yg<=Yin)
    sqrtV={}; pbar={}
    Ds=[0.1,0.5,1.0,2.0]
    etas=[0.70,0.85,1.00,1.20,1.45]
    N=120000
    rows=[]
    print(f"window Y in[{Ycut},{Yin}]  N={N}   PREDICTION: -lnP ~ h*^2/6 (c~0.167), uniform in rho")
    print(f"{'D':>5}{'eta':>6}{'rho=D/l':>8}{'h*':>6}{'P_esc':>9}{'-lnP':>7}{'(-lnP)/h*^2':>12}")
    for D in Ds:
        pb=canard(Yg,D); pbar[D]=pb; sV=np.sqrt(np.maximum(Vfun(Yg,D),0.0)); sqrtV[D]=sV
        for k,eta in enumerate(etas):
            v=variance(Yg,pb,eta); sig=np.sqrt(v)
            hstar=np.min((2*pb/sig)[inner])
            P=sim_escape(Yg,sV,D,eta,inner,N,seed=700+k+37*int(10*D))
            ell=eta**(2/3); rho=D/ell
            nlP=-np.log(P) if P>0 else np.nan
            ratio=nlP/hstar**2 if P>0 else np.nan
            print(f"{D:>5.1f}{eta:>6.2f}{rho:>8.3f}{hstar:>6.2f}{P:>9.5f}{nlP:>7.2f}{ratio:>12.3f}")
            if P*N>=20: rows.append((D,rho,eta,hstar,P,nlP))
    rows=np.array(rows)
    # per-D fit -lnP = c * h*^2 (through origin) + report slope
    print(f"\n{'D':>5}{'rho-range':>14}{'slope c (h*^2)':>16}{'R2':>7}{'npts':>6}")
    fits={}
    for D in Ds:
        sel=rows[rows[:,0]==D]
        if len(sel)<3: continue
        h2=sel[:,3]**2; y=sel[:,5]
        c=np.sum(h2*y)/np.sum(h2*h2); pred=c*h2
        r2=1-np.sum((y-pred)**2)/np.sum((y-y.mean())**2)
        rr=sel[:,1]
        print(f"{D:>5.1f}{rr.min():>7.2f}-{rr.max():<6.2f}{c:>16.4f}{r2:>7.3f}{len(sel):>6d}")
        fits[D]=(c,r2,sel)
    cs=np.array([fits[D][0] for D in fits])
    print(f"\nUNIFORMITY: slope c across rho:  mean={cs.mean():.4f}  "
          f"spread={cs.max()-cs.min():.4f}  rel.spread={(cs.max()-cs.min())/cs.mean():.2f}")
    print(f"(predicted c=1/6={1/6:.3f}; pure-Kramers. Dynamic sweep shifts the constant but not the h*^2 law / uniformity.)")

    # ---- figure ----
    fig,ax=plt.subplots(1,3,figsize=(15,4.3))
    cols=plt.cm.viridis(np.linspace(0.1,0.85,len(Ds)))
    for D,cc in zip(Ds,cols):
        if D not in fits: continue
        c,r2,sel=fits[D]; h2=sel[:,3]**2; y=sel[:,5]; rho=sel[0,1]
        ax[0].plot(h2,y,'o',color=cc,ms=6,label=f"D={D}, rho={sel[:,1].mean():.2f}")
        xx=np.linspace(0,h2.max()*1.05,50); ax[0].plot(xx,c*xx,'-',color=cc,lw=1)
    xx=np.linspace(0,max(rows[:,3])**2*1.05,50)
    ax[0].plot(xx,xx/6,'k--',lw=1.5,label='Kramers 1/6')
    ax[0].set_xlabel(r'$h^{*2}$ (barrier in sd units)$^2$'); ax[0].set_ylabel(r'$-\ln P_{\rm esc}$')
    ax[0].set_title('(A) sub-Gaussian: $-\\ln P$ linear in $h^{*2}$'); ax[0].legend(fontsize=7)
    # panel B: P vs h* log scale (collapse)
    for D,cc in zip(Ds,cols):
        if D not in fits: continue
        _,_,sel=fits[D]; ax[1].semilogy(sel[:,3],sel[:,4],'o-',color=cc,ms=5,label=f"rho={sel[:,1].mean():.2f}")
    ax[1].set_xlabel(r'$h^*=\min\,2\bar p/\sigma$'); ax[1].set_ylabel(r'$P_{\rm esc}$ (log)')
    ax[1].set_title('(B) escape prob vs physical barrier'); ax[1].legend(fontsize=7)
    # panel C: slope c vs rho (uniformity)
    rr=[fits[D][2][:,1].mean() for D in fits]; cc2=[fits[D][0] for D in fits]
    ax[2].plot(rr,cc2,'s-',color='crimson',ms=8)
    ax[2].axhline(cs.mean(),ls=':',color='gray',label=f'mean={cs.mean():.3f}')
    ax[2].set_xlabel(r'$\rho=\Delta/\ell$ (through the merge)'); ax[2].set_ylabel(r'slope $c$')
    ax[2].set_title('(C) UNIFORMITY: $c(\\rho)$ ~ flat'); ax[2].set_xscale('log'); ax[2].legend(fontsize=8)
    ax[2].set_ylim(0,max(cc2)*1.4)
    fig.suptitle('No-early-escape one-sided tail:  P_esc <= C exp(-c h*^2),  c uniform in rho',fontsize=12)
    fig.subplots_adjust(left=0.06,right=0.98,bottom=0.12,top=0.86,wspace=0.27)
    plt.savefig("figures/no_early_escape_test.png",dpi=110)
    print("\nsaved figures/no_early_escape_test.png")

if __name__=="__main__":
    main()
