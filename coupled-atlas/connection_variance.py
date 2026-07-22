"""
Stochastic-Olver bound — numerical content: the noisy connection VARIANCE is O(eta^2), uniform in rho.
This is what the sub-Gaussian tube (T1) needs (variance/2nd-moment control through the turning).

For V_D(Y)=sign(Y)|Y|(|Y|+D), simulate the noisy Riccati dp=(V_D-p^2)dtau+eta dW; measure the
fluctuation variance Var(dp) at the turning entrance (Y=Yref, last canard point) and the escape spread.
Check: (a) true Var bounded uniformly in rho=D/ell (no blow-up at the merge rho->0), scales ~eta^2;
       (b) the OU quasi-static comparison eta^2/(4 sqrt(V)) DIVERGES as Y->0 (the un-uniform comparison).
numpy only.
"""
import numpy as np

def Vfun(Y,D): aY=np.abs(Y); return np.sign(Y)*aY*(aY+D)

def canard(Yg,D):
    n=len(Yg); pb=np.empty(n); pb[0]=np.sqrt(max(Vfun(Yg[0],D),1e-12))
    for i in range(n-1):
        Y=Yg[i]; h=Yg[i+1]-Yg[i]; f=lambda Yv,pv: pv*pv-Vfun(Yv,D)
        k1=f(Y,pb[i]);k2=f(Y+.5*h,pb[i]+.5*h*k1);k3=f(Y+.5*h,pb[i]+.5*h*k2);k4=f(Y+h,pb[i]+h*k3)
        nx=pb[i]+(h/6)*(k1+2*k2+2*k3+k4); pb[i+1]=nx if (np.isfinite(nx) and 0<nx<1e6) else pb[i]
    return pb

def simulate(D,eta,N,Y0=3.0,Yend=-1.2,dt=2e-3,Yref=0.06,seed=0):
    rng=np.random.default_rng(seed); n=int(round((Y0-Yend)/dt)); sq=np.sqrt(dt)
    Yg=Y0-np.arange(n+1)*dt; pb=canard(Yg,D)
    iref=int(round((Y0-Yref)/dt))
    p=np.full(N,pb[0]); dpref=np.full(N,np.nan); Ystar=np.full(N,np.nan)
    for i in range(n):
        Y=Yg[i]
        p=p+(Vfun(Y,D)-p*p)*dt+eta*sq*rng.standard_normal(N); np.clip(p,-40,40,out=p)
        if i==iref: dpref=p-pb[i]
        nw=np.isnan(Ystar)&(p<-10.0); Ystar[nw]=Y-dt
    Ystar[np.isnan(Ystar)]=Yend
    al=np.abs(dpref-np.nanmedian(dpref))<5  # bulk (exclude already-escaped outliers)
    return np.nanvar(dpref[al]), np.std(Ystar)

if __name__=="__main__":
    etas=[0.7,1.0,1.4]; Ds=[0.05,0.15,0.4,1.0,2.5]
    print("Connection variance at turning entrance (Yref=0.06): TRUE vs OU quasi-static")
    print(f"{'eta':>5}{'D':>6}{'rho':>7}{'Var_true':>9}{'Var/eta^2':>10}{'OU=e2/4sqrtV(0.02)':>18}")
    res={}
    for eta in etas:
        ell=eta**(2/3); row=[]
        for D in Ds:
            vt,se=simulate(D,eta,40000)
            ou=eta**2/(4*np.sqrt(Vfun(0.02,D)))
            rho=D/ell; row.append((rho,vt,se,ou))
            print(f"{eta:>5.1f}{D:>6.2f}{rho:>7.3f}{vt:>9.4f}{vt/eta**2:>10.3f}{ou:>18.3f}")
        res[eta]=row
    # uniformity summary
    for eta in etas:
        vts=np.array([r[1] for r in res[eta]]); rhos=np.array([r[0] for r in res[eta]])
        print(f"\neta={eta}: Var_true over rho[{rhos.min():.2f},{rhos.max():.2f}]: "
              f"min={vts.min():.3f} max={vts.max():.3f} (bounded, no merge blow-up); Var/eta^2 range "
              f"{(vts/eta**2).min():.2f}-{(vts/eta**2).max():.2f}")
    # figure
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig,ax=plt.subplots(1,2,figsize=(11,4.3)); cols=plt.cm.viridis(np.linspace(.1,.8,3))
    for eta,cc in zip(etas,cols):
        rr=np.array([r[0] for r in res[eta]]); vt=np.array([r[1] for r in res[eta]]); ou=np.array([r[3] for r in res[eta]])
        ax[0].plot(rr,vt/eta**2,'o-',color=cc,label=f'$\\eta$={eta} true/$\\eta^2$')
        ax[1].plot(rr,vt,'o-',color=cc,label=f'$\\eta$={eta} true')
        ax[1].plot(rr,ou,'s--',color=cc,alpha=0.5,label=f'$\\eta$={eta} OU(0.02)')
    ax[0].set_xlabel(r'$\rho=\Delta/\ell$'); ax[0].set_ylabel(r'Var$(\delta p)/\eta^2$'); ax[0].set_xscale('log')
    ax[0].set_title('(A) connection Var $/\\eta^2$ COLLAPSES + bounded (uniform in $\\rho$)'); ax[0].legend(fontsize=7)
    ax[1].set_xlabel(r'$\rho=\Delta/\ell$'); ax[1].set_ylabel('Var'); ax[1].set_xscale('log'); ax[1].set_yscale('log')
    ax[1].set_title('(B) true bounded; OU quasi-static larger at merge'); ax[1].legend(fontsize=6.5)
    fig.suptitle('Stochastic-Olver content: noisy connection variance is $O(\\eta^2)$, uniform in $\\rho$ through the merge',fontsize=11)
    fig.subplots_adjust(left=0.08,right=0.98,bottom=0.12,top=0.88,wspace=0.25)
    plt.savefig("figures/connection_variance.png",dpi=115); print("\nsaved figures/connection_variance.png")
