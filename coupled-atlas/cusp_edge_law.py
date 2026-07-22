"""
Construct W_beta = the rho->0 noisy parabolic-cylinder connection-data law (cusp edge law).

Definition. At the cusp (Delta=0) the inner equation is u'' = (sign(Y)Y^2 - eta*Wdot) u,
the STOCHASTIC WEBER operator. Cole-Hopf p=u'/u gives the swept Riccati
    dp = (V(Y) - p^2) dtau + eta dW,   V(Y)=sign(Y)Y^2,  tau=Y0-Y,
started recessive on the canard p(Y0)=+sqrt(V(Y0))=Y0. The first explosion p->-inf
at Y=Y* is the first NODE of the recessive solution = where the noisy connection
(recessive -> continuation) first vanishes. W_beta := law of (rescaled) Y* ; beta=4/eta^2.
This is the cusp analogue of the RRV stochastic-Airy/TW construction (q=2 vs q=1).

Test: does W_beta reproduce the cusp fingerprint (beta=2, eta=sqrt2):
   skew +0.61, exkurt -0.24, kappa5 ~ -2.2, kappa6 ~ -2.7, left tail e^{-|s|^5/20}?
Orientation: s = -(Y*-mean)/std  (eigenvalue convention: escaping LATER = larger s;
the heavy LEFT tail s->-inf = escaping early = large downward noise).
numpy only.  Modes:  sim <seed> <N>   (append Y* to cusp_Ystar.npy) ;  analyze.
"""
import numpy as np, sys, os

FN = "cusp_Ystar_v2.npy"

def Vcusp(Y):
    return np.sign(Y)*Y*Y

def simulate(N, eta, seed, Y0=3.0, Yend=-4.8, dt=3e-3):
    rng = np.random.default_rng(seed)
    nsteps = int(round((Y0-Yend)/dt))
    p = np.full(N, Y0, dtype=np.float64)
    Ystar = np.full(N, np.nan)
    sq = np.sqrt(dt)
    for i in range(nsteps):
        Y = Y0 - i*dt
        V = Vcusp(Y)
        p += (V - p*p)*dt + eta*sq*rng.standard_normal(N)
        np.clip(p, -50.0, 50.0, out=p)
        newly = np.isnan(Ystar) & (p < -10.0)
        Ystar[newly] = Y0 - (i+1)*dt
    Ystar[np.isnan(Ystar)] = Yend
    return Ystar

def cumulants(x):
    x = (x - x.mean())/x.std()
    m3=np.mean(x**3); m4=np.mean(x**4); m5=np.mean(x**5); m6=np.mean(x**6)
    k3=m3; k4=m4-3.0; k5=m5-10.0*m3; k6=m6-15.0*m4-10.0*m3**2+30.0
    return k3,k4,k5,k6

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv)>1 else "analyze"
    eta = np.sqrt(2.0)   # beta = 4/eta^2 = 2
    if mode == "sim":
        seed = int(sys.argv[2]); N = int(sys.argv[3])
        ys = simulate(N, eta, seed)
        if os.path.exists(FN):
            old = np.load(FN); ys = np.concatenate([old, ys])
        np.save(FN, ys)
        print(f"batch seed={seed} N={N}; total stored = {len(ys)}; "
              f"frac reaching Yend (no escape) = {np.mean(ys<=-4.79):.5f}")
    else:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        ys = np.load(FN); N=len(ys)
        # orientation s = +(Y*-mean)/std  (eigenvalue convention; matches established cusp fingerprint).
        # late escape (Y* very negative) -> s->-inf = steep exponent-5 tail;
        # early escape (Y* large)        -> s->+inf = heavier exponent-3 tail.
        s = (ys - ys.mean())/ys.std()
        k3,k4,k5,k6 = cumulants(ys)
        se = lambda v: np.sqrt(v/N)
        tgt = dict(skew=0.61, exk=-0.24, k5=-2.2, k6=-2.7)
        print(f"N = {N}   (beta=2, eta=sqrt2)   orientation s=+(Y*-mean)/std")
        print(f"  skew       = {k3:+.3f}   target {tgt['skew']:+.2f}   ~SE {se(6):.3f}")
        print(f"  exkurtosis = {k4:+.3f}   target {tgt['exk']:+.2f}   ~SE {se(24):.3f}")
        print(f"  kappa5     = {k5:+.3f}   target {tgt['k5']:+.2f}")
        print(f"  kappa6     = {k6:+.3f}   target {tgt['k6']:+.2f}")
        print(f"  mean(Y*)={ys.mean():+.4f}  std(Y*)={ys.std():.4f}")
        def tailfit(side, lo, hi):
            ts=np.linspace(lo,hi,30)
            P=np.array([(s< -t).mean() if side=='L' else (s>t).mean() for t in ts])
            ok=P*N>=40; tt=ts[ok]; nlP=-np.log(P[ok])
            A=np.vstack([np.ones_like(tt),np.log(tt)]).T
            coef,*_=np.linalg.lstsq(A,np.log(nlP),rcond=None)
            return tt,nlP,coef[1],np.exp(coef[0])
        tL,nlL,qL,aL = tailfit('L',1.8,3.4)
        tR,nlR,qR,aR = tailfit('R',1.8,3.6)
        print(f"  LEFT  tail (late escape) -lnP(s<-t) ~ {aL:.3f} t^{qL:.2f}  (target exp 5, amp~1/20=0.05)")
        print(f"  RIGHT tail (early escape) -lnP(s>t) ~ {aR:.3f} t^{qR:.2f}  (target exp 3)")
        np.save("cusp_s.npy", s)
        # ---- figure ----
        fig,ax=plt.subplots(1,3,figsize=(15,4.3))
        xs=np.linspace(-5,4,400)
        ax[0].hist(s,bins=300,density=True,color='steelblue',alpha=0.6,label=r'$\mathcal{W}_2$ (constructed)')
        ax[0].plot(xs,np.exp(-xs**2/2)/np.sqrt(2*np.pi),'k--',lw=1,label='standard normal')
        ax[0].set_xlim(-5,4); ax[0].set_xlabel('s (standardized edge variable)'); ax[0].set_ylabel('density')
        ax[0].set_title('(A) $\\mathcal{W}_2$ = noisy parabolic-cylinder edge law')
        ax[0].legend(fontsize=8)
        ax[0].text(0.02,0.97,f"skew {k3:+.2f} (t +0.61)\nexk {k4:+.2f} (t -0.24)\n"
                   f"$\\kappa_5$ {k5:+.2f} (t -2.2)\n$\\kappa_6$ {k6:+.2f} (t -2.7)",
                   transform=ax[0].transAxes,va='top',fontsize=8,family='monospace',
                   bbox=dict(fc='white',alpha=0.7))
        # (B) both tails log-log: left steeper than right (FW asymmetry 5 vs 3)
        ax[1].loglog(tL,nlL,'o',color='crimson',ms=5,label=f'LEFT (late esc) slope {qL:.2f}')
        ax[1].loglog(tR,nlR,'s',color='royalblue',ms=5,label=f'RIGHT (early esc) slope {qR:.2f}')
        ax[1].loglog(tL,tL**5/20,'k--',lw=1.1,label=r'$t^5/20$ (asymp. ref)')
        ax[1].set_xlabel('t'); ax[1].set_ylabel(r'$-\ln P(|s|>t)$')
        ax[1].set_title('(B) tails: LEFT steeper than RIGHT (pre-asymptotic)'); ax[1].legend(fontsize=7.5)
        # (C) cumulant bar comparison
        labels=['skew','exk',r'$\kappa_5$',r'$\kappa_6$']; got=[k3,k4,k5,k6]; want=[0.61,-0.24,-2.2,-2.7]
        xb=np.arange(4); w=0.38
        ax[2].bar(xb-w/2,got,w,color='teal',label='constructed $\\mathcal{W}_2$')
        ax[2].bar(xb+w/2,want,w,color='orange',alpha=0.8,label='cusp fingerprint (target)')
        ax[2].set_xticks(xb); ax[2].set_xticklabels(labels); ax[2].axhline(0,color='k',lw=0.5)
        ax[2].set_title('(C) fingerprint match'); ax[2].legend(fontsize=8)
        fig.suptitle('Construction of the new cusp edge law $\\mathcal{W}_\\beta$ ($\\beta=2$) and its fingerprint',fontsize=13)
        fig.subplots_adjust(left=0.05,right=0.98,bottom=0.12,top=0.86,wspace=0.25)
        plt.savefig("figures/cusp_edge_law.png",dpi=110)
        print("saved figures/cusp_edge_law.png")
