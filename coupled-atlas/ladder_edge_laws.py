"""
Higher-ladder edge laws — escape fingerprints across the catastrophe degree q.

Inner equation u'' = (sign(Y)|Y|^q - eta Wdot) u  (swept multicritical potential).
Riccati p=u'/u:  dp = (V_q(Y) - p^2) dtau + eta dW,  V_q=sign(Y)|Y|^q, tau=Y0-Y,
started recessive p(Y0)=+sqrt(V_q(Y0)). First explosion Y* = the edge variable.
W^(q)_beta := law of standardized Y* (orientation s=+(Y*-mean)/std). beta=4/eta^2=2.

Reframe test (catastrophe q = multicritical/PII-hierarchy index k):
  left  tail exponent  2q+1   (q=1->3 TW, q=2->5 cusp, q=3->7, q=4->9)
  right tail exponent  3q/2
  cumulants: a monotone progression of higher-order-TW laws.
numpy only.  modes: sim <q> <seed> <N> ;  analyze
"""
import numpy as np, sys, os

def Vq(Y, q):
    return np.sign(Y)*np.abs(Y)**q

def simulate(q, N, eta, seed, Y0=3.0, Yend=-4.5, dt=4e-3):
    rng = np.random.default_rng(seed)
    nsteps = int(round((Y0-Yend)/dt)); sq = np.sqrt(dt)
    p = np.full(N, Vq(Y0,q)**0.5, dtype=np.float64)
    Ystar = np.full(N, np.nan)
    for i in range(nsteps):
        Y = Y0 - i*dt
        p += (Vq(Y,q) - p*p)*dt + eta*sq*rng.standard_normal(N)
        np.clip(p, -60.0, 60.0, out=p)
        newly = np.isnan(Ystar) & (p < -12.0)
        Ystar[newly] = Y0 - (i+1)*dt
    Ystar[np.isnan(Ystar)] = Yend
    return Ystar

def cumulants(x):
    x = (x - x.mean())/x.std()
    m3=np.mean(x**3); m4=np.mean(x**4); m5=np.mean(x**5); m6=np.mean(x**6)
    return m3, m4-3.0, m5-10.0*m3, m6-15.0*m4-10.0*m3**2+30.0

def tailfit(s, side, lo, hi, N):
    ts=np.linspace(lo,hi,30)
    P=np.array([(s< -t).mean() if side=='L' else (s>t).mean() for t in ts])
    ok=P*N>=40; tt=ts[ok]
    if tt.size<4: return np.nan
    nlP=-np.log(P[ok]); A=np.vstack([np.ones_like(tt),np.log(tt)]).T
    coef,*_=np.linalg.lstsq(A,np.log(nlP),rcond=None); return coef[1]

if __name__=="__main__":
    mode=sys.argv[1]; eta=np.sqrt(2.0)
    if mode=="sim":
        q=int(sys.argv[2]); seed=int(sys.argv[3]); N=int(sys.argv[4])
        fn=f"ladder_q{q}.npy"
        ys=simulate(q,N,eta,seed)
        if os.path.exists(fn): ys=np.concatenate([np.load(fn),ys])
        np.save(fn,ys)
        print(f"q={q} seed={seed} N={N}: total {len(ys)}; noesc frac {np.mean(ys<=-4.49):.5f}; "
              f"mean(Y*) {ys.mean():+.3f} std {ys.std():.3f}")
    else:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        qs=[1,2,3,4]; rows=[]; dens={}
        print(f"{'q':>3}{'k=q':>5}{'N':>9}  {'skew':>7}{'exk':>7}{'k5':>8}{'k6':>8}  "
              f"{'Ltail':>7}{'2q+1':>5}  {'Rtail':>7}{'3q/2':>5}")
        for q in qs:
            fn=f"ladder_q{q}.npy"
            if not os.path.exists(fn): continue
            ys=np.load(fn); N=len(ys); s=(ys-ys.mean())/ys.std()
            k3,k4,k5,k6=cumulants(ys); qL=tailfit(s,'L',1.8,3.4,N); qR=tailfit(s,'R',1.8,3.6,N)
            print(f"{q:>3}{q:>5}{N:>9}  {k3:>+7.3f}{k4:>+7.3f}{k5:>+8.3f}{k6:>+8.3f}  "
                  f"{qL:>7.2f}{2*q+1:>5}  {qR:>7.2f}{1.5*q:>5.1f}")
            rows.append((q,k3,k4,k5,k6,qL,qR)); dens[q]=s
        R=np.array(rows)
        # figure
        fig,ax=plt.subplots(1,3,figsize=(15,4.3))
        cols=plt.cm.viridis(np.linspace(0.1,0.85,len(dens)))
        xs=np.linspace(-5,4,400)
        for (q,s),cc in zip(dens.items(),cols):
            h,e=np.histogram(s,bins=250,density=True); ctr=0.5*(e[1:]+e[:-1])
            ax[0].plot(ctr,h,color=cc,lw=1.4,label=f'q={q} (k={q})')
        ax[0].plot(xs,np.exp(-xs**2/2)/np.sqrt(2*np.pi),'k--',lw=0.9,label='normal')
        ax[0].set_xlim(-5,4); ax[0].set_xlabel('s (standardized edge)'); ax[0].set_ylabel('density')
        ax[0].set_title('(A) ladder of edge laws $\\mathcal{W}^{(q)}_2$'); ax[0].legend(fontsize=8)
        ax[1].plot(R[:,0],R[:,1],'o-',color='teal',label='skew'); ax[1].plot(R[:,0],R[:,2],'s-',color='orange',label='exk')
        ax[1].plot(R[:,0],R[:,3],'^-',color='crimson',label=r'$\kappa_5$'); ax[1].plot(R[:,0],R[:,4],'d-',color='purple',label=r'$\kappa_6$')
        ax[1].axhline(0,color='k',lw=0.4); ax[1].set_xlabel('catastrophe degree q'); ax[1].set_ylabel('standardized cumulant')
        ax[1].set_title('(B) cumulant progression vs q'); ax[1].legend(fontsize=8); ax[1].set_xticks(qs)
        ax[2].plot(R[:,0],R[:,5],'o',color='crimson',ms=8,label='left tail (fit)')
        ax[2].plot(R[:,0],2*R[:,0]+1,'--',color='crimson',lw=1,label='2q+1 (asymp.)')
        ax[2].plot(R[:,0],R[:,6],'s',color='royalblue',ms=8,label='right tail (fit)')
        ax[2].plot(R[:,0],1.5*R[:,0],'--',color='royalblue',lw=1,label='3q/2 (asymp.)')
        ax[2].set_xlabel('catastrophe degree q'); ax[2].set_ylabel('tail exponent')
        ax[2].set_title('(C) tail exponents (fit = pre-asymptotic)'); ax[2].legend(fontsize=8); ax[2].set_xticks(qs)
        fig.suptitle('The catastrophe ladder of edge laws (q=k = multicritical / higher-order-TW index)',fontsize=12)
        fig.subplots_adjust(left=0.05,right=0.98,bottom=0.12,top=0.87,wspace=0.24)
        plt.savefig("figures/ladder_edge_laws.png",dpi=110); print("saved figures/ladder_edge_laws.png")
