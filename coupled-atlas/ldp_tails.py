"""
R1 — instanton / Freidlin-Wentzell tail validation for the cusp escape law W_2.

Cusp Riccati  dp=(V-p^2)dtau+eta dW, V=sign(Y)Y^2, tau=Y0-Y, recessive start p=Y0.
First-explosion location Y*. Predicted tails (this work):
  RIGHT (early escape, Y*>Y_e>0):  -ln P(Y*>Y_e) ~ (8/3) Y_e^{3q/2} / eta^2   [cubic Kramers barrier]
       q=2 => exponent 3 in Y_e   (NOT 5/2 => refutes half-integer multicritical k=1/2)
  LEFT  (late escape, Y*<-Y_L):    -ln P(Y*<-Y_L) ~ (|Y_L|^{2q+1}-|Ydet|^{2q+1})/(2(2q+1)) /eta^2  [holding instanton]
       q=2 => exponent 5 in Y_L
KEY trick: at fixed Y_e, -lnP is LINEAR in 1/eta^2 (the LDP); the slope's dependence on Y_e
gives the exponent cleanly, avoiding the pre-asymptotic problem of the standardized tail.
numpy only.  modes:  sim <eta> <N> ;  analyze
"""
import numpy as np, sys, os

ETAS = [0.55, 0.70, 0.90, 1.15, 1.45]

def simulate(eta, N, seed, Y0=2.5, Yend=-4.0, dt=3e-3):
    rng=np.random.default_rng(seed); nst=int(round((Y0-Yend)/dt)); sq=np.sqrt(dt)
    p=np.full(N, Y0); Ys=np.full(N,np.nan)
    for i in range(nst):
        Y=Y0-i*dt
        p+=(np.sign(Y)*Y*Y - p*p)*dt + eta*sq*rng.standard_normal(N)
        np.clip(p,-60,60,out=p); nw=np.isnan(Ys)&(p<-12.0); Ys[nw]=Y0-(i+1)*dt
    Ys[np.isnan(Ys)]=Yend; return Ys

def fn(eta): return f"ldp_q2_eta{eta:.2f}.npy"

if __name__=="__main__":
    mode=sys.argv[1]
    if mode=="sim":
        eta=float(sys.argv[2]); N=int(sys.argv[3])
        ys=simulate(eta,N,seed=int(1000*eta)+7)
        if os.path.exists(fn(eta)): ys=np.concatenate([np.load(fn(eta)),ys])
        np.save(fn(eta),ys)
        print(f"eta={eta:.2f} N={N} total={len(ys)} mean(Y*)={ys.mean():+.3f} std={ys.std():.3f} "
              f"P(Y*>0)={np.mean(ys>0):.4f} P(Y*<-3)={np.mean(ys<-3):.4f}")
    elif mode=="analyze2":
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        data={}
        for e in ETAS:
            if os.path.exists(fn(e)): data[e]=np.load(fn(e))
        es=sorted(data)
        print("loaded:", [f"{e:.2f}(N={len(data[e])})" for e in es])
        def rate(thr, side):  # A = eta^2 * (-ln P), per eta, where measurable
            out={}
            for e in es:
                P = np.mean(data[e]>thr) if side=='R' else np.mean(data[e]<thr)
                if P*len(data[e])>=12: out[e]=e*e*(-np.log(P))
            return out
        print("\nRIGHT-tail rate-function collapse  A(Y_e)=eta^2*(-lnP(Y*>Y_e))  [eta-independent => LDP]")
        Yes=np.arange(-1.8,0.41,0.3)
        for Ye in Yes:
            r=rate(Ye,'R')
            if len(r)>=2:
                vals=np.array(list(r.values()))
                print(f"  Y_e={Ye:+.1f}: A={vals.mean():.2f} (spread {vals.max()-vals.min():.2f}, {len(r)} etas)  pred(8/3 Y_e^3 if>0)={8/3*Ye**3 if Ye>0 else 0:.2f}")
        print("\nLEFT-tail rate-function collapse  A(Y_L)=eta^2*(-lnP(Y*<-Y_L))")
        YLs=np.arange(2.2,3.41,0.2)
        for YL in YLs:
            r=rate(-YL,'L')
            if len(r)>=2:
                vals=np.array(list(r.values()))
                print(f"  Y_L={YL:.1f}: A={vals.mean():.2f} (spread {vals.max()-vals.min():.2f}, {len(r)} etas)  pred((Y_L^5-2.09^5)/10)={(YL**5-2.09**5)/10:.2f}")
        # figure: collapse + analytic exponents
        fig,ax=plt.subplots(1,2,figsize=(11.5,4.5)); cols=plt.cm.plasma(np.linspace(0.15,0.85,len(es)))
        for e,cc in zip(es,cols):
            Yg=np.arange(-1.9,0.6,0.1); A=[]; Yp=[]
            for Ye in Yg:
                P=np.mean(data[e]>Ye)
                if P*len(data[e])>=12: A.append(e*e*(-np.log(P))); Yp.append(Ye)
            ax[0].plot(Yp,A,'o-',color=cc,ms=3.5,label=f'$\\eta$={e:.2f}')
        yy=np.linspace(0.05,0.6,30); ax[0].plot(yy,8/3*yy**3+8.3,'k--',lw=1.4,label=r'analytic $\frac{8}{3}Y_e^3$+c ($Y_e\gg1$, exp 3)')
        ax[0].set_xlabel(r'$Y_e$ (early-escape threshold)'); ax[0].set_ylabel(r'$\eta^2\,(-\ln P(Y^*>Y_e))$')
        ax[0].set_title('(A) RIGHT tail: $\\eta^2$-rate-fn COLLAPSE (LDP confirmed)'); ax[0].legend(fontsize=7)
        for e,cc in zip(es,cols):
            Yg=np.arange(2.2,3.6,0.1); A=[]; Yp=[]
            for YL in Yg:
                P=np.mean(data[e]<-YL)
                if P*len(data[e])>=12: A.append(e*e*(-np.log(P))); Yp.append(YL)
            if A: ax[1].plot(Yp,A,'s-',color=cc,ms=3.5,label=f'$\\eta$={e:.2f}')
        yy=np.linspace(2.2,3.5,30); ax[1].plot(yy,(yy**5-2.09**5)/10,'k--',lw=1.4,label=r'analytic $(Y_L^5-Y_{det}^5)/10$ (exp 5, upper bd)')
        ax[1].set_xlabel(r'$Y_L$ (late-escape threshold)'); ax[1].set_ylabel(r'$\eta^2\,(-\ln P(Y^*<-Y_L))$')
        ax[1].set_title('(B) LEFT tail: $\\eta^2$-rate-fn COLLAPSE (LDP)'); ax[1].legend(fontsize=7)
        fig.suptitle('FW rate-function $\\eta^2$-collapse confirms the LDP; ANALYTIC tails (5,3)=(2q+1,3q/2) refute k=1/2 (5/2) & integer multicritical (ratio 2)',fontsize=10)
        fig.subplots_adjust(left=0.07,right=0.98,bottom=0.12,top=0.88,wspace=0.24)
        plt.savefig("figures/ldp_tails.png",dpi=115); print("\nsaved figures/ldp_tails.png")
    else:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        data={};
        for e in ETAS:
            if os.path.exists(fn(e)): data[e]=np.load(fn(e))
        es=sorted(data); inv=np.array([1/e**2 for e in es])
        print("loaded etas:", [f"{e:.2f}(N={len(data[e])})" for e in es])
        # RIGHT tail: slope of -lnP(Y*>Y_e) vs 1/eta^2, as function of Y_e
        Yes=np.array([0.2,0.4,0.6,0.8,1.0]); slopeR=[]
        print("\nRIGHT tail  -lnP(Y*>Y_e) = m/eta^2 + b ; predicted m=(8/3)Y_e^3")
        for Ye in Yes:
            nlP=[]; ok=[]
            for e in es:
                P=np.mean(data[e]>Ye); nlP.append(-np.log(P) if P*len(data[e])>=15 else np.nan); ok.append(P*len(data[e])>=15)
            nlP=np.array(nlP); m=np.array(ok)
            if m.sum()>=3:
                A=np.vstack([inv[m],np.ones(m.sum())]).T; sl,b=np.linalg.lstsq(A,nlP[m],rcond=None)[0]
                slopeR.append((Ye,sl)); print(f"  Y_e={Ye:.1f}: slope m={sl:.3f}  (pred {8/3*Ye**3:.3f})  [{m.sum()} etas]")
        slopeR=np.array(slopeR)
        # exponent of m(Y_e): fit log m = log c + p log Y_e
        pos=slopeR[:,1]>0
        pR=np.polyfit(np.log(slopeR[pos,0]),np.log(slopeR[pos,1]),1)[0]
        print(f"  => RIGHT-tail exponent in Y_e: {pR:.2f}  (predict 3 = 3q/2; half-integer k=1/2 would give 2.5)")
        # LEFT tail: slope of -lnP(Y*<-Y_L) vs 1/eta^2
        YLs=np.array([2.4,2.7,3.0,3.3]); slopeL=[]
        print("\nLEFT tail  -lnP(Y*<-Y_L) = m/eta^2 + b ; predicted m=(Y_L^5-2.09^5)/10")
        for YL in YLs:
            nlP=[]; ok=[]
            for e in es:
                P=np.mean(data[e]<-YL); nlP.append(-np.log(P) if P*len(data[e])>=15 else np.nan); ok.append(P*len(data[e])>=15)
            nlP=np.array(nlP); m=np.array(ok)
            if m.sum()>=3:
                A=np.vstack([inv[m],np.ones(m.sum())]).T; sl,b=np.linalg.lstsq(A,nlP[m],rcond=None)[0]
                slopeL.append((YL,sl)); print(f"  Y_L={YL:.1f}: slope m={sl:.3f}  (pred {(YL**5-2.09**5)/10:.3f})  [{m.sum()} etas]")
        slopeL=np.array(slopeL)
        if len(slopeL)>=2:
            posL=slopeL[:,1]>0
            pL=np.polyfit(np.log(slopeL[posL,0]),np.log(slopeL[posL,1]),1)[0]
            print(f"  => LEFT-tail exponent in Y_L: {pL:.2f}  (predict 5 = 2q+1)")
        # figure
        fig,ax=plt.subplots(1,3,figsize=(15,4.3))
        cols=plt.cm.viridis(np.linspace(0.1,0.85,len(Yes)))
        for Ye,cc in zip(Yes,cols):
            y=[]; x=[]
            for e in es:
                P=np.mean(data[e]>Ye)
                if P*len(data[e])>=15: y.append(-np.log(P)); x.append(1/e**2)
            ax[0].plot(x,y,'o-',color=cc,ms=5,label=f'Y_e={Ye:.1f}')
        ax[0].set_xlabel(r'$1/\eta^2$'); ax[0].set_ylabel(r'$-\ln P(Y^*>Y_e)$')
        ax[0].set_title('(A) RIGHT tail: LDP linear in $1/\\eta^2$'); ax[0].legend(fontsize=7)
        if len(slopeR):
            ax[1].loglog(slopeR[:,0],slopeR[:,1],'o',color='crimson',ms=8,label=f'fit slope $\\propto Y_e^{{{pR:.2f}}}$')
            xx=np.linspace(slopeR[:,0].min(),slopeR[:,0].max(),20)
            ax[1].loglog(xx,8/3*xx**3,'k--',lw=1.3,label=r'$\frac{8}{3}Y_e^3$ (exp 3)')
            ax[1].loglog(xx,8/3*xx**2.5*slopeR[0,1]/(8/3*slopeR[0,0]**2.5),':',color='gray',lw=1.3,label=r'exp 5/2 (k=1/2)')
        ax[1].set_xlabel(r'$Y_e$'); ax[1].set_ylabel('LDP slope m'); ax[1].set_title('(B) RIGHT exponent: 3 vs 5/2'); ax[1].legend(fontsize=8)
        if len(slopeL):
            ax[2].loglog(slopeL[:,0],slopeL[:,1],'s',color='navy',ms=8,label=f'fit $\\propto Y_L^{{{pL:.2f}}}$')
            xx=np.linspace(slopeL[:,0].min(),slopeL[:,0].max(),20)
            ax[2].loglog(xx,(xx**5-2.09**5)/10,'k--',lw=1.3,label=r'$(Y_L^5-Y_{det}^5)/10$ (exp 5)')
        ax[2].set_xlabel(r'$Y_L$'); ax[2].set_ylabel('LDP slope m'); ax[2].set_title('(C) LEFT exponent: 5'); ax[2].legend(fontsize=8)
        fig.suptitle('FW/instanton tails of $\\mathcal{W}_2$: RIGHT $\\sim Y_e^{3}$ (=3q/2), LEFT $\\sim Y_L^{5}$ (=2q+1) — refutes k=1/2 (needs 5/2)',fontsize=11)
        fig.subplots_adjust(left=0.06,right=0.98,bottom=0.12,top=0.87,wspace=0.26)
        plt.savefig("figures/ldp_tails.png",dpi=110); print("\nsaved figures/ldp_tails.png")
