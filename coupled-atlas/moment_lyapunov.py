"""
Stage 1 (Route C, decisive go/no-go) — moment-hierarchy Lyapunov check for the last lemma.

Inner SDE x=dp:  dx = -(2 pbar x + x^2) dtau + eta dW.
Moment ODEs:  d/dtau E[x^k] = -2 pbar k E[x^k] - k E[x^{k+1}] + (eta^2/2) k(k-1) E[x^{k-2}].
Integrate along the REAL swept canard pbar(tau) (from canard(), not frozen), for several Delta (rho).
Lyapunov:  H_m = sum_{k=2}^{2m} (beta^k/k!) E[x^k].  Drift defect (the test, want <=0):
   Dd_m = [drift part of dH_m] + 2 c0 H_m
        = -2 pbar * sum_{k=2}^{2m} beta^k/(k-1)! E[x^k]   (restoring, rate 2 pbar >= 2 c0)
          - sum_{k=2}^{2m} beta^k/(k-1)! E[x^{k+1}]        (cubic coupling = the obstruction)
          + 2 c0 H_m .
The noise source (eta^2 beta^2/2)(1+H_{m-1}) ~ O(eta^2) is added separately.
GO if Dd_m <= O(eta^2) pointwise (incl. turning) with a single Delta-independent beta.
numpy only.
"""
import numpy as np
from math import factorial

def Vfun(Y,D): aY=np.abs(Y); return np.sign(Y)*aY*(aY+D)
def canard(Yg,D):
    n=len(Yg); pb=np.empty(n); pb[0]=np.sqrt(max(Vfun(Yg[0],D),1e-12))
    for i in range(n-1):
        Y=Yg[i]; h=Yg[i+1]-Yg[i]; f=lambda Yv,pv: pv*pv-Vfun(Yv,D)
        k1=f(Y,pb[i]);k2=f(Y+.5*h,pb[i]+.5*h*k1);k3=f(Y+.5*h,pb[i]+.5*h*k2);k4=f(Y+h,pb[i]+h*k3)
        nx=pb[i]+(h/6)*(k1+2*k2+2*k3+k4); pb[i+1]=nx if (np.isfinite(nx) and 0<nx<1e6) else pb[i]
    return pb

K=14; karr=np.arange(K+1)
def dM(M,pb,eta):
    Mp1=np.concatenate([M[1:],[0.0]])          # M[k+1]
    Mm2=np.concatenate([[0.0,0.0],M[:-2]])      # M[k-2]
    d=-2*pb*karr*M - karr*Mp1 + (eta**2/2)*karr*(karr-1)*Mm2
    d[0]=0.0; return d

def run(D,eta,Y0=3.0,Yend=0.04,dt=1e-3):
    Yg=np.arange(Y0,Yend,-dt); pb=canard(Yg,D); n=len(Yg)
    v0=eta**2/(4*pb[0]); M=np.zeros(K+1); M[0]=1.0
    for j in range(1,K//2+1): M[2*j]=factorial(2*j)/(2**j*factorial(j))*v0**j  # Gaussian (2j-1)!! v0^j
    rec=[]  # (Y, pbar, M.copy())
    for i in range(n-1):
        h=dt
        k1=dM(M,pb[i],eta);k2=dM(M+.5*h*k1,pb[i],eta);k3=dM(M+.5*h*k2,pb[i],eta);k4=dM(M+h*k3,pb[i+1],eta)
        M=M+(h/6)*(k1+2*k2+2*k3+k4)
        rec.append((Yg[i+1],pb[i+1],M.copy()))
    return rec

def defect(rec,beta,m,c0,eta):
    w=np.array([beta**k/factorial(k) for k in range(2*m+2)])  # weights up to 2m+1
    maxDd=-1e9; argY=None; Hmax=0
    for (Y,pb,M) in rec:
        H = sum(w[k]*M[k] for k in range(2,2*m+1))
        A = -2*pb*sum((beta**k/factorial(k-1))*M[k] for k in range(2,2*m+1))
        B = -sum((beta**k/factorial(k-1))*M[k+1] for k in range(2,2*m+1))
        Dd = A+B+2*c0*H
        if Dd>maxDd: maxDd=Dd; argY=Y
        Hmax=max(Hmax,abs(H))
    return maxDd, argY, Hmax

if __name__=="__main__":
    eta=0.3; c0=0.69
    Ds=[0.05,0.2,0.5,1.0,2.5]; betas=[0.5,1.0,2.0,3.0]
    print(f"eta={eta}  c0={c0}  noise source ~ eta^2 beta^2/2 = {eta**2/2:.4f}*beta^2")
    print("max_tau Dd_m  (want <= O(eta^2); GO if <=0 with a single Delta-indep beta)\n")
    for m in [2,3,5]:
        print(f"--- m={m} (moments up to x^{2*m}) ---")
        print(f"{'beta':>6}"+"".join(f"  D={D:<5}" for D in Ds)+"   where(worst Y)")
        for beta in betas:
            row=[]; worstY=[]
            for D in Ds:
                rec=run(D,eta); md,aY,Hm=defect(rec,beta,m,c0,eta); row.append(md); worstY.append(aY)
            flag="GO" if max(row)<=2*eta**2*beta**2 else ("ok" if max(row)<=0.05 else "FAIL")
            print(f"{beta:>6.1f}"+"".join(f"{r:>9.4f}" for r in row)+f"   {flag} (Y~{np.mean(worstY):+.2f})")
        print()
    # ---- figure: Dd(Y) along the sweep (incl. turning) + max Dd vs rho ----
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    def traj(D,beta,m):
        rec=run(D,eta); Ys=[];Dd=[]
        for (Y,pb,M) in rec:
            H=sum(beta**k/factorial(k)*M[k] for k in range(2,2*m+1))
            A=-2*pb*sum(beta**k/factorial(k-1)*M[k] for k in range(2,2*m+1))
            B=-sum(beta**k/factorial(k-1)*M[k+1] for k in range(2,2*m+1))
            Ys.append(Y); Dd.append(A+B+2*c0*H)
        return np.array(Ys),np.array(Dd)
    fig,ax=plt.subplots(1,2,figsize=(11,4.3))
    for D,c in [(0.05,'crimson'),(0.5,'teal'),(2.5,'navy')]:
        for m,ls in [(2,'-'),(5,'--')]:
            Ys,Dd=traj(D,1.0,m); ax[0].plot(Ys,Dd,ls,color=c,lw=1.5,label=f'$\\rho$~{D/eta**(2/3):.2f}, m={m}')
    ax[0].axhline(0,color='k',lw=0.6); ax[0].axvline(0,color='gray',ls=':',lw=0.8)
    ax[0].set_xlabel('Y (turning at 0)'); ax[0].set_ylabel('drift defect  Hdot_m + 2 c0 H_m')
    ax[0].set_title('(A) defect < 0 all Y incl. turning (beta=1; m=2 solid, 5 dashed)'); ax[0].legend(fontsize=7); ax[0].set_xlim(2.0,-0.05)
    rhos=[];mx={b:[] for b in betas}
    for D in Ds:
        rhos.append(D/eta**(2/3))
        for beta in betas:
            rec=run(D,eta); md,_,_=defect(rec,beta,5,c0,eta); mx[beta].append(md)
    for beta in betas: ax[1].plot(rhos,mx[beta],'o-',label=f'$\\beta$={beta}')
    ax[1].axhline(0,color='k',lw=0.6); ax[1].set_xscale('log')
    ax[1].set_xlabel(r'$\rho=\Delta/\ell$'); ax[1].set_ylabel('max over tau of defect (m=5)')
    ax[1].set_title('(B) GO: defect < 0 uniformly in rho (single beta)'); ax[1].legend(fontsize=8)
    fig.suptitle('Route C Lyapunov go/no-go: drift defect < 0 through the turning, all m, uniform in rho => GO',fontsize=11)
    fig.subplots_adjust(left=0.08,right=0.98,bottom=0.12,top=0.88,wspace=0.26)
    plt.savefig("figures/moment_lyapunov.png",dpi=115); print("saved figures/moment_lyapunov.png")
