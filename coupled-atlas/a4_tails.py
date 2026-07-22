"""
Program 1, Tier A / step A4 — numerical audit of the two-sided tail bounds that give moment-determinacy.
(See regime-tests/PROGRAM1_TIERA_THEOREM.md.)

Three checks, all cheap and decisive:
 (1) integrated-Riccati identity  eta W_s = q_s - q0 + s^3/3 + int q^2  (the identity the left-tail proof rests on)
 (2) branch-2 pure-Riccati survival from -L over unit time  <=  2 exp(-(L-1)^2/2 eta^2)   [rigorous bound]
 (3) FP left tail -logP(Y*<-s): lies ABOVE the rigorous rate s^5/(72 eta^2), exponent ~5.
numpy + fp_cusp only.
"""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from fp_cusp import solve_fp

# ---------- (1) integrated-Riccati identity ----------
def ident(eta=1.4,N=100000,dt=2e-4,s=3.0,q0=0.5,seed=3):
    rng=np.random.default_rng(seed); n=int(round(s/dt)); sq=np.sqrt(dt)
    q=np.full(N,q0); W=np.zeros(N); Iq2=np.zeros(N); alive=np.ones(N,bool)
    for i in range(n):
        t=i*dt; V=-t*t; dW=sq*rng.standard_normal(N)
        Iq2[alive]+=q[alive]**2*dt; W[alive]+=dW[alive]
        q[alive]+=(V-q[alive]**2)*dt+eta*dW[alive]; alive&=(q>-40)
    lhs=eta*W[alive]; rhs=q[alive]-q0+s**3/3+Iq2[alive]
    lb=q[alive]-q0+s**3/3
    return alive.sum(),np.abs(lhs-rhs).max(),bool(np.all(lhs>=lb-1e-6)),W[alive]

# ---------- (2) branch-2 survival from deep -L ----------
def surv_from_deep(L,eta=1.0,N=60000,dt=4e-4,T=1.0,seed=1):
    rng=np.random.default_rng(seed); n=int(T/dt); sq=np.sqrt(dt)
    q=np.full(N,-L); alive=np.ones(N,bool)
    for i in range(n):
        q[alive]+=(-q[alive]**2)*dt+eta*sq*rng.standard_normal(alive.sum())
        alive&=(q>-1e6)&np.isfinite(q)
    return alive.mean()

# ---------- (3) FP left tail vs rigorous rate ----------
def fp_left(beta,ss):
    eta=2.0/np.sqrt(beta)
    y,F=solve_fp(eta=eta,Y0=3.0,pmin=-13.0,pmax=8.0,dp=0.025,dt=8e-5,tau_max=11.0)
    F=np.clip(F,1e-300,1.0); mlp=[]
    for s in ss:
        i=np.searchsorted(y,-s); mlp.append(-np.log(max(F[i] if 0<i<len(y) else 1e-300,1e-300)))
    return eta,np.array(mlp)

if __name__=="__main__":
    print("=== (1) integrated-Riccati identity (survivors) ===")
    nsv,err,ok,Wsv=ident()
    print(f"  survivors={nsv}, max|err|={err:.1e} (Euler dt), branch-1 lower bound holds: {ok}")
    print(f"  survivor W_s: min={Wsv.min():.2f} mean={Wsv.mean():.2f}  (pushed positive => Gaussian cost)")

    print("\n=== (2) pure-Riccati survival from -L (unit time) vs rigorous bound 2 exp(-(L-1)^2/2eta^2) ===")
    eta=1.0
    for L in [2.,3.,4.,5.]:
        f=surv_from_deep(L,eta); b=2*np.exp(-(L-1)**2/(2*eta**2))
        print(f"  L={L:.0f}: surv={f:.2e}  bound={b:.2e}  {'OK' if f<=b else 'VIOLATION'}")

    print("\n=== (3) FP left tail: bound s^5/(72 eta^2) holds, exponent ~5 ===")
    ss=np.linspace(3.0,7.5,20)
    fig,ax=plt.subplots(1,2,figsize=(11,4.2))
    for beta,col in [(2.0,"#264653"),(4.0,"#e76f51")]:
        eta,mlp=fp_left(beta,ss); rig=ss**5/(72*eta**2); sharp=ss**5/(10*eta**2)
        m=(mlp>5)&(mlp<600); p=np.polyfit(np.log(ss[m]),np.log(mlp[m]),1)[0]
        holds=np.all(mlp[m]>=rig[m])
        print(f"  beta={beta}: clean-window exponent={p:.2f} (pred 5); bound holds: {holds}")
        ax[0].plot(ss,mlp,"o-",color=col,ms=3,label=f"beta={beta}: -logP(FP)")
        ax[0].plot(ss,rig,"--",color=col,lw=1,label=f"beta={beta}: rigorous s^5/72eta^2")
        ax[1].loglog(ss[m],mlp[m],"o-",color=col,ms=3,label=f"beta={beta} (exp {p:.2f})")
        ax[1].loglog(ss[m],sharp[m],":",color=col,lw=1)
    ax[0].set_xlabel("s"); ax[0].set_ylabel("-log P(Y*<-s)"); ax[0].set_ylim(0,650)
    ax[0].set_title("(a) left tail vs rigorous rate (bound below data)"); ax[0].legend(fontsize=7)
    ax[1].set_xlabel("s"); ax[1].set_ylabel("-log P"); ax[1].set_title("(b) log-log: exponent ~5 (dotted = sharp s^5/10eta^2)")
    ax[1].legend(fontsize=7)
    fig.suptitle("A4 left tail: rigorous exp-5 bound (const 1/72) below FP; sharp instanton 1/10 (dotted)",fontsize=10)
    fig.subplots_adjust(left=0.08,right=0.98,top=0.9,bottom=0.12,wspace=0.25)
    fig.savefig("figures/a4_tails.png",dpi=130); print("\nsaved figures/a4_tails.png")
