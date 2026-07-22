"""
Program 2, Route 2b — CONVERGED weak-noise cumulant SERIES off the exact Weber backbone, to test the
resurgent-growth signature WITHOUT the FP boundary-layer contamination.

Method: sample white noise xi; integrate the perturbation hierarchy u_k'' - V u_k = xi u_{k-1} (k=1..K)
with recessive BC. For a ladder of small eta, form the truncated field U(Y)=sum_k eta^k u_k(Y) and find its
first node Y*(eta) per realization (in a window around the base node Y*_0 -- no boundary layer, no rare escape).
Average -> cumulants(eta); fit the eta-series of mean, var/eta^2, skew/eta.  Compare coefficient growth to the
(contaminated) FP extraction  Var=e2*(0.200-0.245 e2+0.572 e2^2), ratios 1.23->2.33.
numpy only, vectorized over realizations.
"""
import numpy as np

def V(Y): return np.sign(Y)*Y*Y

def fields(N, Y0=4.0, Yend=-5.0, h=1.0e-3, K=3, seed=0, Wwin=120):
    rng=np.random.default_rng(seed)
    nY=int(round((Y0-Yend)/h)); Yg=Y0-np.arange(nY+1)*h; dY=-h; sh=np.sqrt(h)
    # deterministic backbone
    u0=np.empty(nY+1); u0p=np.empty(nY+1); u0[0]=1.0; u0p[0]=-Y0
    for i in range(nY):
        u0[i+1]=u0[i]+u0p[i]*dY; u0p[i+1]=u0p[i]+V(Yg[i])*u0[i]*dY
    istar=next(i for i in range(1,nY+1) if Yg[i]<0 and u0[i-1]*u0[i]<0)
    lo,hi=istar-Wwin,istar+Wwin
    # perturbation fields u_1..u_K, store only window slices
    u=[np.zeros(N) for _ in range(K+1)]; up=[np.zeros(N) for _ in range(K+1)]
    store=[np.zeros((hi-lo,N)) for _ in range(K+1)]  # k=0..K
    for i in range(nY):
        Vi=V(Yg[i]); dB=sh*rng.standard_normal(N)
        newp=[None]*(K+1)
        for k in range(1,K+1):
            src=u0[i] if k==1 else u[k-1]
            newp[k]=up[k]+Vi*u[k]*dY+src*dB
        for k in range(1,K+1):
            u[k]=u[k]+up[k]*dY
        for k in range(1,K+1):
            up[k]=newp[k]
        if lo<=i<hi:
            store[0][i-lo]=u0[i]
            for k in range(1,K+1): store[k][i-lo]=u[k]
    Ywin=Yg[lo:hi]
    return Ywin,store,K

def cumulants(Ywin,store,K,eta):
    hi=store[0].shape[0]; N=store[0].shape[1]
    U=np.zeros((hi,N))
    for k in range(K+1): U+=eta**k*store[k]
    # find first sign change down the window (Ywin descending)
    sgn=np.sign(U)
    chg=(sgn[:-1]*sgn[1:]<0)
    Ystar=np.full(N,np.nan)
    # first True per column
    idx=np.where(chg.any(0), chg.argmax(0), -1)
    good=idx>=0
    j=idx[good]; cols=np.where(good)[0]
    Ulo=U[j,cols]; Uhi=U[j+1,cols]; Ylo=Ywin[j]; Yhi=Ywin[j+1]
    f=Ulo/(Ulo-Uhi)
    Ystar[cols]=Ylo+f*(Yhi-Ylo)
    Ys=Ystar[np.isfinite(Ystar)]
    m=Ys.mean(); sd=Ys.std(); sk=np.mean(((Ys-m)/sd)**3)
    return m,sd*sd,sk,len(Ys)

if __name__=="__main__":
    K=3; h=1e-3; N=200000
    Ywin,store,K=fields(N,h=h,K=K,seed=1)
    etas=np.array([0.06,0.09,0.12,0.16,0.20,0.25])
    rows=[]
    print(f"=== converged weak-noise cumulants (K={K} orders, h={h}, N={N}) ===")
    print(f"  {'eta':>6} {'eta^2':>7} {'mean':>9} {'var':>10} {'var/e2':>8} {'skew':>8} {'skew/eta':>8}")
    for eta in etas:
        m,v,sk,n=cumulants(Ywin,store,K,eta)
        rows.append((eta,m,v,sk)); e2=eta*eta
        print(f"  {eta:6.3f} {e2:7.4f} {m:9.4f} {v:10.6f} {v/e2:8.4f} {sk:8.4f} {sk/eta:8.4f}")
    R=np.array(rows); et=R[:,0]; e2=et*et
    cm=np.polyfit(e2,R[:,1],2)[::-1]           # mean = m0 + m1 e2 + m2 e2^2
    cv=np.polyfit(e2,R[:,2]/e2,2)[::-1]        # var/e2 = v0 + v1 e2 + v2 e2^2
    cs=np.polyfit(e2,R[:,3]/et,1)[::-1]        # skew/eta = s0 + s1 e2
    print("\n  CONVERGED perturbative series (eta^2=4/beta):")
    print(f"    <Y*> = {cm[0]:+.3f} {cm[1]:+.3f} e2 {cm[2]:+.3f} e2^2     (backbone -2.188, m1~0.212)")
    print(f"    Var  = e2*( {cv[0]:.4f} {cv[1]:+.4f} e2 {cv[2]:+.4f} e2^2 ) (C_V=0.1339 closed-form)")
    print(f"    skew = eta*( {cs[0]:+.3f} {cs[1]:+.3f} e2 )                (s0~1.18 MC)")
    print(f"\n  coefficient ratios:  var |v1/v0|={abs(cv[1]/cv[0]):.2f}, |v2/v1|={abs(cv[2]/cv[1]):.2f}")
    print(f"  [contaminated FP gave v0=0.200, |v1/v0|=1.23, |v2/v1|=2.33]")
