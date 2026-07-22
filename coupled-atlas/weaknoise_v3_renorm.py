"""
Program 2 -- the RENORMALIZED v3^ren via variance-reduced sampling.

v3 = Var(Y4) + 2Cov(Y3,Y5) + 2Cov(Y2,Y6) + 2Cov(Y1,Y7). The only divergence is the boundary self-contraction
<xi(Y*_0)^2> from (i) Var(Y4) via B4=Y1^3/3, (ii) Cov(Y1,Y7) via R7=Y1^5/5 (both machine-verified, notes sec 8).
Because s0=xi_delta(node) and Y1 are JOINTLY GAUSSIAN, the divergent moments are known in closed form:
    <s0^2 B4^2> = (5/3)<s0^2> v0^3 + 10 c^2 v0^2 ,   <s0^2 Y1 R7> = 3<s0^2> v0^3 + 18 c^2 v0^2 ,  c=<s0 Y1>.
RENORMALIZE = drop the <s0^2>-divergent pieces; keep the finite companions (10 c^2 v0^2, 18 c^2 v0^2). This lets us
sample only the FINITE-variance regular parts:
    Var(Y4)^ren   = <A4^2> + 2<A4 s0 B4> + 10 c^2 v0^2 - <Y4>^2 ,   A4 = Y4 - s0 B4
    Cov(Y1,Y7)^ren= <Y1 (Y7 - s0^2 R7)> + 18 c^2 v0^2
    v3^ren = Var(Y4)^ren + 2<Y3 Y5> + 2(<Y2 Y6>-<Y2><Y6>) + 2 Cov(Y1,Y7)^ren
Then extrapolate delta->0. (Caveat: Y5,Y6 still carry LINEAR s1,s2 whose high delta-derivatives inflate the
variance of Cov(Y3,Y5),Cov(Y2,Y6) -- so precision is bounded; we report honest error + delta-trend.)
"""
import numpy as np
from weaknoise_highorder import make_Yfuncs
import weaknoise_highorder as W
from math import comb
from numpy.polynomial.hermite_e import hermeval
from scipy.ndimage import gaussian_filter1d
from weaknoise_greens import build_backbone, node, V

def sample(delta,MAX,YF,N,h=1e-3,seed=0,Y0=4.0,Yend=-2.6,batch=4000):
    funcs,_,M=YF
    Yg,u0,u0p,ps,psp=build_backbone(Y0=Y0,Yend=Yend,h=h)
    nd=node(Yg,u0,u0p,ps,psp)
    Ystar=nd['Ystar']; ist=nd['istar']; frac=nd['f']; u0ps=nd['u0p_star']
    Vs=[W.Vd(Ystar,d) for d in range(MAX+1)]
    nY=len(Yg)-1; sigpix=delta/h
    half=int(np.ceil(6*delta/h)); iL=ist-1
    lo=max(0,iL-half-2); hi=min(nY+1,iL+half+3)
    dg=Ystar-Yg[lo:hi]; z=dg/delta
    phi0=np.exp(-z*z/2)/(np.sqrt(2*np.pi)*delta)
    Jmax=MAX-2; Kj=[]
    for j in range(Jmax+1):
        c=np.zeros(j+1); c[j]=1.0
        Kj.append(((-1)**j)*delta**(-j)*hermeval(z,c)*phi0)
    rng=np.random.default_rng(seed); Ys=[]; s0=[]
    Nsum=0
    while Nsum<N:
        nb=min(batch,N-Nsum); Nsum+=nb
        dW=np.sqrt(h)*rng.standard_normal((nb,nY+1))
        xi0=gaussian_filter1d(dW,sigma=sigpix,axis=1,mode='constant')/h
        s=[ (dW[:,lo:hi]*Kj[j][None,:]).sum(1) for j in range(Jmax+1) ]
        uu=np.zeros((MAX,nb)); upp=np.zeros((MAX,nb)); snap={}; dY=-h
        for i in range(nY):
            Vi=V(Yg[i]); xi_i=xi0[:,i]
            src=np.empty((MAX,nb)); src[0]=u0[i]; src[1:]=uu[:MAX-1]
            newup=upp+Vi*uu*dY+xi_i*src*dY
            uu=uu+upp*dY; upp=newup
            if i in (ist-1,ist): snap[i]=(uu.copy(),upp.copy())
        lo_s=snap[ist-1]; hi_s=snap[ist]
        Uval=lo_s[0]+frac*(hi_s[0]-lo_s[0]); Uder=lo_s[1]+frac*(hi_s[1]-lo_s[1])
        Umat=np.zeros((MAX+1,MAX+1,nb)); Umat[0,1]=u0ps; Umat[1:,0]=Uval; Umat[1:,1]=Uder
        for k in range(MAX+1):
            for m in range(2,MAX+1):
                p=m-2; tot=np.zeros(nb)
                for j in range(p+1):
                    cc=comb(p,j); tot=tot+cc*Vs[j]*Umat[k,p-j]
                    if k>=1 and j<=Jmax: tot=tot+cc*s[j]*Umat[k-1,p-j]
                Umat[k,m]=tot
        Ufull=np.zeros((MAX+1,M,nb)); Ufull[:, :MAX+1, :]=Umat
        argl=[Ufull[k,m] for k in range(MAX+1) for m in range(M)]
        Ys.append(np.array([f(*argl) for f in funcs])); s0.append(s[0])
    return np.concatenate(Ys,axis=1), np.concatenate(s0)

def v3ren_from(Y,s0):
    Y1,Y2,Y3,Y4,Y5,Y6,Y7=Y
    v0=np.mean(Y1**2); c=np.mean(s0*Y1)
    B4=Y1**3/3; R7=Y1**5/5
    A4=Y4-s0*B4
    VarY4_ren=np.mean(A4**2)+2*np.mean(A4*s0*B4)+10*c**2*v0**2-np.mean(Y4)**2
    Y7mod=Y7-s0**2*R7
    Cov17_ren=np.mean(Y1*Y7mod)+18*c**2*v0**2
    Cov35=np.mean(Y3*Y5)                       # <Y3>=0
    Cov26=np.mean(Y2*Y6)-np.mean(Y2)*np.mean(Y6)
    v3=VarY4_ren+2*Cov35+2*Cov26+2*Cov17_ren
    return dict(v0=v0,c=c,VarY4_ren=VarY4_ren,Cov35=Cov35,Cov26=Cov26,Cov17_ren=Cov17_ren,v3=v3)

if __name__=="__main__":
    MAX=7; YF=make_Yfuncs(MAX)
    print("=== v3^ren via renormalized (variance-reduced) sampling ===")
    print(f"  {'delta':>6} {'v0':>7} {'VarY4^r':>9} {'2Cov35':>9} {'2Cov26':>9} {'2Cov17^r':>9} {'v3^ren':>9} {'+-':>7}")
    res={}
    for delta in [0.09,0.07,0.055,0.045]:
        seeds=[sample(delta,MAX,YF,N=30000,seed=s) for s in range(4)]
        vals=[v3ren_from(Y,s0) for Y,s0 in seeds]
        def M(k): return np.mean([v[k] for v in vals])
        def SE(k): return np.std([v[k] for v in vals])/2
        res[delta]=(M('v3'),SE('v3'))
        print(f"  {delta:>6.3f} {M('v0'):>7.4f} {M('VarY4_ren'):>9.4f} {2*M('Cov35'):>9.4f} {2*M('Cov26'):>9.4f} "
              f"{2*M('Cov17_ren'):>9.4f} {M('v3'):>9.4f} {SE('v3'):>7.4f}",flush=True)
    ds=np.array(sorted(res)); vs=np.array([res[d][0] for d in ds])
    A=np.vstack([np.ones_like(ds),ds,ds**2]).T; cq,*_=np.linalg.lstsq(A,vs,rcond=None)
    Al=np.vstack([np.ones_like(ds),ds]).T; cl,*_=np.linalg.lstsq(Al,vs,rcond=None)
    print(f"\n  delta->0 extrapolation:  linear={cl[0]:+.4f}   quadratic={cq[0]:+.4f}")
    print(f"  => v3^ren ~ {0.5*(cl[0]+cq[0]):+.3f}   (compare v2=0.11; ratio v3/v2={0.5*(cl[0]+cq[0])/0.11:.2f})")
