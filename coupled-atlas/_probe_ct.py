import numpy as np
from weaknoise_highorder import make_Yfuncs
import weaknoise_highorder as W
from _probe_div import run_full  # reuse the full-Y sampler
# extend run_full to also return <s0^2>: quick reimplement returning s0 too
from math import comb
from numpy.polynomial.hermite_e import hermeval
from scipy.ndimage import gaussian_filter1d
from weaknoise_greens import build_backbone, node, V

def run_s(delta,MAX,YF,N,h=1e-3,seed=0,Y0=4.0,Yend=-2.6,batch=4000):
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
    rng=np.random.default_rng(seed); allY=[]; s0sq=0.0; Nsum=0
    while Nsum<N:
        nb=min(batch,N-Nsum); Nsum+=nb
        dW=np.sqrt(h)*rng.standard_normal((nb,nY+1))
        xi0=gaussian_filter1d(dW,sigma=sigpix,axis=1,mode='constant')/h
        s=[ (dW[:,lo:hi]*Kj[j][None,:]).sum(1) for j in range(Jmax+1) ]
        s0sq+=np.sum(s[0]**2)
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
        allY.append(np.array([f(*argl) for f in funcs]))
    Y=np.concatenate(allY,axis=1)
    return Y, s0sq/Nsum

MAX=7; YF=make_Yfuncs(MAX)
v0=0.1339
print("Verify counterterm: predicted Var(Y4)_div = <s0^2>*(5/3)v0^3, 2Cov(Y1,Y7)_div = <s0^2>*6 v0^3")
print(f"  (5/3)v0^3={5/3*v0**3:.5f}   6 v0^3={6*v0**3:.5f}   total (23/3)v0^3={23/3*v0**3:.5f}")
print(f"  {'delta':>6} {'<s0^2>':>9} {'Var(Y4)':>9} {'2Cov17':>9} | {'VarY4/<s2>':>10} {'2Cv17/<s2>':>10}")
rows=[]
for delta in [0.10,0.07,0.05]:
    Ys=[run_s(delta,MAX,YF,N=20000,seed=s)[0] for s in range(3)]
    _,s0sq=run_s(delta,MAX,YF,N=20000,seed=99)
    Y=np.concatenate(Ys,axis=1)
    VarY4=Y[3].var(); Cov17=np.cov(Y[0],Y[6])[0,1]
    rows.append((delta,s0sq,VarY4,2*Cov17))
    print(f"  {delta:>6.3f} {s0sq:>9.3f} {VarY4:>9.4f} {2*Cov17:>9.4f} | {VarY4/s0sq:>10.5f} {2*Cov17/s0sq:>10.5f}")
print("\n  divergence SLOPE via finite differences (d/d<s0^2>):")
for a in range(len(rows)-1):
    d1,s1,V1,C1=rows[a]; d2,s2,V2,C2=rows[a+1]
    print(f"   delta {d1}->{d2}: dVarY4/d<s2>={ (V2-V1)/(s2-s1):.5f} (pred {5/3*v0**3:.5f})  d2Cov17/d<s2>={ (C2-C1)/(s2-s1):.5f} (pred {6*v0**3:.5f})")
