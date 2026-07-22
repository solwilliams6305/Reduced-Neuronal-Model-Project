"""
Program 2, Phase 2 (large-order) -- extend the regularized node-shift recursion to arbitrary order to get
v3, v4 and read off the LARGE-ORDER growth of the weak-noise variance coefficients v_n. This is the decisive
test of the resurgent-trans-series hypothesis: factorial growth v_n ~ n!/A^n (asymptotic, genuine trans-series,
Borel singularity at A) vs finite radius v_n ~ C/R^n (convergent). Var = sum_m v_{m-1} eta^{2m},
        v_{m-1} = Var(Y_m) + 2 sum_{k=1}^{m-1} Cov(Y_k, Y_{2m-k}).
So v_{m-1} needs Y_1..Y_{2m-1}: v1<-Y3, v2<-Y5, v3<-Y7, v4<-Y9.

Node-shift functionals generated symbolically (cached). Fields marched with mollified noise (width delta);
U[k][m]=u_k^{(m)}(Y*_0) built from the field ODE u_k''=V u_k + xi_delta u_{k-1} (Leibniz), with xi_delta^{(j)}(Y*_0)
from Hermite-weighted local sums. delta->0 extrapolation gives the boundary-Wick values. Precision degrades with
order (high xi-derivatives ~ delta^{-(j+1/2)}); we report per-delta values + extrapolation + honest error.
"""
import numpy as np, sympy as sp, os, pickle, time
from math import comb
from numpy.polynomial.hermite_e import hermeval
from scipy.ndimage import gaussian_filter1d
from weaknoise_greens import build_backbone, node, V

def make_Yfuncs(MAX, cache=True):
    M=MAX+1
    U=[[sp.symbols(f'U{k}_{m}') for m in range(M)] for k in range(MAX+1)]
    Y=[0]+[sp.symbols(f'Y{n}') for n in range(1,MAX+1)]
    syms=[U[k][m] for k in range(MAX+1) for m in range(M)]
    fn=f"_yexprs_{MAX}.txt"
    if cache and os.path.exists(fn):
        with open(fn) as f: exprs=[sp.sympify(line.strip()) for line in f if line.strip()]
    else:
        eta=sp.symbols('eta')
        D=sum(eta**k*Y[k] for k in range(1,MAX+1))
        expr=sp.expand(sum(eta**k*U[k][m]*D**m/sp.factorial(m) for k in range(MAX+1) for m in range(M))).subs(U[0][0],0)
        sol={}; exprs=[]
        for n in range(1,MAX+1):
            c=sp.expand(sp.expand(expr).coeff(eta,n).subs(sol))
            sol[Y[n]]=sp.expand(sp.solve(c,Y[n])[0]); exprs.append(sol[Y[n]])
        if cache:
            with open(fn,'w') as f:
                for e in exprs: f.write(sp.srepr(e)+"\n")
    funcs=[sp.lambdify(syms,e.subs(U[0][0],0),'numpy') for e in exprs]
    return (funcs,MAX,M)

def Vd(Y,d): return {0:-Y*Y,1:-2*Y,2:-2.0}.get(d,0.0)

def run(delta, MAX, YF, N=12000, h=1e-3, seed=0, Y0=4.0, Yend=-2.6, batch=4000):
    funcs,_,M=YF
    Yg,u0,u0p,ps,psp=build_backbone(Y0=Y0,Yend=Yend,h=h)
    nd=node(Yg,u0,u0p,ps,psp)
    Ystar=nd['Ystar']; ist=nd['istar']; frac=nd['f']; u0ps=nd['u0p_star']
    Vs=[Vd(Ystar,d) for d in range(MAX+1)]
    nY=len(Yg)-1; sigpix=delta/h
    half=int(np.ceil(6*delta/h)); iL=ist-1
    lo=max(0,iL-half-2); hi=min(nY+1,iL+half+3)
    dg=Ystar-Yg[lo:hi]; z=dg/delta
    phi0=np.exp(-z*z/2)/(np.sqrt(2*np.pi)*delta)
    # xi^{(j)}(Y*_0) kernels: phi^{(j)}(x)=(-1)^j delta^{-j} He_j(x/delta) phi0(x)
    Jmax=MAX-2
    Kj=[]
    for j in range(Jmax+1):
        c=np.zeros(j+1); c[j]=1.0
        Kj.append(((-1)**j)*delta**(-j)*hermeval(z,c)*phi0)
    rng=np.random.default_rng(seed)
    momsum={}; Nsum=0
    while Nsum<N:
        nb=min(batch,N-Nsum); Nsum+=nb
        dW=np.sqrt(h)*rng.standard_normal((nb,nY+1))
        xi0=gaussian_filter1d(dW,sigma=sigpix,axis=1,mode='constant')/h
        s=[ (dW[:,lo:hi]*Kj[j][None,:]).sum(1) for j in range(Jmax+1) ]
        uu=np.zeros((MAX,nb)); upp=np.zeros((MAX,nb)); snap={}
        dY=-h
        for i in range(nY):
            Vi=V(Yg[i]); xi_i=xi0[:,i]
            src=np.empty((MAX,nb)); src[0]=u0[i]; src[1:]=uu[:MAX-1]
            newup=upp+Vi*uu*dY+xi_i*src*dY
            uu=uu+upp*dY; upp=newup
            if i in (ist-1,ist): snap[i]=(uu.copy(),upp.copy())
        lo_s=snap[ist-1]; hi_s=snap[ist]
        Uval=lo_s[0]+frac*(hi_s[0]-lo_s[0]); Uder=lo_s[1]+frac*(hi_s[1]-lo_s[1])
        Umat=np.zeros((MAX+1,MAX+1,nb))
        Umat[0,1]=u0ps; Umat[1:,0]=Uval; Umat[1:,1]=Uder
        for k in range(MAX+1):
            for m in range(2,MAX+1):
                p=m-2; tot=np.zeros(nb)
                for j in range(p+1):
                    cc=comb(p,j); tot=tot+cc*Vs[j]*Umat[k,p-j]
                    if k>=1 and j<=Jmax: tot=tot+cc*s[j]*Umat[k-1,p-j]
                Umat[k,m]=tot
        Ufull=np.zeros((MAX+1,M,nb)); Ufull[:, :MAX+1, :]=Umat
        argl=[Ufull[k,m] for k in range(MAX+1) for m in range(M)]
        Ys=[f(*argl) for f in funcs]   # Y1..Y_MAX
        # accumulate all needed products: means, and cross moments Y_k Y_l
        for k in range(1,MAX+1):
            momsum[('m',k)]=momsum.get(('m',k),0)+np.sum(Ys[k-1])
            for l in range(k,MAX+1):
                momsum[('c',k,l)]=momsum.get(('c',k,l),0)+np.sum(Ys[k-1]*Ys[l-1])
    n=Nsum
    mean={k:momsum[('m',k)]/n for k in range(1,MAX+1)}
    def cov(k,l):
        a,b=min(k,l),max(k,l)
        return momsum[('c',a,b)]/n - mean[k]*mean[l]
    vs={}
    for m in range(1,(MAX+1)//2+1):
        if 2*m-1<=MAX:
            val=cov(m,m)+2*sum(cov(k,2*m-k) for k in range(1,m))
            vs[m-1]=val
    return dict(delta=delta,vs=vs,mean=mean)

if __name__=="__main__":
    MAX=7
    print(f"generating node-shift functionals to order {MAX}...",flush=True)
    YF=make_Yfuncs(MAX)
    print("=== stability probe: v0,v1,v2,v3 at a few delta (h=1e-3, N=12000) ===",flush=True)
    print("  targets v0=0.134 v1=0.110 v2~0.11\n")
    for delta in [0.12,0.09,0.06]:
        t0=time.time(); d=run(delta,MAX,YF,N=12000,h=1e-3,seed=1)
        vs=d['vs']
        print(f"  d={delta:.3f} ({time.time()-t0:.0f}s): "+"  ".join(f"v{i}={vs[i]:+.4f}" for i in sorted(vs)),flush=True)
