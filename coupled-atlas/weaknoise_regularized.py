"""
Program 2, Route 2b -- v2 (O(eta^6) variance coeff) via the REGULARIZED node-shift recursion (Route a/b bridge).

The node-shift recursion is noise-safe only through Y3; Y4 needs u1'''(Y*_0) = V'* u1 + V* u1' + xi(Y*_0) u0'*
and Y5 needs also xi'(Y*_0) -- white noise (and its derivative) at the fixed deterministic node. These are
distributional pathwise, but FINITE in the cumulant moments (boundary-Wick, kickoff sec 4/5bis).

Here we compute them by MOLLIFYING the noise: xi_delta(Y) = \int phi_delta(Y-s) dW(s), phi_delta a Gaussian of
width delta. Then xi_delta and all its Y-derivatives are finite smooth functions, the recursion runs through Y5
per realization, and v2(delta) -> v2 as delta -> 0. This realizes the boundary-Wick value numerically (no need to
guess the boundary 1/2 factor -- the limit supplies the correct constants, including the xi'(Y*_0) piece in Y5).

Node-shift functionals Y1..Y5 are generated symbolically (Taylor inversion of u(Y*_0+D)=0) in terms of
U[k,m]=u_k^{(m)}(Y*_0); the derivatives m>=2 use the field ODE u_k'' = V u_k + xi_delta u_{k-1} differentiated
(Leibniz), evaluating xi_delta^{(j)}(Y*_0) from the mollifier.

Validation ladder: v0 (0.1339), v1 (+0.110) must be recovered as delta->0; then v2 is the new number.
scipy for RK4 backbone; numpy for the vectorized mollified march.
"""
import numpy as np
import sympy as sp
from math import comb
from scipy.ndimage import gaussian_filter1d
from weaknoise_greens import build_backbone, node, V

# ---- symbolic node-shift functionals Y1..Y5 -> lambdified ----
def make_Yfuncs():
    eta=sp.symbols('eta'); K,M=6,7
    U=[[sp.symbols(f'U{k}_{m}') for m in range(M)] for k in range(K)]
    Y=[0]+[sp.symbols(f'Y{n}') for n in range(1,6)]
    D=sum(eta**k*Y[k] for k in range(1,6))
    expr=sum(eta**k*U[k][m]*D**m/sp.factorial(m) for k in range(K) for m in range(M))
    expr=sp.expand(expr).subs(U[0][0],0); sol={}
    for n in range(1,6):
        c=sp.expand(sp.expand(expr).coeff(eta,n).subs(sol))
        sol[Y[n]]=sp.expand(sp.solve(c,Y[n])[0])
    syms=[U[k][m] for k in range(K) for m in range(M)]
    funcs=[sp.lambdify(syms,sol[Y[n]].subs(U[0][0],0),'numpy') for n in range(1,6)]
    return funcs,(K,M)
YFUNCS,(K,M)=make_Yfuncs()

def Vd(Y,d):  # d-th derivative of V(Y)=sign(Y)Y^2 at Y<0 (V=-Y^2)
    return {0:-Y*Y,1:-2*Y,2:-2.0,3:0.0,4:0.0,5:0.0}[d]

def run(delta, N=20000, h=1e-3, seed=0, Y0=4.0, Yend=-2.6, batch=4000):
    # backbone + node (RK4, deterministic)
    Yg,u0,u0p,ps,psp=build_backbone(Y0=Y0,Yend=Yend,h=h)
    nd=node(Yg,u0,u0p,ps,psp)
    Ystar=nd['Ystar']; ist=nd['istar']; frac=nd['f']; u0ps=nd['u0p_star']
    Vs=[Vd(Ystar,d) for d in range(6)]
    nY=len(Yg)-1
    sigpix=delta/h
    # local analytic derivative kernels at the node (window +-6 delta)
    half=int(np.ceil(6*delta/h)); iL=ist-1
    lo=max(0,iL-half-2); hi=min(nY+1,iL+half+3)
    dgrid=Ystar-Yg[lo:hi]
    gg=np.exp(-dgrid**2/(2*delta**2))/(np.sqrt(2*np.pi)*delta); xx=dgrid/delta**2
    Kj=[gg, -xx*gg, (xx*xx-1/delta**2)*gg, (-xx**3+3*xx/delta**2)*gg]
    rng=np.random.default_rng(seed)
    acc={a:0.0 for a in ['Y1sq','Y2','Y2sq','Y1Y3','Y3sq','Y2Y4','Y1Y5','Y4','Ns']}
    done=0
    while done<N:
        nb=min(batch,N-done); done+=nb
        dW=np.sqrt(h)*rng.standard_normal((nb,nY+1))
        # smooth noise on full grid: xi_delta = (Gaussian filter of dW)/h
        xi0=gaussian_filter1d(dW, sigma=sigpix, axis=1, mode='constant')/h
        # xi^{(j)}(Y*_0) via local analytic kernels
        s_node=[ (dW[:,lo:hi]*Kj[j][None,:]).sum(1) for j in range(4) ]
        # march fields u1..u5 stacked: uu[k-1], upp[k-1] shape (5,nb). u_k'' = V u_k + xi0 u_{k-1}
        uu=np.zeros((5,nb)); upp=np.zeros((5,nb))
        dY=-h; snap={}
        for i in range(nY):
            Vi=V(Yg[i]); xi_i=xi0[:,i]
            src=np.empty((5,nb)); src[0]=u0[i]; src[1:]=uu[:4]      # u_{k-1}
            newup=upp+Vi*uu*dY+xi_i*src*dY
            uu=uu+upp*dY; upp=newup
            if i in (ist-1,ist):
                snap[i]=(uu.copy(),upp.copy())
        lo_s=snap[ist-1]; hi_s=snap[ist]
        Uval=lo_s[0]+frac*(hi_s[0]-lo_s[0])   # (5,nb) = U_{k,0}, k=1..5
        Uder=lo_s[1]+frac*(hi_s[1]-lo_s[1])   # (5,nb) = U_{k,1}, k=1..5
        Umat=np.zeros((6,6,nb))
        Umat[0,0]=0.0; Umat[0,1]=u0ps
        Umat[1:,0]=Uval; Umat[1:,1]=Uder
        s=s_node
        for k in range(6):
            for m in range(2,6):
                p=m-2; tot=np.zeros(nb)
                for j in range(p+1):
                    c=comb(p,j)
                    tot=tot+c*Vs[j]*Umat[k,p-j]
                    if k>=1:
                        xij=(s[j] if j<4 else 0.0)
                        tot=tot+c*xij*Umat[k-1,p-j]
                Umat[k,m]=tot
        Ufull=np.zeros((K,M,nb)); Ufull[:, :6, :]=Umat
        argl=[Ufull[k,m] for k in range(K) for m in range(M)]
        Y1,Y2,Y3,Y4,Y5=[f(*argl) for f in YFUNCS]
        acc['Y1sq']+=np.sum(Y1**2); acc['Y2']+=np.sum(Y2); acc['Y2sq']+=np.sum(Y2**2)
        acc['Y1Y3']+=np.sum(Y1*Y3); acc['Y3sq']+=np.sum(Y3**2); acc['Y2Y4']+=np.sum(Y2*Y4)
        acc['Y1Y5']+=np.sum(Y1*Y5); acc['Y4']+=np.sum(Y4); acc['Ns']+=nb
    n=acc['Ns']
    Y1sq=acc['Y1sq']/n; mY2=acc['Y2']/n; Y2sq=acc['Y2sq']/n; Y1Y3=acc['Y1Y3']/n
    Y3sq=acc['Y3sq']/n; Y2Y4=acc['Y2Y4']/n; Y1Y5=acc['Y1Y5']/n; mY4=acc['Y4']/n
    v0=Y1sq; m1=mY2; v1=(Y2sq-mY2**2)+2*Y1Y3
    covY2Y4=Y2Y4-mY2*mY4
    v2=Y3sq+2*covY2Y4+2*Y1Y5
    return dict(delta=delta,v0=v0,m1=m1,v1=v1,Y3sq=Y3sq,covY2Y4=covY2Y4,Y1Y5=Y1Y5,v2=v2,N=int(n))

if __name__=="__main__":
    import time
    print("=== regularized node-shift recursion: v0,v1 (validate) then v2 ===")
    print("  targets: v0=0.1339  v1=+0.110\n")
    for delta in [0.12,0.09,0.06]:
        t0=time.time()
        d=run(delta=delta,N=12000,h=1e-3,seed=1)
        print(f"  delta={delta:.3f} N={d['N']} ({time.time()-t0:.0f}s): v0={d['v0']:.4f} m1={d['m1']:+.4f} v1={d['v1']:+.4f}"
              f"  || Y3^2={d['Y3sq']:.4f} 2Cov24={2*d['covY2Y4']:+.4f} 2<Y1Y5>={2*d['Y1Y5']:+.4f} v2={d['v2']:+.4f}",flush=True)
