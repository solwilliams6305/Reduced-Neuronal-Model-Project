import numpy as np
from weaknoise_highorder import make_Yfuncs
import weaknoise_highorder as W

def run_full(delta,MAX,YF,N=16000,h=1e-3,seed=0,Y0=4.0,Yend=-2.6,batch=4000):
    from math import comb
    from numpy.polynomial.hermite_e import hermeval
    from scipy.ndimage import gaussian_filter1d
    from weaknoise_greens import build_backbone, node, V
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
    rng=np.random.default_rng(seed)
    allY=[]
    Nsum=0
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
        Ys=np.array([f(*argl) for f in funcs])  # (MAX,nb)
        allY.append(Ys)
    Y=np.concatenate(allY,axis=1)  # (MAX, Ntot)
    return Y

MAX=7; YF=make_Yfuncs(MAX)
print("delta   Var(Y2)  Var(Y3)  Var(Y4)  Var(Y5)  |  <Y4>=m2   Cov(Y1,Y5) Cov(Y2,Y4)  v3_pieces_sum")
for delta in [0.10,0.07,0.05,0.035]:
    Ys=[run_full(delta,MAX,YF,N=16000,h=1e-3,seed=s) for s in range(3)]
    Y=np.concatenate(Ys,axis=1)
    def var(m): return Y[m-1].var()
    def cov(a,b): return np.cov(Y[a-1],Y[b-1])[0,1]
    m2=Y[3].mean()
    v3=var(4)+2*cov(3,5)+2*cov(2,6)+2*cov(1,7)
    print(f"{delta:.3f}  {var(2):8.4f} {var(3):8.4f} {var(4):8.4f} {var(5):8.4f}  |  {m2:+8.4f}  {cov(1,5):+9.4f} {cov(2,4):+9.4f}  v3={v3:+.3f}")
