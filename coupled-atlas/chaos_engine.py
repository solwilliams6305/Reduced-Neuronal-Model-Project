"""
Deterministic Wiener-chaos moment engine for the weak-noise coefficients v_n (Program 2).

Every node functional is an element of the Wiener chaos: element = sum_p I_p(h_p), h_p a SYMMETRIC p-kernel.
We store kernels in WEIGHT-ABSORBED form  h_abs[i..] = h[i..] * prod sqrt(w[i]),  so
    <I_p(f),I_p(g)> = p! * sum(f_abs*g_abs)          (plain sum, weights absorbed)
    I_p(f)*I_q(g)   = sum_r C(p,r)C(q,r) r! I_{p+q-2r}( SYM( contract_r(f_abs,g_abs) ) )   (Ito product)
The node fields are Markovian bond-chains (causal Green's function G(s,r)=u0(s)psi(r)-u0(r)psi(s)):
    u_k(node): ordered kernel  G(node,s1) G(s1,s2)...G(s_{k-1},s_k) u0(s_k),  node<s1<...<s_k
    u_k'(node): same with G(node,s1) -> dG/dY(node,s1).
Validate on v0=<Y1^2>=0.1339 and v1=Var(Y2)+2<Y1Y3>=0.110 (boundary-free) before adding boundary noise.
"""
import numpy as np, itertools
from math import comb, factorial
from weaknoise_greens import build_backbone, node, V

# ---------- chaos element ----------
class Ch:
    def __init__(self, terms=None):  # terms: dict order-> absorbed symmetric tensor (order 0 -> scalar array shape ())
        self.t={} if terms is None else {k:v for k,v in terms.items()}
    @staticmethod
    def scalar(c):
        return Ch({0:np.array(float(c))})
    def __add__(a,b):
        o=Ch(a.t)
        for k,v in b.t.items(): o.t[k]=o.t.get(k,0)+v
        return o
    def __rmul__(self,c): return self.__mul__(c)
    def scale(self,c):
        return Ch({k:c*v for k,v in self.t.items()})

def _sym(T):
    p=T.ndim
    if p<2: return T
    s=np.zeros_like(T);
    for perm in itertools.permutations(range(p)):
        s+=np.transpose(T,perm)
    return s/factorial(p)

def cmul(A,B,MAXORD):
    out=Ch()
    for p,f in A.t.items():
        for q,g in B.t.items():
            for r in range(min(p,q)+1):
                m=p+q-2*r
                if m>MAXORD: continue
                coef=comb(p,r)*comb(q,r)*factorial(r)
                if p==0 and q==0:
                    val=coef*(f*g); out.t[0]=out.t.get(0,0)+val; continue
                if r==0:
                    # outer product
                    if p==0: K=f*g
                    elif q==0: K=g*f
                    else: K=np.tensordot(f,g,axes=0)
                else:
                    # contract first r indices of f with first r of g
                    fa=np.arange(r); ga=np.arange(r)
                    K=np.tensordot(f,g,axes=(list(range(r)),list(range(r))))
                    if p==r and q==r: K=np.array(K)  # scalar
                K=_sym(K) if (np.ndim(K)>=2) else K
                out.t[m]=out.t.get(m,0)+coef*K
    return out

def cpow(A,e,MAXORD):
    r=Ch.scalar(1.0)
    for _ in range(e): r=cmul(r,A,MAXORD)
    return r

def inner(A,B):  # <A,B> full expectation (sum over matching orders, p! * sum(prod))
    s=0.0
    for k,f in A.t.items():
        if k in B.t:
            g=B.t[k]
            if k==0: s+=float(f*B.t[0])
            else: s+=factorial(k)*np.sum(f*g)
    return s

# ---------- build node-field kernels ----------
def build(n_grid=60, MAXORD=3, Y0=4.0, Yend=-2.6, hbb=2e-4):
    Yg,u0,u0p,ps,psp=build_backbone(Y0=Y0,Yend=Yend,h=hbb)
    nd=node(Yg,u0,u0p,ps,psp); Yst=nd['Ystar']; u0ps=nd['u0p_star']; c=nd['c']; cp=nd['cp']
    # resample onto n_grid points over [Yst,Y0]
    s=np.linspace(Yst,Y0,n_grid); ds=s[1]-s[0]
    U=np.interp(s,Yg[::-1],u0[::-1]); PSI=np.interp(s,Yg[::-1],ps[::-1])
    w=np.full(n_grid,ds); w[0]*=0.5; w[-1]*=0.5    # trapezoid
    sw=np.sqrt(w)
    G=np.outer(U,PSI)-np.outer(PSI,U)              # G[i,j]=G(s_i,s_j)
    Gnode=-c*U                                      # G(node,s)= -c u0(s)  (node value)
    H=u0ps*PSI-cp*U                                 # dG/dY(node,s)
    Vst=V(Yst)
    # ordered bond chain for u_k(node): node<s1<...<sk
    triu=np.triu(np.ones((n_grid,n_grid)),1)        # i<j
    def chain(head):   # head[i]=G(node,s_i) or H(s_i); returns dict k-> Ch for u_k value (k=1..MAXORD)
        out={}
        # k=1: f1(s)=head(s)*u0(s)
        f1=head*U
        out[1]=Ch({1:(f1*sw)/factorial(1)})
        # build ordered tensors iteratively: T_k[s1..sk]=head[s1]*G[s1,s2]*..*G[s_{k-1},sk]*u0[sk], ordered
        # represent as: prev vector over last index carrying head*prod bonds up to s_{k-1}, then extend
        # do explicit for k=2,3 (needed); general via loop building ordered dense tensor
        Gmask=G*triu   # bond with ordering i<j baked in
        for k in range(2,MAXORD+1):
            # vectorized ordered chain: P[...,i_{s-1},i_s] = P[...,i_{s-1}] * G[i_{s-1},i_s] * (i_{s-1}<i_s)
            P=head.copy()                          # shape (n,), index i1
            for _ in range(2,k+1):
                P=P[...,None]*Gmask                # (...,n,1)*(n,n) -> (...,n,n), ordering baked in
            T=P*U                                  # tail u0(s_k) on last axis
            Th=_sym(T)*_wten(sw,k)                 # symmetrize + absorb weights
            out[k]=Ch({k:Th})
        return out
    Uk0=chain(Gnode)   # u_k(node)
    Uk1=chain(H)       # u_k'(node)
    return dict(s=s,w=w,sw=sw,U=U,PSI=PSI,G=G,Gnode=Gnode,H=H,u0ps=u0ps,c=c,cp=cp,Vst=Vst,
                Uk0=Uk0,Uk1=Uk1,MAXORD=MAXORD,n=n_grid)

def _wten(sw,k):
    T=sw.copy()
    for _ in range(k-1): T=np.multiply.outer(T,sw)
    return T

if __name__=="__main__":
    for n in [40,60,80]:
        B=build(n_grid=n,MAXORD=3)
        u0ps=B['u0ps']; Uk0=B['Uk0']; Uk1=B['Uk1']; Vst=B['Vst']
        # Y1 = -U10/u0ps  (1-chaos)
        U10=Uk0[1]; U11=Uk1[1]; U20=Uk0[2]; U21=Uk1[2]; U30=Uk0[3]
        Y1=U10.scale(-1.0/u0ps)
        v0=inner(Y1,Y1)
        # Y2 = -(U11*Y1 + U20)/u0ps
        Y2=cmul(U11,Y1,3).__add__(U20).scale(-1.0/u0ps)
        m1=Y2.t.get(0,np.array(0.0));
        VarY2=inner(Y2,Y2)-float(m1)**2
        # Y3 = -( (1/6)Vst u0ps Y1^3 + U11*Y2 + (1/2)Vst U10 Y1^2 + U21*Y1 + U30)/u0ps
        Y1sq=cpow(Y1,2,3); Y1cu=cpow(Y1,3,3)
        Y3=(Y1cu.scale(Vst*u0ps/6.0)
            .__add__(cmul(U11,Y2,3))
            .__add__(cmul(U10,Y1sq,3).scale(0.5*Vst))
            .__add__(cmul(U21,Y1,3))
            .__add__(U30)).scale(-1.0/u0ps)
        Y1Y3=inner(Y1,Y3)
        v1=VarY2+2*Y1Y3
        print(f" n={n}: v0={v0:.4f}(t0.1339)  m1={float(m1):+.4f}(t0.212)  Var(Y2)={VarY2:.4f}(t0.0636)  <Y1Y3>={Y1Y3:+.4f}(t0.0231)  v1={v1:+.4f}(t0.110)")
