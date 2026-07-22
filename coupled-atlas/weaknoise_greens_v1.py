"""
Program 2, Route 2b -- DETERMINISTIC (MC-free) computation of the weak-noise coefficients m1, Var(Y2), <Y1 Y3>, v1
via the Green's-function / Wick reduction. Extends weaknoise_greens.py (which validated v0=0.1343 deterministically).

Every field u_k is a k-fold iterated Ito integral against the SAME white noise, with causal kernel
        G(s,r) = u0(s) psi(r) - u0(r) psi(s)     [ L u = F  =>  u(Y)=\int_Y^{Y0} G(Y,s)F(s)ds ].
At the node Y*_0 (u0(Y*_0)=0):  G0(s):=G(Y*_0,s) = -c u0(s),  c=psi(Y*_0);   H(s):=d/dY G(Y,s)|_{Y*_0} = u0'* psi(s) - c' u0(s).
Node-shift functionals (all fields at Y*_0):
  Y1 = alpha \int u0^2 xi,      alpha=c/u0'*        (1st chaos)
  u1* = -u0'* Y1,   u1'* = \int H u0 xi,            u2* = \int\int_{s<r} P2 xi xi,  P2=-c u0(s)G(s,r)u0(r)
  Y2 = -(u1'* Y1 + u2*)/u0'*                        (0th+2nd chaos)
  Y3 = -( (1/6)V* u0'* Y1^3 + u1'* Y2 + (1/2)V* u1* Y1^2 + u2'* Y1 + u3* )/u0'*
Using chaos orthogonality <1st,3rd>=0 (=> <Y1 u3*>=0) the needed moments are finite deterministic integrals:
  m1 = <Y2> = -alpha \int H u0^3 / u0'*
  Var(Y2) = (1/u0'*^2)[ Var(P) + Var(Q) + 2Cov(P,Q) ],   P=u1'* Y1, Q=u2*
     Var(P)=alpha^2[ (\int H^2 u0^2)(\int u0^4) + (\int H u0^3)^2 ]
     Var(Q)=\int\int_{s<r} P2^2
     Cov(P,Q)=\int\int_{s<r} Sp*P2,   Sp(s,r)=alpha( f(s)g(r)+f(r)g(s) ), f=H u0, g=u0^2
  <Y1 Y3> = -(1/u0'*)[ (1/6)V* u0'* (3 v0^2) + <Y1 u1'* Y2> + (1/2)V*(-u0'* 3 v0^2) + <Y1^2 u2'*> ]
     <Y1 u1'* Y2> = -(1/u0'*)[ alpha^2( (\int u0^4)(\int H^2 u0^2) + 2(\int H u0^3)^2 ) + Cov(P,Q) ]
     <Y1^2 u2'*>  = 2 alpha^2 \int\int_{s<r} P2'(s,r) u0(s)^2 u0(r)^2,   P2'=H(s)G(s,r)u0(r)
  v1 = Var(Y2) + 2<Y1 Y3>
Targets (MC, weaknoise_v1.py): m1=+0.2124, Var(Y2)=0.0636, <Y1 Y3>=+0.0231, v1=+0.110.
"""
import numpy as np
from weaknoise_greens import build_backbone, node, V

def deterministic_coeffs(h=2e-4, hq=None, Y0=4.0, Yend=-6.0):
    Yg,u0,u0p,ps,psp=build_backbone(Y0=Y0,Yend=Yend,h=h)
    nd=node(Yg,u0,u0p,ps,psp)
    Ystar=nd['Ystar']; u0ps=nd['u0p_star']; c=nd['c']; cp=nd['cp']
    Vstar=V(Ystar)
    alpha=c/u0ps
    # restrict to s in [Y*_0, Y0], ascending order for quadrature
    mask=Yg>=Ystar
    s=Yg[mask][::-1].copy()
    U=u0[mask][::-1].copy(); PSI=ps[mask][::-1].copy()
    H=u0ps*PSI - cp*U            # H(s)
    # optionally subsample to a coarser quadrature grid hq for the O(n^2) doubles
    if hq is not None:
        step=max(1,int(round(hq/h)))
        s=s[::step]; U=U[::step]; PSI=PSI[::step]; H=H[::step]
    ds=np.gradient(s)  # ~ +h
    # --- single integrals (trapezoid) ---
    def I(fn): return np.trapz(fn,s)
    Iu4=I(U**4); IHu3=I(H*U**3); IH2u2=I(H**2*U**2)
    v0=alpha**2*Iu4
    m1=-alpha*IHu3/u0ps
    # --- Green matrix G(s,r)=U(s)PSI(r)-U(r)PSI(s) on grid ---
    n=len(s)
    G=np.outer(U,PSI)-np.outer(PSI,U)      # G[i,j]=U_i PSI_j - PSI_i U_j = G(s_i,s_j)
    # ordered region s<r  => i<j (upper triangle, i row, j col)
    w=ds  # quad weights ~ h
    W2=np.outer(w,w)
    triu=np.triu(np.ones((n,n)),k=1)       # strict upper i<j
    # P2(s,r)=-c U_i G_ij U_j  for i<j
    P2=-c*np.outer(U,U)*G
    # Sp(s,r)=alpha(f_i g_j + f_j g_i), f=H U, g=U^2
    f=H*U; g=U**2
    Sp=alpha*(np.outer(f,g)+np.outer(g,f))
    # P2'(s,r)=H_i G_ij U_j
    P2p=np.outer(H,np.ones(n))*G*np.outer(np.ones(n),U)  # H_i G_ij U_j
    def dint(M):  # \int\int_{i<j} M_ij w_i w_j
        return np.sum(triu*M*W2)
    VarQ=dint(P2**2)
    CovPQ=dint(Sp*P2)
    VarP=alpha**2*(IH2u2*Iu4 + IHu3**2)
    VarY2=(VarP+VarQ+2*CovPQ)/u0ps**2
    # <Y1 u1'* Y2>
    Y1u1pY2=-(1.0/u0ps)*(alpha**2*(Iu4*IH2u2+2*IHu3**2)+CovPQ)
    # <Y1^2 u2'*> = 2 alpha^2 \int\int_{i<j} P2'_ij U_i^2 U_j^2
    Y12u2p=2*alpha**2*dint(P2p*np.outer(U**2,U**2))
    Y1Y3=-(1.0/u0ps)*((1/6)*Vstar*u0ps*(3*v0**2)
                      + Y1u1pY2
                      + 0.5*Vstar*(-u0ps*3*v0**2)
                      + Y12u2p)
    v1=VarY2+2*Y1Y3
    return dict(Ystar=Ystar,u0ps=u0ps,c=c,alpha=alpha,Vstar=Vstar,
                v0=v0,m1=m1,VarY2=VarY2,VarP=VarP/u0ps**2,VarQ=VarQ/u0ps**2,
                CovPQ=CovPQ/u0ps**2,Y1Y3=Y1Y3,v1=v1,n=n)

if __name__=="__main__":
    print("=== deterministic (MC-free) weak-noise coefficients via Green/Wick ===")
    print("  targets (MC weaknoise_v1): v0=0.1339 m1=+0.2124 Var(Y2)=0.0636 <Y1Y3>=+0.0231 v1=+0.110\n")
    for hq in [8e-3,4e-3,2e-3]:
        d=deterministic_coeffs(h=2e-4,hq=hq)
        print(f"  hq={hq:.0e} (n={d['n']}): v0={d['v0']:.4f}  m1={d['m1']:+.4f}  "
              f"Var(Y2)={d['VarY2']:.4f}  <Y1Y3>={d['Y1Y3']:+.4f}  v1={d['v1']:+.4f}")
    print("\n  breakdown at finest hq:")
    d=deterministic_coeffs(h=2e-4,hq=2e-3)
    print(f"    Var(P)/u0ps^2={d['VarP']:.4f}  Var(Q)/u0ps^2={d['VarQ']:.4f}  2Cov(P,Q)/u0ps^2={2*d['CovPQ']:.4f}")
    print(f"    Vstar={d['Vstar']:.4f}  alpha={d['alpha']:.3e}  c={d['c']:.3e}  u0ps={d['u0ps']:.1f}")
