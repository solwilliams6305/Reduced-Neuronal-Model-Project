"""
Program 2, Route 2b -- DETERMINISTIC Green's-function / Wick reformulation of the weak-noise expansion.

The field perturbations solve L u_k = xi u_{k-1}, L = d^2/dY^2 - V, recessive (zero) data at Y0. With the causal
Green's function of L (zero Cauchy data at Y0),
        u(Y) = \int_Y^{Y0} G(Y,s) F(s) ds  solves  L u = F,
        G(Y,s) = u0(Y) psi(s) - u0(s) psi(Y),   psi(Y0)=0, psi'(Y0)=1  (Wronskian W = u0 psi' - u0' psi = 1),
every u_k is a k-fold iterated Ito integral against the SAME white noise xi. Hence every cumulant coefficient of
the first-node law is a DETERMINISTIC multiple integral (Wick's theorem) of products of G and u0 -- NO Monte Carlo.

This module builds u0, psi, G numerically and validates the reformulation on the leading variance coefficient
        v0 = <Y1^2> = (c/u0*)^2 \int_{Y*_0}^{Y0} u0(s)^4 ds,   c = psi(Y*_0),  u0* = u0'(Y*_0),
which must reproduce the closed-form/MC value 0.1339. This is the on-ramp: once the kernel is validated on v0,
the same deterministic Wick integrals give v1 (check vs +0.110) and then v2 via boundary-Wick for the xi(Y*_0)
term (Route a of the kickoff).

scipy for a clean RK4 backbone; numpy for the quadratures.
"""
import numpy as np

def V(Y): return np.sign(Y)*Y*Y

def build_backbone(Y0=4.0, Yend=-6.0, h=2e-4):
    """RK4-march u0 (recessive: u0(Y0)=1,u0'=-Y0) and psi (psi(Y0)=0,psi'=1) downward. Return grid + fields."""
    nY=int(round((Y0-Yend)/h))
    Yg=Y0-np.arange(nY+1)*h
    u0=np.empty(nY+1); u0p=np.empty(nY+1); ps=np.empty(nY+1); psp=np.empty(nY+1)
    u0[0],u0p[0]=1.0,-Y0
    ps[0],psp[0]=0.0,1.0
    dY=-h
    def rhs(Y,y,yp): return yp, V(Y)*y   # y''=V y
    for i in range(nY):
        Y=Yg[i]
        for (arr,arrp) in ((u0,u0p),(ps,psp)):
            y,yp=arr[i],arrp[i]
            k1y,k1p=rhs(Y,y,yp)
            k2y,k2p=rhs(Y+0.5*dY,y+0.5*dY*k1y,yp+0.5*dY*k1p)
            k3y,k3p=rhs(Y+0.5*dY,y+0.5*dY*k2y,yp+0.5*dY*k2p)
            k4y,k4p=rhs(Y+dY,y+dY*k3y,yp+dY*k3p)
            arr[i+1]=y+dY/6*(k1y+2*k2y+2*k3y+k4y)
            arrp[i+1]=yp+dY/6*(k1p+2*k2p+2*k3p+k4p)
    return Yg,u0,u0p,ps,psp

def node(Yg,u0,u0p,ps,psp):
    istar=next(i for i in range(1,len(Yg)) if Yg[i]<0 and u0[i-1]*u0[i]<0)
    f=u0[istar-1]/(u0[istar-1]-u0[istar])
    Ystar=Yg[istar-1]+f*(Yg[istar]-Yg[istar-1])
    def interp(a): return a[istar-1]+f*(a[istar]-a[istar-1])
    return dict(istar=istar,f=f,Ystar=Ystar,u0p_star=interp(u0p),
                c=interp(ps),cp=interp(psp))

if __name__=="__main__":
    print("=== deterministic Green's-function reformulation: validate v0 (no MC) ===")
    for h in [4e-4,2e-4,1e-4]:
        Yg,u0,u0p,ps,psp=build_backbone(h=h)
        nd=node(Yg,u0,u0p,ps,psp)
        Ystar=nd['Ystar']; u0ps=nd['u0p_star']; c=nd['c']
        # integrate u0^4 over [Y*_0, Y0]; grid is descending, restrict to s>=Ystar
        mask=Yg>=Ystar
        s=Yg[mask][::-1]; integrand=(u0[mask][::-1])**4
        I4=np.trapz(integrand,s) if hasattr(np,'trapz') else np.trapezoid(integrand,s)
        v0=(c/u0ps)**2*I4
        print(f"  h={h:.0e}: Y*_0={Ystar:.5f}  u0'*={u0ps:.4f}  c=psi(Y*_0)={c:.4f}  "
              f"int u0^4={I4:.5f}  ->  v0={v0:.5f}   (target 0.1339)")
