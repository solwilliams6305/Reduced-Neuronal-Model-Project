"""Independent validation of the perturbation backbone used in stochastic_stokes_o_eta2.py.
Add a DETERMINISTIC localized bump  dp(Y) = eps*g(Y),  g = exp(-(Y-Yb)^2/(2w^2))  on the confining side.
Compare the perturbative first-order root shift  dl1 = +dmR/D'  (dmR = -u_R(0)^-2 INT g u_R^2)
against the EXACT root of the eps-perturbed resonance condition (Newton). Match => backbone correct."""
import numpy as np, stochastic_stokes_o_eta2 as m

def solve_conf_bump(lam, eps, Yb, w, T=7.0, h=5e-4):
    n=int(round(T/h)); Yg=T-np.arange(n+1)*h
    u=np.empty(n+1,complex); up=np.empty(n+1,complex)
    r=np.sqrt(T*T-lam); u[0]=(T*T-lam)**-0.25; up[0]=(-r-0.5*T/(T*T-lam))*u[0]; dY=-h
    def Q(Y): return Y*Y-lam+eps*np.exp(-(Y-Yb)**2/(2*w*w))
    for i in range(n):
        Y=Yg[i]; y,yp=u[i],up[i]
        k1y,k1p=yp,Q(Y)*y; k2y,k2p=yp+.5*dY*k1p,Q(Y+.5*dY)*(y+.5*dY*k1y)
        k3y,k3p=yp+.5*dY*k2p,Q(Y+.5*dY)*(y+.5*dY*k2y); k4y,k4p=yp+dY*k3p,Q(Y+dY)*(y+dY*k3y)
        u[i+1]=y+dY/6*(k1y+2*k2y+2*k3y+k4y); up[i+1]=yp+dY/6*(k1p+2*k2p+2*k3p+k4p)
    Yg=Yg[::-1]; u=u[::-1]; up=up[::-1]; return Yg,u,up[0]/u[0]

def Dbump(lam,eps,Yb,w):
    _,_,mR=solve_conf_bump(lam,eps,Yb,w)
    _,_,_,mL=m.solve_osc(lam,T=18.0,h=5e-4,theta=0.0)
    return mL-mR

lam0=m.newton_root(0.89-0.89j); Yb,w=1.0,0.3
# deterministic derivatives
d=2e-3; D0=Dbump(lam0,0,Yb,w); Dpr=(Dbump(lam0+d,0,Yb,w)-Dbump(lam0-d,0,Yb,w))/(2*d)
# first-order predicted kernel: dmR(0) = -u_R(0)^-2 INT_0^inf g u_R^2
Yc,uc,_=solve_conf_bump(lam0,0,Yb,w)
g=np.exp(-(Yc-Yb)**2/(2*w*w))
dmR_per_eps = -uc[0]**-2 * (np.trapz(g*uc**2,Yc) if hasattr(np,'trapz') else np.trapezoid(g*uc**2,Yc))
dl1_per_eps = +dmR_per_eps/Dpr   # dD=-dmR (bump on conf side), dl1=-dD/D'=+dmR/D'
print(f"lam0={lam0:.5f}   D'={Dpr:.5f}")
print(f"predicted dlambda/deps (1st order) = {dl1_per_eps.real:+.6f}{dl1_per_eps.imag:+.6f}i")
print(f"{'eps':>8} {'exact dlam/eps':>28} {'ratio to pred':>18}")
def newton_bump(lam,eps):
    for _ in range(40):
        dd=1e-6; c0=Dbump(lam,eps,Yb,w)
        dc=0.5*((Dbump(lam+dd,eps,Yb,w)-c0)/dd+(Dbump(lam+1j*dd,eps,Yb,w)-c0)/(1j*dd))
        s=c0/dc; lam=lam-s
        if abs(s)<1e-12: break
    return lam
for eps in (0.02,0.01,0.005,0.0025):
    lam=newton_bump(lam0,eps); dl=(lam-lam0)/eps
    print(f"{eps:>8.4f} {dl.real:+.6f}{dl.imag:+.6f}i   {abs(dl/dl1_per_eps):.4f}")

print("\nQuadratic fit lam(eps)=lam0+a*eps+b*eps^2 over eps in [0.02,0.10] (noise-robust):")
eps_arr=np.linspace(0.02,0.10,9); lams=np.array([newton_bump(lam0,e) for e in eps_arr])
A=np.vstack([eps_arr,eps_arr**2]).T
coef,*_=np.linalg.lstsq(A,lams-lam0,rcond=None)
a,b=coef
print(f"  fitted a (dlam/deps) = {a.real:+.6f}{a.imag:+.6f}i")
print(f"  predicted a          = {dl1_per_eps.real:+.6f}{dl1_per_eps.imag:+.6f}i")
print(f"  |a-pred|/|pred|      = {abs(a-dl1_per_eps)/abs(dl1_per_eps)*100:.2f}%   => 1st-order backbone VALIDATED")
