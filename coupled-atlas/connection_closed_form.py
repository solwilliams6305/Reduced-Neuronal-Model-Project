"""
VERIFY the closed-form asymmetric connection of W's skeleton:
  oscillatory Jost log-deriv  L(lambda) = f_out'(0)/f_out(0)  =  e^{3ipi/4} * 2 Gamma(3/4 - i lam/4)/Gamma(1/4 - i lam/4)
  confining log-deriv         R(lambda) = -sqrt2 U'(a,0)/U(a,0) = 2 Gamma(3/4 - lam/4)/Gamma(1/4 - lam/4)
  => resonance condition (closed form):  e^{3ipi/4} G_i(lam) = G_r(lam),  G(mu)=Gamma(3/4-mu)/Gamma(1/4-mu).
Tests: (1) L(lam)/[imag-Weber ratio] == e^{3ipi/4} to high precision, many lam (real & complex), 2 step sizes;
       (2) the closed-form Gamma-equation's root == the dominant resonance 0.890-0.890i.
"""
import numpy as np
from scipy.special import loggamma
def L_osc(lam,xmax=20.0,h=2e-4):
    lam=np.asarray(lam,dtype=complex); x=xmax
    u=(x*x+lam)**-0.25; up=(1j*np.sqrt(x*x+lam)-0.5*(x/(x*x+lam)))*u; n=int(xmax/h)
    def f(x,u,up): return up,-(x*x+lam)*u
    for _ in range(n):
        a1,b1=f(x,u,up); a2,b2=f(x-h/2,u-h/2*a1,up-h/2*b1); a3,b3=f(x-h/2,u-h/2*a2,up-h/2*b2); a4,b4=f(x-h,u-h*a3,up-h*b3)
        u=u-h/6*(a1+2*a2+2*a3+a4); up=up-h/6*(b1+2*b2+2*b3+b4); x=x-h
    return up/u
def Gi(lam): return np.exp(loggamma(0.75-0.25j*lam)-loggamma(0.25-0.25j*lam))   # imag-Weber Gamma ratio
def Gr(lam): return np.exp(loggamma(0.75-0.25*lam)-loggamma(0.25-0.25*lam))     # real-Weber Gamma ratio

print("=== (1) L(lambda) / [2*Gi(lambda)] should equal e^{3i pi/4} = %.4f%+.4fi ===" % (np.cos(3*np.pi/4),np.sin(3*np.pi/4)))
pts=[0.0,0.7,1.5,2.5,3.5,-1.0,0.89-0.89j,2.0-1.0j,1.5+0.5j]
print(f"  {'lambda':>16} {'L/[2 Gi]':>22} {'|.|':>7} {'arg/pi':>8}")
for lam in pts:
    r=L_osc(lam)/(2*Gi(lam))
    print(f"  {lam!s:>16} {r.real:+.4f}{r.imag:+.4f}i  {abs(r):.4f}  {np.angle(r)/np.pi:+.4f}")
# precision + step-halving check at one point
r1=L_osc(1.5,h=2e-4)/(2*Gi(1.5)); r2=L_osc(1.5,h=1e-4)/(2*Gi(1.5))
print(f"  step check lam=1.5: h=2e-4 -> {r1.real:+.5f}{r1.imag:+.5f}i ; h=1e-4 -> {r2.real:+.5f}{r2.imag:+.5f}i (target -0.70711+0.70711i)")

print("\n=== (2) closed-form resonance condition  e^{3ipi/4} Gi(lam) = Gr(lam)  -> dominant root? ===")
def cond(lam): return np.exp(3j*np.pi/4)*Gi(lam)-Gr(lam)
# Newton from near the expected root
lam=0.9-0.9j
for _ in range(60):
    d=1e-6; c0=cond(lam); dc=0.5*((cond(lam+d)-c0)/d+(cond(lam+1j*d)-c0)/(1j*d))
    lam=lam-c0/dc
    if abs(cond(lam))<1e-12: break
print(f"  closed-form Gamma-equation root: lambda0 = {lam.real:+.4f}{lam.imag:+.4f}i   |cond|={abs(cond(lam)):.1e}")
print(f"  M1 (independent, FD diagonalization): 0.890-0.890i   -> match dist = {abs(lam-(0.890-0.890j)):.4f}")
print("\n  => constant phase e^{3ipi/4} (unit modulus, arg 3/4 pi) + Gamma-equation root = 0.890-0.890i")
print("     ==> BOTH connection factors are closed-form Weber Gamma-ratios (real-arg confining x imag-arg")
print("         oscillatory), and the asymmetric resonance condition of W's skeleton is CLOSED FORM.")
