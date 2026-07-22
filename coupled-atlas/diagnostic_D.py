"""
Diagnostic: is the closed-form Wronskian D(lam) correct, and which left branch is outgoing?
Compare D_closed (gamma zero-values) vs D_ODE (numerically integrate the two Jost solutions to Y=0,
take the Wronskian). Two left branches: a_L=-i lam/2, c=sqrt2 e^{-ipi/4} (~e^{+iY^2/2}); and conj.
If a branch's D_closed == D_ODE (up to lam-independent const), that branch's gamma form is correct.
"""
import cmath, numpy as np
_c=[0.99999999999980993,676.5203681218851,-1259.1392167224028,771.32342877765313,
    -176.61502916214059,12.507343278686905,-0.13857109526572012,9.9843695780195716e-6,1.5056327351493116e-7]
def clogG(z):
    z=complex(z)
    if z.real<0.5: return cmath.log(cmath.pi)-cmath.log(cmath.sin(cmath.pi*z))-clogG(1-z)
    z-=1; x=_c[0]
    for i in range(1,9): x+=_c[i]/(z+i)
    t=z+7.5; return 0.5*cmath.log(2*cmath.pi)+(z+0.5)*cmath.log(t)-t+cmath.log(x)
def rgamma(z):
    z=complex(z)
    if z.real>=0.5: return cmath.exp(-clogG(z))
    return cmath.sin(cmath.pi*z)/cmath.pi*cmath.exp(clogG(1-z))   # entire 1/Gamma (0 at poles)
def U0(a):  return cmath.pi**0.5 * 2**(-0.25-a/2) * rgamma(0.75+a/2)   # U(a,0)  DLMF 12.2.6
def Up0(a): return -cmath.pi**0.5 * 2**(0.25-a/2) * rgamma(0.25+a/2)   # U'(a,0) DLMF 12.2.7
def D_closed(lam,branch):
    aR=-lam/2
    if branch==0: aL=-1j*lam/2; c=cmath.sqrt(2)*cmath.exp(-1j*cmath.pi/4)
    else:         aL=+1j*lam/2; c=cmath.sqrt(2)*cmath.exp(+1j*cmath.pi/4)
    uR,upR=U0(aR), cmath.sqrt(2)*Up0(aR)
    uL,upL=U0(aL), c*Up0(aL)
    return uR*upL-upR*uL
def D_ODE(lam,branch,Y0=6.0,dt=1e-3):
    # right recessive from +Y0 to 0
    def integ(Y0s,u,up,Yend):
        Y=Y0s; n=int(round(abs(Y0s-Yend)/dt)); h=-dt if Y0s>Yend else dt
        for _ in range(n):
            def Q(Yv): return (1 if Yv>0 else -1)*Yv*Yv-lam
            k1u,k1p=up,Q(Y)*u
            k2u,k2p=up+.5*h*k1p,Q(Y+.5*h)*(u+.5*h*k1u)
            k3u,k3p=up+.5*h*k2p,Q(Y+.5*h)*(u+.5*h*k2u)
            k4u,k4p=up+h*k3p,Q(Y+h)*(u+h*k3u)
            u+=(h/6)*(k1u+2*k2u+2*k3u+k4u); up+=(h/6)*(k1p+2*k2p+2*k3p+k4p); Y+=h
        return u,up
    uR,upR=integ(Y0,1+0j,-cmath.sqrt(Y0*Y0-lam),0.0)
    s=1 if branch==0 else -1                      # e^{+iY^2/2} (out) vs e^{-iY^2/2}
    uL,upL=integ(-Y0,1+0j,(-1j*s)*(-Y0)*1,0.0)    # u'(-Y0)=i s Y u ; Y=-Y0 => -i s Y0
    return uR*upL-upR*uL
for lam in [1.0+0j,2.0+0j,1.5-1.0j,3.0-1.0j]:
    print(f"lam={lam}:")
    for b in (0,1):
        dc,do=D_closed(lam,b),D_ODE(lam,b)
        print(f"   branch{b}: D_closed/D_ODE = {dc/do:.4f}   (|ratio|={abs(dc/do):.4f}, arg={cmath.phase(dc/do):+.3f})")
