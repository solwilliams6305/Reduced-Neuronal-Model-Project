"""
Decisive test: do D_closed (gamma Wronskian) and D_ODE (direct integration) share ZEROS (resonances)?
Zeros are normalization-free. Compare lowest zeros of each over the complex-lambda grid.
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
    return cmath.sin(cmath.pi*z)/cmath.pi*cmath.exp(clogG(1-z))
def U0(a):  return cmath.pi**0.5*2**(-0.25-a/2)*rgamma(0.75+a/2)
def Up0(a): return -cmath.pi**0.5*2**(0.25-a/2)*rgamma(0.25+a/2)
def Dc(lam,b):
    aR=-lam/2; aL=(-1j if b==0 else 1j)*lam/2; c=cmath.sqrt(2)*cmath.exp((-1j if b==0 else 1j)*cmath.pi/4)
    uR,upR=U0(aR),cmath.sqrt(2)*Up0(aR); uL,upL=U0(aL),c*Up0(aL)
    w=uR*upL-upR*uL; return w/(abs(uR*upL)+abs(upR*uL)+1e-30)      # normalized
# vectorized ODE Wronskian
def integ(lam,Y0s,u,up,Yend,dt=1.5e-3):
    Y=Y0s; n=int(round(abs(Y0s-Yend)/dt)); h=-dt if Y0s>Yend else dt
    for _ in range(n):
        QY=np.sign(Y)*Y*Y-lam; QYh=np.sign(Y+.5*h)*(Y+.5*h)**2-lam; QYf=np.sign(Y+h)*(Y+h)**2-lam
        k1u,k1p=up,QY*u; k2u,k2p=up+.5*h*k1p,QYh*(u+.5*h*k1u)
        k3u,k3p=up+.5*h*k2p,QYh*(u+.5*h*k2u); k4u,k4p=up+h*k3p,QYf*(u+h*k3u)
        u=u+(h/6)*(k1u+2*k2u+2*k3u+k4u); up=up+(h/6)*(k1p+2*k2p+2*k3p+k4p); Y+=h
        m=np.abs(u)>1e8
        if m.any(): u[m]/=1e8; up[m]/=1e8
    return u,up
def Do(lam,b,Y0=6.0):
    uR,upR=integ(lam,Y0,np.ones_like(lam),-np.sqrt(Y0*Y0-lam),0.0)
    s=-1j if b==0 else 1j                       # e^{+iY^2/2}: u'/u=iY -> at -Y0: -iY0
    uL,upL=integ(lam,-Y0,np.ones_like(lam),(s)*(-Y0)*np.ones_like(lam),0.0)
    w=uR*upL-upR*uL; return w/(np.abs(uR*upL)+np.abs(upR*uL)+1e-30)
EE,GG=np.meshgrid(np.linspace(0.3,4.5,43),np.linspace(0.05,2.5,26)); lam=(EE-1j*GG).ravel()
for b in (0,1):
    do=np.abs(Do(lam,b)); dc=np.array([abs(Dc(l,b)) for l in lam])
    print(f"\n--- branch {b} ---")
    for name,arr in [("D_ODE",do),("D_closed",dc)]:
        o=np.argsort(arr); seen=[]
        out=[]
        for idx in o:
            e,g=EE.ravel()[idx],GG.ravel()[idx]
            if all(abs(e-s[0])>0.6 or abs(g-s[1])>0.6 for s in seen):
                seen.append((e,g)); out.append(f"{e:.2f}-{g:.2f}i({arr[idx]:.2f})")
            if len(seen)>=4: break
        print(f"  {name} lowest |.| at: {', '.join(out)}")
