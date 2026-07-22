"""
Stage 1: cusp edge gap determinant. Kernel from the ROTATED recessive Weber solution psi^theta of the
cusp edge ODE  phi'' = sign(z) z^2 phi  (z=t e^{i th}), normalized to U(0,sqrt2 z) ~ e^{-z^2/2}/(2^{1/4} sqrt z).
K^th(x,y)=[psi(x)psi'(y)-psi'(x)psi(y)]/(x-y) (CD form, x,y on the rotated contour).
F(s)=det(1-K^th)|_{(s,T0)}.  Validate standardized cumulants vs FP W: skew +0.60, exk -0.24.
Pipeline already verified on Airy->TW2. numpy only.
"""
import numpy as np
def build_psi(theta,t0=6.0,dt=2e-3,tend=-6.0):
    e4=np.exp(4j*theta); n=int(round((t0-tend)/dt)); ts=t0-np.arange(n+1)*dt
    z0=t0*np.exp(1j*theta)
    psi0=np.exp(-z0*z0/2)/(2**0.25*np.sqrt(z0)); psip0=-t0*np.exp(2j*theta)*psi0
    psi=np.empty(n+1,complex); psit=np.empty(n+1,complex); u,up=psi0,psip0; psi[0]=u; psit[0]=up
    for i in range(n):
        t=ts[i]; h=-dt
        def Q(tv): return np.sign(tv)*tv*tv*e4
        k1u,k1p=up,Q(t)*u; k2u,k2p=up+.5*h*k1p,Q(t+.5*h)*(u+.5*h*k1u)
        k3u,k3p=up+.5*h*k2p,Q(t+.5*h)*(u+.5*h*k2u); k4u,k4p=up+h*k3p,Q(t+h)*(u+h*k3u)
        u+=(h/6)*(k1u+2*k2u+2*k3u+k4u); up+=(h/6)*(k1p+2*k2p+2*k3p+k4p); psi[i+1]=u; psit[i+1]=up
        if abs(u)>1e12: pass
    return ts[::-1],psi[::-1],psit[::-1]
def cumulants(theta):
    tg,psi,psit=build_psi(theta); e4=np.exp(4j*theta)
    def ipv(arr,x): return np.interp(x,tg,arr.real)+1j*np.interp(x,tg,arr.imag)
    def Fdet(s,T0=6.0,n=70):
        t,w=np.polynomial.legendre.leggauss(n); xq=0.5*(T0-s)*t+0.5*(T0+s); wq=0.5*(T0-s)*w*np.exp(1j*theta)
        p=ipv(psi,xq); pt=ipv(psit,xq)                       # psi, psi' at nodes (vectorized)
        dx=xq[:,None]-xq[None,:]; np.fill_diagonal(dx,1.0)
        Km=(p[:,None]*pt[None,:]-pt[:,None]*p[None,:])/dx
        np.fill_diagonal(Km, pt*pt - p*np.sign(xq)*xq*xq*e4*p)  # CD diagonal
        sw=np.sqrt(wq); return np.linalg.det(np.eye(n)-sw[:,None]*Km*sw[None,:])
    ss=np.linspace(-4.5,1.5,46); F=np.array([Fdet(s).real for s in ss])
    F=np.clip(F,0,1); f=np.gradient(F,ss); f=np.clip(f,0,None)
    if np.trapz(f,ss)<1e-6: return None
    f/=np.trapz(f,ss); m=np.trapz(ss*f,ss); v=np.trapz((ss-m)**2*f,ss); sd=v**0.5
    sk=np.trapz((ss-m)**3*f,ss)/sd**3; ek=np.trapz((ss-m)**4*f,ss)/v**2-3
    return m,sd,sk,ek
print("cusp gap determinant det(1-K^theta): standardized cumulants vs FP W (skew +0.60, exk -0.24)")
for th in [0.45,0.55,0.65]:
    r=cumulants(th)
    if r: print(f"  theta={th:.2f}: mean={r[0]:+.3f} std={r[1]:.3f} skew={r[2]:+.3f} exk={r[3]:+.3f}")
    else: print(f"  theta={th:.2f}: degenerate")
