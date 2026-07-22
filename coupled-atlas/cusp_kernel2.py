"""
Stage 1, alternative kernel: the TWO-Jost-solution form K=[u+(x)u-(y)-u+(y)u-(x)]/((x-y)W),
u+ recessive-right, u- decaying-left (both solve phi''=sign(t)t^2 e4 phi on the rotated contour), W=Wronskian.
Check theta-independence + cumulants vs FP W (skew +0.60, exk -0.24). numpy only.
"""
import numpy as np
def integ(theta,t0,u0,up0,tend):
    e4=np.exp(4j*theta); dt=2e-3; n=int(round(abs(t0-tend)/dt)); h=-dt if t0>tend else dt
    ts=t0+np.arange(n+1)*h; u=np.empty(n+1,complex); up=np.empty(n+1,complex); a,b=u0,up0; u[0]=a; up[0]=b
    for i in range(n):
        t=ts[i]
        def Q(tv): return np.sign(tv)*tv*tv*e4
        k1u,k1p=b,Q(t)*a; k2u,k2p=b+.5*h*k1p,Q(t+.5*h)*(a+.5*h*k1u)
        k3u,k3p=b+.5*h*k2p,Q(t+.5*h)*(a+.5*h*k2u); k4u,k4p=b+h*k3p,Q(t+h)*(a+h*k3u)
        a+=(h/6)*(k1u+2*k2u+2*k3u+k4u); b+=(h/6)*(k1p+2*k2p+2*k3p+k4p); u[i+1]=a; up[i+1]=b
    return ts,u,up
def cumulants(theta):
    e2=np.exp(2j*theta)
    tp,up_,upp=integ(theta,6.0,1+0j,-6*e2,-6.0)          # recessive right (full contour)
    tm,um_,ump=integ(theta,-6.0,1+0j,-6j*e2,6.0)         # decaying left (full contour)
    # interpolate both on common grid
    g=np.linspace(-6,6,4000)
    def ipv(ts,arr):
        o=np.argsort(ts.real); return np.interp(g,ts.real[o],arr.real[o])+1j*np.interp(g,ts.real[o],arr.imag[o])
    uP=ipv(tp,up_); uPp=ipv(tp,upp); uM=ipv(tm,um_); uMp=ipv(tm,ump)
    W=uP*uMp-uPp*uM; W=np.median(W.real)+1j*np.median(W.imag)   # Wronskian (const)
    def at(arr,x): return np.interp(x,g,arr.real)+1j*np.interp(x,g,arr.imag)
    def Fdet(s,T0=6.0,n=70):
        t,w=np.polynomial.legendre.leggauss(n); xq=0.5*(T0-s)*t+0.5*(T0+s); wq=0.5*(T0-s)*w*np.exp(1j*theta)
        a=at(uP,xq); b=at(uM,xq)
        dx=xq[:,None]-xq[None,:]; np.fill_diagonal(dx,1.0)
        Km=(a[:,None]*b[None,:]-a[None,:]*b[:,None])/dx/W
        np.fill_diagonal(Km,0.0)
        sw=np.sqrt(wq); return np.linalg.det(np.eye(n)-sw[:,None]*Km*sw[None,:])
    ss=np.linspace(-4.5,1.5,46); F=np.array([Fdet(s).real for s in ss])
    F=np.clip(F,0,1); f=np.gradient(F,ss); f=np.clip(f,0,None)
    if np.trapz(f,ss)<1e-6: return None
    f/=np.trapz(f,ss); m=np.trapz(ss*f,ss); v=np.trapz((ss-m)**2*f,ss); sd=v**0.5
    return m,sd,np.trapz((ss-m)**3*f,ss)/sd**3,np.trapz((ss-m)**4*f,ss)/v**2-3
print("two-Jost kernel det: cumulants vs FP W (skew +0.60, exk -0.24); check theta-independence")
for th in [0.45,0.55,0.65]:
    r=cumulants(th)
    print(f"  theta={th:.2f}: "+("degenerate" if r is None else f"mean={r[0]:+.3f} std={r[1]:.3f} skew={r[2]:+.3f} exk={r[3]:+.3f}"))
