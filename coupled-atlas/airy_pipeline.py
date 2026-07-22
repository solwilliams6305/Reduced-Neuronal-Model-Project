"""
Pipeline verification: reproduce Airy -> Tracy-Widom_2 via Nystrom on the Airy kernel.
K_Airy(x,y)=[Ai(x)Ai'(y)-Ai'(x)Ai(y)]/(x-y);  F_TW2(s)=det(1-K)|_{(s,inf)}.
Ai built by integrating u''=x u from large x with the normalized asymptotic IC (numpy only).
Target: TW2 mean -1.77, std 0.90, skew +0.224, exk +0.093.
"""
import numpy as np
# --- Ai by integration with normalized asymptotic IC ---
x0=8.0; dx=1e-3; xs=np.arange(x0,-7.0,-dx)
z=(2.0/3.0)*x0**1.5
Ai=np.empty(len(xs)); Aip=np.empty(len(xs))
Ai[0]=np.exp(-z)/(2*np.sqrt(np.pi)*x0**0.25)
Aip[0]=-x0**0.25*np.exp(-z)/(2*np.sqrt(np.pi))
u,up=Ai[0],Aip[0]
for i in range(len(xs)-1):
    x=xs[i]; h=-dx
    k1u,k1p=up,x*u; k2u,k2p=up+.5*h*k1p,(x+.5*h)*(u+.5*h*k1u)
    k3u,k3p=up+.5*h*k2p,(x+.5*h)*(u+.5*h*k2u); k4u,k4p=up+h*k3p,(x+h)*(u+h*k3u)
    u+=(h/6)*(k1u+2*k2u+2*k3u+k4u); up+=(h/6)*(k1p+2*k2p+2*k3p+k4p)
    Ai[i+1]=u; Aip[i+1]=up
print(f"Ai(0)={np.interp(0,xs[::-1],Ai[::-1]):.5f} (exact 0.35503);  Ai'(0)={np.interp(0,xs[::-1],Aip[::-1]):.5f} (exact -0.25882)")
def aifun(x): return np.interp(x,xs[::-1],Ai[::-1]), np.interp(x,xs[::-1],Aip[::-1])
def Kairy(x,y):
    ax,apx=aifun(x); ay,apy=aifun(y)
    out=(ax*apy-apx*ay)/(x-y) if abs(x-y)>1e-9 else apx*apx-x*ax*ax
    return out
def Fdet(s,X=8.0,n=90):
    t,w=np.polynomial.legendre.leggauss(n); xq=0.5*(X-s)*t+0.5*(X+s); wq=0.5*(X-s)*w
    K=np.array([[Kairy(xi,xj) for xj in xq] for xi in xq])
    sw=np.sqrt(wq); M=np.eye(n)-sw[:,None]*K*sw[None,:]
    return np.linalg.det(M)
ss=np.linspace(-6,3,46); F=np.array([Fdet(s) for s in ss])
F=np.clip(F,0,1); f=np.gradient(F,ss); f=np.clip(f,0,None); f/=np.trapz(f,ss)
m=np.trapz(ss*f,ss); v=np.trapz((ss-m)**2*f,ss); sd=np.sqrt(v)
sk=np.trapz((ss-m)**3*f,ss)/sd**3; ek=np.trapz((ss-m)**4*f,ss)/v**2-3
print(f"pipeline TW2: mean={m:.3f}(-1.77) std={sd:.3f}(0.90) skew={sk:+.3f}(+0.224) exk={ek:+.3f}(+0.093)")
