import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.special import loggamma
def L_osc(lam,xmax=18.0,h=3e-4):
    lam=np.asarray(lam,dtype=complex); x=xmax
    u=(x*x+lam)**-0.25; up=(1j*np.sqrt(x*x+lam)-0.5*(x/(x*x+lam)))*u; n=int(xmax/h)
    def f(x,u,up): return up,-(x*x+lam)*u
    for _ in range(n):
        a1,b1=f(x,u,up);a2,b2=f(x-h/2,u-h/2*a1,up-h/2*b1);a3,b3=f(x-h/2,u-h/2*a2,up-h/2*b2);a4,b4=f(x-h,u-h*a3,up-h*b3)
        u=u-h/6*(a1+2*a2+2*a3+a4);up=up-h/6*(b1+2*b2+2*b3+b4);x=x-h
    return up/u
def Gi(lam): return np.exp(loggamma(0.75-0.25j*lam)-loggamma(0.25-0.25j*lam))
def Gr(lam): return np.exp(loggamma(0.75-0.25*lam)-loggamma(0.25-0.25*lam))
fig,ax=plt.subplots(1,2,figsize=(11,4.4))
# (a) constant phase e^{3ipi/4}
pts=np.array([0.0,0.7,1.5,2.5,3.5,-1.0,0.89-0.89j,2-1j,1.5+0.5j])
r=np.array([ (L_osc(p)/(2*Gi(p))) for p in pts])
th=np.linspace(0,2*np.pi,200); ax[0].plot(np.cos(th),np.sin(th),color="#ccc",lw=1)
ax[0].plot(r.real,r.imag,"o",color="#264653",ms=8)
ax[0].plot([np.cos(3*np.pi/4)],[np.sin(3*np.pi/4)],"*",color="#e76f51",ms=18,label=r"$e^{3i\pi/4}$")
ax[0].set_aspect("equal"); ax[0].set_xlabel("Re"); ax[0].set_ylabel("Im"); ax[0].legend(fontsize=10)
ax[0].set_title(r"(a) $L(\lambda)/[2\Gamma(3/4-i\lambda/4)/\Gamma(1/4-i\lambda/4)]=e^{3i\pi/4}$"+"\n(9 real & complex $\\lambda$, all coincide)")
# (b) |condition| in complex-lambda plane; root at resonance
re=np.linspace(0.2,2.0,140); im=np.linspace(-1.6,-0.2,120); RE,IM=np.meshgrid(re,im); LAM=RE+1j*IM
C=np.abs(np.exp(3j*np.pi/4)*Gi(LAM)-Gr(LAM))
im0=ax[1].pcolormesh(re,im,np.log10(C+1e-6),cmap="viridis",shading="auto")
ax[1].plot([0.8896],[-0.8896],"*",color="#e76f51",ms=18,label=r"root $0.890-0.890i$")
ax[1].plot([0.890],[-0.890],"o",mfc="none",mec="w",ms=13,label="M1 (FD) resonance")
ax[1].set_xlabel(r"Re $\lambda$"); ax[1].set_ylabel(r"Im $\lambda$"); ax[1].legend(fontsize=8,loc="upper right")
ax[1].set_title(r"(b) $\log_{10}|e^{3i\pi/4}\Gamma(3/4-i\lambda/4)/\Gamma(1/4-i\lambda/4)-\Gamma(3/4-\lambda/4)/\Gamma(1/4-\lambda/4)|$")
fig.colorbar(im0,ax=ax[1],shrink=0.85)
fig.suptitle("Closed-form asymmetric connection of W's skeleton: confining (real Weber) x oscillatory (imaginary Weber) Gamma-ratios",fontsize=10)
fig.subplots_adjust(left=0.08,right=0.98,top=0.85,bottom=0.12,wspace=0.28)
fig.savefig("figures/connection_closed_form.png",dpi=130); print("saved figures/connection_closed_form.png")
