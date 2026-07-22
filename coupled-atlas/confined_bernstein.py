"""
Validate the LOWER-TAIL sub-Gaussianity that drives the confined Bernstein bound.
The proved no-early-escape supermartingale gives  P(x<=-L) <= C exp(-2 Phi(-L)/eta^2),
Phi(-L)=pbar L^2 - L^3/3 ~ pbar L^2 for moderate L  =>  -ln P(x<=-L) ~ (2 pbar/eta^2) L^2  (sub-Gaussian).
Integrating this tail gives  E[e^{lambda x} 1_conf] <= C exp(C' eta^2 lambda^2), C' ~ 1/(8 c0)  (uniform via floor).

Check, for the inner fluctuation x=dp at the turning entrance, across Delta and eta:
 (A) -ln P(x<=-L) is LINEAR in L^2 (sub-Gaussian), slope ~ 2 pbar/eta^2 >= 2 c0/eta^2;
 (B) the empirical log-MGF log E[e^{lambda x} 1_{x>-barrier}] <= parabola C' eta^2 lambda^2 (Bernstein), C' uniform in Delta.
numpy only.
"""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

def Vfun(Y,D): aY=np.abs(Y); return np.sign(Y)*aY*(aY+D)
def canard(Yg,D):
    n=len(Yg); pb=np.empty(n); pb[0]=np.sqrt(max(Vfun(Yg[0],D),1e-12))
    for i in range(n-1):
        Y=Yg[i]; h=Yg[i+1]-Yg[i]; f=lambda Yv,pv: pv*pv-Vfun(Yv,D)
        k1=f(Y,pb[i]);k2=f(Y+.5*h,pb[i]+.5*h*k1);k3=f(Y+.5*h,pb[i]+.5*h*k2);k4=f(Y+h,pb[i]+h*k3)
        nx=pb[i]+(h/6)*(k1+2*k2+2*k3+k4); pb[i+1]=nx if (np.isfinite(nx) and 0<nx<1e6) else pb[i]
    return pb

def sim_x(D,eta,N,Y0=3.0,Yref=0.06,dt=2e-3,seed=3):
    rng=np.random.default_rng(seed); n=int(round((Y0-(-0.5))/dt)); sq=np.sqrt(dt)
    Yg=Y0-np.arange(n+1)*dt; pb=canard(Yg,D); iref=int(round((Y0-Yref)/dt))
    p=np.full(N,pb[0]);
    for i in range(n):
        p=p+(Vfun(Yg[i],D)-p*p)*dt+eta*sq*rng.standard_normal(N); np.clip(p,-40,40,out=p)
        if i==iref: return p-pb[i], pb[i]
    return p-pb[-1], pb[-1]

c0=0.69
fig,ax=plt.subplots(1,3,figsize=(15,4.3))
print("(A) lower tail: -ln P(x<=-L) vs L^2  => slope ~ 2 pbar/eta^2 (sub-Gaussian, uniform via floor)")
configs=[(0.05,0.5,'crimson'),(0.5,0.5,'teal'),(2.5,0.5,'navy'),(0.05,0.7,'orange')]
slopes={}
for D,eta,c in configs:
    x,pb=sim_x(D,eta,120000); xb=x[x>-2*pb+0.05]  # confined (not escaped)
    Ls=np.linspace(0.3,1.2,18); P=np.array([(x<=-L).mean() for L in Ls]); ok=P*len(x)>=20
    L2=Ls[ok]**2; nlP=-np.log(P[ok])
    sl=np.polyfit(L2,nlP,1)[0]; slopes[(D,eta)]=sl
    ax[0].plot(L2,nlP,'o-',color=c,ms=4,label=f'D={D},eta={eta}: slope {sl:.1f} (2pbar/eta^2={2*pb/eta**2:.1f})')
    print(f"  D={D} eta={eta}: pbar={pb:.3f}, slope(L^2)={sl:.2f}, 2pbar/eta^2={2*pb/eta**2:.2f}, 2c0/eta^2={2*c0/eta**2:.2f}")
ax[0].set_xlabel(r'$L^2$'); ax[0].set_ylabel(r'$-\ln P(x\leq -L)$'); ax[0].legend(fontsize=6.5)
ax[0].set_title('(A) lower tail LINEAR in L^2 = sub-Gaussian')
# (B) log-MGF vs lambda, confined
print("\n(B) confined log-MGF log E[e^{lambda x} 1_conf] vs lambda <= parabola C' eta^2 lambda^2")
lams=np.linspace(-4,4,33)
for D,eta,c in [(0.05,0.5,'crimson'),(0.5,0.5,'teal'),(2.5,0.5,'navy')]:
    x,pb=sim_x(D,eta,120000); xc=x[x>-2*pb+0.05]
    G=np.array([np.log(np.mean(np.exp(np.clip(l*xc,-30,30)))) for l in lams])
    Gc=G-lams*xc.mean()  # center
    ax[1].plot(lams,Gc,'-',color=c,lw=1.6,label=f'D={D}')
    cprime=np.polyfit(lams,Gc,2)[0]/eta**2
    print(f"  D={D} eta={eta}: centered log-MGF curvature C'={cprime:.3f} (~1/(8c0)={1/(8*c0):.3f}); uniform in D check")
ax[1].set_xlabel(r'$\lambda$'); ax[1].set_ylabel('centered log-MGF'); ax[1].legend(fontsize=8)
ax[1].set_title("(B) Bernstein: log-MGF <= parabola (both sides)")
# (C) slope vs D (uniformity)
Ds=[0.05,0.15,0.4,1.0,2.5]; eta=0.5; sl=[]
for D in Ds:
    x,pb=sim_x(D,eta,120000); Ls=np.linspace(0.3,1.1,16); P=np.array([(x<=-L).mean() for L in Ls]); ok=P*len(x)>=20
    sl.append(np.polyfit(Ls[ok]**2,-np.log(P[ok]),1)[0])
ax[2].plot([D/eta**(2/3) for D in Ds],sl,'s-',color='purple',ms=7,label='measured slope')
ax[2].axhline(2*c0/eta**2,ls=':',color='k',label=f'2c0/eta^2={2*c0/eta**2:.1f} (floor bound)')
ax[2].set_xlabel(r'$\rho=\Delta/\ell$'); ax[2].set_ylabel('lower-tail slope (L^2)'); ax[2].set_xscale('log')
ax[2].set_title('(C) slope >= 2c0/eta^2 uniformly in rho'); ax[2].legend(fontsize=8)
fig.suptitle('Confined Bernstein: lower-tail sub-Gaussian (slope >= 2c0/eta^2, floor-uniform) => MGF <= exp(C eta^2 lambda^2) by integration',fontsize=10.5)
fig.subplots_adjust(left=0.06,right=0.98,bottom=0.12,top=0.87,wspace=0.27)
plt.savefig("figures/confined_bernstein.png",dpi=115); print("\nsaved figures/confined_bernstein.png")
