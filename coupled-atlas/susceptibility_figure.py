import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
# recompute minimal data
def solve_field(lam,a,Y0=3.0,Yend=-4.0,dt=1e-4,u0=1.0,v0=None):
    def Va(Y): return (Y*Y if Y>=0 else -(1+a)*Y*Y)
    n=int(round((Y0-Yend)/dt));Y=np.empty(n+1);u=np.empty(n+1);v=np.empty(n+1)
    Y[0]=Y0;u[0]=u0;v[0]=(np.sqrt(max(Va(Y0)-lam,1e-9)) if v0 is None else v0)
    for i in range(n):
        Yi=Y[0]-i*dt;q=Va(Yi)-lam;u[i+1]=u[i]+v[i]*dt;v[i+1]=v[i]+q*u[i]*dt;Y[i+1]=Yi-dt
    return Y,u,v
def fz(Y,u,v):
    for i in range(1,len(u)):
        if Y[i]<0 and u[i-1]*u[i]<0:
            f=u[i-1]/(u[i-1]-u[i]);return Y[i-1]+f*(Y[i]-Y[i-1])
    return np.nan
lam=1.0
avals=np.linspace(-0.08,0.08,9); ys=[fz(*solve_field(lam,a)) for a in avals]
dpred=0.1205; y0=fz(*solve_field(lam,0.0))
lam=1.0
def Vq_a(Y,a): return np.where(Y>=0,Y*Y,-(1+a)*Y*Y)
def esc(a,eta,N,Y0=3.0,Yend=-5.0,dt=1.5e-3,seed=1,thr=25.0):
    rng=np.random.default_rng(seed);n=int(round((Y0-Yend)/dt));sq=np.sqrt(dt)
    p=np.full(N,np.sqrt(max(Y0*Y0-lam,1e-9)));Ys=np.full(N,np.nan)
    for i in range(n):
        Y=Y0-i*dt;V=(Y*Y if Y>=0 else -(1+a)*Y*Y)-lam
        p+=(V-p*p)*dt+eta*sq*rng.standard_normal(N);np.clip(p,-60,60,out=p)
        nw=np.isnan(Ys)&(p<-thr);Ys[nw]=Y0-(i+1)*dt
    return Ys[~np.isnan(Ys)]
def sk(x):d=x-x.mean();return np.mean(d**3)/np.mean(d**2)**1.5
aa=np.linspace(-0.25,0.25,7); sks=[sk(esc(a,np.sqrt(2.0),80000,seed=7)) for a in aa]
fig,ax=plt.subplots(1,2,figsize=(10,4.2))
ax[0].plot(avals,ys,"o",color="#264653",ms=6,label="Y*(a) (numeric)")
ax[0].plot(avals,y0+dpred*avals,"r--",lw=1.3,label=f"derived slope +{dpred:.3f}")
ax[0].set_xlabel("symmetry-breaking a");ax[0].set_ylabel("escape location Y*")
ax[0].set_title("(a) location susceptibility dY*/da: derived=+0.121 vs FD=+0.124")
ax[0].legend(fontsize=8)
ax[1].plot(aa,sks,"s-",color="#7a3b8f",ms=6)
sl=np.polyfit(aa,sks,1)[0]
ax[1].set_xlabel("symmetry-breaking a");ax[1].set_ylabel("escape skew")
ax[1].set_title(f"(b) skew susceptibility d(skew)/da=+{sl:.2f} (E meas +0.36)\ndeeper osc. well -> heavier left tail -> +skew")
fig.suptitle("Cusp within-class symmetry-breaking susceptibility: escape location DERIVED (Green's fn), skew mechanism",fontsize=10.5)
fig.subplots_adjust(left=0.08,right=0.97,top=0.85,bottom=0.13,wspace=0.28)
fig.savefig("figures/susceptibility.png",dpi=130);print("saved figures/susceptibility.png")
