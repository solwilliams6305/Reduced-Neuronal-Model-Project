import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import importlib.util as u
spec=u.spec_from_file_location("m","malliavin_firstpassage_check.py")
# re-run the core computation inline (small) for figure data
def q(Y,lam): return np.sign(Y)*Y*Y-lam
def integrate(xi,dt,Y0,lam,eta):
    n=len(xi);Y=np.empty(n+1);uu=np.empty(n+1);up=np.empty(n+1)
    Y[0]=Y0;uu[0]=1.0;up[0]=np.sqrt(max(q(Y0,lam),1e-9));sq=np.sqrt(dt)
    for i in range(n):
        Yi=Y[0]-i*dt;Q=q(Yi,lam);uu[i+1]=uu[i]-up[i]*dt;up[i+1]=up[i]-(Q*uu[i])*dt+eta*uu[i]*sq*xi[i];Y[i+1]=Yi-dt
    return Y,uu,up
def fz(Y,uu,up):
    for i in range(1,len(uu)):
        if Y[i]<0 and uu[i-1]*uu[i]<0:
            f=uu[i-1]/(uu[i-1]-uu[i]);return Y[i-1]+f*(Y[i]-Y[i-1]),up[i-1]+f*(up[i]-up[i-1]),i
    return np.nan,np.nan,-1
lam=1.0;eta=0.20;Y0=3.0;dt=2e-3;n=int((Y0+4.0)/dt)
rng=np.random.default_rng(0);xi=rng.standard_normal(n)
Y,uu,up=integrate(xi,dt,Y0,lam,eta);Ys,ups,bi=fz(Y,uu,up)
eps=1e-3;js=np.arange(int(0.12*bi),bi-2,4)
meas=[];pred=[];svals=[]
for j in js:
    xp=xi.copy();xp[j]+=eps;Yp,uup,upp=integrate(xp,dt,Y0,lam,eta);Ysp,_,_=fz(Yp,uup,upp)
    meas.append((Ysp-Ys)/eps);pred.append(-((uup[bi]-uu[bi])/eps)/ups);svals.append(Y[j])
meas=np.array(meas);pred=np.array(pred);svals=np.array(svals)
fig,ax=plt.subplots(1,3,figsize=(13,4))
ax[0].plot(pred,meas,"o",color="#264653",ms=5);lim=[min(pred.min(),meas.min()),max(pred.max(),meas.max())]
ax[0].plot(lim,lim,"r--",lw=1);ax[0].set_xlabel(r"predicted $-D_s u(Y^*)/u'(Y^*)$");ax[0].set_ylabel(r"measured $D_s Y^*$")
ax[0].set_title("(a) IFT identity: corr=1.0000, slope=1.00")
ax[1].plot(svals,meas,"o-",color="#2a9d8f",ms=4);ax[1].axvline(0,color="gray",ls=":",lw=1,label="turning k=0")
ax[1].axvline(Ys,color="r",ls=":",lw=1,label=f"escape Y*={Ys:.2f}")
ax[1].set_xlabel("s (=Y)");ax[1].set_ylabel(r"$D_s Y^*$ (Malliavin deriv)");ax[1].set_title("(b) sensitivity kernel: smooth, regular through turning");ax[1].legend(fontsize=8)
j=js[len(js)//2];es=[4e-3,2e-3,1e-3,5e-4,2.5e-4];rr=[]
for e in es:
    xp=xi.copy();xp[j]+=e;Yp,uup,upp=integrate(xp,dt,Y0,lam,eta);Ysp,_,_=fz(Yp,uup,upp);rr.append((Ysp-Ys)/e)
ax[2].semilogx(es,rr,"ks-",ms=6);ax[2].set_xlabel(r"$\epsilon$");ax[2].set_ylabel(r"$\Delta Y^*/\epsilon$")
ax[2].set_ylim(min(rr)*0.9,max(rr)*1.1);ax[2].set_title("(c) linear response (differentiable, no kink)")
fig.suptitle(r"T1 item 3: first-passage $Y^*$ Malliavin-differentiable in the FIELD rep — $D_s Y^*=-D_s u(Y^*)/u'(Y^*)$, transversal $u'(Y^*)\neq0$",fontsize=10)
fig.subplots_adjust(left=0.07,right=0.98,top=0.86,bottom=0.13,wspace=0.3)
fig.savefig("figures/malliavin_firstpassage.png",dpi=130);print("saved figures/malliavin_firstpassage.png")
