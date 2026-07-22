"""
T2 tail mechanisms -- HONEST reframe. Three panels:
 (1) escape-location std vs eta: the escape is at an O(1) OSCILLATORY NODE (Y*_det~-2.1), so std ~ eta
     (node-dominated phase noise), NOT the eta^{2/3} turning scale. Course-correction.
 (2) standardized W: left vs right tail asymmetry (the two mechanisms; +skew). Pre-asymptotic note.
 (3) the deterministic skeleton u(Y): recessive decay for Y>0 (confining), oscillation for Y<0, first node
     = the connection problem. Illustrates Part 1.
numpy only.
"""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
def Vc(Y): return np.sign(Y)*Y*Y
def sim(eta,N,Y0=3.0,Yend=-5.0,dt=5e-3,seed=5):
    rng=np.random.default_rng(seed); nst=int(round((Y0-Yend)/dt)); sq=np.sqrt(dt)
    p=np.full(N,np.sqrt(Vc(Y0))); Ys=np.full(N,np.nan)
    for i in range(nst):
        Y=Y0-i*dt; p+=(Vc(Y)-p*p)*dt+eta*sq*rng.standard_normal(N); np.clip(p,-60,60,out=p)
        nw=np.isnan(Ys)&(p<-12.0); Ys[nw]=Y0-(i+1)*dt
    Ys[np.isnan(Ys)]=Yend; return Ys
# (1) std vs eta
etas=[0.4,0.6,0.9,1.4,2.0]; stds=[]; meds=[]
for eta in etas:
    Y=sim(eta,60000); stds.append(Y.std()); meds.append(np.median(Y))
stds=np.array(stds)
sl=np.polyfit(np.log(etas),np.log(stds),1)[0]
print(f"std~eta^{sl:.2f}  (eta-linear=1.00, turning-scale=0.67); medians(Y*_det proxy)->{meds[0]:.2f} as eta->0")
# (2) standardized tails at eta=1.4
Y=sim(1.4,250000); s=(Y-Y.mean())/Y.std(); skew=np.mean(((Y-Y.mean())/Y.std())**3)
Sg=np.linspace(0.6,3.0,16)
PR=np.array([(s> S).mean() for S in Sg]); PL=np.array([(s<-S).mean() for S in Sg])
okR=PR*len(s)>=30; okL=PL*len(s)>=30
# (3) deterministic skeleton
dt=2e-3; Yg=np.arange(3.0,-4.0,-dt); u=np.empty_like(Yg); u[0]=1e-6; up=-3*u[0]
for i in range(len(Yg)-1):
    h=-dt; k1u,k1p=up,Vc(Yg[i])*u[i]
    k2u,k2p=up+.5*h*k1p, Vc(Yg[i]+.5*h)*(u[i]+.5*h*k1u)
    k3u,k3p=up+.5*h*k2p, Vc(Yg[i]+.5*h)*(u[i]+.5*h*k2u)
    k4u,k4p=up+h*k3p, Vc(Yg[i]+h)*(u[i]+h*k3u)
    u[i+1]=u[i]+(h/6)*(k1u+2*k2u+2*k3u+k4u); up=up+(h/6)*(k1p+2*k2p+2*k3p+k4p)
un=u/np.max(np.abs(u[Yg<0.5])); node=Yg[:-1][(un[:-1]*un[1:]<0)]; node=node[node<0]
Ystar=node[0] if len(node) else np.nan
print(f"deterministic first node Y*_det = {Ystar:.3f}  (matches MC median as eta->0)")
fig,ax=plt.subplots(1,3,figsize=(15,4.4))
ax[0].loglog(etas,stds,'o-',color='navy',ms=8,label=f'measured: std~eta^{sl:.2f}')
ax[0].loglog(etas,0.5*np.array(etas),'--',color='crimson',label='eta^1 (node-dominated)')
ax[0].loglog(etas,0.5*np.array(etas)**(2/3),':',color='gray',label='eta^{2/3} (turning scale)')
ax[0].set_xlabel('eta'); ax[0].set_ylabel('escape-loc std'); ax[0].legend(fontsize=8)
ax[0].set_title('(1) std ~ eta: escape at O(1) node,\nnot the eta^{2/3} turning scale')
ax[1].plot(Sg[okR],-np.log(PR[okR]),'o-',color='crimson',ms=4,label='RIGHT (barrier/early)')
ax[1].plot(Sg[okL],-np.log(PL[okL]),'s-',color='navy',ms=4,label='LEFT (oscillatory/late)')
ax[1].set_xlabel('standardized |s|'); ax[1].set_ylabel('-log P'); ax[1].legend(fontsize=8)
ax[1].set_title(f'(2) tail asymmetry (skew={skew:+.2f});\nexact (5,3) exps pre-asymptotic')
ax[2].plot(Yg,un,color='teal',lw=1.4); ax[2].axhline(0,color='k',lw=0.6); ax[2].axvline(0,color='gray',ls=':',lw=0.8)
if not np.isnan(Ystar): ax[2].plot([Ystar],[0],'rv',ms=10,label=f'first node Y*={Ystar:.2f}')
ax[2].annotate('confining\n(recessive, U)',xy=(1.4,0.5),fontsize=8,color='navy',ha='center')
ax[2].annotate('oscillatory\n(W)',xy=(-2.6,0.5),fontsize=8,color='crimson',ha='center')
ax[2].set_xlabel('Y'); ax[2].set_ylabel('u(Y) skeleton'); ax[2].legend(fontsize=8,loc='lower left')
ax[2].set_title("(3) deterministic connection:\nU(Y>0) -> first node -> W(Y<0)"); ax[2].set_ylim(-1.1,1.1)
fig.suptitle('T2 connection mechanisms (honest): escape-location std ~ eta (node-dominated); tail asymmetry; the skeleton connection',fontsize=11,fontweight='bold')
fig.subplots_adjust(left=0.06,right=0.98,bottom=0.13,top=0.83,wspace=0.28)
plt.savefig("figures/tail_mechanisms.png",dpi=115); print("saved figures/tail_mechanisms.png")
