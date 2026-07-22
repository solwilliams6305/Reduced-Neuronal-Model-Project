"""
T2 tail-constant validation (non-confounded, mechanism-driven).
Prediction from (inner scale std ~ eta^{2/3}) x (FW exponents: right 3, left 5):
  RIGHT tail  P(dY* > d) ~ exp(-c_R d^3 / eta^2)        => collapses under eta^2      (barrier, clean)
  LEFT  tail  P(dY* <-d) ~ exp(-c_L d^5 / eta^{10/3})   => collapses under eta^{10/3} (oscillatory, ANOMALOUS)
Cross-check: right should NOT collapse under eta^{10/3}, left should NOT collapse under eta^2.
Cusp V=sign(Y)Y^2. numpy only.
"""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
def Vc(Y): return np.sign(Y)*Y*Y
def sim(eta,N,Y0=3.0,Yend=-5.0,dt=6e-3,seed=5):
    rng=np.random.default_rng(seed); nst=int(round((Y0-Yend)/dt)); sq=np.sqrt(dt)
    p=np.full(N,np.sqrt(Vc(Y0))); Ys=np.full(N,np.nan)
    for i in range(nst):
        Y=Y0-i*dt; p+=(Vc(Y)-p*p)*dt+eta*sq*rng.standard_normal(N); np.clip(p,-60,60,out=p)
        nw=np.isnan(Ys)&(p<-12.0); Ys[nw]=Y0-(i+1)*dt
    Ys[np.isnan(Ys)]=Yend; return Ys
etas=[1.0,1.4,2.0]; cols=['navy','teal','crimson']; N=200000
data={}
for eta in etas:
    Y=sim(eta,N); data[eta]=Y-np.median(Y)
    print(f"eta={eta}: median={np.median(Y):.3f} std={Y.std():.3f} (eta^(2/3)={eta**(2/3):.3f})")
def tail(x, side, dgrid):
    if side=='R': P=np.array([(x> d).mean() for d in dgrid])
    else:         P=np.array([(x<-d).mean() for d in dgrid])
    ok=P*len(x)>=40; return dgrid[ok], -np.log(P[ok])
fig,ax=plt.subplots(2,2,figsize=(12,9))
for eta,c in zip(etas,cols):
    x=data[eta]
    dR=np.linspace(0.3*eta**(2/3),2.2*eta**(2/3),16); dL=np.linspace(0.3*eta**(2/3),2.2*eta**(2/3),16)
    for j,(side,d) in enumerate([('R',dR),('L',dL)]):
        dd,nlp=tail(x,side,d)
        ax[0,j].plot(dd,nlp*eta**2,'o-',color=c,ms=4,label=f'eta={eta}')
        ax[1,j].plot(dd,nlp*eta**(10/3),'s-',color=c,ms=4,label=f'eta={eta}')
        if eta==1.4 and len(dd)>6:  # report local exponent
            sl=np.polyfit(np.log(dd[-7:]),np.log(nlp[-7:]),1)[0]
            print(f"  {side} tail local exponent (eta=1.4) ~ {sl:.2f}  (predicted {'3' if side=='R' else '5'})")
titles=[["RIGHT tail x eta^2  (predict COLLAPSE)","LEFT tail x eta^2  (predict NO collapse)"],
        ["RIGHT tail x eta^{10/3}  (predict NO collapse)","LEFT tail x eta^{10/3}  (predict COLLAPSE)"]]
for a in range(2):
    for b in range(2):
        ax[a,b].set_title(titles[a][b],fontsize=10); ax[a,b].legend(fontsize=8); ax[a,b].grid(alpha=0.3)
        ax[a,b].set_xlabel('dY*  (physical fluctuation)'); ax[a,b].set_ylabel('rescaled  -logP')
fig.suptitle('T2 tail scalings: right = eta^2 (barrier, clean FW) vs left = eta^{10/3} (oscillatory, anomalous)',fontsize=12,fontweight='bold')
fig.subplots_adjust(left=0.07,right=0.98,bottom=0.07,top=0.92,hspace=0.28,wspace=0.2)
plt.savefig("figures/tail_scaling.png",dpi=110); print("saved figures/tail_scaling.png")
