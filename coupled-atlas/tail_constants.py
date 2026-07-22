"""
Validate the analytic tail constants of W against the smooth Fokker-Planck F_beta (non-confounded).
LEFT (late escape, holding instanton):  -ln F(Y*) ~ |Y*|^5/(10 eta^2)   [predict const 1/10 in eta^2-units]
RIGHT (early escape): test exponent & const empirically (barrier vanishes at turning => dynamic, not static).
Key check: eta^2 * (-ln F) vs |Y*|^5 should COLLAPSE across beta onto slope 1/10 (left, large |Y*|).
numpy only.
"""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

betas=[1.0,2.0,4.0,8.0]
fig,ax=plt.subplots(1,2,figsize=(11,4.4))
print("LEFT tail: fit  eta^2*(-lnF) = a*|Y*|^5 + b  ; predict a=1/10=0.100")
cols=plt.cm.viridis(np.linspace(0.1,0.85,4))
for b,cc in zip(betas,cols):
    y,F=np.load(f"fpb_{b:.1f}.npy"); eta=2.0/np.sqrt(b)
    # left tail
    m=(F>2e-4)&(F<0.04)&(y<-2.0)
    Yl=y[m]; L=eta**2*(-np.log(F[m]))
    if m.sum()>=5:
        A=np.vstack([np.abs(Yl)**5, np.ones_like(Yl)]).T
        a,bb=np.linalg.lstsq(A,L,rcond=None)[0]
        # local exponent
        lx=np.log(np.abs(Yl)); ly=np.log(L); ex=np.polyfit(lx,ly,1)[0]
        print(f"  beta={b}: eta^2={eta**2:.2f}  a={a:.4f} (pred 0.100)  b={bb:.2f}  local-exp={ex:.2f} (pred 5)  [{m.sum()} pts, |Y*|<= {-Yl.max():.2f}..{-Yl.min():.2f}]")
        ax[0].plot(np.abs(Yl)**5, L,'o',color=cc,ms=4,label=f'$\\beta$={b:.0f}')
xx=np.linspace(20,np.abs(Yl).max()**5,50); ax[0].plot(xx,xx/10-3.98,'k--',lw=1.4,label=r'naive instanton $|Y^*|^5/10$')
ax[0].set_xlabel(r'$|Y^*|^5$'); ax[0].set_ylabel(r'$\eta^2(-\ln F)$'); ax[0].legend(fontsize=8)
ax[0].set_title('(A) LEFT tail: NO $\\eta^2$-collapse (phase-persistence, anomalous)')

print("\nRIGHT tail: fit eta^2*(-ln(1-F)) vs |Y* - mean|, local exponent & const")
for b,cc in zip(betas,cols):
    y,F=np.load(f"fpb_{b:.1f}.npy"); eta=2.0/np.sqrt(b)
    mu=np.trapz(y*np.gradient(F,y),y)
    m=((1-F)>2e-4)&((1-F)<0.04)&(y>mu)
    Yr=y[m]; R=eta**2*(-np.log(1-F[m]))
    if m.sum()>=5:
        lx=np.log(Yr-mu); ly=np.log(R); ex=np.polyfit(lx,ly,1)[0]
        print(f"  beta={b}: local-exp(Y*-mean)={ex:.2f}  range Y*={Yr.min():.2f}..{Yr.max():.2f} (mean {mu:.2f})  [{m.sum()} pts]")
        ax[1].plot(Yr, R,'o-',color=cc,ms=4,label=f'$\\beta$={b:.0f}')
ax[1].set_xlabel(r'$Y^*$ (early-escape)'); ax[1].set_ylabel(r'$\eta^2(-\ln(1-F))$'); ax[1].legend(fontsize=8)
ax[1].set_title('(B) RIGHT tail: clean $\\eta^2$-COLLAPSE (FW barrier)')
fig.suptitle('Two tails, two mechanisms: RIGHT = confining-barrier (clean FW $1/\\eta^2$); LEFT = oscillatory phase-persistence (anomalous) -- the asymmetry is structural',fontsize=10)
fig.subplots_adjust(left=0.08,right=0.98,bottom=0.12,top=0.88,wspace=0.24)
plt.savefig("figures/tail_constants.png",dpi=115); print("\nsaved figures/tail_constants.png")
