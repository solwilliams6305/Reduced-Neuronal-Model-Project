"""
Accurate resonances of H=-d2/dY2+sign(Y)Y2 via the outgoing-matching Wronskian (direct integration),
+ a clean landscape figure. These supersede the crude-scan lambda0=1.64-1.49i.
"""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
def integ(lam,Y0s,u,up,Yend,dt=1.0e-3):
    Y=Y0s; n=int(round(abs(Y0s-Yend)/dt)); h=-dt if Y0s>Yend else dt
    for _ in range(n):
        QY=np.sign(Y)*Y*Y-lam; QYh=np.sign(Y+.5*h)*(Y+.5*h)**2-lam; QYf=np.sign(Y+h)*(Y+h)**2-lam
        k1u,k1p=up,QY*u; k2u,k2p=up+.5*h*k1p,QYh*(u+.5*h*k1u)
        k3u,k3p=up+.5*h*k2p,QYh*(u+.5*h*k2u); k4u,k4p=up+h*k3p,QYf*(u+h*k3u)
        u=u+(h/6)*(k1u+2*k2u+2*k3u+k4u); up=up+(h/6)*(k1p+2*k2p+2*k3p+k4p); Y+=h
        m=np.abs(u)>1e8
        if m.any(): u[m]/=1e8; up[m]/=1e8
    return u,up
def Dn(lam,Y0=6.0):   # normalized outgoing-matching Wronskian
    uR,upR=integ(lam,Y0,np.ones_like(lam),-np.sqrt(Y0*Y0-lam),0.0)
    uL,upL=integ(lam,-Y0,np.ones_like(lam),(1j)*(-Y0)*np.ones_like(lam),0.0)  # outgoing branch
    w=uR*upL-upR*uL; return w/(np.abs(uR*upL)+np.abs(upR*uL)+1e-30)
EE,GG=np.meshgrid(np.linspace(0.3,4.7,56),np.linspace(0.1,2.2,33)); lam=(EE-1j*GG).ravel()
D=np.abs(Dn(lam)).reshape(EE.shape)
# extract local minima (resonances)
res=[]
for i in range(1,D.shape[0]-1):
    for j in range(1,D.shape[1]-1):
        if D[i,j]<0.07 and D[i,j]==D[i-1:i+2,j-1:j+2].min():
            res.append((EE[i,j],-GG[i,j],D[i,j]))
res.sort(key=lambda t:t[0])
print("resonances (outgoing-matching Wronskian, |D|<0.07 local minima):")
for e,g,d in res: print(f"   lam = {e:+.2f}{g:+.2f}i   |D|={d:.3f}")
fig,ax=plt.subplots(figsize=(7.4,4.8))
pc=ax.pcolormesh(EE[0],-GG[:,0],np.log10(D+1e-3),shading='auto',cmap='viridis_r')
for e,g,d in res: ax.plot(e,g,'r*',ms=15)
ax.axvline(1,ls=':',color='w',lw=1); ax.axvline(3,ls=':',color='w',lw=1)
ax.text(1.03,-0.3,'bare 1',color='w',fontsize=7); ax.text(3.03,-0.3,'bare 3',color='w',fontsize=7)
ax.set_xlabel('Re lambda'); ax.set_ylabel('Im lambda (=-width/2)')
ax.set_title('TRUE resonances (outgoing-matching Wronskian, direct integration)\n'
             'complex => non-self-adjoint; string Re~0.9,2.3,4.1, width ~2')
fig.colorbar(pc,label='log10 |D_normalized|')
fig.subplots_adjust(left=0.1,right=1.0,bottom=0.13,top=0.85)
plt.savefig("figures/resonances_true.png",dpi=115); print("saved figures/resonances_true.png")
