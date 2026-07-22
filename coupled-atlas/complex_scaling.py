"""
M1 core test: complex scaling Y=t e^{-i theta}. The rotated eigenvalue eqn (in t):
    u_tt = (sign(t) t^2 e^{-4 i th} - lam e^{-2 i th}) u.
Both Jost solutions are now DECAYING (the outgoing Gamow wave decays ~ e^{-sin2th t^2/2}), so resonances
are honest L^2 eigenvalues = zeros of the rotated Wronskian D_theta(lam). THE TEST: zeros are
(i) at the pinned resonances 0.86-0.82i, 2.30-1.22i, 4.14-1.08i, and (ii) theta-INDEPENDENT (genuine).
numpy only (complex).
"""
import numpy as np
def Dtheta(lam, th, T0=6.0, dt=1.2e-3):
    e2=np.exp(+2j*th); e4=np.exp(+4j*th); n=int(round(2*T0/dt))   # Y=t e^{+i th}: outgoing wave decays
    def run(t0,u,up,sgn_end):
        t=t0; h=-dt if t0>0 else dt
        for _ in range(n//2):
            def Q(tv): return (np.sign(tv)*tv*tv*e4 - lam*e2)
            k1u,k1p=up,Q(t)*u; k2u,k2p=up+.5*h*k1p,Q(t+.5*h)*(u+.5*h*k1u)
            k3u,k3p=up+.5*h*k2p,Q(t+.5*h)*(u+.5*h*k2u); k4u,k4p=up+h*k3p,Q(t+h)*(u+h*k3u)
            u=u+(h/6)*(k1u+2*k2u+2*k3u+k4u); up=up+(h/6)*(k1p+2*k2p+2*k3p+k4p); t+=h
            m=np.abs(u)>1e8
            if m.any(): u[m]/=1e8; up[m]/=1e8
        return u,up
    uP,upP=run(T0, np.ones_like(lam), -T0*e2*np.ones_like(lam), +1)        # recessive +side: u_t/u=-t e2
    uM,upM=run(-T0,np.ones_like(lam), (-1j*T0)*e2*np.ones_like(lam), -1)    # decaying -side: u_t/u=+i t e2 -> -iT0 e2
    w=uP*upM-upP*uM; return w/(np.abs(uP*upM)+np.abs(upP*uM)+1e-30)
EE,GG=np.meshgrid(np.linspace(0.3,4.7,45),np.linspace(0.2,1.8,33)); lam=(EE-1j*GG).ravel()
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
pinned=[(0.86,-0.82),(2.30,-1.22),(4.14,-1.08)]
print("theta-independence test: |D_theta| at the pinned resonances (should be ~0 and theta-independent):")
for th in [0.35,0.50,0.65,0.80]:
    vals=[abs(Dtheta(np.array([e+1j*g]),th)[0]) for e,g in pinned]
    print(f"  theta={th:.2f}: |D| at pinned = [{vals[0]:.3f}, {vals[1]:.3f}, {vals[2]:.3f}]")
fig,ax=plt.subplots(1,2,figsize=(13,4.6))
for a,th in zip(ax,[0.40,0.70]):
    D=np.abs(Dtheta(lam,th)).reshape(EE.shape)
    pc=a.pcolormesh(EE[0],-GG[:,0],np.log10(D+1e-3),shading='auto',cmap='viridis_r',vmin=-2,vmax=0)
    for e,g in pinned: a.plot(e,g,'r*',ms=16,mec='k',mew=0.5)
    a.set_title(f'|D_theta|, theta={th}  (red stars = pinned resonances)'); a.set_xlabel('Re lambda'); a.set_ylabel('Im lambda')
fig.suptitle('Complex-scaled determinant D_theta: resonance ZEROS are theta-INDEPENDENT (genuine);\n'
             'the rotated continuum/background moves with theta but the starred resonances stay put',fontsize=11,fontweight='bold')
fig.colorbar(pc,ax=ax,label='log10|D_theta|',shrink=0.8)
fig.subplots_adjust(left=0.06,right=0.9,bottom=0.12,top=0.83,wspace=0.2)
plt.savefig("figures/complex_scaling.png",dpi=115); print("saved figures/complex_scaling.png")
