"""
Resonance determinant of the asymmetric operator H = -d^2/dY^2 + sign(Y)Y^2.
The confining (+Y^2) quasi-bound states leak out the inverted (-Y^2) side => COMPLEX resonances
(non-self-adjoint structure). Compute the incoming/outgoing amplitudes of the recessive-at-+inf
solution at -Y0; resonance = incoming amplitude vanishes (purely outgoing, Gamow state).
Demonstrates the complex spectral data that obstructs a real-spectrum (standard) determinant.
numpy only (complex).
"""
import numpy as np
def integrate(lam,Y0=6.0,dt=1.5e-3):
    """lam: complex array. Vectorized; loop only over time steps. Returns A,B amplitude arrays."""
    n=int(round(2*Y0/dt)); Y=Y0
    u=np.ones_like(lam); up=-np.sqrt(Y0*Y0-lam)   # recessive IC at +Y0
    for i in range(n):
        h=-dt
        QY=np.sign(Y)*Y*Y-lam; QYh=np.sign(Y+.5*h)*(Y+.5*h)**2-lam; QYf=np.sign(Y+h)*(Y+h)**2-lam
        k1u,k1p=up, QY*u
        k2u,k2p=up+.5*h*k1p, QYh*(u+.5*h*k1u)
        k3u,k3p=up+.5*h*k2p, QYh*(u+.5*h*k2u)
        k4u,k4p=up+h*k3p, QYf*(u+h*k3u)
        u=u+(h/6)*(k1u+2*k2u+2*k3u+k4u); up=up+(h/6)*(k1p+2*k2p+2*k3p+k4p); Y+=h
        m=np.abs(u)>1e10
        if m.any(): u[m]/=1e10; up[m]/=1e10
    k=np.sqrt(Y0*Y0+lam)                 # wavenumber at -Y0 (V=-Y^2-lam)
    A=(u-1j*up/k)/2; B=(u+1j*up/k)/2     # u = A e^{+iS} + B e^{-iS}
    return A,B
Es=np.linspace(0.2,4.2,41); Gs=np.linspace(0.02,2.2,23)
EE,GG=np.meshgrid(Es,Gs); lam=(EE-1j*GG).ravel()
A,B=integrate(lam)
ratio=np.minimum(np.abs(A),np.abs(B))/(np.abs(A)+np.abs(B)+1e-30)
order=np.argsort(ratio)
print("lowest normalized min(|A|,|B|)/sum (candidate resonances, incoming amplitude ~0):")
for idx in order[:6]:
    which='B' if np.abs(B[idx])<np.abs(A[idx]) else 'A'
    print(f"  lam = {EE.ravel()[idx]:.2f} - {GG.ravel()[idx]:.2f}i   ratio={ratio[idx]:.3f}  (~0: {which})")
E0,G0=EE.ravel()[order[0]],GG.ravel()[order[0]]
Ef=np.linspace(E0-0.18,E0+0.18,19); Gf=np.linspace(max(0.01,G0-0.18),G0+0.18,19)
EEf,GGf=np.meshgrid(Ef,Gf); lamf=(EEf-1j*GGf).ravel(); Af,Bf=integrate(lamf)
rf=np.minimum(np.abs(Af),np.abs(Bf))/(np.abs(Af)+np.abs(Bf)+1e-30); o=np.argsort(rf)
l0E,l0G=EEf.ravel()[o[0]],GGf.ravel()[o[0]]
print(f"\nrefined lowest resonance: lam_0 = {l0E:.3f} - {l0G:.3f}i")
print("  Im<0 (finite width) => COMPLEX resonance => non-self-adjoint spectral data (vs real TW/Weber).")
print("  nearest confining bare level lam=1 (half-oscillator); the shift+width = inverted-side leakage.")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
R=ratio.reshape(GG.shape)
fig,ax=plt.subplots(figsize=(7.2,4.8))
pc=ax.pcolormesh(Es,-Gs,np.log10(R+1e-6),shading='auto',cmap='viridis_r')
# mark the resonance string (local minima along Re axis)
for idx in order[:12]:
    e,g=EE.ravel()[idx],GG.ravel()[idx]
    if ratio[idx]<0.02: ax.plot(e,-g,'rx',ms=9,mew=2)
ax.plot([l0E],[-l0G],'r*',ms=16,label=f'lowest resonance {l0E:.2f}-{l0G:.2f}i')
ax.axvline(1,ls=':',color='w',lw=1); ax.axvline(3,ls=':',color='w',lw=1)
ax.text(1.02,-0.15,'bare lvl 1',color='w',fontsize=7); ax.text(3.02,-0.15,'bare lvl 3',color='w',fontsize=7)
ax.set_xlabel('Re lambda (resonance energy)'); ax.set_ylabel('Im lambda (= -width/2)')
ax.set_title('Complex resonances of H=-d2/dY2+sign(Y)Y2\n(confining levels leak out the inverted side => Im<0)')
fig.colorbar(pc,label='log10 normalized |incoming amp|'); ax.legend(fontsize=8,loc='lower right')
fig.subplots_adjust(left=0.1,right=1.0,bottom=0.13,top=0.86)
plt.savefig("figures/resonance_landscape.png",dpi=115); print("saved figures/resonance_landscape.png")
