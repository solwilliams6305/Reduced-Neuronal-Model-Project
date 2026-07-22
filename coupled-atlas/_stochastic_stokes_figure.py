import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import stochastic_stokes_o_eta2 as m
lam0=m.newton_root(0.89-0.89j); R=m.assemble(lam0); O2=R['Omega2']; pv=R['var_l1']
truevar=1.7835
fig,ax=plt.subplots(1,2,figsize=(13.5,5.3))
# ---- Panel A: lambda-plane
a=ax[0]
a.axhline(0,color='0.8',lw=.7); a.axvline(0,color='0.8',lw=.7)
for ang,txt,c in [(-45,'det. Borel phase -45 deg (lambda0)','tab:blue'),(-54,'measured stochastic ~ -54 deg','tab:red')]:
    r=1.5; a.plot([0,r*np.cos(np.radians(ang))],[0,r*np.sin(np.radians(ang))],'--',color=c,lw=1.2,alpha=.7)
a.plot(lam0.real,lam0.imag,'o',color='tab:blue',ms=11,zorder=5,label=r'$\lambda_0=0.890-0.890i$ (det.)')
# mean shifted roots
for e2,al in [(0.1,.35),(0.2,.55),(0.4,.85)]:
    lam=lam0+e2*O2; a.plot(lam.real,lam.imag,'s',color='tab:green',ms=7,alpha=al,zorder=4)
    a.annotate(f'$\\eta^2$={e2}',(lam.real,lam.imag),textcoords='offset points',xytext=(6,4),fontsize=7,color='tab:green')
lamf=lam0+0.4*O2
a.annotate('',xy=(lamf.real,lamf.imag),xytext=(lam0.real,lam0.imag),
           arrowprops=dict(arrowstyle='->',color='tab:green',lw=1.6))
a.plot([],[],'s',color='tab:green',label=r'$E[\lambda]=\lambda_0+\eta^2\Omega_2$ (mean; stays ~ -45$^\circ$)')
# fluctuation ellipse at eta^2=0.3 (pseudo-var anisotropy)
e2=0.3; th=np.linspace(0,2*np.pi,200)
# pseudo-variance E[dlam^2]=pv*eta^2 gives anisotropy; semi-axes from |true var| +- |pv|
sig=np.sqrt(truevar*e2); axmaj=np.sqrt(max(truevar+abs(pv),0)*e2); axmin=np.sqrt(max(truevar-abs(pv),0)*e2)
phi=np.angle(pv)/2
ell=lam0+(axmaj*np.cos(th)*np.exp(1j*phi)+1j*axmin*np.sin(th)*np.exp(1j*phi))
a.plot(ell.real,ell.imag,':',color='tab:purple',lw=1.4,label=r'fluctuation 1$\sigma$ ($\eta^2$=0.3): rms$\approx$1.34$\eta$')
a.set_xlabel(r'Re $\lambda$'); a.set_ylabel(r'Im $\lambda$'); a.set_title('Noise-dressed connection root: mean barely rotates,\nfluctuation is large & anisotropic')
a.legend(fontsize=7.5,loc='upper left'); a.set_xlim(0.2,1.6); a.set_ylim(-1.6,-0.2); a.set_aspect('equal')
# ---- Panel B: lemma test delta(0) vs Ito
b=ax[1]
hs=np.array([4e-3,2e-3,1e-3,5e-4]); locs=np.array([abs(l) for _,l in m.location_shift_test(lam0,h_list=hs)])
A_R=np.array([abs(m.chaos_pieces(lam0,h=h,theta=0.3,To=18.0)['A_R']) for h in hs])
b.loglog(1/hs,locs,'o-',color='tab:red',label=r'Borel LOCATION shift $E[\delta^2 P]$ (period): $\propto 1/h$ = $\delta(0)$ DIVERGENT')
b.loglog(1/hs,A_R,'s-',color='tab:green',label=r'connection datum $\int u_R^4$ (Ito): FINITE, convergent')
b.loglog(1/hs,locs[0]*(hs[0]/hs),'--',color='0.6',lw=1,label=r'reference slope $\propto 1/h$')
b.set_xlabel('1/h  (grid resolution)'); b.set_ylabel('|O($\\eta^2$) mean shift|')
b.set_title('LEMMA TEST: location shift is $\\delta(0)$-divergent (needs v3 renorm);\nconnection datum is renorm-clean (non-local Ito functional)')
b.legend(fontsize=7.5,loc='center left'); b.grid(True,which='both',alpha=.25)
fig.suptitle('Program 2 route(b) sub-attack 1: O($\\eta^2$) stochastic Stokes constant $\\Omega_2=-0.451+0.352i$ '
             r'(|.|=0.572, arg=142$^\circ\approx$3$\pi$/4)',fontsize=11.5,fontweight='bold')
fig.tight_layout(rect=[0,0,1,0.95]); fig.savefig("figures/stochastic_stokes_o_eta2.png",dpi=120)
print("saved figures/stochastic_stokes_o_eta2.png")
