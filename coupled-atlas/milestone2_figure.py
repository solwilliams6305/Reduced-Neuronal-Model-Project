import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from math import factorial as fac
import milestone2_contamination as MC

fig,ax=plt.subplots(1,2,figsize=(13.5,5.2))
# Panel A: identifiability -- residual vs true phase (flat => can't discriminate)
ns=np.arange(1,6); tgt=MC.REAL[1:6]; rho=1.258
tts=np.arange(35,66,1); curve=[]
for tt in tts:
    th=np.radians(tt); best=1e9
    for rr in np.linspace(1.0,2.4,80):
        Bm=np.vstack([2*rho**(-(ns+1))*np.cos((ns+1)*th)*np.array([fac(k) for k in ns]),
                      2*rho**(-(ns+1))*np.sin((ns+1)*th)*np.array([fac(k) for k in ns]),
                      rr**(-(ns+1))*np.array([fac(k) for k in ns])]).T
        coef,*_=np.linalg.lstsq(Bm,tgt,rcond=None); r=np.sqrt(np.mean((Bm@coef-tgt)**2)); best=min(best,r)
    curve.append(best)
curve=np.array(curve)
a=ax[0]
a.plot(tts,curve/curve.min(),'-',color='tab:blue',lw=2)
a.axvspan(35,56,color='tab:blue',alpha=0.08)
a.axvline(45,color='tab:green',ls='--',lw=1.5,label='det. 45$^\\circ$ (from $\\lambda_0$)')
a.axvline(54,color='tab:red',ls='--',lw=1.5,label='claimed 54$^\\circ$')
a.axhline(1.0,color='0.7',lw=0.7)
a.set_xlabel('true conjugate-pair phase $\\theta_{\\rm true}$ (deg)')
a.set_ylabel('fit residual / min residual')
a.set_title('(A) 6 coefficients CANNOT discriminate the phase:\nresidual is flat over [35,56]$^\\circ$ — 45$^\\circ$ fits as well as 54$^\\circ$')
a.legend(fontsize=9,loc='upper left'); a.set_ylim(0.9,3.0); a.grid(alpha=.2)
# Panel B: contamination bias -- true 45 pair + real pole -> apparent phase
Rs=np.linspace(0,2.5,26)
p22=[];pdar=[]
for R in Rs:
    v=MC.make_ladder(45,rho_pair=1.258,rho_real=1.5,R=R)
    r=MC.borel_pade_phase(v,2,2); p22.append(r[0] if r else np.nan)
    pdar.append(MC.darboux_phase(v)[0])
b=ax[1]
b.plot(Rs,p22,'o-',color='tab:purple',ms=4,label='[2/2] Borel-Pade')
b.plot(Rs,pdar,'s-',color='tab:orange',ms=4,label='Darboux')
b.axhline(45,color='tab:green',ls='--',lw=1.5,label='TRUE phase = 45$^\\circ$')
b.axhspan(53,59,color='tab:red',alpha=0.12,label='claimed 54-59$^\\circ$ band')
b.set_xlabel('real instanton-pole amplitude R (at $|\\zeta|$=1.5)')
b.set_ylabel('APPARENT pair phase reported (deg)')
b.set_title('(B) The known real instanton pole BIASES a true 45$^\\circ$\nup to 54-64$^\\circ$ in exactly the estimators used')
b.legend(fontsize=9,loc='upper left'); b.grid(alpha=.2); b.set_ylim(40,72)
fig.suptitle('Milestone 2: the ~54$^\\circ$-vs-45$^\\circ$ stochastic Borel-phase gap is NOT robustly established from 6 coefficients',
             fontsize=12,fontweight='bold')
fig.tight_layout(rect=[0,0,1,0.95]); fig.savefig("figures/milestone2_phase_identifiability.png",dpi=120)
print("saved figures/milestone2_phase_identifiability.png")
