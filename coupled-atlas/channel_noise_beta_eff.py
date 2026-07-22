import numpy as np
from scipy.integrate import solve_ivp
# Deterministic backbone u_c: sensitivity kernel Phi ~ u_c^2, functional weight w=Phi^2 ~ u_c^4.
# beta_eff = 4 INT w / INT eta^2 w  (Ito isometry, Tier A). Check: (1) where does w peak? (2) does the
# u_c^4-weighted <phi^2> predict the measured variance shift for the +/-10% profile?
def V(Y): return np.sign(Y)*Y*Y
Y0=3.0; Yend=-6.0
sol=solve_ivp(lambda Y,z:[z[1],V(Y)*z[0]],[Y0,Yend],[1.0,np.sqrt(V(Y0))],rtol=1e-11,atol=1e-13,dense_output=True,max_step=1e-3)
Yg=np.linspace(Y0,Yend,60000); uc=sol.sol(Yg)[0]
# first node Y* (Y<0)
i=next(k for k in range(1,len(Yg)) if Yg[k]<0 and uc[k-1]*uc[k]<0)
Ystar=Yg[i-1]-uc[i-1]*(Yg[i]-Yg[i-1])/(uc[i]-uc[i-1])
mask=Yg>=Ystar; s=Yg[mask]; w=uc[mask]**4          # weight Phi^2 ~ u_c^4 on [Y*,Y0]
w=w/np.trapz(w[::-1],s[::-1])                        # normalize (s descending)
# where does w peak / concentrate?
order=np.argsort(s); s_o=s[order]; w_o=w[order]
speak=s_o[np.argmax(w_o)]
# fraction of weight in turning band [-0.5,0.5], confining [0.5,3], oscillatory-landing [Y*,-0.5]
def frac(lo,hi):
    m=(s_o>=lo)&(s_o<=hi); return np.trapz(w_o[m],s_o[m])
print(f"Y*_det={Ystar:.3f}; sensitivity weight w=u_c^4 on [Y*,Y0]:")
print(f"  peak of w at s={speak:.2f}")
print(f"  weight fraction: landing [Y*,-0.5]={frac(Ystar,-0.5):.2f}  turning [-0.5,0.5]={frac(-0.5,0.5):.2f}  confining [0.5,3]={frac(0.5,3):.2f}")
# beta_eff variance formula: predicted Var(profile)/Var(const) = <phi^2>_w for eta=eta0*phi
for amp in [0.10,0.15]:
    phi=1+amp*np.tanh(s_o); pred=np.trapz(phi**2*w_o,s_o)/np.trapz(w_o,s_o)
    print(f"  +/-{int(amp*100)}% tanh: predicted Var(profile)/Var(const) = <phi^2>_w = {pred:.3f}")
print("  [measured FP std: const 0.689 -> +/-10% 0.715 => Var ratio (0.715/0.689)^2 = 1.077; +/-15% 0.730 => 1.123]")
print("\nOFF-FAMILY check (is the channel-noise law any single beta?):")
print("  beta-family: std DECREASES with beta (b1..8: .904/.689/.482/.320), skew INCREASES (.26/.60/.85/.85).")
print("  mild +/-10%: std=0.715 (>0.689 => beta<2) BUT skew=0.701 (>0.601 => beta>2) -- INCONSISTENT.")
print("  => channel noise deforms the law OFF the additive beta-family (same CLASS/tail, not a single beta).")
