"""
Item-2 DERIVATION: the left-tail (persistence) instanton action, artifact-free (ODE BVP, no diffusive PDE).
FW: -log P(Y*<-s) = I(s)/eta^2 with I(s)=min (1/2) int phi^2 dt, phi = pdot - V + p^2 the noise control that
prevents Riccati explosion until depth t=s.  Hamilton: pdot=pi+V-p^2, pidot=2 p pi, action=(1/2)int pi^2.
Claim to verify: I(s) -> s^5/10  (=> exponent 5, q=1, constant 1/10; the measured q~0.72 is then pre-asymptotic).
scipy.solve_bvp. [DERIVED + NUMERIC check].
"""
import numpy as np
from scipy.integrate import solve_bvp

def V(t):  # t=-Y; oscillatory side t>0: V=-t^2; confining t<0: V=+t^2
    return -np.sign(t)*t*t

def instanton_action(s, p0=0.0, M=20.0, t0=1e-3, n=1200):
    # EL system on [t0,s], BC: p(t0)=p0 (incoming, O(1)); p(s)=-M (explosion threshold)
    t=np.linspace(t0,s,n)
    def ode(t,y):
        p,pi=y
        return np.vstack([pi+V(t)-p*p, 2*p*pi])
    def bc(ya,yb):
        return np.array([ya[0]-p0, yb[0]+M])
    # initial guess: p ~ small then dive to -M; pi ~ t^2
    pg=np.where(t<0.9*s, 1.0/np.maximum(t,0.3), -M*(t-0.9*s)/(0.1*s))
    pg=np.clip(pg,-M,3.0); pig=t*t
    sol=solve_bvp(ode,bc,t,np.vstack([pg,pig]),max_nodes=200000,tol=1e-6,verbose=0)
    if not sol.success: return np.nan,sol
    tt=np.linspace(t0,s,4000); pp,pipi=sol.sol(tt)
    S=0.5*np.trapz(pipi**2,tt)
    return S,sol

print("=== instanton action I(s) vs the predicted asymptote s^5/10 ===")
ss=np.array([2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0,13.0])
Is=[]
for s in ss:
    S,_=instanton_action(s)
    Is.append(S)
    print(f"  s={s:5.1f}: I(s)={S:9.3f}   s^5/10={s**5/10:9.3f}   ratio I/(s^5/10)={S/(s**5/10):.3f}")
Is=np.array(Is)
# asymptotic exponent and constant from the largest s
big=ss>=5
pexp=np.polyfit(np.log(ss[big]),np.log(Is[big]),1)[0]
Cbig=Is[-1]/ss[-1]**5
print(f"\n  large-s exponent of I(s) = {pexp:.3f}  (predicted 5)")
print(f"  large-s constant I/s^5   = {Cbig:.4f} (predicted 0.100)")
print(f"  => -logP(Y*<-s) = I(s)/eta^2 -> s^5/(10 eta^2) = beta s^5/40 :  exponent 5, q=1 (DERIVED)")

# p0 / M robustness (the s^5 coefficient must be universal)
print("\n=== robustness of the s^5 coefficient (should be ~0.1 regardless of p0, M) ===")
for p0 in [-0.5,0.0,0.5]:
    S,_=instanton_action(10.0,p0=p0); print(f"  p0={p0:+.1f}: I(10)/10^5={S/1e5:.4f}")
for M in [12.0,20.0,40.0]:
    S,_=instanton_action(10.0,M=M); print(f"  M={M:4.0f}: I(10)/10^5={S/1e5:.4f}")

# effective exponent q_eff over the MEASURED s-range if one (wrongly) fits -logP vs beta at fixed moderate s:
# with -logP = I(s)/eta^2 EXACTLY, q_eff=1 at fixed s.  Show the s-crossover of the local power a(s)=dlogI/dlog s.
print("\n=== pre-asymptotic local power a(s)=dlogI/dlogs (Theta=s^2/2; a_Theta=a_s/2) ===")
for i in range(1,len(ss)-1):
    a=np.log(Is[i+1]/Is[i-1])/np.log(ss[i+1]/ss[i-1])
    print(f"  s={ss[i]:4.1f} (Theta={ss[i]**2/2:4.1f}): a_s={a:.2f} -> a_Theta={a/2:.2f}  (asymptote 5 -> 2.5)")
