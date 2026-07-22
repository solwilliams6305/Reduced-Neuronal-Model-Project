"""
Artifact-free confirmation that the left tail is the instanton law -logP = I(s)/eta^2 (q=1), NOT q=0.72.
MC has NO numerical diffusion (unlike the coarse FP), so the instanton predicts the collapse
    eta^2 * (-log P(Theta>t))  ->  I(sqrt(2 Theta))   (beta-INDEPENDENT).
Test on existing MC (persistence_data.npz: Th_2,Th_4,Th_8) vs a fresh solve_bvp interpolant of I(s).
If the FP's q=0.72 were real, this collapse would FAIL (curves would fan out as beta^{-0.28}).
"""
import numpy as np
from scipy.integrate import solve_bvp
def V(t): return -np.sign(t)*t*t
def Iaction(s,p0=0.0,M=20.0,t0=1e-3,n=1000):
    t=np.linspace(t0,s,n)
    def ode(t,y): p,pi=y; return np.vstack([pi+V(t)-p*p,2*p*pi])
    def bc(ya,yb): return np.array([ya[0]-p0,yb[0]+M])
    pg=np.clip(np.where(t<0.9*s,1.0/np.maximum(t,0.3),-M*(t-0.9*s)/(0.1*s)),-M,3.0)
    sol=solve_bvp(ode,bc,t,np.vstack([pg,t*t]),max_nodes=200000,tol=1e-6)
    if not sol.success: return np.nan
    tt=np.linspace(t0,s,3000); pp,pi=sol.sol(tt); return 0.5*np.trapz(pi**2,tt)
# I(s) interpolant
sgrid=np.arange(2.0,6.6,0.25); Igrid=np.array([Iaction(s) for s in sgrid])
def Iinterp(s): return np.interp(s,sgrid,Igrid)

d=np.load("persistence_data.npz")
print("=== eta^2*(-logP) vs instanton I(sqrt(2Theta)); collapse (beta-indep) => q=1 confirmed ===")
print(f"  {'Theta':>6} {'s':>5} {'I(s)':>8} | " + " ".join(f"b={int(b):<2d}:eta^2*(-logP)" for b in [2,4,8]))
for Th in [3.0,3.5,4.0,4.5,5.0]:
    s=np.sqrt(2*Th); row=f"  {Th:6.1f} {s:5.2f} {Iinterp(s):8.2f} | "
    for b in [2,4,8]:
        Thd=d[f"Th_{b}"]; tot=Thd.size; cnt=(Thd>Th).sum()
        if cnt>=20:
            eta2=4.0/b; row+=f"   {eta2*(-np.log(cnt/tot)):6.2f}      "
        else:
            row+="     (<20)     "
    print(row)
print("\n  If the three eta^2*(-logP) columns agree with each other AND ~ I(s): q=1 (instanton), artifact-free.")
print("  (FP's q=0.72 would instead make them DECREASE with beta as beta^{-0.28}.)")
