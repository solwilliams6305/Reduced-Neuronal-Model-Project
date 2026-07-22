"""
Deterministic asymmetric connection data of W's skeleton, and the c=1 / Weber Gamma-function test.

Skeleton:  u'' = (sign(Y) Y^2 - lambda) u.   Y>0: confining Weber (real, ~e^{-Y^2/2}); Y<0: inverted/oscillatory
Weber. The canard/attracting solution (log-deriv p=+sqrt(V) at Y0) is continued down through the turning to Y<0;
its first zero = Y*_det (the deterministic anchor of W, ~ -2.1) and its oscillatory phase = the CONNECTION PHASE
phi(lambda). Question: is phi(lambda) the closed-form c=1/Weber Gamma-function phase?

c=1 / DLMF-12 candidates for the connection phase (confining U-side ~ 1/Gamma(3/4+a/2),Gamma(1/4+a/2);
oscillatory W-side ~ arg Gamma(1/2 - i lambda/2)), a=-lambda/2:
   phi_pred(lambda) = arg Gamma(1/2 - i*lambda/2) + c0 + c1*lambda   (fit the two constants; test the Gamma shape)
We compute phi by integration for a grid of lambda and check whether the Gamma SHAPE (not a free curve) fits.
scipy available. [NUMERIC + analytic test].
"""
import numpy as np
from scipy.special import loggamma
from scipy.integrate import solve_ivp

def V(Y,lam): return np.sign(Y)*Y*Y - lam

def nodes_and_phase(lam, Y0=6.0, Yend=-9.0):
    # canard IC: p=+sqrt(V)>0 at Y0 (attracting branch).
    V0=V(Y0,lam); u0=1.0; up0=np.sqrt(max(V0,1e-9))
    def rhs(Y,y): return [y[1], V(Y,lam)*y[0]]
    sol=solve_ivp(rhs,[Y0,Yend],[u0,up0],max_step=1.5e-3,rtol=1e-11,atol=1e-13,dense_output=True)
    Yg=np.linspace(Y0,Yend,20000); u,_=sol.sol(Yg)
    nodes=[]
    for i in range(1,len(Yg)):
        if Yg[i]<0 and u[i-1]*u[i]<0:
            f=u[i-1]/(u[i-1]-u[i]); nodes.append(Yg[i-1]+f*(Yg[i]-Yg[i-1]))
    nodes=np.array(nodes)
    Th=nodes*nodes/2.0                       # phase-depth of each node (increasing)
    # nodes obey Theta_n = n*pi - phi  => spacing pi; fit to extract phi robustly (skip 1st node: turning transient)
    n=np.arange(len(Th))
    use=n>=1
    slope,inter=np.polyfit(n[use],Th[use],1)  # slope~pi ; Theta = slope*n + inter
    # connection phase: at integer n, Theta = n*pi - phi  => -phi = inter (mod pi); phi = -inter mod pi
    phi=np.mod(-inter,np.pi)
    return nodes[0] if len(nodes) else np.nan, Th, slope, phi

print("=== ground-truth skeleton: first node + connection phase phi(lambda), FINE grid ===")
lams=np.arange(-2.6,0.91,0.10)     # stay below the lambda=1 confining pole
Y1=[]; phis=[]
for lam in lams:
    n1,Th,slope,phi=nodes_and_phase(lam); Y1.append(n1); phis.append(phi)
Y1=np.array(Y1); phis=np.array(phis)
phi_un=np.unwrap(2*phis)/2         # remove the mod-pi node-lattice wrap
print(f"  lambda=0: Y*_det={Y1[np.argmin(abs(lams))]:+.3f} (target -2.19).  phi range [{phi_un.min():.2f},{phi_un.max():.2f}]")

def argG(z): return np.imag(loggamma(z))
# DECISIVE: de-linearize. phi = linear(WKB) + Stokes-Gamma(nonlinear). Subtract best line from BOTH phi and each
# candidate, then correlate the NONLINEAR residuals — this isolates the connection-coefficient (Gamma) content.
def delin(y):
    c=np.polyfit(lams,y,1); return y-np.polyval(c,lams)
r_phi=delin(phi_un)
cands={
  "argG(1/2 - i lam/2)  [c=1 oscillatory]": np.array([argG(0.5-0.5j*l) for l in lams]),
  "argG(1/4 - i lam/4)"                    : np.array([argG(0.25-0.25j*l) for l in lams]),
  "argG(3/4 - i lam/4)"                    : np.array([argG(0.75-0.25j*l) for l in lams]),
  "argG(1/2 - i lam) [double]"             : np.array([argG(0.5-1.0j*l) for l in lams]),
  "-(pi/2)*lam*ln|lam| WKB-only (null)"    : np.array([-(np.pi/2)*l*np.log(abs(l)+1e-9) for l in lams]),
}
print("\n=== de-linearized residual test (isolates the Stokes-Gamma content) ===")
print("  correlation of nonlinear(phi) with nonlinear(candidate) + slope of the match:")
best=None
for name,g in cands.items():
    rg=delin(g);
    corr=np.corrcoef(r_phi,rg)[0,1]; slope=np.polyfit(rg,r_phi,1)[0]
    amp=np.std(r_phi)/max(np.std(rg),1e-12)
    print(f"   {name:38s}: corr={corr:+.4f}  match-slope={slope:+.3f}  (resid amp ratio {amp:.2f})")
    if best is None or abs(corr)>abs(best[1]): best=(name,corr)
print(f"\n  best nonlinear match: {best[0]} (corr {best[1]:+.4f})")
print("  => a HIGH residual correlation with a specific Gamma (and near-zero with the WKB-only null) is the")
print("     decisive signature that the connection phase carries that c=1/Weber Gamma Stokes coefficient.")
np.savez("connection_gamma.npz",lams=lams,Y1=Y1,phis=phis,phi_un=phi_un,r_phi=r_phi)
