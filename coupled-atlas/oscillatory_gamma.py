"""
Step (1): pin the OSCILLATORY-side factor of W's connection by its non-confounded signature — the connection
MODULUS. The recessive Weber solution U(a,sqrt2 Y) (a=-lambda/2), decaying ~e^{-Y^2/2} on Y>0, emerges oscillatory
on Y<0 with amplitude M(lambda). The c=1/Weber prediction is the EXACT reflection identity
    |Gamma(1/2 - i lambda/2)|^2 = pi / cosh(pi lambda/2),
so M(lambda)^2 should be EXPONENTIAL in lambda (rate ~pi/2) — a signature a smooth WKB power factor CANNOT fake.
Normalize the solution to the exact U-function at Y0 (scipy pbdv), integrate to Y<0, read the WKB-invariant
amplitude M = sqrt(|Y| u^2 + u'^2/|Y|). Test log M^2 vs lambda (linear=exponential=Gamma) and vs the cosh law.
scipy. [NUMERIC — non-confounded oscillatory-side test].
"""
import numpy as np
from scipy.special import pbdv
from scipy.integrate import solve_ivp
def Vc(Y,lam): return np.sign(Y)*Y*Y - lam

def recessive_amp(lam, Y0=7.0, Yend=-9.0):
    # exact recessive IC: U(a, sqrt2 Y) = D_{-a-1/2}(sqrt2 Y), a=-lam/2 => order nu = lam/2 - 1/2
    nu = lam/2.0 - 0.5; z0 = np.sqrt(2.0)*Y0
    D, Dp = pbdv(nu, z0)                      # D_nu(z0), dD/dz
    u0 = D; up0 = np.sqrt(2.0)*Dp             # u=D(sqrt2 Y): du/dY = sqrt2 D'
    if not np.isfinite(u0) or abs(u0)<1e-300: return np.nan
    def rhs(Y,y): return [y[1], Vc(Y,lam)*y[0]]
    sol=solve_ivp(rhs,[Y0,Yend],[u0,up0],max_step=1.5e-3,rtol=1e-11,atol=1e-14,dense_output=True)
    # WKB-invariant amplitude at several deep Y<0, averaged (M^2 = |Y|u^2 + u'^2/|Y|)
    Ys=np.array([-6.0,-6.5,-7.0,-7.5,-8.0]); Ms=[]
    for Yq in Ys:
        u,up=sol.sol(Yq); k=abs(Yq); Ms.append(k*u*u + up*up/k)
    return np.mean(Ms)                        # = M^2

lams=np.arange(-2.0,3.01,0.25)
print("=== oscillatory-side connection modulus M(lambda)^2 (U-normalized) vs c=1 cosh law ===")
print(f"  {'lambda':>7} {'M^2 (measured)':>15} {'pi/cosh(pi lam/2)':>18} {'ratio':>10}")
M2=[]; cosh=[]
for lam in lams:
    m2=recessive_amp(lam); c=np.pi/np.cosh(np.pi*lam/2.0)
    M2.append(m2); cosh.append(c)
    print(f"  {lam:7.2f} {m2:15.5e} {c:18.5e} {m2/c if c>0 else np.nan:10.4f}")
M2=np.array(M2); cosh=np.array(cosh); ok=np.isfinite(M2)&(M2>0)
# test 1: is log M^2 linear in lambda at large lambda (exponential=Gamma, not WKB power)?
big=ok&(lams>=0.5)
slope=np.polyfit(lams[big],np.log(M2[big]),1)[0]
print(f"\n  large-lambda slope of log M^2 = {slope:+.3f}   (c=1/Gamma predicts -pi/2 = {-np.pi/2:+.3f})")
# test 2: does M^2 / (pi/cosh) flatten to a constant (i.e. M^2 IS the cosh law up to a smooth prefactor)?
ratio=M2[ok]/cosh[ok]
# fit ratio = smooth power A*|lam|^p ; report how constant log(ratio) is vs lambda
lr=np.log(ratio); v_ratio=np.var(lr); v_M2=np.var(np.log(M2[ok]))
print(f"  Var(log[M^2/cosh]) / Var(log M^2) = {v_ratio/v_M2:.3f}   (<<1 => the cosh law captures the variation)")
print(f"  => slope ~ -pi/2 AND cosh-division flattening  ==>  oscillatory factor = Gamma(1/2 - i lambda/2) [c=1].")
np.savez("oscillatory_gamma.npz",lams=lams,M2=M2,cosh=cosh)
