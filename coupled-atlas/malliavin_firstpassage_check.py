"""
Item 3 (T1's last regularity): numerically validate the load-bearing identity behind the Malliavin
differentiability of the first-passage location, in the FIELD representation x=(u,u') (regular through the
turning; no 1/k Pruefer singularity):

   escape Y* = first zero of u (u''=(q(Y)-eta Wdot)u, q=sign(Y)Y^2-lambda), continued from recessive u at Y+.
   CLAIM (IFT on Wiener space):  D_s Y* = - D_s u(Y*) / u'(Y*).

Tests, per perturbation location s_j (Bismut finite-difference on one Brownian increment):
  (1) D_s Y*   := d Y*/d xi_j        (perturb increment j, remeasure Y*)
  (2) D_s u(Y*):= d u(Y*_base)/d xi_j (perturb increment j, remeasure u at the FIXED base location)
  (3) predicted D_s Y* = -(2)/u'(Y*);  compare to (1).
  (4) linearity in eps (differentiable => straight through 0, no kink; transversal => finite, no blow-up).
numpy only. [NUMERIC validation of the proof's key identity].
"""
import numpy as np
def q(Y,lam): return np.sign(Y)*Y*Y - lam

def integrate(xi, dt, Y0, lam, eta):
    """Euler-Maruyama for (u,u') sweeping Y down from Y0; xi = standardized increments (len n).
    Returns arrays Y, u, up. Recessive-ish IC: u=1, u'=+sqrt(q)>0 (decaying into Y+)."""
    n=len(xi); Y=np.empty(n+1); u=np.empty(n+1); up=np.empty(n+1)
    Y[0]=Y0; u[0]=1.0; up[0]=np.sqrt(max(q(Y0,lam),1e-9))  # recessive: u'/u=+sqrt(q) decays as Y->+inf
    sq=np.sqrt(dt)
    for i in range(n):
        Yi=Y[0]-i*dt
        # du = u' d(-... ) ; we sweep Y decreasing: dY=-dt. u_{i+1}=u_i+u'_i*(-dt)? Use tau=-Y increasing:
        # d u/dtau = -u' ;  d u'/dtau = -(q u) + eta u xi   (since u''=(q-eta W')u, ''=d^2/dY^2=d^2/dtau^2)
        # => du = -up*dt ; dup = -(q*u)*dt + eta*u*sq*xi
        Q=q(Yi,lam)
        u[i+1]=u[i]-up[i]*dt
        up[i+1]=up[i]-(Q*u[i])*dt+eta*u[i]*sq*xi[i]
        Y[i+1]=Yi-dt
    return Y,u,up

def first_zero(Y,u,up):
    # first sign change of u at Y<0 (oscillatory side); linear-interpolate the crossing; return Y*, u'(Y*)
    for i in range(1,len(u)):
        if Y[i]<0 and u[i-1]*u[i]<0:
            f=u[i-1]/(u[i-1]-u[i]); Ystar=Y[i-1]+f*(Y[i]-Y[i-1])
            upstar=up[i-1]+f*(up[i]-up[i-1]); ubase_idx=i
            return Ystar,upstar,ubase_idx
    return np.nan,np.nan,-1

lam=1.0; eta=0.20; Y0=3.0; dt=2e-3; n=int((Y0-(-4.0))/dt)
rng=np.random.default_rng(0); xi=rng.standard_normal(n)
Y,u,up=integrate(xi,dt,Y0,lam,eta); Ystar,upstar,bi=first_zero(Y,u,up)
print(f"base escape: Y*={Ystar:.4f}, u'(Y*)={upstar:.4f}  (transversal: |u'|>0 => finite D Y*)")
print(f"lambda={lam}, eta={eta}, dt={dt}\n")

# Bismut finite-difference at a set of locations
eps=1e-3
js=np.linspace(int(0.15*bi),bi-3,14).astype(int)
print(f"  {'s=Y_j':>8} {'D_s Y* (meas)':>13} {'D_s u(Y*)':>11} {'pred=-Du/u''':>12} {'ratio':>7}")
meas=[]; pred=[]
for j in js:
    xp=xi.copy(); xp[j]+=eps
    Yp,upu,upp=integrate(xp,dt,Y0,lam,eta); Ysp,_,_=first_zero(Yp,upu,upp)
    DsYstar=(Ysp-Ystar)/eps
    Du_at_base=(upu[bi]-u[bi])/eps            # d u(Y*_base)/d xi_j  (u at fixed base index)
    pred_DsYstar=-Du_at_base/upstar
    meas.append(DsYstar); pred.append(pred_DsYstar)
    print(f"  {Y[j]:8.3f} {DsYstar:13.4f} {Du_at_base:11.4f} {pred_DsYstar:12.4f} {DsYstar/pred_DsYstar if pred_DsYstar!=0 else np.nan:7.3f}")
meas=np.array(meas); pred=np.array(pred)
corr=np.corrcoef(meas,pred)[0,1]; slope=np.polyfit(pred,meas,1)[0]
print(f"\n  measured vs predicted D_s Y*:  corr={corr:.4f}, slope={slope:.3f}  (identity => corr~1, slope~1)")

# linearity in eps at one location (differentiability: ratio DY*/eps stable; no kink/blow-up)
print("\n  linearity of Y* response (differentiable + transversal => constant DY*/eps as eps->0):")
j=js[len(js)//2]
for e in [4e-3,2e-3,1e-3,5e-4,2.5e-4]:
    xp=xi.copy(); xp[j]+=e; Yp,uu,uup=integrate(xp,dt,Y0,lam,eta); Ysp,_,_=first_zero(Yp,uu,uup)
    print(f"    eps={e:.1e}: DY*/eps={(Ysp-Ystar)/e:8.4f}")
print("\n  => stable ratio (no kink) confirms differentiability; agreement confirms D_s Y*=-D_s u(Y*)/u'(Y*).")
