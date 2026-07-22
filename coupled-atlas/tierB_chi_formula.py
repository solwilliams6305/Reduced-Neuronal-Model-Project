"""
Tier B / Lemma L1 groundwork: verify the closed-form susceptibility
    chi = -(1/k*^2) * INT_{Y*}^{Y0} s u(s)^2 ds ,   k* = u'(Y*)
against a direct finite-difference dY*/dDelta, DETERMINISTICALLY (noise off) -- decisive formula check.
Then STOCHASTICALLY (beta=2) sample k* and chi to determine the small-ball law P(k*<eps): is it power-law
(=> L1 needs a real small-ball estimate) or lognormal/faster (=> E[k*^-2]<inf easily)?  This decides L1's difficulty.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def Vd(Y,Delta): a=abs(Y); return np.sign(Y)*a*(a+Delta)

# ---------- deterministic recessive solution + first node ----------
def det_node(Delta, Y0=6.0, Yend=-5.0):
    u0=np.exp(-Y0*Y0/2); up0=-Y0*u0
    sol=solve_ivp(lambda Y,z:[z[1],Vd(Y,Delta)*z[0]],[Y0,Yend],[u0,up0],
                  rtol=1e-11,atol=1e-14,dense_output=True,max_step=0.005)
    # first node Y<0
    Yg=np.linspace(0.0,Yend,60000); ug=sol.sol(Yg)[0]
    i=np.where(ug[:-1]*ug[1:]<0)[0][0]
    Ys=brentq(lambda y:sol.sol(y)[0],Yg[i],Yg[i+1])
    return Ys, sol

def chi_formula(Y0=6.0, Yend=-5.0):
    Ys,sol=det_node(0.0,Y0,Yend)
    kstar=sol.sol(Ys)[1]
    ss=np.linspace(Ys,Y0,80000); us=sol.sol(ss)[0]
    J=np.trapz(ss*us**2, ss)
    return -J/kstar**2, Ys, kstar, J

if __name__=="__main__":
    print("=== (1) DETERMINISTIC check of chi formula vs finite-difference dY*/dDelta ===")
    chi_f,Ys0,kstar,J=chi_formula()
    dd=1e-3
    Yp,_=det_node(+dd); Ym,_=det_node(-dd)
    chi_fd=(Yp-Ym)/(2*dd)
    print(f"  Y*_det={Ys0:.5f}, k*=u'(Y*)={kstar:.4f}, INT s u^2={J:.5e}")
    print(f"  chi (closed form -J/k*^2) = {chi_f:.5f}")
    print(f"  chi (finite-diff dY*/dD)  = {chi_fd:.5f}")
    print(f"  match: {'YES' if abs(chi_f-chi_fd)<1e-3 else 'NO'}  (rel err {abs(chi_f-chi_fd)/abs(chi_fd):.1e})")

    print("\n=== (2) STOCHASTIC small-ball of k* at beta=2: is P(k*<eps) power-law or faster? ===")
    # sample the noisy FIELD (u,u') via the linear SDE, find first node, record k*=|u'(Y*)|
    eta=np.sqrt(2.0); N=400000; Y0=3.0; Yend=-6.0; dt=1.0e-3; sq=np.sqrt(dt)
    rng=np.random.default_rng(7)
    u=np.full(N,np.exp(-Y0*Y0/2)); up=np.full(N,-Y0*u[0])
    kstar=np.full(N,np.nan); Ystar=np.full(N,np.nan)
    n=int(round((Y0-Yend)/dt))
    uprev=u.copy()
    for i in range(n):
        Y=Y0-i*dt
        dW=sq*rng.standard_normal(N)
        unew=u+up*dt
        upnew=up+Vd(Y,0.0)*u*dt+eta*u*dW
        # detect first node: u changed sign (u>0 -> <=0), record where not yet recorded
        cross=np.isnan(kstar)&(u>0)&(unew<=0)
        kstar[cross]=np.abs(upnew[cross]); Ystar[cross]=Y-dt
        u,up=unew,upnew
    ok=~np.isnan(kstar); k=kstar[ok]
    print(f"  nodes found: {ok.sum()}/{N};  k* stats: mean={k.mean():.3f} median={np.median(k):.3f} min={k.min():.4f}")
    # small-ball: P(k*<eps) vs eps  (log-log slope = theta; chi~k^-2 tail alpha=theta/2)
    for eps in [0.5,0.3,0.2,0.12,0.07,0.04]:
        p=np.mean(k<eps); print(f"    P(k*<{eps:.2f})={p:.3e}")
    eps_grid=np.array([0.3,0.2,0.12,0.07,0.04]); P=np.array([np.mean(k<e) for e in eps_grid])
    good=P>1e-4
    if good.sum()>=3:
        theta=np.polyfit(np.log(eps_grid[good]),np.log(P[good]),1)[0]
        print(f"  => P(k*<eps) ~ eps^theta, theta={theta:.2f}  (chi~k*^-2 => tail alpha=theta/2={theta/2:.2f}; measured alpha~2)")
        print(f"     E[k*^-2]<inf iff theta>2: {'YES (L1 holds)' if theta>2 else 'BORDERLINE/NO'}")
    # direct: E[k*^-2] and its stability across halves
    h=len(k)//2
    for q in [1,2,3]:
        m=np.mean(k**(-q)); m1=np.mean(k[:h]**(-q)); m2=np.mean(k[h:]**(-q))
        print(f"    E[k*^-{q}]={m:.3f} (halves {m1:.3f},{m2:.3f}, ratio {m2/m1:.2f})")
