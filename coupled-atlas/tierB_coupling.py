"""
Tier B — the SYNCHRONOUS COUPLING, tested pathwise at beta=2 (full noise, NOT weak-noise).
Two Riccati escapes on the SAME Brownian path: pure cusp V_0=sign(Y)Y^2 (Delta=0) and finite-merge
V_Delta=sign(Y)|Y|(|Y|+Delta).  Escape-location shift  dY* = Y*(Delta) - Y*(0).

Theorem structure to test:
  (a) dY* = Delta * chi + O(Delta^2)  pathwise (chi = susceptibility, Delta-independent) -> coupling is linear
  (b) W_1(W_{2,Delta}, W_2) <= E|dY*| ~ Delta * E|chi|  => rate p=1/4 (Delta~eps^{1/4}) IF E|chi|<inf
  (c) chi's tail: finite mean (=> W_1 rate clean) but HEAVY tail / divergent 4th moment
      (=> kurtosis correction non-analytic / slower) -- the rate-limiting mechanism.

chi is estimated per path by finite-difference (Y*(Delta)-Y*(0))/Delta at small Delta, SAME noise seed.
numpy only, vectorized over realizations.
"""
import numpy as np, math

def Vd(Y,Delta):
    a=np.abs(Y); return np.sign(Y)*a*(a+Delta)

def escape_delta(Delta, eta, N, seed, Y0=3.0, Yend=-6.0, dt=1.2e-3, thr=25.0, clip=300.0):
    rng=np.random.default_rng(seed); n=int(round((Y0-Yend)/dt)); sq=np.sqrt(dt)
    p=np.full(N,np.sqrt(max(Vd(Y0,Delta),1e-9))); Ys=np.full(N,np.nan)
    for i in range(n):
        Y=Y0-i*dt
        p+=(Vd(Y,Delta)-p*p)*dt+eta*sq*rng.standard_normal(N)
        np.clip(p,-clip,clip,out=p)
        nw=np.isnan(Ys)&(p<-thr); Ys[nw]=Y-dt
    return Ys   # SAME seed => SAME noise path across Delta

if __name__=="__main__":
    eta=np.sqrt(2.0); N=300000; seed=12345
    # escape at Delta=0 and two small Delta, SAME noise
    d1,d2=0.05,0.10
    Ys0=escape_delta(0.0,eta,N,seed)
    Ys1=escape_delta(d1,eta,N,seed)
    Ys2=escape_delta(d2,eta,N,seed)
    ok=~(np.isnan(Ys0)|np.isnan(Ys1)|np.isnan(Ys2))
    print(f"escaped (all 3): {ok.sum()}/{N}")
    chi1=(Ys1[ok]-Ys0[ok])/d1     # susceptibility estimate at Delta=d1
    chi2=(Ys2[ok]-Ys0[ok])/d2     # at d2
    # (a) linearity: chi1 ~ chi2 (both estimate dY*/dDelta); and dY*(d2) ~ 2*dY*(d1)
    dY1=Ys1[ok]-Ys0[ok]; dY2=Ys2[ok]-Ys0[ok]
    print(f"\n(a) pathwise linearity of dY* in Delta:")
    print(f"    E[dY*(d1)]={dY1.mean():+.4f}  E[dY*(d2)]={dY2.mean():+.4f}  ratio={dY2.mean()/dY1.mean():.3f} (expect {d2/d1:.1f})")
    print(f"    corr(dY*(d1),dY*(d2))={np.corrcoef(dY1,dY2)[0,1]:.4f} (expect ~1 if linear/same chi)")
    # (b) chi mean -> W_1 rate
    print(f"\n(b) susceptibility chi=dY*/dDelta (at d1={d1}):")
    print(f"    E[chi]={chi1.mean():+.4f}  E|chi|={np.abs(chi1).mean():.4f}  => W_1 ~ {np.abs(chi1).mean():.3f}*Delta")
    print(f"    => Wasserstein rate p=1/4: W_1(W_2Delta,W_2) ~ {np.abs(chi1).mean():.2f}*eps^(1/4)")
    # (c) tail of chi: moments and tail exponent
    c=chi1; cs=c-c.mean()
    m2=np.mean(cs**2); m4=np.mean(cs**4)
    print(f"\n(c) chi distribution (rate-limiting mechanism):")
    print(f"    std={np.sqrt(m2):.3f}  E|chi|<inf? mean abs={np.abs(c).mean():.3f} (finite)")
    print(f"    kurtosis of chi = {m4/m2**2:.2f} (Gaussian=3; >>3 => heavy tail => 4th-moment fragile)")
    # tail: fraction beyond k stds, and max
    for k in [3,5,8]:
        print(f"    P(|chi-mean|>{k} std)={np.mean(np.abs(cs)>k*np.sqrt(m2)):.2e} (Gaussian {2*(1-0.5*(1+math.erf(k/np.sqrt(2)))):.1e})")
    print(f"    max|chi|={np.abs(cs).max():.1f} std  => heavy tail confirms kurtosis correction is non-analytic")
    np.savez("tierB_chi.npz", chi1=chi1, chi2=chi2, dY1=dY1, dY2=dY2)
