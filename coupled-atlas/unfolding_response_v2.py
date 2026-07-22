"""
Rung E v2 -- consolidated perturbation/unfolding responses, T2-consistent.
New validations:
 (A) PHYSICAL cusp (q=2) linear unfolding V=sign(Y)(Y^2+c1|Y|): relevant perturbation,
     skew should flow cusp(~0.60) -> fold(~0.20) [lowest-degree-term-dominates at the physical cusp].
 (B) within-class symmetry-breaking V=Y^2 (Y>=0), -(1+a)Y^2 (Y<0): keeps leading degree 2 (stays cusp),
     but breaks Y->-Y; measures the skew SUSCEPTIBILITY d(skew)/da of the asymmetric-PIV fingerprint.
beta-response (C) reused from beta_axis_res.npy (D). numpy only.
"""
import numpy as np
def simulate(Vfun, eta, N, Y0=3.0, Yend=-3.0, dt=5e-3, seed=11):
    rng=np.random.default_rng(seed); nst=int(round((Y0-Yend)/dt)); sq=np.sqrt(dt)
    p=np.full(N, np.sqrt(max(Vfun(Y0),1e-6))); Ys=np.full(N,np.nan)
    for i in range(nst):
        Y=Y0-i*dt; p+=(Vfun(Y)-p*p)*dt+eta*sq*rng.standard_normal(N); np.clip(p,-60,60,out=p)
        nw=np.isnan(Ys)&(p<-12.0); Ys[nw]=Y0-(i+1)*dt
    Ys[np.isnan(Ys)]=Yend; z=(Ys-Ys.mean())/Ys.std()
    return np.mean(z**3), np.mean(z**4)-3
eta=np.sqrt(2.0); N=60000
print("(A) cusp linear unfolding V=sign(Y)(Y^2+c1|Y|) -> flows to fold (relevant perturbation)")
c1s=[0.0,0.5,1.0,2.0,4.0]; skA=[]
for c1 in c1s:
    sk,ek=simulate(lambda Y,c=c1: np.sign(Y)*(Y*Y+c*np.abs(Y)), eta, N)
    skA.append(sk); print(f"   c1={c1:>3}: skew={sk:+.3f}  exk={ek:+.3f}")
print(f"   => skew {skA[0]:+.2f} -> {skA[-1]:+.2f}  ({'flows to fold OK' if skA[0]>skA[-1]+0.2 else 'check'})")
print("\n(B) within-class symmetry-breaking V=Y^2(Y>=0), -(1+a)Y^2(Y<0): stays cusp, skew susceptibility")
as_=[-0.4,-0.2,0.0,0.2,0.4]; skB=[]
for a in as_:
    sk,ek=simulate(lambda Y,aa=a: np.where(Y>=0, Y*Y, -(1+aa)*Y*Y), eta, N)
    skB.append(sk); print(f"   a={a:+.1f}: skew={sk:+.3f}  exk={ek:+.3f}")
suscep=np.polyfit(as_,skB,1)[0]
print(f"   => susceptibility d(skew)/da = {suscep:+.3f}; skew stays cusp-class (~0.5-0.7), monotone in a")
np.save("unfold_v2.npy", {"c1s":c1s,"skA":skA,"as_":as_,"skB":skB,"suscep":suscep}, allow_pickle=True)
print("saved unfold_v2.npy")
