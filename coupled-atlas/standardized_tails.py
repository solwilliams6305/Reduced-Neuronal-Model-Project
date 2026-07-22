"""
Clean measurement of W's standardized tails (left=persistence/late, right=Kramers/early), eta-dependence.
A persistence/Gamow universal exponent would be eta-INDEPENDENT in the standardized variable.
"""
import numpy as np
def Vc(Y): return np.sign(Y)*Y*Y
def escape(eta,N,Y0=3.0,Yend=-8.0,dt=6e-3,seed=2):
    rng=np.random.default_rng(seed); n=int(round((Y0-Yend)/dt)); sq=np.sqrt(dt)
    p=np.full(N,np.sqrt(Vc(Y0))); Ys=np.full(N,np.nan)
    for i in range(n):
        Y=Y0-i*dt; p+=(Vc(Y)-p*p)*dt+eta*sq*rng.standard_normal(N); np.clip(p,-60,60,out=p)
        nw=np.isnan(Ys)&(p<-12.0); Ys[nw]=Y0-(i+1)*dt
    return Ys[~np.isnan(Ys)]
def tailexp(s,side):
    S=np.linspace(1.0,3.0,16)
    P=np.array([( (s<-S0).mean() if side=='L' else (s>S0).mean()) for S0 in S]); ok=P*len(s)>=30
    # local exponent: -logP ~ S^a  => slope of log(-logP) vs log S
    nlp=-np.log(P[ok]); a=np.polyfit(np.log(S[ok]),np.log(nlp),1)[0]; return a
for eta in [0.7,1.0,1.4]:
    Ys=escape(eta,250000); s=(Ys-Ys.mean())/Ys.std()
    aL=tailexp(s,'L'); aR=tailexp(s,'R'); sk=np.mean(s**3)
    print(f"eta={eta} (beta={4/eta**2:.1f}): skew={sk:+.3f}  LEFT exp={aL:.2f}  RIGHT exp={aR:.2f}  "
          f"({'L lighter' if aL>aR else 'R lighter'})")
print("\nGamow/persistence universal => eta-INDEPENDENT exponents. Kramers right => eta^2 LDP (heavier, lower exp).")
