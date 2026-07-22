"""
A+B decisive check: is W's late-escape (LEFT) tail a Gamow/persistence exponent?
Escape = first node of u; depth = phase Theta = Y*^2/2. Survival S(Theta)=P(Theta* > Theta) = persistence.
Gamow predicts S ~ e^{-2|Im l0| Theta} = e^{-1.64 Theta}, eta-INDEPENDENT (property of the deterministic resonance).
Also: dominant-resonance gives a Rayleigh-in-|Y*| bulk (skew sign opposite to W) => Gamow gives the TAIL not the bulk.
numpy only.
"""
import numpy as np
def Vc(Y): return np.sign(Y)*Y*Y
def escape(eta,N,Y0=3.0,Yend=-8.0,dt=6e-3,seed=1):
    rng=np.random.default_rng(seed); n=int(round((Y0-Yend)/dt)); sq=np.sqrt(dt)
    p=np.full(N,np.sqrt(Vc(Y0))); Ys=np.full(N,np.nan)
    for i in range(n):
        Y=Y0-i*dt; p+=(Vc(Y)-p*p)*dt+eta*sq*rng.standard_normal(N); np.clip(p,-60,60,out=p)
        nw=np.isnan(Ys)&(p<-12.0); Ys[nw]=Y0-(i+1)*dt
    return Ys[~np.isnan(Ys)]
print("LEFT-tail survival S(Theta)=P(Y*^2/2 > Theta); Gamow/persistence predicts slope 2|Im l0|=1.64, eta-indep")
for eta in [0.7,1.0]:
    Ys=escape(eta,200000); Ys=Ys[Ys<0]; Th=Ys*Ys/2.0
    z=(Ys-Ys.mean())/Ys.std(); sk=np.mean(z**3); ek=np.mean(z**4)-3
    Tg=np.linspace(3.0,9.0,16); S=np.array([(Th>t).mean() for t in Tg]); ok=S*len(Th)>=25
    slope=np.polyfit(Tg[ok],-np.log(S[ok]),1)[0]
    print(f"  eta={eta}: W skew={sk:+.3f} exk={ek:+.3f} (FP +0.60/-0.24); LEFT-tail rate in Theta = {slope:.2f}  (Gamow 1.64)")
# dominant-resonance Gamow bulk prediction: f(Y*) ~ |Y*| e^{-1.64 Y*^2/2}=|Y*|e^{-0.82 Y*^2} (Rayleigh in R=|Y*|)
R=np.linspace(0.01,6,6000); fR=R*np.exp(-0.82*R*R); fR/=np.trapz(fR,R)
mR=np.trapz(R*fR,R); vR=np.trapz((R-mR)**2*fR,R); skR=np.trapz((R-mR)**3*fR,R)/vR**1.5
print(f"\ndominant-resonance Gamow bulk = Rayleigh(|Y*|): skew(|Y*|)={skR:+.3f} => skew(Y*=-|Y*|)={-skR:+.3f}")
print("  W bulk skew is +0.60 (opposite sign) => Gamow survival gives the LEFT TAIL, not the bulk (bulk +skew = right/Kramers side).")
