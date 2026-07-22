"""
Decides Piece-2 tractability: how does the persistence (left-tail survival) rate scale with beta?
rate(Theta=Y*^2/2) ~ beta^gamma. gamma=1 => drift-diffusion first-passage (closed-form-able).
gamma=1/2 (or other) => anomalous => needs genuine non-stationary persistence theory.
"""
import numpy as np
def Vc(Y): return np.sign(Y)*Y*Y
def escape(eta,N,Y0=3.0,Yend=-8.0,dt=5e-3,seed=3):
    rng=np.random.default_rng(seed); n=int(round((Y0-Yend)/dt)); sq=np.sqrt(dt)
    p=np.full(N,np.sqrt(Vc(Y0))); Ys=np.full(N,np.nan)
    for i in range(n):
        Y=Y0-i*dt; p+=(Vc(Y)-p*p)*dt+eta*sq*rng.standard_normal(N); np.clip(p,-60,60,out=p)
        nw=np.isnan(Ys)&(p<-12.0); Ys[nw]=Y0-(i+1)*dt
    return Ys[~np.isnan(Ys)]
betas=[2.0,4.0,8.0,16.0]; rates=[]
for b in betas:
    eta=2.0/np.sqrt(b); Ys=escape(eta,250000); Ys=Ys[Ys<0]; Th=Ys*Ys/2
    Tg=np.linspace(3.2,7.0,14); S=np.array([(Th>t).mean() for t in Tg]); ok=S*len(Th)>=40
    r=np.polyfit(Tg[ok],-np.log(S[ok]),1)[0]; rates.append(r)
    print(f"  beta={b:>4.0f} (eta={eta:.2f}): left-tail survival rate in Theta = {r:.2f}")
g=np.polyfit(np.log(betas),np.log(rates),1)[0]
print(f"\nrate ~ beta^{g:.2f}   =>  {'drift-diffusion (gamma~1, Piece-2 closed-form-able)' if abs(g-1)<0.2 else 'ANOMALOUS (gamma!=1) => needs genuine non-stationary persistence theory'}")
