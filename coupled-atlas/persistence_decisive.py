"""
Frontier item 2, decisive pass: is the persistence (left-tail) rate exponent gamma (rate ~ beta^gamma)
genuinely anomalous (~0.6) or an artifact of linear-fitting a curved -logP at low stats?

Strategy (non-confounded):
  1. High-stats escape sampling, dt-converged, generous thr/clip so the deep tail is not distorted.
  2. Measure the FUNCTIONAL FORM of the deep left tail: -log P(Theta>t) vs t. Pure exponential (constant
     hazard) => straight line => a well-defined 'rate'. Curved => 'rate' is ill-defined and any single gamma
     is fit-window-dependent (the artifact hypothesis).
  3. Extract beta-scaling of the rate over a WIDE beta range (2..32) with a consistent tail window in the
     SAME standardized position, and attempt a data collapse.
  4. Report standardized shape exponent a_L(beta) too.
scipy available.  Tags [NUMERIC].
"""
import numpy as np

def Vc(Y): return np.sign(Y)*Y*Y
def escape(eta, N, Y0=3.0, Yend=-9.0, dt=1.2e-3, seed=2, thr=20.0, clip=250.0):
    rng=np.random.default_rng(seed); n=int(round((Y0-Yend)/dt)); sq=np.sqrt(dt)
    p=np.full(N,np.sqrt(Vc(Y0))); Ys=np.full(N,np.nan)
    for i in range(n):
        Y=Y0-i*dt
        p+=(Vc(Y)-p*p)*dt+eta*sq*rng.standard_normal(N); np.clip(p,-clip,clip,out=p)
        nw=np.isnan(Ys)&(p<-thr); Ys[nw]=Y0-(i+1)*dt
    return Ys[~np.isnan(Ys)]

betas=[2.0,4.0,8.0,16.0,32.0]
N=600000
print("=== dt convergence check (beta=8, tail rate over Theta in [3,6]) ===")
for dt in [2.4e-3,1.2e-3,0.6e-3]:
    Ys=escape(2.0/np.sqrt(8.0),200000,dt=dt,seed=5); Th=(Ys[Ys<0]**2)/2
    tg=np.linspace(3.0,6.0,10); S=np.array([(Th>t).mean() for t in tg]); ok=S*len(Th)>=40
    r=np.polyfit(tg[ok],-np.log(S[ok]),1)[0]
    print(f"  dt={dt:.1e}: tail rate={r:.3f}  (n_tail@Theta>3 = {(Th>3).sum()})")

print("\n=== functional form + beta-scaling ===")
data={}
for b in betas:
    eta=2.0/np.sqrt(b); Ys=escape(eta,N,seed=int(b)+1); Ysn=Ys[Ys<0]; Th=(Ysn**2)/2
    data[b]=(Ys,Th)
    # functional form: fit -logP(Theta>t) to linear (rate) and quadratic; report curvature
    tg=np.linspace(2.6,6.5,22); S=np.array([(Th>t).mean() for t in tg]); ok=S*len(Th)>=30
    t2=tg[ok]; y2=-np.log(S[ok])
    lin=np.polyfit(t2,y2,1); quad=np.polyfit(t2,y2,2)
    r2lin=1-np.sum((y2-np.polyval(lin,t2))**2)/np.sum((y2-y2.mean())**2)
    curv=quad[0]  # coefficient of t^2; ~0 => straight (exponential tail)
    # standardized shape exponent a_L
    s=(Ys-Ys.mean())/Ys.std(); Sg=np.linspace(1.0,3.0,16)
    PL=np.array([(s<-S0).mean() for S0 in Sg]); okL=PL*len(s)>=30
    aL=np.polyfit(np.log(Sg[okL]),np.log(-np.log(PL[okL])),1)[0]
    print(f"  beta={b:>4.0f}: rate(lin)={lin[0]:.3f} (R2={r2lin:.4f}) curv={curv:+.4f} | a_L={aL:.2f} | "
          f"mean(Y*)={Ys.mean():+.3f} std={Ys.std():.3f} n_tail(Th>4)={(Th>4).sum()}")

rates=np.array([np.polyfit(np.linspace(2.6,6.5,22)[ (np.array([(data[b][1]>t).mean() for t in np.linspace(2.6,6.5,22)])*len(data[b][1]))>=30 ],
        -np.log(np.array([(data[b][1]>t).mean() for t in np.linspace(2.6,6.5,22)])[ (np.array([(data[b][1]>t).mean() for t in np.linspace(2.6,6.5,22)])*len(data[b][1]))>=30 ]),1)[0] for b in betas])
g=np.polyfit(np.log(betas),np.log(rates),1)[0]
print(f"\n  rate vs beta:  {[f'{r:.2f}' for r in rates]}")
print(f"  gamma (rate~beta^gamma) = {g:.3f}   [gamma=1 drift-diffusion | gamma~0.6 anomalous]")
np.savez("persistence_data.npz", betas=np.array(betas), rates=rates, gamma=g,
         **{f"Th_{int(b)}":data[b][1] for b in betas}, **{f"Ys_{int(b)}":data[b][0] for b in betas})
print("saved persistence_data.npz")
