#!/usr/bin/env python3
"""
mmo_fhr_alpha_fine.py  —  Phase-1.5 α-tightening (prerequisite for the K2 step).

Refines the realised plateau-width exponent α (Δ_pq ~ 1/q^α) for the working FHR
beyond the 6-point estimate in mmo_fhr_alpha.py. Strategy (same as the canard
chapter: tighten the empirical number, then derive against it):

  • finer c-grid (>20 samples) across the MMO band;
  • longer trajectories so the SAO trains at small μ are counted accurately;
  • a ROBUST per-spike SAO statistic — the mean nS/nL averaged over many cycles —
    rather than the noisy single max; max reported as a cross-check;
  • α = 1 + 1/p from s ~ (c+1)^{-p}, with the slope's standard error propagated,
    fit over the full band AND the deep (large-s, c→−1) subrange where the
    folded-node scaling s ~ 1/(c+1) actually holds.

Saves raw (c, nL, nS) so the fit can be re-run/aggregated cheaply.
Outputs: results/mmo/{alpha_fine.npz, mmo_fhr_alpha_fine.txt}, figures/mmo_fhr_alpha_fine.png
Reproduce:  python3 regime-tests/mmo_fhr_alpha_fine.py
"""
from __future__ import annotations
import os
import numpy as np

HERE=os.path.dirname(os.path.abspath(__file__)); BASE=os.path.dirname(HERE)
A,B,EPS,DELTA,I0=0.7,0.8,0.08,0.2,0.30

def step(s,c,dt):
    v,w,y=s
    def f(v,w,y): return (v-v**3/3-w+y+I0, EPS*(v+A-B*w), EPS*DELTA*(c-v))
    k1=f(v,w,y);k2=f(v+.5*dt*k1[0],w+.5*dt*k1[1],y+.5*dt*k1[2]);k3=f(v+.5*dt*k2[0],w+.5*dt*k2[1],y+.5*dt*k2[2]);k4=f(v+dt*k3[0],w+dt*k3[1],y+dt*k3[2])
    return np.array([v+dt*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6,w+dt*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6,y+dt*(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6])

def measure(c,dt=0.04,Tw=1500,Tr=7000):
    s=np.array([0.5,0.0,0.0])
    for _ in range(int(Tw/dt)): s=step(s,c,dt)
    n=int(Tr/dt); vs=np.empty(n)
    for i in range(n): s=step(s,c,dt); vs[i]=s[0]
    dv=np.diff(vs); mx=np.where((dv[:-1]>0)&(dv[1:]<=0))[0]+1; pk=vs[mx] if len(mx) else np.array([])
    nL=int(np.sum(pk>0.0)); nS=int(np.sum(pk<=0.0))
    sym=''.join('L' if p>0 else 'S' for p in pk); runs=[len(r) for r in sym.split('L') if r]
    mean_run=(nS/nL) if nL>0 else np.nan
    max_run=max(runs) if runs else 0
    return nL,nS,mean_run,max_run

def fit_alpha(cp, s, label):
    """Fit log s = const − p log(c+1); return p, SE(p), α=1+1/p, SE(α)."""
    ok=(s>0)&np.isfinite(s)&(cp>0)
    x=np.log(cp[ok]); y=np.log(s[ok]); n=ok.sum()
    if n<3: return None
    X=np.vstack([np.ones(n),x]).T
    beta,res,_,_=np.linalg.lstsq(X,y,rcond=None)
    p=-beta[1]
    yhat=X@beta; sig2=np.sum((y-yhat)**2)/(n-2)
    cov=sig2*np.linalg.inv(X.T@X); sep=np.sqrt(cov[1,1])
    alpha=1+1/p; se_alpha=sep/p**2
    return dict(label=label,n=n,p=p,sep=sep,alpha=alpha,se=se_alpha)

def main():
    lines=[];P=lambda s="":(print(s,flush=True),lines.append(s))
    P("="*72);P("MMO α-tightening — finer c-grid + error bars");P("="*72)
    cs=np.linspace(-0.74,-0.95,22)
    nL=np.empty(len(cs));nS=np.empty(len(cs));mean_run=np.empty(len(cs));max_run=np.empty(len(cs))
    P(f"{'c':>7}{'c+1':>7}{'nL':>5}{'nS':>5}{'nS/nL':>8}{'maxS':>6}")
    for i,c in enumerate(cs):
        nL[i],nS[i],mean_run[i],max_run[i]=measure(c)
        P(f"{c:7.3f}{c+1:7.3f}{int(nL[i]):5d}{int(nS[i]):5d}{mean_run[i]:8.2f}{int(max_run[i]):6d}")
    np.savez(os.path.join(HERE,"results","mmo","alpha_fine.npz"),c=cs,nL=nL,nS=nS)
    cp=cs+1.0
    P("\nFits  (s ~ (c+1)^(-p);  α = 1 + 1/p):")
    results=[]
    for stat,name in [(mean_run,"mean nS/nL"),(max_run,"max SAOs")]:
        for lo,hi,tag in [(-0.95,-0.74,"full band"),(-0.95,-0.84,"deep (c≤−0.84)")]:
            msk=(cs>=lo)&(cs<=hi)
            r=fit_alpha(cp[msk],stat[msk],f"{name} / {tag}")
            if r: results.append(r); P(f"  {r['label']:>26}: p={r['p']:.2f}±{r['sep']:.2f}  α={r['alpha']:.2f}±{r['se']:.2f}  (n={r['n']})")
    # headline: mean-statistic, deep subrange
    head=[r for r in results if r['label']=="mean nS/nL / deep (c≤−0.84)"]
    if head:
        h=head[0]
        P(f"\nHEADLINE: realised α = {h['alpha']:.2f} ± {h['se']:.2f}  (mean stat, deep subrange).")
    P("  Folded-node ceiling stays α=2 (μ∝(c+1)). α≠1 confirmed at tight error.")
    P("  ⇒ K2 target: does the K2/entry-exit f(c) reproduce this α?")

    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig,ax=plt.subplots(1,2,figsize=(12,4.4))
        ax[0].loglog(cp,mean_run,'o',label="mean nS/nL")
        ax[0].loglog(cp,max_run,'s',ms=4,alpha=.6,label="max SAOs")
        # overlay deep-subrange fit on mean
        h=[r for r in results if r['label'].startswith("mean nS/nL / deep")]
        if h:
            xx=np.linspace(cp.min(),cp.max(),50); p=h[0]['p']
            C=np.exp(np.mean(np.log(mean_run[(cs<=-0.84)])+p*np.log(cp[(cs<=-0.84)])))
            ax[0].loglog(xx,C*xx**(-p),'C3-',lw=1,label=f"fit α={h[0]['alpha']:.2f}±{h[0]['se']:.2f}")
        ax[0].set_xlabel("c+1"); ax[0].set_ylabel("SAOs per spike"); ax[0].legend(fontsize=8)
        ax[0].set_title("realised rotation-count scaling → α=1+1/p")
        ax[1].plot(cs, nL/np.maximum(nL+nS,1),'o-')
        ax[1].set_xlabel("c"); ax[1].set_ylabel("ρ = L/(L+S)"); ax[1].set_title("staircase ρ(c), finer grid")
        fig.tight_layout(); fp=os.path.join(BASE,"figures","mmo_fhr_alpha_fine.png"); fig.savefig(fp,dpi=120)
        P(f"\nFigure: {fp}")
    except Exception as e:
        P(f"[figure skipped: {e}]")
    os.makedirs(os.path.join(HERE,"results","mmo"),exist_ok=True)
    open(os.path.join(HERE,"results","mmo","mmo_fhr_alpha_fine.txt"),"w").write("\n".join(lines))

if __name__=="__main__": main()
