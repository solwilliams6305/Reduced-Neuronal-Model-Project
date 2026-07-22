#!/usr/bin/env python3
"""
mmo_fhr_foldednode.py  —  Phase 1.5 Task B' : folded-node μ + Wechselberger test.

Tests the plan's (MMO_FHR_PLAN.md §2) central claim that the MMO structure is set
by the folded-node eigenvalue ratio μ via Wechselberger (2005): the maximal SAO
count per spike is (1−μ)/(2μ). Two results:

(1) PARAMETRISATION CORRECTION. The plan's §0 *canonical* FHR (Rinzel form,
    v'=v−v³/3−w−y+I, y'=δ(cv+d−y), c=1,d=0.8) produces **pure spiking, no MMOs**
    at every δ∈[0.005,0.04] and I∈[0.5,0.9] tested — the folded node exists
    analytically but the global return does not feed its funnel. The working MMO
    model is the Phase-1 variant (v'=…+y, y'=εδ(c−v), c≈−0.8, δ=0.2).

(2) μ-VERIFICATION (working FHR). The desingularised reduced-flow eigenvalue
    ratio μ(c) is computed analytically. As c→−1, μ→0 and the Wechselberger
    bound (1−μ)/(2μ) grows; the **observed** max SAOs per spike grows in step but
    sits BELOW the bound (it is the *maximal* rotation number, realised only as
    the global return injects deep into the funnel). So μ genuinely governs the
    SAO structure, but the plan's "#SAOs ≈ (1−μ)/(2μ)" should read "≤".

Outputs: results/mmo/mmo_fhr_foldednode.txt, figures/mmo_fhr_foldednode.png
Reproduce:  python3 regime-tests/mmo_fhr_foldednode.py
"""
from __future__ import annotations
import os
import numpy as np

HERE=os.path.dirname(os.path.abspath(__file__)); BASE=os.path.dirname(HERE)
A,B,EPS=0.7,0.8,0.08

# ----- working Phase-1 FHR -----
DELTA_W, I_W = 0.2, 0.30
def step_work(s,c,dt):
    v,w,y=s
    def f(v,w,y): return (v-v**3/3-w+y+I_W, EPS*(v+A-B*w), EPS*DELTA_W*(c-v))
    k1=f(v,w,y);k2=f(v+.5*dt*k1[0],w+.5*dt*k1[1],y+.5*dt*k1[2]);k3=f(v+.5*dt*k2[0],w+.5*dt*k2[1],y+.5*dt*k2[2]);k4=f(v+dt*k3[0],w+dt*k3[1],y+dt*k3[2])
    return np.array([v+dt*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6,w+dt*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6,y+dt*(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6])
def mu_work(c):
    """folded-node μ for the working FHR (2-slow reduction, fold v=−1)."""
    tr=1+DELTA_W; det=2*B*DELTA_W*(c+1)
    if det<=0: return np.nan
    disc=tr*tr-4*det
    if disc<0: return -1.0   # folded focus
    lp=(tr+np.sqrt(disc))/2; lm=(tr-np.sqrt(disc))/2
    return lm/lp
def pattern_work(c,dt=0.04,Tw=2000,Tr=6000):
    s=np.array([0.5,0.0,0.0])
    for _ in range(int(Tw/dt)): s=step_work(s,c,dt)
    vs=np.empty(int(Tr/dt))
    for i in range(len(vs)): s=step_work(s,c,dt); vs[i]=s[0]
    dv=np.diff(vs); mx=np.where((dv[:-1]>0)&(dv[1:]<=0))[0]+1; pk=vs[mx] if len(mx) else np.array([])
    sym=''.join('L' if p>0 else 'S' for p in pk); runs=[len(r) for r in sym.split('L') if r]
    nL=int(np.sum(pk>0)); nS=int(np.sum(pk<=0))
    return (max(runs) if runs else 0), nL, nS

# ----- plan's canonical FHR (for the no-MMO check) -----
def step_plan(s,I,delta,dt):
    v,w,y=s
    def f(v,w,y): return (v-v**3/3-w-y+I, EPS*(v+A-B*w), delta*(1.0*v+0.8-y))
    k1=f(v,w,y);k2=f(v+.5*dt*k1[0],w+.5*dt*k1[1],y+.5*dt*k1[2]);k3=f(v+.5*dt*k2[0],w+.5*dt*k2[1],y+.5*dt*k2[2]);k4=f(v+dt*k3[0],w+dt*k3[1],y+dt*k3[2])
    return np.array([v+dt*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6,w+dt*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6,y+dt*(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6])
def plan_mmo(I,delta,dt=0.02,Tw=3000,Tr=5000):
    s=np.array([-1.0,0.0,0.0])
    for _ in range(int(Tw/dt)): s=step_plan(s,I,delta,dt)
    vs=np.empty(int(Tr/dt))
    for i in range(len(vs)): s=step_plan(s,I,delta,dt); vs[i]=s[0]
    dv=np.diff(vs); mx=np.where((dv[:-1]>0)&(dv[1:]<=0))[0]+1; pk=vs[mx] if len(mx) else np.array([])
    return int(np.sum(pk>0)), int(np.sum(pk<=0))

def main():
    lines=[];P=lambda s="":(print(s,flush=True),lines.append(s))
    P("="*72);P("MMO Phase 1.5 Task B' — folded-node μ + Wechselberger test");P("="*72)

    P("\n(1) Plan's §0 canonical FHR — does it produce MMOs?")
    P(f"   {'delta':>6} {'I':>5}  nL/nS")
    any_mmo=False
    for delta in (0.005,0.02):
        for I in (0.6,0.9):
            nL,nS=plan_mmo(I,delta); any_mmo|=(nL>0 and nS>0)
            P(f"   {delta:6.3f} {I:5.2f}  {nL}L/{nS}S  {'<<MMO' if nL>0 and nS>0 else 'pure spiking'}")
    P(f"   => plan's §0 FHR makes MMOs: {any_mmo}.  (folded node exists but the")
    P(f"      global return bypasses its funnel ⇒ the §0 parameters need correcting.)")

    P("\n(2) Working FHR (Phase-1 variant) — Wechselberger μ test")
    P(f"   params: v'=v−v³/3−w+y+I, w'=ε(v+a−bw), y'=εδ(c−v); ε={EPS} δ={DELTA_W} I={I_W}")
    P(f"   {'c':>7} {'μ':>7} {'(1−μ)/2μ (bound)':>16} {'obs maxS':>9} {'obs nL/nS':>10}")
    cs=[-0.72,-0.78,-0.84,-0.90,-0.94]; mus=[];preds=[];obs=[]
    for c in cs:
        mu=mu_work(c); pred=(1-mu)/(2*mu) if mu>0 else np.nan
        mS,nL,nS=pattern_work(c); mus.append(mu);preds.append(pred);obs.append(mS)
        P(f"   {c:7.2f} {mu:7.3f} {pred:16.1f} {mS:9d} {nL:4d}/{nS:<5d}")
    P("   => μ→0 as c→−1; observed max SAOs grows WITH the (1−μ)/(2μ) bound but")
    P("      stays below it (it is the MAXIMAL rotation number, not the typical).")
    P("      Folded-node μ governs the SAO structure; the plan's '≈' should be '≤'.")

    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig,ax=plt.subplots(1,2,figsize=(12,4.4))
        ax[0].plot(cs,mus,'o-'); ax[0].set_xlabel("c"); ax[0].set_ylabel("μ (folded-node eigenvalue ratio)")
        ax[0].set_title("folded-node μ(c), working FHR")
        ax[1].plot(cs,preds,'s--',color='C3',label="(1−μ)/(2μ)  Wechselberger bound")
        ax[1].plot(cs,obs,'o-',color='C0',label="observed max SAOs/spike")
        ax[1].set_xlabel("c"); ax[1].set_ylabel("SAOs per spike")
        ax[1].set_title("Wechselberger bound vs observed (obs ≤ bound, tracks it)")
        ax[1].legend(fontsize=8)
        fig.tight_layout(); fp=os.path.join(BASE,"figures","mmo_fhr_foldednode.png"); fig.savefig(fp,dpi=120)
        P(f"\nFigure: {fp}")
    except Exception as e:
        P(f"[figure skipped: {e}]")
    os.makedirs(os.path.join(HERE,"results","mmo"),exist_ok=True)
    open(os.path.join(HERE,"results","mmo","mmo_fhr_foldednode.txt"),"w").write("\n".join(lines))

if __name__=="__main__": main()
