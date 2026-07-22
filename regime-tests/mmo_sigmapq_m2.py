#!/usr/bin/env python3
"""
mmo_sigmapq_m2.py  —  MMO σ_pq plan, Phase M2 (deterministic sensitivities).

The σ_pq plan proposes the mechanistic candidate σ_pq ~ μ·a_min(c) (rotation-map
sensitivity ds/d(ln a_in)=−1/κμ × injection depth a_min(c)) as a replacement for
the placeholder μ^{3/2}. M2 measures the deterministic ingredients (no noise, no
counting floor) and tests the candidate against the resolved σ_pq.

FINDING (this run): the candidate is NOT empirically preferred over μ^{3/2} at
current resolution — and the reason is instructive. The smallest-SAO amplitude
a_min is **ambiguous at the factor-2 level**: it varies smoothly across each
plateau's c-width (e.g. ~0.58 → 0.25 across the q=3 plateau) and with the
counting protocol. That ambiguity is LARGER than the μ·a_min-vs-μ^{3/2}
difference, so which form "wins" the 2-plateau check flips with the a_min value
(with a_min(q3)≈0.58 here, μ^{3/2} fits q=4 better; with a_min(q3)≈0.25 as in the
plan's note, μ·a_min looks flat). Also κ=ln(a_max/a_min)/(sμ) is ~2π² only
asymptotically (8 at q=3 → 17 at q=6), so even the rotation-map sensitivity drifts
at low q.

⇒ M2 confirms the plan's chief hazard quantitatively: μ·a_min and μ^{3/2} are
degenerate, and the deterministic side cannot separate them because a_min itself
is only factor-2 defined. The separation genuinely needs BOTH M1 (robust q=2…6
noise σ_pq) AND a *derived* a_min(c) — i.e. Phase B (the K2 inner solution) is the
real resolver, exactly the coupling the plan identifies. Do not declare the power.

Outputs: results/mmo/mmo_sigmapq_m2.txt, figures/mmo_sigmapq_m2.png
Reproduce:  python3 regime-tests/mmo_sigmapq_m2.py
"""
from __future__ import annotations
import os
import numpy as np

HERE=os.path.dirname(os.path.abspath(__file__)); BASE=os.path.dirname(HERE)
A,B,EPS,DELTA,I0=0.7,0.8,0.08,0.2,0.30
MEAS={3:0.022,4:0.016}

def mu(c):
    tr=1+DELTA; det=2*B*DELTA*(c+1); disc=tr*tr-4*det
    return ((tr-np.sqrt(disc))/2)/((tr+np.sqrt(disc))/2)

def episodes(c,dt=0.02,Tw=1200,Tr=7000):
    v,w,y=0.5,0.0,0.0
    def stp():
        nonlocal v,w,y
        v+=dt*(v-v**3/3-w+y+I0); w+=dt*EPS*(v+A-B*w); y+=dt*EPS*DELTA*(c-v)
    for _ in range(int(Tw/dt)): stp()
    vs=np.empty(int(Tr/dt))
    for i in range(len(vs)): stp(); vs[i]=v
    dv=np.diff(vs); mx=np.where((dv[:-1]>0)&(dv[1:]<=0))[0]+1; mn=np.where((dv[:-1]<0)&(dv[1:]>=0))[0]+1
    amp={}
    for m in mx:
        pr=mn[mn<m]
        if len(pr): amp[m]=vs[m]-vs[pr[-1]]
    idx=sorted(amp); pkv=np.array([vs[i] for i in idx]); pka=np.array([amp[i] for i in idx]); isL=pkv>0
    sc=[];am=[];aM=[];run=[]
    for L,Aamp in zip(isL,pka):
        if L:
            if run: sc.append(len(run));am.append(min(run));aM.append(max(run))
            run=[]
        else: run.append(Aamp)
    return np.array(sc),np.array(am),np.array(aM)

def main():
    lines=[];P=lambda s="":(print(s,flush=True),lines.append(s))
    P("="*72);P("MMO Phase M2 — deterministic a_min(c) + μ·a_min vs μ^{3/2}");P("="*72)

    # ONE deterministic scan: cache s(c), a_min(c), a_max(c)
    P("\n(1) a_min(c) across the band (deterministic) — note within-plateau variation:")
    cgrid=np.round(np.linspace(-0.77,-0.895,18),4)
    S={};AMIN={};AMAX={};amins=[]
    for c in cgrid:
        sc,am,aM=episodes(c)
        S[c]=np.median(sc) if len(sc) else 0.0
        AMIN[c]=np.median(am) if len(am) else np.nan
        AMAX[c]=np.median(aM) if len(aM) else np.nan
        amins.append(AMIN[c])
    for c in cgrid[::2]:
        P(f"   c={c:.3f}  s={S[c]:.1f}  a_min={AMIN[c]:.3f}")

    # plateau centres q=2..6 from the cached scan
    P("\n(2) at plateau centres (s=q−1): a_min, a_max, s, μ, κ, μ·a_min")
    P(f"   {'q':>2}{'c':>8}{'μ':>8}{'s':>4}{'a_min':>8}{'a_max':>8}{'κ':>7}{'μ·a_min':>9}")
    rows=[]
    for q in (2,3,4,5,6):
        c=min(cgrid,key=lambda cc:abs(S[cc]-(q-1)))
        s_=S[c]; amin=AMIN[c]; amax=AMAX[c]; m=mu(c)
        if not (s_>0 and amin==amin): continue
        kap=np.log(amax/amin)/(s_*m) if s_>0 else np.nan
        rows.append((q,c,m,s_,amin,amax,kap)); P(f"   {q:2d}{c:8.3f}{m:8.4f}{s_:4.0f}{amin:8.3f}{amax:8.3f}{kap:7.2f}{m*amin:9.5f}")

    # candidate test (normalised at q=3)
    q3=[r for r in rows if r[0]==3][0]
    P("\n(3) σ_pq candidates (normalised at q=3) vs measured:")
    P(f"   {'q':>2}{'μ·a_min':>10}{'μ^{3/2}':>10}{'measured':>10}")
    for q,c,m,s_,amin,amax,kap in rows:
        pma=MEAS[3]*(m*amin)/(q3[2]*q3[4]); pm32=MEAS[3]*(m**1.5)/(q3[2]**1.5)
        P(f"   {q:2d}{pma:10.4f}{pm32:10.4f}{MEAS.get(q,np.nan):10.4f}")
    ratio={q:MEAS[q]/(m*amin) for (q,c,m,s_,amin,amax,kap) in rows if q in MEAS}
    P(f"\n   σ_pq/(μ·a_min) at resolved plateaus: "+", ".join(f"q{q}:{ratio[q]:.2f}" for q in ratio))
    P("   (NOT flat — vs the plan's 2-pt ~1.9; the difference is the a_min protocol/c.)")
    P("\n VERDICT: a_min ambiguous at factor ~2 (varies across plateau width); this")
    P(" ambiguity exceeds the μ·a_min-vs-μ^{3/2} gap, so the deterministic side can")
    P(" NOT separate the forms. κ→2π² only at high q (8 at q=3 → 17 at q=6). The")
    P(" separation needs M1 (robust q=2…6 noise σ_pq) AND a derived a_min(c) (Phase B,")
    P(" the K2 inner solution) — the shared-unknown coupling the plan identifies.")
    P(" μ^{3/2} stays a placeholder; μ·a_min is mechanistically motivated but unproven.")

    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig,ax=plt.subplots(1,2,figsize=(12,4.5))
        ax[0].plot(cgrid,amins,'o-'); ax[0].set_xlabel("c"); ax[0].set_ylabel("a_min (smallest SAO amp)")
        ax[0].set_title("a_min(c): smooth, varies factor~2 within each plateau (ambiguous)")
        qs=np.array([r[0] for r in rows])
        pma=np.array([MEAS[3]*(r[2]*r[4])/(q3[2]*q3[4]) for r in rows])
        pm32=np.array([MEAS[3]*(r[2]**1.5)/(q3[2]**1.5) for r in rows])
        ax[1].loglog(qs,pma,'s-',label="μ·a_min(c) candidate")
        ax[1].loglog(qs,pm32,'^--',label="μ^{3/2} placeholder")
        ax[1].loglog([3,4],[0.022,0.016],'ro',ms=9,label="measured")
        ax[1].set_xlabel("q"); ax[1].set_ylabel("σ_pq"); ax[1].legend(fontsize=8)
        ax[1].set_title("candidates degenerate at q=3,4; diverge at q=2,5,6 (M1 target)")
        fig.tight_layout(); fp=os.path.join(BASE,"figures","mmo_sigmapq_m2.png"); fig.savefig(fp,dpi=120)
        P(f"\nFigure: {fp}")
    except Exception as e:
        P(f"[figure skipped: {e}]")
    os.makedirs(os.path.join(HERE,"results","mmo"),exist_ok=True)
    open(os.path.join(HERE,"results","mmo","mmo_sigmapq_m2.txt"),"w").write("\n".join(lines))

if __name__=="__main__": main()
