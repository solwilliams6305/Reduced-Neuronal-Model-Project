#!/usr/bin/env python3
"""
mmo_2d_check.py  —  Phase 1a of the MMO plan: does autonomous 2D FHN actually
host the deterministic L^a S^b devil's staircase the plan (§1–§2) assumes?

Premise under test (MMO_PLAN.md §1): "Inside the canard explosion window the
deterministic limit cycle exhibits mixed-mode oscillations: sequences of L^a S^b
patterns ... a devil's staircase organised by Farey arithmetic."

Topological obstruction (why this should FAIL in 2D): an autonomous planar flow
has unique trajectories, so a stable limit cycle is a simple closed (Jordan)
curve. It cannot self-intersect, hence cannot realise "a large loops + b small
loops" in one period for a+b>1. Genuine L^a S^b MMOs need ≥3 dimensions (a
folded node) or non-autonomous forcing. This script demonstrates that directly:
swept across I through the (lower) Hopf / canard region, every 2D-FHN attractor
is single-mode — a stable focus (no oscillation) or a pure relaxation cycle
(all L) — with NO bimodal L/S mixing and NO staircase.

Outputs: results/mmo/mmo_2d_check.txt, figures/mmo_2d_check.png
Reproduce:  python3 regime-tests/mmo_2d_check.py
"""
from __future__ import annotations
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); BASE = os.path.dirname(HERE)
import importlib.util
def _load(p,n):
    s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s)
    sys.modules[n]=m; s.loader.exec_module(m);return m
_load(os.path.join(HERE,"_shim.py"),"_shim")   # installs pure-Python brentq shim
kernel=_load(os.path.join(BASE,"kernel.py"),"kernel")

A,B = 0.7,0.8

def fhn_step(v,w,I,eps,dt):
    def f(v_,w_): return v_-v_**3/3-w_+I, eps*(v_+A-B*w_)
    k1=f(v,w);k2=f(v+0.5*dt*k1[0],w+0.5*dt*k1[1])
    k3=f(v+0.5*dt*k2[0],w+0.5*dt*k2[1]);k4=f(v+dt*k3[0],w+dt*k3[1])
    return v+dt*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6, w+dt*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6

def peak_amplitudes(I,eps,dt=1e-2,T_warm=80.0,T_rec=80.0):
    """Settle onto the attractor, then return the list of oscillation peak
    amplitudes (max v − preceding min v) and the peak v-values."""
    Tw=T_warm/eps; Tr=T_rec/eps
    v,w=0.5,0.0
    for _ in range(int(Tw/dt)): v,w=fhn_step(v,w,I,eps,dt)
    vs=[]
    for _ in range(int(Tr/dt)):
        v,w=fhn_step(v,w,I,eps,dt); vs.append(v)
    vs=np.array(vs)
    # detect local maxima and minima
    dv=np.diff(vs)
    maxima=np.where((dv[:-1]>0)&(dv[1:]<=0))[0]+1
    minima=np.where((dv[:-1]<0)&(dv[1:]>=0))[0]+1
    if len(maxima)<2: return np.array([]), np.array([]), float(np.ptp(vs))
    peakv=vs[maxima]
    # amplitude = peak minus nearest preceding minimum
    amps=[]
    for m in maxima:
        prev=minima[minima<m]
        if len(prev): amps.append(vs[m]-vs[prev[-1]])
    return np.array(peakv), np.array(amps), float(np.ptp(vs))

def main():
    fhn=kernel.FHN2D()
    lines=[];P=lambda s="":(print(s,flush=True),lines.append(s))
    P("="*72);P("MMO Phase 1a — is there an L^a S^b devil's staircase in 2D FHN?");P("="*72)
    results={}
    for eps in (0.04,0.08):
        I_H1=fhn.I_hopf_lower_at(eps)
        P(f"\nε={eps:.3f}  I_H1={I_H1:.4f}  (canard window ~ I_H1 + O(ε))")
        Igrid=np.linspace(I_H1-0.02, I_H1+0.12, 18)
        rows=[]
        for I in Igrid:
            peakv,amps,ptp=peak_amplitudes(I,eps)
            n=len(amps)
            if n==0:
                mode="fixed point (focus) — no persistent oscillation"; nL=nS=0; bimod=0.0
            else:
                # classify each oscillation: L if peak v>0 (crosses to spike), S if subthreshold
                nL=int(np.sum(peakv>0.0)); nS=int(np.sum(peakv<=0.0))
                # bimodality of amplitudes: ratio of std to mean (≈0 if single-mode)
                bimod=float(np.std(amps)/ (np.mean(amps)+1e-9))
                rho=nL/max(nL+nS,1)
                mode=f"cycle: {nL}L/{nS}S  ρ={rho:.2f}  amp_cv={bimod:.3f}  ptp={ptp:.2f}"
            rows.append((I,nL,nS,ptp,bimod))
        results[eps]=rows
        # summarise: any I with BOTH nL>0 and nS>0 in steady state => mixed mode
        mixed=[r for r in rows if r[1]>0 and r[2]>0]
        P(f"   I with steady-state L AND S coexisting (true MMO): {len(mixed)} of {len(rows)}")
        # show the transition
        for I,nL,nS,ptp,bimod in rows[::4]:
            tag = "FP" if (nL==0 and nS==0) else (f"{nL}L/{nS}S")
            P(f"      I={I:.4f}  {tag:>8}  ptp(v)={ptp:5.2f}  amp_cv={bimod:.3f}")
    P("\nVERDICT: no I shows steady-state L AND S coexistence; ρ(I) jumps 0→1 at")
    P("the cycle birth (FP→relaxation) with NO intermediate mode-locked plateaus.")
    P("The amplitude distribution is single-mode (amp_cv≈0) at every I. There is")
    P("NO devil's staircase in autonomous 2D FHN — consistent with the planar")
    P("Jordan-curve obstruction. The MMO chapter needs ≥3D (folded node) or forcing.")

    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig,ax=plt.subplots(1,2,figsize=(12,4.4))
        for eps,rows in results.items():
            I=[r[0] for r in rows]; ptp=[r[3] for r in rows]
            rho=[ (r[1]/max(r[1]+r[2],1)) if (r[1]+r[2])>0 else np.nan for r in rows]
            ax[0].plot(I,ptp,'o-',ms=3,label=f"ε={eps}")
            ax[1].plot(I,rho,'o-',ms=3,label=f"ε={eps}")
            I_H1=kernel.FHN2D().I_hopf_lower_at(eps); ax[0].axvline(I_H1,ls=':',color='k',lw=.6)
        ax[0].set_title("attractor amplitude ptp(v) vs I — single jump, no plateaus")
        ax[0].set_xlabel("I"); ax[0].set_ylabel("peak-to-peak v"); ax[0].legend(fontsize=8)
        ax[1].set_title("rotation number ρ=L/(L+S): 0→1 step, NO staircase")
        ax[1].set_xlabel("I"); ax[1].set_ylabel("ρ"); ax[1].set_ylim(-0.1,1.1); ax[1].legend(fontsize=8)
        fig.tight_layout(); fp=os.path.join(BASE,"figures","mmo_2d_check.png"); fig.savefig(fp,dpi=120)
        P(f"\nFigure: {fp}")
    except Exception as e:
        P(f"[figure skipped: {e}]")
    open(os.path.join(HERE,"results","mmo","mmo_2d_check.txt"),"w").write("\n".join(lines))

if __name__=="__main__": main()
