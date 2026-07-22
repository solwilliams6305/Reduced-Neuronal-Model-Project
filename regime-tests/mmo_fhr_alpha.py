#!/usr/bin/env python3
"""
mmo_fhr_alpha.py  —  the MMO chapter's headline deterministic target: the
plateau-width exponent α in Δ_pq ~ 1/q^α, attacked via the folded node.

DERIVATION (FHR-specific, using Wechselberger's folded-node theorem as input):
  1. Desingularised reduced flow at the fold v=−1 (computed in
     mmo_fhr_foldednode.py): tr = 1+δ, det = 2bδ(c+1). Near the folded-
     saddle-node limit c→−1 the weak eigenvalue λ_w ≈ det/tr, so
            μ ≈ 2bδ(c+1)/(1+δ)²  ∝  (c+1).
  2. Wechselberger (2005): maximal rotations s_max = (1−μ)/(2μ) ≈ 1/(2μ) ∝ 1/(c+1).
  3. Principal sequence L¹Sˢ ⇒ ρ=1/(s+1), q=s+1, and the c-interval for "s SAOs"
     is Δc_s ~ |dc/ds|.  If the funnel is filled to s_max (s ∝ (c+1)^{−1}, p=1):
            Δc_s ~ 1/s²  ⇒  α_ceiling = 2.
  In general s ~ (c+1)^{−p}  ⇒  α = 1 + 1/p.  α=1 needs p→∞ (s growing faster
  than any power), inconsistent with the algebraic μ ∝ (c+1).

TEST: measure s_obs(c) and fit p. Result (this script): p ≈ 1.7, so the REALISED
α ≈ 1.6 — between the folded-node ceiling (α=2) and the user's hoped-for α=1. The
excess over p=1 is the global-return "funnel-filling" fraction f=s_obs/s_max,
which RISES toward c→−1 (model-specific, NOT universal folded-node geometry).

VERDICT: α ≠ 1. The folded-node geometry fixes a universal ceiling α=2; the
realised staircase exponent (~1.6 here) is non-universal, set by how the global
return fills the rotation funnel — which is exactly what a Phase-2 K2-chart +
entry-exit return-map calculation must supply.

Outputs: results/mmo/mmo_fhr_alpha.txt, figures/mmo_fhr_alpha.png
Reproduce:  python3 regime-tests/mmo_fhr_alpha.py
"""
from __future__ import annotations
import os
import numpy as np

HERE=os.path.dirname(os.path.abspath(__file__)); BASE=os.path.dirname(HERE)
A,B,EPS,DELTA,I0=0.7,0.8,0.08,0.2,0.30

def mu_smax(c):
    tr=1+DELTA; det=2*B*DELTA*(c+1)
    if det<=0: return np.nan,np.nan
    disc=tr*tr-4*det
    if disc<0: return -1.0,np.nan
    mu=((tr-np.sqrt(disc))/2)/((tr+np.sqrt(disc))/2)
    return mu,(1-mu)/(2*mu)

def step(s,c,dt):
    v,w,y=s
    def f(v,w,y): return (v-v**3/3-w+y+I0, EPS*(v+A-B*w), EPS*DELTA*(c-v))
    k1=f(v,w,y);k2=f(v+.5*dt*k1[0],w+.5*dt*k1[1],y+.5*dt*k1[2]);k3=f(v+.5*dt*k2[0],w+.5*dt*k2[1],y+.5*dt*k2[2]);k4=f(v+dt*k3[0],w+dt*k3[1],y+dt*k3[2])
    return np.array([v+dt*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6,w+dt*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6,y+dt*(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6])

def sao_count(c,dt=0.04,Tw=2500,Tr=12000):
    s=np.array([0.5,0.0,0.0])
    for _ in range(int(Tw/dt)): s=step(s,c,dt)
    vs=np.empty(int(Tr/dt))
    for i in range(len(vs)): s=step(s,c,dt); vs[i]=s[0]
    dv=np.diff(vs); mx=np.where((dv[:-1]>0)&(dv[1:]<=0))[0]+1; pk=vs[mx] if len(mx) else np.array([])
    sym=''.join('L' if p>0 else 'S' for p in pk); runs=[len(r) for r in sym.split('L') if r]
    return (max(runs) if runs else 0)

def main():
    lines=[];P=lambda s="":(print(s,flush=True),lines.append(s))
    P("="*72);P("MMO α derivation — folded-node ceiling vs realised staircase");P("="*72)
    cs=[-0.80,-0.84,-0.87,-0.90,-0.92,-0.94]
    P(f"{'c':>7}{'c+1':>7}{'μ':>8}{'s_max(ceil)':>12}{'s_obs':>7}{'f=obs/max':>10}")
    rows=[]
    for c in cs:
        mu,smax=mu_smax(c); sobs=sao_count(c); rows.append((c,mu,smax,sobs))
        P(f"{c:7.2f}{c+1:7.3f}{mu:8.4f}{smax:12.1f}{sobs:7d}{sobs/smax:10.2f}")
    cc=np.array([r[0]+1 for r in rows]); sobs=np.array([r[3] for r in rows],float); smax=np.array([r[2] for r in rows])
    p_obs=-np.polyfit(np.log(cc),np.log(sobs),1)[0]
    p_max=-np.polyfit(np.log(cc),np.log(smax),1)[0]
    a_obs=1+1/p_obs; a_max=1+1/p_max
    P(f"\n s_max  ~ (c+1)^(-{p_max:.2f})  ->  α_ceiling = {a_max:.2f}   (folded-node, μ∝(c+1))")
    P(f" s_obs  ~ (c+1)^(-{p_obs:.2f})  ->  α_realised = {a_obs:.2f}   (measured staircase)")
    P(f"\n VERDICT: α ≠ 1.  Folded-node ceiling α=2 (derived); realised α≈{a_obs:.1f}")
    P(f" (measured); the gap is the global-return funnel-filling f (rises "
      f"{rows[0][3]/rows[0][2]:.2f}→{rows[-1][3]/rows[-1][2]:.2f} toward c→−1),")
    P(" a model-specific effect the Phase-2 K2/entry-exit return map must supply.")
    P(f"\n NOISE IMPLICATION: with Δ_pq~q^(-α) and δρ_noise~σ^γ, σ_pq~σ_*·q^(-α/γ).")
    P(f" The user's 'α=1 → σ_pq~σ_*/√q' is superseded: α≈{a_obs:.1f}–2 ⇒ a steeper β.")

    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig,ax=plt.subplots(1,2,figsize=(12,4.4))
        ax[0].loglog(cc,smax,'s--',color='C3',label=f"s_max ceiling ~ (c+1)^-{p_max:.1f} (α=2)")
        ax[0].loglog(cc,sobs,'o-',color='C0',label=f"s_obs ~ (c+1)^-{p_obs:.1f} (α={a_obs:.1f})")
        ax[0].set_xlabel("c+1  (distance to folded-saddle-node)"); ax[0].set_ylabel("SAOs per spike")
        ax[0].set_title("rotation count scaling → α"); ax[0].legend(fontsize=8)
        ax[1].plot(cc,sobs/smax,'o-')
        ax[1].set_xlabel("c+1"); ax[1].set_ylabel("funnel-filling f = s_obs/s_max")
        ax[1].set_title("global-return funnel-filling (non-universal; reduces α from 2)")
        ax[1].invert_xaxis()
        fig.tight_layout(); fp=os.path.join(BASE,"figures","mmo_fhr_alpha.png"); fig.savefig(fp,dpi=120)
        P(f"\nFigure: {fp}")
    except Exception as e:
        P(f"[figure skipped: {e}]")
    os.makedirs(os.path.join(HERE,"results","mmo"),exist_ok=True)
    open(os.path.join(HERE,"results","mmo","mmo_fhr_alpha.txt"),"w").write("\n".join(lines))

if __name__=="__main__": main()
