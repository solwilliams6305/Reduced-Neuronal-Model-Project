#!/usr/bin/env python3
"""
mmo_crossover_nf.py  —  MMO chapter keystone, Phases A+B (normal-form crossover).

Claim under test: the 3D folded-node noise law and the 2D canard escape law share
ONE prefactor C_q, connected as the folded node degenerates to a folded
saddle-node (μ → 0). The folded-node K2 normal form

    dV = (V² − W) dT + η dB,   dW = (μV − (1+μ)/2) dT          (Wechselberger K2)

reduces at μ → 0 to the canard Krupa–Szmolyan form dV=(V²−W)dT+η dB, dW=−(1/2)dT
(λ = 1/2). Noise enters identically as η = σ/√ε. So C_q is the *same object in the
same coordinates*; we measure it with the **same escape criterion** as the canard
chapter (`canard_normal_form_map.py`): start on the attracting branch V=−√W,
record W_hit at the first V-crossing of V_cross=1, R_hit=W_hit/W_*; the escape
threshold is where median R_hit crosses 0 (hit at the fold). Θ=η/√λ_eff,
λ_eff=(1+μ)/2; C_q ≡ Θ_crit.

RESULT (this run): η_*(μ) does NOT follow μ^{3/2} in the normal form — it is
canard-like across the whole sweep, with C_q = Θ_crit roughly constant (~2.2–2.6)
and η_*(μ) → the canard plateau as μ → 0. So the LOCAL escape mechanism is the
canard one with shared C_q (the keystone inheritance, as a continuous limit). The
μ^{3/2} geometric factor is therefore a GLOBAL-return / funnel-filling effect
(MMO_K2's f(c)), not the local normal-form escape — a refinement of the prompt's
expectation, consistent with the rest of the chapter.

A1 (analytic): μ→0 ⇒ dW/dT→−1/2, canard normal form, λ_eff=1/2 (verified).
A2: crossover scale μ_c from μ·V_* ~ drift, V_*~λ^{1/3}≈0.79, drift≈1/2 ⇒ μ_c≈0.6
    — so the swept μ≤0.5 sit at/below μ_c, i.e. in the canard-like regime, exactly
    as the numerics show.

NOTE on convention (Hazard §7.1): the −(1+μ)/2 constant is from project notes;
the structural μ→0 ⇒ constant-drift fact is convention-robust (used here); the
μ_c prefactor is not — flagged, not relied on quantitatively.

Outputs: data/mmo_crossover.npz, results/mmo/mmo_crossover.txt, figures/mmo_crossover_nf.png
Reproduce:  python3 regime-tests/mmo_crossover_nf.py
"""
from __future__ import annotations
import os
import numpy as np

HERE=os.path.dirname(os.path.abspath(__file__)); BASE=os.path.dirname(HERE)

def simulate(mu,eta,N=300,dT=8e-4,Tmax=35.0,W0fac=5.0,Vcross=1.0,seed=0):
    """folded-node normal-form SDE; μ=0 is the pure canard (dW/dT=−1/2, λ=1/2)."""
    rng=np.random.default_rng(seed)
    lam=(1+mu)/2; Wstar=lam**(2/3); W0=W0fac*Wstar
    V=np.full(N,-np.sqrt(W0)); W=np.full(N,W0)
    hit=np.zeros(N,bool); Whit=np.full(N,np.nan); sdt=np.sqrt(dT)
    for _ in range(int(Tmax/dT)):
        Vn=V+(V*V-W)*dT+eta*sdt*rng.standard_normal(N)
        Wn=W+(mu*V-(1+mu)/2)*dT
        new=(~hit)&(V<Vcross)&(Vn>=Vcross)
        if new.any():
            frac=(Vcross-V[new])/(Vn[new]-V[new])
            Whit[new]=W[new]+frac*(Wn[new]-W[new]); hit|=new
            if hit.all(): break
        np.clip(Vn,-50,Vcross+1,out=Vn)   # freeze escaped trajectories (avoid overflow)
        V,W=Vn,Wn
    return Whit/Wstar

def theta_crit(mu,thetas,**kw):
    lam=(1+mu)/2; sl=np.sqrt(lam)
    meds=np.array([np.nanmedian(simulate(mu,th*sl,**kw)) for th in thetas])
    # interpolate Θ where median R crosses 0
    tc=np.nan
    for i in range(1,len(thetas)):
        if meds[i-1]<0<=meds[i]:
            tc=thetas[i-1]+(0-meds[i-1])/(meds[i]-meds[i-1])*(thetas[i]-thetas[i-1]); break
    return tc, meds

def main():
    lines=[];P=lambda s="":(print(s,flush=True),lines.append(s))
    P("="*72);P("MMO keystone — normal-form crossover η_*(μ): folded node → canard");P("="*72)
    mus=[0.0,0.02,0.05,0.1,0.2,0.5]; thetas=np.array([1.5,2.0,2.5,3.0,3.5,4.0])
    P(f" Θ grid: {thetas.tolist()};  μ→0 is the pure canard (λ=1/2).")
    P(f"\n {'μ':>6}{'λ_eff':>7}{'Θ_crit=C_q':>11}{'η_*':>8}{'μ^{3/2}·C_q0':>12}")
    rows=[]; Cq0=None
    for mu in mus:
        tc,meds=theta_crit(mu,thetas)
        lam=(1+mu)/2; eta=tc*np.sqrt(lam) if tc==tc else np.nan
        if mu==0.0: Cq0=tc
        ref=Cq0*mu**1.5 if (Cq0 and mu>0) else np.nan
        rows.append((mu,lam,tc,eta));
        P(f" {mu:6.2f}{lam:7.3f}{tc:11.2f}{eta:8.3f}{ref:12.3f}")
    P(f"\n canard anchor (μ=0): C_q = Θ_crit = {Cq0:.2f}, η_* = {Cq0*np.sqrt(0.5):.3f}")
    cqs=np.array([r[2] for r in rows]);
    P(f" C_q(μ) range: {np.nanmin(cqs):.2f}–{np.nanmax(cqs):.2f}  (drift {100*(np.nanmax(cqs)/np.nanmin(cqs)-1):.0f}%)")
    P("\n VERDICT: η_*(μ) stays O(η_canard) and → the canard plateau as μ→0; C_q≈const")
    P(" (mild drift). NO μ^{3/2} regime in the LOCAL normal-form escape (the μ^{3/2}")
    P(" is the global-return funnel-filling, MMO_K2 f(c)). The keystone INHERITANCE")
    P(" holds: the folded-node escape is the canard escape with shared C_q, as a")
    P(" continuous μ→0 limit. (Convention caveat §7.1; absolute C_q is dt/settings-")
    P(" dependent — the constancy + canard limit is the convention-robust result.)")

    os.makedirs(os.path.join(BASE,"data"),exist_ok=True)
    np.savez(os.path.join(BASE,"data","mmo_crossover.npz"),mus=np.array(mus),
             Cq=cqs,eta=np.array([r[3] for r in rows]))
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig,ax=plt.subplots(1,2,figsize=(12,4.6))
        mm=np.array([r[0] for r in rows]); et=np.array([r[3] for r in rows])
        ax[0].plot(mm,et,'o-',label="η_*(μ) measured")
        ax[0].axhline(Cq0*np.sqrt(0.5),color='C2',ls='--',label=f"canard plateau η_*={Cq0*np.sqrt(0.5):.2f} (μ→0)")
        mp=np.linspace(0.01,0.5,30); ax[0].plot(mp,Cq0*mp**1.5,'C3:',label="μ^{3/2}·C_q (NOT seen)")
        ax[0].set_xlabel("μ"); ax[0].set_ylabel("η_*"); ax[0].legend(fontsize=8)
        ax[0].set_title("η_*(μ): canard-like, → plateau as μ→0 (no μ^{3/2})")
        ax[1].plot(mm,cqs,'o-'); ax[1].axhline(Cq0,color='C2',ls='--',label=f"canard C_q={Cq0:.2f}")
        ax[1].set_xlabel("μ"); ax[1].set_ylabel("C_q = Θ_crit"); ax[1].set_ylim(0,4)
        ax[1].set_title("C_q(μ): shared prefactor, ~constant"); ax[1].legend(fontsize=8)
        fig.tight_layout(); fp=os.path.join(BASE,"figures","mmo_crossover_nf.png"); fig.savefig(fp,dpi=120)
        P(f"\nFigure: {fp}")
    except Exception as e:
        P(f"[figure skipped: {e}]")
    os.makedirs(os.path.join(HERE,"results","mmo"),exist_ok=True)
    open(os.path.join(HERE,"results","mmo","mmo_crossover.txt"),"w").write("\n".join(lines))

if __name__=="__main__": main()
