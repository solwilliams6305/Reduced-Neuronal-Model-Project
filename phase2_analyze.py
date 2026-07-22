"""Phase 2 analysis: load fold sweep JSONs, test H1.
   (1) per-eps barrier B(h): converges as h->0; knee near eps^{2/3}.
   (2) resolution threshold h_crit(eps) at fixed rel-error tol; fit log-log slope q.
   (3) converged barrier B*(eps) ~ eps  (saddle-avoidance / along-manifold route).
"""
import json, glob, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

OUT="/sessions/confident-pensive-heisenberg/mnt/outputs"
files=sorted(glob.glob(f"{OUT}/fold_eps_*.json"), key=lambda f:-float(f.split("_")[-1][:-5]))
data={}
for f in files:
    d=json.load(open(f)); eps=d["eps"]
    runs=[r for r in d["runs"] if r["nv"]>=20 and r["nw"]>=8]  # drop only absurdly coarse grids
    runs=sorted(runs,key=lambda r:r["h"])
    data[eps]=dict(geom=d["geom"], h=np.array([r["h"] for r in runs]),
                   B=np.array([r["B"] for r in runs]))

TOL=0.25
epss=sorted(data, reverse=True)
hcrit={}; Bstar={}
print(f"{'eps':>8} {'eps^2/3':>9} {'B*(=minh)':>11} {'B*/eps':>8} {'h_crit(25%)':>11} {'hc/eps^2/3':>10}")
for eps in epss:
    h=data[eps]["h"]; B=data[eps]["B"]
    Bs=B[0]; Bstar[eps]=Bs                      # finest-h reference
    e=(B-Bs)/Bs                                  # relative error vs finest
    # first crossing of TOL scanning from fine to coarse
    hc=np.nan
    for k in range(1,len(h)):
        if e[k]>=TOL:
            # log-linear interp in (h, e)
            hc=np.exp(np.interp(TOL,[e[k-1],e[k]],[np.log(h[k-1]),np.log(h[k])])); break
    hcrit[eps]=hc
    e23=eps**(2/3)
    print(f"{eps:8g} {e23:9.4f} {Bs:11.3e} {Bs/eps:8.4f} {hc:11.4f} {hc/e23:10.3f}")

# ---- fit q on the clean cases (eps<=0.01) ----
fiteps=[e for e in epss if e<=0.01 and np.isfinite(hcrit[e])]
le=np.log([e for e in fiteps]); lh=np.log([hcrit[e] for e in fiteps])
q,c=np.polyfit(le,lh,1)
print(f"\nFIT  log h_crit = q*log eps + c :  q = {q:.3f}  (H1 predicts 2/3 = 0.667)")
print(f"     prefactor h_crit ~ {np.exp(c):.3f} * eps^{q:.3f}")
# fit B* ~ eps^p
lbe=np.log([Bstar[e] for e in epss]); lae=np.log(epss)
p,cb=np.polyfit(lae,lbe,1)
print(f"     B*(eps) ~ {np.exp(cb):.3f} * eps^{p:.3f}   (saddle-avoidance route ~ eps^1)")

# ---------------- figures ----------------
plt.figure(figsize=(7,5))
for eps in epss:
    h=data[eps]["h"]; B=data[eps]["B"]
    plt.loglog(h,B,'o-',label=f"ε={eps:g}")
    plt.axvline(eps**(2/3),ls=':',color=plt.gca().lines[-1].get_color(),alpha=.5)
plt.xlabel("mesh spacing h"); plt.ylabel("computed barrier  B(h,ε)")
plt.title("Fold quasipotential barrier vs mesh\n(dotted = ε$^{2/3}$ inner scale; B plateaus once h≲ε$^{2/3}$)")
plt.legend(fontsize=8); plt.grid(True,which='both',alpha=.3); plt.tight_layout()
plt.savefig(f"{OUT}/fig_barrier_vs_h.png",dpi=140)

plt.figure(figsize=(6.5,5))
efit=np.array(fiteps); hfit=np.array([hcrit[e] for e in fiteps])
eout=np.array([e for e in epss if e>0.01 and np.isfinite(hcrit[e])])
hout=np.array([hcrit[e] for e in eout])
plt.loglog(efit,hfit,'ks',ms=9,label="measured $h_{crit}$ (25% error)")
if len(eout): plt.loglog(eout,hout,'s',mfc='none',mec='gray',ms=9,
                         label="ε=0.03 (coarse-grid, excluded from fit)")
xs=np.array([min(epss),max(epss)])
plt.loglog(xs,np.exp(c)*xs**q,'r-',lw=2,label=f"fit slope q={q:.2f}")
anch=np.exp(np.mean(np.log(hfit))-(2/3)*np.mean(np.log(efit)))   # 2/3 line through fit centroid
plt.loglog(xs,anch*xs**(2/3),'b--',lw=2,label="H1 prediction: slope 2/3")
plt.xlabel("ε (timescale separation)"); plt.ylabel("resolution threshold  $h_{crit}$")
plt.title(f"H1 test: mesh needed to resolve the fold\nmeasured exponent q={q:.2f}  vs  predicted 2/3=0.667")
plt.legend(fontsize=9); plt.grid(True,which='both',alpha=.3); plt.tight_layout()
plt.savefig(f"{OUT}/fig_hcrit_vs_eps.png",dpi=140)

plt.figure(figsize=(6.5,5))
plt.loglog(epss,[Bstar[e] for e in epss],'o-',color='purple',label="converged barrier B*(ε)")
plt.loglog(xs,np.exp(cb)*xs**p,'k--',label=f"fit B*~ε$^{{{p:.2f}}}$")
plt.xlabel("ε"); plt.ylabel("converged barrier B*(ε)")
plt.title("Barrier is the along-manifold (saddle-avoidance) route:\nB* ∝ ε, not the frozen (4/3)δ$^{3/2}$=%.3f"%data[epss[0]]["geom"]["dV_frozen"])
plt.legend(fontsize=9); plt.grid(True,which='both',alpha=.3); plt.tight_layout()
plt.savefig(f"{OUT}/fig_Bstar_vs_eps.png",dpi=140)
print("\nsaved figures: fig_barrier_vs_h.png, fig_hcrit_vs_eps.png, fig_Bstar_vs_eps.png")
