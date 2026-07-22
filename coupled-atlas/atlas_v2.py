"""
Consolidated singularity-labeled regime atlas (v2) — T2-consistent.
Two figures: (1) the regime table; (2) the fold/cusp beta-inversion panel.
Organizing principle = CATASTROPHE LADDER (NOT higher-order TW, refuted); cusp = asymmetric PIV-family.
"""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------- (1) the regime table ----------------
# regime, codim, parameter region (coupled FHN), edge marginal (beta-family), beta-trend, process, status
rows = [
 ("generic fold","0","$g$ generic; $\\sigma^2\\!\\sim\\!\\varepsilon$","$\\mathrm{TW}_\\beta$","skew $\\downarrow\\beta$","Airy$_2$ (forced ✓; intrinsic open)","PROVED [A]","#E8EEF7"),
 ("folded node","1","canard window; $\\mu=\\lambda_w/\\lambda_s$","$\\mathrm{TW}_\\beta$","skew $\\downarrow\\beta$","Airy$_2$-type (cited)","PROVED [BlowDown]","#E1F1FB"),
 ("cusp","2","$g\\!\\to\\!g_{\\rm crit}\\!\\propto\\!\\sqrt{\\varepsilon}$","$\\mathcal{W}_\\beta$  (asym. PIV-family)","skew $\\uparrow\\beta$ (inverted)","Weber proc. (forced ✓; intrinsic open)","CONSTRUCTED [B,T2]","#F3EAF6"),
 ("swallowtail","3","+1 tuned param","$\\mathcal{W}^{(3)}_\\beta$","—","— (open)","NUMERICAL [C]","#E7F4EE"),
 ("butterfly","4","+2 tuned params","$\\mathcal{W}^{(4)}_\\beta$","—","— (open)","NUMERICAL [C]","#FBF0E5"),
]
headers=["singularity","codim","parameter region (coupled FHN)","edge marginal ($\\beta$-family)","$\\beta$-trend","limiting process","status"]
colw=[0.11,0.05,0.20,0.21,0.11,0.21,0.13]
x0=np.concatenate([[0.005],0.005+np.cumsum(colw)])

fig,ax=plt.subplots(figsize=(16,4.6)); ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off')
ytop=0.82; rh=0.135; hh=0.11
for j,h in enumerate(headers):
    ax.text(x0[j]+colw[j]/2,ytop+hh/2,h,ha='center',va='center',fontsize=8.6,fontweight='bold')
ax.plot([0.003,0.997],[ytop,ytop],color='k',lw=1.2)
for i,r in enumerate(rows):
    yb=ytop-(i+1)*rh
    ax.add_patch(plt.Rectangle((0.003,yb),0.994,rh,facecolor=r[7],edgecolor='none',zorder=0))
    for j in range(7):
        fw='bold' if j in (0,3) else 'normal'; fs=8.7 if j not in (2,3,5) else 7.7
        ax.text(x0[j]+colw[j]/2,yb+rh/2,r[j],ha='center',va='center',fontsize=fs,fontweight=fw)
# physical/abstract divider after row 3
ax.plot([0.003,0.997],[ytop-3*rh,ytop-3*rh],color='crimson',lw=1.4,ls='--')
ax.text(0.012,ytop-1.5*rh,"PHYSICAL\n(coupled FHN,\n2 params)",rotation=90,ha='center',va='center',fontsize=6.8,color='navy',fontweight='bold')
ax.text(0.012,ytop-4*rh,"ABSTRACT\nladder",rotation=90,ha='center',va='center',fontsize=6.8,color='saddlebrown',fontweight='bold')
ax.text(0.5,0.965,"Singularity-labeled regime atlas — coupled FitzHugh–Nagumo (v2, T2-consistent)",ha='center',fontsize=13.5,fontweight='bold')
ax.text(0.5,0.04,"Organizing principle: CATASTROPHE LADDER (codim $q$). Edge tails $(2q{+}1,\\,3q/2)$; $\\beta=4/\\eta^2$.  "
        "Cusp = ASYMMETRIC PIV-family isomonodromy — NOT standard PIV, NOT higher-order TW (ratio 5/3$\\neq$2, refuted).  "
        "Coupling: $\\Delta(g)=2\\sqrt{-2g/3}$, $g_{\\rm crit}\\!\\propto\\!\\sqrt{\\varepsilon}$.",
        ha='center',fontsize=8.3,style='italic',bbox=dict(boxstyle='round',fc='#FFF7E6',ec='#C9A24B'))
plt.savefig("figures/regime_atlas_v2.png",dpi=115,bbox_inches='tight'); print("saved figures/regime_atlas_v2.png")

# ---------------- (2) the beta-inversion panel ----------------
fold,cusp=np.load("beta_axis_res.npy"); betas=[1,2,4]
TWskew=[0.293,0.224,0.165]
fig2,ax2=plt.subplots(figsize=(6.2,4.4))
ax2.plot(betas,fold,'o-',color='navy',ms=8,label='fold q=1 (escape-loc sim)')
ax2.plot(betas,TWskew,'^--',color='navy',ms=7,alpha=0.5,label='$\\mathrm{TW}_\\beta$ ref (proved class)')
ax2.plot(betas,cusp,'s-',color='crimson',ms=8,label='cusp q=2 ($\\mathcal{W}_\\beta$)')
ax2.set_xlabel(r'$\beta=4/\eta^2$'); ax2.set_ylabel('escape-location skewness'); ax2.set_xticks(betas)
ax2.set_title('$\\beta$-inversion: fold TW$_\\beta$ (skew $\\downarrow$) vs cusp $\\mathcal{W}_\\beta$ (skew $\\uparrow$)')
ax2.legend(fontsize=8.5); ax2.grid(alpha=0.3)
fig2.subplots_adjust(left=0.13,right=0.96,bottom=0.12,top=0.91); plt.savefig("figures/beta_inversion.png",dpi=115); print("saved figures/beta_inversion.png")
