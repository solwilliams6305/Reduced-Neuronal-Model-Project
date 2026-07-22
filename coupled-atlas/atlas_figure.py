"""
D — the singularity-labeled regime atlas (the 'periodic table' of the catastrophe ladder),
organized by the higher-order-TW / multicritical (Painleve-II-hierarchy) reframe.
Pure rendering of the synthesized data. numpy/matplotlib only.
"""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

# columns: q, catastrophe, V_q(Y), inner stochastic operator, edge law, skew, exk, (L,R) tail exp, process, status
rows = [
 ("1","fold","sign(Y)·Y","stochastic Airy\n$-\\partial^2+Y+\\frac{2}{\\sqrt{\\beta}}\\dot W$",
   "$\\mathrm{TW}_\\beta$","+0.20","+0.02","3 , 1.5","Airy$_2$","CITED (RRV)"),
 ("2","cusp","sign(Y)·Y$^2$","stochastic Weber\n$-\\partial^2+\\mathrm{sgn}(Y)Y^2+\\frac{2}{\\sqrt{\\beta}}\\dot W$",
   "$\\mathcal{W}_\\beta$ (NEW)","+0.61","−0.24","5 , 3","Weber proc.","CONSTRUCTED"),
 ("3","swallowtail","sign(Y)·|Y|$^3$","stoch. higher-Airy (k=3)",
   "$\\mathcal{W}^{(3)}_\\beta$","+0.96","+0.17","7 , 4.5","—","NEW (this work)"),
 ("4","butterfly","sign(Y)·Y$^4$","stoch. higher-Airy (k=4)",
   "$\\mathcal{W}^{(4)}_\\beta$","+1.22","+0.71","9 , 6","—","NEW (this work)"),
]
headers = ["q=k","catastrophe","normal form  $V_q$","stochastic inner operator",
           "marginal edge law","skew","exk","tail exp\n(L , R)","limiting\nprocess","status"]
colw = [0.05,0.12,0.12,0.22,0.13,0.06,0.06,0.08,0.09,0.12]
x0 = np.concatenate([[0.01],0.01+np.cumsum(colw)])
rowcol = ["#E8EEF7","#F3EAF6","#E7F4EE","#FBF0E5"]

fig,ax=plt.subplots(figsize=(15.5,5.6)); ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off')
ytop=0.86; rh=0.165; hh=0.10
# header
for j,h in enumerate(headers):
    ax.text(x0[j]+colw[j]/2, ytop+hh/2, h, ha='center', va='center', fontsize=8.4, fontweight='bold')
ax.plot([0.005,0.995],[ytop,ytop],color='k',lw=1.2)
# rows
for i,r in enumerate(rows):
    yb=ytop-(i+1)*rh
    ax.add_patch(plt.Rectangle((0.005,yb),0.99,rh,facecolor=rowcol[i],edgecolor='none',zorder=0))
    for j,val in enumerate(r):
        fw='bold' if j in (0,4) else 'normal'
        fs=8.8 if j!=3 else 7.6
        ax.text(x0[j]+colw[j]/2, yb+rh/2, val, ha='center', va='center', fontsize=fs, fontweight=fw)
ax.plot([0.005,0.995],[ytop-4*rh,ytop-4*rh],color='k',lw=1.0)
# unfolding arrow (down the ladder = relevant perturbation)
arr=FancyArrowPatch((0.045,ytop-3.55*rh),(0.045,ytop-0.55*rh),arrowstyle='-|>',mutation_scale=16,
                    color='crimson',lw=2.0,connectionstyle="arc3,rad=-0.25")
ax.add_patch(arr)
ax.text(0.005,ytop-2*rh,"unfolding\n(relevant\nperturbation)\nflows UP",rotation=90,ha='center',va='center',
        fontsize=7.2,color='crimson',fontweight='bold')
# title + reframe banner
ax.text(0.5,0.975,"Regime atlas — the catastrophe ladder of noise-induced escape",
        ha='center',fontsize=14,fontweight='bold')
ax.text(0.5,0.045,"Reframe: q = multicritical index k of the higher-order Tracy–Widom (Painlevé-II hierarchy) family.  "
        "Fold = k=1 (KPZ/TW).  Tail exponents 2q+1 (left), 3q/2 (right).  "
        "Generic perturbations flow the class DOWN to the fold; rung q needs codim q−1 tuning.",
        ha='center',fontsize=8.6,style='italic',
        bbox=dict(boxstyle='round',fc='#FFF7E6',ec='#C9A24B'))
plt.savefig("figures/regime_atlas.png",dpi=115,bbox_inches='tight'); print("saved figures/regime_atlas.png")
