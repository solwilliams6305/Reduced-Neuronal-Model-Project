import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
d=np.load("ladder_field_data.npz")
qs=[1,2,3,4]; names={1:"fold",2:"cusp",3:"swallowtail",4:"butterfly"}
cols={1:"#e76f51",2:"#264653",3:"#2a9d8f",4:"#8a5a83"}
ameas={1:1.03,2:1.53,3:1.85,4:2.07}
fig,ax=plt.subplots(1,3,figsize=(14,4.4))

# (a) measured vs predicted exponent
a=ax[0]
ap=[3*q/(q+2) for q in qs]; am=[ameas[q] for q in qs]
a.plot(ap,am,"o",ms=9,color="#264653")
for q in qs: a.annotate(f"q={q}\n{names[q]}",(3*q/(q+2),ameas[q]),fontsize=7,ha="left",va="top",xytext=(4,-2),textcoords="offset points")
a.plot([0.9,2.1],[0.9,2.1],"r--",lw=1)
a.set_xlabel(r"predicted $a=3q/(q+2)$"); a.set_ylabel(r"measured jitter exponent"); a.set_title("(a) unified exponent, q=1..4 (~2%)")

# (b) jitter curves
a=ax[1]
for q in qs:
    D=d[f"D_{q}"]; V=d[f"V_{q}"]; a.loglog(D,V,"o-",color=cols[q],ms=5,label=f"q={q} ({names[q]}): {ameas[q]:.2f}")
a.set_xlabel(r"phase-depth $\Theta_q$"); a.set_ylabel("Var(spacing)"); a.set_title(r"(b) $\mathrm{Var(spacing)}\propto\Theta^{-3q/(q+2)}$"); a.legend(fontsize=7)

# (c) the class map: a(q)=3q/(q+2), boundary a=1 at q=1
a=ax[2]
qg=np.linspace(0.3,4.5,200); a.plot(qg,3*qg/(qg+2),"k-",lw=1.8)
a.scatter(qs,[3*q/(q+2) for q in qs],color=[cols[q] for q in qs],s=60,zorder=5)
a.axhline(1,color="gray",ls="--",lw=1)
a.axhspan(1,2.3,color="#2a9d8f",alpha=0.08); a.axhspan(0.4,1,color="#e76f51",alpha=0.08)
a.text(3.3,1.9,"class I\n(lattice, bounded)",fontsize=8,color="#1d6a5f")
a.text(0.5,0.62,"class III",fontsize=8,color="#a33")
a.annotate("q=1 fold:\nclass II (GUE/Airy)\na=1 boundary",(1,1),fontsize=7.5,ha="left",va="center",
           xytext=(20,-24),textcoords="offset points",arrowprops=dict(arrowstyle="->",lw=0.8))
for q in qs: a.annotate(names[q],(q,3*q/(q+2)),fontsize=7,xytext=(3,6),textcoords="offset points")
a.set_xlabel("catastrophe degree q (A$_{q+1}$)"); a.set_ylabel(r"jitter exponent $a=3q/(q+2)$")
a.set_title("(c) fold q=1 = the class II$\\leftrightarrow$I boundary"); a.set_ylim(0.4,2.3); a.set_xlim(0.3,4.5)

fig.suptitle(r"Catastrophe-ladder node processes unified: jitter $\propto\Theta^{-3q/(q+2)}$; fold (q=1) is the critical GUE(II)$\to$lattice(I) edge",fontsize=11)
fig.subplots_adjust(left=0.06,right=0.98,top=0.88,bottom=0.12,wspace=0.27)
fig.savefig("figures/ladder_processes.png",dpi=130); print("saved figures/ladder_processes.png")
