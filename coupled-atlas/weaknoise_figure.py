"""
Program 2, Route 2b — summary figure: Step 0 reconciliation + corrected resurgence picture.
(a) variance coefficient v0: converged perturbative + closed-form + MC all agree at 0.134; FP drifts to 0.20.
(b) beta=2 overshoot: leading/2-term weak-noise partial sums straddle the true W_2 cumulants (non-perturbative).
Values established in PROGRAM2_ROUTE2B_NOTES.md (weaknoise_v1.py, stochastic_scope.py, transseries_orders.py).
"""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

fig,ax=plt.subplots(1,2,figsize=(11.5,4.4))

# (a) Var/eta^2 vs eta^2 : three methods
a=ax[0]
# converged perturbative (var/e2 at small eta, from weaknoise_series clean window + leading v0,v1)
v0,v1=0.1337,0.110
e2=np.linspace(0,0.6,50)
a.plot(e2, v0+v1*e2, "-", color="#2a9d8f", lw=2, label=f"converged: $v_0$+$v_1\\eta^2$ = {v0:.3f}+{v1:.2f}$\\eta^2$")
a.axhline(0.1339,color="#264653",ls="--",lw=1.3,label="closed-form $C_V$=0.1339")
# MC points (stochastic_scope): Var/eta^2 vs eta^2
mc_e2=np.array([0.10,0.15,0.20,0.30,0.50])**2
mc_r=np.array([0.1352,0.1367,0.1387,0.1449,0.1663])
a.plot(mc_e2,mc_r,"o",color="#e9c46a",ms=7,mec="k",label="small-$\\eta$ MC")
# FP points (transseries_orders): var/e2
fp_e2=np.array([0.05,0.08,0.11,0.15,0.20,0.26,0.33,0.42])
fp_var=np.array([0.00990,0.01441,0.01916,0.02587,0.03489,0.04659,0.06135,0.08180])
a.plot(fp_e2,fp_var/fp_e2,"s",color="#e76f51",ms=6,label="FP (contaminated → 0.20)")
a.set_xlabel(r"$\eta^2=4/\beta$"); a.set_ylabel(r"Var$(Y^*)/\eta^2$")
a.set_title("(a) Step 0: leading variance coeff pinned at 0.134\n(FP small-$\\eta$ artifact retracted)")
a.legend(fontsize=7.5,loc="upper left"); a.set_ylim(0.09,0.22); a.set_xlim(0,0.55)

# (b) beta=2 overshoot: partial sums vs truth
a=ax[1]
eta2=2.0
# variance partial sums / eta^2 (compare to true Var/eta^2)
S=[v0, v0+v1*eta2]
true_var_r=0.47/eta2
a.plot([0,1],S,"o-",color="#2a9d8f",ms=9,lw=2,label="Var/$\\eta^2$ partial sums (β=2)")
a.axhline(true_var_r,color="#264653",ls="--",lw=1.5,label=f"true Var/$\\eta^2$={true_var_r:.2f}")
a.annotate("term 2 > term 1\n⇒ β=2 outside\nperturbative regime",xy=(1,S[1]),xytext=(0.35,0.42),
           fontsize=8,arrowprops=dict(arrowstyle="->",lw=1))
a.set_xticks([0,1]); a.set_xticklabels(["$S_0=v_0$","$S_1=v_0+v_1\\eta^2$"])
a.set_ylabel(r"Var$(Y^*)/\eta^2$ at $\beta=2$")
a.set_title("(b) β=2 is non-perturbative (overshoot)\nskew: $s_0\\eta$=1.64 vs true 0.61 (×2.7)")
a.legend(fontsize=8,loc="lower right"); a.set_ylim(0.1,0.42)

fig.suptitle("Route 2b Step 0: converged coefficients ($C_V$=0.134, $s_0$=1.17, $v_1$=+0.110); FP artifacts retracted",fontsize=10.5)
fig.subplots_adjust(left=0.07,right=0.98,top=0.86,bottom=0.13,wspace=0.24)
fig.savefig("figures/weaknoise_sectors.png",dpi=130)
print("saved figures/weaknoise_sectors.png")
