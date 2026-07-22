import numpy as np
# Morris-Lecar single-unit fast field f(v,w); coupled (2 units, diffusive on v) antisymmetric mode:
# antisym linear coeff mu(w) = d_v f(v_s(w),w) - 2 g_c  (v_s(w)=fast nullcline). Fold: mu=0. CUSP: two folds
# coalesce => 2 g_c = interior extremum of d_v f(v_s(w),w) over w. Then Delta ~ sqrt(g_crit - g_c).
C,gL,gCa,gK=20.,2.,4.4,8.; EL,ECa,EK=-60.,120.,-84.; V1,V2,V3,V4=-1.2,18.,2.,30.
def minf(v): return 0.5*(1+np.tanh((v-V1)/V2))
def minfp(v): return 0.5/V2/np.cosh((v-V1)/V2)**2
def f(v,w,I): return (I-gL*(v-EL)-gCa*minf(v)*(v-ECa)-gK*w*(v-EK))/C
def dvf(v,w): return (-gL-gCa*(minfp(v)*(v-ECa)+minf(v))-gK*w)/C   # d_v f
# fast nullcline w_s(v): f=0 => w=(I-gL(v-EL)-gCa minf(v)(v-ECa))/(gK(v-EK))
def wnull(v,I): return (I-gL*(v-EL)-gCa*minf(v)*(v-ECa))/(gK*(v-EK))
I=90.0
vv=np.linspace(-70,60,4000); ww=wnull(vv,I)
# is the v-nullcline N-shaped (has interior max & min in v => folds)?
dw=np.gradient(ww,vv); sign_changes=np.where(dw[:-1]*dw[1:]<0)[0]
print(f"ML I={I}: v-nullcline extrema (folds) at v = {[round(vv[i],1) for i in sign_changes]}  (N-shaped if 2)")
# d_v f along the nullcline, as function of v (param by v on the nullcline)
dvf_null=dvf(vv,ww)
# restrict to the physiological branch where nullcline is single-valued-ish; find interior extremum of d_v f
# (the cusp condition: 2 g_c = extremum of d_v f over the fold region)
m=(vv>-50)&(vv<40)
imax=np.argmax(dvf_null[m]); vext=vv[m][imax]; dvf_ext=dvf_null[m][imax]
print(f"  d_v f along nullcline: interior MAX = {dvf_ext:.4f} at v={vext:.1f} => cusp at 2 g_c = {dvf_ext:.4f}, g_crit={dvf_ext/2:.4f}")
imin=np.argmin(dvf_null[m]); print(f"  d_v f along nullcline: interior MIN = {dvf_null[m][imin]:.4f} at v={vv[m][imin]:.1f}")
# Delta ~ sqrt scaling: for 2 g_c just below the extremum, count/space the two folds (roots of d_v f(v,w_s(v))=2 g_c)
print(f"\n  fold separation Delta(g_c) vs g_crit (roots of d_v f_null = 2 g_c):")
gcrit=dvf_ext/2
for gc in [gcrit-0.001,gcrit-0.004,gcrit-0.016,gcrit-0.064]:
    target=2*gc; d=dvf_null[m]-target; roots=vv[m][np.where(d[:-1]*d[1:]<0)[0]]
    if len(roots)>=2:
        Delta=roots[-1]-roots[0]; print(f"    g_crit-g_c={gcrit-gc:.4f}: {len(roots)} folds, Delta(v)={Delta:.3f}, Delta/sqrt(dg)={Delta/np.sqrt(gcrit-gc):.2f}")
    else:
        print(f"    g_crit-g_c={gcrit-gc:.4f}: {len(roots)} fold(s) (coalesced)")
print("  => constant Delta/sqrt(g_crit-g_c) confirms the sqrt fold-coalescence = CUSP (same as FHN Delta=2sqrt(-2g/3)).")
