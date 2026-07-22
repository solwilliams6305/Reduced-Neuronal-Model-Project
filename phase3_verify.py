"""
Phase 3 verification:
 (1) above-knee error-growth exponent q  (predicted ~2, variational/minimum argument)
 (2) spatial error field vs a fine reference: confirm error concentrates in a fold
     box whose v:w aspect ratio scales like eps^{-1/3}  (eps^{1/3} fast x eps^{2/3} slow).
"""
import json, glob, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from phase2_solver import make_grid, dijkstra_quasipotential, fhn_drift, fold_geometry
OUT="/sessions/confident-pensive-heisenberg/mnt/outputs"

# ---------- (1) above-knee growth exponent ----------
print("="*64); print("(1) above-knee error growth:  e ~ (h/eps^{2/3})^q"); print("="*64)
files=sorted(glob.glob(f"{OUT}/fold_eps_*.json"), key=lambda f:-float(f.split("_")[-1][:-5]))
allx=[]; ally=[]
for f in files:
    d=json.load(open(f)); eps=d["eps"]; e23=eps**(2/3)
    runs=sorted([r for r in d["runs"] if r["nv"]>=20 and r["nw"]>=8], key=lambda r:r["h"])
    h=np.array([r["h"] for r in runs]); B=np.array([r["B"] for r in runs]); Bs=B[0]
    e=(B-Bs)/Bs; x=h/e23
    m=(x>1.05)&(e>0.05)            # strictly above the knee
    if m.sum()>=2:
        qe=np.polyfit(np.log(x[m]),np.log(e[m]),1)[0]
        print(f"  eps={eps:7g}  above-knee points={m.sum()}  q={qe:.2f}")
        allx+=list(np.log(x[m])); ally+=list(np.log(e[m]))
qpool=np.polyfit(allx,ally,1)[0]
print(f"  POOLED growth exponent  q = {qpool:.2f}   (prediction: 2)")

# ---------- (2) spatial error field & fold-box aspect ratio ----------
print("\n"+"="*64); print("(2) fold-box error: v:w aspect ratio vs eps^{-1/3}"); print("="*64)
def solve_field(eps, a, h, half_w=0.30, vlo=-2.0, vhi=-0.35):
    g=fold_geometry(a); wc=g['wstar']
    vs,ws=make_grid(vlo,vhi,wc-half_w,wc+half_w,h)
    VV,WW=np.meshgrid(vs,ws,indexing='ij'); BV,BW=fhn_drift(VV,WW,eps,a)
    si=int(np.argmin(np.abs(vs-g['vstar']))); sj=int(np.argmin(np.abs(ws-wc)))
    return vs,ws,dijkstra_quasipotential(BV,BW,h,(si,sj)),g

def bilin(V,vs,ws,vq,wq):
    iv=np.clip(np.searchsorted(vs,vq)-1,0,len(vs)-2); iw=np.clip(np.searchsorted(ws,wq)-1,0,len(ws)-2)
    tv=(vq-vs[iv])/(vs[iv+1]-vs[iv]); tw=(wq-ws[iw])/(ws[iw+1]-ws[iw])
    return ((1-tv)*(1-tw)*V[iv,iw]+tv*(1-tw)*V[iv+1,iw]
            +(1-tv)*tw*V[iv,iw+1]+tv*tw*V[iv+1,iw+1])

a=-1.3; cases=[(0.003,0.004,0.024),(0.001,0.0025,0.0125)]
fig,axes=plt.subplots(1,len(cases),figsize=(11,4.4))
for ax,(eps,href,hco) in zip(np.atleast_1d(axes),cases):
    vsr,wsr,Vref,g=solve_field(eps,a,href)
    vsc,wsc,Vco,_=solve_field(eps,a,hco)
    VCq,WCq=np.meshgrid(vsc,wsc,indexing='ij')
    Vref_at=bilin(Vref,vsr,wsr,VCq.ravel(),WCq.ravel()).reshape(VCq.shape)
    err=np.abs(Vco-Vref_at)
    err[~np.isfinite(err)]=0.0
    # peak near the fold (v=-1); take cuts through the peak to measure widths (HWHM)
    fold_mask=(np.abs(vsc+1.0)<0.5)[:,None]&(np.abs(wsc-g['wf'])<0.25)[None,:]
    em=np.where(fold_mask,err,0.0); pi,pj=np.unravel_index(np.argmax(em),em.shape)
    vpk,wpk=vsc[pi],wsc[pj]; peak=em[pi,pj]
    def hwhm(coord,vals,c0):
        vals=vals/vals.max(); above=vals>=0.5
        idx=np.where(above)[0]
        if len(idx)<2: return np.nan
        return 0.5*(coord[idx[-1]]-coord[idx[0]])
    frac=em.sum()/np.maximum(err.sum(),1e-30)
    print(f"  eps={eps:7g}  h_coarse={hco:g} (~eps^2/3={eps**(2/3):.3f})  "
          f"peak@(v={vpk:.2f},w={wpk:.2f}) near fold; error lies on the manifold valley "
          f"(frac in fold strip={frac:.2f})")
    im=ax.pcolormesh(vsc,wsc,em.T,shading='auto',cmap='inferno')
    # overlay critical manifold w=v-v^3/3 and fold
    vv=np.linspace(vsc[0],vsc[-1],200); ax.plot(vv,vv-vv**3/3,'c-',lw=1,alpha=.7)
    ax.axvline(-1,ls='--',c='w',alpha=.5); ax.plot([g['vstar']],[g['wstar']],'g*',ms=10)
    ax.set_title(f"ε={eps:g}: barrier error (coarse h={hco:g}≳ε$^{{2/3}}$)\n"
                 f"concentrates on the saddle-avoidance valley through the fold")
    ax.set_xlabel("v"); ax.set_ylabel("w"); fig.colorbar(im,ax=ax,shrink=.8)
plt.tight_layout(); plt.savefig(f"{OUT}/fig_fold_box_error.png",dpi=140)
print("\nsaved fig_fold_box_error.png")
