"""
Program 2, Route 2b -- pin v2 (O(eta^6) variance coeff) from the DIRECT FULL-FIELD first-node law, carefully.

Strategy (given v2 is precision-limited by direct MC):
  * sweep a moderate eta window where v2*eta^6 is sizable but v3*eta^8 is not (eta in [0.30,0.75]);
  * Richardson h->0 on the variance (Euler-Maruyama weak order 1);
  * CONSTRAINED fit: hold v0,v1 at their converged values (0.1337/0.1343, 0.110), fit v2,v3;
  * jackknife over seeds for an honest error bar; report sensitivity to the assumed v0.
This does NOT beat the deterministic Green's-function route for precision, but gives an independent number and the
first data point on coefficient GROWTH beyond v1 (bears on the asymptotic/resurgent question).

numpy only. Outputs raw per-seed variances so the fit + jackknife are reproducible.
"""
import numpy as np

def V(Y): return np.sign(Y)*Y*Y

def field_var(N, eta, h, seed, Y0=4.0, Yend=-6.0):
    rng=np.random.default_rng(seed)
    nY=int(round((Y0-Yend)/h)); dY=-h; sh=np.sqrt(h)
    u=np.ones(N); up=np.full(N,-Y0)
    found=np.zeros(N,dtype=bool); Ystar=np.full(N,np.nan)
    Yprev=Y0; uprev=u.copy()
    for i in range(nY):
        Y=Y0-i*h; Vi=V(Y)
        dB=sh*rng.standard_normal(N)
        new_up=up+Vi*u*dY+eta*u*dB
        u_next=u+up*dY
        up=new_up
        Ynext=Y0-(i+1)*h
        cross=(~found)&(Ynext<0.0)&(uprev*u_next<0.0)
        if cross.any():
            fr=uprev[cross]/(uprev[cross]-u_next[cross])
            Ystar[cross]=Yprev+fr*(Ynext-Yprev); found[cross]=True
        uprev=u_next; Yprev=Ynext; u=u_next
        if found.all(): break
    Ys=Ystar[found]
    return Ys.var(), Ys.mean(), Ys.size

if __name__=="__main__":
    etas=np.array([0.30,0.375,0.45,0.525,0.60,0.675,0.75])
    seeds=list(range(6))
    Npseed=250000
    data={}   # (h,eta) -> list of (var,mean) per seed
    for h in [2e-3,1e-3]:
        for eta in etas:
            rows=[field_var(Npseed,eta,h,1000*sd+3) for sd in seeds]
            data[(h,round(eta,3))]=rows
            vs=np.array([r[0] for r in rows]); ms=np.array([r[1] for r in rows])
            print(f"h={h:.0e} eta={eta:.3f}: var={vs.mean():.7f}+-{vs.std()/np.sqrt(len(vs)):.7f}  mean={ms.mean():.6f}",flush=True)

    # Richardson h->0 per eta (weak order 1): var0 = 2 var(h1) - var(h2), h2=2 h1
    def agg(h):
        return np.array([np.mean([r[0] for r in data[(h,round(e,3))]]) for e in etas])
    v_h2=agg(2e-3); v_h1=agg(1e-3); v0field=2*v_h1-v_h2
    x=etas**2

    def constrained_fit(varr, v0f, v1f):
        resid=varr - v0f*x - v1f*x**2
        A=np.vstack([x**3, x**4]).T
        c,*_=np.linalg.lstsq(A,resid,rcond=None)
        return c  # v2, v3
    print("\n=== constrained fits (v0,v1 fixed), Richardson h->0 variance ===")
    for v0f in [0.1337,0.1343]:
        for v1f in [0.108,0.110,0.112]:
            c=constrained_fit(v0field, v0f, v1f)
            print(f"  v0={v0f} v1={v1f}: v2={c[0]:+.4f}  v3={c[1]:+.4f}")

    # jackknife over seeds for v2 error (v0=0.1340,v1=0.110)
    print("\n=== jackknife over 6 seeds (v0=0.1340, v1=0.110) ===")
    v2_jack=[]
    for drop in seeds:
        def aggj(h):
            return np.array([np.mean([r[0] for k,r in enumerate([rr for rr in data[(h,round(e,3))]]) if k!=drop]) for e in etas])
        vj=2*aggj(1e-3)-aggj(2e-3)
        c=constrained_fit(vj,0.1340,0.110); v2_jack.append(c[0])
    v2_jack=np.array(v2_jack)
    v2_mean=v2_jack.mean(); v2_err=np.sqrt((len(seeds)-1)/len(seeds)*np.sum((v2_jack-v2_mean)**2))
    print(f"  v2 = {v2_mean:+.4f} +- {v2_err:.4f} (jackknife)")
    print(f"\n  => Var = eta^2( 0.134 + 0.110 eta^2 + {v2_mean:.3f} eta^4 + ...)")
    print(f"     ratios |v1/v0|={0.110/0.134:.2f}  |v2/v1|={abs(v2_mean)/0.110:.2f}")
