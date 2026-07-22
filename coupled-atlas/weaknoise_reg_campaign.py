"""Regularized-recursion campaign: (1) field-march h-convergence, (2) delta->0 extrapolation of v0,v1,v2.
Credibility test: v0,v1 must extrapolate to the pinned 0.1339,+0.110; then v2 is read off the same extrapolation."""
import numpy as np, time
from weaknoise_regularized import run

def sweep(deltas, h, N, seeds=(1,2,3)):
    rows=[]
    for delta in deltas:
        accs=[run(delta=delta,N=N,h=h,seed=s) for s in seeds]
        def M(k): return np.mean([a[k] for a in accs])
        def S(k): return np.std([a[k] for a in accs])/np.sqrt(len(seeds))
        rows.append(dict(delta=delta,v0=M('v0'),v0e=S('v0'),m1=M('m1'),v1=M('v1'),v1e=S('v1'),
                         Y3sq=M('Y3sq'),cov=M('covY2Y4'),Y1Y5=M('Y1Y5'),v2=M('v2'),v2e=S('v2')))
        r=rows[-1]
        print(f"  d={delta:.3f} h={h:.0e}: v0={r['v0']:.4f}({r['v0e']:.4f}) m1={r['m1']:+.4f} "
              f"v1={r['v1']:+.4f}({r['v1e']:.4f})  v2={r['v2']:+.4f}({r['v2e']:.4f})  "
              f"[Y3^2={r['Y3sq']:.3f} 2cov={2*r['cov']:+.3f} 2Y1Y5={2*r['Y1Y5']:+.3f}]",flush=True)
    return rows

if __name__=="__main__":
    print("=== (1) field-march h-convergence at delta=0.09, N=12000x3 ===",flush=True)
    for h in [1.5e-3,1e-3,6e-4]:
        t0=time.time(); sweep([0.09],h=h,N=12000); print(f"     ({time.time()-t0:.0f}s)",flush=True)

    print("\n=== (2) delta-sweep at h=8e-4, N=18000x3 ===",flush=True)
    rows=sweep([0.14,0.11,0.08,0.06,0.045],h=8e-4,N=18000)

    d=np.array([r['delta'] for r in rows])
    for key in ['v0','v1','v2']:
        y=np.array([r[key] for r in rows])
        # fit y = a + b delta + c delta^2  (allow O(delta) boundary + O(delta^2))
        A=np.vstack([np.ones_like(d),d,d**2]).T
        c,*_=np.linalg.lstsq(A,y,rcond=None)
        # also linear in delta^2 only (if boundary term absent)
        A2=np.vstack([np.ones_like(d),d**2]).T
        c2,*_=np.linalg.lstsq(A2,y,rcond=None)
        print(f"  {key}: extrap(a+b d+c d^2)={c[0]:+.4f}   extrap(a+c d^2)={c2[0]:+.4f}")
    print("\n  targets: v0=0.1339 v1=+0.110 ; v2 = the deliverable")
