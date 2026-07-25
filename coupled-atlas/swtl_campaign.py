"""
Swallowtail (q=3) production campaign: v3..v6 across grids, with Richardson extrapolation.
Runs each (k,n) as an ISOLATED subprocess (swtl_production.py 3 k n) so caches never cross-contaminate.
Checkpoints after every job to swtl_results.json, so partial results (v3,v4) survive even if v5/v6 run
for hours. Ordered cheap->expensive.  Launch:  python3 swtl_campaign.py
"""
import subprocess, json, re, os, sys, time

HERE=os.path.dirname(os.path.abspath(__file__))
RESULTS=os.path.join(HERE,'swtl_results.json')
# (k, [grids])  -- cheap first; grids chosen for n->inf Richardson extrapolation
SCHEDULE=[
    (3,[24,30,36]),
    (4,[24,30,36]),
    (5,[16,20,24]),
    (6,[10,12,14]),
]

def load():
    if os.path.exists(RESULTS):
        with open(RESULTS) as f: return json.load(f)
    return {}

def save(d):
    with open(RESULTS,'w') as f: json.dump(d,f,indent=2)

def run_job(k,n):
    t0=time.time()
    p=subprocess.run([sys.executable,os.path.join(HERE,'swtl_production.py'),'3',str(k),str(n)],
                     capture_output=True,text=True,cwd=HERE)
    out=p.stdout+p.stderr
    m=re.search(rf'v{k}\(q=3, n={n}\)\s*=\s*([-+0-9.eE]+)',out)
    val=float(m.group(1)) if m else None
    return val, time.time()-t0, out

def richardson(ns, vs):
    # fit v(n) = v_inf + a/n + b/n^2  (or linear if only 2 points)
    import numpy as np
    ns=np.array(ns,float); vs=np.array(vs,float)
    if len(ns)>=3:
        A=np.vstack([np.ones_like(ns),1/ns,1/ns**2]).T
    else:
        A=np.vstack([np.ones_like(ns),1/ns]).T
    c,*_=np.linalg.lstsq(A,vs,rcond=None); return float(c[0])

def main():
    res=load()
    print(f"=== swallowtail v3..v6 campaign (results -> {RESULTS}) ===",flush=True)
    for k,grids in SCHEDULE:
        kk=str(k); res.setdefault(kk,{})
        for n in grids:
            if str(n) in res[kk]:
                print(f"  v{k}(n={n}) cached = {res[kk][str(n)]:+.6f}",flush=True); continue
            print(f"  running v{k}(n={n}) ...",flush=True)
            val,dt,out=run_job(k,n)
            if val is None:
                print(f"    v{k}(n={n}) FAILED after {dt:.0f}s; tail:\n"+"\n".join(out.splitlines()[-6:]),flush=True)
                continue
            res[kk][str(n)]=val; save(res)
            print(f"    v{k}(n={n}) = {val:+.6f}  [{dt:.0f}s]  (checkpointed)",flush=True)
        # extrapolate this k if we have >=2 grids
        # exclude the 'extrap_ninf' key we write below, else int() raises on every resume
        got=sorted((int(n),v) for n,v in res[kk].items() if n!='extrap_ninf')
        if len(got)>=2:
            ns=[g[0] for g in got]; vs=[g[1] for g in got]
            vinf=richardson(ns,vs); res[kk]['extrap_ninf']=vinf; save(res)
            print(f"  ==> v{k}(n->inf) Richardson = {vinf:+.5f}   from {list(zip(ns,vs))}",flush=True)
    print("\n=== CAMPAIGN COMPLETE ===",flush=True)
    for k in ['3','4','5','6']:
        if k in res and 'extrap_ninf' in res[k]:
            print(f"  v{k} = {res[k]['extrap_ninf']:+.5f}",flush=True)

if __name__=="__main__":
    main()
