"""
SWALLOWTAIL (q=3) production-engine driver for v3..v6, to the cusp paper's standard.

Reuses the validated cusp Wick/transfer engine (chaos_engine + chaos_diagram + chaos_transfer),
monkeypatching the THREE q-dependent points:
  (1) backbone potential + recessive IC  (chaos_engine.V, chaos_engine.build_backbone)
  (2) node V-derivative table Vd_j       (chaos_diagram._get_sub)
  (3) isolated _ybase caches             (chaos_diagram.Ybase, q-tagged filenames)
The node-shift functionals (_yexprs), vk_terms combinatorics, and the moment engine are q-independent.

GATE: q=2 must reproduce the certified cusp anchor v6(n=10)=80.94197 (raw, pre-extrapolation) before
q=3 is trusted.  Usage:  python3 swtl_production.py <q> <k> <n>
"""
import os, sys, time, numpy as np, sympy as sp
os.environ.setdefault('CT_NUMBA', '1')
import chaos_engine as CE
import chaos_diagram as CD
import chaos_transfer as CT
from math import comb, factorial
from _v6_driver import vk_terms, _yfile

Q = int(sys.argv[1]) if len(sys.argv) > 1 else 3
K = int(sys.argv[2]) if len(sys.argv) > 2 else 3
N = int(sys.argv[3]) if len(sys.argv) > 3 else 24

# ---------- (1) backbone: V_q and recessive IC = -sqrt(V) = -|Y0|^{q/2} ----------
def Vq(Y): return np.sign(Y)*np.abs(Y)**Q
def build_backbone_q(Y0=4.0, Yend=-2.6, h=2e-4):
    nY=int(round((Y0-Yend)/h)); Yg=Y0-np.arange(nY+1)*h
    u0=np.empty(nY+1); u0p=np.empty(nY+1); ps=np.empty(nY+1); psp=np.empty(nY+1)
    u0[0],u0p[0]=1.0,-(Y0**(Q/2.0)); ps[0],psp[0]=0.0,1.0; dY=-h
    def rhs(Y,y,yp): return yp, Vq(Y)*y
    for i in range(nY):
        Y=Yg[i]
        for (arr,arrp) in ((u0,u0p),(ps,psp)):
            y,yp=arr[i],arrp[i]
            k1y,k1p=rhs(Y,y,yp); k2y,k2p=rhs(Y+0.5*dY,y+0.5*dY*k1y,yp+0.5*dY*k1p)
            k3y,k3p=rhs(Y+0.5*dY,y+0.5*dY*k2y,yp+0.5*dY*k2p); k4y,k4p=rhs(Y+dY,y+dY*k3y,yp+dY*k3p)
            arr[i+1]=y+dY/6*(k1y+2*k2y+2*k3y+k4y); arrp[i+1]=yp+dY/6*(k1p+2*k2p+2*k3p+k4p)
    return Yg,u0,u0p,ps,psp
CE.V = Vq
CE.build_backbone = build_backbone_q

# ---------- (2) node V-derivative table: V=(-1)^{q+1}Y^q near node (Y<0) ----------
def _get_sub_q(Us,s,Vs,u0ps,Vst):
    key='sub'
    if key in CD._UBCACHE: return CD._UBCACHE[key]
    R=16; Ub={}
    for k in range(R):
        Ub[(k,0)]=sp.Integer(0) if k==0 else Us[k][0]; Ub[(k,1)]=Us[k][1]
    for m in range(2,R):
        for k in range(R):
            p=m-2; tot=0
            for j in range(p+1):
                c=comb(p,j); tot+=c*Vs[j]*Ub[(k,p-j)]
                if k>=1 and j<14: tot+=c*s[j]*Ub[(k-1,p-j)]
            Ub[(k,m)]=sp.expand(tot)
    sub={Us[k][m]:Ub[(k,m)] for k in range(R) for m in range(2,R)}
    sub[Us[0][0]]=0; sub[Us[0][1]]=sp.Float(u0ps)
    Yst=CD.GEO['s'][0]
    Vvals={}
    for j in range(14):
        Vvals[Vs[j]] = Vst if j==0 else (((-1)**(Q+1))*factorial(Q)/factorial(Q-j))*(Yst**(Q-j)) if j<=Q else 0.0
    CD._UBCACHE[key]=(sub,Vvals); return sub,Vvals
CD._get_sub = _get_sub_q

# ---------- (3) q-tagged _ybase cache to avoid cusp collision ----------
import pickle
_orig_parse=CD._parse_term
def Ybase_q(n, Us,s,Vs,Ye, u0ps, Vst):
    cf=f'_ybase_q{Q}_{n}.pkl'
    if os.path.exists(cf):
        with open(cf,'rb') as f: return pickle.load(f)
    sub,Vvals=CD._get_sub(Us,s,Vs,u0ps,Vst)
    if Ye is None:
        # order-aware: CD.load_Y() defaults to _yexprs_11.txt (11 orders), which silently
        # IndexErrors for v6 (needs Y up to idx 13) and v7 (idx 15).  _yfile picks the
        # smallest available table that covers the requested order.
        _,_,_,Ye=CD.load_Y(_yfile(n))
    from collections import defaultdict
    terms=[]
    for t in sp.Add.make_args(Ye[n-1]):
        e=sp.expand(t.subs(sub).subs(Vvals))
        for tt in sp.Add.make_args(e):
            pt=_orig_parse(tt,u0ps)
            if pt: terms.append(pt)
    agg=defaultdict(float)
    for c,at in terms: agg[tuple(sorted(at))]+=c
    terms=[(c,list(k)) for k,c in agg.items() if c!=0.0]
    with open(cf,'wb') as f: pickle.dump(terms,f)
    return terms
CD.Ybase = Ybase_q

# ---------- assemble v_K at grid n=N ----------
if __name__=="__main__":
    t0=time.time()
    CD.setup(n_grid=N, MAXORD=CT.DENSE_MAXORD); CT.clear_all(); CD._UBCACHE.clear()
    u0ps=CD.GEO['u0ps']; Vst=CD.GEO['Vst']; Yst=CD.GEO['s'][0]
    print(f"[q={Q}] node Yst={Yst:.6f}  u0ps={u0ps:.4g}  Vst=V(node)={Vst:.5f}  (Yst^{Q}={Yst**Q:.5f})",flush=True)
    Us=[[sp.symbols(f'U{a}_{m}') for m in range(16)] for a in range(16)]
    s=[sp.symbols(f's{j}') for j in range(14)]; Vs=[sp.symbols(f'Vd{j}') for j in range(14)]
    idxs=sorted(set(p for pr in vk_terms(K) for p in pr[:2]))
    Y={idx: CD.Ybase(idx, Us,s,Vs,None, u0ps, Vst) for idx in idxs}
    print(f"[{time.time()-t0:.0f}s] Y{idxs} built",flush=True)
    mom=CT.moment_hybrid
    def pm(A,B):
        tot=0.0
        for ca,aa in A:
            for cb,ab in B:
                if ca==0.0 or cb==0.0: continue
                tot+=ca*cb*mom(aa+ab)
        return tot
    def mean(A): return sum(c*mom(at) for c,at in A)
    tot=0.0
    for a,b,coef in vk_terms(K):
        t1=time.time()
        val=(pm(Y[a],Y[a])-mean(Y[a])**2) if a==b else (pm(Y[a],Y[b])-mean(Y[a])*mean(Y[b]))
        tot+=coef*val
        print(f"  ({a},{b})x{coef}: {val:+.5f}  [{time.time()-t0:.0f}s]",flush=True)
    print(f"\n[{time.time()-t0:.0f}s] v{K}(q={Q}, n={N}) = {tot:+.6f}",flush=True)
    anchors={(2,6,10):80.94197,(2,7,10):-178.60840}
    if (Q,K,N) in anchors:
        ac=anchors[(Q,K,N)]; print(f"  GATE anchor {ac:+.5f}  |diff|={abs(tot-ac):.2e}  {'PASS' if abs(tot-ac)<1e-3 else 'FAIL'}",flush=True)
