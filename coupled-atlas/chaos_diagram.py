"""
Diagram/pairing Wiener-chaos moment engine (replaces the truncated-product assembly).

Y_n is a polynomial in BASE atoms {U_{k,0},U_{k,1},s_j}: U_{k,.} is a clean k-chaos Markovian-chain kernel,
s_j = xi^(j)(node) a boundary leg.  A moment <prod atoms> = sum over CONTRACTION MULTIGRAPHS M (edges pair legs
across atoms; each atom i has degree = its chaos order; boundary atoms degree 1; NO intra-atom legs, NO
boundary-boundary edges [renormalization drop]):
      <prod> = sum_M [ prod_i order_i! / prod_{i<j} M_ij! ] * TN(M)
TN(M): apply beval (1/2 d^j/dY^j at node) for each boundary edge (removes a bulk leg), then einsum-contract the
remaining bulk-bulk edges.  Never builds a tensor above the largest single atom (<=5 legs). Exact (no truncation).
Validate v0,v1,v2 (0.134,0.110,0.11) then run v3,v4.
"""
import numpy as np, itertools, sympy as sp
from math import comb, factorial
from chaos_engine import build
GEO={}
def setup(**kw):
    GEO.clear(); GEO.update(build(**kw)); GEO['h']=GEO['s'][1]-GEO['s'][0]; _MCACHE.clear()

def beval(T, leg, j):
    n=GEO['n']; sw=GEO['sw']; h=GEO['h']
    T=np.moveaxis(T,leg,0)/sw.reshape((n,)+(1,)*(T.ndim-1))
    for _ in range(j): T=np.gradient(T,h,axis=0)
    return GEO.get('bfac',0.5)*T[0]   # boundary-Wick endpoint factor (0.5=symmetric/Strat; tunable)

# ---- atoms: ('U',k,m) -> kernel tensor (order k); ('s',j) -> boundary order j ----
def atom_kernel(a):
    kind=a[0]
    if kind=='U':
        _,k,m=a
        base=GEO['Uk0'][k] if m==0 else GEO['Uk1'][k]
        return base.t[k]   # absorbed symmetric k-tensor
    return None
def atom_order(a):
    return a[1] if a[0]=='U' else 1

def _enum_multigraphs(orders, is_bnd):
    """yield symmetric edge-count dicts {(i,j):m} with deg(i)=orders[i], no i-i, no bnd-bnd edge, bnd deg 1."""
    N=len(orders); pairs=[(i,j) for i in range(N) for j in range(i+1,N) if not(is_bnd[i] and is_bnd[j])]
    res=[]
    def rec(pi, rem):
        if pi==len(pairs):
            if all(r==0 for r in rem): res.append({})
            return
        i,j=pairs[pi]; maxm=min(rem[i],rem[j])
        for m in range(maxm+1):
            rem[i]-=m; rem[j]-=m
            for sub in _rec_collect(pi+1,rem,pairs):
                d=dict(sub);
                if m>0: d[(i,j)]=m
                res.append(d)
            rem[i]+=m; rem[j]+=m
        return
    # simpler explicit recursion returning lists
    N=len(orders)
    def rec2(pi,rem):
        if pi==len(pairs):
            return [dict()] if all(r==0 for r in rem) else []
        # prune: for each atom a, max degree it can still get = sum of rem[partner] over remaining pairs; must >= rem[a]
        avail=[0]*N
        for qi in range(pi,len(pairs)):
            u,v=pairs[qi]; avail[u]+=rem[v]; avail[v]+=rem[u]
        for a in range(N):
            if rem[a]>avail[a]: return []
        i,j=pairs[pi]; out=[]; maxm=min(rem[i],rem[j])
        for m in range(maxm+1):
            rem2=list(rem); rem2[i]-=m; rem2[j]-=m
            for sub in rec2(pi+1,rem2):
                if m>0: sub=dict(sub); sub[(i,j)]=m
                out.append(sub)
        return out
    return rec2(0,list(orders))

_MCACHE={}
def moment(atoms):
    """<prod atoms> full contraction (memoized by canonical atom multiset)."""
    if not atoms: return 1.0
    key=tuple(sorted(atoms))
    if key in _MCACHE: return _MCACHE[key]
    orders=[atom_order(a) for a in atoms]; is_bnd=[a[0]=='s' for a in atoms]
    if sum(orders)%2!=0: _MCACHE[key]=0.0; return 0.0
    # connectivity prune: a bulk atom (order k) needs k pair-slots to OTHER atoms; if isolated-impossible ->0
    tot=0.0
    for M in _enum_multigraphs(orders, is_bnd):
        mult=1.0
        for i in range(len(atoms)): mult*=factorial(orders[i])
        for e,m in M.items(): mult/=factorial(m)
        tot+=mult*_tn(atoms,M)
    _MCACHE[key]=tot
    return tot
def clear_cache(): _MCACHE.clear()

def _tn(atoms, M):
    N=len(atoms)
    # tensors for bulk atoms (copy); boundary atoms handled via edges
    tens={i:(atom_kernel(atoms[i]).copy() if atoms[i][0]=='U' else None) for i in range(N)}
    # apply boundary edges: for each edge (i,j) where one is boundary
    bulk_edges=[]
    for (i,j),m in M.items():
        bi=atoms[i][0]=='s'; bj=atoms[j][0]=='s'
        if bi and bj: return 0.0
        if bi or bj:
            b,bulk=(i,j) if bi else (j,i)
            jorder=atoms[b][1]
            # m should be 1 (boundary degree 1)
            for _ in range(m):
                tens[bulk]=beval(tens[bulk],0,jorder)   # remove a leg
        else:
            bulk_edges.append(((i,j),m))
    # now einsum the bulk-bulk edges among remaining bulk tensors
    # assign letters
    letters=iter('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')
    leg_labels={i:[] for i in range(N) if atoms[i][0]=='U'}
    for (i,j),m in bulk_edges:
        for _ in range(m):
            L=next(letters); leg_labels[i].append(L); leg_labels[j].append(L)
    # build einsum operands
    ops=[]; subs=[]
    for i in leg_labels:
        T=tens[i]
        if T is None: T=np.array(1.0)
        if T.ndim!=len(leg_labels[i]):
            # mismatch => this multigraph is inconsistent for this atom (shouldn't happen)
            return 0.0
        ops.append(T); subs.append(''.join(leg_labels[i]))
    if not ops: return 1.0
    if len(ops)==1: return float(ops[0].sum() if subs[0]=='' else np.einsum(subs[0]+'->',ops[0]))
    return float(np.einsum(','.join(subs)+'->', *ops, optimize='greedy'))

# ---- build Y_n as polynomial in base atoms ----
def load_Y(fname='_yexprs_11.txt'):
    Us=[[sp.symbols(f'U{k}_{m}') for m in range(16)] for k in range(16)]
    s=[sp.symbols(f's{j}') for j in range(14)]; Vs=[sp.symbols(f'Vd{j}') for j in range(14)]
    with open(fname) as f: Ye=[sp.sympify(l.strip()) for l in f if l.strip()]
    return Us,s,Vs,Ye

import pickle, os
_UBCACHE={}
def _get_sub(Us,s,Vs,u0ps,Vst):
    key='sub'
    if key in _UBCACHE: return _UBCACHE[key]
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
    Vvals={Vs[0]:Vst, Vs[1]:-2*GEO['s'][0], Vs[2]:-2.0}
    for j in range(3,14): Vvals[Vs[j]]=0.0
    _UBCACHE[key]=(sub,Vvals); return sub,Vvals

def _parse_term(term,u0ps):
    coef,mono=term.as_coeff_Mul(); atoms=[]; coef=float(coef)
    for fac in sp.Mul.make_args(mono):
        b,ex=fac.as_base_exp(); nm=str(b); ex=int(ex)
        if nm.startswith('U'):
            k=int(nm[1:nm.index('_')]); m=int(nm[nm.index('_')+1:])
            if k==0:
                if m==0: coef=0.0
                else: coef*=u0ps**ex
            else: atoms+=[('U',k,m)]*ex
        elif nm.startswith('s'): atoms+=[('s',int(nm[1:]))]*ex
        else: return None
    return (coef,atoms) if coef!=0.0 else None

def Ybase(n, Us,s,Vs,Ye, u0ps, Vst):
    """return list of (coef_float,[atoms]) for Y_n in base atoms. Cached to disk (n-independent of grid)."""
    cf=f'_ybase_{n}.pkl'
    if os.path.exists(cf):
        with open(cf,'rb') as f: return pickle.load(f)
    sub,Vvals=_get_sub(Us,s,Vs,u0ps,Vst)
    terms=[]
    # term-by-term substitution+expand to avoid a giant intermediate
    for t in sp.Add.make_args(Ye[n-1]):
        e=sp.expand(t.subs(sub).subs(Vvals))
        for tt in sp.Add.make_args(e):
            pt=_parse_term(tt,u0ps)
            if pt: terms.append(pt)
    # merge duplicate atom-multisets
    from collections import defaultdict
    agg=defaultdict(float)
    for c,at in terms: agg[tuple(sorted(at))]+=c
    terms=[(c,list(k)) for k,c in agg.items() if c!=0.0]   # keep tiny-but-real coeffs (~1/u0ps^k, u0ps~2e4)
    with open(cf,'wb') as f: pickle.dump(terms,f)
    return terms

def prod_moment(A,B):
    """<Y_a Y_b> = sum over monomial pairs coef*coef*moment(atoms_a+atoms_b)."""
    tot=0.0
    for ca,aa in A:
        for cb,ab in B:
            tot+=ca*cb*moment(aa+ab)
    return tot
def mean(A): return sum(c*moment(at) for c,at in A)

if __name__=="__main__":
    import time,sys
    what=sys.argv[1] if len(sys.argv)>1 else 'v2'
    def Cov(A,Bx): return prod_moment(A,Bx)-mean(A)*mean(Bx)
    if what=='v2':
        for n in [30,40,50]:
            setup(n_grid=n, MAXORD=3); u0ps=GEO['u0ps']; Vst=GEO['Vst']
            Us,s,Vs,Ye=load_Y(); Y={k:Ybase(k,Us,s,Vs,Ye,u0ps,Vst) for k in range(1,6)}
            v0=prod_moment(Y[1],Y[1]); m1=mean(Y[2]); v1=(prod_moment(Y[2],Y[2])-m1**2)+2*prod_moment(Y[1],Y[3])
            v2=prod_moment(Y[3],Y[3])+2*Cov(Y[2],Y[4])+2*Cov(Y[1],Y[5])
            print(f" n={n}: v0={v0:.4f} v1={v1:+.4f} v2={v2:+.4f}",flush=True)
    if what=='v3':
        for n in [30,40]:
            t0=time.time(); setup(n_grid=n, MAXORD=4); u0ps=GEO['u0ps']; Vst=GEO['Vst']
            Us,s,Vs,Ye=load_Y(); Y={k:Ybase(k,Us,s,Vs,Ye,u0ps,Vst) for k in range(1,8)}
            m4=mean(Y[4]); VarY4=prod_moment(Y[4],Y[4])-m4**2
            v3=VarY4+2*Cov(Y[3],Y[5])+2*Cov(Y[2],Y[6])+2*Cov(Y[1],Y[7])
            print(f" n={n}: Var(Y4)={VarY4:.4f} v3={v3:+.4f}  (target ~-0.3, sign flip)  [{time.time()-t0:.0f}s]",flush=True)
    if what=='v4':
        for n in [24,30]:
            t0=time.time(); setup(n_grid=n, MAXORD=5); u0ps=GEO['u0ps']; Vst=GEO['Vst']
            Us,s,Vs,Ye=load_Y(); Y={k:Ybase(k,Us,s,Vs,Ye,u0ps,Vst) for k in range(1,10)}
            m5=mean(Y[5]); VarY5=prod_moment(Y[5],Y[5])-m5**2
            v4=VarY5+2*Cov(Y[4],Y[6])+2*Cov(Y[3],Y[7])+2*Cov(Y[2],Y[8])+2*Cov(Y[1],Y[9])
            print(f" n={n}: Var(Y5)={VarY5:.4f} v4={v4:+.4f}  (factorial~-4 vs convergent~-0.06)  [{time.time()-t0:.0f}s]",flush=True)
