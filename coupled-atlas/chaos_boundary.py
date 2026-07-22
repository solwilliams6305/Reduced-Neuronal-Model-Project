"""
Boundary-enabled Wiener-chaos engine (extends chaos_engine.py) -> renormalized v2, v3, v4.

A chaos element is a dict keyed by (bulk_p, bnd_tuple) -> absorbed symmetric tensor of shape (n,)*bulk_p,
where bnd_tuple is a sorted tuple of boundary-noise orders j (each s_j = xi^(j)(node), a 1-chaos "boundary leg").
RENORMALIZATION PRESCRIPTION (drops the divergent boundary self-contractions, keeps finite boundary evals):
  * boundary--bulk contraction  <s_j , a>  ->  BEVAL_j[a] = (1/2) a^(j)(node)   [KEPT]
  * boundary--boundary contraction <s_i,s_j> -> DROPPED                          [renormalization]
Node fields carry a u0(node)=0 factor in the node-adjacent slot, so BEVAL is finite (no leftover divergence).
Validate on v2 (~0.11) before trusting v3,v4.
"""
import numpy as np, itertools
from math import comb, factorial
import chaos_engine as CE
from chaos_engine import build, _sym, _wten

class B:
    def __init__(s,t=None): s.t={} if t is None else dict(t)
    @staticmethod
    def scalar(c): return B({(0,()):np.array(float(c))})
    def add(a,b):
        o=B(a.t)
        for k,v in b.t.items(): o.t[k]=o.t.get(k,0)+v
        return o
    def scale(s,c): return B({k:c*v for k,v in s.t.items()})

# global geometry (set by setup)
GEO={}
def setup(n_grid=40, MAXORD=5, **kw):
    Bd=build(n_grid=n_grid,MAXORD=MAXORD,**kw)
    GEO.clear(); GEO.update(Bd)
    n=Bd['n']; sw=Bd['sw']; w=Bd['w']
    GEO['inv_sw']=1.0/sw
    return Bd

def beval(tensor, leg, j):
    """contract bulk `leg` of absorbed `tensor` with boundary order j: (1/2) d^j/dY^j (kernel)|node, on that leg.
       kernel = tensor/sqrt(w) on that leg; node = grid index 0; derivative via finite differences along the grid."""
    n=GEO['n']; sw=GEO['sw']; s=GEO['s']; h=s[1]-s[0]
    # move leg to axis 0
    T=np.moveaxis(tensor,leg,0)
    # un-absorb that leg: divide by sw along axis0
    T=T/ sw.reshape((n,)+(1,)*(T.ndim-1))
    # j-th derivative along axis0 (one-sided at node = index 0), then value at node
    D=T
    for _ in range(j):
        D=np.gradient(D,h,axis=0)
    val=0.5*D[0]   # (1/2) a^(j)(node); remaining legs still absorbed
    return val   # shape (n,)*(p-1)

def bexpect(A):
    """<A> : sum of fully-contracted terms. A term (p,bnd)->T contributes only if it can fully contract:
       all bulk legs pair among themselves? No -- <I_p>=0 for p>=1. Only (0,()) survives. So <A>=A.t.get((0,()),0)."""
    return float(A.t.get((0,()),np.array(0.0)))

def bmul(A,Bx,MAXTOT):
    """ordinary product of chaos elements with boundary legs, renormalized (drop bnd-bnd contractions).
       Contract r bulk-bulk pairs and any number of bulk-boundary pairs. Enumerate contractions."""
    out=B()
    for (pa,ba),f in A.t.items():
        for (pb,bb),g in Bx.t.items():
            # choose rb = # bulk-bulk contractions (0..min(pa,pb))
            for rb in range(min(pa,pb)+1):
                # choose how many boundary legs of A contract with bulk legs of B (mA) and vice versa (mB)
                for mA in range(len(ba)+1):
                    for mB in range(len(bb)+1):
                        # need enough bulk legs left after bulk-bulk contraction to absorb boundary contractions
                        if mB> pa-rb or mA> pb-rb: continue
                        # resulting bulk order:
                        p_out=(pa-rb-mB)+(pb-rb-mA)
                        bnd_out=len(ba)-mA+len(bb)-mB
                        if p_out+bnd_out>MAXTOT: continue
                        # combinatorial multiplicity & the tensor contraction:
                        term=_contract(pa,ba,f,pb,bb,g,rb,mA,mB)
                        if term is None: continue
                        coef,T,bnds=term
                        key=(p_out,tuple(sorted(bnds)))
                        out.t[key]=out.t.get(key,0)+coef*T
    return out

def _contract(pa,ba,f,pb,bb,g,rb,mA,mB):
    # bulk-bulk: contract first rb bulk legs of f with first rb of g (symmetric kernels)
    n=GEO['n']
    # multiplicity for choosing which legs: C(pa,rb)C(pb,rb)rb!  (bulk-bulk, symmetric)
    coef=comb(pa,rb)*comb(pb,rb)*factorial(rb)
    # boundary of A (mA of them) contract with bulk legs of B; choose which boundary orders and which B bulk legs
    # for symmetric kernels the choice of B-leg is symmetric -> multiplicity C(pb-rb, mA)*mA! * (choose mA of ba)
    # similarly mB of bb contract with A bulk legs
    # We'll handle small numbers by explicit loops over which boundary orders.
    baL=list(ba); bbL=list(bb)
    from itertools import combinations
    total=None; acc=0.0
    # iterate over subsets of A-boundary of size mA and B-boundary of size mB
    for aset in combinations(range(len(baL)),mA):
        for bset in combinations(range(len(bbL)),mB):
            jA=[baL[i] for i in aset]; jB=[bbL[i] for i in bset]
            remA=[baL[i] for i in range(len(baL)) if i not in aset]
            remB=[bbL[i] for i in range(len(bbL)) if i not in bset]
            # contract: first rb bulk of f with first rb of g; then mA boundary(jA) with next mA bulk of g;
            # then mB boundary(jB) with next mB bulk of f. Do bulk-bulk via tensordot, boundary via beval.
            F=f.copy(); G=g.copy()
            # bulk-bulk
            if rb>0:
                K=np.tensordot(F,G,axes=(list(range(rb)),list(range(rb))))
            else:
                K=np.tensordot(F,G,axes=0) if (F.ndim>0 and G.ndim>0) else (F*G if F.ndim==0 or G.ndim==0 else F*G)
            # after tensordot: axes = [F-remaining (pa-rb)] + [G-remaining (pb-rb)]
            nFa=pa-rb; nGa=pb-rb
            # boundary jA contract with G-remaining legs (which are axes nFa..nFa+nGa-1)
            axoff=nFa
            for jj in jA:
                K=beval(K, axoff, jj)  # removes that axis
                # after removal, axes shift; G-remaining count decreases; keep axoff pointing at first G axis
                nGa-=1
            # boundary jB contract with F-remaining legs (axes 0..nFa-1)
            for jj in jB:
                K=beval(K, 0, jj); nFa-=1
            mult= factorial(mA)*comb(pb-rb,mA) * factorial(mB)*comb(pa-rb-mB if False else pa-rb, mB)
            # accumulate with remaining boundary legs
            if total is None:
                total=(np.zeros_like(K), remA+remB)
            acc_term=coef*mult*K
            total=(total[0]+acc_term, remA+remB)
    if total is None: return None
    T,bnds=total
    T=_sym(T) if np.ndim(T)>=2 else T
    return 1.0, T, bnds

def bpow(A,e,MT):
    r=B.scalar(1.0)
    for _ in range(e): r=bmul(r,A,MT)
    return r

if __name__=="__main__":
    # validate v0,v1 through the boundary engine (no boundary noise used yet) to check bmul matches
    for n in [40,60]:
        setup(n_grid=n,MAXORD=3)
        u0ps=GEO['u0ps']; Vst=GEO['Vst']
        def elt(ch):  # wrap a CE.Ch (bulk only) into B
            return B({(k,()):v for k,v in ch.t.items()})
        U10=elt(GEO['Uk0'][1]); U11=elt(GEO['Uk1'][1]); U20=elt(GEO['Uk0'][2]); U21=elt(GEO['Uk1'][2]); U30=elt(GEO['Uk0'][3])
        Y1=U10.scale(-1/u0ps)
        v0=bexpect(bmul(Y1,Y1,6))
        Y2=bmul(U11,Y1,6).add(U20).scale(-1/u0ps)
        m1=bexpect(Y2)
        Y2c=Y2.add(B.scalar(-m1))
        VarY2=bexpect(bmul(Y2c,Y2c,6))
        Y1sq=bpow(Y1,2,6); Y1cu=bpow(Y1,3,6)
        Y3=(Y1cu.scale(Vst*u0ps/6).add(bmul(U11,Y2,6)).add(bmul(U10,Y1sq,6).scale(0.5*Vst))
            .add(bmul(U21,Y1,6)).add(U30)).scale(-1/u0ps)
        Y1Y3=bexpect(bmul(Y1,Y3,6))
        print(f" n={n}: v0={v0:.4f}(.1339) m1={m1:+.4f}(.212) Var(Y2)={VarY2:.4f}(.0636) <Y1Y3>={Y1Y3:+.4f}(.0231) v1={VarY2+2*Y1Y3:+.4f}(.110)")
