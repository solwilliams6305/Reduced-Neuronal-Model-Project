"""Validate the boundary chaos engine on v2 (~0.11): build U_{k,m} with boundary noise s_j, evaluate Y3,Y4,Y5,
compute v2 = <Y3^2> + 2Cov(Y2,Y4) + 2Cov(Y1,Y5).  Truncate chaos at order 3 (suffices for v2)."""
import numpy as np, sympy as sp
from math import comb
import chaos_boundary as CB
from chaos_boundary import B, bmul, bexpect, setup, GEO

def Vd(Y,d): return {0:-Y*Y,1:-2*Y,2:-2.0}.get(d,0.0)

def cov(A,Bx,MT):
    mA=bexpect(A); mB=bexpect(Bx)
    return bexpect(bmul(A,Bx,MT))-mA*mB

def run(n=40, MT=3):
    setup(n_grid=n,MAXORD=MT)
    u0ps=GEO['u0ps']; Yst=GEO['s'][0]
    Vs=[Vd(Yst,d) for d in range(6)]
    def elt(ch): return B({(k,()):v for k,v in ch.t.items()})
    # base U_{k,0},U_{k,1}
    Umat={}
    for k in range(1,MT+1):
        Umat[(k,0)]=elt(GEO['Uk0'][k]); Umat[(k,1)]=elt(GEO['Uk1'][k])
    Umat[(0,0)]=B.scalar(0.0); Umat[(0,1)]=B.scalar(u0ps)
    s=[B({(0,(j,)):np.array(1.0)}) for j in range(6)]   # boundary noises s_j
    # recursion U_{k,m}=sum_j C(m-2,j)[V^(j)U_{k,m-2-j} + s_j U_{k-1,m-2-j}]
    for m in range(2,MT+1):
        for k in range(0,MT+1):
            acc=B(); p=m-2
            for j in range(p+1):
                c=comb(p,j)
                if (k,p-j) in Umat: acc=acc.add(Umat[(k,p-j)].scale(c*Vs[j]))
                if k>=1 and (k-1,p-j) in Umat:
                    acc=acc.add(bmul(s[j],Umat[(k-1,p-j)],MT).scale(c))
            Umat[(k,m)]=acc
    # load Y-formulas
    Usym=[[sp.symbols(f'U{k}_{mm}') for mm in range(10)] for k in range(10)]
    with open('_yexprs_9.txt') as f:
        Yexpr=[sp.sympify(l.strip()) for l in f if l.strip()]
    def evalY(n_):
        e=Yexpr[n_-1].subs(Usym[0][1],u0ps).subs(Usym[0][0],0)
        poly=sp.expand(e)
        res=B()
        for term in sp.Add.make_args(poly):
            coef,mono=term.as_coeff_Mul()
            el=B.scalar(float(coef))
            for fac in sp.Mul.make_args(mono):
                base,exp=fac.as_base_exp()
                # base is some U{k}_{m}
                nm=str(base)
                k=int(nm[1:nm.index('_')]); mm=int(nm[nm.index('_')+1:])
                atom=Umat.get((k,mm),B())   # missing (k>MT) => pure high chaos, zero at this truncation
                for _ in range(int(exp)):
                    el=bmul(el,atom,MT)
            res=res.add(el)
        return res
    Y1=evalY(1); Y2=evalY(2); Y3=evalY(3); Y4=evalY(4); Y5=evalY(5)
    Y3sq=bexpect(bmul(Y3,Y3,MT))
    Cov24=cov(Y2,Y4,MT); Cov15=cov(Y1,Y5,MT)
    v2=Y3sq+2*Cov24+2*Cov15
    return dict(v0=bexpect(bmul(Y1,Y1,MT)),m1=bexpect(Y2),Y3sq=Y3sq,Cov24=Cov24,Cov15=Cov15,v2=v2)

if __name__=="__main__":
    print("boundary engine -> v2 (target ~0.11; decomposition <Y3^2>~0.145, 2Cov24~0.008, 2Cov15~-0.044)")
    for n in [30,40,50]:
        d=run(n=n,MT=3)
        print(f"  n={n}: v0={d['v0']:.4f} m1={d['m1']:+.4f} | <Y3^2>={d['Y3sq']:.4f} 2Cov24={2*d['Cov24']:+.4f} 2Cov15={2*d['Cov15']:+.4f}  v2={d['v2']:+.4f}")
