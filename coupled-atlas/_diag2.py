import numpy as np, sympy as sp
from math import comb
import chaos_boundary as CB
from chaos_boundary import B, bmul, bexpect, setup, GEO
def Vd(Y,d): return {0:-Y*Y,1:-2*Y,2:-2.0}.get(d,0.0)
def cov(A,Bx,MT): return bexpect(bmul(A,Bx,MT))-bexpect(A)*bexpect(Bx)

Usym=[[sp.symbols(f'U{k}_{mm}') for mm in range(10)] for k in range(10)]
with open('_yexprs_9.txt') as f: Yexpr=[sp.sympify(l.strip()) for l in f if l.strip()]

def build_Y(n=40, MT=3, zero_s=False):
    setup(n_grid=n,MAXORD=MT); u0ps=GEO['u0ps']; Yst=GEO['s'][0]; Vs=[Vd(Yst,d) for d in range(6)]
    def elt(ch): return B({(k,()):v for k,v in ch.t.items()})
    Umat={}
    for k in range(1,MT+1): Umat[(k,0)]=elt(GEO['Uk0'][k]); Umat[(k,1)]=elt(GEO['Uk1'][k])
    Umat[(0,0)]=B.scalar(0.0); Umat[(0,1)]=B.scalar(u0ps)
    s=[B() if zero_s else B({(0,(j,)):np.array(1.0)}) for j in range(6)]
    for m in range(2,MT+1):
        for k in range(0,MT+1):
            acc=B(); p=m-2
            for j in range(p+1):
                cc=comb(p,j)
                if (k,p-j) in Umat: acc=acc.add(Umat[(k,p-j)].scale(cc*Vs[j]))
                if k>=1 and (k-1,p-j) in Umat: acc=acc.add(bmul(s[j],Umat[(k-1,p-j)],MT).scale(cc))
            Umat[(k,m)]=acc
    def evalY(n_):
        e=sp.expand(Yexpr[n_-1].subs(Usym[0][1],u0ps).subs(Usym[0][0],0)); res=B()
        for term in sp.Add.make_args(e):
            coef,mono=term.as_coeff_Mul(); el=B.scalar(float(coef))
            for fac in sp.Mul.make_args(mono):
                base,exp=fac.as_base_exp(); nm=str(base)
                k=int(nm[1:nm.index('_')]); mm=int(nm[nm.index('_')+1:])
                atom=Umat.get((k,mm),B())
                for _ in range(int(exp)): el=bmul(el,atom,MT)
            res=res.add(el)
        return res
    return {n_:evalY(n_) for n_ in range(1,6)}

for zs in [True,False]:
    Y=build_Y(zero_s=zs); MT=3
    Y3sq=bexpect(bmul(Y[3],Y[3],MT)); C24=cov(Y[2],Y[4],MT); C15=cov(Y[1],Y[5],MT)
    tag="s=0 (reg only)" if zs else "s=ON (full) "
    print(f" {tag}: <Y3^2>={Y3sq:.4f}  2Cov24={2*C24:+.4f}  2Cov15={2*C15:+.4f}  v2={Y3sq+2*C24+2*C15:+.4f}")

print("\n--- MT sensitivity (reg, s=0) and m2=<Y4> ---")
for MT in [3,4,5]:
    Y=build_Y(n=36,MT=MT,zero_s=True)
    C24=cov(Y[2],Y[4],MT)
    m2=bexpect(Y[4])
    print(f"  MT={MT}: 2Cov24_reg={2*C24:+.4f}   m2=<Y4>={m2:+.4f}  (mollified m2~+0.078)")
