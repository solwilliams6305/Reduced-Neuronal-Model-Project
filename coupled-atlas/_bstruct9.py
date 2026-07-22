import sympy as sp, time
from sympy import binomial
MAX=9; M=MAX+1
U=[[sp.symbols(f'U{k}_{m}') for m in range(M)] for k in range(MAX+1)]
with open('_yexprs_9.txt') as f:
    Yexpr=[sp.sympify(l.strip()) for l in f if l.strip()]
s=[sp.symbols(f's{j}') for j in range(8)]; Vs=[sp.symbols(f'Vd{j}') for j in range(8)]
Ub={}
for k in range(MAX+1):
    Ub[(k,0)]=0 if k==0 else U[k][0]; Ub[(k,1)]=U[k][1]
for m in range(2,MAX+1):
    for k in range(MAX+1):
        p=m-2; tot=0
        for j in range(p+1):
            c=binomial(p,j); tot+=c*Vs[j]*Ub[(k,p-j)]
            if k>=1 and j<8: tot+=c*s[j]*Ub[(k-1,p-j)]
        Ub[(k,m)]=sp.expand(tot)
subs={U[k][m]:Ub[(k,m)] for k in range(MAX+1) for m in range(2,M)}; subs[U[0][0]]=0
Y1s=-U[1][0]/U[0][1]
t0=time.time()
def smono(name):
    print(f'  {name}s',flush=True)
for n in range(4,10):
    Yb=sp.expand(Yexpr[n-1].subs(subs))
    poly=sp.Poly(Yb,*s)
    quad={}; cub={}
    for monom,coeff in poly.terms():
        d=sum(monom)
        if d==2: quad[monom]=coeff
        elif d==3: cub[monom]=coeff
    def lbl(m): return '*'.join('s'+str(i) for i,e in enumerate(m) for _ in range(e))
    print(f"Y{n}: quad-s: {[lbl(m) for m in quad] or 'none'} ;  cubic-s: {[lbl(m) for m in cub] or 'none'}")
    # simplify each quadratic coeff in terms of Y1 = -U1_0/U0_1
    for m,c in quad.items():
        cs=sp.simplify(c)
        # try express as multiple of Y1^p
        print(f"     coeff[{lbl(m)}] = {cs}")
print(f'({time.time()-t0:.0f}s)')
