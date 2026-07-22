import sympy as sp
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
Y1=-U[1][0]/U[0][1]
def linc(n,j):
    Yb=sp.expand(Yexpr[n-1].subs(subs))
    c=Yb.coeff(s[j],1)
    for k in range(8):
        if k!=j: c=c.subs(s[k],0)
    return sp.simplify(c)
# linear s-coeffs needed for v4's Var(Y5) [s0,s1], Cov(Y4,Y6) [Y4:s0, Y6:s0]
print('B4 (Y4 s0):', sp.simplify(linc(4,0)), ' = Y1^3/3?', sp.simplify(linc(4,0)-Y1**3/3)==0)
print('B50(Y5 s0):', sp.simplify(linc(5,0)))
print('B51(Y5 s1):', sp.simplify(linc(5,1)))
print('B60(Y6 s0):', sp.simplify(linc(6,0)))
