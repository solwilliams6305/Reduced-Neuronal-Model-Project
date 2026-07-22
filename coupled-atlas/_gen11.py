import sympy as sp, time
from sympy import binomial
MAX=11; M=MAX+1
t0=time.time()
U=[[sp.symbols(f'U{k}_{m}') for m in range(M)] for k in range(MAX+1)]
Y=[0]+[sp.symbols(f'Y{n}') for n in range(1,MAX+1)]
eta=sp.symbols('eta')
D=sum(eta**k*Y[k] for k in range(1,MAX+1))
expr=sp.expand(sum(eta**k*U[k][m]*D**m/sp.factorial(m) for k in range(MAX+1) for m in range(M))).subs(U[0][0],0)
sol={}; exprs=[]
for n in range(1,MAX+1):
    c=sp.expand(sp.expand(expr).coeff(eta,n).subs(sol))
    yn=sp.expand(sp.solve(c,Y[n])[0]); sol[Y[n]]=yn; exprs.append(yn)
    print(f'  Y{n}: {len(sp.Add.make_args(yn))} terms  [{time.time()-t0:.0f}s]',flush=True)
with open('_yexprs_11.txt','w') as f:
    for e in exprs: f.write(sp.srepr(e)+"\n")
print(f'DONE order 11 in {time.time()-t0:.0f}s',flush=True)
