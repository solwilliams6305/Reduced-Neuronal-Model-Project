import sympy as sp, time
from sympy import Rational
MAX=13; M=MAX+1
t0=time.time()
U=[[sp.symbols(f'U{k}_{m}') for m in range(M)] for k in range(MAX+1)]
Yv=[0]+[sp.symbols(f'Y{n}') for n in range(1,MAX+1)]
# truncated power series in eta as list of length M (coeff of eta^p), coefficients are sympy exprs
def smul(a,b):
    c=[sp.Integer(0)]*M
    for i in range(M):
        if a[i]==0: continue
        for j in range(M-i):
            if b[j]!=0: c[i+j]=c[i+j]+a[i]*b[j]
    return c
def sadd(a,b): return [a[i]+b[i] for i in range(M)]
# D = sum_{j>=1} eta^j Y_j
D=[sp.Integer(0)]*M
for j in range(1,MAX+1): D[j]=Yv[j]
# precompute D^m for m=0..MAX (truncated)
Dpow=[[sp.Integer(1)]+[sp.Integer(0)]*(M-1)]  # D^0
for m in range(1,MAX+1):
    Dpow.append(smul(Dpow[-1],D))
    print(f'  D^{m} built [{time.time()-t0:.0f}s]',flush=True)
# u = sum_k eta^k sum_m U[k][m] D^m / m!
uS=[sp.Integer(0)]*M
for k in range(MAX+1):
    for m in range(M):
        coef=Rational(1,sp.factorial(m))
        term=[coef*U[k][m]*Dpow[m][p] for p in range(M)]
        # multiply by eta^k = shift by k
        for p in range(M-1,-1,-1):
            uS[p]= uS[p] + (term[p-k] if p-k>=0 else 0)
uS=[sp.expand(x.subs(U[0][0],0)) for x in uS]
print(f'  u series built [{time.time()-t0:.0f}s]',flush=True)
# solve order by order
sol={}; exprs=[]
U01=U[0][1]
for n in range(1,MAX+1):
    eqn=sp.expand(uS[n].subs(sol))
    yn=sp.expand(sp.solve(eqn,Yv[n])[0])
    sol[Yv[n]]=yn; exprs.append(yn)
    print(f'  Y{n}: {len(sp.Add.make_args(yn))} terms [{time.time()-t0:.0f}s]',flush=True)
with open('_yexprs_13.txt','w') as f:
    for e in exprs: f.write(sp.srepr(e)+"\n")
print(f'DONE order 11 [{time.time()-t0:.0f}s]',flush=True)
