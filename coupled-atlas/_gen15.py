"""Extend the Y-functional ladder to order 15 (needed for v7): reuse Y1..Y13 from
_yexprs_13.txt, solve orders 14 and 15 with TERM-BY-TERM substitution+expand (the trick that
avoids the giant-intermediate OOM that killed a naive order-12 run).  Writes _yexprs_15.txt."""
import sympy as sp, time, resource
from sympy import Rational

MAX = 15; M = MAX + 1
t0 = time.time()
U = [[sp.symbols(f'U{k}_{m}') for m in range(M)] for k in range(MAX + 1)]
Yv = [0] + [sp.symbols(f'Y{n}') for n in range(1, MAX + 1)]

def rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30

def smul(a, b):
    c = [sp.Integer(0)] * M
    for i in range(M):
        if a[i] == 0: continue
        for j in range(M - i):
            if b[j] != 0: c[i + j] = c[i + j] + a[i] * b[j]
    return c

# D = sum_{j>=1} eta^j Y_j  (symbolic placeholders -- cheap; blowup only happens at subs time)
D = [sp.Integer(0)] * M
for j in range(1, MAX + 1): D[j] = Yv[j]
Dpow = [[sp.Integer(1)] + [sp.Integer(0)] * (M - 1)]
for m in range(1, MAX + 1):
    Dpow.append(smul(Dpow[-1], D))
    print(f'  D^{m} built [{time.time()-t0:.0f}s, {rss_gb():.2f}GB]', flush=True)

# u = sum_k eta^k sum_m U[k][m] D^m / m!
uS = [sp.Integer(0)] * M
for k in range(MAX + 1):
    for m in range(M):
        coef = Rational(1, sp.factorial(m))
        term = [coef * U[k][m] * Dpow[m][p] for p in range(M)]
        for p in range(M - 1, -1, -1):
            uS[p] = uS[p] + (term[p - k] if p - k >= 0 else 0)
uS = [sp.expand(x.subs(U[0][0], 0)) for x in uS]
print(f'  u series built [{time.time()-t0:.0f}s, {rss_gb():.2f}GB]', flush=True)

# load Y1..Y13 (already solved)
with open('_yexprs_13.txt') as f:
    exprs = [sp.sympify(l.strip()) for l in f if l.strip()]
assert len(exprs) == 13
sol = {Yv[j]: exprs[j - 1] for j in range(1, 14)}
print(f'  Y1..Y13 loaded [{time.time()-t0:.0f}s, {rss_gb():.2f}GB]', flush=True)

U01 = U[0][1]
for n in [14, 15]:
    # uS[n] contains Yv[n] ONLY via the k=0,m=1 term U01*Yv[n]; everything else involves Y_{<n}.
    # Solve U01*Y_n + rest = 0  =>  Y_n = -rest/U01, substituting term by term.
    pieces = []
    tlist = sp.Add.make_args(uS[n])
    for i, t in enumerate(tlist):
        if t.has(Yv[n]):
            assert sp.expand(t - U01 * Yv[n]) == 0, f"Y{n} appears nonlinearly: {t}"
            continue
        e = sp.expand(t.subs(sol) / U01)
        pieces.append(-e)
        if (i + 1) % 50 == 0:
            print(f'    Y{n}: term {i+1}/{len(tlist)} [{time.time()-t0:.0f}s, {rss_gb():.2f}GB]', flush=True)
    yn = sp.Add(*pieces)
    yn = sp.expand(yn)     # final merge of like monomials (pieces already expanded -> cheap-ish)
    sol[Yv[n]] = yn
    exprs.append(yn)
    print(f'  Y{n}: {len(sp.Add.make_args(yn))} terms [{time.time()-t0:.0f}s, {rss_gb():.2f}GB]', flush=True)

with open('_yexprs_15.txt', 'w') as f:
    for e in exprs: f.write(sp.srepr(e) + "\n")
print(f'DONE order 15 [{time.time()-t0:.0f}s, {rss_gb():.2f}GB]', flush=True)
