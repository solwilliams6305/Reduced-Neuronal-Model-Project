"""Hastings-McLeod PII s->-inf asymptotic coefficients + median-resummation Stokes structure
(route R1 for the TW_beta complex Stokes constant; the tail inherits this via Theorem T1).

q'' = 2q^3 + s q,  q ~ sqrt(-s/2) sum_k a_k (-s)^{-3k}  (s->-inf).
Result: a_0=1, a_{k>=1}<0 (NON-alternating) => positive-axis Borel singularity => imaginary
lateral ambiguity => median-resummation (complex) Stokes constant.  Growth a_k ~ -C Gamma(2k-1/2)
(9/8)^k, C~0.139~1/(4 sqrt pi); action A=2 sqrt2/3.
"""
import numpy as np, json
import sympy as sp
from scipy.special import gammaln

def hm_coeffs(N=10):
    u = sp.symbols('u', positive=True)                # u = 1/z = -1/s
    A = list(sp.symbols('a0:%d' % N))
    F = sum(A[k]*u**(3*k) for k in range(N))
    q = sp.sqrt(1/(2*u))*F
    qz = -u**2*sp.diff(q, u); qzz = -u**2*sp.diff(qz, u)
    E = sp.expand((qzz - (2*q**3 - (1/u)*q))*sp.sqrt(2*u))
    ser = sp.series(E, u, 0, 3*N).removeO()
    sol = {A[0]: 1}
    for k in range(1, N):
        expr = sp.expand(ser.subs(sol))
        for t, c in sp.collect(expr, u, evaluate=False).items():
            cc = sp.expand(c)
            if cc.has(A[k]) and not any(cc.has(A[j]) for j in range(k+1, N)):
                s = sp.solve(cc, A[k])
                if s:
                    sol[A[k]] = sp.simplify(s[0]); break
    return [sol[A[k]] for k in range(N)]

if __name__ == "__main__":
    exact = hm_coeffs(10)
    a = np.array([float(x) for x in exact])
    print("a_k (exact):", [str(x) for x in exact[:5]], "...")
    print("signs:", ['+' if v > 0 else '-' for v in a], " (non-alternating for k>=1)")
    k = np.arange(len(a)); r = a[1:]/a[:-1]
    print("ratio/(4k^2) -> 9/8:", np.round(r[1:]/(4*k[1:-1]**2), 4))
    y = np.log(np.abs(a[1:])) - k[1:]*np.log(9/8)
    from scipy.optimize import brentq
    g = brentq(lambda gg: (y-gammaln(2*k[1:]+gg))[-1]-(y-gammaln(2*k[1:]+gg))[-2], -3, 3)
    C = np.exp(np.mean((y-gammaln(2*k[1:]+g))[-3:]))
    print(f"growth: a_k ~ -{C:.4f} Gamma(2k+{g:.3f}) (9/8)^k;  action A=2sqrt2/3={2*np.sqrt(2)/3:.4f}")
    print(f"        C={C:.4f} vs 1/(4 sqrt pi)={1/(4*np.sqrt(np.pi)):.4f}")
    json.dump(list(a), open('_hm_coeffs.json', 'w'))
