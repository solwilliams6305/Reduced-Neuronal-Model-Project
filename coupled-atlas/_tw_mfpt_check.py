"""Minimal check: reduced escape rate R(g) closed form vs direct MFPT double integral.
Finite cutoffs: the MFPT weight e^{V(y)/D}e^{-V(z)/D} is dominated by the barrier y=-1, well z=+1;
integrand negligible beyond |.|=3, and the inner lower-limit blow-up (V->-inf as z->-inf) must be
capped to avoid overflow."""
import mpmath as mp, sys
mp.mp.dps = 30
L = mp.mpf(3)      # finite cutoff (|V/D| already huge here)

def qt(J):
    out = [mp.mpf(1)]
    for j in range(1, J + 1):
        num = mp.mpf(1)
        for i in range(6):
            num *= (6*j - i)
        den = mp.mpf(576)*(3*j)*(3*j-1)*(3*j-2)*(2*j)*(2*j-1)
        out.append(out[-1]*num/den)
    return out
QT = qt(60)
V = lambda t: t**3/3 - t

def R_direct(g):
    D = mp.mpf(2)/g
    def inner(y):
        pts = [y, 1, L] if y < 1 else [y, L]
        return mp.quad(lambda z: mp.e**(-V(z)/D), pts)
    tau = (1/D)*mp.quad(lambda y: mp.e**(V(y)/D)*inner(y), [-L, -1, 1])
    return 1/tau

def R_formula(g):
    D = mp.mpf(2)/g
    terms = [abs(QT[j]*D**j) for j in range(len(QT))]
    K = min(range(1, len(terms)), key=lambda k: terms[k])
    Q = mp.fsum(QT[j]*D**j for j in range(K+1))
    return (1/mp.pi)*mp.e**(-mp.mpf(2)/3*g)/Q**2, K

print("  g     R_direct(MFPT)      R_formula          rel.diff     K*", flush=True)
for g in (mp.mpf(10), mp.mpf(16)):
    Rd = R_direct(g)
    Rf, K = R_formula(g)
    print(f"  {int(g):>2}  {mp.nstr(Rd,12):>18} {mp.nstr(Rf,12):>18}  {mp.nstr(abs(Rd-Rf)/Rd,3):>10}  {K}", flush=True)
