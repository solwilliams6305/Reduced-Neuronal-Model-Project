"""
Closed-form resonance determinant D(lambda) of H=-d^2/dY^2+sign(Y)Y^2, verified against the ODE resonances.
Derivation (Wronskian of the two Jost solutions at Y=0, DLMF 12.2.6-7 zero-values):
  D(lam) ∝  e^{-i pi/4}/(G(3/4-lam/4)G(1/4-i lam/4)) - 1/(G(1/4-lam/4)G(3/4-i lam/4))     [1/G entire]
  ZEROS  <=>  R(lam) = e^{-i pi/4} R(i lam),   R(w)=Gamma(3/4-w/4)/Gamma(1/4-w/4).
Confining ratio R(lam) (real arg) glued to oscillatory R(i lam) (imag arg) by interface phase e^{-i pi/4}.
"""
import cmath
_c=[0.99999999999980993,676.5203681218851,-1259.1392167224028,771.32342877765313,
    -176.61502916214059,12.507343278686905,-0.13857109526572012,9.9843695780195716e-6,1.5056327351493116e-7]
def clogG(z):
    z=complex(z)
    if z.real<0.5: return cmath.log(cmath.pi)-cmath.log(cmath.sin(cmath.pi*z))-clogG(1-z)
    z-=1; x=_c[0]
    for i in range(1,9): x+=_c[i]/(z+i)
    t=z+7.5; return 0.5*cmath.log(2*cmath.pi)+(z+0.5)*cmath.log(t)-t+cmath.log(x)
def G(z): return cmath.exp(clogG(z))
print("gamma checks: G(1)=%.6f  G(.5)=%.6f(=sqrtpi %.6f)  G(5)=%.4f"%(G(1).real,G(0.5).real,cmath.pi**.5,G(5).real))
E=cmath.exp(-1j*cmath.pi/4)
def R(w): return cmath.exp(clogG(0.75-w/4)-clogG(0.25-w/4))   # ratio via log-diff (no overflow)
def F(lam): return R(lam)-E*R(1j*lam)          # resonance: F=0
def newton(l0,it=80):
    l=complex(l0)
    for _ in range(it):
        h=1e-6; d=(F(l+h)-F(l-h))/(2*h)
        if abs(d)<1e-14: break
        step=F(l)/d
        if abs(step)>0.4: step=step/abs(step)*0.4   # damp
        l=l-step
        if abs(step)<1e-11: break
    return l
print("\nclosed-form resonances (Newton from the ODE-found values):")
ode=[1.64-1.49j, 3.6-1.5j]
for g0 in ode:
    z=newton(g0)
    print(f"  start {g0:+.2f} -> D-zero lam = {z.real:+.4f}{z.imag:+.4f}i   |F|={abs(F(z)):.2e}")
# independent coarse scan for the lowest zeros (no seeding from ODE)
print("\nindependent grid scan (|F| minima, lower half-plane):")
import numpy as np
best=[]
for E_ in np.linspace(0.5,4.5,81):
    for Gg in np.linspace(0.3,2.4,43):
        lam=E_-1j*Gg; best.append((abs(F(lam)),E_,Gg))
best.sort()
seen=[]
for v,e,g in best:
    if all(abs(e-s[0])>0.5 or abs(g-s[1])>0.5 for s in seen):
        z=newton(e-1j*g); seen.append((e,g))
        print(f"   min near {e:.2f}-{g:.2f}i -> refined lam={z.real:+.4f}{z.imag:+.4f}i  |F|={abs(F(z)):.1e}")
    if len(seen)>=4: break
