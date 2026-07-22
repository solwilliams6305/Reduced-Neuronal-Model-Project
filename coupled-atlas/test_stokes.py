"""
Diagnose the missing ingredient. The naive Wronskian gives condition
  e^{-ipi/4}(1-S) R(i lam) = (1+S) R(lam),   R(w)=G(3/4-w/4)/G(1/4-w/4),
where S is the inverted-side Stokes multiplier of the outgoing U-combination (naive form had S=0).
Solve for the S REQUIRED to make D=0 at the true (ODE) resonances, and compare to the PREDICTED
Stokes factor S_pred = sqrt(2pi)/Gamma(1/2 - i lam/2) (the oscillatory-side data). If they match
(up to a fixed phase), the closed form is D=0 <=> that condition -> M3 essentially closed.
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
def R(w): return cmath.exp(clogG(0.75-w/4)-clogG(0.25-w/4))
E=cmath.exp(-1j*cmath.pi/4)
res=[0.86-0.82j, 2.30-1.22j, 4.14-1.08j]   # ACCURATE ODE resonances
print("at each ODE resonance lam:  S_required (for naive-form D=0)  vs  S_pred=sqrt(2pi)/Gamma(1/2-i lam/2)")
for lam in res:
    Sreq=(E*R(1j*lam)-R(lam))/(E*R(1j*lam)+R(lam))
    Spred=cmath.sqrt(2*cmath.pi)/G(0.5-1j*lam/2)
    print(f"  lam={lam.real:+.2f}{lam.imag:+.2f}i: "
          f"S_req=({Sreq.real:+.3f}{Sreq.imag:+.3f}i,|{abs(Sreq):.3f}|,arg{cmath.phase(Sreq):+.2f})  "
          f"S_pred=(|{abs(Spred):.3f}|,arg{cmath.phase(Spred):+.2f})  ratio|.|={abs(Sreq)/abs(Spred):.3f}")
print("\n(if |S_req|/|S_pred| ~ const and arg differ by a const => same Stokes structure, closed form fixed up to that phase)")
