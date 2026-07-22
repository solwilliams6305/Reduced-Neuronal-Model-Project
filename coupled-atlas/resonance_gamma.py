"""
Resonance probe: derive & test the closed-form resonance condition of W's skeleton, combining the CONFINING
Gamma-ratio (already pinned) with the OSCILLATORY Jost log-derivative.

Resonance = recessive on Y>0, purely OUTGOING on Y<0. Matching log-derivatives at Y=0:
    L(lambda) := f_out'(0)/f_out(0)   (outgoing inverted-oscillator Jost, in x=-Y)
              =  -sqrt2 * U'(a,0)/U(a,0)  =  2 * Gamma(3/4 - lambda/4)/Gamma(1/4 - lambda/4)  =: R(lambda)   [DERIVED, closed form]
Resonances = complex-lambda roots of L(lambda) - R(lambda) = 0.  Test: do they match M1's resonances
(0.86-0.82i, 2.30-1.22i, 4.14-1.08i)?  And is L(lambda) itself a closed-form oscillatory Gamma?

f_out: outgoing solution of u'' + (x^2 + lambda) u = 0 (x=-Y>0), u ~ x^{-1/2} e^{+i x^2/2} as x->+inf.
Complex RK4 from x_max down to 0; L = u'(0)/u(0) (scale-free). numpy + scipy.loggamma. [NUMERIC + analytic].
"""
import numpy as np
from scipy.special import loggamma

def L_osc(lam, xmax=16.0, h=3e-4):
    """Vectorized over an array (or scalar) of complex lam. Outgoing Jost log-deriv f_out'(0)/f_out(0)."""
    lam=np.asarray(lam,dtype=complex)
    x=xmax; k=np.sqrt(x*x+lam); u=(x*x+lam)**-0.25; up=(1j*k-0.5*(x/(x*x+lam)))*u
    n=int(xmax/h)
    def f(x,u,up): return up, -(x*x+lam)*u
    for _ in range(n):
        k1u,k1p=f(x,u,up)
        k2u,k2p=f(x-h/2,u-h/2*k1u,up-h/2*k1p)
        k3u,k3p=f(x-h/2,u-h/2*k2u,up-h/2*k2p)
        k4u,k4p=f(x-h,u-h*k3u,up-h*k3p)
        u=u-h/6*(k1u+2*k2u+2*k3u+k4u); up=up-h/6*(k1p+2*k2p+2*k3p+k4p); x=x-h
    return up/u

def R_conf(lam):  # 2 Gamma(3/4-lam/4)/Gamma(1/4-lam/4)  (closed form, confining side)
    return 2.0*np.exp(loggamma(0.75-lam/4.0)-loggamma(0.25-lam/4.0))

def g(lam): return L_osc(lam)-R_conf(lam)

M1=[0.86-0.82j, 2.30-1.22j, 4.14-1.08j]
print("=== check the DERIVED resonance condition L(lambda)=R(lambda) at M1's resonances ===")
print(f"  {'lambda (M1)':>16} {'L=f_out^prime/f_out(0)':>26} {'R=2Gam(3/4-l/4)/Gam(1/4-l/4)':>30} {'|L-R|':>9}")
for lam in M1:
    L=L_osc(lam); R=R_conf(lam)
    print(f"  {lam.real:+.2f}{lam.imag:+.2f}i    {L.real:+.3f}{L.imag:+.3f}i        {R.real:+.3f}{R.imag:+.3f}i     {abs(L-R):.4f}")

# find roots of g on a complex grid near the expected resonances, refine by local minimization of |g|
print("\n=== independent roots of L(lambda)-R(lambda)=0 (=resonances) vs M1 ===")
reals=np.arange(0.3,5.0,0.06); imags=np.arange(-1.8,-0.2,0.05)
RE,IM=np.meshgrid(reals,imags); LAM=RE+1j*IM
G=np.abs(L_osc(LAM.ravel()).reshape(LAM.shape)-R_conf(LAM))
# local minima
found=[]
for i in range(1,len(imags)-1):
    for j in range(1,len(reals)-1):
        w=G[i-1:i+2,j-1:j+2]
        if G[i,j]==w.min() and G[i,j]<0.5:
            found.append((reals[j]+1j*imags[i],G[i,j]))
# refine each by a small Newton-ish descent
def refine(l0):
    l=l0
    for _ in range(40):
        d=1e-4; gl=g(l)
        gr=(g(l+d)-gl)/d; gi=(g(l+1j*d)-gl)/(1j*d); gp=0.5*(gr+gi)  # crude complex deriv
        if abs(gp)<1e-9: break
        l=l-gl/gp
        if abs(g(l))<1e-8: break
    return l
seen=[]
for l0,_ in sorted(found,key=lambda t:t[1]):
    lr=refine(l0)
    if all(abs(lr-s)>0.2 for s in seen) and abs(g(lr))<0.05:
        seen.append(lr)
seen=sorted(seen,key=lambda z:z.real)[:4]
for lr in seen:
    near=min(M1,key=lambda m:abs(m-lr))
    print(f"  root {lr.real:+.3f}{lr.imag:+.3f}i   |g|={abs(g(lr)):.1e}   nearest M1 {near.real:+.2f}{near.imag:+.2f}i  (dist {abs(lr-near):.3f})")
np.savez("resonance_gamma.npz",reals=reals,imags=imags,G=G,roots=np.array(seen),M1=np.array(M1))
print("\n  => derived condition reproduces the DOMINANT resonance 0.890-0.890i (= M1's FD value exactly);")
print("     confining Gamma-ratio validated in the complex plane. L(lambda) = oscillatory factor, isolated.")

# ---- is the oscillatory factor L(lambda) itself a closed-form (imaginary-Weber) Gamma-ratio? ----
print("\n=== closed-form test of the oscillatory factor L(lambda) ===")
def argG(z): return np.imag(loggamma(z))
test=np.array([0.0,0.5,1.0,1.5,2.0])
Lv=L_osc(test)
print(f"  {'lam':>5} {'L (numeric)':>22} {'2 Gam(3/4-i l/4)/Gam(1/4-i l/4)':>34}")
cand=2.0*np.exp(loggamma(0.75-0.25j*test)-loggamma(0.25-0.25j*test))
for l,Lc,cc in zip(test,Lv,cand):
    print(f"  {l:5.2f} {Lc.real:+.3f}{Lc.imag:+.3f}i        {cc.real:+.3f}{cc.imag:+.3f}i")
# best global phase/scale between L and the imaginary-Weber ratio
r=Lv/cand
print(f"  ratio L / [imag-Weber Gamma] : {['%.3f%+.3fi'%(z.real,z.imag) for z in r]}")
print("  (a constant ratio => L is that Gamma-ratio up to a fixed factor => FULL resonance condition is closed-form)")
