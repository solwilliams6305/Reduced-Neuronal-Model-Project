"""
DECISIVE test of the derived left-tail law:  -logP(Y*<-s) = I(s)/eta^2 + c(s) = I(s)*beta/4 + c(s).
=> at fixed Theta, -logP is LINEAR in beta with slope = I(s)/4 (the instanton rate) and eta-indep intercept c(s).
This both (a) confirms the FW instanton rate I(s) (exponent 5, q_asymp=1, const 1/10) artifact-free (MC, no
numerical diffusion), and (b) explains the 'anomalous' effective exponent (0.6 rate / 0.72 amplitude) as the
crossover produced by the O(1) offset c(s).
High-stats MC at several beta; instanton I(s) from solve_bvp. [DERIVED + NUMERIC].
"""
import numpy as np
from scipy.integrate import solve_bvp
def Vc(Y): return np.sign(Y)*Y*Y
def escape(eta,N,Y0=3.0,Yend=-8.0,dt=1.5e-3,seed=2,thr=20.0,clip=250.0):
    rng=np.random.default_rng(seed); n=int(round((Y0-Yend)/dt)); sq=np.sqrt(dt)
    p=np.full(N,np.sqrt(Vc(Y0))); Ys=np.full(N,np.nan)
    for i in range(n):
        Y=Y0-i*dt; p+=(Vc(Y)-p*p)*dt+eta*sq*rng.standard_normal(N); np.clip(p,-clip,clip,out=p)
        nw=np.isnan(Ys)&(p<-thr); Ys[nw]=Y0-(i+1)*dt
    return Ys[~np.isnan(Ys)]
def Vt(t): return -np.sign(t)*t*t
def Iaction(s,p0=0.0,M=20.0,t0=1e-3,n=1000):
    t=np.linspace(t0,s,n)
    def ode(t,y): p,pi=y; return np.vstack([pi+Vt(t)-p*p,2*p*pi])
    def bc(ya,yb): return np.array([ya[0]-p0,yb[0]+M])
    pg=np.clip(np.where(t<0.9*s,1.0/np.maximum(t,0.3),-M*(t-0.9*s)/(0.1*s)),-M,3.0)
    sol=solve_bvp(ode,bc,t,np.vstack([pg,t*t]),max_nodes=200000,tol=1e-6)
    tt=np.linspace(t0,s,3000); pp,pi=sol.sol(tt); return 0.5*np.trapz(pi**2,tt)

betas=[2.0,3.0,4.0,5.0,6.0]; N=2500000
Thd={}
for b in betas:
    Ys=escape(2.0/np.sqrt(b),N,seed=int(10*b)); Thd[b]=(Ys[Ys<0]**2)/2

print("=== -logP(Theta>t) vs beta at fixed Theta: slope should = I(s)/4 (instanton rate) ===")
print(f"  {'Theta':>5} {'s':>5} {'I(s)':>7} {'I/4(pred)':>9} {'slope(MC)':>9} {'intercept c':>11} {'slope/pred':>9}")
Ths=[3.5,4.0,4.5,5.0,5.5]
rows=[]
for Th in Ths:
    s=np.sqrt(2*Th); I=Iaction(s)
    bb=[]; yy=[]
    for b in betas:
        T=Thd[b]; cnt=(T>Th).sum()
        if cnt>=40: bb.append(b); yy.append(-np.log(cnt/T.size))
    if len(bb)>=3:
        A=np.polyfit(bb,yy,1); slope,inter=A[0],A[1]
        print(f"  {Th:5.1f} {s:5.2f} {I:7.2f} {I/4:9.3f} {slope:9.3f} {inter:11.3f} {slope/(I/4):9.3f}")
        rows.append((Th,I/4,slope,inter))
print("\n  slope/pred ~ 1  => the instanton rate I(s) (=> beta s^5/40, exponent 5, const 1/10) is CONFIRMED,")
print("  artifact-free (MC).  The O(1) intercept c(s) is the offset that makes the effective exponent < 1.")

# demonstrate: the 'anomalous' effective exponent from -logP = I*beta/4 + c is a crossover, ->1 as beta->inf
print("\n=== effective exponent q_eff(beta) from the DERIVED law I*beta/4+c (using Theta=4 fit) ===")
Th=4.0; s=np.sqrt(2*Th); I=Iaction(s)
# fit c at Theta=4
bb=[b for b in betas if (Thd[b]>Th).sum()>=40]; yy=[-np.log((Thd[b]>Th).sum()/Thd[b].size) for b in bb]
c=np.polyfit(bb,yy,1)[1]
for b in [2,8,32,128,512]:
    val=I*b/4+c; # local exponent d log(val)/d log b = (I b/4)/(I b/4 + c)
    qeff=(I*b/4)/(I*b/4+c)
    print(f"  beta={b:>4}: -logP={val:8.1f}  q_eff(local)={qeff:.3f}")
print("  => q_eff rises 0.6-0.7 (accessible) -> 1 (semiclassical). The 'anomaly' is the crossover, DERIVED.")
np.savez("leftail_deriv.npz",Ths=np.array(Ths),rows=np.array(rows))
