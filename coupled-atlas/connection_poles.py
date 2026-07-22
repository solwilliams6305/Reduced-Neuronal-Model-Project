"""
DECISIVE non-confounded test of the c=1/Weber Gamma connection: the confining side u''=(Y^2-lambda)u is an
inverted-well with bound states at lambda=1,3,5,... (odd integers) = the POLES of Gamma(1/4-lambda/4) &
Gamma(3/4-lambda/4). A true Gamma connection makes the oscillatory-side connection phase JUMP by pi at each of
these (and adds a Y>0 node); a smooth WKB/linear phase CANNOT. So: scan lambda across 0->6, track the number of
Y>0 nodes of the continued solution and the first oscillatory node, and check the jumps land on odd integers.
scipy. [NUMERIC — decisive discriminator].
"""
import numpy as np
from scipy.integrate import solve_ivp
def V(Y,lam): return np.sign(Y)*Y*Y - lam
def scan(lam, Y0=6.5, Yend=-8.0):
    V0=V(Y0,lam); def_=None
    def rhs(Y,y): return [y[1], V(Y,lam)*y[0]]
    sol=solve_ivp(rhs,[Y0,Yend],[1.0,np.sqrt(max(V0,1e-9))],max_step=1.5e-3,rtol=1e-11,atol=1e-13,dense_output=True)
    Yg=np.linspace(Y0,Yend,24000); u,_=sol.sol(Yg)
    npos=0; nodes_neg=[]
    for i in range(1,len(Yg)):
        if u[i-1]*u[i]<0:
            Yc=Yg[i-1]+ (u[i-1]/(u[i-1]-u[i]))*(Yg[i]-Yg[i-1])
            if Yc>0: npos+=1
            else: nodes_neg.append(Yc)
    Th1=(nodes_neg[0]**2/2) if nodes_neg else np.nan
    return npos, Th1

print("=== bound-state / pole scan: do Y>0 nodes increment at odd integers (Gamma poles)? ===")
print(f"  {'lambda':>7} {'#(Y>0 nodes)':>12} {'Theta1 (1st osc node depth)':>28}")
lams=np.concatenate([np.arange(0.0,8.05,0.25), np.array([2.9,2.95,2.99,3.01,6.9,6.95,6.99,7.01,7.05])])
lams=np.sort(lams)
prev=None; jumps=[]
for lam in lams:
    npos,Th1=scan(lam)
    mark=""
    if prev is not None and npos!=prev: mark=f"  <-- jump {prev}->{npos}"; jumps.append(lam)
    print(f"  {lam:7.2f} {npos:12d} {Th1:28.3f}{mark}")
    prev=npos
print(f"\n  node-increment (pole) locations: {[f'{j:.2f}' for j in jumps]}")
print(f"  predicted Gamma poles (confining bound states) = odd integers 1,3,5 ...")
if jumps:
    near=[round(2*round((j-1)/2)+1) for j in jumps]
    dev=np.mean([abs(j-n) for j,n in zip(jumps,near)])
    print(f"  nearest odd integers: {near}; mean deviation {dev:.2f}  => {'MATCH (Gamma pole structure confirmed)' if dev<0.25 else 'no clean match'}")
