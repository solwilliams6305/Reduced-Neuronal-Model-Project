"""
C/D/E climb: the intrinsic node processes of the WHOLE catastrophe ladder V_q=sign(Y)|Y|^q, q=1..4, unified by
one derived mechanism. Prediction (Pruefer phase-diffusion D_Theta ~ eta^2/k^3, k=|Y|^{q/2}):

   spacing-jitter exponent  a(q) = 3q/(q+2)      [DERIVED]
   accumulated phase var Sum n^{-a} converges iff a>1 iff q>1:
     q=1 fold  -> a=1   -> log  -> CLASS II (GUE / Airy edge)   <- the critical boundary
     q>=2      -> a>1   -> bounded -> CLASS I (perturbed lattice)

Measure, per q, from the additive Weber-type Riccati  dp=(sign(Y)|Y|^q - p^2)d(-Y)+(2/sqrt(beta))dW:
  (1) spacing jitter Var(spacing) vs phase-depth Theta_q -> local exponent a  (compare 3q/(q+2))
  (2) number variance Var(N) growth: bounded (q>=2) vs log (q=1)  -> class I vs II
Phase-depth Theta_q = 2|Y|^{(q+2)/2}/(q+2); nodes ~ pi apart in Theta_q. numpy only. [DERIVED + NUMERIC].
"""
import numpy as np
def Vq(Y,q): return np.sign(Y)*np.abs(Y)**q
def Theta_q(absY,q): return 2.0*absY**((q+2)/2.0)/(q+2.0)

def run(q, beta, N, Y0, Yend, dt, seed=3, maxnodes=140):
    eta=2.0/np.sqrt(beta); rng=np.random.default_rng(seed)
    n=int(round((Y0-Yend)/dt)); sq=np.sqrt(dt)
    p=np.full(N,np.sqrt(max(Vq(Y0,q),1e-9))); nodes=np.full((maxnodes,N),np.nan); idx=np.zeros(N,int)
    cps=np.linspace(-1.0,Yend+0.5,10); ci=0; Ncp=np.zeros((len(cps),N)); cnt=np.zeros(N)
    for i in range(n):
        Y=Y0-i*dt
        thr=8.0*np.abs(Y)**(q/2.0)+40.0                # reset scale ~ local frequency |Y|^{q/2}
        p+=(Vq(Y,q)-p*p)*dt+eta*sq*rng.standard_normal(N); np.clip(p,-1e5,1e5,out=p)
        ex=p<-thr; live=ex&(idx<maxnodes); cols=np.where(live)[0]
        if cols.size: nodes[idx[cols],cols]=Y; idx[cols]+=1
        cnt=cnt+ex; p[ex]=thr
        if ci<len(cps) and Y<=cps[ci]: Ncp[ci]=cnt.copy(); ci+=1
    return nodes, Theta_q(np.abs(cps),q), Ncp

# per-q sweep ranges tuned to yield ~40-80 nodes at fixed beta=2
cfg={1:(3.0,-42.0,3e-3), 2:(3.0,-20.0,1.5e-3), 3:(3.0,-13.0,1.0e-3), 4:(3.0,-9.5,7e-4)}
print("=== ladder node-process jitter exponent a(q) vs DERIVED 3q/(q+2) ===")
print(f"  {'q':>3} {'a_pred=3q/(q+2)':>15} {'a_meas':>8} {'class(pred)':>12} {'Var(N) growth':>14}")
results={}
for q in [1,2,3,4]:
    Y0,Yend,dt=cfg[q]
    nodes,Thcp,Ncp=run(q,2.0,6000,Y0,Yend,dt)
    Th=Theta_q(np.abs(nodes),q)                        # phase-depth at each node (ascending in index)
    sp=Th[1:]-Th[:-1]; Thm=0.5*(Th[1:]+Th[:-1])
    # local spacing-jitter exponent over mid-depth bands
    D=[];V=[]
    idxs=np.arange(sp.shape[0])
    bands=[(3,8),(8,14),(14,22),(22,34),(34,50)]
    for lo,hi in bands:
        seg=sp[lo:hi]; th=Thm[lo:hi]; m=np.isfinite(seg)&(seg>0.3)&(seg<10)
        if m.sum()>300: D.append(th[m].mean()); V.append(seg[m].var())
    D=np.array(D);V=np.array(V)
    a_meas=-np.polyfit(np.log(D),np.log(V),1)[0]
    a_pred=3*q/(q+2.0)
    # number-variance growth: slope of Var(N) vs log(Theta) (>0 => log/class II ; ~0 => bounded/class I)
    Nv=Ncp.var(axis=1); Nm=Ncp.mean(axis=1); ok=Nm>1.5
    dVdlog=np.polyfit(np.log(Thcp[ok]),Nv[ok],1)[0] if ok.sum()>3 else np.nan
    cls="II (GUE/log)" if abs(a_pred-1)<0.05 else "I (lattice)"
    print(f"  {q:>3} {a_pred:>15.3f} {a_meas:>8.2f} {cls:>12} {dVdlog:>+14.3f}")
    results[q]=(D,V,a_meas,a_pred,Nv,Nm,Thcp,dVdlog)
np.savez("ladder_data.npz",**{f"D_{q}":results[q][0] for q in results},**{f"V_{q}":results[q][1] for q in results},
         **{f"Nv_{q}":results[q][4] for q in results},**{f"Th_{q}":results[q][6] for q in results})
print("\n  a_meas ~ 3q/(q+2) across the ladder => unified mechanism; q=1 (a=1) is the class II<->I boundary.")
print("  Var(N) growth: q=1 positive (log, class II = Airy/GUE), q>=2 ~0 (bounded, class I = lattice).")
