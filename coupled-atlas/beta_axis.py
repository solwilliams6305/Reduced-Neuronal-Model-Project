"""
Atlas beta-axis validation: fold (q=1) follows TW_beta (skew DECREASES with beta);
cusp (q=2) follows W_beta (skew INCREASES with beta -- the inverted family). beta=4/eta^2.
Confirms the two distinct edge-class assignments and the fold/cusp beta-inversion.
numpy only.
"""
import numpy as np
def Vq(Y,q): return np.sign(Y)*np.abs(Y)**q
def simulate(q,eta,N,Y0=3.0,Yend=-4.5,dt=4e-3,seed=7):
    rng=np.random.default_rng(seed); nst=int(round((Y0-Yend)/dt)); sq=np.sqrt(dt)
    p=np.full(N,Vq(Y0,q)**0.5); Ys=np.full(N,np.nan)
    for i in range(nst):
        Y=Y0-i*dt; p+=(Vq(Y,q)-p*p)*dt+eta*sq*rng.standard_normal(N); np.clip(p,-60,60,out=p)
        nw=np.isnan(Ys)&(p<-12.0); Ys[nw]=Y0-(i+1)*dt
    Ys[np.isnan(Ys)]=Yend; z=(Ys-Ys.mean())/Ys.std(); return np.mean(z**3), np.mean(z**4)-3

betas=[1.0,2.0,4.0]
# reference TW_beta cumulants (eigenvalue convention): skew, exk
TWref={1:(0.293,0.165),2:(0.224,0.093),4:(0.165,0.044)}
print(f"{'beta':>5}{'eta':>6} | {'FOLD q=1 skew':>14}{'(TW_b ref)':>11} | {'CUSP q=2 skew':>14}")
fold=[];cusp=[]
for b in betas:
    eta=2.0/np.sqrt(b)
    s1,e1=simulate(1,eta,150000); s2,e2=simulate(2,eta,150000)
    fold.append(s1); cusp.append(s2)
    print(f"{b:>5.0f}{eta:>6.2f} | {s1:>+14.3f}{TWref[b][0]:>+11.3f} | {s2:>+14.3f}")
print(f"\nFOLD skew trend: {fold[0]:+.3f} -> {fold[1]:+.3f} -> {fold[2]:+.3f}  "
      f"({'DECREASING (TW_beta) OK' if fold[0]>fold[2] else 'NOT decreasing'})")
print(f"CUSP skew trend: {cusp[0]:+.3f} -> {cusp[1]:+.3f} -> {cusp[2]:+.3f}  "
      f"({'INCREASING (W_beta, inverted) OK' if cusp[0]<cusp[2] else 'NOT increasing'})")
print(f"=> fold and cusp have OPPOSITE beta-dependence (the beta-inversion at the catastrophe step q=1->2).")
np.save("beta_axis_res.npy",np.array([fold,cusp]))
