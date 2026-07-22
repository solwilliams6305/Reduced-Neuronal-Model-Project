import numpy as np
def solve_fp_eta(eta_fn, Y0=3.0, pmin=-8.0, pmax=7.0, dp=0.02, dt=4e-5, tau_max=7.5, deep=False):
    if deep: pmin,pmax,dp,dt,tau_max=-13.0,8.0,0.025,8e-5,11.0
    pc=np.arange(pmin+dp/2,pmax,dp); pf=np.arange(pmin,pmax+dp/2,dp); pf2=pf*pf
    p0=np.sqrt(max(np.sign(Y0)*Y0*Y0,1e-9)); rho=np.exp(-0.5*((pc-p0)/0.30)**2); rho/=rho.sum()*dp
    nst=int(round(tau_max/dt)); taus=np.empty(nst+1); S=np.empty(nst+1); taus[0]=0; S[0]=rho.sum()*dp
    for k in range(nst):
        Y=Y0-k*dt; V=np.sign(Y)*Y*Y; D=eta_fn(Y)**2/2.0; drift=V-pf2
        rl=np.concatenate(([0.0],rho)); rr=np.concatenate((rho,[0.0]))
        J=np.where(drift>0,rl,rr)*drift-D*(rr-rl)/dp; J[-1]=0.0
        rho=rho-dt*(J[1:]-J[:-1])/dp; np.maximum(rho,0.0,out=rho); taus[k+1]=(k+1)*dt; S[k+1]=rho.sum()*dp
    y=(Y0-taus)[::-1]; F=np.clip(S[::-1],1e-300,1.0); return y,F
def cums(y,F):
    f=np.gradient(np.clip(F,0,1),y); f=np.clip(f,0,None); f/=np.trapz(f,y)
    m1=np.trapz(y*f,y); c=y-m1; v=np.trapz(c*c*f,y); sd=np.sqrt(v); z=c/sd
    return m1,sd,np.trapz(z**3*f,y),np.trapz(z**4*f,y)-3
s2=np.sqrt(2.0)
print("=== PIECE 2: strong-noise (beta=2) deformation — is the ODD part of eta(Y) the 2nd parameter? ===")
print("  beta-family (std,skew): (.904,.26)(.689,.60)(.482,.85)(.320,.85)")
profs={"const":lambda Y:s2,
       "ODD +0.15 tanh(Y)":lambda Y:s2*(1+0.15*np.tanh(Y)),
       "ODD -0.15 tanh(Y)":lambda Y:s2*(1-0.15*np.tanh(Y)),
       "EVEN +0.15(tanh|Y|-.5)":lambda Y:s2*(1+0.15*(np.tanh(abs(Y))-0.5)),
       "EVEN +0.15(|Y|/2-.5)":lambda Y:s2*(1+0.15*(abs(Y)/2-0.5))}
for name,fn in profs.items():
    y,F=solve_fp_eta(fn); m,sd,sk,ek=cums(y,F)
    print(f"  {name:26}: std={sd:.3f} skew={sk:+.3f} exk={ek:+.3f}")
print("  => if ODD profiles move OFF the (std,skew) beta-curve but EVEN stay ON it (~rescaled beta_eff),")
print("     the 2nd (off-family) parameter is the ODD moment of eta(Y). [odd potential sign(Y)Y^2 couples to odd noise]")

print("\n=== PIECE 1: tail RATE picks up eta at the ESCAPE DEPTH (instanton), not the turning? ===")
# two profiles agreeing in bulk (|Y|<2.5) but differing DEEP (|Y|>2.5): should share bulk, differ in tail rate
prof_deep={"const":lambda Y:s2,
           "x1.4 DEEP (|Y|>2.5)":lambda Y: 1.4*s2 if abs(Y)>2.5 else s2}
ss=np.array([3.0,3.5,4.0,4.5,5.0,5.5,6.0])
for name,fn in prof_deep.items():
    y,F=solve_fp_eta(fn,deep=True); mlp=[]
    for s in ss:
        i=np.searchsorted(y,-s); mlp.append(-np.log(max(F[i] if 0<i<len(y) else 1e-300,1e-300)))
    print(f"  {name:22}: -logP at s={list(ss)} = {[round(x,1) for x in mlp]}")
print("  => if 'x1.4 deep' has SMALLER -logP at large s (ratio ~ 1/1.4^2=0.51), tail rate ~ 1/eta(depth)^2 (instanton). ")
