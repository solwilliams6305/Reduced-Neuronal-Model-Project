import numpy as np
# fp_cusp with STATE-DEPENDENT noise eta(Y) (channel-noise model): D(Y)=eta(Y)^2/2, Y=Y0-tau.
# Riccati dp=(sign(Y)Y^2 - p^2)dtau + eta(Y)dW ; eta depends on the SWEEP Y (not p) => no Ito/Strat drift here.
def solve_fp_eta(eta_fn, Y0=3.0, pmin=-8.0, pmax=7.0, dp=0.02, dt=4e-5, tau_max=7.5):
    pc=np.arange(pmin+dp/2,pmax,dp); pf=np.arange(pmin,pmax+dp/2,dp); pf2=pf*pf
    p0=np.sqrt(max(np.sign(Y0)*Y0*Y0,1e-9))
    rho=np.exp(-0.5*((pc-p0)/0.30)**2); rho/=rho.sum()*dp
    nst=int(round(tau_max/dt)); taus=np.empty(nst+1); S=np.empty(nst+1); taus[0]=0; S[0]=rho.sum()*dp
    for k in range(nst):
        Y=Y0-k*dt; V=np.sign(Y)*Y*Y; D=eta_fn(Y)**2/2.0; drift=V-pf2
        rl=np.concatenate(([0.0],rho)); rr=np.concatenate((rho,[0.0]))
        J=np.where(drift>0,rl,rr)*drift - D*(rr-rl)/dp; J[-1]=0.0
        rho=rho-dt*(J[1:]-J[:-1])/dp; np.maximum(rho,0.0,out=rho)
        taus[k+1]=(k+1)*dt; S[k+1]=rho.sum()*dp
    y=(Y0-taus)[::-1]; F=S[::-1]
    F=np.clip(F,0,1); f=np.gradient(F,y); f=np.clip(f,0,None); f/=np.trapz(f,y)
    m1=np.trapz(y*f,y); c=y-m1; v=np.trapz(c**2*f,y); sd=np.sqrt(v); z=c/sd
    sk=np.trapz(z**3*f,y); ek=np.trapz(z**4*f,y)-3
    return m1,sd,sk,ek

s2=np.sqrt(2.0)  # eta for beta=2
print("REFERENCE constant-eta beta-family (eta^2=4/beta):")
for b in [1,2,4,8]:
    eta=2/np.sqrt(b); m,sd,sk,ek=solve_fp_eta(lambda Y,e=eta:e)
    print(f"  beta={b}: eta={eta:.3f}  mean={m:+.3f} std={sd:.3f} skew={sk:+.3f} exk={ek:+.3f}")

# escape region ~ [-2.5,-1] (typical Y*). Test LOCALIZATION: change eta OUTSIDE vs INSIDE this region.
def band(Y,inside,outside,lo=-2.5,hi=-1.0): return inside if lo<=Y<=hi else outside
print("\nSTATE-DEPENDENT eta(Y) profiles (baseline eta=sqrt2 in escape region [-2.5,-1]):")
cases={
 "A const sqrt2 (=beta2)": lambda Y: s2,
 "B noise x2 OUTSIDE escape": lambda Y: s2 if -2.5<=Y<=-1 else 2*s2,
 "C noise x0.5 OUTSIDE escape": lambda Y: s2 if -2.5<=Y<=-1 else 0.5*s2,
 "D noise HALVED INSIDE escape": lambda Y: 0.5*s2 if -2.5<=Y<=-1 else s2,
 "E gradient eta=sqrt2(1+0.2Y)": lambda Y: s2*(1+0.2*Y),
 "F degenerate: eta small near turning Y~0": lambda Y: s2*np.tanh(abs(Y)/0.8),
}
print(f"  {'case':30} {'mean':>7} {'std':>6} {'skew':>7} {'exk':>7}")
for name,fn in cases.items():
    m,sd,sk,ek=solve_fp_eta(fn)
    print(f"  {name:30} {m:+7.3f} {sd:6.3f} {sk:+7.3f} {ek:+7.3f}")
print("\nDECISIVE: if B (eta changed OUTSIDE escape) ~ A => localization/robustness (only escape-region noise matters).")
print("          D (eta halved INSIDE) should ~ beta=8 (skew/exk of sharper law). Skew sign + / exk sign - = class markers.")
