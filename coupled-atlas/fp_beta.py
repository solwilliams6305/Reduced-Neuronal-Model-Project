"""
beta-family consistency test for the PIV sigma-form of the cusp law W_beta.

Solve W_beta via the escape Fokker-Planck PDE for several beta=4/eta^2 (smooth, MC-free).
Natural variable s = (Y* - b)/eta^p : origin b = deterministic edge -2.09 (spectral edge),
scale eta^p with p from the operator balance (predict p=2/5); the test SCANS p and confirms it.
PIV sigma-form:  (sigma'')^2 = A (s sigma'-sigma)^2 + B sigma'^3 + [per-beta: D_b sigma'^2 + E_b sigma'],
A,B SHARED across beta (the rigid PIV structure), D_b,E_b = the monodromy theta(beta) freedom.
Multiple distinct curves (skew varies with beta) must lie on ONE (s sig'-sig)^2+sig'^3 surface
=> breaks the 1-curve confound. Compare joint R^2 to a non-Painleve NULL (skew-normal family).
numpy only.  modes:  solve <beta> ;  test
"""
import numpy as np, sys, os

def solve_fp(eta, Y0=3.0, pmin=-8.0, pmax=7.0, dp=0.025, dt=1.0e-4, tau_max=8.5):
    D=eta*eta/2.0
    pc=np.arange(pmin+dp/2,pmax,dp); pf=np.arange(pmin,pmax+dp/2,dp); pf2=pf*pf
    p0=np.sqrt(np.sign(Y0)*Y0*Y0)
    rho=np.exp(-0.5*((pc-p0)/0.30)**2); rho/=rho.sum()*dp
    nst=int(round(tau_max/dt)); taus=np.empty(nst+1); S=np.empty(nst+1); taus[0]=0; S[0]=1.0
    for k in range(nst):
        Y=Y0-k*dt; V=np.sign(Y)*Y*Y; drift=V-pf2
        rl=np.concatenate(([0.0],rho)); rr=np.concatenate((rho,[0.0]))
        J=np.where(drift>0,rl,rr)*drift - D*(rr-rl)/dp; J[-1]=0.0
        rho=rho-dt*(J[1:]-J[:-1])/dp; np.maximum(rho,0,out=rho)
        taus[k+1]=(k+1)*dt; S[k+1]=rho.sum()*dp
    return (Y0-taus)[::-1], S[::-1]

def cumulants(y,F):
    f=np.gradient(F,y); f=np.clip(f,0,None); f/=np.trapz(f,y)
    m1=np.trapz(y*f,y); c=y-m1; sd=np.sqrt(np.trapz(c*c*f,y)); z=c/sd
    return m1,sd,np.trapz(z**3*f,y),np.trapz(z**4*f,y)-3

def derivs(y,F,b,sc):
    s=(y-b)/sc; m=(F>0.05)&(F<0.95); sg=s[m]; lF=np.log(F[m])
    i=np.argsort(sg); sg,lF=sg[i],lF[i]
    cc=np.polyfit(sg,lF,8); P=np.poly1d(cc)
    se=np.linspace(sg.min()+0.1,sg.max()-0.1,40)
    return se,P.deriv(1)(se),P.deriv(2)(se),P.deriv(3)(se)

def joint_R2(curves, b_det, p_exp):
    # curves: list of (eta, y, F, mean). returns joint R^2 of shared-A,B PIV form.
    X1=[];X2=[];Y=[];blocks=[]; nb=len(curves)
    for i,(eta,y,F,mu) in enumerate(curves):
        b = b_det if b_det is not None else mu
        se,sig,sp,spp=derivs(y,F,b,eta**p_exp)
        X1.append((se*sp-sig)**2); X2.append(sp**3); Y.append(spp**2)
        blocks.append((i,len(se),sp**2,sp))
    X1=np.concatenate(X1);X2=np.concatenate(X2);Y=np.concatenate(Y); n=len(Y)
    cols=[X1,X2]
    off=0; lens=[len(c[2]) for c in blocks]
    starts=np.cumsum([0]+lens)
    for i,(bi,L,d,e) in enumerate(blocks):
        cD=np.zeros(n); cE=np.zeros(n); cD[starts[i]:starts[i+1]]=d; cE[starts[i]:starts[i+1]]=e
        cols+= [cD,cE]
    Xj=np.column_stack(cols); sc=np.maximum(np.abs(Xj).max(0),1e-12); Xn=Xj/sc
    cf,*_=np.linalg.lstsq(Xn,Y,rcond=None); R2=1-np.sum((Y-Xn@cf)**2)/np.sum((Y-Y.mean())**2)
    return R2

# ---- skew-normal null family (smooth, non-Painleve), matched skew range ----
def skewnormal_F(alpha, xg):
    phi=np.exp(-xg**2/2)/np.sqrt(2*np.pi)
    Phi=np.cumsum(phi)*(xg[1]-xg[0]); Phi/=Phi[-1]
    Phia=np.interp(alpha*xg, xg, Phi)
    f=2*phi*Phia; F=np.cumsum(f)*(xg[1]-xg[0]); F/=F[-1]
    return F

if __name__=="__main__":
    mode=sys.argv[1]
    if mode=="solve":
        beta=float(sys.argv[2]); eta=2.0/np.sqrt(beta)
        y,F=solve_fp(eta); m1,sd,sk,ek=cumulants(y,F)
        np.save(f"fpb_{beta:.1f}.npy",np.vstack([y,F]))
        print(f"beta={beta:.1f} eta={eta:.3f}: mean={m1:+.3f} std={sd:.3f} skew={sk:+.3f} exk={ek:+.3f}")
    else:
        betas=[1.0,2.0,4.0,8.0]; curves=[]
        print("loaded W_beta (skew varies => genuine family):")
        for b in betas:
            if not os.path.exists(f"fpb_{b:.1f}.npy"): continue
            y,F=np.load(f"fpb_{b:.1f}.npy"); eta=2.0/np.sqrt(b); m1,sd,sk,ek=cumulants(y,F)
            curves.append((eta,y,F,m1)); print(f"  beta={b}: eta={eta:.2f} mean={m1:+.3f} std={sd:.3f} skew={sk:+.3f}")
        # confirm natural scale: std vs eta
        et=np.array([c[0] for c in curves]); sds=np.array([cumulants(c[1],c[2])[1] for c in curves])
        pscale=np.polyfit(np.log(et),np.log(sds),1)[0]
        print(f"\n  std(W_beta) ~ eta^{pscale:.3f}  (operator-balance prediction: 2/5=0.40)")
        # scan p (scale exponent) and origin, report joint shared-A,B R^2
        print("\n  JOINT shared-A,B PIV sigma-form R^2 (origin=deterministic -2.09):")
        best=(-9,None)
        for p in [0.30,0.40,0.50,0.667]:
            r=joint_R2(curves,-2.09,p); print(f"    scale eta^{p:.3f}: joint R^2 = {r:+.4f}")
            if r>best[0]: best=(r,p)
        rmean=joint_R2(curves,None,best[1])
        print(f"    (origin=per-beta mean, best p): joint R^2 = {rmean:+.4f}")
        # NULL: skew-normal family with matched skews
        xg=np.linspace(-8,8,4000)
        sknull=[]
        for (eta,y,F,mu),al in zip(curves,[0.0,2.0,4.0,8.0]):
            Fn=skewnormal_F(al,xg); sknull.append((eta,xg.copy(),Fn,np.trapz(xg*np.gradient(Fn,xg),xg)))
        rnull=max(joint_R2(sknull,None,p) for p in [0.30,0.40,0.50,0.667])
        print(f"\n  NULL (skew-normal family), best joint R^2 = {rnull:+.4f}")
        print(f"\n  VERDICT: cusp joint R^2 = {best[0]:+.4f} (best p={best[1]}) vs null {rnull:+.4f}")
        print("    cusp>>null with a sensible scale => shared PIV sigma-form (Painleve-IV consistent).")
        print("    cusp ~ null  => no shared sigma-form; W is new beyond PIV (tails (5,3) already place it outside soft-edge family).")
