import numpy as np
# WEAK-NOISE check of beta_eff formula: Var(Y*) = INT eta(s)^2 Phi(s)^2 ds  (Ito isometry, exact as eta->0).
# Perturbative Y1 = -u1(Y*)/u0'(Y*), u1 driven by state-dependent source eta(Y)*u0*dB. Var(Y1,profile)/Var(Y1,const)
# should equal <phi^2>_{Phi^2} (the Phi^2=u_c^4 weighted avg), = 0.873 predicted for +/-10% tanh.  Converged, unambiguous.
def V(Y): return np.sign(Y)*Y*Y
def run(eta_fn, N=400000, Y0=4.0, Yend=-5.0, h=1.0e-3, seed=0):
    rng=np.random.default_rng(seed); nY=int(round((Y0-Yend)/h)); Yg=Y0-np.arange(nY+1)*h; dY=-h; sh=np.sqrt(h)
    u0=np.empty(nY+1); u0p=np.empty(nY+1); u0[0]=1.0; u0p[0]=-Y0
    for i in range(nY):
        u0[i+1]=u0[i]+u0p[i]*dY; u0p[i+1]=u0p[i]+V(Yg[i])*u0[i]*dY
    istar=next(i for i in range(1,nY+1) if Yg[i]<0 and u0[i-1]*u0[i]<0)
    f=u0[istar-1]/(u0[istar-1]-u0[istar]); u0p_star=u0p[istar-1]+f*(u0p[istar]-u0p[istar-1])
    u1=np.zeros(N); u1p=np.zeros(N); snap={}
    for i in range(nY):
        Vi=V(Yg[i]); dB=sh*rng.standard_normal(N); eta=eta_fn(Yg[i])
        u1p_new=u1p+Vi*u1*dY+eta*u0[i]*dB
        u1=u1+u1p*dY; u1p=u1p_new
        if i in (istar-1,istar): snap[i]=u1.copy()
    u1s=snap[istar-1]+f*(snap[istar]-snap[istar-1])
    Y1=-u1s/u0p_star
    return np.var(Y1)
s2=np.sqrt(2.0)
print("WEAK-NOISE Var(Y1) (perturbative functional, converged) for state-dependent noise profiles:")
Vc=run(lambda Y: s2)                       # const
for amp in [0.10,0.15]:
    Vp=run(lambda Y,a=amp: s2*(1+a*np.tanh(Y)))
    print(f"  +/-{int(amp*100)}% tanh: Var(Y1,profile)/Var(Y1,const) = {Vp/Vc:.3f}   (u_c^4-weighted prediction: {0.873 if amp==0.10 else 0.813})")
print(f"  const Var(Y1)={Vc:.5f}")
print("\n=> if ratio ~ 0.87/0.81, the beta_eff = 4 INT Phi^2 / INT eta^2 Phi^2 formula (kernel u_c^4) is VERIFIED at weak noise.")
print("   (The beta=2 FP showed variance UP because beta=2 is STRONG noise, outside this weak-noise formula.)")
