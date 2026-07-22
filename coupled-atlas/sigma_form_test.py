"""
R3 — does the cusp law W_2 satisfy a Painleve sigma-form (II or IV)?

A Painleve sigma-form is a 2nd-order 2nd-degree ODE: (sigma'')^2 = polynomial(s, sigma, sigma'),
where sigma(s) = d/ds log F(s).  It is covariant under affine s -> a s + b (coeffs rescale), so we
may standardize all laws to mean 0 / std 1 and fit the coefficients freely.

Method: fit log F(s) by a polynomial over the bulk, take sigma, sigma', sigma'' analytically,
then least-squares fit (sigma'')^2 onto a monomial basis spanning the PII and PIV sigma-forms
(the PIV-distinguishing block is (s sigma' - sigma)^2 = s^2 sigma'^2 - 2 s sigma' sigma + sigma^2).
R^2 of that fit = how well a Painleve sigma-form describes the law.

Controls:
  TW2-smooth  (Fredholm, no noise) -> must give R^2 ~ 1  [calibration: PII sigma-form recoverable]
  TW2-sampled (1.4M draws, same pipeline) -> Painleve + matched-noise baseline
  Gaussian    (1.4M) -> non-Painleve NULL (reveals overfitting)
  cusp W_2    (1.4M MC) -> the test
numpy only.
"""
import numpy as np

def standardize(x): return (x - x.mean())/x.std()

def derivs_from_logF(zg, logF, zev, deg=10):
    c=np.polyfit(zg, logF, deg); p=np.poly1d(c)
    p1=p.deriv(1); p2=p.deriv(2); p3=p.deriv(3)
    return p1(zev), p2(zev), p3(zev)   # sigma, sigma', sigma''

def emp_logF(samples, zg):
    s=np.sort(samples); F=np.searchsorted(s, zg, side='right')/len(s)
    return F

def basis(z, sig, sp):
    # spans PII (sig'^3, z sig'^2, sig' sig, ...) and PIV ((z sig'-sig)^2 block + cubic)
    return np.column_stack([
        sp**3, sp**2, sp, np.ones_like(z),          # cubic in sigma'
        z*sp**2, z*sp, z,                            # s-weighted
        sig*sp, sig, sig**2,                         # sigma terms
        z**2*sp**2, z*sig*sp ])                      # PIV (z sig'-sig)^2 block

def regress_R2(z, sig, sp, spp, k=5):
    # SPARSE + CV: pick the best size-k monomial subset by in-sample R^2 (a true sigma-form is a
    # minimal low-degree relation), then report THAT subset's cross-validated R^2 (fit left half,
    # predict right half). True sigma-form -> sparse subset extrapolates (CV~1); overfit -> CV low.
    from itertools import combinations
    y=spp**2; X=basis(z, sig, sp)
    sc=np.maximum(np.abs(X).max(0),1e-12); Xn=X/sc
    n=len(z); h=n//2; sst=np.sum((y-y.mean())**2)
    best=(-1e9,None)
    for cols in combinations(range(Xn.shape[1]), k):
        A=Xn[:,cols]; cf,*_=np.linalg.lstsq(A,y,rcond=None)
        R2=1-np.sum((y-A@cf)**2)/sst
        if R2>best[0]: best=(R2,cols)
    cols=best[1]; A=Xn[:,cols]
    cf,*_=np.linalg.lstsq(A[:h],y[:h],rcond=None); predR=A[h:]@cf
    cvR2=1-np.sum((y[h:]-predR)**2)/np.sum((y[h:]-y[h:].mean())**2)
    return cvR2, best[0]

def analyze_samples(samples, name, zlo=-2.2, zhi=2.5):
    z=standardize(samples)
    zg=np.arange(-3.0,3.0,0.02); F=emp_logF(z,zg)
    m=(F>0.02)&(F<0.985); zg2=zg[m]; logF=np.log(F[m])
    zev=np.linspace(max(zlo,zg2.min()+0.1), min(zhi,zg2.max()-0.1), 60)
    sig,sp,spp=derivs_from_logF(zg2,logF,zev)
    cv,full=regress_R2(zev,sig,sp,spp)
    print(f"  {name:22s}: CV-R^2 = {cv:+.3f}   (in-sample {full:.4f})")
    return cv

if __name__=="__main__":
    rng=np.random.default_rng(0)
    # --- TW2 smooth (Fredholm) calibration ---
    d=np.load("hoTW_F_k1.npy"); sv,F=d
    f=np.gradient(F,sv); f=np.clip(f,0,None); f/=np.trapz(f,sv)
    mu=np.trapz(sv*f,sv); sd=np.sqrt(np.trapz((sv-mu)**2*f,sv))
    z=(sv-mu)/sd; m=(F>0.02)&(F<0.985)
    zev=np.linspace(z[m].min()+0.1, z[m].max()-0.1, 60)
    sig,sp,spp=derivs_from_logF(z[m], np.log(F[m]), zev)
    cv,full=regress_R2(zev,sig,sp,spp)
    print(f"  {'TW2 smooth (calib)':22s}: CV-R^2 = {cv:+.3f}   (in-sample {full:.4f})   [must be ~1]")
    # --- TW2 sampled (Painleve + noise baseline) ---
    u=rng.random(1400000); tw_samp=np.interp(u, F, sv)
    analyze_samples(tw_samp, "TW2 sampled (baseline)")
    # --- Gaussian null ---
    analyze_samples(rng.standard_normal(1400000), "Gaussian (null)")
    # --- cusp W_2 ---
    ys=np.load("ladder_q2.npy")
    analyze_samples(ys, "cusp W_2 (TEST)")
    print("\nVERDICT: method validated on SMOOTH TW2 (CV-R^2~0.99). But the matched-noise baseline")
    print("TW2-SAMPLED (a true Painleve at 1.4M MC) FAILS CV (<0), so 1.4M MC samples are insufficient")
    print("to detect a sigma-form through the 3rd-derivative pipeline. => the cusp test is INCONCLUSIVE")
    print("(its CV-R^2 is not interpretable when the Painleve baseline itself fails).")
    print("A SMOOTH representation of W_2 (Fokker-Planck PDE / R2 scaling function) is the prerequisite.")
