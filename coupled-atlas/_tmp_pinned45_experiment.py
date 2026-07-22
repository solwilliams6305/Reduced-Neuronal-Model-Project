import sys, numpy as np
from math import pi
from scipy.special import gamma as Gamma
from scipy.optimize import least_squares

# The real sharpened 7-coefficient ladder v0..v6
VREAL = np.array([0.134, 0.111, 0.104, -0.030, -0.451, -1.19, -1.90])

# ---------- models ----------
def pair_term(C, zeta, theta, phi, alpha, n):
    return 2*C*zeta**(-(n+1.0))*Gamma(n+1.0+alpha)*np.cos((n+1.0)*theta-phi)

def real_term(Cr, zr, ar, n):
    return Cr*zr**(-(n+1.0))*Gamma(n+1.0+ar)

# ---------- generic robust fitters ----------
def fit_single_pair(v, multistart_theta=None, lean=False):
    """Free-theta single complex-pair Darboux fit (5 params). Returns (params, maxrelresid)."""
    ns = np.arange(len(v), dtype=float); y = np.asarray(v, float)
    scale = np.abs(y)+1e-3
    def resid(p): return (pair_term(p[0],p[1],p[2],p[3],p[4],ns)-y)/scale
    lo=[1e-6,0.3,np.radians(5),-2*pi,-3.0]; hi=[50.0,6.0,np.radians(120),2*pi,3.0]
    if multistart_theta is not None: thetas=multistart_theta
    elif lean: thetas=np.radians(np.arange(25,100,5.0))
    else: thetas=np.radians(np.arange(20,110,2.0))
    zs = [1.6,2.2,3.0] if lean else [1.2,1.6,2.0,2.4,3.0]
    as_ = [0.0] if lean else [-0.5,0.0,0.5]
    best=None
    for th0 in thetas:
        for z0 in zs:
            for a0 in as_:
                p0=[0.1,z0,th0,0.5,a0]
                try: r=least_squares(resid,p0,bounds=(lo,hi),xtol=1e-14,ftol=1e-14,max_nfev=2000)
                except ValueError: continue
                if best is None or r.cost<best.cost: best=r
    p=best.x; return p, float(np.max(np.abs(best.fun)))

def fit_pinned_pair_plus_real(v, theta_deg=45.0, ar_mode='eq_a'):
    """Pair theta PINNED + real instanton. ar_mode: 'eq_a' (ar=alpha), 'zero' (ar=0), 'free'.
    Params fitted: C,zeta,phi,alpha,Cr,zr[,ar]. Returns (fullparams_dict, maxrelresid)."""
    ns=np.arange(len(v),dtype=float); y=np.asarray(v,float); scale=np.abs(y)+1e-3
    th=np.radians(theta_deg)
    def unpack(p):
        C,zeta,phi,alpha,Cr,zr = p[:6]
        if ar_mode=='eq_a': ar=alpha
        elif ar_mode=='zero': ar=0.0
        else: ar=p[6]
        return C,zeta,phi,alpha,Cr,zr,ar
    def resid(p):
        C,zeta,phi,alpha,Cr,zr,ar=unpack(p)
        m=pair_term(C,zeta,th,phi,alpha,ns)+real_term(Cr,zr,ar,ns)
        return (m-y)/scale
    lo=[1e-6,0.3,-2*pi,-3.0,-50,0.3]; hi=[50,6.0,2*pi,3.0,50,6.0]
    if ar_mode=='free': lo=lo+[-3.0]; hi=hi+[3.0]
    best=None
    for z0 in [1.4,1.8,2.2,2.6]:
        for zr0 in [1.5,1.9,2.4]:
            for crsign in [+1,-1]:
                for a0 in [-0.5,0.0,0.5]:
                    p0=[0.1,z0,0.5,a0,crsign*0.1,zr0]
                    if ar_mode=='free': p0=p0+[0.0]
                    try: r=least_squares(resid,p0,bounds=(lo,hi),xtol=1e-15,ftol=1e-15,max_nfev=6000)
                    except ValueError: continue
                    if best is None or r.cost<best.cost: best=r
    C,zeta,phi,alpha,Cr,zr,ar=unpack(best.x)
    d=dict(C=C,zeta=zeta,theta_deg=theta_deg,phi=phi,alpha=alpha,Cr=Cr,zr=zr,ar=ar)
    return d, float(np.max(np.abs(best.fun))), best

def predict_pinned(d, N):
    ns=np.arange(N,dtype=float)
    th=np.radians(d['theta_deg'])
    return pair_term(d['C'],d['zeta'],th,d['phi'],d['alpha'],ns)+real_term(d['Cr'],d['zr'],d['ar'],ns)

def predict_single(p,N):
    ns=np.arange(N,dtype=float); return pair_term(p[0],p[1],p[2],p[3],p[4],ns)

# =========================================================================
print("="*80); print("TASK 1: fit REAL ladder with pinned-45 pair + real vs unpinned single-pair")
print("="*80)
print("real ladder v0..v6:", [f"{x:+.3f}" for x in VREAL])

# reference: unpinned single pair (free theta), and theta-scan at 45 & 50
psp, rsp = fit_single_pair(VREAL)
print(f"\n[unpinned single-pair, theta FREE]  theta={np.degrees(psp[2]):.1f}  |zeta|={psp[1]:.3f}  "
      f"C={psp[0]:.4f}  phi={np.degrees(psp[3]):.1f}  alpha={psp[4]:+.3f}  maxrelresid={rsp:.4f}")
# pinned single-pair at 45 and 50 (theta multistart pinned)
for thd in (45.0,50.0):
    p,r=fit_single_pair(VREAL, multistart_theta=[np.radians(thd)])
    # force theta stay: refit with theta bound tight
    print(f"[single-pair seeded {thd:g}]  ended theta={np.degrees(p[2]):.1f}  maxrelresid={r:.4f}")

print("\n[pinned-45 pair + real instanton]:")
for mode in ('eq_a','zero','free'):
    d,r,_=fit_pinned_pair_plus_real(VREAL, 45.0, mode)
    nfit = 6 if mode!='free' else 7
    print(f"  ar_mode={mode:5s} ({nfit} params): pair |zeta|={d['zeta']:.3f} C={d['C']:.4f} "
          f"phi={np.degrees(d['phi']):.1f} alpha={d['alpha']:+.3f} | real zr={d['zr']:.3f} "
          f"Cr={d['Cr']:+.4f} ar={d['ar']:+.3f}  maxrelresid={r:.4f}")
    if mode=='eq_a':
        pv=predict_pinned(d,9)
        print(f"    -> predicts v7={pv[7]:+.3f}, v8={pv[8]:+.3f}")
        d_keep=d

# single-pair v7/v8 predictions
pv_sp=predict_single(psp,9)
print(f"\n[unpinned single-pair] predicts v7={pv_sp[7]:+.3f}, v8={pv_sp[8]:+.3f} (theta={np.degrees(psp[2]):.1f})")

# =========================================================================
# Choose a "true" complex pair consistent with the real ladder for synthetic generation.
# Use the unpinned single-pair fit's C, phi, alpha as the pair template (physically motivated);
# we then vary theta_true, pair radius zeta, and the real-pole (Cr, zr).
def make_synth(theta_true_deg, zeta, Cr, zr, C=None, phi=None, alpha=None, ar=None, N=7):
    C   = psp[0] if C   is None else C
    phi = psp[3] if phi is None else phi
    alpha = psp[4] if alpha is None else alpha
    ar  = alpha if ar is None else ar
    ns=np.arange(N,dtype=float); th=np.radians(theta_true_deg)
    return pair_term(C,zeta,th,phi,alpha,ns)+real_term(Cr,zr,ar,ns)

print("\n"+"="*80)
print("TASK 2: synthetic bias map -- theta_true=45, fit with FREE-theta single pair")
print("="*80)
Cr_step1 = d_keep['Cr']  # from pinned-45 eq_a fit
print(f"(step-1 pinned-45 real-pole strength Cr={Cr_step1:+.4f}; grid spans +-factor 3 around it)")
zeta_grid = [1.8, 2.0, 2.3, 2.6]
Cr_grid = sorted(set([round(x,4) for x in
            [Cr_step1/3, Cr_step1, Cr_step1*3, -Cr_step1, 0.05, 0.15, 0.3, -0.1, -0.3]]))
zr_grid = [1.5, 1.9, 2.4]
print(f"pair radius zeta in {zeta_grid}; real Cr in {Cr_grid}; zr in {zr_grid}")
print(f"\n{'zeta':>5} {'zr':>5} {'Cr':>8} | {'fit_theta':>9} {'fit_zeta':>8} {'resid':>7}  flag")
hits=[]
for zeta in zeta_grid:
    for zr in zr_grid:
        for Cr in Cr_grid:
            vs = make_synth(45.0, zeta, Cr, zr)
            p,r = fit_single_pair(vs, lean=True)
            ft=np.degrees(p[2])
            flag=""
            if 47.5<=ft<=51.0 and 0.09<=r<=0.16: flag="<== masquerade (theta~48-50 & resid~0.12)"
            if flag: hits.append((zeta,zr,Cr,ft,r))
            print(f"{zeta:5.2f} {zr:5.2f} {Cr:+8.4f} | {ft:9.1f} {p[1]:8.3f} {r:7.4f}  {flag}")
print(f"\n# masquerade hits (fit_theta in 47.5-51 AND resid 0.09-0.16): {len(hits)}")
for h in hits: print(f"   zeta={h[0]} zr={h[1]} Cr={h[2]:+.4f} -> fit_theta={h[3]:.1f} resid={h[4]:.4f}")

print("\n"+"="*80)
print("TASK 3: reverse check -- theta_true=50 + real pole, fit with FREE-theta single pair")
print("="*80)
print(f"\n{'zeta':>5} {'zr':>5} {'Cr':>8} | {'fit_theta':>9} {'fit_zeta':>8} {'resid':>7}")
for zeta in [2.0, 2.3]:
    for zr in zr_grid:
        for Cr in [Cr_step1, Cr_step1*3, -Cr_step1, 0.15, -0.15]:
            vs = make_synth(50.0, zeta, round(Cr,4), zr)
            p,r = fit_single_pair(vs, lean=True)
            print(f"{zeta:5.2f} {zr:5.2f} {round(Cr,4):+8.4f} | {np.degrees(p[2]):9.1f} {p[1]:8.3f} {r:7.4f}")
# clean no-real-pole controls
print("\n[controls: pure pair, NO real pole -- single-pair should recover theta_true exactly]")
for tt in (45.0, 48.0, 50.0):
    for zeta in (2.0, 2.3):
        vs=make_synth(tt, zeta, 0.0, 1.9)
        p,r=fit_single_pair(vs, lean=True)
        print(f"  theta_true={tt:g} zeta={zeta}: fit_theta={np.degrees(p[2]):.2f} fit_zeta={p[1]:.3f} resid={r:.5f}")
