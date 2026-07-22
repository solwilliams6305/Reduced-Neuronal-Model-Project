"""Pre-registration of v7/v8 predictions + discrimination + nonparametric checks.
Pure fitting on the exact ladder v0..v6 = [0.134,0.111,0.104,-0.030,-0.451,-1.19,-1.90].
Model convention (matches _borel_analysis.darboux_model): v_n ~ 2C zeta^{-(n+1)} Gamma(n+1+alpha) cos((n+1)theta - phi).
So order index k = n+1 runs 1..7 for v0..v6.
"""
import numpy as np
from math import factorial, pi
from scipy.special import gamma as Gamma
from scipy.optimize import least_squares

V = np.array([0.134, 0.111, 0.104, -0.030, -0.451, -1.19, -1.90])   # v0..v6
NS = np.arange(len(V), dtype=float)

# ---------- flexible darboux fitters ----------
def pair_model(p, n):
    C, zeta, theta, phi, alpha = p
    return 2*C*zeta**(-(n+1.0))*Gamma(n+1.0+alpha)*np.cos((n+1.0)*theta - phi)

def fit_pair(v, theta_fix=None, alpha_fix=None, nmin=0):
    ns = np.arange(len(v), dtype=float)[nmin:]
    y = np.asarray(v, float)[nmin:]
    scale = np.abs(y) + 1e-3
    def unpack(p):
        # order of free params: C, zeta, [theta], phi, [alpha]
        i = 0
        C = p[i]; i+=1
        zeta = p[i]; i+=1
        if theta_fix is None: theta = p[i]; i+=1
        else: theta = np.radians(theta_fix)
        phi = p[i]; i+=1
        if alpha_fix is None: alpha = p[i]; i+=1
        else: alpha = alpha_fix
        return np.array([C,zeta,theta,phi,alpha])
    def resid(p):
        return (pair_model(unpack(p), ns) - y)/scale
    best=None
    for th0 in np.radians([40,45,48,50,55,60]):
        for z0 in [1.2,1.5,1.9,2.3]:
            p0=[0.1, z0]
            lo=[1e-6,0.3]; hi=[50.0,6.0]
            if theta_fix is None:
                p0.append(th0); lo.append(np.radians(5)); hi.append(np.radians(120))
            p0.append(0.5); lo.append(-2*pi); hi.append(2*pi)
            if alpha_fix is None:
                p0.append(0.0); lo.append(-1.5); hi.append(1.5)
            try:
                r=least_squares(resid,p0,bounds=(lo,hi),xtol=1e-14,ftol=1e-14,max_nfev=5000)
            except ValueError:
                continue
            if best is None or r.cost<best.cost:
                best=r; bp=unpack(r.x)
    res=np.max(np.abs(best.fun))
    return bp, res

def predict_pair(p, ks):
    # ks are n-values (v_n). returns model at those n
    return pair_model(p, np.asarray(ks,float))

# ---------- two-singularity (pair + real), theta pinned ----------
def twosing_model(p, n):
    C,zeta,theta,phi,alpha,Cr,zr,ar = p
    pair=2*C*zeta**(-(n+1.0))*Gamma(n+1.0+alpha)*np.cos((n+1.0)*theta-phi)
    real=Cr*zr**(-(n+1.0))*Gamma(n+1.0+ar)
    return pair+real

def fit_twosing(v, theta_fix=45.0, alpha_fix=0.0, ar_fix=0.0, nmin=0):
    ns=np.arange(len(v),dtype=float)[nmin:]
    y=np.asarray(v,float)[nmin:]
    scale=np.abs(y)+1e-3
    def unpack(p):
        i=0
        C=p[i];i+=1; zeta=p[i];i+=1
        theta=np.radians(theta_fix)
        phi=p[i];i+=1
        alpha=alpha_fix
        Cr=p[i];i+=1; zr=p[i];i+=1
        ar=ar_fix
        return np.array([C,zeta,theta,phi,alpha,Cr,zr,ar])
    def resid(p):
        return (twosing_model(unpack(p),ns)-y)/scale
    best=None
    for z0 in [1.5,1.9,2.3]:
        for sr in [+1,-1]:
            for zr0 in [1.5,1.9,2.5]:
                p0=[0.1,z0,0.5,sr*0.05,zr0]
                lo=[1e-6,0.3,-2*pi,-10,0.3]
                hi=[50,6,2*pi,10,6]
                try:
                    r=least_squares(resid,p0,bounds=(lo,hi),xtol=1e-14,ftol=1e-14,max_nfev=5000)
                except ValueError:
                    continue
                if best is None or r.cost<best.cost:
                    best=r; bp=unpack(r.x)
    res=np.max(np.abs(best.fun))
    return bp, res

# ================= TASK 1: pre-register v7,v8 with uncertainty =================
print("="*80)
print("TASK 1: v7/v8 predictions per hypothesis (with v6+/-0.15, v5+/-0.01 refit spread)")
print("="*80)

def perturbations():
    out=[]
    for dv6 in (0.0, +0.15, -0.15):
        for dv5 in (0.0, +0.01, -0.01):
            w=V.copy(); w[6]+=dv6; w[5]+=dv5
            out.append(((dv6,dv5), w))
    return out

def summarize(name, fitfn):
    v7s=[]; v8s=[]; thetas=[]; resids=[]; zetas=[]
    for (d,w) in perturbations():
        p,res = fitfn(w)
        pr = pair_model(p, np.array([7.0,8.0])) if len(p)==5 else twosing_model(p, np.array([7.0,8.0]))
        v7s.append(pr[0]); v8s.append(pr[1]); resids.append(res)
        thetas.append(np.degrees(p[2])); zetas.append(p[1])
    v7s=np.array(v7s); v8s=np.array(v8s)
    # central = nominal (first entry, d=(0,0))
    print(f"\n[{name}]")
    print(f"  central: theta={thetas[0]:.1f}deg zeta={zetas[0]:.3f} resid={resids[0]:.3f}")
    print(f"  v7 = {v7s[0]:+.2f}  spread over perturb [{v7s.min():+.2f},{v7s.max():+.2f}] "
          f"(+/-{0.5*(v7s.max()-v7s.min()):.2f})")
    print(f"  v8 = {v8s[0]:+.2f}  spread over perturb [{v8s.min():+.2f},{v8s.max():+.2f}] "
          f"(+/-{0.5*(v8s.max()-v8s.min()):.2f})")
    return v7s[0], 0.5*(v7s.max()-v7s.min()), v8s[0], 0.5*(v8s.max()-v8s.min())

# H1: 45 pinned + real pole (two-sing)
h1 = summarize("H1: theta=45 pinned + real instanton (two-sing, alpha=ar=0)",
               lambda w: fit_twosing(w, theta_fix=45.0))
# H2: free single pair (~48-50)
h2 = summarize("H2: free single pair (alpha free)",
               lambda w: fit_pair(w, theta_fix=None, alpha_fix=None))
# H3: 55 pinned single pair
h3 = summarize("H3: theta=55 pinned single pair (alpha free)",
               lambda w: fit_pair(w, theta_fix=55.0, alpha_fix=None))
# extra: 50 pinned and 48 pinned single pair for discrimination
h50 = summarize("H2b: theta=50 pinned single pair (alpha free)",
               lambda w: fit_pair(w, theta_fix=50.0, alpha_fix=None))
h45 = summarize("H1b: theta=45 pinned single pair (alpha free)",
               lambda w: fit_pair(w, theta_fix=45.0, alpha_fix=None))

# ================= TASK 2: discrimination power of v7 =================
print("\n"+"="*80)
print("TASK 2: v7 precision to discriminate 45 vs 50 at 3:1 residual odds")
print("="*80)
# predictions from single-pair pinned fits
v7_45, e45, _, _ = h45
v7_50, e50, _, _ = h50
gap = abs(v7_50 - v7_45)
v7_48 = h2[0]  # free fit lands ~48.6
print(f"  single-pair pinned predictions: v7(45)={v7_45:+.2f}, v7(48free)={v7_48:+.2f}, v7(50)={v7_50:+.2f}, gap(45-50)={gap:.2f}, gap(45-48)={abs(v7_48-v7_45):.2f}")
# expected grid error on v7: v6 was +/-0.15; envelope ratio v7/v6 ~ |v6|/|v5|*... estimate
env_ratio = abs(V[6])/abs(V[5])   # ~1.6 but that's suppressed(node). Use factorial-ish growth
# per notes: assume proportional to v6's, scaled by envelope ratio ~7x
for sigma7 in [0.3, 1.0, 2.0, 3.0]:
    # residual-odds: model as chi over the single v7 point. odds ~ exp(-(dev^2)/(2 sigma^2)) ratio
    # If truth=45: resid_45 ~ 0, resid_50 = gap. 3:1 odds in likelihood requires gap/sigma >= sqrt(2 ln 3)
    thresh = np.sqrt(2*np.log(3))
    decisive = gap/sigma7 >= thresh
    print(f"  sigma7={sigma7:.1f}: gap/sigma={gap/sigma7:.2f} vs 3:1-threshold {thresh:.2f} -> "
          f"{'DECISIVE' if decisive else 'not decisive'}")
print(f"  => 45-vs-50: need sigma7 <= gap/{np.sqrt(2*np.log(3)):.2f} = {gap/np.sqrt(2*np.log(3)):.2f} for 3:1 on a single v7 datum")
gap48 = abs(v7_48-v7_45)
print(f"  => 45-vs-48: gap={gap48:.2f}, need sigma7 <= {gap48/np.sqrt(2*np.log(3)):.2f}; "
      f"at sigma7=1: gap/sig={gap48/1.0:.2f} ({'decisive' if gap48/1.0>=1.48 else 'NOT'}); "
      f"at sigma7=3: gap/sig={gap48/3.0:.2f} ({'decisive' if gap48/3.0>=1.48 else 'NOT'})")

# ================= TASK 3a: sign-run / node region map =================
print("\n"+"="*80)
print("TASK 3a: nonparametric (theta,phi) region for sign pattern +,+,+,-,-,-,- with v3 & v6 near-node")
print("="*80)
target_signs = np.array([+1,+1,+1,-1,-1,-1,-1])
ks = np.arange(1,8)   # k=n+1 = 1..7
def region_scan():
    ths = np.arange(30,75.01,0.5)
    phs = np.arange(-180,180.01,1.0)
    allowed=[]      # sign pattern only
    allowed_node=[] # sign pattern + v3 node + v6 node
    for th in ths:
        thr=np.radians(th)
        for ph in phs:
            phr=np.radians(ph)
            c = np.cos(ks*thr - phr)
            if np.all(np.sign(c)==target_signs):
                allowed.append((th,ph))
                # near-node conditions: |cos(4theta-phi)| small (v3, k=4) AND
                # |cos(7theta-phi)| suppressed vs neighbor (v6, k=7)
                cnode3 = abs(c[3])   # k=4 -> index 3
                cnode6 = abs(c[6])   # k=7 -> index 6
                # require v3 node: |c(k=4)| < 0.35 ; v6 suppressed: |c(k=7)| < |c(k=6)|
                if cnode3 < 0.35 and cnode6 < abs(c[5]):
                    allowed_node.append((th,ph))
    return np.array(allowed), np.array(allowed_node)

allowed, allowed_node = region_scan()
if len(allowed):
    print(f"  sign-pattern-only allowed theta range: [{allowed[:,0].min():.1f}, {allowed[:,0].max():.1f}]deg "
          f"({len(allowed)} grid cells)")
else:
    print("  no theta admits the sign pattern (check)")
if len(allowed_node):
    ths=allowed_node[:,0]
    print(f"  sign+node(v3<0.35, v6<neighbor) allowed theta: [{ths.min():.1f},{ths.max():.1f}]deg "
          f"({len(allowed_node)} cells)")
    for probe in [45,48,50,55,60]:
        hit = np.any(np.abs(allowed_node[:,0]-probe)<0.6)
        print(f"    theta={probe}deg in node-region? {'YES' if hit else 'no'}")
else:
    print("  node-region empty")

# --- AMPLITUDE-FREE envelope-consistency test (alpha=0): b_n/cos((n+1)th-phi) must be clean geometric ---
# g_n = b_n / cos((n+1)theta - phi) = 2C zeta^{-(n+1)} : all same sign, log|g_n| linear in n.
print("  --- amplitude-free envelope-consistency (b_n/cos must be single-sign geometric) ---")
bb = np.array([V[n]/factorial(n) for n in range(len(V))])
def envelope_scatter(th_deg):
    thr=np.radians(th_deg)
    best=None
    for ph in np.arange(-180,180.01,0.5):
        phr=np.radians(ph)
        c=np.cos(ks*thr-phr)
        if not np.all(np.sign(c)==target_signs):
            continue
        # exclude near-node terms (|c|<0.25) from the geometric fit to avoid blow-up
        mask=np.abs(c)>=0.25
        if mask.sum()<4:
            continue
        g=bb[mask]/c[mask]
        if not np.all(g>0):       # 2C zeta^{-(n+1)} must be single-signed
            continue
        x=ks[mask].astype(float); lg=np.log(g)
        A=np.vstack([np.ones_like(x),x]).T
        sol,_,_,_=np.linalg.lstsq(A,lg,rcond=None)
        resid=np.max(np.abs(A@sol-lg))     # max log-residual of geometric fit
        zeta_env=np.exp(-sol[1])           # slope = -log zeta
        if best is None or resid<best[0]:
            best=(resid,ph,zeta_env,mask.sum())
    return best
print("    theta  best-phi  logResid(geom)  zeta_env  #pts   verdict(<0.15 clean)")
for th in [40,42,44,45,46,48,50,52,54,55,57,60]:
    r=envelope_scatter(th)
    if r is None:
        print(f"    {th:5.1f}   ----     (no phi: sign or single-sign envelope fails)  EXCLUDED")
    else:
        resid,ph,ze,npts=r
        verdict='clean' if resid<0.15 else ('marginal' if resid<0.4 else 'poor')
        print(f"    {th:5.1f}   {ph:6.1f}    {resid:8.3f}       {ze:6.3f}   {npts:d}    {verdict}")

# tighter node on v3 only (the strong constraint), report phi at 45/48/55
print("  --- phi window per theta (sign+v3-node) ---")
for probe in [45,48,50,55,60]:
    thr=np.radians(probe)
    phs=np.arange(-180,180.01,0.5); good=[]
    for ph in phs:
        c=np.cos(ks*thr-np.radians(ph))
        if np.all(np.sign(c)==target_signs) and abs(c[3])<0.35:
            good.append(ph)
    if good:
        print(f"    theta={probe}: phi in [{min(good):.1f},{max(good):.1f}]deg, min|cos(4th-phi)| feasible")
    else:
        print(f"    theta={probe}: NO phi gives sign+v3node")

# ================= TASK 3b: 3-term recursion phase extraction =================
print("\n"+"="*80)
print("TASK 3b: Hadamard/3-term recursion phase extraction on b_n = v_n/n!")
print("="*80)
b = np.array([V[n]/factorial(n) for n in range(len(V))])
print(f"  b_n = v_n/n!: {['%+.4e'%x for x in b]}")
# For a single complex pair b_n ~ A r^n cos(n*theta+psi), satisfies b_{n+1} = 2 (cos theta / rho) b_n - (1/rho^2) b_{n-1}
# i.e. 3-term recursion b_{n+1} + p b_n + q b_{n-1} = 0 with q=1/rho^2>0, p=-2 cos theta/rho.
# Solve least squares for (p,q) over interior n, extract rho=1/sqrt(q), theta=acos(-p/(2 sqrt(q))).
def three_term(bv):
    # rows: b_{n+1} = -p b_n - q b_{n-1}
    n0=1
    rows=[]; rhs=[]
    for n in range(1,len(bv)-1):
        rows.append([bv[n], bv[n-1]])   # [-p,-q] multiply
        rhs.append(-bv[n+1])
    A=np.array(rows); rr=np.array(rhs)
    sol,res,rank,sv=np.linalg.lstsq(A,rr,rcond=None)
    p,q=sol   # A@[p,q] = rhs? Let's define: b_{n+1} = -p b_n - q b_{n-1} => p*b_n + q*b_{n-1} = -b_{n+1}
    cond = sv[0]/sv[-1] if sv[-1]>0 else np.inf
    return p,q,cond,sv
p,q,cond,sv = three_term(b)
print(f"  3-term LS over all interior points: p={p:+.4f} q={q:+.4f} cond={cond:.2e} sv={sv}")
if q>0:
    root_mag=np.sqrt(q)          # |char root| = growth rate of b_n
    zeta_borel=1.0/root_mag      # Borel radius |t_sing| = 1/growth-rate
    arg=-p/(2*np.sqrt(q))
    if abs(arg)<=1:
        theta=np.degrees(np.arccos(arg))
        print(f"    => b_n growth-rate |root|={root_mag:.3f}, Borel radius |zeta|={zeta_borel:.3f}, theta={theta:.1f}deg")
    else:
        print(f"    => q>0 but |arg|={abs(arg):.2f}>1: no real theta (over/under-damped) -> real-axis dominance")
else:
    print(f"    => q<=0: recursion not a clean complex pair (q={q:.4f})")
# condition-number comparison: 5-coeff (v0..v4) vs 7-coeff (v0..v6)
b5=b[:5]; p5,q5,c5,sv5=three_term(b5)
print(f"  5-coeff (v0..v4) 3-term: p={p5:+.4f} q={q5:+.4f} cond={c5:.2e}")
print(f"  7-coeff (v0..v6) 3-term: cond={cond:.2e}  (lower cond = better conditioned)")

# also do the 2-term ratio (Hadamard) diagnostic
print("  --- Hadamard ratio b_{n+1}/b_n ---")
with np.errstate(divide='ignore', invalid='ignore'):
    rat=b[1:]/b[:-1]
print(f"    ratios: {['%+.3f'%x for x in rat]}  (sign flips mark oscillation; |ratio|~|zeta| envelope)")
