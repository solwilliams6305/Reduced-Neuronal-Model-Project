"""
Frontier item 1 (per HANDOFF_TO_FABLE5): the exact 2-point structure of the intrinsic Weber node process.

We measure, directly from the canonical additive-noise Weber Riccati (the defining process of W_beta),
    dp = (sign(Y) Y^2 - p^2) d(-Y) + (2/sqrt(beta)) dW,   p ~ +|Y| recessive at Y0,
node = explosion (p<-thr, reset +thr).  Node phase-depth Theta = Y^2/2.

Deliverables:
  (1) mean node spacing vs depth (sanity: should be ~pi and depth-independent  -> stationary pi-lattice in Theta)
  (2) structure factor S(k) via a Hann-windowed estimator; small-k power-law exponent alpha  (THE class arbiter):
          S ~ k^alpha, alpha>1  <=>  Var(N) bounded            = class I  (stronger than GUE)  [current claim]
          S ~ k       (alpha=1) <=>  Var(N) ~ log L            = class II (GUE / Airy edge)
          S ~ k^a<1            <=>  Var(N) ~ L^{1-a}           = class III
     Bragg peak near k = 2pi/pi = 2 confirms the pi-lattice.
  (3) pair correlation g2(r): hard core at 0 (repulsion), oscillations at period pi decaying to 1.
  (4) displacement-variance test: jitter u_n = Theta_n - (n*sbar + off); running Var(u) vs depth
          saturates (class I) or grows ~log (class II).  Independent cross-check of (2).

scipy available here.  Tags [NUMERIC]/[DERIVED]; decisive, non-confounded (exponent, not bulk-cumulant).
"""
import numpy as np

def Vc(Y): return np.sign(Y)*Y*Y

def sample_nodes(beta, N, Y0=3.0, Yend=-28.0, dt=1.5e-3, seed=11, maxnodes=160):
    """Vectorized over N realizations; record up to maxnodes node-depths Y each. thr scaled ~|Y| to keep
    the reset phase-loss ~2/C rad small and depth-uniform (avoids the fixed-thr spacing bias at depth)."""
    eta = 2.0/np.sqrt(beta)
    rng = np.random.default_rng(seed)
    n = int(round((Y0-Yend)/dt)); sq = np.sqrt(dt)
    p = np.full(N, np.sqrt(max(Vc(Y0),1e-9)))
    nodes = np.full((maxnodes, N), np.nan); idx = np.zeros(N, int)
    for i in range(n):
        Y = Y0 - i*dt
        thr = 20.0*abs(Y) + 50.0                      # C=20 => reset phase-loss ~0.1 rad (3% of pi), depth-uniform
        p += (Vc(Y) - p*p)*dt + eta*sq*rng.standard_normal(N)
        np.clip(p, -1e4, 1e4, out=p)
        ex = p < -thr
        live = ex & (idx < maxnodes); cols = np.where(live)[0]
        if cols.size:
            nodes[idx[cols], cols] = Y; idx[cols] += 1
        p[ex] = thr
    return nodes   # (maxnodes, N) node depths Y (descending), nan-padded

def theta_of(nodes):
    return 0.5*nodes*nodes   # Theta = Y^2/2, ascending with node index

# ---------------------------------------------------------------- run
BETA = 2.0
nodes = sample_nodes(BETA, N=6000)
Th = theta_of(nodes)                      # (maxnodes, N)

# (1) spacing vs depth ----------------------------------------------------
sp = Th[1:] - Th[:-1]                      # consecutive phase gaps, per realization
print("=== (1) node spacing vs depth (expect ~pi=3.1416, depth-independent) ===")
for lo,hi in [(0,5),(5,15),(15,30),(30,60),(60,110)]:
    seg = sp[lo:hi]; seg = seg[np.isfinite(seg)]; seg = seg[(seg>0.5)&(seg<8)]
    if seg.size>50:
        print(f"  node index [{lo:3d},{hi:3d}):  mean spacing={seg.mean():.3f}  CV={seg.std()/seg.mean():.3f}  (n={seg.size})")
allsp = sp[np.isfinite(sp)]; allsp=allsp[(allsp>0.5)&(allsp<8)]
sbar = allsp.mean(); rho = 1.0/sbar
print(f"  overall: mean spacing sbar={sbar:.4f}  intensity rho=1/sbar={rho:.4f}")

# (2) structure factor S(k) ----------------------------------------------
# Bulk window in Theta; Hann taper suppresses the k=0 forward-peak leakage. Normalize S->1 at large k.
Tlo, Thi = 8.0, 360.0; L = Thi - Tlo
kgrid = np.linspace(0.02, 3.0, 600)
num = np.zeros_like(kgrid); wsum = 0.0
for r in range(Th.shape[1]):
    t = Th[:, r]; t = t[np.isfinite(t)]; t = t[(t>=Tlo)&(t<=Thi)]
    if t.size < 5: continue
    x = t - Tlo
    w = np.sin(np.pi*x/L)**2                # Hann weight, zero at window ends
    A = w @ np.exp(-1j*np.outer(x, kgrid))  # sum_j w_j e^{-i k x_j}, shape (len k,)
    num += np.abs(A)**2
    wsum += np.sum(w*w)
S = num / wsum                              # -> 1 at large k (Poisson floor); ->0 small k if hyperuniform

# small-k exponent: fit log S ~ alpha log k over a band above the leakage floor (k>~4pi/L) and below Bragg
kfit_lo, kfit_hi = 4*np.pi/L*1.5, 0.9
band = (kgrid>kfit_lo)&(kgrid<kfit_hi)
alpha, c0 = np.polyfit(np.log(kgrid[band]), np.log(S[band]), 1)
# Bragg peak location
bragg_region = (kgrid>1.5)&(kgrid<2.6)
kbragg = kgrid[bragg_region][np.argmax(S[bragg_region])]
print("\n=== (2) structure factor S(k) ===")
print(f"  leakage floor ~4pi/L = {4*np.pi/L:.3f}; fit band k in [{kfit_lo:.3f},{kfit_hi:.3f}]")
print(f"  small-k power:  S(k) ~ k^alpha,  alpha = {alpha:+.2f}")
print(f"    => alpha>1 class I (bounded, >GUE) | alpha~1 class II (GUE log) | alpha<1 class III")
print(f"  Bragg peak at k={kbragg:.3f}  (pi-lattice predicts 2pi/sbar={2*np.pi/sbar:.3f}); S_peak={S[bragg_region].max():.2f}")
for kk in [0.05,0.1,0.2,0.4,0.8,1.5,2.0,2.5]:
    j=np.argmin(abs(kgrid-kk)); print(f"    S({kgrid[j]:.2f})={S[j]:.3f}")

# (3) pair correlation g2(r) ---------------------------------------------
rmax=6*np.pi; nb=240; edges=np.linspace(0,rmax,nb+1); cent=0.5*(edges[1:]+edges[:-1]); dr=edges[1]-edges[0]
hist=np.zeros(nb); Npairs_norm=0.0
for r in range(Th.shape[1]):
    t=Th[:,r]; t=t[np.isfinite(t)]; t=t[(t>=Tlo)&(t<=Thi)]
    if t.size<5: continue
    d=np.abs(t[:,None]-t[None,:]); d=d[np.triu_indices(t.size,1)]
    hist+=np.histogram(d,bins=edges)[0]
    Npairs_norm += t.size
# g2(r) = (counts in shell)/(rho * dr * (#points)) , averaged; factor 2 since we took unordered pairs once
g2 = hist / (rho*dr*Npairs_norm)  * 2.0
print("\n=== (3) pair correlation g2(r) (hard core g2(0+)->0; ->1 at large r) ===")
for rr in [0.3,1.0,np.pi,2*np.pi,3*np.pi,4*np.pi]:
    j=np.argmin(abs(cent-rr)); print(f"    g2(r={cent[j]:.2f})={g2[j]:.3f}")

# (4) displacement variance vs depth -------------------------------------
# jitter u_n = Theta_n - (off + n*sbar); running Var over realizations at each node index
print("\n=== (4) displacement variance vs depth (saturates=class I; grows ~log=class II) ===")
idxn=np.arange(Th.shape[0])
# per-realization offset removed via robust median so we isolate fluctuation, not drift
U=Th - (np.nanmedian(Th - idxn[:,None]*sbar, axis=0)[None,:] + idxn[:,None]*sbar)
for ni in [5,10,20,40,80,120]:
    row=U[ni]; row=row[np.isfinite(row)]
    if row.size>200: print(f"    node {ni:3d} (Theta~{ni*sbar:6.1f}): Var(displacement)={row.var():.4f}  (n={row.size})")

np.savez("node_sk_data.npz", kgrid=kgrid, S=S, cent=cent, g2=g2, sbar=sbar, alpha=alpha)
print("\nsaved node_sk_data.npz")
