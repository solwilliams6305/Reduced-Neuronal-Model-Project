"""
Careful re-check of g_crit ~ -0.58 sqrt(eps) in the fingerprint onset of the genuine coupled FHN.
Prior tension: a run with eps=0.008 gave the OPPOSITE onset trend. Hypotheses: (i) eps=0.008 is too slow
(few spikes in fixed T -> noisy), (ii) outer corrections (~eps) drive the sub-Gaussianity, not the inner cusp.
Test: THREE eps in the reliable fast-spiking range (0.012, 0.018, 0.027), FINE g-grid, high M & long T; measure
the spread-amplification onset AND the exk sign-flip onset; check g_onset/sqrt(eps) ~ const (-0.58).
numpy only.
"""
import numpy as np
def sweep(g_vals, eps, c=0.99, sigma=0.02, M=3500, T=3600.0, dt=0.02, warmup=1000.0, seed=0):
    ng=len(g_vals); G=np.asarray(g_vals,float)[:,None]; rng=np.random.default_rng(seed); sdt=np.sqrt(dt)
    v=np.tile([0.2,-0.2],(ng,M,1)).astype(float); w=np.zeros((ng,M,2))
    def step(v,w,noise):
        coup=G[:,:,None]*(v[:,:,::-1]-v)
        return v+(-v**3+3*v-w+coup)*dt+noise, w+eps*(v-c)*dt
    for _ in range(int(warmup/dt)): v,w=step(v,w,0.0)
    on_up=np.zeros((ng,M),bool); wmax=np.full((ng,M),-9.0); gidx=np.tile(np.arange(ng)[:,None],(1,M)); pw=[]; pg=[]
    for _ in range(int(T/dt)):
        noise=sigma*sdt*rng.standard_normal((ng,M,2)); v,w=step(v,w,noise); v1,w1=v[:,:,0],w[:,:,0]
        on_up|=(v1>1.2); wmax=np.where(on_up,np.maximum(wmax,w1),wmax); esc=on_up&(v1<0.5)
        if esc.any(): pw.append(wmax[esc].copy()); pg.append(gidx[esc].copy())
        on_up=np.where(esc,False,on_up); wmax=np.where(esc,-9.0,wmax)
    return (np.concatenate(pw) if pw else np.array([])),(np.concatenate(pg) if pg else np.array([],int))
def phys(x):
    if x.size==0: return x
    m=np.median(x); mad=np.median(np.abs(x-m))+1e-12; return x[np.abs(x-m)<10*mad]
def cross(gv,y,tgt):
    gv=np.array(gv); y=np.array(y); o=np.argsort(gv); gv,y=gv[o],y[o]
    for i in range(len(gv)-1):
        if np.isfinite(y[i]) and np.isfinite(y[i+1]) and (y[i]-tgt)*(y[i+1]-tgt)<=0:
            return gv[i]+(tgt-y[i])*(gv[i+1]-gv[i])/(y[i+1]-y[i])
    return np.nan

gvals=[0.03,0.0,-0.03,-0.05,-0.07,-0.09,-0.11,-0.13,-0.16,-0.20]
epss=[0.012,0.018,0.027]
print("careful sqrt(eps) re-check — single-unit peel-off, fine grid, high stats")
res={}
for eps in epss:
    pw,pg=sweep(gvals,eps,seed=int(1e5*eps))
    gg=[];sd=[];ek=[]
    for j,g in enumerate(gvals):
        x=phys(pw[pg==j])
        if x.size>300:
            d=x-x.mean(); var=np.mean(d**2); gg.append(g); sd.append(np.sqrt(var)); ek.append(np.mean(d**4)/var**2-3)
    sd=np.array(sd); mid=0.5*(sd.min()+sd.max())
    g_spread=cross(gg,sd,mid); g_exk=cross(gg,ek,0.0)
    res[eps]=(g_spread,g_exk)
    print(f"  eps={eps:.3f} (sqrt={np.sqrt(eps):.3f}): spread-onset g={g_spread:+.3f} (g/sqrt(eps)={g_spread/np.sqrt(eps):+.2f}); "
          f"exk-flip g={g_exk:+.3f} (g/sqrt(eps)={g_exk/np.sqrt(eps):+.2f})  [n_events ok, {len(gg)} g-pts]")
print("\n  target: g_crit/sqrt(eps) ~ -0.58 (constant across eps) => sqrt(eps) scaling holds")
rs=[res[e][0]/np.sqrt(e) for e in epss if np.isfinite(res[e][0])]
re=[res[e][1]/np.sqrt(e) for e in epss if np.isfinite(res[e][1])]
if rs: print(f"  spread-onset g/sqrt(eps): {[f'{x:+.2f}' for x in rs]}  (const? spread={np.std(rs):.2f})")
if re: print(f"  exk-flip   g/sqrt(eps): {[f'{x:+.2f}' for x in re]}  (const? spread={np.std(re):.2f})")
print("  => small std across eps AND value near -0.58 => sqrt(eps) CONFIRMED; else the fingerprint-onset does NOT scale as sqrt(eps).")
