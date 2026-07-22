"""
Construct W_beta as the rho->0 noisy parabolic-cylinder connection data, and test
whether it reproduces the cusp fingerprint.

STRUCTURAL DEFINITION. The rho->0 (pure cusp) inner equation is the NOISY
parabolic-cylinder (Weber) equation  u'' = (sign(Y)Y^2 - eta*dW) u.  The recessive
solution at Y->+inf, propagated through the turning, has a first node Y* (= first
explosion of the Cole-Hopf Riccati dp=(sign(Y)Y^2 - p^2)dtau + eta dW).  The
location of that node is set by the (noisy) recessive->oscillatory CONNECTION DATA
of the parabolic-cylinder parametrix.  Define
        W_beta := law of the standardized first-node location,  beta = 4/eta^2.
This is the candidate new T2 cusp edge law.  We compute its standardized cumulants
and tails and compare to the independently-measured cusp fingerprint:
        skew +0.61, exk -0.24, k5 ~ -2.2, k6 ~ -2.7, left-tail exp 5 (e^{-|s|^5/20}).
(TW_2 for contrast: skew +0.224, exk +0.093 -- a DIFFERENT, lighter law.)
numpy only, float32 for speed.
"""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

def first_node(eta, N, Y0=4.0, Yend=-3.5, dt=2e-3, Pthr=10.0, seed=0):
    """first-explosion location Y* of dp=(sign(Y)Y^2 - p^2)dtau + eta dW, start on canard.
    Yend=-3.5 ensures ~99.98% of paths have escaped (no left-tail truncation)."""
    rng = np.random.default_rng(seed)
    nsteps = int((Y0-Yend)/dt)
    p = np.full(N, np.float32(np.sqrt(Y0*Y0)), dtype=np.float32)
    Ystar = np.full(N, np.nan, dtype=np.float32)
    alive = np.ones(N, dtype=bool)
    sq = np.float32(eta*np.sqrt(dt)); dtf = np.float32(dt)
    for i in range(nsteps):
        Y = Y0 - i*dt
        V = np.float32(np.sign(Y)*Y*Y)
        p += (V - p*p)*dtf + sq*rng.standard_normal(N).astype(np.float32)
        np.clip(p, -60.0, 60.0, out=p)
        newly = alive & (p < -Pthr)
        Ystar[newly] = Y
        alive &= ~newly
        if Y < -0.6 and alive.sum() < 0.0003*N:   # ~all exploded
            break
    return Ystar[~np.isnan(Ystar)].astype(np.float64)

def cumulants(x):
    m = x.mean(); d = x - m
    mu = [ (d**k).mean() for k in range(2,7) ]   # mu2..mu6
    mu2,mu3,mu4,mu5,mu6 = mu
    k2=mu2; k3=mu3; k4=mu4-3*mu2**2; k5=mu5-10*mu3*mu2
    k6=mu6-15*mu4*mu2-10*mu3**2+30*mu2**3
    return dict(n=len(x), mean=m, sd=np.sqrt(k2),
                skew=k3/k2**1.5, exk=k4/k2**2, k5n=k5/k2**2.5, k6n=k6/k2**3)

def run_beta(eta, batches, Npb, seed0):
    xs=[]
    for b in range(batches):
        xs.append(first_node(eta, Npb, seed=seed0+b))
    return np.concatenate(xs)

def main():
    fp = dict(skew=0.61, exk=-0.24, k5n=-2.2, k6n=-2.7)   # cusp fingerprint (beta=2)
    print("Constructing W_beta = noisy parabolic-cylinder first-node (connection-data) law\n")
    # main: beta=2 (eta=sqrt2). 3 batches, also used for stability.
    batch=[first_node(np.sqrt(2.0), 120000, seed=100+b) for b in range(3)]
    x2 = np.concatenate(batch)
    c2 = cumulants(x2)
    flip = np.sign(c2['skew']) != np.sign(fp['skew'])   # law defined up to reflection
    if flip: x2=-x2; c2=cumulants(x2)
    print(f"beta=2 (eta=sqrt2): N={c2['n']}   [oriented {'(reflected)' if flip else '(as-is)'}]")
    print(f"  {'':12}{'W_2 (this)':>12}{'cusp fp':>10}{'TW_2':>9}")
    print(f"  {'skew':12}{c2['skew']:>12.3f}{fp['skew']:>10.2f}{0.224:>9.3f}")
    print(f"  {'exc kurt':12}{c2['exk']:>12.3f}{fp['exk']:>10.2f}{0.093:>9.3f}")
    print(f"  {'kappa5/s^5':12}{c2['k5n']:>12.3f}{fp['k5n']:>10.2f}{'-':>9}")
    print(f"  {'kappa6/s^6':12}{c2['k6n']:>12.3f}{fp['k6n']:>10.2f}{'-':>9}")
    sub=[cumulants(b if not flip else -b) for b in batch]
    print(f"  stability (3 batches N=1.3e5): skew {[round(s['skew'],2) for s in sub]}, "
          f"exk {[round(s['exk'],2) for s in sub]}, k5 {[round(s['k5n'],1) for s in sub]}, "
          f"k6 {[round(s['k6n'],1) for s in sub]}")

    print("\nbeta-dependence (skew, exk):")
    for eta,bet in [(2.0,1.0),(1.0,4.0)]:
        cb = cumulants((lambda z: -z if flip else z)(first_node(eta,120000,seed=300+int(10*eta))))
        print(f"  beta={bet:>3.0f} (eta={eta:.3f}):  skew={cb['skew']:+.3f}  exk={cb['exk']:+.3f}")
    print(f"  beta=  2 (eta=1.414):  skew={c2['skew']:+.3f}  exk={c2['exk']:+.3f}")

    # tails of W_2 (standardized)
    s = (x2 - x2.mean())/x2.std()
    print("\ntail exponents of standardized W_2  (FW predict: left 2q+1=5, right 3q/2=3):")
    def tail_exp(s, side):
        ss = -s if side=='left' else s
        ss = ss[ss>0]; ss.sort()
        # survival P(>u); fit ln(-ln P) vs ln u on the tail decade
        u = ss[int(0.90*len(ss)):int(0.999*len(ss))]
        P = 1 - np.arange(len(ss))[int(0.90*len(ss)):int(0.999*len(ss))]/len(ss)
        good = P>0
        A = np.polyfit(np.log(u[good]), np.log(-np.log(P[good])), 1)
        return A[0]
    aL = tail_exp(s,'left'); aR = tail_exp(s,'right')
    print(f"  left-tail exponent  ~ {aL:.2f}  (target 5; super-Gaussian, hard to pin)")
    print(f"  right-tail exponent ~ {aR:.2f}  (target 3)")

    # ---- figure ----
    fig,ax=plt.subplots(1,3,figsize=(15,4.3))
    # (A) density vs normal
    ax[0].hist(s, bins=120, density=True, color='teal', alpha=0.6, label=r'$\mathcal{W}_2$ (constructed)')
    g=np.linspace(-5,5,200); ax[0].plot(g, np.exp(-g*g/2)/np.sqrt(2*np.pi),'k--',label='standard normal')
    ax[0].set_xlim(-5,5); ax[0].set_xlabel('standardized first-node location $s$')
    ax[0].set_ylabel('density'); ax[0].set_title(f'(A) $\\mathcal{{W}}_2$: skew {c2["skew"]:+.2f}, exk {c2["exk"]:+.2f}')
    ax[0].legend(fontsize=8)
    # (B) tails
    for side,col in [('left','crimson'),('right','navy')]:
        ss=(-s if side=='left' else s); ss=ss[ss>0]; ss.sort()
        P=1-np.arange(len(ss))/len(ss); m=(P>1e-5)&(ss>0.5)
        ax[1].plot(np.log(ss[m]), np.log(-np.log(P[m])), '.', color=col, ms=2, label=f'{side} tail')
    for al,lab,col in [(5,'slope 5 (left target)','crimson'),(3,'slope 3 (right target)','navy')]:
        xx=np.linspace(-0.3,1.4,10); ax[1].plot(xx, al*xx + ( -0.5 if al==5 else -1.5),'--',color=col,lw=1,label=lab)
    ax[1].set_xlabel(r'$\ln|s|$'); ax[1].set_ylabel(r'$\ln(-\ln P_{\rm survival})$')
    ax[1].set_title('(B) tail exponents (slope = exponent)'); ax[1].legend(fontsize=7)
    # (C) cumulant comparison
    labels=['skew','exk','k5n','k6n']; W=[c2[k] for k in labels]; F=[fp[k] for k in labels]
    xpos=np.arange(4); w=0.35
    ax[2].bar(xpos-w/2, W, w, color='teal', label=r'$\mathcal{W}_2$ constructed')
    ax[2].bar(xpos+w/2, F, w, color='orange', label='cusp fingerprint')
    ax[2].axhline(0,color='k',lw=0.6); ax[2].set_xticks(xpos)
    ax[2].set_xticklabels([r'skew',r'exk',r'$\kappa_5$',r'$\kappa_6$'])
    ax[2].set_title('(C) fingerprint match'); ax[2].legend(fontsize=8)
    fig.suptitle('Construction of the new cusp edge law W_beta (noisy parabolic-cylinder connection data)',fontsize=12)
    fig.subplots_adjust(left=0.06,right=0.98,bottom=0.12,top=0.86,wspace=0.27)
    plt.savefig("figures/w_beta_construction.png",dpi=110)
    print("\nsaved figures/w_beta_construction.png")

if __name__=="__main__":
    main()
