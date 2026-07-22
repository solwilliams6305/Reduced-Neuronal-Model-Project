"""
Program 2, Route 2b -- the NEXT variance coefficient v2 (O(eta^6)) by the DIRECT FULL-FIELD first-node method.

This sidesteps the node-shift recursion wall (the pointwise xi(Y*_0) term needed for Y4/Y5). Instead of expanding
the field in powers of eta, we march the FULL noisy field
        u'' = (V(Y) + eta * dW/dY) u ,   V(Y)=sign(Y)Y^2 ,   u(Y0)=1, u'(Y0)=-Y0,
downward in Y (Euler-Maruyama), and read off Y* = first node with Y<0. The first-node law's cumulants in eta give
        <Y*>  = Y*_0 + m1 eta^2 + m2 eta^4 + ...
        Var   =        v0 eta^2 + v1 eta^4 + v2 eta^6 + ...
        skew  =        s0 eta   + s1 eta^3 + ...
With v0=0.1337, v1=+0.110, m1=+0.212 already PINNED (converged direct-functional, weaknoise_v1.py), we extract
v2 (and m2, s1) by subtracting the known lower orders and fitting the residual across several eta.

The Euler march is itself a regularization at scale h (white noise mollified to correlation length ~h); we
Richardson-extrapolate h->0. This is Route (b) of the kickoff (regularized full field) in its cleanest form --
no explicit mollifier, the grid IS the regulator -- and it is the INDEPENDENT TARGET the boundary-Wick (Route a)
computation of v2 must reproduce.

numpy only, vectorized over realizations.
"""
import numpy as np

def V(Y): return np.sign(Y)*Y*Y

def firstnode_moments(N, eta, h, seed, Y0=4.0, Yend=-6.0):
    """March the full noisy field, return (mean,var,skew, frac_with_node) of the first node Y*<0."""
    rng=np.random.default_rng(seed)
    nY=int(round((Y0-Yend)/h)); dY=-h; sh=np.sqrt(h)
    u=np.ones(N); up=np.full(N,-Y0)
    found=np.zeros(N,dtype=bool); Ystar=np.full(N,np.nan)
    Yprev=Y0; uprev=u.copy()
    for i in range(nY):
        Y=Y0-i*h; Vi=V(Y)
        dB=sh*rng.standard_normal(N)
        new_up=up+Vi*u*dY+eta*u*dB
        u_next=u+up*dY
        up=new_up;
        Ynext=Y0-(i+1)*h
        # detect first node with Ynext<0: sign change of u across the step, not yet found
        cross=(~found)&(Ynext<0.0)&(uprev*u_next<0.0)
        if cross.any():
            # linear interpolation of the zero between (Yprev,uprev) and (Ynext,u_next)
            frac=uprev[cross]/(uprev[cross]-u_next[cross])
            Ystar[cross]=Yprev+frac*(Ynext-Yprev)
            found[cross]=True
        uprev=u_next; Yprev=Ynext
        u=u_next
        if found.all(): break
    Ys=Ystar[found]
    m=Ys.mean(); c=Ys-m
    var=np.mean(c**2); sk=np.mean(c**3)/var**1.5
    return m,var,sk,found.mean(),Ys.size

if __name__=="__main__":
    v0,v1,m1=0.1337,0.110,0.2124
    Ystar0=-2.18762
    print("=== DIRECT FULL-FIELD first-node law: extract v2 (O(eta^6) variance) ===")
    print(f"  known (converged, weaknoise_v1): v0={v0}, v1={v1}, m1={m1}, Y*_0={Ystar0}")
    etas=[0.35,0.45,0.55,0.65,0.75]
    for h in [2e-3,1e-3]:
        print(f"\n  --- h={h:.0e} ---")
        print(f"  {'eta':>5} {'mean':>10} {'var':>12} {'skew':>8} {'v2_est':>10} {'m2_est':>9} {'found%':>7}")
        for eta in etas:
            # average over 4 seeds, N per seed
            Ms=[];Vs=[];Sk=[];Fr=[]
            for sd in range(4):
                m,var,sk,fr,n=firstnode_moments(N=250000,eta=eta,h=h,seed=100*sd+7)
                Ms.append(m);Vs.append(var);Sk.append(sk);Fr.append(fr)
            m=np.mean(Ms);var=np.mean(Vs);sk=np.mean(Sk);fr=np.mean(Fr)
            e2=eta*eta
            # residual variance beyond known orders -> v2
            v2=(var - v0*e2 - v1*e2*e2)/(e2**3)
            m2=((m-Ystar0) - m1*e2)/(e2*e2)
            print(f"  {eta:>5.2f} {m:>10.5f} {var:>12.6f} {sk:>8.4f} {v2:>10.4f} {m2:>9.4f} {100*fr:>6.2f}")
