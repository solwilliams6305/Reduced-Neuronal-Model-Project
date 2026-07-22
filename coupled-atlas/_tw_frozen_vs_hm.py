"""
_tw_frozen_vs_hm.py -- the exact relation between the FROZEN escape Stokes data
(S0=1/pi, index gamma'=0) and the FULL Hastings-McLeod backbone (C=sqrt(2/(3pi^3)),
index gamma=-1/2).  Purpose: characterize precisely what the "x-extension" does, using
the two series already in hand, before claiming to "derive" the sqrt(2/(3pi)) factor.

Dictionary (matched two-instanton exponentials): g = beta a^{3/2} = sqrt2 (-s)^{3/2},
so (-s)^{-3} = 2/g^2.  The HM series q ~ sqrt(-s/2) sum_k a_k (-s)^{-3k}, re-expanded in g,
lands at EVEN orders g^{-2k} with coefficient a_k 2^k.  Compare to the frozen loop series
F(g)=sum_m d_m g^{-m} at m=2k.

Growth (both established):
   frozen : d_{2k} ~ -(1/pi) Gamma(2k) (9/4)^k          (gamma'=0, base 3/2 => (3/2)^{2k}=(9/4)^k)
   HM     : a_k 2^k ~ -C Gamma(2k-1/2) (9/4)^k           (gamma=-1/2, base 9/8 => (9/8)^k 2^k=(9/4)^k)
Same action (9/4)^k, but the INDEX differs by 1/2 and the CONSTANT by C/(1/pi)=pi C=sqrt(2/(3pi)).
So the ratio  r_k := d_{2k}/(a_k 2^k) ~ (1/(pi C)) Gamma(2k)/Gamma(2k-1/2) ~ sqrt(3pi/2)*sqrt(2k):
   r_k / sqrt(2k) -> 1/(pi C) = sqrt(3pi/2) = 1/sqrt(2/(3pi)).
This makes the "x-extension" precise: it shifts the resurgent index by -1/2 AND rescales the
constant by sqrt(2/(3pi)); the two come bundled (a determinant/collective-mode signature),
NOT a stand-alone multiplicative factor.
"""
import mpmath as mp, json, os
mp.mp.dps = 50
HERE = os.path.dirname(os.path.abspath(__file__))
C_HM = mp.sqrt(mp.mpf(2)/(3*mp.pi**3))


def qt_coeffs(J):
    out = [mp.mpf(1)]
    for j in range(1, J + 1):
        num = mp.mpf(1)
        for i in range(6):
            num *= (6*j - i)
        den = mp.mpf(576)*(3*j)*(3*j-1)*(3*j-2)*(2*j)*(2*j-1)
        out.append(out[-1]*num/den)
    return out


def frozen_d(M):
    qt = qt_coeffs(M)
    s = [mp.mpf(0)]*(M+1); s[0] = mp.mpf(1)
    for m in range(1, M+1):
        s[m] = -mp.fsum(qt[k]*s[m-k] for k in range(1, m+1))
    e = [mp.fsum(s[k]*s[m-k] for k in range(m+1)) for m in range(M+1)]
    return [e[m]*mp.mpf(2)**m for m in range(M+1)]     # d_m


def main():
    M = 160
    d = frozen_d(M)
    a = json.load(open(os.path.join(HERE, "_hm_coeffs_exact.json")))["a_exact"]
    def ak(k):
        n, den = a[k].split("/"); return mp.mpf(int(n))/mp.mpf(int(den))

    print("# frozen d_{2k}  vs  HM a_k*2^k   (same base (9/4)^k; index & constant differ)")
    print(f"  {'k':>3} {'d_{2k}':>16} {'a_k*2^k':>16} {'r_k=ratio':>14} {'r_k/sqrt(2k)':>14}")
    Kmax = min(len(a)-1, M//2)
    rs = []
    for k in list(range(2, 8)) + list(range(Kmax-4, Kmax+1)):
        if k > Kmax: continue
        ak2 = ak(k)*mp.mpf(2)**k
        r = d[2*k]/ak2
        rs.append((k, r/mp.sqrt(2*k)))
        print(f"  {k:>3} {mp.nstr(d[2*k],10):>16} {mp.nstr(ak2,10):>16} "
              f"{mp.nstr(r,10):>14} {mp.nstr(r/mp.sqrt(2*k),10):>14}")

    # extrapolate r_k/sqrt(2k) -> 1/(pi C) = sqrt(3pi/2)
    # (use a few of the largest-k values with a simple 1/k Richardson)
    xs = [mp.mpf(1)/k for k, _ in rs[-5:]]
    ys = [v for _, v in rs[-5:]]
    # Neville to 0
    P = list(ys)
    for j in range(1, len(xs)):
        for i in range(len(xs)-1, j-1, -1):
            P[i] = (xs[i]*P[i-1] - xs[i-j]*P[i])/(xs[i]-xs[i-j])
    lim = P[-1]
    target = 1/(mp.pi*C_HM)
    print(f"\n  r_k/sqrt(2k) -> {mp.nstr(lim,14)}")
    print(f"  1/(pi C)      = {mp.nstr(target,14)}  = sqrt(3pi/2) = {mp.nstr(mp.sqrt(3*mp.pi/2),14)}")
    print(f"  = 1/sqrt(2/(3pi)) : sqrt(2/(3pi)) = {mp.nstr(mp.sqrt(mp.mpf(2)/(3*mp.pi)),12)}")
    print(f"  diff = {mp.nstr(abs(lim-target),3)}")
    print()
    print("  CONCLUSION: frozen (S0=1/pi, gamma'=0)  -->  full HM (C=sqrt(2/(3pi^3)), gamma=-1/2)")
    print("  is a joint (constant x sqrt(2/(3pi)), index -1/2) shift, NOT a lone multiplicative factor.")
    print("  The x-direction fluctuation (one extra continuous mode) supplies BOTH; C is already")
    print("  fixed rigorously by the T1 Painleve-II reduction, so no separate derivation of C is needed.")


if __name__ == "__main__":
    main()
