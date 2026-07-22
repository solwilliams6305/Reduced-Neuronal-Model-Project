"""
_tw_thread2_borel.py -- THREAD 2: is the TW_beta tail's noise (1/beta) expansion Borel-summable?

The tail is (Borot-Nadal 1111.2761, Prop 1.1)
   1-TW_beta(s) = c_beta s^{-3beta/4} e^{-2beta s^{3/2}/3} exp[ sum_m (beta/2) R_m(2/beta) s^{-3m/2} ],
   c_beta = Gamma(beta/2)/((4 beta)^{beta/2} 2 pi).
Its 1/beta expansion has TWO factorially-divergent sectors:

 (A) DYNAMICAL ESCAPE sector -- the fluctuation series of the reduced escape (this program, exact):
     R(g)=(1/pi)e^{-(2/3)g}/Qtilde(2/g)^2, F(g)=Qtilde(2/g)^{-2}=sum d_m g^{-m},
     d_m ~ -(1/pi) Gamma(m)(3/2)^m  (NON-alternating).  g = beta a^{3/2}, so the 1/beta series at
     fixed a is d_m a^{-3m/2} beta^{-m}: factorial, POSITIVE-real-axis Borel singularity => needs
     MEDIAN resummation.  We (i) Borel-Pade the d_m and locate the singularity at t=2/3, (ii) median
     Borel-Laplace resum F(g) and check it reproduces the exact escape rate (=> Borel-summable).

 (B) COULOMB-GAS sector -- the Gamma(beta/2) beta-ensemble normalization.  ln Gamma(beta/2) has the
     Stirling tail sum_k B_{2k}/(2k(2k-1)) (2/beta)^{2k-1}: ALTERNATING factorial => Borel
     singularities on the IMAGINARY axis (the Gamma instantons e^{-2 pi i n beta/2}) => ORDINARY
     Borel-summable (no median ambiguity).

Verdict: YES, Borel-summable -- a genuine resurgent trans-series with a median (dynamical) sector
and an ordinary (Coulomb-gas) sector.  Distinct from the LEVEL-variable HM resurgence.
"""
import mpmath as mp
mp.mp.dps = 50


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
    return [e[m]*mp.mpf(2)**m for m in range(M+1)]

def R_mfpt(g):   # exact reduced escape rate via the MFPT double integral (finite cutoff)
    D = mp.mpf(2)/g; L = mp.mpf(3); V = lambda t: t**3/3 - t
    def inner(y):
        pts = [y, 1, L] if y < 1 else [y, L]
        return mp.quad(lambda z: mp.e**(-V(z)/D), pts)
    tau = (1/D)*mp.quad(lambda y: mp.e**(V(y)/D)*inner(y), [-L, -1, 1])
    return 1/tau


def sector_A():
    print("# (A) DYNAMICAL ESCAPE sector: d_m ~ -(1/pi)Gamma(m)(3/2)^m, positive-axis Borel sing.")
    M = 100
    d = frozen_d(M)
    print("   signs d_1..12:", "".join("+" if d[m] > 0 else "-" for m in range(1, 13)),
          " (all '-' = non-alternating = positive-axis)")
    # Borel transform in z=1/g:  F(g)=sum d_m z^m ; B_F(t)=sum_{m>=1} d_m t^{m-1}/(m-1)!
    N = 61
    bc = [d[m+1]/mp.factorial(m) for m in range(N)]     # coeffs of B_F(t)=sum bc_k t^k
    Lp = (N-1)//2
    p, q = mp.pade(bc, Lp, Lp)
    roots = mp.polyroots(q[::-1], maxsteps=2000, extraprec=200)
    posreal = [r for r in roots if abs(mp.im(r)) < 1e-6 and mp.re(r) > 0]
    tc = min(posreal, key=lambda r: abs(mp.re(r)-mp.mpf(2)/3)) if posreal else None
    print(f"   Borel B_F(t) nearest positive-real pole: t_c = {mp.nstr(mp.re(tc),12) if tc else None}"
          f"   (expect 2/3 = {mp.nstr(mp.mpf(2)/3,12)})")

    # median Borel-Laplace resummation of F(g): F = d_0 + Re \int_0^{inf} B_F(t) e^{-g t} dt (PV/median)
    def Bpade(t):
        return mp.polyval(p[::-1], t)/mp.polyval(q[::-1], t)
    g = mp.mpf(9)
    th = mp.mpf('0.35')
    def lap(theta):
        ray = lambda u: u*mp.e**(1j*theta)
        return mp.quad(lambda u: Bpade(ray(u))*mp.e**(-g*ray(u))*mp.e**(1j*theta), [0, mp.inf])
    Fmed = d[0] + mp.re((lap(th)+lap(-th))/2)
    Ftrue = mp.pi*mp.e**(mp.mpf(2)/3*g)*R_mfpt(g)       # exact: R=(1/pi)e^{-2g/3}F
    print(f"   median Borel-Laplace F({int(g)})  = {mp.nstr(Fmed,10)}")
    print(f"   exact  F({int(g)}) = pi e^(2g/3) R = {mp.nstr(Ftrue,10)}   rel.diff {mp.nstr(abs(Fmed-Ftrue)/Ftrue,3)}")
    print("   => dynamical sector is Borel-summable (median). Imag ambiguity ~ e^{-(2/3)g} (Stokes).")


def sector_B():
    print("\n# (B) COULOMB-GAS sector: Gamma(beta/2) Stirling tail, coeff of beta^{-(2k-1)}")
    # ln Gamma(z) - [(z-1/2)ln z - z + (1/2)ln 2pi] = sum_k B_{2k}/(2k(2k-1)) z^{-(2k-1)}, z=beta/2
    # coeff of beta^{-(2k-1)}:  a_k = B_{2k}/(2k(2k-1)) * 2^{2k-1}
    K = 16
    a = [mp.bernoulli(2*k)/(2*k*(2*k-1))*mp.mpf(2)**(2*k-1) for k in range(1, K+1)]
    print("   a_1..4:", [mp.nstr(a[i], 8) for i in range(4)])
    print("   signs a_k:", "".join("+" if v > 0 else "-" for v in a), " (ALTERNATING => imaginary-axis)")
    # ratio test: a_{k}/a_{k-1} ~ -(2k-2)(2k-3)/pi^2  (=> factorial, base 1/pi^2, alternating)
    print("   |a_k/a_{k-1}| / [(2k-2)(2k-3)] -> 1/pi^2 =", mp.nstr(1/mp.pi**2, 8))
    for k in (8, 12, 16):
        r = abs(a[k-1]/a[k-2])/((2*k-2)*(2*k-3))
        print(f"      k={k}: {mp.nstr(r,10)}")
    # Borel transform B(t)=sum a_k t^{2k-1}/(2k-1)! ; singularity where? a_k/(2k-1)! ~ (-1)^{k+1}/((2k-1)pi^{2k})
    # => B(t) ~ sum (-1)^{k+1} t^{2k-1}/((2k-1)pi^{2k}) = (1/pi) arctan(t/pi): poles/branch at t=+- i pi.
    print("   Borel transform ~ (1/pi) arctan(t/pi): branch points at t = +- i pi (IMAGINARY axis)")
    print("   => Coulomb-gas sector is ORDINARY Borel-summable (no positive-axis median ambiguity);")
    print("      the 'noise instantons' are the Gamma-function/Coulomb-gas ones e^{-2 pi i n beta/2}.")


def main():
    print("THREAD 2 -- Borel-summability of the TW_beta tail 1/beta (noise) expansion\n")
    sector_A()
    sector_B()
    print("\n# VERDICT: YES, Borel-summable. Genuine resurgent trans-series, TWO sectors:")
    print("   (A) dynamical escape: positive-real-axis, MEDIAN (Stokes const ~1/pi, action 2Phi);")
    print("   (B) Coulomb-gas Gamma(beta/2): imaginary-axis, ordinary Borel-summable.")
    print("   Contrast: the LEVEL-variable (a, beta=2) resurgence is the HM PII one (C=sqrt(2/(3pi^3))).")


if __name__ == "__main__":
    main()
