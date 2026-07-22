"""
DEPTH — test whether the cusp law W_2 is the k=2 higher-order Tracy-Widom distribution.

Generalized Airy:  Ai_k(x) = (1/pi) int_0^inf cos( t^(2k+1)/(2k+1) + x t ) dt   (k=1 = standard Airy).
Higher-order Airy kernel:  K_k(x,y) = int_0^inf Ai_k(x+u) Ai_k(y+u) du.
Higher-order TW (beta=2):  F_k(s) = det(I - K_k)|_{L^2(s,inf)}   (Bornemann-Nystrom / Gauss-Legendre).
k=1 -> standard TW2 (GUE): mean -1.7711, std 0.9018, skew +0.2241, exkurt +0.0934  (validation anchor).
Then k=2 -> compare cumulants to the cusp W_2: skew +0.61, exk -0.24.
numpy only (leggauss for quadrature; np.linalg.det for the Fredholm determinant).
"""
import numpy as np

# ---------- generalized Airy on a grid ----------
def airy_k_grid(k, xs, T=16.0, Nt=60000):
    """Ai_k(x) = (1/pi) int_0^T cos(t^(2k+1)/(2k+1)+x t) dt, Hann-tapered tail to kill truncation ringing."""
    t = np.linspace(0.0, T, Nt)
    dt = t[1]-t[0]
    taper = np.ones_like(t)
    cut = int(0.8*Nt)
    taper[cut:] = 0.5*(1+np.cos(np.pi*(t[cut:]-t[cut])/(T-t[cut])))
    p = t**(2*k+1)/(2*k+1)
    # phase = p + x*t ; integrate cos(phase) over t for each x
    out = np.empty_like(xs)
    for i, x in enumerate(xs):
        out[i] = np.sum(np.cos(p + x*t)*taper)*dt/np.pi
    return out

def make_airy_interp(k, xlo=-10.0, xhi=20.0, dx=0.01):
    xg = np.arange(xlo, xhi, dx)
    # shorter, well-sampled T for higher k (t^(2k+1) phase oscillates faster; tail cancels anyway)
    T = 16.0 if k == 1 else 7.0
    Nt = 70000
    ag = airy_k_grid(k, xg, T=T, Nt=Nt)
    return xg, ag

def Aik(xq, xg, ag):
    return np.interp(xq, xg, ag, left=0.0, right=0.0)

# ---------- kernel and Fredholm determinant ----------
def kernel_matrix(xnodes, xg, ag, U=18.0, Nu=900):
    u = np.linspace(0.0, U, Nu); du = u[1]-u[0]
    # A[i,m] = Ai_k(xnodes[i] + u[m])
    A = Aik(xnodes[:,None] + u[None,:], xg, ag)
    return A @ (A.T) * du   # K[i,j] = sum_m A[i,m]A[j,m] du

def F_of_s(s, xg, ag, Smax=8.0, Ngl=48):
    xi, wl = np.polynomial.legendre.leggauss(Ngl)
    x = s + (Smax - s)*(xi+1)/2.0
    w = (Smax - s)/2.0*wl
    K = kernel_matrix(x, xg, ag)
    sw = np.sqrt(np.abs(w))
    M = np.eye(len(x)) - (sw[:,None]*sw[None,:])*K
    sign, logdet = np.linalg.slogdet(M)
    return sign*np.exp(logdet)

def cumulants_from_cdf(svals, F):
    F = np.clip(F, 0.0, 1.0)
    f = np.gradient(F, svals)
    f = np.clip(f, 0, None); f /= np.trapz(f, svals)
    m1 = np.trapz(svals*f, svals)
    c = svals - m1
    m2 = np.trapz(c**2*f, svals); sd = np.sqrt(m2)
    z = c/sd
    mu3 = np.trapz(z**3*f, svals); mu4 = np.trapz(z**4*f, svals)
    return m1, sd, mu3, mu4-3.0

if __name__ == "__main__":
    import sys
    k = int(sys.argv[1]) if len(sys.argv)>1 else 1
    # validate generalized Airy for k=1 against known Airy values
    if k == 1:
        xg, ag = make_airy_interp(1)
        for xv, ref in [(0.0,0.355028),(-1.0,0.535561),(2.0,0.034924),(-2.0,0.227407)]:
            print(f"  Ai(x={xv:+.0f}) = {Aik(np.array([xv]),xg,ag)[0]:+.5f}   ref {ref:+.5f}")
    else:
        xg, ag = make_airy_interp(k)
    # Fredholm determinant F_k(s) over a grid -> cumulants
    svals = np.arange(-9.0, 6.01, 0.10)
    F = np.array([F_of_s(s, xg, ag) for s in svals])
    m1, sd, sk, exk = cumulants_from_cdf(svals, F)
    print(f"\nk={k}  higher-order TW (beta=2):")
    print(f"  mean={m1:+.4f}  std={sd:.4f}  skew={sk:+.4f}  exkurt={exk:+.4f}")
    if k==1: print("  TW2 ref:   mean -1.7711  std 0.9018  skew +0.2241  exkurt +0.0934")
    if k==2: print("  cusp W_2:  skew +0.61  exk -0.24  (target if reframe holds)")
    # left-tail exponent of F_k:  -ln F(s) ~ |s|^p  as s->-inf  (TW: p=3=2k+1; multicritical order)
    Fc = np.clip(F, 1e-12, 1.0)
    sel = (Fc > 1e-6) & (Fc < 0.05) & (svals < m1)
    if sel.sum() >= 4:
        t = np.log(-np.log(Fc[sel])); xlg = np.log(-svals[sel])
        p = np.polyfit(xlg, t, 1)[0]
        print(f"  left-tail exponent fit  -lnF ~ |s|^{p:.2f}   (order 2k+1 = {2*k+1})")
    np.save(f"hoTW_F_k{k}.npy", np.vstack([svals, F]))
