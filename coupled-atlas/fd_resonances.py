"""
Independent validation of the complex-scaled resonances: finite-difference diagonalization of
H_theta = -e^{-2i th} d^2/dt^2 + e^{2i th} sign(t) t^2  on a rotated contour grid.
Genuine resonances appear as theta-INDEPENDENT complex eigenvalues (Im<0). Compare to the Jost-Wronskian
values 0.86-0.82i, 2.30-1.22i, 4.14-1.08i. numpy.linalg.eig (independent of the Wronskian construction).
"""
import numpy as np
def Heig(th, L=9.0, N=900):
    t=np.linspace(-L,L,N); h=t[1]-t[0]
    kin=-np.exp(-2j*th)/h**2
    H=np.zeros((N,N),complex)
    np.fill_diagonal(H, -2*kin + np.exp(2j*th)*np.sign(t)*t*t)
    idx=np.arange(N-1); H[idx,idx+1]=kin; H[idx+1,idx]=kin
    w=np.linalg.eigvals(H)
    w=w[(w.imag<-0.2)&(w.imag>-2.2)&(w.real>0.2)&(w.real<5.0)]
    return w[np.argsort(w.real)]
pin=np.array([0.86-0.82j,2.30-1.22j,4.14-1.08j])
for th in [0.40,0.55]:
    w=Heig(th)
    # match each pinned to nearest FD eigenvalue
    print(f"theta={th:.2f}: FD eigenvalues (Im<0) near pinned resonances:")
    for p in pin:
        k=np.argmin(np.abs(w-p));
        print(f"   pinned {p.real:.2f}{p.imag:+.2f}i  ->  FD {w[k].real:+.3f}{w[k].imag:+.3f}i   (dist {abs(w[k]-p):.3f})")
