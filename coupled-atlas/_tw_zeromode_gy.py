"""
_tw_zeromode_gy.py -- derive the x-extension factor sqrt(2/(3pi)) as the Gelfand-Yaglom
zero-mode (collective-coordinate) Jacobian of the escape instanton.

Reduced cubic-barrier escape (PROGRAM2_TWBETA_S_OF_A_NOTES.md): dP=-(P^2-1)dtau+2 eps_eff dW,
V(P)=P^3/3 - P (well P=+1, barrier P=-1, height Delta V = V(-1)-V(+1) = 4/3 = 2 A_g, A_g=2/3).

Instanton = uphill trajectory dP/dtau = V'(P) = P^2-1, i.e. the TANH KINK
    P*(tau) = -tanh(tau)        (P: +1 at tau=-inf  ->  -1 at tau=+inf),
translation zero mode  psi0 = dP*/dtau = -sech^2(tau).

Gelfand-Yaglom fact: a fluctuation determinant with a bosonic zero mode psi0 is regularized as
det'(M) with the zero eigenvalue removed; trading the zero mode for its collective coordinate
introduces the Jacobian  ||psi0||/sqrt(2 pi hbar).  Its hbar-independent (Stokes-constant) part is
    J_zm = sqrt( ||psi0||^2 / (2 pi) ),      ||psi0||^2 = \int psi0^2 dtau = \int sech^4 = 4/3 = Delta V
so  J_zm = sqrt( (4/3)/(2 pi) ) = sqrt(2/(3 pi)) = sqrt(A_g/pi).
The accompanying hbar^{1/2} becomes the -1/2 large-order index shift (one continuous zero mode).

Hence the full level-variable (Hastings-McLeod) Stokes constant factorizes as
    C = S0 * J_zm = (1/pi) * sqrt(2/(3 pi)) = sqrt(2/(3 pi^3)),
with S0 = 1/pi the transverse (Kramers, P-direction) determinant [the frozen escape RATE, which
factors the zero mode out], and J_zm the escape-instanton zero-mode Jacobian [which the level-
variable resurgence counts back in].
"""
import mpmath as mp
mp.mp.dps = 40

A_g   = mp.mpf(2)/3
C_HM  = mp.sqrt(mp.mpf(2)/(3*mp.pi**3))
S0    = 1/mp.pi

# ---- 1. tanh kink solves the instanton ODE dP/dtau = P^2 - 1 ----
P  = lambda t: -mp.tanh(t)
dP = lambda t: -mp.sech(t)**2
res = lambda t: dP(t) - (P(t)**2 - 1)
print("# 1. instanton P*=-tanh solves dP/dtau=P^2-1 :",
      "max|resid|", mp.nstr(max(abs(res(t)) for t in (-2,-1,0,1,2)), 3))

# ---- 2. barrier height Delta V and the virial identity ||psi0||^2 = int V' dP = Delta V ----
V   = lambda p: p**3/3 - p
dV  = mp.quad(lambda p: (p**2 - 1), [1, -1])     # \int_{+1}^{-1} V'(P) dP = V(-1)-V(+1)
print("\n# 2. Delta V and zero-mode norm")
print("   Delta V = V(-1)-V(+1)        =", mp.nstr(V(-1)-V(1), 12), " (=4/3)")
print("   int_{+1}^{-1} V'(P) dP       =", mp.nstr(dV, 12), " (virial: = Delta V)")
norm2 = mp.quad(lambda t: mp.sech(t)**4, [-mp.inf, 0, mp.inf])   # ||psi0||^2 = int sech^4
print("   ||psi0||^2 = int sech^4 dtau =", mp.nstr(norm2, 12), " (= 4/3 = Delta V = 2 A_g)")
print("   checks:  4/3 =", mp.nstr(mp.mpf(4)/3, 12),
      "  |norm2-4/3| =", mp.nstr(abs(norm2 - mp.mpf(4)/3), 3))

# ---- 3. the Gelfand-Yaglom zero-mode Jacobian ----
J_zm = mp.sqrt(norm2/(2*mp.pi))
print("\n# 3. zero-mode Jacobian  J_zm = sqrt(||psi0||^2/(2 pi))")
print("   J_zm            =", mp.nstr(J_zm, 20))
print("   sqrt(2/(3 pi))  =", mp.nstr(mp.sqrt(mp.mpf(2)/(3*mp.pi)), 20))
print("   sqrt(A_g/pi)    =", mp.nstr(mp.sqrt(A_g/mp.pi), 20))
print("   diff            =", mp.nstr(abs(J_zm - mp.sqrt(mp.mpf(2)/(3*mp.pi))), 3))

# ---- 4. assemble the full Stokes constant ----
print("\n# 4. C = S0 * J_zm")
C_assembled = S0 * J_zm
print("   S0 (frozen, transverse Kramers) =", mp.nstr(S0, 16), " (=1/pi)")
print("   J_zm (escape zero mode)         =", mp.nstr(J_zm, 16))
print("   C = S0*J_zm                     =", mp.nstr(C_assembled, 20))
print("   C_HM = sqrt(2/(3 pi^3)) [Dunne] =", mp.nstr(C_HM, 20))
print("   |C_assembled - C_HM|            =", mp.nstr(abs(C_assembled - C_HM), 3))
print("\n   => sqrt(2/(3pi)) DERIVED as the GY zero-mode Jacobian sqrt(||psi0||^2/2pi),")
print("      ||psi0||^2 = int sech^4 = 4/3 = Delta V (virial); index shift -1/2 = one zero mode.")
