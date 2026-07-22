"""
MILESTONE 2, Mechanism A (Burekovic-Schafer-Grauer, Handoff-B sub-attack 2):
COMPLEX INSTANTON + ONE-LOOP GAUSSIAN FLUCTUATION DETERMINANT.

Question: the bare real far-tail Freidlin-Wentzell instanton S=s^5/10 gives arg 0. Does a COMPLEX
saddle of the weak-noise problem, DRESSED by its one-loop fluctuation determinant (Gelfand-Yaglom /
functional determinant of the second-variation operator, including the zero mode), produce a
Borel-pair phase near +-54 deg (vs the deterministic connection root arg -45 deg)?

We work with the connection/resonance problem of the Weber skeleton
        u'' = (sign(Y) Y^2 - lambda) u ,
whose deterministic connection root is lambda0 = 0.8896 - 0.8896 i (arg -45 deg). The Borel singularity
of the weak-noise VARIANCE series Var(Y*) = eta^2 (v0 + v1 eta^2 + ...) sits at zeta ~ 1.2-1.8 e^{+-i(54-63) deg}.

STRUCTURE of this script (each piece is an independent, certifiable number):
  A) The COMPLEX WKB (Voros) action between the relevant complex turning-point pair, arg of e^{iP}.
     P(lambda) = oint sqrt(lambda - sign(Y)Y^2) dY. This is the DETERMINISTIC bare-saddle action;
     it controls the connection factor. Reproduce lambda0's phase from it (sanity), then the action
     evaluated between the confining turning points Y=+-sqrt(lambda) (closed-form: pi*lambda/2).
  B) The ONE-LOOP fluctuation determinant around the saddle. Two independent computations:
     (B1) The Gelfand-Yaglom functional determinant of the second-variation operator
          -d^2/dY^2 + (sign(Y)Y^2 - lambda) on the saddle contour, with the zero mode removed.
     (B2) The WKB (Van Vleck / DDJ) one-loop prefactor phase = the Maslov/turning-point phase
          accumulated by the fluctuation determinant = -(1/2) sum over turning points of the
          Stokes rotation, PLUS the amplitude/period phase d(action)/dlambda.
  C) NET effective Borel-pair phase = arg[ saddle action contribution ] combined with the
     one-loop determinant phase, and whether the determinant supplies the ~9-18 deg 45->54 rotation.

numpy/scipy only. Every printed number comes from code run here.
"""
import numpy as np
from scipy.special import loggamma
from scipy.integrate import quad

I = 1j

# ================================================================================================
#  0. The deterministic connection root lambda0 (closed-form Weber Gamma-equation), for reference.
# ================================================================================================
def Gi(lam): return np.exp(loggamma(0.75-0.25j*lam)-loggamma(0.25-0.25j*lam))
def Gr(lam): return np.exp(loggamma(0.75-0.25*lam)-loggamma(0.25-0.25*lam))
def cond(lam): return np.exp(3j*np.pi/4)*Gi(lam)-Gr(lam)
def newton(f, lam, n=80, d=1e-7):
    for _ in range(n):
        c0 = f(lam); dc = 0.5*((f(lam+d)-c0)/d + (f(lam+1j*d)-c0)/(1j*d))
        step = c0/dc; lam = lam - step
        if abs(step) < 1e-13: break
    return lam
LAM0 = newton(cond, 0.9-0.9j)

# ================================================================================================
#  A. The COMPLEX WKB / Voros action.
#     The skeleton has an ASYMMETRIC potential Q(Y)=sign(Y)Y^2. The connection between the
#     recessive (Y>0) and oscillatory (Y<0) sides is controlled by the WKB action integral
#     I(lambda) = int sqrt(Q(Y)-lambda) dY across the classically forbidden (barrier) region and
#     around the complex turning points. We compute several natural action objects and their args.
# ================================================================================================

def action_confining_between_TPs(lam):
    """CONFINING side barrier action between the pair of turning points Y=+-sqrt(lam).
    S = int_{-sqrt lam}^{+sqrt lam} sqrt(lam - Y^2) dY.  Closed form = pi*lam/2 (semicircle).
    This is the 'bulk' classical action for the Y>0 Weber problem analytically continued through
    Y=0; it is the leading exponent of the connection / of the variance's Borel singularity."""
    a = np.sqrt(lam)
    # numeric on straight complex segment -sqrt(lam) -> +sqrt(lam), sqrt branch chosen positive imag.
    t = np.linspace(-1, 1, 4001)
    Y = a * t
    integrand = np.sqrt(lam - Y**2 + 0j)
    S_num = np.trapz(integrand, Y)
    S_closed = np.pi*lam/2
    return S_closed, S_num

def action_oscillatory_TPs(lam):
    """OSCILLATORY side turning points: Q=-Y^2, so -Y^2=lam -> Y=+-i sqrt(lam).
    The relevant action between the two complex conjugate turning points +-i sqrt(lam):
    S_osc = int_{-i sqrt lam}^{+i sqrt lam} sqrt(-Y^2 - lam) dY  along the imaginary axis.
    Substituting Y=i y: sqrt(-(-y^2)-lam)=sqrt(y^2-lam), dY=i dy -> S = i int_{-sqrt lam}^{sqrt lam} sqrt(y^2-lam) dy
    = i * (i pi lam/2) = -pi lam/2  (the mirror of the confining action, with an extra i)."""
    a = np.sqrt(lam)
    t = np.linspace(-1, 1, 4001)
    Y = 1j*a*t            # from -i sqrt(lam) to +i sqrt(lam)
    integrand = np.sqrt(-Y**2 - lam + 0j)
    S_num = np.trapz(integrand, Y)
    S_closed = -np.pi*lam/2
    return S_closed, S_num

def full_voros_period(lam):
    """The full periodic Voros action / classical period that appears in the exact-WKB quantization
    of the ASYMMETRIC well: oint over the closed classical orbit that threads BOTH the confining
    turning points (Y=+-sqrt lam) AND the oscillatory ones. For the Weber connection the relevant
    combination is the difference of the two side-actions (the connection is a RATIO of the two Voros
    symbols). Return the confining action, oscillatory action, and their sum/difference args."""
    Sc, _ = action_confining_between_TPs(lam)
    So, _ = action_oscillatory_TPs(lam)
    return Sc, So

# ================================================================================================
#  B1. GELFAND-YAGLOM functional determinant of the second-variation operator.
#     The one-loop fluctuation determinant around the instanton is det(-d^2/dY^2 + Q(Y) - lambda),
#     with Q(Y)=sign(Y)Y^2. Gelfand-Yaglom: det[ -d^2 + W(Y) ] on [a,b] with Dirichlet BC equals
#     (up to a Y-independent normalization) psi(b), where psi solves -psi'' + W psi = 0, psi(a)=0,
#     psi'(a)=1. The DETERMINANT'S PHASE (arg of the GY endpoint value) is the one-loop phase.
#     We compute it on the complex saddle contour and track its phase as a function of lambda near LAM0.
# ================================================================================================

def gelfand_yaglom(lam, Ya=-8.0, Yb=8.0, n=8000, contour_rot=0.0):
    """Solve the GY IVP  psi'' = (sign(Y)Y^2 - lambda) psi , psi(Ya)=0, psi'(Ya)=1, on a straight
    complex contour from Ya to Yb (optionally rotated by contour_rot to make oscillatory tails decay).
    Returns psi(Yb) = the Gelfand-Yaglom determinant (relative to the free operator normalization).
    NOTE: with sign(Y)Y^2 the two half-lines differ; we keep the real contour split at 0."""
    e = np.exp(1j*contour_rot)
    # contour Y = s*e, s real from Ya to Yb (so sign(Y) uses sign(s))
    s = np.linspace(Ya, Yb, n+1)
    ds = (s[1]-s[0])*e
    psi = 0.0+0j; dpsi = 1.0+0j
    out = np.empty(n+1, complex); out[0]=psi
    def W(sv):
        Y = sv*e
        return np.sign(sv)*Y*Y - lam
    for i in range(n):
        sv = s[i]
        k1 = dpsi;                 l1 = W(sv)*psi
        k2 = dpsi+0.5*ds*l1;       l2 = W(sv+0.5*(s[1]-s[0]))*(psi+0.5*ds*k1)
        k3 = dpsi+0.5*ds*l2;       l3 = W(sv+0.5*(s[1]-s[0]))*(psi+0.5*ds*k2)
        k4 = dpsi+ds*l3;           l4 = W(sv+(s[1]-s[0]))*(psi+ds*k3)
        psi  = psi  + ds/6*(k1+2*k2+2*k3+k4)
        dpsi = dpsi + ds/6*(l1+2*l2+2*l3+l4)
        out[i+1]=psi
    return psi, dpsi, s, out

# ================================================================================================
#  B2. WKB one-loop prefactor phase (Van Vleck / DDJ).
#     For a saddle S(lambda), the one-loop prefactor is  1/sqrt(S''(lambda))  (Van Vleck), and its
#     phase adds to the saddle. Additionally each simple turning point crossed contributes a Stokes
#     rotation e^{-i pi/6} (Airy) or e^{+-i pi/4} (Maslov). We compute S(lambda), S'(lambda),
#     S''(lambda) numerically for the confining barrier action and read off the prefactor phase.
# ================================================================================================

def barrier_action_and_derivs(lam, dl=1e-4):
    """S(lam)=pi*lam/2 confining semicircle action (closed form). Its derivatives are trivial but we
    ALSO compute the OSCILLATORY-side action and the connection-relevant combination numerically."""
    def Sc(l): return np.pi*l/2
    S   = Sc(lam)
    Sp  = (Sc(lam+dl)-Sc(lam-dl))/(2*dl)
    Spp = (Sc(lam+dl)-2*Sc(lam)+Sc(lam-dl))/dl**2
    return S, Sp, Spp


if __name__ == "__main__":
    def a(z): return np.degrees(np.angle(z))
    def sci(z): return f"{z.real:+.5f}{z.imag:+.5f}i"
    print("="*90)
    print("MILESTONE 2, Mechanism A: complex instanton + one-loop fluctuation determinant")
    print("="*90)
    print(f"\n[ref] deterministic connection root lambda0 = {sci(LAM0)}  |.|={abs(LAM0):.4f}  arg={a(LAM0):+.2f} deg")
    print(f"      target Borel-pair phase: ~ +-54 to 63 deg (measured), radius |zeta| ~ 1.2-1.8")

    print("\n" + "-"*90)
    print("A. Complex WKB / Voros action (bare saddle)")
    print("-"*90)
    Sc_cl, Sc_num = action_confining_between_TPs(LAM0)
    So_cl, So_num = action_oscillatory_TPs(LAM0)
    print(f"  Confining barrier action S_c = pi*lam0/2 (semicircle between Y=+-sqrt lam0):")
    print(f"     closed form = {sci(Sc_cl)}   |.|={abs(Sc_cl):.4f}  arg={a(Sc_cl):+.2f} deg")
    print(f"     numeric     = {sci(Sc_num)}  (check)")
    print(f"  Oscillatory action S_o (between Y=+-i sqrt lam0):")
    print(f"     closed form = {sci(So_cl)}   |.|={abs(So_cl):.4f}  arg={a(So_cl):+.2f} deg")
    print(f"     numeric     = {sci(So_num)}  (check)")
    print(f"  arg(S_c) = arg(lam0) + 90 deg? {a(Sc_cl):+.2f} vs {a(LAM0)+90:+.2f}")

    print("\n" + "-"*90)
    print("B1. Gelfand-Yaglom functional determinant of -d^2/dY^2 + (sign(Y)Y^2 - lambda)")
    print("-"*90)
    for rot in (0.0, 0.15, 0.30):
        psiB, dpsiB, s, out = gelfand_yaglom(LAM0, Ya=-7.0, Yb=7.0, n=12000, contour_rot=rot)
        print(f"  contour_rot={rot:.2f}: GY det psi(Yb) = {sci(psiB)}  |.|={abs(psiB):.3e}  arg={a(psiB):+.2f} deg")

    print("\n" + "-"*90)
    print("B2. WKB one-loop prefactor (Van Vleck 1/sqrt(S'')) phase")
    print("-"*90)
    S, Sp, Spp = barrier_action_and_derivs(LAM0)
    print(f"  S(lam0)={sci(S)}  S'={sci(Sp)}  S''={sci(Spp)}")
    if Spp != 0:
        pref = 1.0/np.sqrt(Spp)
        print(f"  Van Vleck prefactor 1/sqrt(S'') = {sci(pref)}  arg={a(pref):+.2f} deg  (S''=0 for linear-in-lam action!)")
    else:
        print("  S'' = 0 (action linear in lambda) -> no Van Vleck rotation from this action")
