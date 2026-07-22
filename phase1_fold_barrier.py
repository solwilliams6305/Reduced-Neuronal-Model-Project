"""
Phase 1 — Closed-form fold barrier, non-degenerate (isotropic) noise.
Independent symbolic re-derivation by TWO routes, checked against the stored
project constants (FOLD_ESCAPE_PREFACTORS.md): K_fold = 4/3, A0 = 1/pi,
escape rate ~ (sqrt(delta)/pi) exp(-(8/3) delta^{3/2}/sigma^2).

Route A  (HJ / gradient + Kramers): the frozen fast-subsystem potential barrier.
Route B  (Freidlin-Wentzell instanton): minimum-action escape path; action = quasipotential.

Convention (matches Cameron / the project files):
  SDE  dx = b(x) dt + sigma dW ,   rate ~ exp(-V/sigma^2),
  Hamiltonian H(x,p) = 1/2 |p|^2 + b.p ,  V solves H(x, grad V)=0.
  For a gradient fast field b = -U', the quasipotential V = 2U, so the
  escape *action* (exponent) is  Delta V = 2 Delta U.
"""
import sympy as sp

x, V, W, p = sp.symbols('x V W p', real=True)
delta = sp.symbols('delta', positive=True)   # distance to fold > 0
sqrt = sp.sqrt
ok = True
def check(name, got, want):
    global ok
    eq = sp.simplify(got - want) == 0
    ok = ok and eq
    print(f"  [{'PASS' if eq else 'FAIL'}]  {name}:  got {sp.nsimplify(got)}   want {want}")

print("="*72)
print("ROUTE A — frozen fast-subsystem potential barrier (HJ/gradient side)")
print("="*72)
# Fast subsystem at frozen slow value delta>0 (distance to fold):
#   x' = b(x) = x^2 - delta = -U'(x)   =>   U(x) = delta*x - x^3/3
U = delta*x - x**3/3
Up = sp.diff(U, x)
print("  U(x) = delta*x - x^3/3 ,  U'(x) =", Up)
# critical points x = -+ sqrt(delta):  well at x=-sqrt(delta), saddle at x=+sqrt(delta)
xwell, xsaddle = -sqrt(delta), +sqrt(delta)
Upp = sp.diff(U, x, 2)
print("  U''(well)   =", sp.simplify(Upp.subs(x, xwell)),  "(>0 => minimum, well)")
print("  U''(saddle) =", sp.simplify(Upp.subs(x, xsaddle)),"(<0 => maximum, saddle)")
dU = sp.simplify(U.subs(x, xsaddle) - U.subs(x, xwell))
check("barrier  Delta U", dU, sp.Rational(4,3)*delta**sp.Rational(3,2))
check("K_fold (= DeltaU / delta^{3/2})", sp.simplify(dU/delta**sp.Rational(3,2)), sp.Rational(4,3))

# Kramers prefactor for the cubic barrier:  A0 = sqrt(U''_min |U''_max|)/(2 pi) / sqrt(delta)
kmin = sp.Abs(Upp.subs(x, xwell)); kmax = sp.Abs(Upp.subs(x, xsaddle))
A0 = sp.simplify(sqrt(kmin*kmax)/(2*sp.pi) / sqrt(delta))   # strip the sqrt(delta) scaling
check("Kramers prefactor A0", A0, 1/sp.pi)

print()
print("="*72)
print("ROUTE B — Freidlin-Wentzell instanton (minimum-action escape path)")
print("="*72)
# Hamiltonian H(x,p) = 1/2 p^2 + b(x) p ,  b(x) = x^2 - delta.
b = x**2 - delta
H = sp.Rational(1,2)*p**2 + b*p
# zero-energy manifold H=0 -> p=0 (deterministic relaxation) or p = -2 b (fluctuation path)
p_fluct = sp.solve(H, p)
print("  H=0  =>  p in", p_fluct, " (p=0 relaxation; p=-2b is the uphill instanton)")
p_inst = sp.Rational(-2)*b              # = 2(delta - x^2) >= 0 on (-sqrt d, sqrt d)
# instanton action S = integral p dx along the escape path well -> saddle
S = sp.integrate(p_inst, (x, xwell, xsaddle))
S = sp.simplify(S)
check("instanton action S = ∫ p dx", S, sp.Rational(8,3)*delta**sp.Rational(3,2))
check("S = 2 * DeltaU  (action = quasipotential exponent)", sp.simplify(S - 2*dU), 0)

print()
print("="*72)
print("CROSS-CHECK against stored project constants (FOLD_ESCAPE_PREFACTORS.md)")
print("="*72)
print("  stored:  K_fold = 4/3 ,  A0 = 1/pi ,  rate ~ (sqrt(delta)/pi) exp(-(8/3) delta^{3/2}/sigma^2)")
rate_exponent = sp.simplify(S)   # exp(-S/sigma^2) = exp(-(8/3) delta^{3/2}/sigma^2)
check("rate exponent coefficient (S)", rate_exponent, sp.Rational(8,3)*delta**sp.Rational(3,2))
print()
print("BOTH ROUTES AGREE." if ok else "DISCREPANCY — investigate.")
print("Summary:  DeltaU = (4/3) delta^{3/2};  quasipotential barrier DeltaV = S = (8/3) delta^{3/2};")
print("          escape rate lambda(delta) = (sqrt(delta)/pi) exp(-(8/3) delta^{3/2}/sigma^2).")
import sys; sys.exit(0 if ok else 1)
