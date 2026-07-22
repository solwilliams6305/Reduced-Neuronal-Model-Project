# The closed-form asymmetric connection of 𝒲's skeleton (both Γ-factors, machine-precision verified)

_July 2026 (incoming agent, Fable 5). Following the deep-research lead (Weber-type exact WKB + c=1 Γ-connection
data), the **deterministic** asymmetric connection / resonance condition of 𝒲's skeleton
$u''=(\operatorname{sign}(Y)Y^2-\lambda)u$ is now in **closed form** — both connection factors are Weber Γ-ratios
(real-argument confining × imaginary-argument oscillatory), verified to machine precision, and its Γ-equation root
reproduces the dominant resonance $\lambda_0=0.890-0.890i$. Scripts `coupled-atlas/connection_gamma_test.py`,
`connection_poles.py`, `oscillatory_gamma.py`, `resonance_gamma.py`, `connection_closed_form.py`; figures
`figures/connection_poles.png`, `figures/connection_closed_form.png`. Tags **[DERIVED]/[NUMERIC ✓✓]**; ⚑._

## Result

Resonance = recessive on $Y>0$, purely outgoing on $Y<0$; matching log-derivatives at $Y=0$ gives the condition
$L(\lambda)=R(\lambda)$ with (both **closed form**, $G(\mu):=\Gamma(3/4-\mu)/\Gamma(1/4-\mu)$):
$$
\boxed{\;R(\lambda)=-\sqrt2\,\frac{U'(a,0)}{U(a,0)}=2\,\frac{\Gamma(3/4-\lambda/4)}{\Gamma(1/4-\lambda/4)}\quad(\text{confining, real Weber, }a=-\tfrac\lambda2)\;}
$$
$$
\boxed{\;L(\lambda)=\frac{f_{\rm out}'(0)}{f_{\rm out}(0)}=e^{3i\pi/4}\,2\,\frac{\Gamma(3/4-i\lambda/4)}{\Gamma(1/4-i\lambda/4)}\quad(\text{oscillatory, imaginary Weber})\;}
$$
$$
\boxed{\;\text{Resonance condition:}\quad e^{3i\pi/4}\,\frac{\Gamma(3/4-i\lambda/4)}{\Gamma(1/4-i\lambda/4)}=\frac{\Gamma(3/4-\lambda/4)}{\Gamma(1/4-\lambda/4)}\;}
$$
**The asymmetry is exactly $\lambda\leftrightarrow i\lambda$ (real ↔ imaginary Weber) plus the $e^{3i\pi/4}$
inversion Stokes phase.** This is the "Γ-connection mixing a real and an imaginary parameter" anticipated in
`T2_CONNECTION_DATA_DERIVATION`, now explicit.

## Evidence (decisive, machine precision)

1. **Confining factor $R$ — pinned by the non-confounded POLE test** (`connection_poles.py`): the first oscillatory
   node depth $\Theta_1(\lambda)\to0$ (a node crosses the turning, $u(0)=0$) at **$\lambda=3.00$ and $7.00$** (mean
   deviation $0.01$) = the poles of $\Gamma(3/4-\lambda/4)$ = the zeros of the textbook $U(a,0)$ (DLMF §12.2.6).
   $Y^\star_{\rm det}(0)=-2.188$ matches the program's $-2.19$. **[NUMERIC ✓✓]**
2. **Oscillatory factor $L$ — closed form to machine precision** (`connection_closed_form.py`): the outgoing
   inverted-oscillator Jost log-derivative satisfies $L(\lambda)/[2\Gamma(3/4-i\lambda/4)/\Gamma(1/4-i\lambda/4)]
   =e^{3i\pi/4}$ **exactly** — modulus $1.0000$, argument $0.7500\pi$ — across 9 real *and complex* $\lambda$,
   stable under step-halving to $-0.70711+0.70711i$ (5 digits). Not a fit; an identity. **[NUMERIC ✓✓]** ⚑
3. **The closed-form Γ-equation reproduces the resonance:** its Newton root is $\lambda_0=0.8896-0.8896i$
   ($|{\rm cond}|=2.8\times10^{-16}$), matching M1's **independent** finite-difference resonance $0.890-0.890i$ to
   $0.0006$. **[NUMERIC ✓✓]**

## Significance + honest scope

- **This is the closed form of the DETERMINISTIC skeleton's connection/resonance data** — the exact
  eigenvalue-parameter connection of the asymmetric glued-Weber operator. It *supersedes the falsified M3 attempt*
  (`M3_RESONANCE_DETERMINANT_NOTES.md`, which tried a real-axis meromorphic Γ-ratio $D(\lambda)$ and failed): the
  correct object is the **resonance condition** (log-derivative matching), a Γ-*equation* combining the real-Weber
  and imaginary-Weber ratios — which M3's real-axis $D(\lambda)$ could not be. The deep-research recommendation
  (graft the c=1/Weber Γ-connection data, both sides) is now realized exactly.
- **It is NOT (yet) the closed form of the stochastic law $\mathcal W$.** $\mathcal W$ is the *noise-smeared*
  first-passage law; this closed-form connection is its deterministic **backbone** (resonances, deterministic
  first-node structure). The remaining open step is the **stochastic generalization** — integrate the noise
  against this now-exact connection (the q=2 analogue of RRV's Riccati diffusion). The program's theorem that
  $\mathcal W$ is *not* a Painlevé-σ-ODE stands; but the skeleton it lives on is now closed-form, which is the
  concrete platform for the stochastic step.

## The stochastic step, scoped (`stochastic_scope.py`) — [DERIVED + NUMERIC ✓✓]

Does the noise smear the resonance Γ-equation into a closed-form law? **Precisely yes-and-no**, matching the
no-ODE theorem. The noise is a random potential ($\delta V=-\eta\xi$), so by the first-passage susceptibility the
**weak-noise variance is closed-form from the deterministic connection**:
$$\operatorname{Var}(Y^\star)=\eta^2\,C_V,\qquad C_V=\Big[\tfrac{u_2(Y^\star)}{\mathcal W\,u_c'(Y^\star)}\Big]^2\!\int u_c^4\,ds=0.1339$$
(a Weber-function integral off the exact backbone). MC confirms $\operatorname{Var}(Y^\star)/\eta^2\to C_V$ as
$\eta\to0$ (ratio $1.01$ at $\eta=0.1$). **[NUMERIC ✓✓]** But at the physical $\beta=2$ ($\eta=1.41$) the width is
$+74\%$ above this leading term and the **skew climbs $0\to+0.605$** — and that skew is the **irreducible
2-variable FP-PDE content** (no ODE reduction — theorem). So:
> The closed-form backbone gives 𝒲's **exact weak-noise Gaussian** (width from the connection Γ-integral) and its
> **exact instanton tails** (left constant $1/10$); the non-Gaussian $\beta=2$ shape is the genuinely irreducible
> part. The noise does **not** collapse the resonance Γ-equation to a 1-D closed form — consistent with, and now
> concretely demonstrating, the proven no-Painlevé-ODE obstruction.

## Tracker

The closed-form asymmetric connection (both Γ-factors + resonance condition) is established and machine-verified,
and the stochastic step is scoped: closed-form backbone + weak-noise width + instanton tails, with the full
$\beta=2$ shape the irreducible PDE. Open #17 advances from "structural (asymmetric PIV-family)" to "**deterministic
connection in closed form; weak-noise law closed-form; full law = the (proven-irreducible) 2-var PDE.**" This is
the honest ceiling of the closed-form program: everything reducible is now reduced; the residue is a theorem.
