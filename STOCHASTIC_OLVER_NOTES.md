# The stochastic-Olver uniform-connection bound — the load-bearing estimate that closes T1

_June 2026. Figure `coupled-atlas/figures/connection_variance.png`; script `connection_variance.py`.
The single estimate the whole T1 tube reduces to. Partial rigorous progress + the validated uniformity +
the precise residual. Tags [CITED]/[DERIVED]/[NUMERIC]/[OPEN]; ⚑ = load-bearing._

## 1. Precise statement

Inner equation $u''=(V_\Delta(Y)-\lambda-\eta\dot W)u$, $V_\Delta=\operatorname{sign}(Y)|Y|(|Y|+\Delta)$, with
the recessive solution at $Y\to+\infty$ continued through the turning to the oscillatory region $Y<0$. Let
$\Theta_\Delta(\eta)$ be the **noisy Prüfer connection phase** accumulated across the turning window
$[Y_-,Y_+]$ (entrance recessive, exit oscillatory), and $\delta\Theta=\Theta_\Delta(\eta)-\Theta_\Delta(0)$
its noise-induced fluctuation about the deterministic Weber connection.

> **Stochastic-Olver bound (target).** There are constants $C_m<\infty$, **independent of $\Delta$**, with
> $$\mathbb E\big[(\delta\Theta)^{2m}\big]\ \le\ C_m\,\eta^{2m}\qquad(m=1,2,\dots),\quad\text{uniformly in }
> \rho=\Delta/\ell\in(0,\infty).$$
> Equivalently: the noisy connection coefficient is a **sub-Gaussian** perturbation of the deterministic
> (Olver/Weber) connection, with **width $O(\eta)$ uniform through the merge**. The norm is the scalar
> $L^{2m}(\mathbb P)$ moment of the connection phase; "uniform in $\rho$" = constants do **not** degenerate
> as the two turning structures coalesce ($\Delta\to0$, simple-Airy $\to$ quadratic-Weber).

The $m=1$ case (the **variance**) is what the sub-Gaussian tube (T1) needs at leading order; all $m$ give
the full Gaussian tube.

## 2. Proof skeleton

**Step 0 [CITED].** Olver's deterministic uniform parabolic-cylinder connection: for $u''=(V_\Delta-\lambda)u$,
the recessive solution's continuation through the turning is given by a **parabolic-cylinder (Weber)
parametrix with explicit error bounds, uniform as the turning structure varies** (Olver, *Asymptotics and
Special Functions* Ch. 11; the two-coalescing-turning-points / PC case), after a change of variable reducing
$V_\Delta$ to the Olver normal form. This fixes $\Theta_\Delta(0)$ and the deterministic phase, uniform in $\Delta$.

**Step 1 [DERIVED].** Prüfer phase SDE in the oscillatory region ($k=\sqrt{|\lambda-V_\Delta|}$, $p=k\cot\theta$):
$$\theta'=k+\frac{k'}{2k}\sin2\theta+\frac{\eta\sin^2\theta}{k}\,\dot W .$$
Linearising, the noise-induced fluctuation is $\displaystyle\delta\theta(Y)=\eta\!\int\frac{\sin^2\theta_{\rm det}}{k}\,dW+(\text{nonlinear, Itô})$.

**Step 2 [DERIVED] — outer variance, bounded at both ends.**
$$\operatorname{Var}(\delta\Theta_{\rm out})=\eta^2\!\int_{\rm out}\frac{\mathbb E[\sin^4\theta]}{k^2}\,dY
\approx\tfrac{3}{8}\,\eta^2\!\int_{\rm out}\frac{dY}{|\lambda-V_\Delta|},\quad\text{cut at the turning width }w_t.$$
- **Airy end ($\rho\gg1$, simple turning slope $s\!\approx\!\Delta$):** $\int\sim\frac1s\ln(L/w_t)$,
  $w_t\sim s^{-1/3}$ ⇒ $\operatorname{Var}\sim\eta^2\frac{\ln\Delta}{\Delta}\to0$.
- **Weber/cusp end ($\rho\to0$, quadratic turning):** $k^2\approx Y^2$, $\int_{w_t}^L\frac{dY}{Y^2}\sim1/w_t$,
  $w_t\sim O(1)$ (cusp WKB width) ⇒ $\operatorname{Var}\sim\eta^2\cdot O(1)$.
- ⇒ the **outer** variance is $O(\eta^2)$ at both ends; bounded uniformly in $\rho$ (modulo a possible mild log).

**Step 3 [OPEN ⚑ — the single hardest sub-step].** **Inner PC region** $|Y-Y_t|\lesssim w_t$, where the
WKB/Prüfer picture fails and Olver's PC parametrix governs. Here the noise coefficient $\eta\sin^2\theta/k$
is singular ($k\to0$), and the white-noise perturbation must be carried through the **PC parametrix** and the
connection-coefficient fluctuation bounded **uniformly in $\Delta$** as the turning degenerates. This is the
load-bearing step: Olver's bounds are deterministic and for a fixed turning structure; the upgrade needs
(i) the stochastic perturbation of the PC parametrix (a singular-noise perturbation of $U(a,z)$), and
(ii) uniformity as simple→quadratic ($\Delta\to0$). **Not done.**

**Where deterministic uniformity survives / fails.** It **survives** in the outer region (Step 2: the WKB
phase-diffusion variance is bounded, the $1/k^2$ singularity integrable away from $w_t$). It is **at risk**
in the inner PC region (Step 3): the noise is singular exactly where Olver's parametrix replaces WKB, and the
uniform-in-$\Delta$ control there is the open content.

## 3. Numerical validation [NUMERIC]

Measured the connection variance $\operatorname{Var}(\delta p)$ at the turning entrance for $V_\Delta$,
$\eta\in\{0.7,1,1.4\}$, $\rho\in[0.04,3.2]$:
- **$O(\eta^2)$, uniform in $\rho$:** $\operatorname{Var}/\eta^2$ stays in $[0.20,0.47]$ across the whole
  $\rho$ range — **bounded, no blow-up at the merge** ($\rho\to0$); the three $\eta$ curves collapse
  (η²-scaling). **[NUMERIC ✓ — the $m=1$ content of the bound]**
- **The OU quasi-static comparison diverges** at the merge ($\eta^2/4\sqrt{V}\to\infty$ as $V\to0$: 3.3, 6.7,
  13.1 at the smallest $\rho$) — confirming the uniformity is **non-trivial** and genuinely requires the
  Weber/PC structure (not the frozen-OU one). **[NUMERIC ✓]**

So the **variance** content of the stochastic-Olver bound is confirmed; what numerics cannot supply is the
*rigorous* inner-PC control (Step 3) or the higher-moment (sub-Gaussian) constants.

## 4. How close this brings T1 to closed

| T1 component | status |
|---|---|
| approach tube (C1), both sides, uniform in ρ | **PROVED** |
| re-entry closed (global supermartingale) | **PROVED** |
| deterministic post-turning parametrix (Weber/Olver) | **CITED** |
| noise carried through = Prüfer phase-diffusion | **DERIVED** |
| **connection variance $O(\eta^2)$ uniform in ρ** | **DERIVED (outer) + NUMERIC (full)** |
| **inner-PC noise control, uniform in Δ** (Step 3) | **OPEN ⚑ — the single residual** |
| sub-Gaussian higher moments $C_m\eta^{2m}$ | **OPEN** (expected from phase-diffusion structure) |

**T1 is reduced to one named, well-localised estimate.** The approach is proved, re-entry closed, the
deterministic connection cited, the noise mechanism derived, and the **key uniformity (the $O(\eta^2)$
connection variance) is numerically confirmed**. What remains for *full closure* is the rigorous **inner-PC
noise control** (Step 3) — carrying the singular white-noise perturbation through Olver's PC parametrix with
constants uniform as the turning degenerates — plus the routine upgrade from variance to sub-Gaussian
moments. **The residual, if it does not fully close, is exactly that single inner-PC estimate**; everything
around it is in hand or validated.

## Net

The stochastic-Olver bound is now a sharply-posed, half-proved object: outer phase-diffusion variance
**derived and bounded uniformly**, the full $O(\eta^2)$-uniform connection variance **numerically confirmed**
(with the OU comparison shown to diverge, so the uniformity is real), and the deterministic backbone **cited**
(Olver). T1 hangs on the **single load-bearing inner-PC noise estimate** — the cleanest and most localised
statement of what closing the tube requires. That estimate — a stochastic perturbation of the
parabolic-cylinder parametrix, uniform through the turning degeneration — is the precise next theorem.
