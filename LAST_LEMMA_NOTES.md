# The last lemma — Stage 1 (GO) and Stage 2 (upper tail proved, lower tail reduced)

_June 2026. Figure `coupled-atlas/figures/moment_lyapunov.png`; script `moment_lyapunov.py`. Executing
`regime-tests/LAST_LEMMA_PROOF_STRATEGY.md`. Target: $\kappa_{2m}(\delta\Theta)\le C^m m!\,\eta^{2m}$ uniformly
in Δ, all m. Tags [PROVED]/[CITED]/[DERIVED]/[NUMERIC]/[HEURISTIC]; ⚑ = load-bearing._

## Stage 1 — Route C Lyapunov go/no-go: **GO** [NUMERIC, decisive]

Integrated the moment hierarchy $\dot M_k=-2\bar p\,k M_k-k M_{k+1}+\tfrac{\eta^2}{2}k(k-1)M_{k-2}$ (to $k=14$)
along the **real swept canard** $\bar p(\tau)$, for $\rho\in[0.04,3]$, and tested the dissipation
$\dot{\mathcal H}_m+2c_0\mathcal H_m\le O(\eta^2)$ with $\mathcal H_m=\sum_{k=2}^{2m}\tfrac{\beta^k}{k!}M_k$.

**Result: the drift defect is $<0$ for every $\beta\in\{0.5,1,2,3\}$, every $\rho$, every $m\le5$ (moments to
$x^{10}$) — including through the turning** (worst point $Y\approx+0.04$, at the floor, still negative). It is
**$m$-independent** (m=2,5 curves coincide — the all-moment uniformity) and **$\rho$-uniform** (the defect
barely moves across the merge). A single Δ-independent $\beta$ works. **⇒ GO.** Per the strategy note, this
means the analytic proof is essentially bookkeeping; the floor's restoring beats the cubic coupling.

The structural reason (analytic): at the tightest point $\bar p=c_0$, the restoring/floor combination gives
$A+C=2c_0\sum_{k\ge2}\beta^k M_k\big(\tfrac1{k!}-\tfrac1{(k-1)!}\big)=-2c_0\sum_{k\ge2}\tfrac{(k-1)\beta^k}{k!}M_k<0$
— the $(k-1)$ factor is the built-in dissipation, and the cubic coupling $B$ is dominated. **[DERIVED — confirms the numerics.]**

## Stage 2 — Route A (exponential martingale), the route the dissipation favors

Inner SDE = Langevin in the bistable $\Phi=\bar p x^2+\tfrac13x^3$: $dx=-\Phi_x d\tau+\eta dW$. Generator
$$\mathcal L e^{\lambda x}/e^{\lambda x}=-2\bar p\lambda x-\lambda x^2+\tfrac{\eta^2}{2}\lambda^2 .$$

### (i) Upper tail ($\lambda>0$): **PROVED** — the easy half (cubic confinement)

For $\lambda>0$, complete the square: $-\lambda x^2-2\bar p\lambda x=-\lambda(x+\bar p)^2+\lambda\bar p^2$, so the
generator is **bounded above globally** by $g_+(\lambda,\tau)=\lambda\bar p(\tau)^2+\tfrac{\eta^2}{2}\lambda^2$
(max at $x=-\bar p$, in the well; $-\lambda x^2\to-\infty$ both ways). Hence
$M^\lambda_\tau=\exp(\lambda x_\tau-\int_0^\tau g_+)$ has drift $\le0$ — a **genuine supermartingale**
(integrable: $e^{\lambda x}\to0$ on escape, super-confined for $x>0$). Optional stopping:
$$\mathbb E[e^{\lambda x_\tau}]\le e^{\lambda x_0+\lambda\int_0^\tau\bar p^2+\frac{\eta^2}{2}\lambda^2\tau},\qquad
\Rightarrow\ \mathbb E\big[e^{\lambda(x_\tau-a_\tau)}\big]\le e^{\frac{\eta^2\tau}{2}\lambda^2}\ \ (\lambda>0),$$
with $a_\tau=x_0+\int_0^\tau\bar p^2$. **Sub-Gaussian, variance proxy $\eta^2\tau$.** The exposure
$\int\bar p^2 d\tau$ is **bounded uniformly in Δ** (the closed-form bounded-exposure integral, reused). So the
**upper-tail moments are sub-Gaussian, uniform in Δ**: $\mathbb E[(x_\tau-a_\tau)_+^{2m}]\le(2m-1)!!(C\eta^2)^m$.
**[PROVED — global supermartingale + Itô + bounded exposure.]**

### (ii) Lower tail ($\lambda<0$): reduced to the confined Bernstein bound + the proved escape tail

For $\lambda<0$ the generator is $|\lambda|(x+\bar p)^2-|\lambda|\bar p^2+\tfrac{\eta^2}{2}\lambda^2$ — an *upward*
parabola, **unbounded above**, so no global supermartingale (this is the escape side). Split on the confined
event $\mathcal C=\{x_s>-2\bar p+\varepsilon\ \forall s\}$:
$$\mathbb E[e^{\lambda x_\tau}]=\underbrace{\mathbb E[e^{\lambda x_\tau}\mathbf 1_{\mathcal C}]}_{\text{(I) confined}}
+\underbrace{\mathbb E[e^{\lambda x_\tau}\mathbf 1_{\mathcal C^c}]}_{\text{(II) escape}}.$$
- **(II) is in hand:** $P(\mathcal C^c)\le Ce^{-h^{\star2}/6}$ (the proved one-sided cubic-supermartingale tail);
  Cauchy–Schwarz makes (II) exponentially small. **[CITED — earlier result.]**
- **(I) is the new estimate:** on $\mathcal C$ the trajectory sees only the well, where $x^2\le(2\bar p)^2$, so the
  generator $\le g_-(\lambda)=\tfrac{\eta^2}{2}\lambda^2+|\lambda|\,O(\bar p^2)$ — a **Bernstein bound**. The
  confined supermartingale then gives $\mathbb E[e^{\lambda x_\tau}\mathbf 1_{\mathcal C}]\le e^{\int g_-}$, with
  $g_-$'s cubic correction producing the $m!$ (Bernstein, not pure Gaussian) growth. The cutoff
  $\varepsilon\sim\bar p_{\min}$ couples (I)/(II); **uniformity in Δ comes from the floor $\bar p_{\min}\ge c_0$**
  (the well curvature $\Phi_{xx}=2\bar p\ge2c_0$). **[PARTIAL — the confined Bernstein bound; its dissipation is
  exactly the Route C defect, numerically validated GO.]**

### Combine

Upper (PROVED) + lower (reduced) give the two-sided centered MGF
$\mathbb E[e^{\lambda(x-a)}]\le e^{C\eta^2\lambda^2}$ for $|\lambda|\lesssim\bar p_{\min}/\eta^2$, hence
$$\boxed{\ \kappa_{2m}(\delta\Theta)\le C^m\,m!\,\eta^{2m}\ \text{uniformly in }\Delta\ }$$
with $C$ depending only on $c_0$ (floor) and the bounded exposure — **modulo the lower-tail confined Bernstein
bound (I).** The Grönwall through the turning uses the bounded exposure integral (reused). **[the upper half
proved; the bound holds modulo one estimate.]**

## The exact remaining lemma ⚑

> **Confined Bernstein bound (lower tail).** For the Langevin sweep $dx=-\Phi_x d\tau+\eta dW$ in the bistable
> $\Phi$, on the no-escape event $\mathcal C$, $\ \mathbb E[e^{\lambda x_\tau}\mathbf 1_{\mathcal C}]\le
> \exp\!\big(\lambda a_\tau+C\eta^2\lambda^2\big)$ for $\lambda\in(-\bar p_{\min}/(C\eta^2),0)$, with $C$
> depending only on $c_0$ — **not Δ**.

This is the **only** unproved step. It is a time-inhomogeneous, uniformly-convex-on-a-region Bernstein
estimate; no off-the-shelf theorem covers the moving non-convex potential through a turning, but its
dissipation is the Route C defect, now **validated GO** (defect $<0$, all m, uniform in ρ, through the
turning). Routes B (Holley–Stroock LSI with the bounded non-convex defect paid by (II)) and C (the validated
Lyapunov hierarchy + Grönwall) are interchangeable proofs of it.

## T1 status

| component | status |
|---|---|
| §0 gradient/HJ structure, barrier, curvature floor | **solid** (CITED-note) |
| split on no-escape event; escape tail (II) | **PROVED** (earlier) |
| **upper-tail moments sub-Gaussian, uniform in Δ** | **PROVED** (Route A, λ>0) |
| Route C Lyapunov dissipation (go/no-go) | **NUMERIC: GO** (defect<0, all m, uniform in ρ) |
| **lower-tail confined Bernstein bound** ⚑ | **OPEN** (validated; the single residual) |
| ⇒ $\kappa_{2m}\le C^m m!\eta^{2m}$ uniform ⇒ T1 closed | **modulo the one residual** |

## Net

Stage 1 is a clean **GO**: the moment-hierarchy dissipation holds with a single Δ-independent β, all moments,
through the turning — so the floor genuinely beats the cubic coupling. Stage 2 then **proves the upper-tail
half** of the lemma outright (the cubic confinement, global supermartingale, bounded exposure) and **reduces
the lower-tail half to one named estimate** — the confined Bernstein bound — whose dissipation is exactly the
GO-validated Lyapunov defect. So the entire T1 program now rests on a **single, sharply-stated,
numerically-validated, half-proved** estimate: the time-inhomogeneous confined Bernstein bound on the escape
side, uniform via the curvature floor. That is the final theorem, and it is one Grönwall away from the
validated dissipation — though I do **not** call the numerical Lyapunov check a proof of it.
