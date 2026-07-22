# Proof attempt — the one-sided "no early escape" tail bound (Tier-3 first sub-goal)

_June 2026. Figure `coupled-atlas/figures/no_early_escape_test.png`; script `no_early_escape_test.py`.
Companion to `UNIFORM_WEBER_TUBE_THEOREM_SCOPING.md` (Tier 3 = inner nonlinear δp² control). This is the
tractable first lemma of the load-bearing step: a **one-sided** bound, the direction in which the cubic
nonlinearity has a definite sign. Tags [PROVED]/[DERIVED]/[NUMERIC]/[CITED]/[OPEN]/[LOAD-BEARING]._

---

## 1. The precise one-sided statement

Inner Riccati SDE, $dp=(V_\Delta(Y)-p^2)\,d\tau+\eta\,dW$, $Y=Y_0-\tau$, canard $\bar p_\Delta(Y)>0$ on the
attracting branch ($Y>0$). Fluctuation $x:=p-\bar p_\Delta$ obeys, **exactly [PROVED]**,
$$dx=\underbrace{-\big(2\bar p_\Delta(\tau)\,x+x^2\big)}_{=:\,b(x,\tau)}\,d\tau+\eta\,dW,
\qquad b=-x(2\bar p_\Delta+x).$$
Bistable cubic drift: stable zero $x=0$, **unstable zero $x=-2\bar p_\Delta$** (the repelling branch =
escape barrier). Early escape = $x$ reaches $-L(\tau)$, $L(\tau)\le 2\bar p_\Delta(\tau)$, **before the
deterministic peel-off** (before $Y$ reaches the turning $0$, i.e. $\tau<T$).

> **Lemma (no early escape).** _For the barrier $L(\tau)=h\,\sigma_\Delta(Y(\tau))\le 2\bar p_\Delta$,
> with $\sigma_\Delta^2$ the (OU) fluctuation variance, and $h^\star:=\min_{0\le\tau\le T}2\bar p_\Delta/\sigma_\Delta$,_
> $$\mathbb P\big(\exists\,\tau<T:\ x(\tau)\le -L(\tau)\big)\ \le\ C\,e^{-c\,h^{\star 2}},$$
> _with **$c=\tfrac16$ a universal constant (independent of $\Delta$ and $\eta$)**, and $C=C(\Delta)$ a
> prefactor. "Uniform in $\rho=\Delta/\ell$" means: the **rate $c$ is ρ-independent outright**, and the
> bound is non-vacuous uniformly provided $\inf_\rho\bar p_{\Delta,\min}\ge c_0>0$ (so $h^\star$ is bounded
> below) and $C(\Delta)$ is bounded across the merge._

**h-scaling.** $\sigma_\Delta^2=\eta^2/(4\bar p)$ ⇒ $h^\star=4\bar p_{\min}^{3/2}/\eta\propto1/\eta$: the
small-noise limit drives $h^\star\to\infty$, and the bound is the small-noise escape rate.

**Why one-sided is the tractable sub-goal.** Downward ($x\to$ escape) the nonlinearity $-x^2$ has a
**definite sign** and the drift is *restoring on the entire pre-escape strip* $x\in(-2\bar p,0)$ (there
$b=-x(2\bar p+x)>0$). The barrier is a genuine cubic potential well, so the escape is a clean
large-deviation event. The full two-sided tube also needs the post-turning region (where $\dot{\bar p}>0$,
see §4) — deferred.

---

## 2. Candidate proof strategies, compared

| strategy | controls O(1) inner δp² **uniformly in ρ**? | obstacle |
|---|---|---|
| **(A) PC-canard exponential supermartingale / cubic scale function** ($U=\Phi$, the bistable potential) | **YES** — the full nonlinear drift sits *in the exponent* (barrier $=\tfrac43\bar p^3$ exact); rate $c=\tfrac16$ is a pure number, ρ-free | prefactor $C=e^{2\int\bar p\,d\tau}$ uniformity (exposure); needs $\dot{\bar p}\le0$ (pre-turning) |
| (B) coupling / comparison to the parabolic-cylinder (linear OU) process | **NO** — $-x^2$ makes $b\le-2\bar p x$, so SDE-comparison only gives a *reduced-rate* linear OU $a_-=2\bar p-L$; the rate **degenerates** ($a_-\to0$) at the full barrier $L\to2\bar p$ | loses the cubic enhancement; valid only for $L\ll\bar p$ |
| (C) one-sided Berglund–Gentz exponential martingale | **NO** — B–G bootstraps $x^2$ as *small* ($|x|\le h\sigma\ll\bar p$); fails when $h\sigma\sim\bar p$ (the inner O(1) regime) | exactly the merge degeneration; works outer (Tier 2) only |

**Verdict.** (A) is the strategy that handles the O(1) term *by confinement, not smallness* — it uses the
exact cubic drift in the Lyapunov exponent. (B) and (C) both linearize/perturb and therefore degenerate at
the merge — they are the very failure the build diagnosed. **Push (A).**

---

## 3. Strategy A, pushed through (the supermartingale)

**Construction [PROVED].** Let $\Phi(x,\tau)=\bar p_\Delta(\tau)x^2+\tfrac13x^3$ (the frozen bistable
potential; $\Phi_x=-b$), and
$$M_\tau=\exp\!\Big(\psi(\tau)+\tfrac{2}{\eta^2}\Phi(x_\tau,\tau)\Big),\qquad
\dot\psi(\tau)=-\big(2\bar p_\Delta(\tau)+2x_+\big),\ \ \psi(0)=0,$$
on the strip $x\in[-L(\tau),x_+]$ ($x_+$ a fixed upper level). A direct Itô computation gives the
**exact identity**
$$\frac{\mathcal L M}{M}=\dot\psi+\frac{2}{\eta^2}\,\dot{\bar p}_\Delta\,x^2+\big(2\bar p_\Delta+2x\big),$$
where the $O(\eta^{-2})$ potential terms **cancel identically** because $U=\Phi$ solves
$b\,\Phi_x+\Phi_x^2=0$ (Hamilton–Jacobi). On the pre-turning window $\dot{\bar p}_\Delta\le0$ (§4), so the
middle term is $\le0$; and on the strip $x\le x_+$, the choice of $\dot\psi$ gives
$\dot\psi+2\bar p+2x=2(x-x_+)\le0$. Hence $\mathcal L M\le0$: **$M$ is a bounded supermartingale on the
strip.** [PROVED — elementary Itô + the HJ cancellation.]

**The bound [PROVED, modulo §4 inputs].** $M_0=1$ (start $x=0$, $\Phi=0$, $\psi=0$). Optional stopping at
$\theta\wedge T$ ($\theta=$ exit from $[-L,x_+]$): $\mathbb E\,M_{\theta\wedge T}\le1$. On the escape event
$E$, $x=-L$, so $M\ge\exp(\psi(T)+\tfrac{2}{\eta^2}\Phi_\ast)$ with $\Phi_\ast=\min_\tau\Phi(-L(\tau),\tau)$.
Therefore
$$\mathbb P(E)\ \le\ \underbrace{e^{\int_0^T(2\bar p_\Delta+2x_+)\,d\tau}}_{=:C}\ \exp\!\Big(-\tfrac{2}{\eta^2}\Phi_\ast\Big).$$
At the **full barrier** $L=2\bar p_\Delta$: $\Phi(-2\bar p)=\tfrac43\bar p^3$, minimized at $\bar p_{\min}$,
so $\tfrac{2}{\eta^2}\Phi_\ast=\tfrac{8\bar p_{\min}^3}{3\eta^2}=\tfrac16 h^{\star2}$ (using
$h^{\star2}=16\bar p_{\min}^3/\eta^2$). This **is the Lemma with $c=\tfrac16$.** $\;\blacksquare$ (modulo §4)

**What this achieves.** The O(1) nonlinearity is handled **exactly**: the cubic $-x^3/3$ enters
$\Phi$ and *reduces* the barrier from the Gaussian $\bar p L^2$ to $\bar p L^2-\tfrac13L^3$ — captured in
closed form, no smallness. The **rate $c=\tfrac16$ is a universal number**, so "ρ-uniform rate" is proved
outright; all ρ-dependence is pushed into $h^\star(\rho)$ and $C(\rho)$.

**Where it is complete vs blocked.**
- **[COMPLETE]** the supermartingale property, the optional-stopping bound, and the explicit rate
  $c=\tfrac16$ at the full barrier — for the pre-turning window, conditional on the deterministic facts
  below.
- **[OPEN — G1, prefactor uniformity, LOAD-BEARING for a *uniform* $C$]** $C=e^{2\int_0^T\bar p_\Delta d\tau}$
  is $\eta$-independent (so it does not touch the exponential rate) but must be bounded **uniformly in ρ**;
  this is the "exposure" $\int_0^T\bar p_\Delta\,d\tau$ through the merge. Likely $O(1)$ (the canard's area
  over the inner window) but unproven uniformly.
- **[OPEN — G2, cited input]** $\inf_\rho\bar p_{\Delta,\min}\ge c_0>0$ — from Olver's uniform
  parabolic-cylinder asymptotics (scoping §3(ii) backbone). Makes $h^\star$ bounded below ⇒ bound
  non-vacuous uniformly.
- **[OPEN — G3, minor]** $\dot{\bar p}_\Delta\le0$ on the pre-turning window for the *exact* PC canard,
  uniform in ρ. Expected (the canard descends monotonically to the turning); needs a line.

---

## 4. The role of "pre-turning" (why this is the *one-sided* sub-goal)

$\dot{\bar p}_\Delta\le0$ holds while the canard descends toward the turning ($Y\downarrow0$). There the
helpful term $\tfrac{2}{\eta^2}\dot{\bar p}x^2\le0$ can be dropped. **Past** the turning $\dot{\bar p}>0$
and this term flips sign — the construction needs modification (a genuinely different argument). The
"early escape (before the deterministic peel-off)" event lives entirely in the pre-turning window, which is
exactly why it is the tractable first sub-goal; the two-sided/full-tube estimate must also handle the
post-turning region. [DERIVED — the sign of $\dot{\bar p}$ pins down the scope.]

---

## 5. Numerical validation

Physical early escape (crossing the repelling branch $p<-\sqrt{V_\Delta}$ in the inner window
$Y\in[0.25,2]$), sweeping $\eta\in[0.70,1.45]$ to sweep $h^\star=\min 2\bar p/\sigma$, at four
$\Delta$ spanning the merge. $N=1.2\times10^5$ paths.

| Δ | ρ-range | fitted slope $c$ ($-\ln P=c\,h^{\star2}$) | $R^2$ |
|---|---|---|---|
| 0.1 | 0.08–0.13 | **0.230** | 0.958 |
| 0.5 | 0.39–0.63 | **0.231** | 0.968 |
| 1.0 | 0.78–1.11 | **0.240** | 0.976 |
| 2.0 | 1.56–2.00 | **0.243** | 0.992 |

- **Sub-Gaussian confirmed:** $-\ln P_{\rm esc}$ is **linear in $h^{\star2}$** ($R^2=0.96$–$0.99$); a pure
  $e^{-c h^{\star2}}$ law (not $e^{-ch}$ nor $e^{-ch^3}$). **[NUMERIC]**
- **Rate uniform in ρ through the merge:** $c=0.230\to0.243$ as $\rho:0.1\to2.0$ — **relative spread 6%**.
  **[NUMERIC]**
- **Consistent with the proved bound:** measured $c\approx0.236>\tfrac16=0.167$, i.e. the truth is *rarer*
  than the supermartingale bound (as it must be — the bound is not tight because the helpful
  $\dot{\bar p}x^2$ term was dropped and the dynamic sweep further suppresses escape). The bound
  $P\le Ce^{-h^{\star2}/6}$ holds with room to spare. **[NUMERIC ✓ vs PROVED]**

_(An earlier run measuring the running-min of $x/\sigma$ over the full window was contaminated by
extreme-value statistics of the long outer pre-history and is **not** the physical escape; the
repelling-branch crossing in the inner window is the correct observable.)_

---

## 6. Status

| component | status |
|---|---|
| fluctuation cubic SDE $dx=-(2\bar p x+x^2)d\tau+\eta dW$; bistable, barrier $\tfrac43\bar p^3$ | **[PROVED]** |
| supermartingale $M=e^{\psi+2\Phi/\eta^2}$, HJ cancellation, $\mathcal LM\le0$ (pre-turning) | **[PROVED]** |
| one-sided bound $P\le Ce^{-h^{\star2}/6}$, **rate $c=\tfrac16$ ρ-independent** | **[PROVED]** (modulo G1–G3) |
| sub-Gaussian $e^{-ch^{\star2}}$, $c$ uniform in ρ (6% spread), bound non-vacuous | **[NUMERIC]** |
| prefactor $C(\rho)=e^{2\int\bar p}$ uniform in ρ (exposure) | **[OPEN — G1, load-bearing for uniform $C$]** |
| $\bar p_{\min}(\rho)\ge c_0>0$ (Olver uniform PC) | **[OPEN — G2, cited backbone]** |
| $\dot{\bar p}_\Delta\le0$ exact canard, pre-turning, uniform | **[OPEN — G3, minor]** |
| post-turning region / two-sided full tube | **[OPEN — beyond this sub-goal]** |

**Net.** The first Tier-3 sub-goal is **essentially proved**: a rigorous one-sided sub-Gaussian
no-early-escape bound, with the O(1) inner nonlinearity controlled **by the cubic scale function
(confinement), not by smallness**, and a **ρ-independent rate $c=\tfrac16$** — numerically confirmed
($c\approx0.236$, 6% spread across the merge, the proved bound holding with margin). What remains is the
prefactor's uniform-in-ρ control (the exposure integral, G1) and the two cited/minor deterministic inputs
(G2, G3). This is the concrete entry point into the load-bearing Tier-3 estimate; the post-turning and
two-sided extensions remain open.
