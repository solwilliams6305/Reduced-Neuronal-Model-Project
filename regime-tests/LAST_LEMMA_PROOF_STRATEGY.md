# Proof strategy for the last lemma
### The two-sided, all-moment nonlinear inner cumulant bound (closes T1)

**Target.** For the inner fluctuation $x=\delta p$ accumulated through the turning,
$$\kappa_{2m}(\delta\Theta)\ \le\ C^{m}\,m!\,\eta^{2m}\qquad\text{uniformly in }\Delta,\ \text{all }m.$$
(The variance $m=1$ is already proved; this is the two-sided, all-moment upgrade. Note $C^m m!$ growth is an *entire-MGF / sub-Gaussian-tailed* bound — leptokurtic core allowed, light tails required.)

---

## 0. The structural observation that makes it tractable

The exact inner SDE is **a gradient (Langevin) diffusion in the bistable potential**:
$$dx=-\big(2\bar p\,x+x^2\big)\,d\tau+\eta\,dW=-\Phi_x\,d\tau+\eta\,dW,\qquad \Phi(x,\tau)=\bar p(\tau)\,x^2+\tfrac13x^3 .$$
This is the *same* $\Phi$ whose exponential martingale $e^{(2/\eta^2)\Phi}$ closed the one-sided bound. Three facts fall straight out and organize everything:

- **Barrier.** $\Phi_x=0$ at $x=0$ (well) and $x=-2\bar p$ (saddle). Barrier height $\Delta\Phi=\tfrac43\bar p^3$, so the escape exponent is $2\Delta\Phi/\eta^2=\tfrac83\bar p^3/\eta^2$ — exactly the proved $h^{\star2}/6$. The lower side is the escape side.
- **Curvature floor.** $\Phi_{xx}=2\bar p+2x$. At the well, $\Phi_{xx}(0)=2\bar p\ge 2\bar p_{\min}(\Delta)\ge 2c_0>0$ by the **Olver floor** (cited, validated $c_0\approx0.69$). This is the uniform restoring that survives the merge — the linchpin for uniformity in $\Delta$.
- **HJ identity.** $\mathcal L\,e^{(2/\eta^2)\Phi}=e^{(2/\eta^2)\Phi}\big[(2/\eta^2)\partial_\tau\Phi+\Phi_{xx}\big]$, with $\partial_\tau\Phi=\dot{\bar p}\,x^2\le0$ pre-turning. The $O(\eta^{-2})$ terms cancel because the drift is $-\nabla\Phi$; only the $O(1)$ curvature term remains. This is the engine.

**Asymmetry of difficulty.** The two sides are *not* equally hard:
- **Upper side $x>0$:** $\Phi\sim\tfrac13x^3\to+\infty$ — super-quadratic confinement, strongest exactly where you'd worry. Upper moments are controlled *for free* by the cubic growth (the upper tail is lighter than Gaussian). Easy.
- **Lower side $x<0$ (toward the barrier):** non-convex, the leptokurtic source. This is the only real work, and it is already half-done.

---

## 1. The plan: split on the no-escape event, then sub-Gaussian-confine the rest

Write the cumulants through the decomposition
$$\mathbb E[e^{\lambda\delta\Theta}]=\underbrace{\mathbb E[e^{\lambda\delta\Theta}\,\mathbf 1_{\text{confined}}]}_{\text{(I) curvature floor}}+\underbrace{\mathbb E[e^{\lambda\delta\Theta}\,\mathbf 1_{\text{escape}}]}_{\text{(II) proved tail}},$$
where *confined* $=\{x_\tau> -2\bar p+\varepsilon\ \forall\tau\}$ (stays out of the saddle neighborhood).

**(II) is already controlled.** $P(\text{escape})\le C\,e^{-h^{\star2}/6}$ from the one-sided cubic-supermartingale bound. On a Cauchy–Schwarz split, (II) contributes $\le \|e^{\lambda\delta\Theta}\|_2\,P(\text{escape})^{1/2}$ — exponentially small, harmless for every fixed $\lambda$ in the target range. `[in hand]`

**(I) is the genuinely new estimate**, and it is a *confined* problem: on the no-escape event the trajectory only ever sees the well, where $\Phi$ is uniformly convex with curvature $\ge 2c_0$. There the diffusion is a time-inhomogeneous Langevin process in a uniformly convex potential, and sub-Gaussian cumulants are the expected output.

---

## 2. How to prove (I): three interchangeable routes, pick the cleanest

**Route A — Exponential martingale + optional stopping (most in-keeping with what worked).**
Generalize the one-sided $M=e^{(2/\eta^2)\Phi}$ to a **two-parameter tilt** $M^{\lambda}_\tau=\exp\!\big(\lambda x_\tau-\!\int_0^\tau g(\lambda,s)\,ds\big)$ and choose $g$ so $M^\lambda$ is a supermartingale on the confined event. From the SDE, $\mathcal L e^{\lambda x}=e^{\lambda x}[-\lambda(2\bar p x+x^2)+\tfrac{\eta^2}{2}\lambda^2]$; the $-\lambda x^2$ term is the obstruction. On the confined region bound $-\lambda x^2\le \lambda\,(2\bar p+\varepsilon)|x|$ (barrier cutoff) and complete the square against the $-2\bar p\lambda x$ drift using the floor. This yields $g(\lambda)\le \tfrac{\eta^2}{2}\lambda^2(1+O(\lambda/\bar p_{\min}))$, i.e. a Bernstein MGF $\le e^{C\eta^2\lambda^2}$ for $|\lambda|\lesssim \bar p_{\min}/\eta^2$ — the cumulant bound, uniform in $\Delta$ via $\bar p_{\min}\ge c_0$. The cubic correction is what produces the $m!$ (not pure Gaussian) growth.

**Route B — Bakry–Émery / log-Sobolev on the well (cleanest citation).**
On $\{x>-2\bar p+\varepsilon\}$, $\nabla^2\Phi=\Phi_{xx}\ge 2c_0>0$, so the *frozen-time* Langevin measure $e^{-2\Phi/\eta^2}$ satisfies a log-Sobolev inequality with constant $\ge 2c_0$, giving sub-Gaussian concentration with variance proxy $\eta^2/(2c_0)$ — manifestly $\Delta$-uniform. The work is the **time-inhomogeneous upgrade**: the process is a transient sweep, not stationary, so you need the moving-potential version (a Grönwall/entropy-dissipation argument with the LSI constant bounded below uniformly along the sweep). This is the one citation-adjacent route; Holley–Stroock handles the bounded non-convex defect near the barrier as a perturbation, paid for by the escape tail (II).

**Route C — Moment-hierarchy Lyapunov function (most elementary, most checkable).**
The moments obey $\frac{d}{d\tau}\mathbb E[x^k]=-2\bar p\,k\,\mathbb E[x^k]-k\,\mathbb E[x^{k+1}]+\tfrac{\eta^2}{2}k(k-1)\mathbb E[x^{k-2}]$. The $-k\,\mathbb E[x^{k+1}]$ coupling is the non-closure (leptokurtosis). Close it from above with the Lyapunov functional $\mathcal H_m(\tau)=\sum_{k\le 2m} a_k\mathbb E[x^k]$, $a_k=\beta^k/k!$, and show $\dot{\mathcal H}_m\le -2c_0\,\mathcal H_m+(\text{source }\sim\eta^2)$ on the confined event — the floor gives the negative feedback, the cubic term is dominated by the $\beta$-geometric weights. Grönwall then yields $\mathbb E[x^{2m}]\le C^m m!\,\eta^{2m}$ directly. This route is the easiest to **validate numerically step-by-step** (the $a_k$ and $\beta$ are explicit) before committing to the analysis.

All three lean on the **same two ingredients**: the curvature floor $\bar p_{\min}\ge c_0$ (uniformity) and the escape tail (the non-convex defect). That convergence is the sign the lemma is really one estimate.

---

## 3. The single load-bearing new step, named

Everything cited or already-proved reduces the lemma to:

> **Uniform confined cumulant bound.** For the Langevin sweep $dx=-\Phi_x d\tau+\eta dW$ in the bistable $\Phi$, *conditioned on no barrier escape*, the cumulants satisfy $\kappa_{2m}\le C^m m!\,\eta^{2m}$ with $C$ depending only on $c_0$ (the curvature floor) — **not** on $\Delta$.

This is a time-inhomogeneous, uniformly-convex-on-a-region concentration statement. It is *new* (no off-the-shelf theorem covers the moving non-convex potential with a uniform floor through a turning), but it is **not** a new black box: it is the all-moment, two-sided sibling of the one-sided bound already proved, and Routes A/C make the cubic correction explicit.

---

## 4. Risks / where it could resist

- **Time-inhomogeneity through the turning.** The LSI/martingale constants must be bounded *uniformly along the sweep*, including the instant of the turning. The floor $c_0$ is what should guarantee this, but the rigorous Grönwall needs $\int \bar p_{\min}(\tau)\,d\tau$ controlled — the same exposure integral already shown bounded (the closed-form $\tfrac1\Delta\ln(\cdots)$). Reuse it.
- **The barrier cutoff $\varepsilon$.** Choosing $\varepsilon$ couples (I) and (II): smaller $\varepsilon$ tightens the convexity floor but loosens the escape tail. Optimize $\varepsilon\sim\bar p_{\min}$.
- **$m$-uniformity of $C$.** Route C's $\beta$ must be chosen $m$-independent; check the geometric weights actually dominate the cubic coupling for all $m$ (a finite computation per $k$).

---

## 5. Concrete first move (cheap, decisive)

Run **Route C numerically as a finite hierarchy**: integrate the moment ODEs up to $k=12$ along the real swept $\bar p(\tau)$, form $\mathcal H_m$ with trial $\beta$, and check $\dot{\mathcal H}_m\le -2c_0\mathcal H_m+O(\eta^2)$ holds pointwise across $\rho\in[0.04,3]$. If the Lyapunov inequality holds numerically with a $\Delta$-independent $\beta$, the analytic proof is essentially bookkeeping; if it fails at the turning, it localizes the exact obstruction. Either outcome is decisive and costs one script.

---

**Honest status.** §0 (gradient structure, floor, HJ identity) is solid. The split §1 is rigorous given the proved escape tail. The three routes §2 are strategies, not theorems — Route A is most self-contained, Route C most checkable, Route B most citable. §3 is the one genuinely new estimate, and it inherits the curvature floor and exposure bound already in hand. Tag everything on execution; do not call the numerical Lyapunov check a proof of §3.
