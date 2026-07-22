# Scoping note — the final theorem: the uniform parabolic-cylinder tube estimate with noise

_June 2026. Scoping only (statement + proof plan), **not** a proof. Companion to
`WEBER_TUBE_BUILD_NOTES.md` (the validated covariance splice) and
`regime-tests/RH_DIRECTION_NOVEL_ANGLES.md` (§0/§1). This is the single open object that closes **T1**
(the uniform tube through the merge) and whose connection data is the new **T2** cusp edge law._

Tags used throughout: **[PROVED]** (rigorous, in hand or immediate); **[CITED]** (in the literature,
applies here after a stated reduction); **[HEURISTIC]** (supported by the build's numerics or formal
asymptotics, not proved); **[OPEN]** (the actual gap); **[LOAD-BEARING]** (an unproved step the whole
theorem rests on).

---

## 0. Setup and notation

After the Kristiansen–Pedersen blow-up of the coupled-FHN antisymmetric canard, the central
(parabolic-cylinder) chart carries the **inner Riccati SDE**

$$dp \;=\; \big(V_\Delta(Y)-p^2\big)\,d\tau \;+\; \eta\,dW_\tau,\qquad Y=Y_0-\tau,\quad \tau\in[0,T],$$

with the swept double-fold potential
$$V_\Delta(Y)=\operatorname{sign}(Y)\,|Y|\,\big(|Y|+\Delta\big),\qquad \Delta=\Delta(g)=2\sqrt{-2g/3}\ \ (\text{branch separation}).$$

- $\Delta=0$ ⇒ pure cusp $V_0=\operatorname{sign}(Y)Y^2$; $\Delta>0$ ⇒ two simple folds at $Y=0,-\Delta$.
- $\eta$ = effective inner noise (the rescaled $\sigma$); escape = **first explosion** $p\to-\infty$ (first node of $u$, Cole–Hopf $p=u'/u$, $u''=(V_\Delta-\eta\dot W)u$).
- **Canard** $\bar p_\Delta(Y)$: the deterministic ($\eta=0$) solution tracking the attracting branch $+\sqrt{V_\Delta}$, finite and $>0$ through the turning. Existence/properties **[CITED — K–P + deterministic Olver/Weber connection]**.
- **Fluctuation** $\delta p=p-\bar p_\Delta$; its linear part solves $d\,\delta p=-2\bar p_\Delta\,\delta p\,d\tau+\eta\,dW$, with variance $v_\Delta$ from $dv/d\tau=-4\bar p_\Delta v+\eta^2$ **[PROVED]**.
- **Inner scale** $\ell=\ell(\eta)$: defined by the breakdown balance $v_{\rm qs}(\ell)\sim V_\Delta(\ell)$ (linearized tube reaches the branch separation), giving $\ell\sim\eta^{2/3}$ for the cusp **[DERIVED; exact constant deferred to §3(ii)]**. The uniformity parameter is the dimensionless ratio $\boxed{\rho:=\Delta/\ell}\in[0,\infty)$.
- **Spliced envelope** $\sigma_\Delta(Y)^2:=v_\Delta(Y)$: $v_{\rm qs}=\eta^2/4\sqrt{V_\Delta}$ in the outer zone $|Y|\gg\ell$, the **parabolic-cylinder** variance $v_{\rm Weber}$ in the inner zone $|Y|\lesssim\ell$, matched in the overlap. (This is the object validated to 15–27% in the build.)
- **Peel-off (exit) time** $\tau^\star=\inf\{\tau:\,p(\tau)\le -M\}$ for a fixed large $M$; the **escape level** is $Y^\star=Y(\tau^\star)$, suitably rescaled.

---

## 1. Theorem statement (TARGET — OPEN)

> **Theorem (uniform parabolic-cylinder tube with noise).** _There exist constants $h_0,c,C>0$ and
> $\eta_0>0$, **all independent of $\Delta$ and $\eta$**, such that for every $\eta\in(0,\eta_0]$,
> every $\Delta\ge 0$, and every $h\in[1,h_0]$:_
>
> **(C1) Uniform confinement.**
> $$\mathbb P\!\left(\ \sup_{0\le\tau<\tau^\star}\ \frac{|\delta p(\tau)|}{\sigma_\Delta(Y(\tau))}\ \ge\ h\ \right)\ \le\ C\,(1+T)\,e^{-c\,h^2},$$
> _and the constants $c,C$ are **uniform in the crossover** $\rho=\Delta/\ell\in[0,\infty)$._
>
> **(C2) Endpoint laws.** _The rescaled escape level $Y^\star$ converges in distribution to:_
> - _$\rho\to\infty$ (noisy-Airy): the Tracy–Widom law $\mathrm{TW}_\beta$, $\beta=4/\eta^2$;_
> - _$\rho\to0$ (noisy-Weber): a law $\mathcal W_\beta$ given by the **connection (Stokes) data of the parabolic-cylinder parametrix** $\bar u$ — the candidate **new cusp edge law**._
>
> **(C3) Uniform crossover.** _The family $\{\mathcal L(Y^\star\mid\rho)\}_{\rho\in[0,\infty]}$ is tight and the convergence in (C2) is uniform on compact $\rho$-ranges; the limit interpolates $\mathrm{TW}_\beta\leftrightarrow\mathcal W_\beta$ monotonically in the cumulant fingerprint._

**What it buys.** (C1) is the **T1** statement (a Berglund–Gentz tube that does **not** degenerate at the
merge, because the comparison sd is the spliced $\sigma_\Delta$, finite through the turning, rather than
the divergent $v_{\rm qs}$). (C2)(ρ→0) **is** the construction of **T2**: the new cusp law is *defined*
as the connection data of $\bar u$, and the theorem asserts it is the actual escape limit. (C3) is the
"uniform in $\Delta/\ell$ through the merge" clause that the falsifier campaign isolated as the crux.

**Honest caveat on the statement.** The exponent in $\ell\sim\eta^{2/3}$ and the precise normalization of
$Y^\star$ are part of §3(ii); the **shape** of the theorem (uniform sub-Gaussian confinement around the
parabolic-cylinder envelope + connection-data limit) is the robust content.

---

## 2. Proof skeleton

Five steps; the OU→Weber comparison enters **rigorously at Step B**.

**Step A — Reduction (in hand). [CITED]** Blow-up to the inner Riccati and existence of the canard
$\bar p_\Delta$ with the stated outer/inner asymptotics. K–P supply the deterministic geometry; the only
new bookkeeping is carrying the additive noise $\eta\,dW$ through the rescaling (linear in the chart, so
clean).

**Step B — The linear propagator via the parabolic-cylinder parametrix. [CITED backbone + NEW glue]**
The linearized flow has propagator $\Phi(\tau,s)=\exp\!\big(-2\!\int_s^\tau \bar p_\Delta\big)$, and
$\bar p_\Delta=\bar u'/\bar u$ with $\bar u$ the **parabolic-cylinder canard**. The covariance is
$C_{\rm lin}(\tau,\tau')=\eta^2\!\int_0^{\min}\Phi(\tau,s)\Phi(\tau',s)\,ds$, whose Gaussian envelope is
exactly the spliced $\sigma_\Delta$. **This is where OU→Weber is made rigorous:** the *outer* asymptotic
of $\bar u$ ($|Y|\gg\ell$) is the WKB/OU form $\bar p_\Delta\approx\sqrt{V_\Delta}$ (Airy at a single
fold); the *inner* form ($|Y|\lesssim\ell$) is the Weber form; **Olver's uniform two-turning-point
asymptotics glue them with explicit error bounds**, uniformly as the folds coalesce. The covariance
splice validated in the build (rel-err 27%) is the numerical shadow of this step.

**Step C — Linear (Gaussian) tube. [NEW, tractable — Tier 1]** Given Step B, a Berglund–Gentz / Doob
maximal-inequality argument yields (C1) **for the linear process** $\delta p_{\rm lin}$, with constants
controlled by $\sup\Phi$ and the envelope — both uniform in $\rho$ *iff* Step B's bounds are uniform.

**Step D — Nonlinear closure. [NEW, LOAD-BEARING — Tier 3]** Reinstate the $-\delta p^2$ term. A
stopping-time bootstrap shows that, up to $\tau^\star$, $\delta p$ stays within $h\,\sigma_\Delta$ and the
nonlinear term is dominated — **provided** the inner-region estimate of §3(i) holds uniformly in $\rho$.
This is the genuinely hard step (the $\sim$20% residual lives here).

**Step E — Exit law = connection data. [NEW — Tier 4]** On the tube, the first explosion is governed by
the parabolic-cylinder parametrix's recessive→oscillatory connection; the noisy connection coefficient's
law is $\mathcal W_\beta$ (ρ→0) and $\mathrm{TW}_\beta$ (ρ→∞). Assemble A–E with uniform constants ⇒ (C1)–(C3).

---

## 3. The two sub-pieces, scoped separately

### (i) Control of the nonlinear $\delta p^2$ correction — the ~20% residual

**What is needed.** A bound that, on the event $\{|\delta p|\le h\sigma_\Delta\}$ up to $\tau^\star$, the
accumulated nonlinear drift $\int \delta p^2$ does not push the path out of the tube before the genuine
peel-off, **with constants uniform in $\rho$**.

**Is it a perturbation of the linear tube?** **In the outer zone, yes; in the inner zone, no.**
**[LOAD-BEARING / HEURISTIC]**
- *Outer* ($|Y|\gg\ell$): $\bar p_\Delta\sim\sqrt{V_\Delta}\gg\sigma_\Delta$, so $\delta p^2$ is smaller
  than the linear restoring $2\bar p_\Delta\delta p$ by the factor $h\sigma_\Delta/\bar p_\Delta\ll1$ —
  a standard bootstrap perturbation (Tier 2). **[CITED-style, expected PROVED]**
- *Inner* ($|Y|\lesssim\ell$): $\bar p_\Delta=O(1)$ (the Weber value) and $\sigma_\Delta=O(1)$, so
  $\delta p^2$ is **comparable** to the linear term — *not* a perturbation. Worse, $\delta p^2$ is exactly
  the term that *drives* the Riccati explosion, so the escape itself cannot be treated perturbatively.
  The nonlinear correction must be controlled **using the parabolic-cylinder structure itself** (the
  confinement of $\bar u$), not by smallness. **This is its own argument, and it is the crux.**

**Verdict.** Sub-piece (i) needs a dedicated inner estimate (a self-consistent / Lyapunov argument on the
parabolic-cylinder operator), **not** a perturbation of the linear tube. The build's 15–27% gap is the
quantitative signature that this term is first-order, not negligible.

### (ii) Parabolic-cylinder connection-coefficient bounds, uniform in $\Delta/\ell$

**What uniformity is required.** Error bounds on the recessive→dominant connection of $u''=(V_\Delta-\cdots)u$
that are **uniform as the two turning points coalesce** ($\Delta\to0$, i.e. $\rho\to0$) **and** down to
the inner scale, so that Step B's propagator constants do not blow up at the merge.

**Does Olver supply it?** **Partly — deterministic backbone yes, stochastic uniformity no.**
- **[CITED — supplies]** Olver's uniform asymptotics for a 2nd-order ODE with two coalescing turning
  points (the parabolic-cylinder comparison, DLMF §"two coalescing turning points") give exactly this
  uniformity **for the deterministic equation**, with computable error bounds in the coalescence
  parameter. After a reduction of $V_\Delta$ to Olver's normal form $\tfrac14 W^2-a$ (a smooth change of
  variable; **routine but must be checked uniform in $\Delta$ [OPEN-minor]**), this is in hand.
- **[OPEN — gap]** Two upgrades are *not* in Olver: **(a)** the additive **noise** $\eta\,dW$ (the
  connection coefficients become random; one needs the *distribution/moment* bounds of the noisy
  connection, uniform in $\rho$); **(b)** propagating the deterministic error bounds through the
  **covariance** $C_{\rm lin}$ uniformly. (a) is the genuinely new analytic content; (b) is technical
  but expected tractable once (a) is set.

**Verdict.** Olver gives the deterministic uniform parabolic-cylinder connection (the backbone of Step B);
the **stochastic** uniform connection (noisy coefficients + covariance propagation, uniform in $\rho$) is
new and is the second load-bearing gap.

---

## 4. Citable vs genuinely new

**Citable (machinery that applies after a stated reduction):**
- **Berglund–Gentz** sample-path tube estimates for slow–fast SDEs near uniformly-hyperbolic slow
  manifolds, and the single-fold (Airy) case — gives Steps C and the *outer* part of D. **[CITED]**
- **Olver** uniform asymptotics for two coalescing turning points (parabolic-cylinder comparison, error
  bounds) — the deterministic backbone of Step B and sub-piece (ii). **[CITED]**
- **Kristiansen–Pedersen** deterministic cusp/folded-node blow-up — Step A geometry, canard existence.
  **[CITED]**
- **Ramírez–Rider–Virág** (stochastic Airy ⇒ $\mathrm{TW}_\beta$) — the $\rho\to\infty$ endpoint (C2). **[CITED]**

**Genuinely new technical content:**
1. **The stochastic parabolic-cylinder comparison** — the noisy Weber tube and its covariance splice as a
   *uniform* object (Step B's stochastic glue + the covariance propagation, §3(ii)(a–b)).
2. **The inner nonlinear $\delta p^2$ control, uniform in $\rho$** (§3(i), inner zone) — where the cusp
   escape is first-order, not a perturbation of the linear tube. **The crux.**
3. **The uniform-in-$\Delta/\ell$ assembly** — gluing single-fold (Airy/B–G) and merged (Weber) regimes
   into one estimate with $\rho$-independent constants (C3).
4. **The identification** $\mathcal W_\beta$ = law of the noisy parabolic-cylinder connection coefficient
   (C2, ρ→0) — the *construction* of the new T2 cusp edge law.

---

## 5. Obstacles and a tiered attack plan (most tractable first)

**Tier 1 — Linear uniform tube (most tractable). [expected PROVED]**
Prove (C1) for $\delta p_{\rm lin}$ with $\rho$-uniform constants, using B–G + Olver's deterministic
uniform connection for the propagator. *Obstacle:* the reduction of $V_\Delta$ to Olver normal form
uniform in $\Delta$. *Payoff:* the Gaussian skeleton of the whole theorem; directly cross-checkable
against the validated covariance splice. **Do this first.**

**Tier 2 — Outer nonlinear closure. [expected PROVED]**
Bootstrap the $\delta p^2$ term in $|Y|\gg\ell$ (perturbative, since $\bar p_\Delta$ large). *Obstacle:*
none beyond standard B–G; mainly bookkeeping uniform in $\Delta$.

**Tier 3 — Inner nonlinear control (LOAD-BEARING). [OPEN]**
The dedicated parabolic-cylinder Lyapunov/self-consistent estimate for $\delta p^2$ in $|Y|\lesssim\ell$,
uniform in $\rho$ (§3(i) inner). *Obstacle:* $\delta p^2$ is first-order and drives the escape — needs the
confinement of $\bar u$ itself. **This is the real theorem; everything else is scaffolding.** Tractable
sub-goal: first prove a *one-sided* (upper) tail bound (escape not too early), which suffices for the
peel-off law's right tail, before the full two-sided tube.

**Tier 4 — Stochastic connection + exit law. [OPEN]**
The noisy parabolic-cylinder connection coefficients (moments uniform in $\rho$) and the identification of
their law with $\mathcal W_\beta$ (Step E, C2 ρ→0). *Obstacle:* §3(ii)(a) — no cited stochastic version;
this defines T2. Tractable sub-goal: compute the *first two cumulants* of the noisy connection
perturbatively in $\eta$ and check against the cusp fingerprint (skew $+0.61$, exkurt $-0.24$).

**Suggested order:** Tier 1 → Tier 2 → Tier 4-cumulants (cheap, high-information cross-check) → Tier 3
(the hard core) → Tier 4-full. Tiers 1–2 are largely assembly of cited machinery; Tiers 3 and 4(a) carry
the genuinely new content.

---

## 6. Status summary

| component | role | status |
|---|---|---|
| inner Riccati + canard + fluctuation/variance ODE | Step A | **[PROVED]** / K–P **[CITED]** |
| spliced envelope $\sigma_\Delta$ finite through merge | comparison object | **[DERIVED]** + **[NUMERIC]** (build, 15–27%) |
| linear propagator via PC parametrix (OU↔Weber glue) | Step B | backbone **[CITED Olver]**, stochastic glue **[OPEN]** |
| linear uniform tube (C1 for $\delta p_{\rm lin}$) | Tier 1 | **[expected PROVED]** (B–G + Olver) |
| outer nonlinear closure | Tier 2 | **[expected PROVED]** |
| **inner nonlinear $\delta p^2$ control, uniform in $\rho$** | Tier 3 | **[OPEN, LOAD-BEARING]** |
| stochastic PC connection bounds uniform in $\rho$ | §3(ii)(a), Tier 4 | **[OPEN]** |
| endpoint laws $\mathrm{TW}_\beta$ (ρ→∞), $\mathcal W_\beta$ (ρ→0) | C2 | ρ→∞ **[CITED RRV]**; ρ→0 = new, **[OPEN]** = def. of T2 |
| full theorem (C1–C3) | — | **[OPEN]** — reduced to Tiers 3 & 4(a) |

**One-line scope.** Tiers 1–2 assemble cited Berglund–Gentz + Olver machinery; the theorem then rests on
exactly **two** new estimates — the **inner nonlinear $\delta p^2$ control** (§3(i)) and the **stochastic
parabolic-cylinder connection bounds** (§3(ii)(a)), both **uniform in $\Delta/\ell$**. The first closes
T1; the second defines and delivers the new T2 cusp law.
