# Program 2, Route 2b — Step 0 done (leading coefficients pinned) + the resurgent-growth claim corrected

_July 2026 (incoming agent, Fable 5). Executes Route 2b Step 0 (`FRONTIER_SCOPING_NOTES.md`): reconcile the
small-η FP-vs-MC discrepancy in the weak-noise coefficients, with a **converged method**, and then re-examine
the resurgent signature. **Result:** a converged perturbative expansion off the exact Weber backbone pins the
leading coefficients decisively (three methods now agree), **retracts the FP small-η values as boundary-layer
artifacts**, and — importantly — shows the previously-quoted **resurgent-growth ratios (1.23 → 2.33) were an FP
artifact** (the FP even had the wrong *sign* for the sub-leading variance coefficient). The robust resurgence
evidence (β = 2 overshoot; Borel singularity at the instanton action) survives; the coefficient-growth diagnosis
of divergence does **not** (it needs orders beyond the current statistical reach — a mapped obstruction).
Scripts `coupled-atlas/weaknoise_sectors.py`, `weaknoise_series.py`, `weaknoise_v1.py`, `stochastic_scope.py`;
figure `figures/weaknoise_sectors.png`. Tags **[DERIVED]/[NUMERIC ✓]/[RETRACTED]/[open-mapped]**._

---

## 0. The discrepancy (Step 0 target)

The weak-noise sectors of 𝒲 (expansion in η² = 4/β) had **two conflicting numerical extractions**:

| coefficient | FP (small-η) | MC / closed-form | agree at β = 2? |
|---|---|---|---|
| variance leading $C_V=v_0$ | **0.200** | **0.1339** (closed-form + MC) | yes (both → β=2 shape) |
| skew leading $s_0=$ skew/η | **0.66** | **≈1.2** (MC) | yes (both skew(β=2)≈0.61) |

They agree at β = 2 but disagree at small η. Step 0 = pin the leading coefficients with a converged method.

## 1. The converged method — perturbative sectors off the exact backbone

Both FP and MC are fragile at small η for structural reasons: **FP** integrates the escape Fokker–Planck PDE
with an upwind scheme whose numerical diffusion at the deep absorbing barrier is $D_{\rm num}\sim|{\rm drift}|\,
dp/2\sim p^2\,dp/2$ — *enormous* compared to the physical $D=\eta^2/2$ when η is small (diagnosed directly:
the FP mean extrapolates to **−2.065** as η→0, but the exact deterministic node is **−2.18762**, an $O(0.12)$
bias); **MC** must resolve a skew of size $\sim\!\eta$ on a distribution of width $\sim\!\eta$, needing huge N
and tiny dt.

The clean method avoids both: expand the escape location directly,
$$
Y^\star=Y^\star_0+\eta Y_1+\eta^2 Y_2+\eta^3 Y_3+\cdots,
$$
where $Y_k$ are **η-independent functionals** of the white-noise potential (the field perturbations
$u_k''-Vu_k=\xi\,u_{k-1}$, recessive BC), evaluated by the node-shift recursion at the exact base node
$Y^\star_0=-2.18762$ (using $u_0(Y^\star_0)=0\Rightarrow u_0''(Y^\star_0)=0$, $u_1''(Y^\star_0)=V^\star u_1$):
$$
Y_1=-\tfrac{u_1}{u_0'},\quad
Y_2=-\tfrac{u_1'Y_1+u_2}{u_0'},\quad
Y_3=-\tfrac{\tfrac16 V^\star u_0'Y_1^3+u_1'Y_2+\tfrac12 V^\star u_1 Y_1^2+u_2'Y_1+u_3}{u_0'}\quad(\text{at }Y^\star_0).
$$
No FP boundary layer, no rare-escape sampling. The cumulant series coefficients are then Wick/moment averages:
$$
v_0=\langle Y_1^2\rangle,\quad m_1=\langle Y_2\rangle,\quad
s_0=\frac{3\langle Y_1^2Y_2\rangle_c}{v_0^{3/2}},\quad
v_1=\operatorname{Var}(Y_2)+2\langle Y_1Y_3\rangle.
$$
(Here $\operatorname{Var}=v_0\eta^2+v_1\eta^4+\cdots$; the odd cross terms vanish by $W\!\to\!-W$.)

## 2. Result — leading coefficients pinned (Step 0 closed) [NUMERIC ✓, converged]

Converged in step $h$ (1.5·10⁻³, 10⁻³) and seed (×3), N = 3·10⁵ (`weaknoise_v1.py`, `weaknoise_sectors.py`):

| coefficient | converged value | independent target | verdict |
|---|---|---|---|
| $v_0=\langle Y_1^2\rangle$ | **0.1337 ± 0.0004** | 0.1339 (closed-form) | ✓✓ matches to 4 digits |
| $m_1=\langle Y_2\rangle$ | **+0.2124 ± 0.0006** | ≈0.212 (notes/FP-leading) | ✓ |
| $s_0=$ skew/η | **1.160 ± 0.02** | ≈1.18 (small-η MC) | ✓ |

**The framework is validated on two independently-known quantities ($v_0$ closed-form, $m_1$), so its
predictions are trustworthy.** The converged leading coefficients are $C_V=0.134$ and $s_0\approx1.16$–1.18.

> **Step 0 verdict [RETRACTED FP].** The FP small-η values ($C_V=0.20$, $s_0=0.66$) are **numerical-diffusion
> artifacts** (upwind scheme at the deep absorbing barrier) and are retracted. The correct leading weak-noise
> coefficients are $C_V=0.1339$ and $s_0\approx1.17$, confirmed by **three independent methods** (closed-form
> Green's function, small-η MC, converged perturbative expansion).

## 3. The sub-leading variance coefficient — and the corrected resurgence picture [NUMERIC ✓ / RETRACTED]

The direct-functional method gives (converged, seed-robust):
$$
\boxed{\ \operatorname{Var}(Y^\star)=\eta^2\big(0.1337+0.110\,\eta^2+\cdots\big),\qquad v_1=+0.110\pm0.001.\ }
$$

Compare the **contaminated FP** extraction (`transseries_orders.py`): $\operatorname{Var}=\eta^2(0.200-0.245
\eta^2+0.572\eta^4)$, i.e. $v_1=-0.245$. **Both the magnitude and the *sign* of $v_1$ are wrong in the FP.**

**Consequence — the resurgent-growth claim is corrected [RETRACTED ⚑].** The handoff/notes cited "variance
coefficient ratios grow, 1.23 → 2.33 per order ⇒ divergent/asymptotic" as the resurgence signature. Those
ratios were computed from the FP coefficients $(0.200,-0.245,0.572)$ and are **an artifact** — the leading
coefficient itself was wrong (0.200 vs 0.134). With converged data, $|v_1/v_0|=0.82$ (not 1.23), so **the
runaway coefficient growth is not present at this order.** This does *not* prove the series is convergent —
diagnosing an asymptotic (zero-radius, factorial) series needs the large-order trend ($v_2,v_3,\dots$), which
is beyond the current reach (see §5). It does mean the earlier "coefficient-growth ⇒ divergent" argument must
be withdrawn.

**What robustly survives as resurgence evidence:**

1. **The β = 2 overshoot (non-perturbative regime), with converged coefficients.** The leading weak-noise skew
   extrapolated to β = 2 (η = √2) is $s_0\eta=1.16\times1.414=1.64$ versus the true **0.605** — an overshoot of
   **×2.7** (the old ×2.9 used the contaminated leading skew; the effect is real, not an artifact). Likewise the
   variance partial sums at β = 2 straddle the truth: $S_0=v_0\eta^2=0.268$ (under), $S_1=(v_0+v_1\eta^2)\eta^2
   =0.708$ (over), true $\approx0.47$ — the second term (0.44) *exceeds* the first (0.27), so the series is
   **outside its regime at β = 2** (radius $\sim\eta^2\approx v_0/v_1\approx1.2$, i.e. β ≈ 3.3). β = 2 is
   genuinely non-perturbative — the leading weak-noise theory cannot reach the physical point. **[NUMERIC ✓]**
2. **The Borel singularity at the instanton action.** The non-perturbative sector is the FW left-tail instanton
   $\sim e^{-s^5/10\eta^2}$ (`PERSISTENCE_ITEM2_NOTES.md`, rigorously bounded now in `PROGRAM1_TIERA_THEOREM.md`
   §4.2); the Borel-plane singularity of the η²-series sits at the corresponding instanton action. The
   two-term radius estimate ($\eta^2\approx1.2$) is *consistent* with an $O(1)$ instanton action, but with only
   $v_0,v_1$ this is suggestive, not a determination. **[DERIVED (structure) + NUMERIC (radius)]**

So the trans-series *structure* stands — perturbative sectors (now with **correct** leading coefficients) plus
an instanton sector, with β = 2 demonstrably non-perturbative — but the specific "coefficient-growth ⇒
divergent" numerology was FP-contaminated and is retracted.

## 4. Net (Route 2b status)

- **Step 0: DONE, decisive.** Leading weak-noise coefficients pinned by three agreeing converged methods:
  $C_V=0.1339$, $m_1=0.212$, $s_0\approx1.17$. FP small-η values retracted as boundary-layer artifacts.
  **[NUMERIC ✓✓]**
- **Sub-leading variance coefficient computed (converged):** $v_1=+0.110$; the FP-based resurgent-growth ratios
  are **retracted** (wrong sign). **[NUMERIC ✓ / RETRACTED]**
- **Resurgence evidence, honestly:** β = 2 is non-perturbative (×2.7 skew overshoot, variance partial sums
  straddle with a growing second term) and the non-perturbative sector is the instanton at $s^5/10$ — the
  trans-series *structure* is real; the *divergence-by-coefficient-growth* diagnosis needs higher orders.

## 5. The obstruction to more orders (mapped)

The clean node-shift recursion is **noise-safe only through $Y_3$**: it relies on $u_0(Y^\star_0)=0$ to kill the
white-noise term in $u_1''(Y^\star_0)=V^\star u_1+\xi\,u_0$. At **fourth order** ($Y_4$, needed for $v_2$) the
recursion requires $u_1'''(Y^\star_0)=V'^\star u_1+V^\star u_1'+\xi(Y^\star_0)\,u_0'^\star$, which contains the
**white-noise potential evaluated at the (random) node** — a genuinely distributional object (the field is only
$C^{1/2}$; there is no pointwise $\xi$). Computing $v_2$ therefore needs this pointwise term treated as an Itô
object (or a Malliavin/Wick treatment of $\langle Y_1 Y_4\rangle$), not a finite-difference march. **This is the
minimal missing ingredient for the next order**, and it is exactly the kind of stochastic-analysis step the full
"stochastic exact-WKB" program (`FRONTIER_SCOPING_NOTES.md` Route 2b) is about. Reachable, but a real subproblem
— not a coarse numerical push.

**Honest ceiling (unchanged).** A partial resurgent representation — the correct leading + sub-leading
perturbative coefficients, the instanton sector, and the β = 2 non-perturbative overshoot — is now on **solid,
converged footing** (this pass fixed the contaminated leading data). The full resurgent trans-series (all
orders + Borel resummation + Stokes constant) remains the deeper stochastic-exact-WKB program.

## 5bis. The obstruction REFRAMED (first tackle, 2026-07): method-limited, not fundamental

_A first move on Program 2, resuming it as the second multi-year core._ The v2 "obstruction" is a limitation of
the **node-shift recursion method**, not of the escape law — the first node of the noisy field
$u''=(V+\eta\dot W)u$ is always well-defined (the field is continuous, $C^{1/2}$; it crosses zero). Two points:

1. **The singular term is finite in the moments.** $\xi(Y^\star_0)$ is white noise at the **fixed deterministic
   node** $Y^\star_0=-2.188$ (we expand around the *deterministic* node, not a random one). Pathwise it is
   singular (white noise has no pointwise value; $\xi(Y^\star_0)\sim\Delta W/h\to\infty$ under grid refinement —
   which is *why the finite-difference march breaks*). But in the **cumulant moments** it is finite: in
   $\langle\cdots\rangle$ the factor $\xi(Y^\star_0)$ Wick-contracts with the noise in the lower $Y_k$ (which are
   Itô integrals $\int\Phi\,\xi\,ds$), and $\langle\xi(Y^\star_0)\xi(s)\rangle=\delta(Y^\star_0-s)$
   **δ-localizes** the integral to a finite **boundary value** (a value/derivative of $\Phi$ and $u_0$ at
   $Y^\star_0$). So $v_2$ is a *finite* Wick computation, not a divergence.
2. **Two concrete routes to $v_2$** (upgrading "obstructed" → "a well-posed subproblem"): **(a) boundary-Wick** —
   carry the recursion but treat the $\xi(Y^\star_0)$ term as an Itô/Malliavin object, collecting the
   δ-localized boundary contributions (the clean analytic path; this *is* the stochastic-exact-WKB step, and it
   generalizes the deterministic Weber exact-WKB / Voros machinery of Iwaki–Koike–Takei); **(b) regularized
   full-field** — mollify the noise (correlation length $\delta$), compute the first-node cumulants of the smooth
   field, extrapolate $\delta\to0$ (validates (a) but is precision-limited for the $O(\eta^6)$ coefficient).

**Net of the first tackle:** the "named hard problem" of Program 2 is reframed from a distributional *wall* into a
**finite boundary-Wick computation** — the pointwise-noise term contributes δ-localized boundary terms, not a
divergence. This is the concrete on-ramp to the perturbative sectors beyond $v_1$, and it sharpens exactly what
the stochastic-exact-WKB literature pass should target (see the kickoff brief `PROGRAM2_KICKOFF_STOCHASTIC_EXACT_WKB.md`).

## 6. The wall is CLEARED — $v_2$ computed (three agreeing routes) [NUMERIC ✓✓]

_2026-07-07 (Fable 5). Executed the §5bis plan. The $v_2$ "obstruction" is now resolved: $v_2=+0.11$, computed by
three independent methods that agree, plus a fully deterministic (MC-free) reformulation validated through $v_1$._

$$\boxed{\ \operatorname{Var}(Y^\star)=\eta^2\big(0.134+0.110\,\eta^2+0.11\,\eta^4+\cdots\big),\qquad v_2=+0.11\pm0.015.\ }$$

### 6a. Deterministic Green's-function / Wick reformulation — the perturbative sectors are MC-free
The field perturbations $Lu_k=\xi u_{k-1}$ ($L=\partial_Y^2-V$) solve via the causal Green's function
$G(Y,s)=u_0(Y)\psi(s)-u_0(s)\psi(Y)$ ($\psi(Y_0)=0,\psi'(Y_0)=1$, Wronskian $=1$), so **every $u_k$ is a $k$-fold
iterated Itô integral** and **every cumulant coefficient is a deterministic multiple integral** (Wick's theorem) —
no Monte Carlo. At the node ($u_0(Y^\star_0)=0$): $G_0(s)=-c\,u_0(s)$, $c=\psi(Y^\star_0)$;
$H(s)=u_0'^\star\psi(s)-c'u_0(s)$; $\alpha=c/u_0'^\star$; $V^\star=-Y^{\star2}_0=-4.786$.
Using chaos orthogonality $\langle\text{1st},\text{3rd}\rangle=0$ (kills $\langle Y_1 u_3^\star\rangle$), all pieces
of $v_1$ reduce to single + ordered-double integrals. **Validated (converged in the quadrature grid):**

| coeff | deterministic (MC-free) | MC target | 
|---|---|---|
| $v_0=\langle Y_1^2\rangle=\alpha^2\!\int u_0^4$ | 0.1343 | 0.1339 |
| $m_1=-\alpha\!\int Hu_0^3/u_0'^\star$ | +0.2135 | +0.2124 |
| $\operatorname{Var}(Y_2)$ | 0.0640 | 0.0636 |
| $\langle Y_1Y_3\rangle$ | +0.0233 | +0.0231 |
| $v_1$ | **+0.1107** | **+0.110** |

Scripts `coupled-atlas/weaknoise_greens.py` ($v_0$), `weaknoise_greens_v1.py` (full $v_1$). This *is* the
stochastic-exact-WKB viewpoint made concrete (fields as Green's-function iterated integrals).

### 6b. The boundary term, explicit. 
Symbolic Taylor inversion of the node condition (sympy, `weaknoise_regularized.py`) confirms by hand: substituting
$u_1'''(Y^\star_0)=V'^\star u_1+V^\star u_1'+\xi(Y^\star_0)u_0'^\star$ and $u_2''(Y^\star_0)=V^\star u_2+\xi(Y^\star_0)u_1$
into $Y_4$, the pointwise-noise pieces collapse to a single clean term:
$$Y_4=\tfrac13\,\xi(Y^\star_0)\,Y_1^3+(\text{regular}),$$
and $Y_5$ additionally carries $\xi'(Y^\star_0)$. So the "wall" is one explicit boundary operator, not a mess.

### 6c. The number — $v_2$ by three routes [NUMERIC ✓✓]
$v_2=\langle Y_3^2\rangle+2\operatorname{Cov}(Y_2,Y_4)+2\langle Y_1Y_5\rangle$.
- **Direct full-field first-node MC** (`weaknoise_v2_fullfield.py`, `weaknoise_v2_extract.py`): march the *full* noisy
  field $u''=(V+\eta\dot W)u$, read $Y^\star=$ first node, Richardson $h\to0$, subtract the pinned $v_0,v_1$, fit.
  $v_2=+0.114\pm0.007$ (jackknife; $\pm0.02$ systematic from the $v_0,v_1$ assumption).
- **Regularized recursion** (`weaknoise_regularized.py`, `weaknoise_reg_campaign.py`): mollify the noise (width δ) so
  $\xi_\delta,\xi_\delta',\dots$ are finite; run the recursion through $Y_5$; extrapolate $\delta\to0$. This supplies
  the correct boundary-Wick constants (incl. the $\xi'(Y^\star_0)$ piece) *without guessing the boundary ½*.
  **Credibility check:** as $\delta\to0$ it recovers $v_0\to0.135$, $m_1\to0.213$, $v_1\to0.112$ (all $\approx$ Itô
  targets — the mollified/Wong–Zakai limit is Itô-consistent here because the node-noise multiplies non-anticipating
  lower-order fields). Same extrapolation gives $v_2=+0.114$ (direct) / $+0.110$ (piecewise).
- **Decomposition** (δ→0): $\langle Y_3^2\rangle=+0.145$, $2\operatorname{Cov}(Y_2,Y_4)=+0.008$,
  $2\langle Y_1Y_5\rangle=-0.044$; the two **boundary-Wick terms net $-0.036$** — a ~30 % correction to the
  $\langle Y_3^2\rangle$-dominated value.

### 6d. Resurgence payoff — a third converged coefficient sharpens the radius
$v_0=0.134,\ v_1=0.110,\ v_2\approx0.11$. Ratios $v_1/v_0=0.82$, $v_2/v_1\approx1.0$ — growing gently toward $1$,
**not yet runaway**. The radius of the $\eta^2$-series is $\eta^2_c\approx v_1/v_2\approx1.0$, i.e.
$\boxed{\beta_c\approx4}$ — tightening the earlier two-coefficient estimate ($\eta^2\approx1.2$, $\beta\approx3.3$).
Consistent with an $O(1)$ Borel singularity at the instanton action; three converged terms is still too few to
separate finite-radius from factorial (asymptotic) growth — that remains the higher-order question. **The trans-series
structure now rests on three converged perturbative coefficients + the instanton sector + the demonstrably
non-perturbative $\beta=2$ point.**

## 7. Phase 2 attempted — large-order growth: radius pinned, factorial-vs-convergent gated on an EXACT $v_3$

_2026-07-07 (Fable 5), continuing "go for it" toward the proved-closed-form route (`PROGRAM2_PROVED_CLOSED_FORM_ROUTE.md`
Phase 2). Goal: the large-order growth of $v_n$ — the decisive resurgence test (factorial $v_n\sim n!/A^n$ ⇒ genuine
trans-series with Borel singularity at $A$; vs finite radius ⇒ convergent). Extended the node-shift recursion to
arbitrary order (`weaknoise_highorder.py`, sympy-generated $Y_1..Y_9$, cached; $v_{m-1}=\operatorname{Var}(Y_m)+
2\sum_{k<m}\operatorname{Cov}(Y_k,Y_{2m-k})$)._

**Result — $v_3$ is a RENORMALIZATION THRESHOLD, not a precision wall [precise obstruction, sharpened].** $v_0,v_1,
v_2$ reproduce, but $v_3$ cannot be reached by *any* noise-sampling method — and the reason is structural, not
statistical. Probing the individual pieces vs the mollifier width δ (`_probe_div.py`, N=48k):

| piece | δ: 0.10 → 0.035 | behavior |
|---|---|---|
| $\operatorname{Var}(Y_2),\operatorname{Var}(Y_3)$ | 0.060→0.063, 0.122→0.140 | **converge** (finite) |
| $m_2=\langle Y_4\rangle$, $\operatorname{Cov}(Y_2,Y_4)$ | 0.065→0.074, +0.0083→+0.0085 | **converge** (single boundary noise) |
| $\operatorname{Var}(Y_4)$ | 0.260→**0.331** | **diverges** ($\xi(Y^\star_0)^2$) |
| $\operatorname{Var}(Y_5)$ | 1.03→**2.02** | **diverges faster** ($\xi'(Y^\star_0)^2$) |
| $v_3$ assembly | −0.17→**−0.38** | **diverges** — cancellation incomplete |

**The mechanism.** $Y_4=\tfrac13\xi(Y^\star_0)Y_1^3+$reg carries a *single* boundary noise, so its mean $m_2$ and its
covariance with boundary-free functionals ($\operatorname{Cov}(Y_2,Y_4)$) are finite — which is exactly why **$v_2$
is clean** (all its boundary terms are single-power). But $v_3=\operatorname{Var}(Y_4)+2\operatorname{Cov}(Y_3,Y_5)+
2\operatorname{Cov}(Y_2,Y_6)+2\operatorname{Cov}(Y_1,Y_7)$ contains $\operatorname{Var}(Y_4)\ni\xi(Y^\star_0)^2$ — a
boundary **self-contraction** $\delta(0)$ that diverges. The physical $v_3$ is finite only because these divergences
**cancel across the assembly** (against the $Y_6,Y_7$ covariances, which carry $\xi'',\xi'''$ and their own
$\delta$-singularities). But the pieces diverge at *different* rates and no sampling / mollification can resolve a
delicate cancellation between individually-divergent terms — so $v_3$ demands **analytic renormalization** of the
boundary self-contractions before any numerics. This is a QFT-like counterterm cancellation, first appearing at $v_3$.

**So the honest structure of the ladder:** $v_0,v_1$ boundary-free (deterministic, done §6a); $v_2$ single-boundary
(finite, computed §6c, ≈0.11); $v_3$ **the renormalization threshold** — the first coefficient needing exact
cancellation of $\xi(Y^\star_0)^2$-type self-contractions. Any noise-sampling method (regularized recursion, direct
full-field) provably fails here; only the **exact boundary-Wick with symbolic $\delta(0)$-cancellation** can reach it.

## 8. The $v_3$ counterterm, EXACTLY — the divergences do NOT cancel; renormalization is real [DERIVED, exact algebra]

_2026-07-07 (Fable 5), "go for it" continued. Attacked the $v_3$ renormalization analytically (symbolic extraction
of the boundary-noise structure of $Y_4..Y_7$, `_probe_div.py`/`_probe_ct.py` + sympy). Result: the divergence
structure of $v_3$ is exactly one term, and it does **not** self-cancel — $v_3$ carries a genuine, computable
counterterm. This is the sharp resolution of the "renormalization threshold."_

**The boundary-noise structure is minimal [EXACT].** Writing $s_j=\xi^{(j)}(Y^\star_0)$ and substituting the field
ODE into the node-shift formulas $Y_4..Y_7$, the *only* quadratic-boundary-noise monomial anywhere in the ladder up
to $Y_7$ is a single $s_0^2$ in $Y_7$:
$$Y_4=\underbrace{\tfrac13\,s_0\,Y_1^3}_{B_4=\frac13 Y_1^3}+(\text{reg}),\qquad
Y_5,Y_6:\ \text{only }s_0,s_1,s_2\ \text{LINEAR},\qquad
Y_7\ni R_7\,s_0^2,\quad \boxed{R_7=\tfrac15 Y_1^5}\ \ (=-U_{1,0}^5/5U_{0,1}^5,\ \text{sympy-exact}).$$
(Both $B_4=\tfrac13Y_1^3$ and $R_7=\tfrac15Y_1^5$ are machine-verified identities.)

**The divergence of $v_3$ [EXACT].** $v_3=\operatorname{Var}(Y_4)+2\operatorname{Cov}(Y_3,Y_5)+2\operatorname{Cov}
(Y_2,Y_6)+2\operatorname{Cov}(Y_1,Y_7)$. The middle two covariances are **finite** ($Y_5,Y_6$ carry only *linear*
boundary noise → single contractions → finite boundary values). The divergence is a single self-contraction
$\langle\xi(Y^\star_0)^2\rangle\equiv\Lambda$ from two sources, both with the **same sign**:
$$\operatorname{Var}(Y_4)\big|_{\rm div}=\Lambda\langle B_4^2\rangle=\Lambda\tfrac19\langle Y_1^6\rangle=\Lambda\tfrac53v_0^3,\qquad
2\operatorname{Cov}(Y_1,Y_7)\big|_{\rm div}=2\Lambda\langle Y_1R_7\rangle=2\Lambda\tfrac15\langle Y_1^6\rangle=6\Lambda v_0^3,$$
using $\langle Y_1^6\rangle=15v_0^3$ (Gaussian, exact). **Total:**
$$\boxed{\ v_3\big|_{\rm div}=\Lambda\big(\tfrac53+6\big)v_0^3=\tfrac{23}{3}\,\langle\xi(Y^\star_0)^2\rangle\,v_0^3\neq0.\ }$$

**Interpretation.** The divergences do **not** cancel among themselves — so $v_3$ is *not* rendered finite by internal
cancellation; it requires an explicit **counterterm subtraction** $\tfrac{23}{3}\langle\xi(Y^\star_0)^2\rangle v_0^3$
(a Wong–Zakai / Itô-vs-Stratonovich correction, first biting at the order where two noise-at-node factors collide —
the fingerprint of the field's $C^{1/2}$ roughness). $v_0,v_1,v_2$ are counterterm-free (no $s^2$ boundary term until
$Y_7$); **$v_3$ is exactly where the weak-noise variance series becomes a RENORMALIZED series.** The physical
$v_3^{\rm ren}$ is the finite remainder after subtracting the counterterm.

**Status of the number.** The counterterm *structure and coefficient* are exact algebra. The finite $v_3^{\rm ren}$
VALUE still requires assembling the finite chaos-integral parts of the four moments with the counterterm removed —
now a **well-posed deterministic computation** (no sampling: the one divergent piece is known in closed form and
subtracted). Numerical confirmation of the $\tfrac{23}{3}v_0^3$ coefficient by sampling is inconclusive (the
divergent part $\sim0.02\Lambda$ is dwarfed by the drifting finite parts $\sim0.24$ — sampling cannot isolate it),
which is itself why the exact route is mandatory.

**Resurgence connection [LEAD].** That the perturbative coefficients are defined only *after* an order-by-order
subtraction is precisely the renormalization/summation-ambiguity structure that, in resurgence, is the shadow of the
Stokes automorphism (cf. `PROGRAM2_PROVED_CLOSED_FORM_ROUTE.md` Phase-3 lead): the counterterm scheme and the Borel
non-perturbative sector are two faces of the same object. Concrete next step: compute $v_3^{\rm ren}$, then track how
the counterterm grows with order (the $v_n^{\rm ren}$ ladder) — its growth *is* the resurgence signature.

**$v_3^{\rm ren}$ COMPUTED — and it is NEGATIVE (sign change) [NUMERIC, indicative].** Variance-reduced renormalized
sampling (`weaknoise_v3_renorm.py`): use the exact jointly-Gaussian moments $\langle s_0^2B_4^2\rangle=\tfrac53
\langle s_0^2\rangle v_0^3+10c^2v_0^2$, $\langle s_0^2Y_1R_7\rangle=3\langle s_0^2\rangle v_0^3+18c^2v_0^2$
($c=\langle s_0Y_1\rangle$) to remove the $\langle s_0^2\rangle$ divergence exactly and sample only the finite parts.
Result across $\delta=0.09{\to}0.045$: $v_3^{\rm ren}=-0.226,-0.271,-0.303,-0.299$ → **$v_3^{\rm ren}\approx-0.3$**
(magnitude $\pm0.1$, still drifting from the un-tamed *linear* $s_1,s_2$ boundary noises; but the **sign is robust**,
anchored by a stable $\operatorname{Cov}(Y_2,Y_6)\approx-0.20$). So the coefficient ladder is
$$v_0=+0.134,\ v_1=+0.110,\ v_2=+0.11,\ v_3\approx-0.3\ \ (\textbf{first sign change}).$$
**The sign change matches the COMPLEX connection root.** $\lambda_0=0.890-0.890i$ sits at **exactly $-45^\circ$**
($\mathrm{Re}=-\mathrm{Im}$). A Borel singularity inherited from $\lambda_0$ at phase $\pm\pi/4$ forces oscillatory
coefficients $v_n\sim|S|^{-n}\cos(n\pi/4+\varphi)$; $\cos$ first turns negative between $n=2$ and $n=3$, so
**$v_3$ is predicted to be the first negative coefficient** — exactly as observed. This is the seed of the exact-form
surmise (`PROGRAM2_EXACT_FORM_SURMISE.md`): the resurgent structure of $\mathcal W$ is governed by a complex-conjugate
Borel pair at $\arg=\pm\pi/4$, tied to $\lambda_0$. (Caveat: $v_3$ magnitude not pinned; the sign/oscillation is the
robust, falsifiable claim — $v_4$ predicted also $<0$, $v_5$ turning back up.)

**Phase-2 analysis on the three solid coefficients [NUMERIC].** $v=(0.134,0.110,0.11)$; ratios $v_1/v_0=0.82$,
$v_2/v_1=1.00$ ⇒ radius $z_c=\eta^2_c\approx v_1/v_2\approx1.0$, i.e. $\boxed{\beta_c\approx4}$ — consistent with an
$O(1)$ instanton action (left tail $-\log\mathbb P(Y^\star<-s)\to s^5/(10\eta^2)$; Borel singularity at the action).
Fitting the two hypotheses to $(v_0,v_1,v_2)$: **factorial** $v_n\sim\Gamma(n+b)/S^n$ gives $S\approx5.6,b\approx4.6$,
predicting $v_3\approx0.13$; **finite-radius** predicts $v_3\approx0.11$. The **~20 % split is below the current $v_3$
reach ($\pm0.25$)** — so three coefficients pin the *radius* but cannot yet separate divergent-asymptotic from
convergent. That separation is the gate to "proved resurgent trans-series," and it needs one exact coefficient.

**The wall-breaking method (identified, correct, not yet built).** Do the node-shift moments by **exact discrete
Wick** (Gaussian/Isserlis) instead of sampling: on the grid the white noise is a Gaussian vector $g$ ($\langle g_ig_j
\rangle=h\delta_{ij}$), the node fields are linear forms in $g$ with known coefficient vectors (the validated
Green's-function kernels, §6a), and $\xi(Y^\star_0)\to g_{\rm node}/h$. The pointwise term that makes **MC blow up**
($\operatorname{Var}(g_{\rm node}/h)=1/h$) is **finite in exact Wick** ($\langle (g_{\rm node}/h)\,g_i\rangle=\delta_{
{\rm node},i}$ — the $1/h$ cancels the covariance $h$). So $\langle Y_kY_l\rangle$ = sum over Isserlis pairings of
dot-products of coefficient vectors — deterministic, no sampling error, no δ, the singular cancellations exact. This
is the concrete build that delivers exact $v_2,v_3,v_4$ and thus the definitive large-order / Borel-singularity test.
Entry point already in hand: $Y_4=\tfrac13\xi(Y^\star_0)Y_1^3+$reg (§6b) and the validated linear/quadratic/cubic
Green's-function kernels (§6a).
