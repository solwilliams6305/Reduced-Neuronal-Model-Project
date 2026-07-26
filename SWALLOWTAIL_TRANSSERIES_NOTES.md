# The swallowtail (A₄, q=3) trans-series and its operator

_July 2026. Built to the 𝒲 (cusp) paper's standard, to the depth the methods reach.
Engine `coupled-atlas/weaknoise_q.py` (q-parameterized, validated against the cusp).
Tags: [PROVED] / [NUMERIC] / [OPEN]._

## The operator that describes it

The swallowtail escape law is the first-node / first-explosion law of the **stochastic
cubic-turning-point operator**
$$u'' = \big(\operatorname{sign}(Y)\,|Y|^3 - \eta\,\dot W\big)u \;=\; (Y^3-\eta\dot W)u,\qquad \beta=4/\eta^2,$$
the q=3 member of the swept multicritical family $u''=(\operatorname{sign}(Y)|Y|^q-\eta\dot W)u$.
Structural facts:

- **Deterministic backbone = order-1/5 Bessel.** $u''=Y^3u$ solves in
  $\sqrt{|Y|}\,\mathcal C_{1/5}\!\big(\tfrac25|Y|^{5/2}\big)$; the recessive solution is the
  **symmetric** combination $J_{1/5}+J_{-1/5}$, and the peel-off levels are its zeros:
  $-2.046707,\,-2.856553,\,-3.419883,\,-3.870609,\,-4.253976$ (matched to $10^{-14}$). **[PROVED]**
- **Z₅ (5-fold) Stokes structure.** The rotation $Y\to e^{2\pi i/5}Y$ maps $\lambda\to e^{4\pi i/5}\lambda$
  (from $e^{-2i\phi}=e^{3i\phi}\Rightarrow\phi=2\pi/5$), so the cubic operator has five Stokes sectors —
  versus the cusp's four (Weber, Z₄). **[PROVED]**
- **Provably non-integrable, a fortiori.** Both cusp obstructions apply and are stronger: the operator
  is unbounded below ($Y^3<0$ for $Y<0$ ⇒ no Fredholm/soft-edge determinant), and the Bessel-1/5
  skeleton is an isomonodromy fixed point further from any Painlevé flow (⇒ no $\sigma$-form). So the
  only closed form is a resurgent trans-series — as for 𝒲. **[PROVED]**
- **The resonance/connection root is a cubic-oscillator connection datum — NOT elementary.** For the
  cusp, the deterministic connection was the closed-form Weber Gamma-equation (root $\lambda_0$, phase
  exactly $-\pi/4$). The swallowtail analog is the connection problem for $u''=(Y^3-\lambda)u$, i.e. the
  **cubic anharmonic oscillator** (Bender–Wu / Voros) — its connection coefficients are genuinely
  transcendental, with no Gamma-function closed form. **This is the key structural finding: the
  swallowtail is one level *less* integrable than the cusp.** The resonance root is complex (from the
  unbounded-below side) but carries no clean closed phase. **[NUMERIC/OPEN]**

## The trans-series (weak-noise variance coefficients)

$\operatorname{Var}(Y^\star)=\eta^2(v_0+v_1\eta^2+v_2\eta^4+\cdots)$, node $Y^\star_0=-2.046707$.
Computed by the q-parameterized Wiener-chaos engine (`weaknoise_q.py`), **validated by reproducing the
cusp's $v_0,v_1$** at q=2 (v0→0.1346, v1→0.1049 as δ→0; targets 0.134, 0.110):

| coeff | value | status |
|---|---|---|
| $v_0$ | **0.04953187** | [PROVED] exact one-loop $v_0=\int_{Y^\star_0}^\infty\varphi^4/\varphi'(Y^\star_0)^4$; engine δ→0 confirms 0.0495 |
| $v_1$ | $\approx 0.067$ | [NUMERIC] engine, δ→0 extrapolated |
| $v_2$ | $\approx 0.118$ | [NUMERIC] engine, δ→0 extrapolated |
| $v_3$ | $\approx 0.057$ | [SUPERSEDED — rough δ-sensitive value; production engine gives 0.177, see below] |

**Homogeneous production ladder (2026-07-24).** The table above splices two engines. For any
Domb–Sykes/Borel work the ladder must be on *one* normalization, so $v_0,v_1,v_2$ were recomputed on
the **production** engine (`swtl_lowrungs.py`, same driver as the campaign's $v_3\ldots v_6$), grids
$n=24,30,36$:

| coeff | production $n\to\infty$ | band over estimators | grid $n=24,30,36$ |
|---|---|---|---|
| $v_0$ | **+0.04969** | — (grid-independent) | 0.049690, 0.049690, 0.049690 |
| $v_1$ | +0.0632 | [0.0632, 0.0652] | 0.065973, 0.065537, 0.065215 |
| $v_2$ | +0.1108 | [0.1095, 0.1177] | 0.121837, 0.119317, 0.117723 |
| $v_3$ | +0.1770 | [0.1732, 0.1828] | 0.187562, 0.184549, 0.182793 |
| $v_4$ | +0.053 | [0.043, 0.086] ⚠ | 0.048035, 0.045638, 0.046804 |

- **Normalization verified:** production $v_0=0.049690$ vs the exact $0.04953187$ → ratio **1.0032**
  (grid-level); $v_1,v_2$ also track the independent weak-noise values (0.067, 0.118). So splicing
  $v_0..v_2$ with $v_3..v_6$ is legitimate.
- ⚠ **$v_4$ is poorly determined.** Its grid sequence is **non-monotone** (0.0480, 0.0456, 0.0468), so
  the campaign's quadratic-in-$1/n$ Richardson (which reports $+0.0858$, nearly double every raw
  value) is fitting curvature to noise and should **not** be used. Estimator spread 0.043–0.086.
- **The $v_3\to v_4$ collapse** (0.177 → 0.053) is the ladder turning over — the swallowtail analogue
  of the cusp's $v_2\to v_3$ sign change (0.104 → −0.030), one rung later. $v_5$ is the decisive
  coefficient. [NUMERIC]

- **Sign pattern $+,+,+,+$** (contrast the cusp's $+,+,+,-$): no early negative run.
- **Factorially divergent / asymptotic** (growing $v_0<v_1<v_2$: 0.0495, 0.067, 0.118), confirming the
  resurgent-trans-series hypothesis. Borel radius $\eta^2_c\approx0.55$ (tighter than the cusp's ≈1.2 —
  the steeper catastrophe is more non-perturbative). **[NUMERIC]**

## The β=2 ground truth — independently verified (2026-07-24)

The target $\operatorname{Var}=0.328$ had been *asserted* in these notes with no traceable
computation in the repo. It is now reproduced from the MC-free FP-PDE (`swtl_groundtruth.py`, the
q-generalization of `fp_cusp.solve_fp`: $V=\operatorname{sign}(Y)|Y|^q$, $p_0=Y_0^{q/2}$):

| q | $\operatorname{Var}(\beta{=}2)$ | $f(2)=\operatorname{Var}/\eta^2$ | skew | ex. kurt |
|---|---|---|---|---|
| 2 (gate) | 0.4740 / 0.4715 | **0.2370 / 0.2358** | +0.601 | −0.243 |
| **3** | **0.3324 / 0.3309** | **0.1662 / 0.1654** | **+0.945** | **+0.150** |

(two resolutions, $dp=0.02/0.01$). The q=2 row is a **regression gate** and reproduces the cusp
ground truth (0.237, +0.607, −0.237). So **Var(β=2) ≈ 0.331** — the asserted 0.328 is confirmed to
1%, and the resummation target is $f(2)\approx0.165$. **[NUMERIC ✓]**

- **New: the swallowtail β=2 fingerprint differs qualitatively from the cusp's.** Skew
  **+0.945** (vs cusp +0.601) and excess kurtosis **+0.150** — *positive*, where the cusp's is
  *negative* (−0.243). The steeper catastrophe is both more skewed and heavier-tailed; the cusp's
  sub-Gaussian bulk signature does not carry over to q=3. **[NUMERIC]**

## Borel plane

- **Real Freidlin–Wentzell instanton (closed-form anchor).** Left-tail rate
  $-\log P(Y^\star<-s)\to I(s)/\eta^2$ with
  $$I(s)=\frac{|s|^{2q+1}}{2(2q{+}1)}\;\xrightarrow{\;q=3\;}\;\frac{|s|^7}{14};$$
  right tail $\to\tfrac43|s|^{4.5}$. The real Borel singularity sits on $\mathbb R_+$ at the reduced
  escape action, exponent **7** (vs the cusp's 5). **[PROVED, exponents; DERIVED, constants]**
  - **CORRECTION (2026-07-24).** This note previously quoted $|s|^{2q+1}/[4(2q{+}1)]=|s|^7/28$ — a
    factor **2 too small**. At $q=2$ that formula gives $s^5/20$, contradicting the *numerically
    verified* cusp value $s^5/10$ (`instanton_action.py`, robust in $p_0,M$) and the 𝒲 paper's
    eq. (instanton). Analytic reason: on the oscillatory side $V=-t^q<0$ the deterministic Riccati
    $\dot p=V-p^2$ explodes; holding it off to depth $s$ forces $p\approx0$, hence $\pi\approx t^q$ and
    $I(s)=\tfrac12\int_0^s t^{2q}dt=s^{2q+1}/[2(2q{+}1)]$ — which reproduces $1/10$ at $q=2$.
    Re-verified by the same BVP machinery at both $q$ (`instanton_action_q.py`): the measured
    $I/s^{2q+1}$ → $0.0930$ (q=2, pred $0.1$) and $0.0701$ (q=3, pred $0.07143$), i.e. ratio-to-notes
    $\to1.86$ and $1.96$. **The anchor is $s^7/14$.**
- **Complex pair** (the cusp had a conjugate pair at $\theta\approx50^\circ$, $|\zeta|\approx1.9$): would
  require the fuller coefficient ladder to locate; not resolved with four coefficients. **[OPEN]**
  → **Now resolved in modulus, NOT in phase** with seven coefficients — see the campaign results below.

## CAMPAIGN RESULT: the 7-coefficient ladder (2026-07-24)

The campaign completed ($v_3..v_6$; $v_0..v_2$ recomputed on the same engine). Analysis in
`coupled-atlas/swtl_borel.py`, whose machinery is **validated** by reproducing the 𝒲 paper's
published cusp resummation to 3 decimals (0.2244/0.2337/**0.2347**/0.2337 vs paper
0.224/0.234/**0.235**/0.234 at $x=1,1.5,2,2.25$).

| $k$ | $v_k$ (best) | honest band | grid |
|---|---|---|---|
| 0 | **+0.04969** | — | 24,30,36 (flat) |
| 1 | +0.0632 | [0.0632, 0.0652] | 24,30,36 |
| 2 | +0.1108 | [0.1095, 0.1177] | 24,30,36 |
| 3 | +0.1770 | [0.1732, 0.1828] | 24,30,36 |
| 4 | +0.031 | [0.031, 0.044] | 24,30,36,42,48 |
| 5 | −1.212 | [−1.538, −1.000] | 16,20,24 |
| 6 | −4.71 | [−11.76, −4.71] ⚠ | 10,12,14 |

**(1) Sign pattern $+,+,+,+,+,-,-$** (cusp: $+,+,+,-,-,-,-$). **Factorial divergence confirmed**:
the envelope collapses to a node at $|v_4|\approx0.03$, then explodes $|v_5|\approx1.2$,
$|v_6|\approx5\text{–}12$. This is *exactly* the cusp's stated fingerprint — "small $|v_3|$ (cosine
node), large $|v_4|$ (anti-node), growing envelope" — **shifted one rung later**, which is itself the
signature of a *smaller* Borel phase $\theta$ (the cosine takes longer to reach its first zero).

**(2) Domb–Sykes** on $b_k=v_k/k!$ gives a **negative** intercept ($1/\zeta=-2.66$): no real positive
singularity dominates the large-order behaviour, i.e. the growth is pair-driven. The fitted "$\zeta$"
is not literal.

**(3) Borel plane.** Complex pair present and its **modulus is stable**: $|\zeta|\approx1.1\text{–}1.5$
(vs the cusp's 1.9). Real positive poles appear at 1.24, 1.31, 1.54, 3.04, 3.46; with the corrected
anchor $I(s)=s^7/14$ the median $|\zeta_r|=1.54$ corresponds to escape depth $s\approx1.55$.
**But $\theta$ is NOT determined** — it drifts *monotonically* with the number of coefficients:

| $K$ (uses $v_0..v_K$) | 3 | 4 | 5 | 6 |
|---|---|---|---|---|
| $\theta$ | 34.4° | 36.4° | 49.0° | 60.4° |

A 26° monotone drift with no sign of a plateau. This is precisely the failure mode the 𝒲 paper had
to **retract** (its 5/6-coefficient "$\theta\approx54$–$63^\circ$" was a small-sample artifact). We
therefore quote **no $\theta$** for the swallowtail. **[OPEN]**

**(4) Median-Borel resummation — DOES NOT reconstruct the ground truth.** Against the verified
$f(2)=0.1655$:

| $x$ | $\beta$ | naive (7-term) | median resummation | truth |
|---|---|---|---|---|
| 1.0 | 4.0 | −5.5 | 0.148 | — |
| 1.5 | 2.67 | −61.7 | 0.147 | — |
| **2.0** | **2.0** | **−337.8** | **0.126** | **0.1655** |

24% short, and — decisively — the spread across Padé orders/contour angles is **255% of the median**
(IQR 0.045–0.365). The same code reproduces the cusp to 1%. Worse, the central-ladder value (0.126)
lies *outside* the range spanned by all eight extrapolation-band corners (0.011–0.047): the
approximant is not behaving monotonically in the data, which is the clearest single statement that
this is **not** a reconstruction. **[NEGATIVE RESULT, honest]**

### Why q=3 is harder than q=2, and what would actually fix it

- $\beta{=}2$ is *relatively deeper* for the swallowtail: $x/|\zeta|\approx2/1.2\approx1.7$ versus the
  cusp's $2/1.9\approx1.05$. Same seven coefficients, further outside the disc.
- **The binding constraint is grid convergence, not coefficient count.** $v_6$ is known only to a
  factor $\approx2.5$ ($-11.8$ to $-4.7$) from grids $n=10,12,14$, which still move 74% between the
  first two points; $v_5$ rests on $n=16,20,24$. The cusp quoted $v_6=-1.90\pm0.15$ (8%). Getting the
  swallowtail Borel plane to the 𝒲 standard needs **higher-$n$ grids at $k=5,6$** (and ideally $v_7$,
  for which `_yexprs_15.txt` already exists), not more rungs at the current resolution.
- ⚠ **The campaign's own Richardson values are unusable at the top of the ladder** and should not be
  read off `swtl_results.json`: its quadratic-in-$1/n$ fit reports $v_5=-1.962$ and $v_6=-103.07$,
  both *below every grid point* while both sequences are rising. `swtl_borel.py` screens estimators
  for wrong-side extrapolation, bad fit residual, and sign flips, and reports bands instead.

**Net.** Structure — divergence, the node/anti-node fingerprint, a complex pair of modulus
$\approx1.2$, a real instanton at $s^7/14$ — is established. Quantitative Borel data ($\theta$) and
the end-to-end resummation are **not**, and the earlier note that these were "compute-bound
extensions, not obstructions of principle" survives, but the compute needed is *finer grids*, not
*more coefficients*.

### The instanton-ANCHORED Darboux fit — and why q=3 fails where the cusp succeeds

Raw Borel–Padé poles mix the complex pair with the real FW singularity, and that contamination is
what biased the cusp's early (retracted) $\theta$. The proper tool pins the real term at a
**profiled** $z_r$ (scanned, not fitted) and extracts the pair at each assumed instanton position:
$$v_n \sim 2C|\zeta|^{-(n+1)}\Gamma(n{+}1{+}\alpha)\cos((n{+}1)\theta-\varphi)\;+\;C_r z_r^{-(n+1)}\Gamma(n{+}1{+}\alpha).$$
Implemented as `anchored_darboux()` in `swtl_borel.py`. **Validated on the cusp**: it returns
$\theta=50.0^\circ$, $|\zeta|=1.85$ — the paper's $50^\circ\pm2$, $|\zeta|\approx1.9$ — with a
$z_r$-spread of only **1.5°** at fixed $\alpha$.

Applied to the swallowtail (with the coarse $v_5,v_6$) the same fit gives a $z_r$-spread of
**18.5°** at $\alpha=0$ — *while the residuals stay small* (0.03–0.08). That combination is the
decisive diagnostic: the model fits fine, but **many $(z_r,\theta)$ pairs fit equally well.** The
q=3 $\theta$ is genuinely **unidentifiable**, not merely badly fitted — the same "6-coeff
unidentifiable + real-pole bias" pathology recorded for the cusp in
[[program2-milestone2-borel-phase]]. **[NUMERIC, honest negative]**

| | cusp q=2 | swallowtail q=3 |
|---|---|---|
| $\theta$ ($\alpha{=}0$) | 50.0° | 48.0° (meaningless — see spread) |
| $z_r$-spread at fixed $\alpha$ | **1.5°** | **18.5°** |
| $|\zeta|$ | 1.85 | 1.22 |
| max fit residual | 0.11–0.12 | 0.03–0.08 |

## Fine-grid campaign for $v_5,v_6$ — COMPLETE, and it REFUTES the diagnosis above (2026-07-24)

All six jobs landed (54 min, 6-way parallel). **The finer grids did not help the Borel analysis at
all** — which retracts the "the limiter is grid convergence, not rung count" claim made above.

| $k$ | grid | values | verdict |
|---|---|---|---|
| 5 | 16,20,24,28,32,36 | −1.7944, −1.6030, −1.5379, −1.5016, −1.4991, **−1.4961** | **CONVERGED** |
| 6 | 10,12,14,16,18,20 | −22.453, −12.939, −11.763, −10.136, −9.564, **−9.182** | still moving |

- **$v_5=-1.496\pm0.005$ — converged.** Gaps +0.191, +0.065, +0.036, **+0.0025, +0.0030**: the last
  two are at the numerical noise floor and *non-monotone*, so the sequence has plateaued and any
  further $1/n$ extrapolation is fitting noise. Note this is a **24% revision** from the $-1.21$ the
  coarse grids implied — the $n=16$–$24$ points were pre-asymptotic, so every $1/n$ extrapolation off
  them over-reached. (`swtl_borel.py` now carries a plateau detector for exactly this.)
- **$v_6\in[-9.18,-5.34]$ — improved (factor 2.5 → 1.7) but NOT converged.** Its tail gaps (+0.57,
  +0.38) are still tapering at the normal $\sim0.8\times$ per step.
- **Yet the resummation is no better:** $f(2)=0.035$ with a converged $v_5$, versus 0.126 with the
  coarse one. **The failure is not a grid-convergence artifact.**

### Protocol correction: "median" means lateral, not multi-order

The 𝒲 paper's recipe is the **diagonal** Padé (for 7 coefficients, $[3/3]$), with *median* denoting
$\tfrac12(f_++f_-)$ — the two **lateral** Borel sums — not a median over Padé orders. Taking the
latter is a different and much worse statistic: on the *cusp* it returns 0.136 at $x=2$ against the
published 0.235 and truth 0.237. `swtl_borel.py` now leads with the diagonal and reports the
multi-order median only as a secondary robustness check.

Apples-to-apples, diagonal $[3/3]$, median of lateral Borel sums:

| $x$ | cusp $[3/3]$ | cusp truth | swallowtail $[3/3]$ | swallowtail truth |
|---|---|---|---|---|
| 1.0 | 0.2244 | 0.232 | 0.0666 | — |
| 1.5 | 0.2337 | 0.239 | 0.0478 | — |
| **2.0** | **0.2347** | **0.237** | **0.0354** | **0.1655** |

The swallowtail is not merely low (79%): it **decreases** in $x$ where the truth must *rise* from
$v_0=0.0497$ to 0.1655. The approximant has the wrong qualitative shape, and it is completely
insensitive to the contour angle ($\varphi=20^\circ$–$45^\circ$ give identical values), so this is
not a contour artifact.

### What the failure is NOT (two hypotheses killed)

1. **Not relative depth.** $\beta{=}2$ sits at $x/|\zeta|\approx1.8$ for the swallowtail vs 1.05 for
   the cusp, so "further outside the disc" was the natural explanation. **Refuted**: the cusp's own
   $[3/3]$ holds 0.2347 → 0.2196 out to $x/|\zeta|=2.11$, a 6% drift — it does not collapse. The
   swallowtail is $\sim4\times$ low at *every* $x$, including $x=1.16$ where its relative depth
   matches the cusp's $\beta{=}2$ exactly.
2. **Not a bad ladder.** Tested by comparing the series to FP at small $x$ — but see the caveat
   below; the test came back **inconclusive-for-the-ladder and damning-for-FP** instead.

### ✅ RESOLVED (2026-07-25): the "FP bias" is first-order upwind diffusion, not a bias

**The caveat below is retracted.** The FP-PDE is *not* biased — it is first-order accurate in $dp$,
and the discrepancy is **numerical diffusion from the upwind advection**, $D_{\rm num}\sim|{\rm
drift}|\,dp/2$. That is *independent of $\eta$*, so it swamps the physical $D=\eta^2/2$ at small
$\eta$ — which is exactly why the error looked like a small-$\eta$ bias. Decisive test at $q=2$,
$\eta=0.3$ (trusted cusp series = 0.14477):

| $dp$ | 0.04 | 0.02 | 0.01 | 0.005 | 0.0025 |
|---|---|---|---|---|---|
| $f$ | 0.23041 | 0.18795 | 0.16674 | 0.15614 | 0.15085 |
| err | +59.2% | +29.8% | +15.2% | +7.9% | +4.2% |

The error **halves as $dp$ halves** — clean $O(dp)$ convergence. Richardson in $dp$ on the last two
points gives $f\to0.1456$ vs the series' 0.14477: **agreement to 0.5%**. So the ladder and the FP
agree once discretization is removed; both are sound.

**Why this matters — it opens a route that needs no new rungs.** With $dp$-extrapolation (or a
second-order/limiter scheme) the ground truth becomes accurate enough to *subtract the known
perturbative part* and read the non-perturbative remainder directly, giving the Borel data
($A$, and $\theta$ via its oscillation in $x$) from the observable rather than from the ladder's
large-order behaviour. Window estimate for $q=3$: at $x\approx0.3$–0.4 the optimal-truncation
ambiguity is 0.4–0.9% of $f$ while the non-perturbative term $\sim e^{-A/x}$ is 1.8–5% — a
signal-to-ambiguity of $\sim4$–5×, workable if FP is pushed to $\sim0.1\%$ (one or two more $dp$
halvings plus Richardson). **This is far cheaper than $v_8$ and is the recommended next attack.**
**[NUMERIC, route newly opened]**

### ⚠ SUPERSEDED caveat (kept for the record): "the FP-PDE is biased high at small η"

The $\eta$-sweep table earlier in this file is **only trustworthy near $\beta=2$**. Control at $q=2$,
where the cusp series is independently trusted:

| $\eta$ | $x$ | FP $f$ | trusted cusp series | |
|---|---|---|---|---|
| 0.3 | 0.09 | 0.1667 | **0.1448** (converged) | FP **+15%** |
| 0.4 | 0.16 | 0.1666 | 0.1539 | FP +8% |
| $\sqrt2$ | 2.0 | 0.2370 | — | matches paper ✓ |

The same upward bias appears at $q=3$ (FP 0.0670 vs series 0.0564 at $x=0.09$, +19%). It is **not**
the initial condition — rescaling the IC width with $\eta$ changes the answer by $<10^{-5}$. Since
the bias is present at $q=2$ where the coefficients are known-good, **it is an FP artifact, and it
does not invalidate the $q=3$ ladder.** The $\beta=2$ value (Var $=0.331$) still passes its $q=2$
gate and stands. **[NUMERIC, open artifact]**

## $v_7$ ATTEMPTED — computed, but WALLED (2026-07-25)

Acting on the "needs more rungs" diagnosis, $v_7$ was computed at $q=3$. It required building
`_ybase_q3_14/15.pkl` from `_yexprs_15.txt` (**82 min**), then the moment assembly:

| $n$ | $v_7$ | assembly time |
|---|---|---|
| 10 | −20.617957 | 1828 s |
| 12 | −45.417927 | 20239 s |

**Neither value is usable**, for a structural reason that also retro-explains the $v_6$ mess:

> **RESOLUTION RULE: a grid of $n$ points cannot resolve $Y_{\rm idx}$ when ${\rm idx}>n$.**
> Such points are not merely noisy, they are systematically wrong.

| $v_k$ | max chaos idx | grids | behaviour |
|---|---|---|---|
| $v_4$ | 9 | 24–48 (all $\gg9$) | smooth |
| $v_5$ | 11 | 16–36 (all $>11$) | smooth, plateaus |
| $v_6$ | 13 | **10, 12** $<13$ | **−22.45, −12.94 — 74% swing** |
| $v_6$ | 13 | 14–20 $\ge13$ | −11.76, −10.14, −9.56, −9.18 smooth |
| $v_7$ | 15 | **10, 12** $<15$ | **−20.62, −45.42 — moves the WRONG way** |

$v_7$'s two points straddle nothing: unlike $v_6$ they move *away* from zero, and a two-point $1/n$
fit gives $-169$, which is meaningless. `swtl_borel.py` now **auto-drops** any grid with
$n<\max{\rm idx}$ (this removed $v_6$'s $n=10,12$ and excluded $v_7$ entirely); $v_6$'s central value
is unchanged at $-5.34$ since it already rested on the $n\ge16$ tail.

**The wall.** Assembly cost measured across the two points scales as $n^{13.2}$ (11.1× for
$10\to12$). Projections to the first *usable* resolution:

| $n$ | projected assembly | usable? |
|---|---|---|
| 14 | 43 h | no ($<15$) |
| **16** | **10.4 days** | **first usable** |
| 18 | 49 days | — |
| 20 | 197 days | — |

So a proper $v_7$ grid $(16,18,20)$ is $\sim250$ days of compute. **$v_7$ is out of reach with this
engine** — matching the cusp's own experience (`[[program2-v6-theta48]]`: "v7 walled by assembly
scale"). $v_8$ is further out still: it needs idx 17, and `_yexprs_17.txt` does not exist (only
orders 7,9,11,13,15 are built). ⚠ Note `_yfile(17)` silently returns `_yexprs_13.txt`, so a naive
$v_8$ run would `IndexError` rather than report the missing table — the same trap that broke $v_6$.
**[NUMERIC, hard wall]**

### $v_7$ LANDED (2026-07-26): $v_7(n{=}16)=-25.7215$

The $n=16$ run completed in **61027 s ≈ 17 h** — the first *usable* resolution ($n\ge$ idx 15; the
$n=10,12$ points stay below threshold and are auto-dropped). One usable grid point, so no
extrapolation: $v_7=-25.72$ with no band. Sign pattern is now $+,+,+,+,+,-,-,-$.

**The 8th coefficient transforms the resummation.** Against the verified truth $f(2)=0.1655$:

| ladder | approximant | $f(1)$ | $f(1.5)$ | $f(2)$ | error |
|---|---|---|---|---|---|
| $v_0..v_6$ | $[3/3]$ | 0.0666 | 0.0478 | 0.0354 | **79% low** |
| $v_0..v_7$ | $[3/4]$ | 0.1563 | 0.1913 | **0.2177** | **32% high** |

Two qualitative gains, not just a smaller number: the truth is now **bracketed** (79% low → 32%
high), and $f(x)$ **increases** with $x$ as it must (the 7-coefficient result *decreased*, which was
flatly wrong — $f$ has to rise from $v_0=0.0497$ to 0.1655). This is what a slowly-converging Padé
sequence over a *correct* representation looks like, and it retires the worry that the
representation itself was at fault. It also confirms the earlier blind prediction (a scan over
plausible $v_7$ gave 0.22–0.25; the actual $v_7=-25.7$ gives 0.2177). **[NUMERIC]**

$v_7$ does **not** rescue $\theta$ from the joint fit — that stays at 38.1° (stable across all four
skew grids) — but it does stabilise Borel–Padé's $\theta$, which forced the calibration analysis
above. **A second usable point ($n=18$, ~35 h) would give $v_7$ a band; it is the single highest-value
remaining computation.**

**What $v_7$ was predicted to buy, before it existed.** Scanning a plausible $v_7$ through the near-diagonal
$[4/3]$ Padé: $f(2)$ lands at 0.22–0.25 for $v_7\in[-20,-8]$, versus 0.035 from the 7-coefficient
$[3/3]$ and truth 0.1655. So the eighth coefficient **flips the error from 79% low to ~40% high** —
the truth becomes *bracketed*, and $\theta$ settles to 51–56°, $|\zeta|\approx1.1$. That
order-to-order oscillation is the signature of a slowly-converging (not broken) Padé sequence, which
supports the "more rungs" reading — but it cannot be cashed in. **[NUMERIC, indicative only]**

## Does the swallowtail belong to the same ladder family? — STRUCTURALLY YES (2026-07-25)

Worth stating precisely, because the resummation failure is easy to misread as evidence against
membership. It is not: it is a numerical method running out of coefficients.

**Confirmed family structure at $q=3$** (each the $q$-general pattern evaluated at $q=3$):

| family law | general form | $q=2$ (cusp) | $q=3$ (swallowtail) |
|---|---|---|---|
| backbone skeleton | Bessel-$1/(q{+}2)$ | Weber / $1/4$ | **Bessel-$1/5$** ✓ |
| Stokes structure | $\mathbb Z_{q+2}$ | $\mathbb Z_4$ | **$\mathbb Z_5$** ✓ |
| exact one-loop | $v_0^{(q)}=\int\varphi^4/\varphi'(\text{node})^4$ | 0.1343 | **0.04953** ✓ |
| tail exponent | $2q+1$ | 5 | **7** ✓ |
| divergence | Gevrey-1 | ✓ | **✓** |
| Borel architecture | pair + real instanton, $\cos(k\theta-\varphi)$ envelope | node at $k{=}3$ | **node at $k{=}4$** ✓ |

**The instanton constant is now a verified FIVE-member family law**, not a two-point coincidence
(`instanton_action_q.py`, BVP at $s=10$; convergence to the asymptote is from below, faster at larger
$q$):

| $q$ | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| predicted $\tfrac1{2(2q+1)}$ | 0.16667 | 0.10000 | 0.07143 | 0.05556 | 0.04545 |
| measured | 0.13987 | 0.09297 | 0.07010 | 0.05538 | 0.04543 |
| ratio | 0.839 | 0.930 | 0.981 | 0.997 | 0.9995 |

$$\boxed{I(s)=\frac{s^{2q+1}}{2(2q+1)}}\quad\text{[NUMERIC, }q=1..5\text{]}$$

**What is NOT confirmed** is the *quantitative resurgent data*: $\theta(q)$ (unidentifiable at 7
coefficients) and the end-to-end representation test (resummation reproducing $\operatorname{Var}$).
So the open question is not "is it in the family" but "**does the family have a $\theta(q)$ law**" —
which needs $\theta$ at two $q$, and we have it at one. Note the weak positive evidence that the
representation itself is sound: including $v_7$ moves $f(2)$ from 0.035 to 0.22–0.25, *bracketing*
the truth 0.1655 — the signature of an unconverged approximant over a correct representation, not a
wrong one.

## THE SKEW LADDER WORKS — θ is finally pinned (2026-07-25, provisional)

After the conformal map and the direct extraction both failed *because they needed $(A,\theta)$ as
input*, the third route — a **second observable** — is not circular: $\kappa_3$'s late orders are
governed by the *same* Borel singularities with different amplitudes, so it adds independent
constraints on the same unknowns using only coefficients we can compute.

**The $\kappa_3$ ladder** ($\kappa_3=\eta^4\sum_j t_j\eta^{2j}$, `coupled-atlas/swtl_kappa3.py`,
grids $n=12,16,20,24$, extrapolated with the same screens as the variance ladder):

$$t_{0:5}=(+0.0169,\ +0.0566,\ +0.1689,\ +0.4045,\ +0.5182,\ -5.225)$$

It is **far cheaper than the variance ladder** — $t_j$ needs chaos index only $2+2j$ versus $v_k$'s
$2k+1$, so the whole 6-coefficient ladder costs ~2 min at $n=12$ against $v_7$'s ~10 days.

**Structural corroboration, before any fitting.** Within each catastrophe *both* observables first
turn negative at the **same rung** — cusp: $v_3$ and $t_3$; swallowtail: $v_5$ and $t_5$ — which is
exactly what a shared $\cos(k\theta-\varphi)$ envelope demands. The swallowtail's node is 2 rungs
later, so its $\theta$ must be **smaller**. This is the first independent confirmation that the
Borel geometry inferred from the variance ladder is physical and not a Padé artifact.

**The joint fit** (`swtl_joint.py`; shared $|\zeta|,\theta$, per-ladder $C,\varphi,\alpha$; 13 data
points, 6–8 parameters) is **stable in $n$**, which no previous method was:

| skew grid $n$ | 12 | 16 | 20 | 24 |
|---|---|---|---|---|
| $\theta$ | 37.2° | 37.9° | 38.0° | **38.0°** |
| $|\zeta|$ | 1.452 | 1.469 | 1.471 | 1.470 |

Four independent lines agree:

| method | $\theta$ |
|---|---|
| joint fit, raw | 37.4–38.0° |
| joint fit, bias-corrected (cusp gate runs +4.7°: 54.7 vs published 50) | **≈33°** |
| node position (flip at rung 5 vs cusp's 3, $\sim50°\times\tfrac46$) | 33.3° |
| low-$K$ Borel–Padé | 35–38° |

$$\boxed{\theta_{\rm swallowtail}\approx33\text{–}42^\circ\ (\text{best }\sim39^\circ),\qquad |\zeta|\approx1.5}$$

**and $\theta=50^\circ$ (the cusp value) is excluded** — the swallowtail's Borel phase is genuinely
smaller, not inherited. **[NUMERIC]**

**REVISED after $v_7$ landed (see below).** The first reading here was "33–38°", from the joint fit
plus a node-position argument. Two corrections: (i) the node argument ($\theta\propto1/k_{\rm node}$)
tacitly assumes the two catastrophes share the phase $\varphi$, which is not justified — a different
$\varphi$ moves the node at fixed $\theta$, so it is **not** a reliable estimator and is demoted to a
consistency check; (ii) with $v_7$, Borel–Padé $\theta$ stabilises at 50.6° ($K=6$: 50.7, $K=7$:
50.6), which at first looked like a flat contradiction of the joint fit's 38°.

**The methods reconcile under matched-$K$ differential calibration against the cusp** (true
$\theta=50°$). What is robust is the *offset*, not the absolute value — Borel–Padé reads $+8.8°$ high
on the cusp itself:

| method (matched $K=6$) | cusp | swallowtail | offset | ⇒ $\theta_{\rm swtl}$ |
|---|---|---|---|---|
| Borel–Padé poles | 58.8° | 50.7° | **−8.1°** | 41.9° |
| Darboux, variance-only | 49.8° | 39.2° | **−10.6°** | 39.4° |
| joint variance+skew | 54.7° | 37.6° | **−17.1°** | 32.9° |

Every method puts the swallowtail **8–17° below** the cusp; the sign and rough magnitude of that
differential are method-independent even though the absolute values disagree by 9°. So
**$\theta_{\rm swtl}<\theta_{\rm cusp}$ is the solid claim**, with the absolute value 33–42°.
The best-calibrated single estimator is Darboux variance-only (bias just $-0.2°$ on the cusp) →
**39.4°**.

⚠ **Honest limits.** (i) The joint fit's max relative residual is 0.8–0.9 (the cusp's is 0.7–1.1),
so the shared-$(|\zeta|,\theta)$ model describes the data only roughly; the *stability* in $n$, not
the fit quality, is what carries the result. (ii) The cusp gate passes but does **not tighten** —
variance-alone gives 49.8° at half the residual (0.129), so the skew ladder's contribution is
$n$-stability and the node argument, not extra precision. (iii) The bias correction is crude.

**Why the grid stops at $n=24$.** $n=28$ was launched and then deliberately abandoned. $\theta$ held
at 38.0° across $n{=}20\to24$ **while $t_5$ itself moved 5%** ($-5.65\to-5.39$) — so $\theta$ is
demonstrably insensitive to the top rung, which is the only thing a further grid point would sharpen.
The residual uncertainty is entirely systematic (model residual, the $+4.7°$ gate bias, the 2.4°
$\alpha$-spread) at $\pm3$–5°, none of which a finer grid touches. Continuing would also have starved
the $v_7$ run: load was 25–39 on 10 cores with the $\kappa_3$ pool holding 8 of them.

## Route 1 (read the non-perturbative sector off the observable) — TRIED AND REFUTED (2026-07-25)

Idea: since $v_8$ is walled, get $(A,\theta)$ from the *physics* instead of the ladder. The
trans-series predicts a remainder $R(x)=f_{\rm exact}-f_{\rm pert}$ that decays as
$e^{-A\cos\theta/x}$ and **oscillates** in $1/x$ with frequency $A\sin\theta$. Fit $R$, read off
$A,\theta$. Implemented in `coupled-atlas/swtl_nonpert.py`.

**Genuine by-product, worth keeping: a second-order FP solver.** The old first-order upwind scheme
carries numerical diffusion $\sim|{\rm drift}|\,dp/2$; a MUSCL/van-Leer TVD reconstruction removes
it. At $q=2,\eta=0.3$ against the converged series 0.144770:

| $dp$ | 0.04 | 0.02 | 0.01 |
|---|---|---|---|
| upwind | +58.78% | +29.64% | +15.08% |
| **MUSCL/van Leer** | **−0.01%** | **+0.11%** | **+0.31%** |

A ~5000× accuracy gain at the same resolution. (The mild *growth* of error with finer $dp$ indicates
a time-step/positivity-clipping floor, so Richardson in $dp$ is **not** valid here — use the direct
value at moderate $dp$.) This solver is reusable for any ground-truth work. **[NUMERIC ✓]**

**But the extraction itself fails its own cusp gate**, where $A=1.9,\theta=50^\circ$ are known:

| perturbative part subtracted | implied $A\cos\theta$ | true | sign pattern of $R$ |
|---|---|---|---|
| optimally truncated sum | 1.526 (25% off) | 1.221 | − − − − − |
| median Borel sum | 0.731 (40% off) | 1.221 | + + + + + + |

**The decisive tell is the absence of oscillation.** A pair at $50^\circ$ must flip the sign of $R$
across $x\in[0.35,1.1]$ (the phase $A\sin\theta/x$ sweeps ~2.5 rad); $R$ is monotone in both
variants. So $R$ is not the pair term — it is method error. With only 7 coefficients the optimal
truncation sits at $K=3$ for every $x$, so the "truncation ambiguity" is really the large, monotone,
omitted tail $v_4x^4+v_5x^5+v_6x^6$; and after Borel subtraction the cusp residual is 0.4–1.3% of
$f$, i.e. the same size as the Padé error and the solver error. For $q=3$ the residual is 5.6%→49%,
transparently the Padé failure rather than physics.

**The obstruction is circular:** extracting the non-perturbative data requires an accurate
perturbative sum, and an accurate perturbative sum requires the Borel data. **[NEGATIVE RESULT,
validated against the cusp]**

## Conformal-map resummation — TRIED AND REFUTED (2026-07-25)

Since $v_7$ is compute-walled, the natural move was a *method* upgrade needing no new coefficients:
conformal-map Borel resummation, the standard fix for evaluating at $x=2$ when the nearest Borel
singularity sits at $|\zeta|\approx1.1$. Implemented in `coupled-atlas/swtl_conformal.py` with the
conjugate-pair map
$$\sigma(u)=1-2\tfrac{u}{A}\cos\theta+\tfrac{u^2}{A^2},\qquad
w(u)=\frac{1-\sqrt\sigma}{1+\sqrt\sigma},\qquad
u(w)=A\Big[\cos\theta-\sqrt{\sigma(w)-\sin^2\theta}\Big],$$
which correctly sends both branch points to $|w|=1$ (verified: $\sigma(\zeta)=0$ exactly).

**It fails the cusp validation gate**, so it was never applied at $q=3$:

| $x$ | plain $[3/3]$ | conformal+Padé | conformal+truncated | truth |
|---|---|---|---|---|
| 1.0 | 0.2244 | 0.1935 | 0.1980 | 0.232 |
| **2.0** | **0.2347** | **0.1687** (29% off) | 0.3034 | **0.237** |
| 2.25 | 0.2337 | 0.1623 | 0.3724 | 0.234 |

**Why — and it is structural, not an implementation bug.** The map is built for *two* singularities,
but the Borel plane has *three*: the conjugate pair **and** the real FW instanton. A real positive
$z_r$ maps to
$|w|=0.13,\,0.12,\,0.08,\,0.01,\,0.22$ for $z_r=1.2,1.5,1.9,2.5,3.8$ — i.e. **strictly inside the
unit disc**, so the mapped series has a singularity in its own domain of convergence and cannot
converge there. This is precisely the "three singularities compete" situation the 𝒲 paper flags
(Rem. inheritance): for the swallowtail the pair ($|\zeta|\approx1.1$–1.5) and the real instanton
($z_r\approx1.24$–1.54) have *comparable moduli*, so no two-singularity map separates them.

**Subtracting the instanton first does not rescue it.** Scanning $(z_r,\gamma,C_r)$ over 63
combinations on the cusp, $f(2)$ swings over $[-0.025,+0.396]$ with **no stable plateau**; the best
value (0.2333) is only obtained by picking the combination nearest the known answer — fitting to the
answer, not a method. Untuned plain $[3/3]$ beats it (0.2347).

**The real lesson.** The conformal map needs the Borel singularity data $(A,\theta,z_r,C_r,\gamma)$
as *input*, and the seven-coefficient ladder does not determine it — $\theta$ is unidentifiable
(20.8° anchor spread at small residual), and $C_r$ is softer still. **So the method upgrade fails for
a data reason, not a numerical one: you cannot build the map without the very information the ladder
is too short to supply.** This closes the "improve the method instead of the engine" route.
**[NEGATIVE RESULT, validated against the cusp]**

### Current best diagnosis

With $v_5$ exact and $v_6$ good to $\sim1.7\times$, the surviving explanation is the **$v_4$
near-node**: $v_4=+0.036$ sits almost exactly on a zero of the $\cos(k\theta-\varphi)$ envelope, so
the ladder is $0.177,\,0.036,\,-1.496,\,-5.339$ — successive ratios $0.21,\,-41,\,3.6$. A $[3/3]$
Padé built on a sequence with a near-zero interior coefficient is badly conditioned, whereas the
cusp's smoothly-growing $-0.030,-0.451,-1.19,-1.90$ (ratios 15, 2.6, 1.6) is not. **This says the
swallowtail needs MORE RUNGS — $v_7,v_8$ — to see past the node, which is the opposite of the
"finer grids" prescription retracted above.** `_yexprs_15.txt` already supports $v_7$ (idx 15).
**[NUMERIC, best current reading]**

## Fine-grid campaign for $v_5,v_6$ — launch record (2026-07-24)

Acting on the diagnosis above (the limiter is grid convergence, not rung count),
`coupled-atlas/swtl_fine.py` extends both grids: **$v_5$ at $n=28,32,36$** and **$v_6$ at
$n=16,18,20$**, run **6-way parallel** (all `_ybase_q3_*.pkl` idx 1–13 are prebuilt, jobs are
single-threaded, and the two tracks use disjoint $n$ so the `_ckpt_{k}_{n}` caches cannot collide).
Log `swtl_fine.log`; results merged into `swtl_results.json` by the parent **only at the end**, so
mid-flight the JSON does not yet contain them.

All six completed in 54 min wall-clock. Per-job cost: $v_5$ 14.5/25.8/39.1 min at $n=28/32/36$;
$v_6$ 20.5/32.9/54.0 min at $n=16/18/20$ — consistent with the observed $\sim n^{4.8}$ scaling.
Results in the table at the top of this section.

## What reaches the full 𝒲 standard, and what doesn't

Reached: the operator identification, the exact skeleton (Bessel-1/5) and its Z₅ structure, the exact
$v_0$, reliable $v_1,v_2$, the factorial-divergence signature, and the closed-form real instanton.
**Not** reached: (i) the full 7-coefficient ladder $v_3,\dots,v_6$ — these are δ-sensitive and need the
production grid-extrapolation engine (`chaos_transfer.py`-class), as the cusp did; (ii) the end-to-end
median-Borel resummation to $\beta{=}2$ — with only four coefficients and a tighter radius (β=2 sits at
$x/A\approx3.6$, more non-perturbative than the cusp's ≈1.7), Borel–Padé scatters and cannot reconstruct
the ground truth $\operatorname{Var}=0.328$. Both are compute-bound extensions, not obstructions of
principle. See [[two-ladders-genuine-a4]], [[theorem-grade-constants-transcendental]].

> **STATUS UPDATE (2026-07-24).** (i) is now **DONE** — the 7-coefficient ladder landed. (ii) is
> **attempted and still negative**: with all seven coefficients the median resummation reaches
> $f(2)=0.126$ against the (now independently verified) truth 0.1655, still scattered. The
> diagnosis has sharpened though — the limiter is **grid convergence of $v_5,v_6$**, not the number
> of coefficients. See "CAMPAIGN RESULT" below.

## Production campaign (v3–v6) — set up and launched (2026-07-24)

The full 7-coefficient ladder is now running on the production Wick/transfer engine, adapted for q=3:
- **Driver**: `coupled-atlas/swtl_production.py` (monkeypatches the 3 q-points: backbone V=Y³ + recessive
  IC −Y0^{3/2}; node V-derivative table Vd0=Yst³, Vd1=3Yst², Vd2=6Yst, Vd3=6; q-tagged `_ybase_q3_*.pkl`
  caches). Reuses chaos_engine/chaos_diagram/chaos_transfer + `_v6_driver.vk_terms` + q-independent `_yexprs`.
- **VALIDATION GATE PASSED**: the same driver at q=2 reproduces the pristine cusp engine term-by-term,
  v3(n=24)=−0.043077 (identical: +0.33104, −0.06320, −0.12394, +0.00009). So q=3 is trustworthy.
- **Campaign**: `coupled-atlas/swtl_campaign.py` (bg PID 94766), grids v3:[24,30,36] v4:[24,30,36]
  v5:[16,20,24] v6:[10,12,14], Richardson n→∞ extrapolation, checkpoints to
  `coupled-atlas/swtl_results.json`, log `swtl_campaign.log`. v3/v4 in hours; v5/v6 grind (v6 MAXORD=7 was
  the cusp's wall). NOT harness-tracked — check `swtl_results.json`. Combine with the real instanton
  (s⁷/28) for the Borel/Stokes analysis once coefficients land.
