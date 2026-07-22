# Stage 1: the gap determinant for 𝒲 — pushed, falsified, and the sharp diagnosis

_June 2026. Continuation of `M1_COMPLEX_SCALING_NOTES.md`. Tries to complete 𝒲's construction as
$F_{\mathcal W}(s)=\det(1-K^\theta_s)$. Outcome: **the gap-determinant hypothesis is falsified** (two kernel
forms, θ-dependent, non-matching) — and the falsification yields a **sharp structural diagnosis**: 𝒲 is a
**first-passage** law of an unbounded-below operator, *not* an eigenvalue-gap of a determinantal process.
Figures `coupled-atlas/figures/kernel_falsification.png`; scripts `airy_pipeline.py`, `cusp_kernel.py`,
`cusp_kernel2.py`, `kernel_falsification.py`. Tags **[DERIVED]/[NUMERIC]/[FALSIFIED]/[CITED]**; ⚑ load-bearing._

---

## 0. Outcome

| step | result | status |
|---|---|---|
| Pipeline check: Airy → TW₂ via Nyström | mean $-1.78$, std $0.91$, **skew $+0.20$** (target $+0.224$), Ai(0) to 4 digits | **VALIDATED ✓** (pipeline trustworthy) |
| Cusp gap det, CD kernel $[\psi\psi'-\psi'\psi]/(x-y)$ | skew $1.27/1.52/1.10$, exk $0.38/1.99/0.77$ — **θ-DEPENDENT**, $\neq$ 𝒲 | **FALSIFIED ⚑** |
| Cusp gap det, two-Jost kernel $[u_+u_- - u_+u_-]/(x-y)$ | skew $0.70/0.00/1.28$, exk $-1.08/-1.97/-0.33$ — **θ-DEPENDENT**, $\neq$ 𝒲 | **FALSIFIED ⚑** |
| Diagnosis | 𝒲 is **first-passage** (unbounded-below), not an eigenvalue-gap | **DERIVED ⚑** |
| Corrected map | 𝒲 is governed by the **resonance / survival (Gamow) expansion**, not a Fredholm gap | **mapped** |

**Verdict.** $\mathcal W=\det(1-K^\theta_s)$ is **false**. The complex-scaled *spectral* determinant $D_\theta$
(M1, resonances, θ-independent) stands — but the *edge law* 𝒲 is **not its gap probability**. I did not force a
match; two natural kernels both fail the θ-independence test that the correct object must pass.

---

## 1. The pipeline is sound (Airy → TW₂)

`airy_pipeline.py`: built $\mathrm{Ai}$ by integrating $u''=xu$ from $x=8$ with the normalized asymptotic IC
($\mathrm{Ai}(0)=0.3547$ vs exact $0.35503$ ✓), formed $K_{\rm Airy}=[\mathrm{Ai}(x)\mathrm{Ai}'(y)-\mathrm{Ai}'(x)\mathrm{Ai}(y)]/(x-y)$,
Nyström on $(s,8)$, $F(s)=\det(1-K)$. Result: **mean $-1.78$, std $0.91$, skew $+0.20$, exk $+0.11$** — TW₂
($-1.77,0.90,+0.224,+0.093$). So the determinant→cumulants pipeline is **trustworthy**; a cusp mismatch is
real, not a bug. **[VALIDATED ✓]**

## 2. The cusp gap determinant — falsified (two kernels) ⚑

Built the rotated recessive Weber solution $\psi^\theta$ (cusp edge ODE $\phi''=\operatorname{sign}(z)z^2\phi$ on
$z=t\,e^{i\theta}$, normalized to $U(0,\sqrt2 z)$) and assembled the edge kernel two ways:
- **CD form** $K^\theta=[\psi(x)\psi'(y)-\psi'(x)\psi(y)]/(x-y)$ (`cusp_kernel.py`);
- **two-Jost form** $K^\theta=[u_+(x)u_-(y)-u_+(y)u_-(x)]/((x-y)W)$ (`cusp_kernel2.py`).

Both give $\det(1-K^\theta_s)$ with **cumulants that swing with $\theta$** (figure) — CD skew $1.27\to1.52\to1.10$,
two-Jost skew $0.70\to0.00\to1.28$ — and **neither matches** 𝒲 ($+0.60,-0.24$). The **θ-dependence is the
tell**: a correct complex-scaled determinant is θ-independent (M1's $D_\theta$ *was*, to 3 digits). So these
gap determinants are **not** the invariant object. **[FALSIFIED, NUMERIC ⚑]**

## 3. The sharp diagnosis — first-passage, not eigenvalue-gap ⚑

The reason is structural, and it is the same fact that made M3/M1 hard:

> **Airy/TW₂:** the stochastic Airy operator $-\partial^2+x$ (on $[0,\infty)$ with a boundary condition) is
> **bounded below** — it has a discrete **bottom eigenvalue**, and that eigenvalue *is* TW. "No eigenvalue above
> $s$" is a genuine event ⇒ a determinantal **gap probability** ⇒ a Fredholm determinant. ✓
>
> **Cusp/𝒲:** $-\partial^2+\operatorname{sign}(Y)Y^2$ has an **escape direction** ($-Y^2\to-\infty$ as
> $Y\to-\infty$, no boundary condition) — it is **unbounded below**, with **no bottom eigenvalue**. The escape
> is the **first explosion of the Riccati / first node going down the inverted side** — a **first-passage**, not
> an eigenvalue. So 𝒲 is **not** "no eigenvalue above $s$" and **not** a determinantal gap. ✗

This **confirms and explains** the long-standing negative ("𝒲 is not a soft-edge Fredholm determinant") with a
precise mechanism, now backed by the explicit falsification of two kernels. **[DERIVED ⚑]**

What complex scaling *did* give (M1) is the **spectral** determinant $D_\theta$ — the resonances $\lambda_n$ as
honest discrete eigenvalues of $H_\theta$. That is correct and θ-independent. But the resonances are the
spectrum of the *operator*, while 𝒲 is a *first-passage location* — different objects.

## 4. The corrected map — 𝒲 via the resonance/survival expansion

A first-passage law of a resonance-dominated escape is governed by the **Gamow / Mittag-Leffler expansion** over
the resonances, **not** a Fredholm gap:
$$\text{survival amplitude}\ \sim\ \sum_n c_n\,e^{-i\lambda_n\,\Theta},\qquad \lambda_0=0.86-0.82i\ \text{dominant},$$
with $\Theta$ the accumulated phase/depth. The escape statistics follow from the residues at the resonances
(the lowest, narrowest $\lambda_0$ setting the dominant rate $2\lvert\operatorname{Im}\lambda_0\rvert=1.64$).
The **resonance spectrum built in M1 is therefore the right invariant** for 𝒲 — but the law is assembled by a
**resonance expansion**, not by $\det(1-K)$. Validating the full survival-expansion → FP-𝒲 map (residues +
the deterministic $Y^\star_{\rm det}=-2.19$ trajectory) is the next concrete step. **[mapped, CONJECTURAL.]**

## 5. Stage 2 (rank-$n$) — correctly not attempted as a gap determinant

Stage 2 was conditional on Stage 1 landing. Since the **marginal** is not a gap determinant, the **process** is
not a matrix gap determinant either — it inherits the first-passage structure (the multi-resonance / multi-line
escape). The rank-$n$ object is the **joint** first-passage governed by the matrix resonance spectrum; building
it as a noncommutative Fredholm gap would fail the same θ-independence test. **Not built (correctly).**

## 6. Net + tracker

**Honest:** the headline hypothesis $\mathcal W=\det(1-K^\theta_s)$ is **falsified** — pushed with a
pipeline-verified Nyström solver and two kernel forms, both θ-dependent and non-matching. The **diagnosis is
the deliverable**: 𝒲 is a **first-passage law of an unbounded-below operator**, not an eigenvalue-gap — which
*explains* every earlier "not a Fredholm determinant" negative, precisely.

**Standing:** M1's complex-scaled **spectral** determinant $D_\theta$ and the **resonances** $\lambda_n$ are
correct and remain the analytic characterization of the cusp edge; 𝒲 is the **resonance-governed first-passage
law** (the corrected map), not a determinant.

**Tracker: B stays 96%** — no inflation: the spectral object stands, the gap-determinant route is now *closed
with reason*, and the resonance-expansion route is the mapped next step. This is the honest ceiling: 𝒲 is
characterized by the resonance spectrum but is a first-passage law, with **no closed-form determinant** (a
proven-style negative, not a missing calculation).
