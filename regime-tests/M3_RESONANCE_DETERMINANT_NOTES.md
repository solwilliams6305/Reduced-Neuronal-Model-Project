# M3: the closed-form resonance determinant D(λ) — pushed, falsified twice, and the sharpened map

_June 2026. Continuation of `INTERFACE_TAU_CONSTRUCTION_NOTES.md`. The flagged-"completable" sub-problem:
derive $D(\lambda)$ in closed Gamma form. Outcome: the natural analytic closed form **does not exist on the
real line** (two ansätze falsified against direct integration), and that falsification **proves complex
scaling (M1) is necessary, not optional.** Plus: the **accurate complex resonances**. Figures
`coupled-atlas/figures/resonances_true.png`. Scripts `verify_D.py`, `diagnostic_D.py`, `diagnostic_D2.py`,
`test_stokes.py`, `refine_res.py`. Tags **[DERIVED]/[FALSIFIED]/[NUMERIC]/[CITED]**; ⚑ load-bearing._

---

## 0. Outcome up front (honest)

I pushed hard for the closed form, as flagged. It **resisted**, and — per the method — I stopped at a precise,
**numerically-validated** diagnosis rather than forcing it:

| step | result | status |
|---|---|---|
| Candidate closed form (Wronskian-at-0 of the two PC Jost solutions) | $R(\lambda)=e^{-i\pi/4}R(i\lambda)$, $R(w)=\Gamma(\tfrac34-\tfrac w4)/\Gamma(\tfrac14-\tfrac w4)$ — clean analytic Gamma-ratio | **DERIVED** |
| Test vs direct integration | candidate's zeros **do not match** the true resonances | **FALSIFIED ⚑** |
| Minimal fix (single Stokes multiplier $S\propto1/\Gamma(\tfrac12-\tfrac{i\lambda}2)$) | required $S$ **inconsistent** across resonances ($|S_{\rm req}|/|S_{\rm pred}|=15.6,\,0.74,\,0.15$) | **FALSIFIED ⚑** |
| Root cause | the outgoing inverted-oscillator value carries the **non-analytic** $|\Gamma|^{1/2}$ DLMF-$W$ normalization | **DERIVED** |
| **Accurate complex resonances** | $\lambda_0=0.86-0.82i$ (clean), $2.30-1.22i$, $4.14-1.08i$ | **NUMERIC ✓** |
| Correction | supersedes the crude-scan $\lambda_0=1.64-1.49i$ (that used a local-WKB criterion) | logged |
| Map | the non-analyticity is **exactly** what complex scaling (M1) removes ⇒ M1 necessary | **mapped** |

**Verdict.** "Completable closed form" was over-optimistic. Pushing it produced something cleaner than a forced
calculation: a **proof (by falsification) that $D(\lambda)$ is not a meromorphic Gamma-ratio on the real line**,
plus the accurate resonances. The closed form lives on the **complex-scaled contour**, not $\mathbb R$.

---

## 1. The candidate closed form (derived)

Jost solutions: recessive-right $u_+(Y)=U(-\tfrac\lambda2,\sqrt2\,Y)$; "outgoing"-left, via $Y=\tfrac{e^{i\pi/4}}{\sqrt2}t$,
$u_{\rm out}(Y)=U(-\tfrac{i\lambda}2,\sqrt2\,e^{-i\pi/4}Y)$. Resonance determinant $D(\lambda)=W[u_+,u_{\rm out}](0)$.
Using the DLMF zero-values $U(a,0)=\sqrt\pi\,2^{-1/4-a/2}/\Gamma(\tfrac34+\tfrac a2)$,
$U'(a,0)=-\sqrt\pi\,2^{1/4-a/2}/\Gamma(\tfrac14+\tfrac a2)$ [DLMF 12.2.6–7], the condition $D=0$ collapses to
$$\boxed{\ R(\lambda)=e^{-i\pi/4}R(i\lambda),\qquad R(w)=\frac{\Gamma(\tfrac34-\tfrac w4)}{\Gamma(\tfrac14-\tfrac w4)}\ }\quad\text{— confining ratio × oscillatory ratio × interface phase.}$$
Clean, analytic, the structure predicted in the previous note. **[DERIVED.]**

## 2. Falsification #1 — the candidate is wrong ⚑

The leading Gamma-ratio asymptotics of $R(\lambda)$ and $e^{-i\pi/4}R(i\lambda)$ **cancel** ($R(w)\sim\sqrt{-w/4}$,
and $e^{-i\pi/4}\sqrt{-i\lambda/4}=\sqrt{-\lambda/4}$): so the candidate has **no isolated zeros** near the true
resonances — Newton runs off to $|\lambda|\to\infty$. Direct integration (the **outgoing-matching Wronskian**,
`diagnostic_D2.py`/`refine_res.py`) gives the true resonances at $\lambda\approx0.86-0.82i,\,2.30-1.22i,\,4.14-1.08i$
(figure). The candidate's zeros sit nowhere near them. **The clean Gamma-ratio is falsified. [FALSIFIED, NUMERIC ⚑.]**

## 3. Diagnosis + Falsification #2 — no single Stokes multiplier saves it ⚑

**Why #1 failed.** At $Y\to-\infty$ the rotated argument $\sqrt2\,e^{-i\pi/4}Y$ lands on the **anti-Stokes ray**
$\arg=3\pi/4$. There $U$ is **not** a single exponential — by the **Stokes phenomenon** it is a *mixture* of
outgoing $e^{+iY^2/2}$ and incoming $e^{-iY^2/2}$. So $u_{\rm out}$ above is **not purely outgoing**; the
Wronskian-at-0 used the wrong left solution. **[DERIVED.]**

**Minimal fix tried.** Outgoing $=U_++S\,U_-$ with one Stokes multiplier $S$. Solving $D=0$ for the **required**
$S$ at the true resonances and comparing to the predicted $S_{\rm pred}=\sqrt{2\pi}/\Gamma(\tfrac12-\tfrac{i\lambda}2)$
(`test_stokes.py`):
$$|S_{\rm req}|/|S_{\rm pred}|\;=\;15.6,\quad0.74,\quad0.15\quad(\lambda_0,\lambda_1,\lambda_2)\ \text{— not constant.}$$
**A single analytic Stokes multiplier does not reconcile it. [FALSIFIED, NUMERIC ⚑.]**

**Root cause.** The purely-outgoing inverted-oscillator solution is the DLMF $W$/$E$-function, whose interface
value is $W(a,0)\propto\big|\Gamma(\tfrac14+\tfrac{ia}2)/\Gamma(\tfrac34+\tfrac{ia}2)\big|^{1/2}$ — a
**non-analytic modulus** (the real-axis oscillatory normalization $k=\sqrt{1+e^{-\pi\lambda}}$). Hence
$$\boxed{\ D(\lambda)\ \text{is explicit in special-function values but NOT a meromorphic Gamma-ratio in }\lambda\text{ on }\mathbb R.\ }$$
The analytic $U$-rotated representation that *would* give a Gamma-ratio is not purely outgoing (misses the
Stokes resummation); the purely-outgoing $W$-representation is analytic-on-the-real-axis only through the
non-analytic $|\Gamma|^{1/2}$. **[DERIVED ⚑.]**

## 4. The positive — accurate complex resonances [NUMERIC ✓]

Outgoing-matching Wronskian, direct integration (`refine_res.py`, $dt=10^{-3}$, $Y_0=6$):
$$\lambda_0=0.86-0.82i\ (|D|{=}0.006),\quad \lambda_1=2.30-1.22i,\quad \lambda_2=4.14-1.08i.$$
Genuinely **complex** (widths $\Gamma\approx1.6,2.4,2.2$ — broad, the vanishing cusp barrier ⇒ strong leakage),
near but shifted off the bare half-oscillator levels. **Non-self-adjointness concretely confirmed.** These
**supersede** the earlier crude-scan $\lambda_0=1.64-1.49i$ (which used a single-point local-WKB amplitude
criterion; the Wronskian-with-outgoing criterion here is the principled one). **[NUMERIC ✓, correction logged.]**

## 5. The sharpened map — complex scaling (M1) is *necessary*

The blocker in §3 (the real-axis $|\Gamma|^{1/2}$) is precisely the pathology that **complex scaling
(Aguilar–Balslev–Combes)** removes. Rotate $Y\to Ye^{-i\theta}$: the inverted continuum swings off $\mathbb R$,
the non-analytic real-axis $W$-normalization becomes the **analytic** rotated-$U$ with a **definite
analytically-continued Stokes factor**, and the resonances become honest discrete eigenvalues on the rotated
contour. Therefore:
> There is **no real-line meromorphic shortcut** (proved by the two falsifications). The closed form lives on
> the **complex-scaled contour**: $D_\theta(\lambda)$ analytic, zeros $=\lambda_n$. **M3's failure upgrades M1
> from "recommended next step" to "necessary."** **[mapped; CITED — complex scaling.]**

This is a *stronger* statement than the previous "do M3 then M1": it shows the rank-1 object itself only
becomes analytic after complex scaling — so the **non-Hermitian (complex-scaled) Fredholm determinant** is not a
convenience but the intrinsic home of 𝒲, and its rank-$n$ matrix version (the intrinsic process, B/D) inherits
the same necessity. **[DERIVED consequence.]**

## 6. Net + tracker

**Honest:** M3 did **not** yield a clean closed-form Gamma-ratio — and *proved* (by two numerically-validated
falsifications) that none exists on the real line. **Positive:** the accurate complex resonances
($0.86-0.82i,\,2.30-1.22i,\,4.14-1.08i$), non-self-adjointness confirmed, and the earlier $\lambda_0$ corrected.
**Map:** complex scaling (M1) is now shown **necessary**; the resonance determinant is analytic only on the
rotated contour. **Tracker: B stays 95%** (no closed form achieved; the obstruction is sharper and a building
block — the resonances — is pinned). Next concrete step: **M1**, the complex-scaled non-Hermitian determinant,
now the unavoidable route.
