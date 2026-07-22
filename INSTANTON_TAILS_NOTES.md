# R1 — Instanton / Freidlin–Wentzell tails of 𝒲₂: the cusp is genuinely new

_June 2026. Figure `coupled-atlas/figures/ldp_tails.png`; script `ldp_tails.py`. The decisive route from
`NOVEL_IDENTIFICATION_ROUTES.md`: compute the large-deviation tails of the cusp escape law and settle
"genuinely new" vs "half-integer (k=½) multicritical." Tags [ANALYTIC]/[NUMERIC]/[CONCLUSION]._

## The two tails, derived analytically

Cusp Riccati $dp=(V-p^2)d\tau+\eta dW$, $V=\operatorname{sign}(Y)|Y|^q$, escape location $Y^\star$;
orientation $s=+(Y^\star-\text{mean})/\text{std}$.

**Right tail (early escape, $s\to+\infty$) — Kramers action over the proven cubic barrier.**
For $Y_e>0$ the canard sits at $\bar p=\sqrt V=Y_e^{q/2}$; the fluctuation has the bistable cubic potential
$\Phi(\delta p)=\bar p\,\delta p^2+\tfrac13\delta p^3$ with barrier $\Phi(-2\bar p)=\tfrac43\bar p^3$
(rigorously established in `NO_EARLY_ESCAPE_NOTES.md`). The escape (Kramers/FW) action is twice the barrier:
$$-\ln P(Y^\star>Y_e)\ \xrightarrow{\ Y_e\gg1\ }\ \frac{2}{\eta^2}\cdot\tfrac43\bar p^3=\frac{8}{3\eta^2}\,Y_e^{3q/2}
\quad\Rightarrow\quad \textbf{right exponent }=\tfrac{3q}{2}.$$
[ANALYTIC — standard 1D escape-action over the proven cubic potential.]

**Left tail (late escape, $s\to-\infty$) — holding instanton.**
Past the deterministic blow-up ($Y<Y_{\rm det}$, where $V=-|Y|^q$), survival requires the noise to oppose
the runaway. The optimal static hold ($p\approx0$) costs FW action
$$-\ln P(Y^\star<-Y_L)\ \approx\ \frac{1}{2\eta^2}\int^{Y_L}|Y|^{2q}\,dY=\frac{Y_L^{2q+1}}{2(2q+1)\eta^2}
\quad\Rightarrow\quad \textbf{left exponent }=2q+1.$$
[ANALYTIC — the $p=0$ hold is the static optimum; it is an *upper bound* on the true instanton action, so
the coefficient $1/(2(2q+1))$ is an overestimate, but the exponent $2q+1$ is robust.]

**Validation at q=1 (fold):** right $3q/2=3/2$, left $2q+1=3$ — **exactly the Tracy–Widom tail exponents**
($e^{-\frac{4}{3}s^{3/2}}$ right, $e^{-\frac{1}{12}|s|^3}$ left). The derivations reproduce TW. ✓

## The decisive conclusion

| | left exp | right exp | left/right ratio |
|---|---|---|---|
| **cusp 𝒲₂ (q=2)** | **5** = 2q+1 | **3** = 3q/2 | **5/3** |
| integer higher-order TW (k) | 4k+3 | (4k+3)/2 | **2** (all k) |
| half-integer k=½ candidate | 5 | **5/2** | 2 |
| fold/TW (q=1) | 3 | 3/2 | 2 |

- The cusp right exponent is **3, not 5/2** ⇒ **𝒲₂ is NOT the half-integer (k=½) multicritical edge.**
- The cusp tail ratio is **5/3, not 2** ⇒ **𝒲₂ is NOT any integer higher-order TW.**
- The fold (q=1) is the *only* member shared by the two families (ratio 2). For q≥2 the catastrophe ladder
  is a **distinct family** with q-dependent tail ratio $(4q+2)/(3q)$.

**𝒲₂ is genuinely new** — definitively, on the tail exponents. This closes the "is it multicritical TW?"
question (no, on both integer and half-integer counts).

## Numerical support (honest scope)

Simulated the cusp escape at $\eta\in\{0.55,0.7,0.9,1.15,1.45\}$ (~4×10⁵ each). The $\eta^2$-rescaled rate
functions $\eta^2(-\ln P)$ **approximately collapse** (improving as $\eta\to0$) — confirming the FW/LDP
$1/\eta^2$ structure. **But** the *reachable* range ($Y_e\lesssim0.4$, $Y_L\lesssim3.3$) is **near-typical**,
not the asymptotic-tail regime where the exponents 3, 5 live (those need $Y_e\gg1$). So:
- **[NUMERIC]** the LDP rate-function structure (η²-collapse) is confirmed; the right rate function is convex
  (consistent with the cubic onset), the left is steep (consistent with a high exponent, with the $p=0$-hold
  coefficient overshooting as predicted).
- **[ANALYTIC]** the exponents (3, 5) themselves rest on the Kramers-barrier / holding-instanton derivations
  above — which reproduce TW exactly at q=1. The Monte-Carlo cannot reach the asymptotic tails (rare events),
  so it supports the framework, not the exponent values.

## Status and next

| item | status |
|---|---|
| right tail exponent 3q/2 (cusp 3) | **ANALYTIC** (Kramers cubic barrier; TW-validated at q=1) |
| left tail exponent 2q+1 (cusp 5) | **ANALYTIC** (holding instanton; TW-validated) |
| 𝒲₂ ≠ integer higher-order TW (ratio 5/3≠2) | **CONCLUSION** |
| 𝒲₂ ≠ half-integer k=½ edge (right 3≠5/2) | **CONCLUSION** |
| FW/LDP η²-collapse | **NUMERIC** (confirmed, near-typical range) |
| left-tail coefficient (exact) | **OPEN** (needs the true instanton BVP, not the p=0 bound) |
| full closed-form identity of 𝒲₂ | **OPEN** → routes R3 (Painlevé-IV σ-form) and R2 (dynamic-pitchfork scaling function) |

**Net.** R1 decisively establishes what 𝒲₂ is **not** (no multicritical TW, integer or half-integer) and
pins its tail exponents (5, 3) analytically, reproducing TW at the fold. The cusp escape law is a genuinely
new universality class with tail ratio $(4q+2)/(3q)$. Its positive *identification* (a closed form) now
points to **Painlevé-IV** (the Weber/parabolic-cylinder structure — R3) and the **dynamic-pitchfork
weak-noise scaling function** (R2), which the tail data above will constrain.
