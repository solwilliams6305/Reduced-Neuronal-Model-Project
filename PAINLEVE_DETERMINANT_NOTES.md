# Notes — the Painlevé-IV determinant test (and the recurring cubic)

_June 2026. Figure `coupled-atlas/figures/painleve_determinant_test.png`; scripts
`painleve_determinant_test.py`, `subordination_probe.py`. Tags: [VALIDATED]/[REFUTE]/[SUPPORTED]/[OPEN];
construction status [CITED]/[RECONSTRUCTED]._

Goal: confirm or refute **cusp law = TW⁽²⁾₂ = a Painlevé-IV (parabolic-cylinder) Fredholm determinant**,
testing against the tighter fingerprint (cubic QQ-map, κ₅≈−2.2, κ₆≈−2.7), not just skew/kurtosis.

## Pipeline, validated [CITED + VALIDATED]

Bornemann–Nyström Fredholm determinant $F(s)=\det(1-K)_{L^2(s,\infty)}$. With the **Airy** kernel
$K_{\rm Ai}=\int_0^\infty\!\mathrm{Ai}(x{+}t)\mathrm{Ai}(y{+}t)\,dt$ it reproduces **TW₂**: $\mathrm{Ai}(0)=0.3547$
(exact 0.3550); determinant mean −1.775 (−1.771), std 0.906 (0.902), skew +0.207 (0.224), exkurt +0.117
(0.093); max kernel eigenvalue 0.998 ≤ 1. **Pipeline trustworthy.**

## The naive parabolic-cylinder determinant — REFUTED [RECONSTRUCTED → REFUTE]

I built the natural cusp analogue $K_C=\int_0^\infty\!Ci(x{+}t)Ci(y{+}t)\,dt$, $Ci$ = decaying solution of
the cusp turning $u''=\mathrm{sign}(y)y^2u$ (the parabolic-cylinder analogue of Ai; WKB-normalized,
**reconstructed, not cited**). It is **not a valid determinantal kernel**:
- **max eigenvalue ≫ 1** (≈ 3.7–15 across scales, vs Airy's ≤1) ⇒ $\det(1-K)$ leaves $[0,1]$;
- the resulting "$F(s)$" is **wildly non-monotone** (min ΔF = −1) — not a CDF.

**Why (structural).** $\int\!Ai\,Ai$ is the restriction of a *projection* because of the shift-completeness
$\int \mathrm{Ai}(z{+}t)\mathrm{Ai}(z{+}s)\,dz=\delta(t{-}s)$ — which holds for the **linear** (Airy) turning.
For the **quadratic** cusp turning the completeness fails, the kernel is not a projection, and the
eigenvalue bound is violated. So **the cusp law is NOT the soft-edge $\int\!\phi\phi$ determinant**, and
that determinant **cannot be reconstructed from the soft-edge analogy.**

## What this does and does not refute

- **REFUTED:** the *specific soft-edge-kernel* construction of TW⁽²⁾₂.
- **NOT refuted, and SUPPORTED:** that **Painlevé IV** governs the cusp via the *dynamical / Riemann–
  Hilbert* route. Positive evidence:
  1. the exact left tail $e^{-|s|^5/20}$ (Freidlin–Wentzell), the parabolic-cylinder/PIV signature;
  2. the **recurring cubic deformation** — both the TW→cusp QQ ($Q\approx Q_{\rm TW}+0.072Q^2-0.035Q^3$)
     *and* the β-family QQ are cubic (see below). PIV's nonlinearity is cubic ($\tfrac32 w^3+4xw^2$);
     a cubic quantile deformation is exactly what a PIV-vs-PII (Airy) relation would produce.
- **OPEN [needs CITED RH machinery]:** a genuine PIV determinant for the cusp must come from the
  Riemann–Hilbert / isomonodromy structure (the PIV τ-function), **not** the $\int\!\phi\phi$ kernel.
  That construction is research-grade and is not reproduced here.

## Subordination / TW⁽q⁾_β stress-test [NUMERIC]

The β-family (cusp, q=2): skew +0.25/+0.62/+0.86 and exkurt −0.80/−0.24/+0.66 for β=1/2/4. The QQ
between β-laws is a **cubic transform** (QQ β=4 vs β=1: cubic residual **0.017** vs affine 0.162). So the
whole β-family is a **one-parameter (cubic) transform family** — *subordination-consistent*: one "bare"
law smoothly reparametrized by β. It is **not** a pure location-scale subordination (the standardized
shape genuinely changes). The recurrence of the **same cubic** in both the q-direction (TW→cusp) and the
β-direction is the strongest structural hint that TW⁽q⁾_β is a low-order (PIV-type) deformation family.

## Verdict

| claim | status |
|---|---|
| Nyström pipeline; Airy → TW₂ | **[VALIDATED]** (CITED kernel) |
| cusp = soft-edge $\int\!Ci\,Ci$ determinant | **[REFUTE]** (not a valid kernel; projection property fails) |
| cusp governed by Painlevé IV (dynamical/RH route) | **[SUPPORTED, not proven]** (|s|⁵ tail + recurring cubic QQ) |
| TW⁽q⁾_β is a low-order (cubic) deformation family | **[NUMERIC]** (q- and β-QQ both cubic) |
| β-family = subordination of one bare law | **[NUMERIC: consistent]** (cubic transform, not pure scale) |
| genuine PIV τ-function determinant for the cusp | **[OPEN]** (RH machinery; not reconstructible from soft edge) |

**Bottom line.** The cusp law is *not* a soft-edge Fredholm determinant — that route is closed, cleanly,
for a structural reason (no projection property at a quadratic turning). But the Painlevé-IV identity is
*supported* by the |s|⁵ tail and the recurring cubic deformation in both the q- and β-directions; pinning
it requires the Riemann–Hilbert/τ-function construction, which is the standing open problem (rung F).
