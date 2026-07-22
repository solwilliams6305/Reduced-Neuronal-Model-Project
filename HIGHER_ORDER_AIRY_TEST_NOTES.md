# Depth — testing the reframe: is the cusp law 𝒲₂ the k=2 higher-order Tracy–Widom?

_June 2026. Figure `coupled-atlas/figures/hoTW_compare.png`; script `higher_order_airy.py`,
`hoTW_compare.py`. The honest test of the "multicritical / higher-order-TW" reframe proposed for the cusp
law. Result: **inconclusive — the reframe is neither confirmed nor cleanly refuted**, because the naive
higher-Airy kernel I built is the wrong object for k≥2. Tags [PROVED/validated]/[NUMERIC]/[OPEN]._

## What was built and validated

Generalized Airy $\mathrm{Ai}_k(x)=\frac1\pi\int_0^\infty\cos\!\big(\tfrac{t^{2k+1}}{2k+1}+xt\big)dt$; higher-order
Airy kernel $K_k(x,y)=\int_0^\infty\mathrm{Ai}_k(x{+}u)\mathrm{Ai}_k(y{+}u)\,du$; higher-order TW
$F_k(s)=\det(I-K_k)|_{L^2(s,\infty)}$ via Bornemann–Nyström (Gauss–Legendre, `np.linalg.slogdet`).

**k=1 validation [PASSED]:** $\mathrm{Ai}_1$ matches tabulated Airy values to 5 digits; the Fredholm
determinant reproduces **TW₂ (GUE)**: mean −1.76 (ref −1.771), std 0.903 (0.902), skew **+0.225** (+0.224),
exkurt **+0.090** (+0.093), left-tail exponent **2.97** (=3). The machinery is correct.

## The k=2 test [INCONCLUSIVE — wrong kernel]

The naive k=2 kernel (cosine-integral $\mathrm{Ai}_2$, phase $t^5/5$) gives:
- skew **−0.00** (symmetric), exkurt −0.05, std 1.11, left-tail exponent **2.42**.

This matches **neither** the cusp $\mathcal W_2$ (skew **+0.61**, exk −0.24) **nor** the expected order-2
multicritical edge (whose left exponent should be **5**, steeper than TW's 3 — but the fit gives 2.4,
*shallower* than k=1). So the naive object is symmetric, wide, and not even a multicritical edge.

**Diagnosis [the honest reason].** The k=1 cosine integral is the unique recessive Airy solution, but the
k=2 "higher-Airy" ODE $y''''=-xy$ has a **two-dimensional decaying space** as $x\to+\infty$. The real
cosine integral is the *wrong* combination — it is not the recessive parametrix that defines the
multicritical edge kernel. The correct k≥2 kernel requires the proper Riemann–Hilbert parametrix (the
Stokes structure of the $t^5$ phase has multiple rays) or, equivalently, the Painlevé-II-hierarchy
representation $F_2(s)=\exp(-\int_s^\infty(x-s)\,q_2(x)^2dx)$ with $q_2$ the 2nd hierarchy member. Building
that is a substantially harder special-function computation, not a one-line generalization.

## Honest conclusion

- **The reframe is UNRESOLVED.** This computation does **not** confirm that $\mathcal W_2$ is the k=2
  higher-order TW, and does **not** cleanly refute it either — because the kernel I could build naively is
  demonstrably the wrong object (it fails to even be an order-2 multicritical edge). No dressing-up: the
  elegant reframe remains a conjecture; the direct test was thwarted by the special-function subtlety.
- **What stands:** the Fredholm-determinant machinery is built and validated (k=1=TW₂), so it is ready to
  *correctly* test the reframe once the proper higher-Airy recessive solution / PII-hierarchy $q_2$ is in
  hand. That is the concrete next step.
- **$\mathcal W_2$'s analytic identity remains OPEN/CONJECTURAL**, exactly as the ledger stated. The cusp
  law is still genuinely new (consistent with the falsifier campaign); we have not yet found its closed
  form, and the most natural candidate (naive higher-Airy) is now known to be the wrong target.

## Status

| item | status |
|---|---|
| higher-order Airy Fredholm-determinant machinery | **[VALIDATED]** (k=1 = TW₂, incl. tail exp 3) |
| naive k=2 kernel = cusp 𝒲₂ ? | **NO** — symmetric, not multicritical (wrong solution) |
| reframe (𝒲₂ = k=2 higher-order TW) | **[OPEN]** — needs the correct higher-Airy recessive / PII-hierarchy $q_2$ |
| 𝒲₂ analytic closed form | **[CONJECTURAL]** — unchanged |
