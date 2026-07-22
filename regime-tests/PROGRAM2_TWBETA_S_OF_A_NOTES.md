# PROGRAM2 — TW_β halo: the a-dependence of the Stokes constant S(a) is ALGEBRAIC

**Date:** 2026-07-22. **Executes:** `PROGRAM2_HANDOFF_C_TW_BETA_HALO.md` open task **#2**, the
remaining half — the tail-specific *a-dependence* of `S(a)` beyond the a-independent backbone
`C=√(2/(3π³))` (fixed earlier today, see `PROGRAM2_TWBETA_STOKES_CONSTANT_NOTES.md`).
**Memory:** `program2-tw-beta-halo.md`.

---

## Result (one line)

The a-dependence of the tail Stokes constant `S(a)` carries **no new transcendental content**: under
an exact scaling symmetry of the backward-Kolmogorov PDE, `a` enters the fluctuation expansion only
through the single variable `g := β a^{3/2}`, so `S(a) = S₀ · a^{p}` with `S₀` an a-independent pure
number and `p` elementary. For the dominant boundary-escape (reduced) problem the constant is, in
**exact closed form,**

$$\boxed{\;S_0=\tfrac1\pi\;=\;0.3183098861837907\ldots\;,\qquad d_m\sim-\tfrac1\pi\,\Gamma(m)\,\bigl(\tfrac32\bigr)^{m}\;(\gamma'=0,\ \rho=\tfrac32\ \text{exact}).\;}$$

This settles the qualitative question the handoff posed — the "a-dependence" is eikonal/Jacobian, not a
new special function — and it is the natural completion of the "no α(β) / eikonal prefactor" theme.

---

## 1. The exact scaling symmetry (the mechanism) — verified symbolically

The tail is the explosion probability of the Riccati SDE `dp=((x+a)−p²)dx+2ε dB`; its no-explosion
probability solves the exact backward-Kolmogorov PDE (O1 note eq. 66)
`u_x + ((x+a)−p²)u_p + 2ε² u_pp = 0`. Under

$$p=\sqrt a\,P,\qquad x+a=a\,X\ \ (X\in[1,\infty)),$$

`sympy` confirms the PDE becomes (after dividing by `a`)
`a^{-3/2}u_X + (X−P²)u_P + 2ε² a^{-3/2}u_{PP}=0`, i.e. the transport part is **a-independent** and the
diffusion coefficient is `2ε_eff²` with

$$\varepsilon_{\rm eff}^2=\varepsilon^2 a^{-3/2}=\frac1{\beta a^{3/2}}=\frac1g,\qquad g:=\beta a^{3/2}.$$

**Consequences (rigorous).** The reduced fluctuation series depends on `a` only through `g`, so its
`β^{-n}` coefficient scales as `c_n(a) ∝ a^{-3n/2}` (up to a single n-independent prefactor; the
`β^{-n}` must pair with `a^{-3n/2}` to build `g^{-n}`), the numerical part being a-independent. Hence
the 1/β-Borel singularity sits at the instanton action `Φ(a)=(2/3)a^{3/2}` (equivalently the reduced
1/g-series has its Borel singularity at `2/3`), and the a-dependence of `S(a)` is purely the elementary
Jacobian: `S(a)=S₀·a^{p}`, **no new transcendental structure in `a`.**

## 2. The reduced escape rate in closed form

The dominant contribution is boundary escape at `x=0` over the cubic barrier `U(p)=p³/3−ap`.
Rescaling `p=√a P` gives the a-independent SDE `dP=−(P²−1)dτ+2ε_eff dW`, `V(P)=P³/3−P` (well at
`P=+1`, barrier at `P=−1`, `ΔV=4/3`). The mean-first-passage double integral factorizes to all orders
in `1/g` (the two Laplace peaks at `P=∓1` decouple), giving the **exact reduced escape rate**

$$R(g)=\frac1\pi\,e^{-\frac23 g}\,\frac1{\tilde Q(2/g)^{2}},\qquad
\tilde Q(D)=\sum_{j\ge0}\tilde q_j D^{j},\quad \tilde q_j=\frac{(6j)!}{576^{\,j}\,(3j)!\,(2j)!}\ \ (\text{exact}),$$

with `q̃_1=5/48`, and (Stirling) `q̃_j ~ (2π)^{-1}Γ(j)(3/4)^j` — factorial, Borel singularity at `D=4/3`
(⇔ the two-instanton `e^{-(4/3)g}=e^{-2Φ}`). The physical escape rate is `r(a,ε)=√a·R(g)` (exact by the
same substitution in the MFPT integral), exhibiting the `√a` Jacobian explicitly.

## 3. The reduced Stokes constant S₀ = 1/π — certified

The reduced fluctuation series `F(g)=Q̃(2/g)^{-2}=Σ_m d_m g^{-m}` (so `R=(1/π)e^{-(2/3)g}F`) has
`d_0=1, d_1=−5/12, d_2=−0.5382, …`, **all `d_m<0`** (non-sign-alternating — the same
positive-axis-Borel / median-resummation signature as Hastings–McLeod). Computed to `m=170` (exact
`q̃_j` via running product; `Q̃^{-2}` by series reciprocal-square) and Richardson-extrapolated:

| quantity | value | status |
|---|---|---|
| base `ρ` | `3/2` | exact (all windows) |
| index `γ'` | `0` | exact (~20 digits) |
| **`S₀`** | **`1/π = 0.31830988618379067…`** | **`S₀−1/π = −4·10⁻²⁵` (25 digits)** |

So `d_m ∼ −(1/π)\,Γ(m)(3/2)^m = −(1/π)(m−1)!\,(3/2)^m`.

## 4. Independent validation — direct MFPT integral

The closed-form `R(g)` was checked against a direct high-precision mean-first-passage double integral
`τ=(1/D)∫_{-∞}^{1}dy\,e^{V(y)/D}∫_y^{∞}dz\,e^{-V(z)/D}`, `D=2/g` (optimal-truncated `Q̃`):

| g | R_direct (MFPT integral) | R_formula (closed form) | rel. diff |
|---|---|---|---|
| 10 | 0.00038489754 | 0.00038455511 | 8.9×10⁻⁴ |
| 16 | 7.2060945×10⁻⁶ | 7.2059764×10⁻⁶ | 1.6×10⁻⁵ |

Agreement sharpens with `g` (8.9×10⁻⁴ → 1.6×10⁻⁵), as it must for an asymptotic formula (optimal
truncation of the divergent `Q̃` series). Confirms `R(g)` and hence the reduced `Q̃`/`d_m`/`S₀=1/π`.

## Relation to the backbone C, and what this means

The **full-tail** transcendental backbone (by Theorem T1, the rigorous Lax-pair map to PII) is
`C=√(2/(3π³))`. The **reduced** (frozen boundary-escape) constant is `1/π`. They factor as

$$C=\sqrt{\tfrac{2}{3\pi^{3}}}=\frac1\pi\cdot\sqrt{\tfrac{2}{3\pi}}\;=\;S_0\cdot\sqrt{\tfrac{2}{3\pi}}.$$

Interpretation: the frozen escape (this note) captures the **scaling mechanism** — i.e. that the
a-dependence of `S(a)` is algebraic — and gives the clean reduced constant `1/π`; the **x-extension**
of the instanton dresses `1/π → C`, contributing the `√(2/(3π))`.

> **REFINED 2026-07-22 (see `PROGRAM2_TWBETA_XEXTENSION_FACTOR_NOTES.md`):** the frozen→full change is
> NOT a lone multiplicative factor — it is a **joint shift** `(S₀=1/π, γ'=0) → (C, γ=−½)`. The constant
> part `√(2/(3π)) = √(A_g/π)` (`A_g=2/3` the action) is the **x-translation zero-mode Jacobian**, bundled
> with the `−½` index shift (one continuous zero mode). Confirmed to 9 digits; matches Dunne 2511.15528
> eq. 2.16, which writes the HM Stokes constant in exactly this `(1/π)√(2/(3π))` factored form (our `a_n`
> reproduce its `−17/72` subleading to 16 digits). NB the resurgence is **level-variable at β=2**, not
> the noise `1/β` (CLDS treats `1/β` as ordinary cumulants — no Stokes structure).

## Impact on the handoff / paper

- Handoff task **#2 essentially closed**: the a-dependence of `S(a)` is **algebraic** (established);
  the only residual is the elementary `x`-extension factor above, now a bounded finite computation
  rather than an open "is there hidden a-structure" question.
- Suitable as a Remark in `TWbeta_Resurgent_skeleton.tex` (T2c): "the a-dependence of `S(a)` is
  eikonal; under `p=√a P` the tail reduces to an a-independent cubic-barrier escape with Stokes
  constant `1/π`, and `C=(1/π)√(2/(3π))`."

## Files

- `coupled-atlas/_tw_stokes_a_dependence.py` — PDE-scaling check (sympy), exact `q̃_j`, reduced series
  `d_m`, Richardson extraction (`ρ=3/2, γ'=0, S₀=1/π`), MFPT cross-check (`MFPT=1` env to run part 3).

## Gotchas (this session)

- **Neville extrapolation is ill-conditioned at large windows.** Windows `W≳40` at dps 50 blew up
  (γ' → 10¹⁶, poisoning S₀); moderate windows `W≤30` at dps 80 are stable and clean. The growth is
  **single-factorial `Γ(m)`** (base 3/2), NOT the double-factorial `Γ(2k)` of the HM `(−s)^{-3}` series
  — a wrong (Γ(2m)) ansatz is what first made the ratio estimator diverge.
- Exact `(6j)!` fractions to `j=170` are far too slow (`(1440)!` ~ 3900-digit ints, O(M³)); use a
  running-product `mpf` recursion for `q̃_j` instead.
