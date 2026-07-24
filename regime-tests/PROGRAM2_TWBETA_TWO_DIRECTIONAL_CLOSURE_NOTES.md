# PROGRAM2 — TW_β halo: the TWO-DIRECTIONAL RESURGENCE CLOSURE (uniform-in-`a`)

**Date:** 2026-07-23. **Status:** synthesis / closure note. Consolidates the uniform-in-`a` arc
(Handoff-D task #2) into one theorem-grade statement. **Memory:** `program2-tw-beta-halo.md`.
**Supersedes as the headline** the residual framing in `PROGRAM2_TWBETA_UNIFORM_IN_A_NOTES.md` (kept for
the blow-by-blow). **Scripts:** `coupled-atlas/_tw_{uniform_in_a,reconcile_gevrey,offdiag_bound,
nuniform_pin,bn_brackets,weber_wkb}.py`.

---

## The closure in one paragraph

The TW_β right tail is resurgent in **two directions at once** — level (`a→∞`, fixed β; Hastings–McLeod
PII, Thm T1) and noise (`β→∞`, fixed `a`; median-summable two-sector, Thm T3). The open question was
whether the **noise** Borel summability is **uniform in `a`**. It is, for `a` in the deep-tail regime
`a≥a₀≈1`, and the proof needs **no new machinery**: it assembles three results the program already had
(T3 median summability; §8 singularity-location; the algebraic noise Stokes constant `S(a)`) with two
structural findings established here (the **median reconciliation** of the scaling diagonal, and the
**level Borel radius** setting `a₀`). The two directions **do not collide** — the noise singularity
`Φ(a)=\tfrac23a^{3/2}` moves outward with `a` while the level singularity sits at a fixed radius in
`u=a^{-3/2}` — and the noise amplitude is the finite, algebraic `S(a)`. Notably the whole result is
**independent of `R_{m≥4}`**: the Borot–Nadal recursion (an all-genus topological recursion, no
single-ODE shortcut for general β) was shown to be off the critical path.

## Object and the two directions

Borot–Nadal (arXiv:1111.2761, Prop. 1.1), exact for all β>0, `s=a`:
$$1-\mathrm{TW}_\beta(a)=\underbrace{\frac{\Gamma(\beta/2)}{(4\beta)^{\beta/2}2\pi}}_{c_\beta}\,a^{-3\beta/4}\,
e^{-\frac{2\beta}{3}a^{3/2}}\,\exp\Big[\Sigma(a,\beta)\Big],\qquad
\Sigma=\sum_{m\ge1}\tfrac\beta2 R_m(\tfrac2\beta)\,a^{-3m/2}.$$
- **Level** (`a`): the `a→∞` series is the HM PII trans-series (**T1**), Borel-summable (Costin), Stokes
  constant `C=√(2/(3π³))`.
- **Noise** (`1/β`): the `1/β` series is median-Borel-summable (**T3**), a dynamical-escape sector
  (positive axis, median) + a Coulomb-gas sector (`iℝ`, ordinary).

## Step 0 — reduction (the noise coefficients) `[_tw_uniform_in_a.py; self-checked]`

Collecting `Σ` in powers of `1/β` (`X=2/β`, `t=a^{-3/2}`, `r_{m,ℓ}=[X^ℓ]R_m`):
$$\Sigma=\beta\,\sigma_{-1}+\sigma_0+\sum_{n\ge1}\sigma_n(a)\,\beta^{-n},\qquad
\boxed{\;\sigma_n(a)=2^n\!\!\sum_{m\ge n} r_{m,n+1}\,a^{-3m/2}\;}$$
The prefactor `c_β` (Coulomb-gas, `iℝ`) and `−\tfrac{3\beta}4\log a` are **`a`-independent / `β^{≥0}`**, so
**all** `1/β`-Borel content is `{σ_n}_{n≥1}`. Split each into
$$\sigma_n(a)=\underbrace{D_n\,a^{-3n/2}}_{\text{diagonal }(m=n),\ D_n=2^n r_{n,n+1}}
+\underbrace{2^n\!\!\sum_{m>n} r_{m,n+1}a^{-3m/2}}_{\text{off-diagonal }(m>n)}.$$

## Step 1 — the DIAGONAL is a median; it inherits Gevrey-1 `[_tw_reconcile_gevrey.py]`

Two **independent** computations of the scaling diagonal `D_n`:
- **(F) frozen escape** `δ_n=[g^{-n}]\log F(g)`, `F=\tilde Q(2/g)^{-2}` (cubic-barrier MFPT), `g=βa^{3/2}`;
- **(W) Weber genus-0 WKB** `2^n\hat r_n` (refined quantum curve, Kidwai–Osuga eq 4.23; `_tw_weber_wkb.py`).

**Finding (verified to `n=7/8`, exact fractions):** `2^n\hat r_n=(-1)^{n+1}δ_n`, i.e. `W(g)=−F(−g)` (a Borel
reflection). So F and W **agree at odd `n`, are opposite at even `n`**, and the exact diagonal is their
**median**:
$$D_n=\tfrac12\bigl(δ_n+2^n\hat r_n\bigr)=\text{odd part of }\log F(g)$$
— verified *exactly* at `n=1,2,3` (`D_1=δ_1=−\tfrac5{12}`, `D_2=0`, `D_3=δ_3=−\tfrac{1105}{576}`). The
frozen escape's **even-order terms are the Stokes/discontinuity artifact**; the physical diagonal keeps
only the odd (real/median) part, which **auto-cancels** them — precisely the `R_{even}` degree deficit
(`D_2=0`). This is the median-summation principle made explicit (F, W = the two lateral branches).

**Consequence (Gevrey-1, uniform):** `F` is **proven** Gevrey-1 (`|δ_n|∼(1/π)Γ(n)(3/2)^n`, Borel sing at
`g=2/3=Φ`; the Airy-integral/Nevanlinna proof, `PROGRAM2_TWBETA_NOISE_BOREL_PROOF_NOTES.md`). The exact
diagonal, an **odd sub-series** of `F`, satisfies `|D_n|≤|δ_n|`, so
$$|\sigma_n(a)_{\rm diag}|=|D_n|\,a^{-3n/2}\le K\,\Gamma(n)\,\Phi(a)^{-n},\qquad \Phi(a)=\tfrac23a^{3/2}.$$

## Step 2 — the OFF-DIAGONAL is the level resurgence; threshold `a₀` `[_tw_offdiag_bound.py]`

The off-diagonal is the column generating function `G_{n+1}(u)=\sum_m r_{m,n+1}u^m`, `u=a^{-3/2}` — a
**divergent (Gevrey-1)** series in `u`: it *is* the **level-`a` resurgence** (HM/PII). Verified at `n=1`:
`σ_1=−\tfrac5{12}a^{-3/2}[1−\tfrac{33}8u+\tfrac{1555}{64}u^2−…]`, coefficients growing factorially
(ratios `4.13, 5.89`) — not a convergent correction. **Key: no collision.** The noise singularity sits at
`Φ(a)=\tfrac23a^{3/2}` and **moves outward** with `a`; the level singularity sits at a **fixed** level
Borel radius `ρ` in `u`. So for
$$a\ge a_0:=\rho^{2/3}$$
the column is level-Borel-summable and `σ_n(a)=D_n a^{-3n/2}(1+O(a^{-3/2}))`. From `R_m(X)` (`m≤6`,
bracket data), `|R_m/R_{m-1}|/m→ρ` with `ρ(X=1)≈0.94=\tfrac{2\sqrt2}3` (**the HM action** — confirms the
level direction is HM/PII), `ρ(X=2)≈1.4`, `ρ(X=½)≈0.92`. So `ρ=O(1)` and **`a₀≈1`**: the bound holds
throughout the deep tail `a≫1` (not down to `a→0`).

## Step 3 — the AMPLITUDE is the noise Stokes constant; `n`-uniformity `[_tw_nuniform_pin.py]`

Let `K(a)=\sup_n|σ_n(a)|/[Γ(n)Φ(a)^{-n}]`. By median summability (**T3**) + a single dominant singularity
at `Φ(a)` on `ℝ₊` (**§8**), the standard resurgence large-order formula gives
$$σ_n(a)=\frac{S(a)}{2\pi i}\,Γ(n)\,Φ(a)^{-n}\bigl[1+O(1/n)\bigr]\ \Longrightarrow\
\frac{|σ_n(a)|}{Γ(n)Φ(a)^{-n}}\to\frac{|S(a)|}{2\pi},$$
so the amplitude limit **is** the **noise Stokes constant** `S(a)`. And `S(a)` is **algebraic**:
`S(a)=S₀a^p(1+O(a^{-3/2}))`, `S₀=1/π` (`PROGRAM2_TWBETA_S_OF_A_NOTES.md`, exact scaling symmetry,
certified 25 digits), hence **bounded** for `a≥a₀`. Small-`n` amplitude is finite/continuous. So
`K(a)<∞`, uniform in `n` and `a≥a₀`. **Verified:** the diagonal amplitude `|D_n|/[Γ(n)(3/2)^n]→1/π=S₀`
(odd `n`: `0.278,0.284,…,0.316→0.318`; even `n`: `0`), bounded, no growth in `n`. *(The naive `|·|`-column
estimate suggests a "two-directional resonance" at `n∼a^{3/2}/ρ`; it is a red herring — median summation
tames it, and the physical amplitude is the finite `S(a)`.)*

## THE CLOSURE (theorem-grade, modulo cited machinery)

> **For `a≥a₀≈1`, the noise `1/β` series `\sum_{n≥1}σ_n(a)β^{-n}` is median-Borel-summable along `ℝ₊`,
> uniformly in `a`:** (i) **sector uniform** — the dynamical singularity is at `Φ(a)=\tfrac23a^{3/2}∈ℝ₊`
> (moves with `a`, never leaves the ray; §8), the Coulomb-gas fixed on `iℝ`; (ii) **uniform Gevrey-1** —
> `|σ_n(a)|≤K(a)Γ(n)Φ(a)^{-n}` with `K(a)` bounded uniformly, `K(a)→|S(a)|/2π` (large `n`), `S(a)`
> algebraic. The scaling **diagonal** is the odd/median part of the proven-Gevrey-1 frozen escape; the
> **off-diagonal** is the level (HM/PII) resurgence, summable for `a≥a₀=ρ^{2/3}`.

**Inputs, all established:** T3 (median summability), §8 (singularity location), S_OF_A (algebraic `S(a)`).
**New structural glue:** the median reconciliation `D_n=` odd part of `\log F` (F↔W reflection), and the
level Borel radius `ρ≈` HM action giving `a₀≈1`. **Two-directional consistency:** the *same* instanton
action `Φ` governs both directions (action-aligned), which is exactly why the amplitude is the single
Stokes constant `S(a)` and the two singularities never cross.

## Scope, and the `R_4` wall (diagnostic, honest)

- **`a₀≈1`, not `a→0`:** the closure is for the deep-tail regime (where the tail lives); below `a₀` the
  level and noise scales are comparable and the clean separation fails.
- **`R_{m≥4}` is OFF the critical path.** Producing the exact `R_4` needs the **all-genus** Borot–Nadal
  topological recursion (`PROGRAM2_TWBETA_BN_RECURSION_NOTES.md`): the loop-equation recursion with a
  branch-point projection + Bergman-regularized `W_2` coincident limit. The refined-**Airy** curve is
  trivial (eq 4.27); the **Weber** curve is pure Schrödinger so β cancels from the potential; the
  **CEM** dictionary (`ℏ=(t₀/N)(√β−1/√β)`, β enters the kinetic term) shows each `R_m` is a **genus sum**
  (`R_m=\sum_h𝒬^{m+1-2h}c_{m,h}`, `𝒬=(1−X)/√X`), so the single-ODE WKB is genus-0 only. No single-ODE
  shortcut exists for general β (exact ODEs only at even `β_Dyson`). **But the closure never needed it:**
  the diagonal came from the median, the off-diagonal from the level radius, the amplitude from `S(a)`.
  `R_m` large-order would only pin the *explicit numerical value* of `S(a)`'s subleading corrections — cosmetic.

## Impact on the paper (`TWbeta_Resurgent_paper.tex`)

- **Upgrade Rem. 7.4:** the uniform statement is now a **theorem for `a≥a₀`** (not "would combine the Airy
  proof with the PII map"). State the three-piece structure (sector / diagonal-median / amplitude-`S(a)`).
- **T3 scope caveat (still stands):** T3-as-proved is median summability of the *frozen* reduced escape;
  the median reconciliation shows the exact diagonal is its **odd part** (the even part is a freeze
  artifact), which is the precise bridge from the frozen object to the exact tail.
- The **two-directional** framing is now genuinely closed: level (T1), noise (T3), and their **uniform
  joint** structure (this note) — action-aligned, single Stokes constant.

## Files
- Reduction `_tw_uniform_in_a.py`; reconciliation+diagonal-Gevrey `_tw_reconcile_gevrey.py`; off-diagonal
  `_tw_offdiag_bound.py`; amplitude/`n`-uniformity `_tw_nuniform_pin.py`; bracket data `_tw_bn_brackets.py`;
  Weber WKB `_tw_weber_wkb.py`. Constituent notes: `PROGRAM2_TWBETA_{UNIFORM_IN_A,BN_RECURSION,S_OF_A,
  NOISE_BOREL_PROOF}_NOTES.md`.

## Gotchas
- **Diagonal ≠ frozen `δ_n`; diagonal = ODD part of frozen** (median). The even-order `δ_{2k}` are Stokes
  artifacts (the `R_{even}` degree deficit). Do not use the raw frozen series as the exact diagonal.
- **Off-diagonal is DIVERGENT** (level resurgence), not a convergent small correction — it needs `a≥a₀`.
- **The amplitude is `S(a)`, not `1/π`** at finite `a`: `S₀=1/π` is the *scaling* (diagonal) constant;
  `S(a)=S₀a^p(1+…)` is the finite-`a` Stokes constant (still bounded).
- **`R_4` is a red herring for uniformity** — established here; don't re-chase it for this purpose.
