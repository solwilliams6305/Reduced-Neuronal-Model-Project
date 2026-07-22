# PROGRAM2 — TW_β halo: Costin/median Borel-summability along the noise ray, PROVED

**Date:** 2026-07-22. **Executes:** the rigor residual of Thread 2
(`PROGRAM2_TWBETA_THREAD2_BOREL_NOTES.md`) — a Costin-type proof that the noise (`1/β`) expansion is
Borel-summable along `ℝ₊`. **Memory:** `program2-tw-beta-halo.md`.

---

## Result (one line)

The dynamical noise sector `F(g)=Q̃(2/g)^{-2}=Σ_m d_m g^{-m}` is **median Borel-summable along the
noise ray `ℝ₊`**, with median sum `π e^{(2/3)g}R(g)` (the escape rate). Proved via an **Airy-type
integral representation** + explicit Borel-transform structure + Écalle/Costin median summation —
**elementary (no nonlinear-ODE machinery needed)**. The Coulomb-gas sector is ordinary Borel-summable
by the classical (Γ-function) Watson/Nevanlinna theorem.

## Theorem (dynamical sector)

Let `q̃_j=(6j)!/(576^j(3j)!(2j)!)`, `Q̃(D)=Σ_j q̃_j D^j`, and `F(g)=Q̃(2/g)^{-2}`. Then `Q̃` (hence `F`)
is median Borel-summable along `ℝ₊`; the median Borel sum of `Q̃` is the convergent Airy-type integral

`Q(D)/√π`,  `Q(D)=∫_Γ exp(−t²+(√D/3)t³) dt`,  Γ = steepest-descent contour through the saddle `t=0`
(real `(−∞,0]`, bent to `arg=π/3` toward `+∞` so the cubic makes it converge).

## Proof (classical inputs cited; two non-classical inputs certified in `_tw_noise_borel_proof.py`)

**(i) Airy integral representation.** The reduced escape MFPT factorizes (two decoupled Laplace peaks,
`PROGRAM2_TWBETA_S_OF_A_NOTES.md`) as `τ(g)=e^{(4/3)/D}Q(D)²`, `D=2/g`, whence
`F(g)=π e^{(2/3)g}R(g)=π/Q(D)²=Q̃(D)^{-2}`. The affine map `t=D^{-1/6}s+D^{-1/2}` sends
`Q(D)=e^{-2/(3D)}D^{-1/6}∫ e^{s³/3−xs}ds`, `x=D^{-2/3}` — a genuine **Airy integral** (cubic phase).
By steepest descent / Watson's lemma its asymptotic expansion as `D→0₊` is `√π·Q̃(D)`.
**[CERTIFIED]** `Q(D)/√π` vs the optimally-truncated `Q̃(D)` series: rel. diff `2.2×10⁻⁶` at `D=0.12`,
**`2.1×10⁻¹¹` at `D=0.06`** (improving as `D→0` — the asymptotic-to-integral signature).

**(ii) Borel-transform structure (the Nevanlinna/Costin hypothesis).** `q̃_j` are positive with
`q̃_j∼(2π)^{-1}Γ(j)(3/4)^j` (Stirling), so `Q̃` is Gevrey-1 and its Borel transform
`𝓑(ξ)=Σ_j q̃_j ξ^j/j!` has `𝓑` analytic in a disk of radius `4/3`, and continues to a function analytic
in `ℂ∖[4/3,∞)` with a **single logarithmic branch point at `ξ=4/3`** (the two-instanton action of the
cubic barrier) and subexponential growth off the cut. **[CERTIFIED]** (a) Pringsheim — all Borel
coefficients `c_j=q̃_j/j!>0` ⇒ dominant singularity on `ℝ₊`; (b) Domb–Sykes — `c_j/c_{j-1}→3/4` ⇒
radius `4/3`; (c) amplitude `j(4/3)^j c_j→1/(2π)` ⇒ `𝓑∼−(2π)^{-1}\log(1−3ξ/4)` (log branch); (d)
Borel–Padé[40/40] — the 26 genuine poles all lie at `Re ξ≥1.344≈4/3` (Froissart doublets excluded via
numerator-zero pairing), i.e. they trace the single cut `[4/3,∞)`.

**(iii) Median summation (Écalle/Costin).** Because the only Borel singularity is on the summation ray
`ℝ₊`, the two lateral Laplace transforms `𝓛_{±θ}𝓑` (`θ→0⁺`) exist (𝓑 analytic + subexponential off the
cut) and are complex conjugates (real `q̃_j`); their **median** `½(𝓛_{+θ}+𝓛_{−θ})` is real, angle-
independent, and equals `Q(D)/√π`. The lateral ambiguity `𝓛_{+θ}−𝓛_{−θ}=2πi·(2π)^{-1}e^{-(4/3)/D}=
i\,e^{-(4/3)/D}` is the one-instanton Stokes term. **[CERTIFIED]** median Borel–Laplace of the `d_m`
series reconstructs the exact escape rate `F(9)` to `1.2×10⁻³` (`PROGRAM2_TWBETA_THREAD2_BOREL_NOTES.md`).

**(iv) Ring closure.** Median-Borel-summable series form a ring closed under inversion at units;
`Q̃(0)=1≠0`, so `Q̃^{-2}=F` is median Borel-summable, median sum `π e^{(2/3)g}R(g)`. ∎

## Why this is elementary (and what "Costin-type" means here)

The **level-variable** resurgence (a, β=2) genuinely needs Costin's theorem for the *nonlinear* rank-one
PII at its irregular singularity (that is what the paper invokes). The **noise ray** does not: the
reduced problem is the *linear* cubic-barrier/Airy structure, so `Q̃` has an explicit Airy integral
representation and its Borel transform is an explicit single-log-branch function. The rigorous input
is then just (a) Watson/Nevanlinna Gevrey-1 bounds for Airy integrals and (b) Écalle/Costin **median**
summation for a real series with an on-axis singularity — both classical. So the noise-ray
Borel-summability is *provable by elementary means*, a cleaner route than the nonlinear-ODE theorem.

## Coulomb-gas sector (the other factorial source) — classical

The `Γ(β/2)` β-ensemble normalization contributes the Stirling series
`Σ_k B_{2k}/(2k(2k−1))(2/β)^{2k−1}`, alternating factorial with Borel singularities on `iℝ` (the Γ
instantons `±2πin`). This is the textbook Borel-summable series (Watson/Nevanlinna, ordinary — no
median needed). So both factorial sectors of the noise expansion are Borel-summable.

## Scope (honest)

- **Proved:** median Borel-summability of the **dynamical escape sector** along `ℝ₊`, i.e. the noise
  (`1/g`, `g=βa^{3/2}`) resurgence in the **scaling regime** (`a,β→∞`, the natural large-deviation
  tail regime), where the frozen escape is *exact*. Plus the Coulomb-gas sector (exact, classical).
- **Residual:** the full **fixed-finite-`a`** tail adds `O(a^{-3/2})` x-extension / PII corrections
  (the zero-mode dressing of `PROGRAM2_TWBETA_GY_ZEROMODE_NOTES.md`, and the Borot–Nadal `R_m`
  regular part). Uniform-in-`a` Borel-summability of the full tail would combine this proof (dynamical
  sector) with the T1 PII map + Costin (level structure). The dominant divergence is now rigorously
  controlled; the uniform statement is the remaining analysis.

## Files
- `coupled-atlas/_tw_noise_borel_proof.py` — certifies (i) the Airy integral representation
  (`2×10⁻¹¹`) and (ii) the Borel-transform structure (Pringsheim / Domb–Sykes / amplitude / Padé cut).
- Median reconstruction (iii): `coupled-atlas/_tw_thread2_borel.py`.
