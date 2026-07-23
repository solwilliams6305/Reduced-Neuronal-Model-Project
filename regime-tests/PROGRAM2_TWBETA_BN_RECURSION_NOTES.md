# PROGRAM2 — TW_β halo: porting the Borot–Nadal recursion for R_{m≥4}

**Date:** 2026-07-22. **Executes:** the residual of `PROGRAM2_TWBETA_UNIFORM_IN_A_NOTES.md` /
HANDOFF_D task #2 — get `R_4` onward (the exact tail noise coefficients at large order) to close the
uniform-in-`a` Gevrey bound. **Memory:** `program2-tw-beta-halo.md`. **Script:**
`coupled-atlas/_tw_bn_brackets.py`. **Source:** Borot–Nadal arXiv:1111.2761 (full PDF extracted).

---

## Headline

The Borot–Nadal `R_m` recursion is **not printed** in the paper — it is a **β-deformed topological
recursion** (their §2.3) that must be reconstructed. I extracted the full machinery from the PDF and
got two clean, immediately usable results, but the faithful port hits **two genuine technical
subtleties** (a branch-point projection and a `W_2` coincident-limit regularization) that make it a real
reconstruction, not a one-line port; and `R_4` specifically needs `ω_1^{[5]}`, which the paper does not
tabulate (it prints the integrals only through `ω_1^{[4]} → R_3`). Status below is honest and staged.

## What is now CLEAN and usable (no recursion needed)

1. **The Coulomb-gas sector is exactly the Bernoulli/Stirling term — a-independent.** Borot–Nadal
   Prop. 2.1 (eq. 2-32) splits the `N^{-m}` exponent coefficient into
   $$-\frac{B_{m+1}}{m(m+1)}\Big(\tfrac2\beta\Big)^{m}\ -\ \beta\!\int_{\alpha(a)}^{\infty}\!\omega_1^{[m+1]}(\alpha')\,d\alpha'.$$
   The first piece has **no `a`-dependence** and (with `B_{m+1}=0` for `m` even) is exactly the
   `Γ(β/2)` Stirling tail `Σ_k B_{2k}/((2k-1)2k)(2/β)^{2k-1}`. **This independently confirms the
   `PROGRAM2_TWBETA_UNIFORM_IN_A_NOTES.md` decomposition**: Coulomb-gas = a-independent (Stirling,
   imaginary-axis), all `a`-dependence (hence all uniformity content) sits in the `∫ω_1^{[m+1]}`
   dynamical piece.

2. **Exact β=1,2,4 tail brackets to `O(s^{-21/2})` (page 5) → `R_m` at three X-points, m=1..6.**
   `log(bracket_β) = Σ_m (β/2) R_m(2/β) s^{-3m/2}`, so the GOE/GUE/GSE brackets give `R_m(2)`,
   `R_m(1)`, `R_m(1/2)`. **Validated** (`_tw_bn_brackets.py`): the paper's `R_1,R_2,R_3` reproduce all
   three brackets exactly. **New data** (m=4,5,6), e.g.
   `R_4(2)=68061/1024`, `R_4(1)=2905/128`, `R_4(1/2)=801555/32768`;
   `R_5(1)=-5218675/49152`; `R_6(1)=7177625/12288`. These are exact finite-β anchors — but **three
   X-points cannot fix `R_m` (deg ≤ m+1) for m≥3**, and in particular give **not** the diagonal
   `r_{m,m+1}` (the `X→∞` top coefficient the uniformity diagonal needs).

## The extracted machinery (for the reconstruction)

Gaussian βE, `V=x²/2t`, edge `a*=2√t`. Uniformize `x=√t(z+1/z)`, `Y=V'-2W_1^{[0]}=t^{-1}√(x²-4t)`,
modified correlators `ω_n^{[k]}(z_i)=Π x'(z_i)·W_n^{[k]}`. Charge position `a=√t(α+1/α)`, `|α|>1`.
- **1-pt recursion (eq. 2-19)**, Gaussian (`Q_1^{[0]}=1/t`, `Q_1^{[1]}=-1/t`, else 0):
  `Y·W_1^{[k]} = Q_1^{[k]} + W_2^{[k-2]}(x,x) + Σ_{k'=1}^{k-1} W_1^{[k']}W_1^{[k-k']} + (1-2/β)∂_x W_1^{[k-1]} + 2(W_1^{[k-1]}(x)-W_1^{[k-1]}(a))/(x-a)`.
- **n-pt recursion (eq. 2-20)** for `W_n^{[k]}`, n≥2 (needed for the `W_2^{[k-2]}(x,x)` term).
- **Given checkpoints:** `ω_1^{[0]}=z^{-1}-z^{-3}`, `ω_1^{[1]}` (2-26), `ω_2^{[0]}` (2-27),
  `ω_1^{[2]}` (2-28), `ω_2^{[1]}` (2-29), `ω_3^{[0]}` (2-30).
- **Assembly (Prop. 2.1, eq. 2-32)** with Selberg prefactor; **explicit** `∫ω_1^{[2]},∫ω_1^{[3]},∫ω_1^{[4]}`
  (eqs. 2-33, 2-34, and the `ω_1^{[4]}` block) — these give `R_1,R_2,R_3`. Edge-scale `a=2√t+N^{-2/3}s`
  (`α→1`, `α²-1 ≈ 2√(a-2)`), leading singular term of `-β∫ω_1^{[m+1]}` `+` the Bernoulli term `→` `(β/2)R_m`.

## Implementation status (the blocker, precisely)

`_tw_uniform_in_a.py` companion probe (scratchpad `bn_probe*.py`): the 1-pt recursion is coded and
- **reproduces `ω_1^{[1]}` up to a clean branch-point term:** naive `RHS/Y` matches eq. 2-26 except the
  residues at `z=±1` are short by `b/2=1/β` each (i.e. missing `+(b/2)(1/(z-1)+1/(z+1))`). This is the
  **β-deformation branch projection** — the `(1-2/β)∂_x` term needs a quantum correction at the
  branch points that naive `Y`-division drops.
- **`k=2` has an additional `W_2` coincident-limit subtlety:** `ω_2^{[0]}=(2/β)/(z_1z_2-1)²` is the
  *regularized* (Bergman-subtracted) 2-point; `W_2^{[0]}(x,x)` in eq. 2-19 is **not** the naive
  `z_1=z_2` value — the coincident limit needs the Eynard-formalism regularization. Until both are
  fixed, `k≥2` is not validated.

Neither subtlety is exotic (both are standard in Chekhov–Eynard–Marchal β-TR), but both need the
precise §2.3 prescription, which the paper states abstractly ("invert `(V'-2W_1^{[0]})Id+O`") and does
not spell out. This is why `R_4` is not delivered here.

## Path to `R_4` onward (concrete)

1. Fix the two subtleties: (a) branch projection `= + Σ_{±1} (residue rule)`; (b) `W_2` coincident
   regularization `W_2^{[0]}(x,x) = lim[ω_2^{[0]}/(x'x') + Bergman-reg]`. **Validate** against the
   given `ω_1^{[1,2]}, ω_2^{[0,1]}, ω_3^{[0]}` (2-26..2-30) and then `R_1,R_2,R_3` via the explicit
   integrals 2-33/2-34 + the edge extraction.
2. Run the recursion to `ω_1^{[5]}` (needs `ω_n^{[k]}`, `n+k≤6`), edge-extract via Prop. 2.1 → `R_4`;
   continue to `R_5,R_6,…`. **Cross-check** each against the `_tw_bn_brackets.py` values at X=2,1,1/2.
3. Feed the diagonal `r_{n,n+1}` (top-degree coeffs) into the uniformity Gevrey bound
   (`PROGRAM2_TWBETA_UNIFORM_IN_A_NOTES.md`, residual (i)) and the off-diagonal into residual (ii).

Alternative to (1)–(2): find a reference/CAS package with the β-deformed TR for the Gaussian/Airy
curve already implemented (Marchal; Chekhov–Eynard) and specialize.

## Files
- `coupled-atlas/_tw_bn_brackets.py` — β=1,2,4 bracket data (m≤6), validates `R_1,R_2,R_3`, tabulates
  `R_4,R_5,R_6` at X=2,1,1/2.
- BN PDF extracted this session (23pp); machinery transcribed above with eq. numbers.

## Gotchas
- **The recursion is not in the paper.** "The first few `R_m` are:" (1-18..1-20) + "a method to find
  them recursively" (§1.3) = the β-TR of §2.3; only `R_1,R_2,R_3` and `∫ω_1^{[≤4]}` are printed.
- **Naive `RHS/Y` is wrong at the branch points** by `O(1/β)` (the `(1-2/β)∂_x` quantum term).
- **Naive coincident `W_2^{[0]}(x,x)` is wrong** — needs Bergman regularization.
- **Three finite-β points don't fix `R_{m≥3}`;** they anchor/validate but miss the `X→∞` diagonal.

---

## UPDATE — the "ready β-TR / Airy curve" route (Kidwai–Osuga arXiv:2204.12431)

Tried the **refined topological recursion** (a clean, geometric β-TR; refinement parameter
`𝒬=√β−1/√β`, Ω-background `ε₁=ℏ√β`, `ε₂=−ℏ/√β`, so `ε₁+ε₂=ℏ𝒬`, `ε₁ε₂=−ℏ²`). Its **quantum curve**
is a 2nd-order ODE `(ε₁²∂ₓ²+q ε₁∂ₓ+r)ψ^{TR}=0` — the exact resummation of all `ω_{g,n}` (so `ln ψ^{TR}`
IS the full tail exponent, no missing pieces, **no projection/coincident subtleties**). Two results:

- **The refined AIRY curve is TRIVIAL — a dead end for TW_β.** For `y²=x` (`x=z², y=z`), eq. (4.27) is
  `(ε₁² d²/dx² − x)ψ^{TR}=0` — the **plain Airy equation** with `ℏ→ε₁=ℏ√β`. No `μ,ν` params, no β-content
  beyond the trivial `ε₁` rescaling. So the local Airy curve carries **none** of the TW_β tail's
  `R_m(2/β)` β-dependence. (Makes sense: TW_β lives in the Airy *kernel*/PII Fredholm determinant, not
  the Airy *curve* wavefunction.)
- **TW_β lives on the WEBER curve — explicit clean ODE.** Weber `y²=¼x²−m∞` is the Gaussian-βE
  semicircle (`x=√m∞(z+1/z)`, `y=½√m∞(z−1/z)`; at `m∞=1`, `Y=½√(x²−4)`). Its refined quantum curve
  (eq. 4.23) is a **β-deformed parabolic-cylinder ODE**:
  $$\Big(\varepsilon_1^2\tfrac{d^2}{dx^2} + C + \tfrac{4m_\infty-x^2}{4}\Big)\psi^{TR}(x)=0,\qquad
    C=\tfrac{\varepsilon_2(\nu_{\infty+}-\nu_{\infty-})-(\varepsilon_1+\varepsilon_2)\mu_{\infty+}}{2}$$
  (`C` is a *constant*, linear in `ε₁,ε₂`). Leading WKB `S₀=∫√(x²/4−m∞)dx ~ (2/3)(x−2√m∞)^{3/2}`
  reproduces the `(2/3)βs^{3/2}` rate — the route is structurally correct.

**Why this is the clean path to `R_{m≥4}` (recommended):** WKB of a *single* ODE is a plain recursion
(`S_k` from `S_{k-1}`), sidestepping the loop-recursion's branch-projection and `W_2` coincident
subtleties entirely. Remaining steps (a focused computation, not a package call): (1) WKB-expand
`W=ψ'/ψ`; (2) fix the identification to the Borot–Nadal resolvent (map `ε₁,K,m∞,μ,ν` ↔ `N,β`), pinned by
matching the **exact** `ω_1^{[1]}` (eq. 2-26) / `∫ω_1^{[2]}` (eq. 2-33); (3) edge-extract via Prop. 2.1;
(4) generate `R_4,R_5,…` and **cross-check against `_tw_bn_brackets.py`** (X=2,1,1/2). The Weber ODE
(4.23) is the concrete object to WKB. **No off-the-shelf code emits `R_m`** — the refined-TR packages
give the formalism/curve, not the TW_β tail coefficients.

---

## UPDATE 2 — the constant-K Weber ODE FAILS; the Chekhov–Eynard–Marchal dictionary (the fix)

**Ran the Weber WKB** (`coupled-atlas/_tw_weber_wkb.py`). **Structural obstruction:** Kidwai–Osuga
eq. 4.23 is a **pure Schrödinger** ODE (no `d/dx` term — Weber has no finite poles), so β can enter
only through the constant `K`. Full edge double-scaling shows the `K`-dependence **cancels**: with
`K=1+q₁X` one gets `R_1 = −15/8` *identically*, independent of `q₁`. And the pure WKB `w_k` (leading
edge coeffs, `−1,−1/4,5/32,−15/64,1105/2048,…`) reproduce the exact top coeff `r_{m,m+1}` only at
non-deficit orders (R_1,R_3) and **alternate** in sign, whereas the true diagonal `2^n r_{n,n+1}` is
non-alternating. So a constant-K, β-in-the-potential ODE cannot make TW_β's β-dependence.

**The Chekhov–Eynard–Marchal dictionary (arXiv:1009.6007 §8) — where β actually enters.**
Measure `Z=∫dλ |Δ(λ)|^{2β} ∏ e^{−(N√β/t₀)V(λ_i)}` (their **β = Dyson/2**; β=1 = Hermitian). Then:
- **Quantum curve** `((ℏ∂)² − U(x))ψ(x)=0` (a Schrödinger ODE), with
- **`ℏ = (t₀/N)(√β − 1/√β)`** (eq. 8.3) — β enters through the **KINETIC** term ℏ, not the potential.
  `ℏ=0` ↔ β=1 (Hermitian); `ℏ→−ℏ` ↔ `β→1/β`.
- Resolvents `W_k = β^{k/2}⟨Σ…⟩_c` (eq. 8.4); loop eq. `W₁²+(√β−1/√β)W₁'+W₂=(N/t₀)(V'W₁−P₁)` (eq. 8.11).
- Gaussian `U(x)=x²/4 − t₀ + O(ℏ)` (semicircle + quantum correction; edge at `2√t₀`).

**So the fix is:** the WKB parameter is `ℏ ∝ (√β−1/√β)/N` (β-dependent), NOT `1/N`. My Weber attempt
used `ε₁=1/N` and put β in `K` — wrong slot. In CEM the genus-0 resolvent is *already* quantum
(Riccati `ω²+ℏω'=V'ω−P₁`, eq. 8.12), and β rides in `ℏ`.

**Honest caveat (why this is not yet a plug-in).** Converting the CEM curve to Borot–Nadal's
**polynomial** `R_m(X)` (X=2/β_Dyson=1/β_CEM) needs the full normalization bookkeeping: (i) `ℏ^m ∝ 𝒬^m`
with `𝒬=√β−1/√β=(1−X)/√X` is **non-polynomial** in X on its own; (ii) the `β^{k/2}` resolvent
normalization and the `N√β/t₀` measure factor must combine with `𝒬^m` to restore a polynomial (and to
fix the rate — the naive `1/𝒬` leading factor does NOT match `−(4/3X)s^{3/2}`). This bookkeeping is the
remaining derivation. **Dictionary = found; assembly of `R_4` = a careful (not mechanical) next step.**
