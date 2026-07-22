# HANDOFF D — TW_β halo: continuing after the two-directional resurgence + full paper

> **STATUS 2026-07-22.** The TW_β "halo result" is **substantially complete on resurgence structure**
> and consolidated into a full paper: **`TWbeta_Resurgent_paper.tex` (repo root, 8pp, builds clean,
> proofread).** This handoff **supersedes `PROGRAM2_HANDOFF_C_TW_BETA_HALO.md`** (which accreted six
> `cont.` banners — read it only for the blow-by-blow history). Remaining work = (a) paper→submission
> polish, (b) four scoped math refinements, (c) the multi-year frontier. **Memory:**
> `program2-tw-beta-halo.md`. **Under git** (2 commits; local only, no remote).

_Self-contained handoff for a fresh Claude Code session. Companion to the cusp paper
`W_ResurgentTransseries_paper.tex` (Program 2's `q=2`, numerical). Reads: this file, then the paper
`TWbeta_Resurgent_paper.tex`, then the `PROGRAM2_TWBETA_*.md` notes as needed._

## The result in one paragraph (the load-bearing picture — don't lose it)

The TW_β right tail (`𝒯_β(a)=P(Λ₀<−a)`, stochastic-Airy bottom eigenvalue Λ₀, weak noise ε=1/√β) is
the first random-operator edge law placed in a Borel-summable resurgent trans-series, and it is
resurgent in **TWO independent directions**:

- **(L) LEVEL** (`a→∞`, fixed β): via the affine **Lax map** `s=−2^{1/3}(x+a)`, `q=−2^{−1/3}p⋆`, the
  fluctuation (Jacobi) operator `𝒥=−∂_x²+6p⋆²−2(x+a)` **is** the exact PII linearization
  `μ²(−∂_s²+6q²+s)`, and the instanton is the **α=1/2 Airy special solution of PII**. So the `a→∞`
  expansion is the Hastings–McLeod trans-series; Borel-summable by **Costin**; Stokes constant
  **`C=√(2/(3π³))=(1/π)√(2/(3π))=0.146632271193848478…`** (28 digits, = Dunne arXiv:2511.15528
  eq.2.16 incl. its `−17/72` subleading to 16 digits; **refutes the guess `1/(4√π)`**).
- **(N) NOISE** (`β→∞`, fixed a): the `1/β` loop expansion is **median Borel-summable**, a **two-sector**
  trans-series — a **dynamical escape** sector (positive-real-axis, median; proved *elementarily* via
  an **Airy-type integral representation** of the reduced escape — no nonlinear-ODE theorem) and a
  **Coulomb-gas `Γ(β/2)`** sector (imaginary-axis, ordinary Borel).

The **a-dependence of the Stokes constant is ALGEBRAIC**: `C = (1/π) × √(2/(3π))`, where `1/π` is the
transverse Kramers determinant and `√(2/(3π))=√(‖ψ₀‖²/2π)` is the **Gelfand–Yaglom zero-mode Jacobian**
of the tanh-kink escape instanton `P⋆=−tanh τ` (`‖ψ₀‖²=∫sech⁴=4/3=ΔV=2A_g`, `A_g=2/3`).

## What is DONE — do NOT rebuild (all certified; in the paper + notes)

1. **Exact Stokes constant `C=√(2/(3π³))`** — `PROGRAM2_TWBETA_STOKES_CONSTANT_NOTES.md`,
   `coupled-atlas/_hm_stokes_exact.py` (+ `_hm_stokes_validate.py`, `_hm_coeffs_exact.json`). Exact
   rational recursion for the HM `s→−∞` coeffs `a_k` (`a_1=−1/8, a_2=−73/128, …`) to N=120 +
   Richardson: `γ=−1/2` (25d), `C` (28d); matched to Dunne + Borel/ODE cross-checks. Refutes `1/(4√π)`.
2. **a-dependence is algebraic** — `PROGRAM2_TWBETA_S_OF_A_NOTES.md`, `_tw_stokes_a_dependence.py`.
   Scaling symmetry `p=√a·P` (sympy-verified) ⇒ `a` enters only through `g=βa^{3/2}`; reduced
   cubic-barrier escape `R(g)=(1/π)e^{−(2/3)g}/Q̃(2/g)²`, `Q̃_j=(6j)!/(576^j(3j)!(2j)!)`, reduced Stokes
   constant `S₀=1/π` (`d_m∼−(1/π)Γ(m)(3/2)^m`).
3. **`√(2/(3π))` = tanh-kink GY zero-mode Jacobian** — `PROGRAM2_TWBETA_GY_ZEROMODE_NOTES.md`,
   `_tw_zeromode_gy.py`. `‖ψ₀‖²=∫sech⁴=4/3=ΔV=2A_g` (virial); `√(‖ψ₀‖²/2π)=√(2/(3π))`; assembles to
   `C=(1/π)√(2/(3π))=√(2/(3π³))`. Also `_tw_frozen_vs_hm.py`: the frozen→full shift is joint
   `(S₀=1/π,γ'=0)→(C,γ=−½)`, confirmed 9 digits.
4. **Thread 2 — the noise `1/β` expansion is Borel-summable** (two sectors) —
   `PROGRAM2_TWBETA_THREAD2_BOREL_NOTES.md`, `_tw_thread2_borel.py`. Dynamical: positive-axis Borel
   sing at `t_c=2/3=A_g` (reduced action Φ; its median ambiguity is the two-instanton `e^{−2βΦ}`);
   median Borel–Laplace reconstructs `R(g)` to `1.2e−3`. Coulomb-gas: `Γ(β/2)` Stirling, alternating,
   imaginary-axis `±iπ`.
5. **Costin/median proof of the noise ray** (dynamical sector, scaling regime) —
   `PROGRAM2_TWBETA_NOISE_BOREL_PROOF_NOTES.md`, `_tw_noise_borel_proof.py`. Airy integral rep
   `Q(D)=∫e^{−t²+(√D/3)t³}dt` (certified vs series to `2e−11`); Borel transform analytic in
   `ℂ∖[4/3,∞)`, single log branch at `4/3` (Pringsheim + Domb–Sykes `→3/4` + amplitude `→1/2π` + Padé);
   Nevanlinna/Écalle–Costin median + ring closure ⇒ Borel-summable.
6. **T2a location-fixing, T2b `Ω₁=−1.12±0.01`** (self-adjoint SAO ⇒ real; two evaluations) — in the
   paper §8; original `O5*.tex` notes.

## OPEN TASKS (priority order)

1. **[highest value] Paper → submission.** `TWbeta_Resurgent_paper.tex` is a complete, proofread draft
   (numerics independently re-verified; abstract "first" claim qualified "to our knowledge"). Remaining:
   (i) a final human read of the abstract + §9 positioning for tone; (ii) pick a venue (CMP / J.Stat.Phys
   / Nonlinearity / a resurgence-friendly venue); (iii) the 4 refinements below are tagged `[to prove]`/
   `[frontier]` in-text and can ship as remarks. Gotcha: the classical citations (TW1994, FN1980, JMU1981,
   FIKN, Clarkson) were web-verified this program (2 title patches applied) — don't re-paraphrase.
2. **[tractable rigor] Uniform-in-`a` noise Borel-summability — REDUCED 2026-07-22**
   (`PROGRAM2_TWBETA_UNIFORM_IN_A_NOTES.md`, `coupled-atlas/_tw_uniform_in_a.py`). Worked from the
   **exact** Borot–Nadal tail (1111.2761 Prop 1.1; `R_m(2/β)` are IN the paper for `m≤3`, validated at
   β=2 vs GUE: `R_1(1)=−35/24`, `R_2(1)=35/16`). **Done:** (i) reduction — the Coulomb-gas prefactor and
   `−(3β/4)log a` are a-independent/`β^{≥0}`, so ALL `1/β`-Borel content is the exact noise coeffs
   `σ_n(a)=2^nΣ_{m≥n}r_{m,n+1}a^{-3m/2}` (`r_{m,ℓ}=[X^ℓ]R_m`, `X=2/β`; self-checked); (ii) **sector
   uniformity is structural** — the dynamical singularity is at the *real* WKB action `Φ(a)∈ℝ₊` (moves
   with `a`, never leaves the ray), Coulomb-gas fixed on `iℝ` ⇒ median sector uniform for `a≥a₀`.
   **FINDING:** the frozen boundary-escape series is **NOT** the exact noise series beyond `n=1` — the
   exact diagonal `2^n r_{n,n+1}` matches frozen `δ_n` at `n=1,3` but **not `n=2`** (`R_2` has degree
   deficit `2<m+1=3` ⇒ exact 2-loop scaling coeff **=0**, frozen `−5/8`). So T3's Airy proof covers the
   *frozen* object; the exact large-order must be read off BN. **RESIDUAL (one estimate):** a uniform
   Gevrey-1 bound `|σ_n(a)|≤KΓ(n)Φ(a)^{-n}` on the exact `σ_n` — needs the **BN loop-equation recursion**
   (§1.3 of 1111.2761) for `R_{m≥4}` (large-order `r_{n,n+1}` growth + off-diagonal `O(a^{-3/2})`
   uniformity; the `n=1` off-diagonal is `σ_1=(−5/12)a^{-3/2}[1−(33/8)a^{-3/2}+(1555/64)a^{-3}+…]`).
   **Next actor:** derive/port the BN recursion, generate `R_4..R_~20`, extract `r_{n,n+1}` growth and
   confirm the on-axis singularity numerically. This is now the natural completion of Thread 2.
   **BN-RECURSION PORT — STARTED 2026-07-22** (`PROGRAM2_TWBETA_BN_RECURSION_NOTES.md`,
   `coupled-atlas/_tw_bn_brackets.py`; full 1111.2761 PDF extracted). The recursion is **not printed** —
   it is a **β-deformed topological recursion** (§2.3): loop eqs 2-8/2-9, Gaussian 1-pt recursion 2-19 in
   uniformization `x=√t(z+1/z)`, assembly Prop 2.1 eq 2-32. **Clean wins:** (a) Coulomb-gas sector = the
   Bernoulli term `−B_{m+1}/(m(m+1))(2/β)^m` (Prop 2.1) — **a-independent**, = the `Γ(β/2)` Stirling
   tail, independently confirming the uniform-note decomposition; (b) page-5 β=1,2,4 brackets to
   `O(s^{-21/2})` give exact `R_m` at `X=2,1,1/2` for m≤6 (validated `R_1..R_3`; new `R_4,R_5,R_6` at 3
   points — e.g. `R_4(1)=2905/128`). **Blocker:** faithful port needs two standard β-TR subtleties the
   paper states only abstractly — a **branch-point projection** (naive `RHS/Y` is short by `+b/2` residue
   at `z=±1`; verified at k=1 vs eq 2-26) and a **`W_2` coincident-limit (Bergman) regularization** (k=2).
   `R_4` needs `ω_1^{[5]}` (paper prints `∫ω_1^{[≤4]}→R_3` only). 3 finite-β points can't fix `R_{m≥4}`
   (deg ≤ m+1) nor the `X→∞` diagonal. **To finish:** fix the two subtleties, validate vs 2-26..2-30 +
   `R_1..R_3`, run to `ω_1^{[5+]}`, edge-extract via Prop 2.1, cross-check vs the bracket data.
3. **[structural] `α=1/2` vs `α=0` family assignment** (paper Rem.~in §5). The tail *saddle* is the
   α=1/2 Airy solution; the β=2 *distribution* is α=0 (HM). Pin the correspondence via the PII
   coalescence / Airy-solution structure. Does NOT block T1 (Costin's class covers all α); sharpens
   *which* PII. Self-contained. See `O1_step1b_laxpair.tex` Q1/Q2.
4. **[novel] Coulomb-gas sector → β-ensemble instanton physics.** The noise expansion's *second* sector
   (the `Γ(β/2)` imaginary-axis instantons `e^{−2πinβ/2}`) is only characterized, not connected to the
   underlying Coulomb-gas / β-ensemble instanton structure. That linkage looks genuinely new and is the
   least-explored piece — potentially its own short paper.
5. **[technical] Rough-potential rigor** of the isomonodromy rigidity (T2a): white noise really
   preserving the rank-3/2 irregular singularity (the Voros periods over the Airy spectral curve being
   noise-independent). Controlled by the O5 renormalization; make rigorous.
6. **[validation] Numerical end-to-end check.** Direct Monte-Carlo / ODE of `𝒯_β(a)` at finite β vs the
   full asymptotic + resummation picture — the analogue of the cusp's `fp_cusp` ground-truth check
   (`coupled-atlas/_w_groundtruth.py` pattern). Confidence-building; moderate effort.
7. **[far frontier, multi-year] Conjecture** (paper §4): a DIRECT exact-WKB / Voros resurgence of the
   white-noise operator `𝒜_β` itself — the whole PII route (T1–T3) is its "shadow". Nikolaev
   arXiv:2410.17224 is the deterministic backbone to noise-dress; no random-operator exact-WKB framework
   exists. Deep-research REFUTED the naive single-thimble transfer — reduce to the finite-dim ODE first.

## Files

- **Paper (primary):** `TWbeta_Resurgent_paper.tex` (8pp, two-directions framing, T1/T2/T3).
  `TWbeta_Resurgent_skeleton.tex` = earlier 6pp version, **superseded** (kept for history).
- **Notes (`regime-tests/`):** `PROGRAM2_TWBETA_{STOKES_CONSTANT,S_OF_A,XEXTENSION_FACTOR,GY_ZEROMODE,
  THREAD2_BOREL,NOISE_BOREL_PROOF}_NOTES.md`.
- **Scripts (`coupled-atlas/`):** `_hm_stokes_exact.py`, `_hm_stokes_validate.py`, `_hm_coeffs_exact.json`,
  `_tw_stokes_a_dependence.py`, `_tw_mfpt_check.py`, `_tw_frozen_vs_hm.py`, `_tw_zeromode_gy.py`,
  `_tw_thread2_borel.py`, `_tw_noise_borel_proof.py`; falsifier `_o1_falsifier.py`.
- **Working O-notes (fuller derivations):** `O1_reduction_note.tex` (reduction + transport hierarchy),
  `O1_step1b_laxpair.tex` (the Lax map), `O5*.tex` (location, Ω₁, Stokes branch).

## Gotchas (learned this program — don't re-discover)

- **LEVEL vs NOISE — the central pitfall.** Costin's nonlinear-ODE theorem governs the **level** (a)
  direction (PII in `a`), NOT the `1/β` series. The noise ray is proved by the **elementary Airy
  integral + Nevanlinna/median** route. The superseded skeleton conflated them; the paper fixes it.
- **The resurgence is level-variable at β=2** (HM PII, Dunne/Cleri–Dunne). CLDS 2510.14433 treat the
  `1/β` corrections as ordinary **cumulants** (typical/bulk, `C₂≈1.6697=`the program's `s₂`) — that is
  NOT the tail, so no contradiction with Thread 2.
- **Borel-singularity naming:** the dynamical sector's sing sits at `t_c=2/3=A_g=Φ` (reduced action);
  its *median ambiguity* is the two-instanton `e^{−2βΦ}`. Do NOT call the singularity "2Φ".
- **Numerics:** `float(Fraction)` overflows for `k≳50` (use `num/den` strings → mpf); Neville
  extrapolation is ill-conditioned at large windows (use `W≤30` at high dps); compute `Q̃_j` by a
  running mpf product, not exact `(6j)!` (astronomical); median Borel–Laplace via rotated ±θ rays;
  Airy-norm bug: `‖Ai(·−Λ)‖²=Ai'(−Λ)²` via `scipy airy(−Λ)[1]²`, not `ai_zeros`.

## Success criteria

The paper is submittable once (task 1) the tone pass is done and a venue chosen; tasks 2–5 either land
or ship as the already-written remarks; task 6 is optional corroboration; task 7 is the mapped frontier.
The headline is the **two-directional resurgence** (level HM-PII + noise two-sector) of a random-operator
edge law — with every constant certified and the noise ray proved by elementary means.
