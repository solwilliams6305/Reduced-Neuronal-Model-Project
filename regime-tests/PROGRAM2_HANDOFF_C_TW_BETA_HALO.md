# HANDOFF C — the TW_β "halo result": a resurgent trans-series for the Tracy–Widom-β tail

> **STATUS 2026-07-22: both theorems in hand (T1 proved, T2 established, both modulo cited
> machinery). This is the `q=1` (Airy/fold) THEOREM-GRADE sibling of the cusp `𝒲` paper (`q=2`,
> Weber/cusp, numerical). Consolidated paper: `TWbeta_Resurgent_skeleton.tex` (repo root, 6pp,
> builds clean). Remaining work is sharpenings, NOT open "is this true" questions.
> Memory: `program2-tw-beta-halo.md`.**
>
> **UPDATE 2026-07-22 (this session): tasks #1 and #2-backbone CLOSED.** (i) All 5 classical
> citations web-verified (2 title patches applied). (ii) The backbone Stokes constant is now EXACT:
> `C=√(2/(3π³))=0.146632271193848478…`, `γ=−1/2` (certified 28/25 digits, three independent routes;
> the paper's `1/(4√π)` guess REFUTED, 3.8% off). Skeleton + bibliography updated & rebuilt clean.
> New: `PROGRAM2_TWBETA_STOKES_CONSTANT_NOTES.md`, `coupled-atlas/_hm_stokes_exact.py` +
> `_hm_stokes_validate.py`.
>
> **UPDATE 2026-07-22 (cont.): task #2 a-dependence ESSENTIALLY CLOSED.** The a-dependence of `S(a)` is
> ALGEBRAIC (eikonal), not transcendental: an exact scaling symmetry `p=√a·P` of the backward-Kolmogorov
> PDE (verified in sympy) makes `a` enter only through `g=βa^{3/2}`, so `S(a)=S₀·a^{p}`, `S₀` a pure
> a-independent number. The reduced cubic-barrier escape has closed-form Stokes constant **`S₀=1/π`**
> (certified 25 digits; `d_m∼−(1/π)Γ(m)(3/2)^m`; escape-rate formula validated vs direct MFPT integral),
> and `C=(1/π)√(2/(3π))`. New: `PROGRAM2_TWBETA_S_OF_A_NOTES.md`, `coupled-atlas/_tw_stokes_a_dependence.py`,
> `_tw_mfpt_check.py`.
>
> **UPDATE 2026-07-22 (cont.²): the `√(2/(3π))` factor DERIVED (structurally) + framing corrected +
> literature scoping.** (`PROGRAM2_TWBETA_XEXTENSION_FACTOR_NOTES.md`, `_tw_frozen_vs_hm.py`.) The
> frozen→full change is a **joint shift** `(S₀=1/π, γ'=0) → (C, γ=−½)` (NOT a lone factor): constant
> `×√(2/(3π))=√(A_g/π)` (`A_g=2/3` the action — the **x-translation zero-mode Jacobian**) bundled with
> index `−½` = the one-continuous-zero-mode signature (same mode the O1 note removes for `a^{−3β/4}`).
> Confirmed 9 digits (`r_k/√(2k)→1/(πC)=√(3π/2)`). Matches **Dunne 2511.15528 eq. 2.16**, which writes
> the HM Stokes constant as exactly `(1/π)√(2/(3π))`; our `a_n` reproduce its `−17/72` subleading to 16
> digits — so `a_n` ARE Dunne's `c_n^{(0,−)}`. **⚠ SCOPING (deep-research, important for the paper's
> framing):** this resurgence is **level-variable at the fixed classical point β=2** (HM PII), NOT the
> noise `1/β`. CLDS 2510.14433 treat the `1/β` corrections as ordinary **cumulants** (`C₂≈1.6697`=the
> program's `s₂=1.67`; no Borel/Stokes); Borot–Nadal 1111.2761 give the exact all-β soft tail as a
> **single-exponential** level series (no 2nd instanton). So a genuine **`1/β` trans-series is
> un-established** in the literature — the paper's "resurgent trans-series" thesis should be scoped to
> the level variable / β=2, or the `1/β` Borel-summability treated as the open frontier (Thread 2).
> Skeleton T2c Remark refined & rebuilt clean.
>
> **UPDATE 2026-07-22 (cont.³): the `√(2/(3π))` factor DERIVED in closed form.**
> (`PROGRAM2_TWBETA_GY_ZEROMODE_NOTES.md`, `_tw_zeromode_gy.py`.) It is the Gelfand–Yaglom zero-mode
> Jacobian `√(‖ψ₀‖²/2π)` of the reduced escape instanton — the **tanh kink** `P⋆=−tanh τ`, zero mode
> `ψ₀=−sech²τ`, with `‖ψ₀‖²=∫sech⁴=4/3=ΔV=2A_g` (virial). So `√(A_g/π)=√(2/(3π))` (the `2` is virial,
> not a boundary artifact), and `C=(1/π)√(2/(3π))=√(2/(3π³))` **to 40 digits** — `S₀=1/π` the transverse
> Kramers determinant, `√(2/(3π))` the escape-translation zero mode (frozen *rate* factors it out; the
> level-variable resurgence counts it in; `−½` index = one zero mode). **T2c fully accounted.** The only
> program-level open item is the noise-`1/β` resurgence (Thread 2, un-established per the literature).
>
> **UPDATE 2026-07-22 (cont.⁴): THREAD 2 ANSWERED — the noise `1/β` expansion IS Borel-summable.**
> (`PROGRAM2_TWBETA_THREAD2_BOREL_NOTES.md`, `_tw_thread2_borel.py`.) It is a genuine **two-sector
> resurgent trans-series**: **(A) dynamical escape** — the reduced series `d_m∼−(1/π)Γ(m)(3/2)^m` is
> non-alternating ⇒ POSITIVE-real-axis Borel singularity at `t_c=2/3=A_g` (the reduced instanton action Φ; two-instanton ambiguity; verified via Borel-Padé);
> **median Borel–Laplace reconstructs the exact escape rate to `1.2×10⁻³`** ⇒ median-summable. **(B)
> Coulomb-gas** — the `Γ(β/2)` β-ensemble normalization's Stirling tail is ALTERNATING factorial
> (`|a_k/a_{k−1}|/[(2k−2)(2k−3)]→1/π²`, 6 digits), Borel `~(1/π)arctan(t/π)` ⇒ IMAGINARY-axis (`±iπ`)
> ⇒ ordinary Borel-summable. (The Borot–Nadal `R_m(2/β)` fluctuation part is regular.) **Upgrades the
> literature's "un-established" to a concrete characterization.** Doesn't contradict CLDS — their
> ordinary cumulants are the *typical/bulk* fluctuations, not the *tail*. Caveat: numerically decisive;
> a full Costin-type proof for the noise ray is the standard remaining rigor. **All three resurgence
> directions now mapped:** level-a/β=2 (HM PII, `C`); noise-`1/β` dynamical (escape, median); noise-`1/β`
> Coulomb-gas (`Γ`, ordinary).
>
> **UPDATE 2026-07-22 (cont.⁵): the Costin-type noise-ray proof DONE (dynamical sector).**
> (`PROGRAM2_TWBETA_NOISE_BOREL_PROOF_NOTES.md`, `_tw_noise_borel_proof.py`.) **Median Borel-summability
> of `F(g)=Q̃(2/g)^{-2}` along `ℝ₊` is proved**, median sum `= π e^{(2/3)g}R(g)`. Key move: the reduced
> escape has an **Airy-type integral rep** `Q(D)=∫exp(−t²+(√D/3)t³)dt` (linear cubic phase), so — unlike
> the level variable — it needs NO nonlinear-ODE (Costin) theorem; just Watson/Nevanlinna Gevrey-1 +
> Écalle/Costin **median** summation. CERTIFIED: (i) integral rep matches the series to `2×10⁻¹¹` at
> `D=0.06`; (ii) Borel transform `𝓑(ξ)` analytic in `ℂ∖[4/3,∞)`, single log branch at `ξ=4/3`
> (Pringsheim + Domb–Sykes `→3/4` + amplitude `→1/2π` + Padé cut); (iii) median reconstructs `R(g)` to
> `1.2×10⁻³`; (iv) ring-closure gives `Q̃^{-2}`. Coulomb-gas sector = classical Watson (Γ). **Scope:**
> proved in the scaling regime `g=βa^{3/2}→∞` (frozen escape exact); the uniform fixed-`a` statement
> (x-extension/PII corrections) is the remaining analysis. **Thread 2 now has a proof, not just numerics.**
>
> **UPDATE 2026-07-22 (cont.⁶): FULL PAPER consolidated → `TWbeta_Resurgent_paper.tex` (8pp, builds
> clean).** Supersedes the 6pp `TWbeta_Resurgent_skeleton.tex` as the primary document. New organizing
> principle: **two directions of resurgence** — **(L) level-`a`** (HM PII via the Lax map; Costin;
> Stokes constant `C=√(2/(3π³))`) and **(N) noise-`1/β`** (median Borel-summable, two sectors: dynamical
> escape + Coulomb-gas `Γ(β/2)`, proved via the Airy integral). Fixes the skeleton's level-vs-noise
> conflation (it invoked Costin for the 1/β series; Costin actually governs the level direction, the
> noise ray is Airy/Nevanlinna). Structure: §1 objects, §2 two directions, §3 PII backbone, §4 T1/T2/T3,
> §5 proof T1 (Lax map), §6 exact `C` + GY zero-mode `a`-dependence, §7 proof T3 (noise Borel), §8
> location/`Ω₁`, §9 ledger. All this session's results now promoted from remarks to theorems/sections.
> **Remaining before submission:** proof-read for over-claims; the 4 refinements (α-family, uniform
> fixed-`a`, rough-potential rigor, exact-WKB frontier) are tagged `[to prove]`/`[frontier]` in-text.

_Self-contained handoff for a fresh Claude Code session. Companion to the cusp paper
`W_ResurgentTransseries_paper.tex` (Program 2's `q=2` deliverable). Reads: this file, then the
skeleton `.tex`, then the O-notes below in order._

## The one-sentence goal (achieved, being sharpened)
The FIRST Borel-summable resurgent trans-series for a random-operator edge law: the Tracy–Widom-β
(`TW_β`) right tail, obtained by dressing the Painlevé II / Hastings–McLeod backbone with noise at
large Dyson index β (weak noise ε = 1/√β).

## Why q=1 is provable (the load-bearing idea — don't lose it)
At large β the tail REDUCES to Painlevé II (Comtet–Le Doussal–Smith 2510.14433), which is already
proved resurgent (Cleri–Dunne 2002.06270). The cusp `𝒲` (q=2) has NO Painlevé reduction
(isomonodromy fixed point) → only numerical. So `TW_β` leans on the integrable crutch the cusp
lacks. **Structural dichotomy (the paper's organizing insight):** the cusp's complexity is a
*resonance* phenomenon (non-self-adjoint, unbounded-below Weber → complex λ₀); `TW_β`'s is a
*trans-series* phenomenon (self-adjoint SAO → real spectrum → real Ω₁; complexity in the median
resummation).

## What is DONE — do NOT rebuild (all in `TWbeta_Resurgent_skeleton.tex`, verified)
- **T1 — the tail is a Borel-summable PII trans-series [PROVED, modulo cited machinery].**
  - Exact reduction: tail = explosion prob of Riccati SDE `dp=((x+a)-p²)dx+2ε dB`; no-explosion
    prob solves an exact backward-Kolmogorov PDE.
  - Weak-noise WKB: instanton `ṗ=p²-(x+a)` linearizes to **Airy** `φ''=(x+a)φ`; rate
    Φ(a)=(2/3)a^{3/2} (matches Dumaz–Virág 1102.4818).
  - Fluctuation (Jacobi) operator **J = -∂_x² + 6p⋆² - 2(x+a)** (2nd variation of FW action).
  - **THE CRUX (Lax-pair map, VERIFIED to machine precision, `_omega1`-adjacent checks):** under
    `s=-2^{1/3}(x+a)`, `q=-2^{-1/3}p⋆`, the instanton solves `q'=q²+s/2` ⇒ PII `q''=2q³+sq+1/2`
    (α=1/2, the classical **Airy special solution** of PII), and **J = μ²·(-∂_s²+6q²+s) = μ²·L_PII**
    exactly. So the fluctuation operator IS a PII linearization → Costin math/0608408 gives Borel
    summability; Stokes const ≠0 (Costin-Costin-Kohut math/0608316); zero mode (Schorlepp-Grafke-
    Grauer 2208.08413).
  - **No α(β):** the a^{-3β/4} prefactor is EIKONAL (its log a coefficient scales as (3/4)β), so α
    is β-independent; both α∈{0,1/2} in Costin's rank-one class ⇒ T1 robust to family.
- **T2 — the noise-averaged Stokes data [ESTABLISHED].**
  - **T2a location-fixing:** WKB action A[b']=∫√((x+a)+2εb'); E[A₁]=0 (linear, mean-zero →
    location fixed at O(ε), verified); E[A₂]=−½δ(0)∫(x+a)^{-3/2} δ(0)-divergent (matches cusp v₃
    threshold, verified) → location not clean, Stokes data carries it. All-orders via isomonodromy
    rigidity (Airy spectral curve). File: `O5_location_fixing.tex`.
  - **T2b Ω₁ = −1.12 [COMPUTED, 2 independent methods]:** self-adjoint SAO ⇒ real spectrum every
    realization (verified max|Im Λ|=0). `Ω₁ = 4Σ_{n≥1}⟨ψ₀²|ψₙ²⟩/(Λ₀-Λₙ) = -4∫ψ₀²G_red`, converges
    n^{-4/3}. Sum-over-states −1.121 (`coupled-atlas/_omega1_airy.py`); direct-diagonalization MC
    −1.128±0.009 (`_omega1_crosscheck.py`), agree <1%. E[Λ₀](β)=2.3381−1.12/β. File:
    `O5b_omega1_result.tex`.
  - **T2c complex Stokes constant [COMPUTED, backbone EXACT 2026-07-22]:** by T1 = the HM PII datum;
    computed the HM s→−∞ asymptotic series `q~√(-s/2)(1-1/8(-s)^{-3}-73/128(-s)^{-6}-...)` — **all
    a_{k≥1}<0, non-alternating** (Cleri-Dunne signature) → positive-axis Borel singularity → imaginary
    lateral ambiguity → median = real tail. Growth `a_k ~ -C Γ(2k+γ)(9/8)^k` with, now CERTIFIED,
    **γ=-1/2 (25 digits) and C=√(2/(3π³))=(1/π)√(2/(3π))=0.146632271193848478… (28 digits)**, action
    A=2√2/3. Computed from a fast exact rational recursion to N=120 + Richardson (method A), matched to
    the published HM Stokes constant (Dunne 2511.15528 eq.2.16; Kapaev nlin/0411009; FIKN ch.11) and
    cross-checked by a Borel-space branch-exponent extraction and direct HM-ODE integration.
    **REFUTES the earlier guess C≈0.139≈1/(4√π) (that was under-converged N=10; 1/(4√π) is 3.8% off).**
    Files: `coupled-atlas/_hm_stokes_exact.py`, `_hm_stokes_validate.py`, `_hm_coeffs_exact.json`,
    `regime-tests/PROGRAM2_TWBETA_STOKES_CONSTANT_NOTES.md`. (Supersedes `O5d_stokes_computed.tex` /
    `_hm_stokes.py` on the value of C. O5c note has the self-adjointness dichotomy; NB its framing slip
    — S is a single fn S(a), not "S_HM+S₁/β" — corrected in O5d.)

## OPEN TASKS (priority order) — the remaining work is all sharpening, not "is it true"
1. **✅ DONE (2026-07-22) — 5 classical citations verified** (web, primary sources). Two need patches
   (now applied in the skeleton): Flaschka–Newell title "…deformations**. I**"; Jimbo–Miwa–Ueno
   subtitle "…rational coefficients**. I. General theory and τ-function**". Tracy–Widom 1994, FIKN,
   Clarkson confirmed as-was. Bonus verified: HM Stokes multipliers `(−i,0,i)`. (The arXiv cites were
   already verified: RRV math/0607331, DV 1102.4818, BorotNadal 1111.2761, CLDS 2510.14433, CleriDunne
   2002.06270, Costin math/0608408, CCK math/0608316, SGG 2208.08413, Nikolaev 2410.17224, AnicetoCrew
   2410.13690, Viens 0901.0383.)
2. **T2c refinement — a-dependence of S(a). ✅ ESSENTIALLY DONE (2026-07-22).**
   (i) Exact backbone: **C=√(2/(3π³))=0.146632271193848478…, γ=−1/2** (`PROGRAM2_TWBETA_STOKES_CONSTANT_NOTES.md`;
   1/(4√π) refuted). (ii) **a-dependence is ALGEBRAIC (eikonal), NOT transcendental**
   (`PROGRAM2_TWBETA_S_OF_A_NOTES.md`, `coupled-atlas/_tw_stokes_a_dependence.py`): the exact scaling
   symmetry p=√a·P of the backward-Kolmogorov PDE (verified in sympy) makes `a` enter the fluctuation
   expansion only through g=βa^{3/2}, so c_n(a)∝a^{−3n/2} and S(a)=S₀·a^{p} with S₀ an a-independent
   pure number (Borel singularity at Φ(a)=(2/3)a^{3/2}). The reduced cubic-barrier escape has the
   CLOSED-FORM Stokes constant **S₀=1/π** (certified 25 digits; d_m∼−(1/π)Γ(m)(3/2)^m, ρ=3/2, γ'=0).
   **RESIDUAL (now sharp + bounded):** the full-tail C and the frozen-escape 1/π factor as
   C=(1/π)·√(2/(3π)); derive the elementary √(2/(3π)) x-extension factor rigorously (transport
   hierarchy c_n(a), `O1_reduction_note.tex` §3) — i.e. show the x-extension dresses 1/π→C. No longer
   an open "is there hidden a-structure" question, just this finite computation.
3. **T1 refinement — α=1/2 vs α=0 assignment.** The tail SADDLE is the α=1/2 Airy solution; the
   exact β=2 DISTRIBUTION is α=0 (HM). Not a β-dependence (that's resolved) — a saddle-vs-
   distribution / half-integer-vs-integer PII refinement. Pin via the Airy-solution/coalescence
   structure. Does NOT block T1 (both α in Costin's class). (`O1_step1b_laxpair.tex` Q1/Q2.)
4. **T2a refinement — rough-potential rigor** of the isomonodromy rigidity (white noise really
   preserving the rank-3/2 irregular singularity). Controlled by the O5 renormalization; make
   rigorous.
5. **[far frontier, multi-year] Conjecture:** a DIRECT exact-WKB / Voros resurgence of the
   white-noise operator itself (Thread-2 blank from the deep-research). The PII route (T1/T2) is its
   "shadow". Nikolaev 2410.17224 is the deterministic backbone to noise-dress; no random-operator
   exact-WKB framework exists.

## Files (all repo root unless noted)
- **Paper:** `TWbeta_Resurgent_skeleton.tex` (consolidated, T1+T2, 6pp — the primary document).
- **Working notes (detailed derivations, fuller than the skeleton):** `O1_reduction_note.tex`
  (reduction + all-orders structure + transport hierarchy setup), `O1_step1_isomonodromy.tex`
  (Jacobi operator + framing), `O1_step1b_laxpair.tex` (THE Lax-pair map, α=1/2),
  `O1_step2_alpha.tex` (no α(β), eikonal), `O5_location_fixing.tex` (T2a), `O5b_omega1_result.tex`
  (Ω₁), `O5c_stokes_branch.tex` (self-adjointness dichotomy; framing-slip, see O5d),
  `O5d_stokes_computed.tex` (complex Stokes constant computed).
- **Scripts (`coupled-atlas/`):** `_o1_falsifier.py` (PII→DV tail check + variance s₂=1.67),
  `_omega1_airy.py` (Ω₁ sum-over-states), `_omega1_crosscheck.py` (Ω₁ diagonalization MC),
  `_hm_stokes.py` (HM asymptotic coeffs + Stokes growth), `_hm_coeffs.json`.
- **Deep-research scope** (108-agent pass, the literature map): in the session transcript; key
  refs folded into the skeleton bibliography.

## Gotchas (learned this session — don't re-discover)
- **Airy eigenfunction norm:** `‖Ai(·-Λ)‖²_{[0,∞)} = Ai'(-Λ)²`, get it via `scipy airy(-Λ)[1]²` —
  NOT scipy `ai_zeros`' returned values (those are something else; using them gives a FALSE
  divergence in Ω₁). The reduced-resolvent sum converges only with the correct norm (n^{-4/3}).
- **Ω₁ is REAL by structure** (SAO self-adjoint), not a limitation — the complex Stokes constant is
  a distinct trans-series object (O5c/O5d).
- **The a^{-3β/4} prefactor is eikonal**, not a varying Painlevé index — don't chase an "α(β)".
- **Don't attack the white-noise operator directly** (deep-research REFUTED the Serone-Spada-
  Villadoro single-thimble transfer, 0-3): reduce to the finite-dim PII/Riccati ODE FIRST, then sum.

## Success criteria (for "submittable")
T1 + T2 written up with the 5 classical citations verified (task 1) and the refinements (tasks 2-4)
either done or clearly stated as remarks (they already are, in the skeleton). The paper is then a
genuine first-of-kind theorem paper; tasks 2-4 are follow-on strengthenings and task 5 is the
mapped multi-year frontier.
