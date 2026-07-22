# Derivation audit — working identities that lack a first-principles derivation

**Purpose.** A project-wide inventory of the constants, laws, and mechanisms the
chapters *use* but have **not** derived from first principles — i.e. the values
that are measured, fitted, cited from the literature, computed numerically, or
rest on a heuristic. This is the "Path-A frontier" of the whole project in one
place. It is the honest answer to "what is asserted vs what is proven."

**The pattern, stated once.** Across every chapter the *laws and exponents* are
analytically derived (blow-up + accumulated-variance + matched-asymptotic
arguments); the *O(1) constants* are almost all measured. The gaps below are
overwhelmingly **constants**, plus a few cited theorems and one refuted heuristic.

**Status tags:** **MEAS** = fitted/measured constant, no derivation · **CITE** =
from the literature, not derived in-project · **NUM** = computed numerically, no
closed form · **HEUR** = scaling/heuristic argument, not rigorous (and has been
wrong) · **PART** = mechanism/form derived, value not.

---

## The audit

| # | Identity | Chapter | Role | Status | What's missing for a first-order derivation |
|---|---|---|---|---|---|
| 1 | **C_q = √(4π ln2) ≈ 2.95 (NF)** | Canard (→ MMO, VdP) | the universal noise prefactor σ_* = C_q√ε·(geom) | **DERIVED** (leading order; `DERIVATION_STRATEGIES.md` §1) — *verified independently*: barrier (4/3)W^{3/2}, hazard collapse H=η²/(4πλ), C_q=√(4π·H_crit) all check out | **validation not yet run:** the measured 2.8 used the *arbitrary* V_cross=1 detector — re-measure with the *intrinsic* separatrix (V=+√W) and test the quantile law C_q(p)=√(−4π ln(1−p)). The 2.95→2.8 gap is the adiabatic/finite-barrier Kramers correction (sign argued, not proven). Closes VdP universal half (#16). |
| 2 | **κ ≈ 2π²** | MMO | folded-node rotation map `ln R = κμ` | **reduced to one integral; prefactor OPEN** (`MMO_K2_ROUTE_AB.md`) | Route B (*verified*): antidamping reduction gives `κ=2π² ⇔ ∫_funnel g dτ = π²`, and **rules out** the conservative WKB envelope (`κ_env~(2/3)ln(1/μ)` — wrong scaling, 10× small; re-derived ✓). The order `c₀=1/μ−1` checks via WKB=s_max. The remaining `2π²` = the parabolic-cylinder **modulus connection coefficient** (DLMF §12.14) — same computation as `a_min` (#3). *Caveats:* the Liouville potential has a sign slip (`Q=Ω²−g²+g'`, not `−g'` — leading-order-immaterial but it lives in the O(μ) drift term); coefficients g, Ω, a↔μ unverified vs Wechselberger §4. **Decisive next step: measure ∫g numerically.** |
| 3 | **a_min(c)** | MMO | smallest-SAO amplitude; sets f(c), α, and the global σ_pq μ-law | **NUM** | closed form from the **Krupa–Wechselberger 2010** transition-regime inner solution (µ < ε^{1/2}); the standard O(ε^{(1−µ)/2}) sector scaling is an upper bound only |
| 4 | **global μ^{3/2} in σ_pq** | MMO | the c/q-dependence of the dissolution threshold | **HEUR/empirical** | the local "sector-crossing" derivation was **refuted** (crossover); the literature sector-packing route `µ·ε^{(1−µ)/2}` is found but needs KW 2010. The power may not even be 3/2 (M2 degeneracy) |
| 5 | **γ ≈ 1** | MMO | rotation-spread exponent δρ ~ σ^γ | **HEUR**, unmeasured | both a *direct measurement* (std(ρ)~σ^γ) and a first-principles value; currently a Path-B guess, not to be back-inferred from α/β |
| 6 | **realised α ≈ 1.55** | MMO | the devil's-staircase exponent | **MEAS** (composition) | rests on measured κ (#2) + numerical a_min (#3); a closed-form a_min ⇒ closed-form α |
| 7 | **s_max = (1−µ)/(2µ)**, secondary-canard count | MMO | the SAO ceiling, α_ceiling = 2 | **CITE** (Wechselberger 2005) | not derived in-project — and applied *below* its stated µ ≫ ε^{1/2} regime (§2.5) |
| 8 | **c ≈ 1.55** (`A_mid = √(c/π²)`) | Tonic | the mid-tonic jitter constant | **PART** | mechanism (fold blow-up of adjoint) + inner eq `−2V` **are** derived; but the value overshoots/drifts (`c_recon ≈ 1.9→4.5`), the bulk shoulder (~36% of ∮Z_v²) is numerical. Needs the coupled O(ε^{1/3}) BVP + matched bulk integral |
| 9 | **σ_onset ≈ 0.015–0.02** | Tonic | onset of non-perturbative fold escape | **MEAS** | the Freidlin–Wentzell *cycle*-escape action (the tonic instance of σ_*); flagged PhD-scope |
| 10 | **instanton commitment −1.06 rad** | Resonator | spike-direction selection | **NUM** (gMAM) | a closed-form minimum-action escape path; only a numerical value exists |
| 11 | **measured B-exponent ≈ δ^{0.5}** | Resonator | the barrier-vs-deficit slope | **PART/explained** | the δ² *corner law* IS derived from the subcritical-Hopf normal form; the measured pre-asymptotic exponent is floor-/global-dominated — needs the full 2D quasipotential beyond the Hopf normal form |
| 12 | **C = K_fold/A_0** | Excitable | the σ_crit = √(ε/C) prefactor | **PART** | stated "computable from geometry, no free parameters," but K_fold, A_0 are not actually computed to a validated number |
| 13 | **τ_v = b/([1−b(1−V_FP²)]·\|1−V_FP²\|) ≈ 0.95** | Excitable | OU timescale in the w_escape model | **DERIVED** (leading order; `DERIVATION_STRATEGIES.md` §13) — *verified independently*: the 2×2 Lyapunov solution Var(δw)=a₂₁²σ²/(2·detJ·\|trJ\|) is exact (re-derived from scratch) | the old 1-D `1/\|1−V²\|=1.73` was the wrong leading order; the 2-D Lyapunov transfer gives ≈0.95, bracketing the measured 1.16 from below (residual ~18% = nonlinear covariance). **Check:** the number hinges on V_FP=−1.2563 (confirm vs kernel.py — if V_FP≈−1.2 it brackets from *above*); test Var∝εσ² and the τ_v(ε) trend. |
| 14 | **explosion washout ~ ε^{3/4}** | Canard | small-cycle noise threshold | **CITE** (BGK 2012) | not derived in-project; experimentally inaccessible in FHN |
| 15 | **ε^{−1/4} Hopf-edge prefactor** | Tonic | edge divergence of A(I,ε) | derived but **unvalidated** | the (I−I_H)^{−1/2} part is derived + confirmed; the ε^{−1/4} is derived but numerically inaccessible — no empirical check |
| 16 | **c_vdp ≈ 0.78** (C_q^vdp now via #1) | Cross-model (VdP) | universality prefactors | **MEAS** (C_q half **DERIVED** via #1) | VdP shares the identical blow-up `dX/dT=X²−Y`, so its normal-form C_q is the *same* √(4π ln2) (#1) — universal half closed. Only `c_vdp` (the tonic-style A_mid prefactor) still needs the VdP inner BVP (cf. #8). |

---

## For contrast — what *is* first-principles (so the gaps are visible against it)

Derived and not in the table: the canard exponent `ε^{1/2}λ^{1/2}` (accumulated
Brownian variance); the tonic form `A_mid = √(c/π²)` and the corrected inner
equation `−2V`; the tonic edge `(I−I_H)^{−1/2}` (Green's-theorem topology); the
resonator first Lyapunov coefficient `a_3 = +0.268` / `l_1 = +0.97` and the `δ²`
corner law (Kuznetsov Hopf formula + amplitude potential); the MMO `µ ∝ (c+1)`
(Jacobian) and `α_ceiling = 2`; the **crossover** `µ → 0` reduction of the
folded-node normal form to the canard form; the excitable `ε^{1/2}` asymptotic
exponent. **The skeleton is derived; the constants are measured.**

---

## Reading the audit — where the leverage is

1. **C_q (#1) is the single highest-leverage derivation.** It is the *same*
   constant in canard, MMO (via the crossover), and VdP. Deriving it once — the
   prefactor of the blow-up escape integral — upgrades three chapters at once.
2. **κ = 2π² (#2) is the most self-contained.** The value almost certainly has a
   clean origin in the folded-node normal form; every step is numerically
   checkable against the measured κ. Best candidate for a short, standalone,
   genuinely-yours derivation.
3. **a_min(c) (#3) is a two-for-one.** It is the shared unknown behind the
   deterministic α (#6) and the global noise μ-law (#4); one closed form (KW 2010
   transition regime) closes both — but it is the hardest (special-function
   inner solution, below the standard regime).
4. **The honest tiering for the thesis:** #1, #2, #8 are *constants* of derived
   laws — closable, high-value, "asymptotic theory" register. #3, #4, #10 need
   special-function inner solutions — genuinely PhD/paper-scope. #7, #14 are
   *cited* results to attribute, not re-derive. #13 is a *flag* (a derived value
   that doesn't match — worth a line of honesty).

**One-line summary:** nothing in the project rests on an *unstated* assumption —
every gap above is already flagged "open/measured/cited" in its chapter; this
table just collects them. The work is sound *as scaling theory validated by
numerics*; the table is exactly the list of derivations that would move it to
*asymptotic theory*.
