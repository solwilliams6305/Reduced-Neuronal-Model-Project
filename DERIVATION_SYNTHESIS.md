# Derivation synthesis — the lay of the land

*A capstone over the `regime-tests/` derivation work. Start here; it points to the
individual documents. Companion to `regime-tests/DERIVATION_AUDIT.md`, which lists
the 16 gaps this work addresses.*

> **Revision note.** Novelty claims in an earlier draft were too strong and were
> corrected after a literature check and an internal-consistency review (see
> **Prior art and honest positioning** below). In particular the circle-map
> "criticality" narrative is withdrawn — it contradicted the project's own
> Phase-1.5 retraction — and the entire noise-MMO cluster (#3/#4/#6) overlaps the
> published Berglund–Gentz–Kuehn programme and carries **no** novelty claim
> pending a full read of it.

---

## The one-paragraph picture

The audit found 16 places where a scaling **law** was derived but its leading
**constant** was only measured. The work since then shows those constants come in
exactly two kinds: ones that **collapse to a clean closed form** (the
fold-escape/quasipotential family) and ones that are **irreducibly global** (least-
action paths, global returns) and stay numerical. Knowing which is which — and
*why* — is the main result. Three chapters' noise thresholds turned out to be one
formula; two more chapters turned out to share a single Airy function; the hardest
corner (folded-node MMOs) collapsed from four constants to one integral; and two
genuine boundaries were proven to be numerical for understood reasons.

---

## Prior art and honest positioning (read this before claiming novelty)

This work sits in a **populated field**, and the closest prior art is close to home:

- **Berglund, Gentz, Kuehn — *From random Poincaré maps to stochastic mixed-mode-
  oscillation patterns*** (arXiv:1312.6353). Noise on 1-fast/2-slow folded-node
  systems; a Markov chain on returns to a Poincaré section; **explicitly the
  interplay between the SAO count and the global return, a saturation phenomenon
  near the folded node, and conditions for when noise increases/decreases the SAO
  count.** That is the qualitative core of the whole `#3/#4/#6` noise-MMO cluster
  here. Kuehn co-authors with Popović. **Nothing in the noise-MMO cluster should
  be called novel until this is read in full**; at best the `Λ(μ)` consolidation
  is a repackaging within their framework.
- **Sacré & Franci — *Singularly perturbed phase response curves*** (arXiv:1506.02010),
  and Sacré's 2013 PhD thesis (Liège). FHN PRCs in the singular limit (same
  `a=0.7, b=0.8`), via a geometric finite-PRC / isochron approach. Their thesis is
  that the *infinitesimal* PRC breaks near the singular limit — exactly the fold
  regime the §B Airy duality addresses. Position the duality as the complementary
  inner-adjoint structure, and check the 2013 thesis for `iPRC = Ai²` before
  claiming it.
- **Berglund & Gentz — *Noise-Induced Phenomena in Slow–Fast Dynamical Systems***.
  Canonical prior art for the fold-escape family (§A); the `√(4π ln2)` collapse is
  a clean closed form *within* this framework, not a new phenomenon.
- **Desroches et al. — *Mixed-Mode Oscillations with Multiple Time Scales*** (SIAM
  Review); Wechselberger 2005; Krupa–Szmolyan 2001 — the folded-node/canard
  machinery every inner result rests on.

**The honest framing for the whole body of work:** clean closed forms and
consolidations *within* known frameworks (Freidlin–Wentzell / Berglund–Gentz–
Kuehn, Wechselberger, Krupa–Szmolyan) — a legitimate contribution, but not new
phenomena, and the noise-MMO side must be positioned against BGK first.

---

## The landscape

### A. The fold-escape family — closed (leading order)

The biggest consolidation: the **canard** noise threshold (#1), the **excitable**
threshold (#12), and the **Van der Pol** cross-check (#16) are the *same*
Freidlin–Wentzell fold-escape problem, and collapse to one formula

```
σ_crit = √( 4π · D · ln2 ) · √ε ,     D = deterministic slow-drift speed AT the fold,
```

with `D = λ = b(I−I_fold)` for the canard and `D = g = |a−1−b(I−2/3)|` for the
excitable; the `√(4π ln2)` is universal. The canard constant is parameter-free,
`C_q = √(4π ln2) ≈ 2.95` (measured ~2.8); the excitable barrier constants come out
of the cubic geometry exactly, `K_fold = 4/3`, `A₀ = 1/π`. The mechanism is the
exact Riccati-cubic barrier `(4/3)W^{3/2}` plus a hazard integral whose Kramers
prefactor cancels the Jacobian and collapses to `η²/(4πλ)`. **Status:** clean
leading order, with named caveats (the Kramers prefactor is not asymptotically
exact near the dominant region; an O(λ) adiabatic correction). → `FOLD_ESCAPE_PREFACTORS.md`.

### B. The tonic–canard duality — an organizing lemma (prior art unchecked), structure closed

The tonic spiking "jitter" constant (#8) and the canard escape are two faces of
**one Airy function** of the fold. The Riccati substitution `V = −u'/u` turns the
fold normal form into the Airy equation; then

```
canard velocity      V   = −u'/u      (Airy log-derivative)
tonic inner iPRC      Z̃_v =  u²  = C·Ai(−ξ)²   (the square)
```

so the iPRC peaks at the first Airy extremum and amplifies across the fold by
`exp(ΔΦ/λ)` — the *same* quasipotential as `C_q`. The constant `c ≈ 1.55` itself
still needs a matching step (an Airy fourth-moment with cutoffs, which is *why*
`c` is model-specific), but the inner solution and the cross-chapter unification
are exact. **It is a one-to-two-step consequence of two classical facts** — the
Riccati→Airy linearization of the canard, and "the adjoint of a scalar
log-derivative equation is its square" — so a folded-node expert would call it the
natural structure, *not* a surprise. It lives in Sacré & Franci's active area
(FHN PRCs in the singular limit; their point is the infinitesimal PRC breaks
exactly here); the `iPRC = Ai²` identity must be checked against Sacré's 2013
thesis before any novelty claim. **Status:** clean organizing lemma; the prize
`c ≈ 1.55` still open. → `TONIC_CMID_AIRY.md`.

### C. The MMO cluster consolidation (circle-map "criticality" withdrawn; overlaps BGK)

The folded-node MMO side (smallest oscillation `a_min` #3, staircase exponent `α`
#6, noise-dissolution `σ_pq` #4). What genuinely survives:

1. **Exponential asymptotics** gave `a_min = a_max · e^{−π² f}` — the inner
   physics is one number (`π² = κ/2`) plus a universal floor `e^{−π²} ≈ 5×10⁻⁵`,
   and the "mysterious exponent" being chased was a red herring (`a_min`'s whole
   `c`-dependence is the global funnel-filling `f`). → `MMO_AMIN_TIER1.md`.
2. **The consolidation:** everything reduces to one passage integral
   `Λ(μ) = ln(a_max/a_min)` — `a_min = a_max e^{−Λ}`, `f = Λ/π²`, `s_obs = Λ/(κμ)`.
   That algebraic collapse is the keepable result.

**Two corrections to the earlier draft.** (i) The circle-map **"criticality"
narrative is withdrawn**: it contradicts the project's own Phase-1.5 retraction
(`MMO_PHASE1_5.md`: *"the FHR staircase is **not** a critical circle map — its
plateaus are organised by the folded node's secondary canards"*), and the measured
`α ≈ 1.55` sits *below* the critical-circle-map range `[2,3]`, which is evidence
*against* criticality. (ii) `α = 1 + 1/(1+r)` is just the existing `α = 1 + 1/p`
relabelled (`p = 1+r`) — not new structure. (iii) The whole noise-MMO content —
SAO-count vs global return, folded-node saturation, noise changing the SAO count —
is the published **Berglund–Gentz–Kuehn** programme (arXiv:1312.6353); the `Λ(μ)`
reduction is at most a repackaging within it. **Status:** the `Λ(μ)` reduction
stands as a consolidation; no novelty claim pending a BGK read.
→ `MMO_AMIN_TIER1.md` (solid), `MMO_AMIN_TIER2.md` / `MMO_CIRCLEMAP.md` (read with
this caveat).

### D. A fix, not a gap — `τ_v` (#13)

The one derived value that *contradicted* measurement (`τ_v`: derived 1.73 vs
measured 1.16). Resolved: the old value was a 1-D estimate; the correct object is
the 2-D Lyapunov noise-transfer of the linearized fixed point, `≈ 0.95`, which
brackets the measurement (the fixed point is a spiral, so it decorrelates faster
than the 1-D guess). **Status:** closed-form leading order. → `DERIVATION_STRATEGIES.md` §#13.

### E. The two honest boundaries

- **Resonator commitment angle `−1.06 rad` (#10).** Proven to be an irreducibly
  nonlinear least-action escape path: the *linear* instanton is an angle-less
  unstable focus (so it carries no preferred direction — which itself derives two
  of that chapter's empirical findings). The number is genuinely numerical.
  → `RESONATOR_INSTANTON_ATTEMPT.md`.
- **`κ = 2π²` (#2).** Mechanism derived (the funnel is a logarithmic spiral;
  growth is anti-damping, *not* the conservative envelope — the envelope is ~10×
  too small). But the exact constant is a delicate parabolic-cylinder connection
  coefficient (execution, not insight), and the numerical test showed `κ` actually
  **drifts** — it is `2π²(1+O(μ))`, not a sharp constant. The numerics also
  confirmed the rotation half cleanly (phase exact to ~1%). → `MMO_K2_ROUTE_AB.md`,
  `MMO_K2_GINT_MEASUREMENT.md` (+ `mmo_k2_gint.py`).

---

## What ranks as most defensible (revised after the prior-art check)

Reordered honestly. Earlier this section led with the Airy duality and the
circle-map "criticality"; both were overstated and are demoted/withdrawn below.

1. **The `C_q` collapse + fold-escape unification (#1/#12/#16).** A clean,
   parameter-free closed number (`C_q = √(4π ln2)`) and a real unification of three
   thresholds into `σ_crit = √(4π D ln2 ε)` — *within* the standard Freidlin–
   Wentzell framework (Berglund–Gentz is the prior art). High value, modest
   novelty, no overlap problem beyond the framework. The most defensible single
   result, and the one I'd lead with.
2. **The `κ` negative result (#2).** Ruling out the conservative envelope (it is
   ~10× too small *and* the wrong scaling), so the rotation growth is anti-damping.
   Honest and useful — the mechanism, not the constant (and the constant `2π²`
   actually drifts).
3. **The `a_min = a_max e^{−π²f}` "red-herring" dissolution (#3).** Showing there
   is no separate inner exponent to chase. Useful clarification — but the
   underlying SAO-count/global-return split is BGK's territory.
4. **The meta-pattern (next section).** The most transferable insight here.

**Demoted / withdrawn:**

- The **Airy duality (#8)** — a one-to-two-step consequence of two classical facts
  (Riccati→Airy; the adjoint of a log-derivative equation is its square). Elegant
  and unifying for *these two chapters*, but the natural structure to an expert,
  and it lives in Sacré & Franci's active area (check Sacré's 2013 thesis). An
  organizing lemma, **not** a surprising headline; `c ≈ 1.55` is still open.
- The **circle-map "criticality" narrative (#3/#4/#6)** — **withdrawn**: it
  contradicts the project's own Phase-1.5 retraction, the measured `α ≈ 1.55` is
  *below* the critical range `[2,3]` (evidence against criticality), `α=1+1/(1+r)`
  is the existing `α=1+1/p` relabelled, and the whole noise-MMO content overlaps
  the published BGK programme. Only the `Λ(μ)` consolidation survives.
- The **`τ_v` 2-D transfer fix (#13)** — correct and useful, but least novel: it is
  using the right standard tool (the Lyapunov equation) instead of the wrong one.

---

## The meta-pattern

Every idea that *worked* had the same shape: **replace a hand-wavy placeholder
with the exact object, which then collapses cleanly.** The fold-escape barrier,
the `τ_v` Lyapunov transfer, the Airy duality — all that move. The things that
*resist* (the resonator angle, the folded-node global return `Λ`) are exactly the
ones where no such collapsing object exists: they are global nonlinear quantities,
which can be *organized* (the circle map) but not reduced to a universal number.
So "needs authentic creativity / stays open" ≈ "the find-the-exact-object trick
fails here."

---

## Status at a glance

| # | Item | Result | Status |
|---|---|---|---|
| 1 | canard `C_q` | `√(4π ln2) ≈ 2.95` | leading order, certified |
| 12 | excitable `C` | `K_fold=4/3, A₀=1/π`; `σ_crit=√(4π g ln2 ε)` | leading order |
| 16 | VdP universal half | same `√(4π ln2)` | rides #1 |
| 13 | `τ_v` | 2-D Lyapunov `≈0.95` (fixes 1.73→1.16) | leading order |
| 8 | tonic `c` | iPRC `= Ai(−ξ)²`, canard `=−u'/u` (duality) | structure exact; number open |
| 2 | `κ=2π²` | mechanism derived; drifts `2π²(1+O(μ))` | mechanism done; constant = PCF |
| 3 | `a_min` | `= a_max e^{−π²f}`; reduces to `Λ(μ)` | reduced to one integral |
| 6 | `α` | `= 1+1/(1+r)`, `Λ∼μ^{−r}`; data 1.55–1.59 | reduced; exponent transition-regime |
| 4 | `σ_pq` μ-law | circle-map noise-unlocking; `q`-law `q^{−α/γ}` | q-law derived; μ-law = map |
| 10 | resonator angle | nonlinear heteroclinic; linear = angle-less focus | boundary (numerical) |

---

## Index of documents (`regime-tests/`)

- `DERIVATION_AUDIT.md` — the 16-gap inventory (the starting point).
- `DERIVATION_STRATEGIES.md` — strategies for `C_q` (#1), `κ` (#2), `τ_v` (#13).
- `FOLD_ESCAPE_PREFACTORS.md` — `C_q` certified; excitable `C` derived; the fold-escape unification.
- `TONIC_CMID_AIRY.md` — the tonic–canard Airy duality (most novel).
- `MMO_K2_ROUTE_AB.md` — `κ`: antidamping reduction (Route B) + PCF connection setup (Route A).
- `MMO_K2_GINT_MEASUREMENT.md` (+ `mmo_k2_gint.py`) — numerical test of the `κ` reduction (phase exact; `κ` drifts).
- `RESONATOR_INSTANTON_ATTEMPT.md` — `−1.06 rad` (#10): the nonlinear-heteroclinic boundary.
- `MMO_AMIN_MULAW.md` — `a_min` (#3) & `σ_pq` (#4): transition-regime diagnosis, the keystone structure.
- `MMO_AMIN_BRAINSTORM.md` — 15 cross-disciplinary angles on `a_min`.
- `MMO_AMIN_TIER1.md` — exponential asymptotics: `a_min = a_max e^{−π²f}`.
- `MMO_AMIN_TIER2.md` — circle-map / Arnold-tongue framing of the global side.
- `MMO_CIRCLEMAP.md` — the explicit return circle map and the cluster reduction.

---

## The remaining frontier

Small and specific:

1. **One MMO integral `Λ(μ) = ln(a_max/a_min)`** — the global-return slow-passage
   contraction. Its transition-regime scaling `r` sets `α` (#6), `f`, and `a_min`
   (#3) together. *Decisive next step:* sweep `ε` to test whether `Λ` depends on
   `μ` alone or on `μ/ε^{1/2}`; that pins `r` and `α`.
2. **The `κ` connection coefficient (#2)** — the exact `2π²` (μ→0 limit), a
   parabolic-cylinder execution best done from the DLMF/Wechselberger text.
3. **The resonator angle (#10)** — accept as numerical, or characterize the
   normal-form Hamiltonian heteroclinic semi-analytically.

Everything else is either closed at leading order or a write-up/rigour step
(the exact non-Kramers prefactor for the fold-escape family; the matching for the
tonic `c`).
