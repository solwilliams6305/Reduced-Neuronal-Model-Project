# MMO — deriving the *global* μ-dependence of σ_pq

> **⚠ Phase M2 executed — deterministic (`mmo_sigmapq_m2.py`, `MMO_NOISE.md` §4.2).**
> Measured a_min(c) and the sensitivities with NO noise (no counting floor). Result
> refines this plan: **the μ·a_min candidate is NOT empirically preferred over
> μ^{3/2} at current resolution.** With the plateau-*centre* a_min, σ_pq/(μ·a_min) =
> 0.87, 1.70 at q=3,4 (NOT the ~1.9 flat of §1's 2-pt note), and μ·a_min
> *underpredicts* the measured σ_pq(q=4) by ~2× while μ^{3/2} matches it. The cause
> is the chief hazard §6.1/§6.2, now quantified: **a_min is factor-~2 ambiguous** —
> it varies smoothly across each plateau's c-width (≈0.54→0.25 across q=3) and with
> the counting protocol, and that ambiguity **exceeds** the μ·a_min-vs-μ^{3/2}
> difference. Also κ=ln(a_max/a_min)/(sμ) → 2π² only asymptotically (8 at q=3 → 18
> at q=6). **Conclusion:** the deterministic side cannot separate the forms; the
> candidates are genuinely degenerate at q=3,4 and only diverge at q=2,5,6 (the M1
> target). Separation requires BOTH M1 (robust q=2…6 noise σ_pq, the counting-floor
> study) AND a **derived** a_min(c) — i.e. **Phase B is the real resolver**, exactly
> the shared-unknown coupling this plan identifies. "Do not declare the power"
> stands; μ^{3/2} remains a placeholder, μ·a_min mechanistically motivated but unproven.

**Goal.** The crossover (`MMO_CROSSOVER.md`) showed the **local** noise escape in
the folded-node funnel is canard-like (μ-independent prefactor C_q, no μ^{3/2}).
So the μ-dependence in the measured `σ_pq ≈ C_q√ε·μ(q)^{3/2}` is a **global**
property — how the global return + funnel-filling map the plateau index q onto
the injection geometry. This plan derives it from the deterministic funnel,
turning the empirical "μ^{3/2}" into a mechanism (and, very likely, into the
deterministic Path-A quantity `a_min(c)`).

**Why it matters.** It is the last non-trivial open piece of the noise side, and
it appears to be the **same unknown** as the deterministic Path-A item
(closed-form `a_min(c)`). If so, one calculation closes two open ends and the
noise side reaches the same first-principles standard as the deterministic side.

---

## 1. The target (and how soft the current "μ^{3/2}" is)

The dissolution threshold was matched with `σ_pq ≈ C_q√ε·μ(c_pq)^{3/2}`, C_q ≈ 8–10,
giving a q-slope β ≈ 1.13 = measured 1.11. **Two honest caveats:**

- The "μ^{3/2}" was *assumed* (inherited from the now-refuted local argument) and
  fit the q-slope; it is **not** independently pinned. From the two cleanly
  resolved plateaus the direct σ_pq-vs-μ power is ~2.4, not 1.5 — the q-slope
  match survives only because μ(q) is not a clean power of q. So the real target
  is "derive the σ_pq(μ) law," and the answer may **not** be exactly μ^{3/2}.
- A mechanistic candidate, `σ_pq ~ μ·a_min(c)` (§2), is *motivated but not
  empirically preferred*. An edge-`a_min` 2-point note first suggested
  `σ_pq/(μ·a_min) ≈ 1.9` (flat), but **Phase M2 (banner) overturned it**: with the
  plateau-*centre* `a_min` the ratio is 0.87, 1.70, and the factor-~2 `a_min`
  ambiguity exceeds the `μ·a_min`-vs-`μ^{3/2}` gap, so the forms are **degenerate**
  at q = 3, 4 (μ^{3/2} actually matches q=4 slightly better there). The target
  stands: derive σ_pq(μ) — which the deterministic side alone cannot pin (M2).

**Deliverable:** the deterministic law σ_pq(μ, c) — its functional form, the
μ/c-power, and the O(1) prefactor's relation to the canard-escape C_q — validated
against a *better-resolved* σ_pq sweep.

---

## 2. The decomposition (the candidate derivation chain)

A plateau ρ = 1/q (s = q−1 SAOs) dissolves when noise makes the realised SAO
count fluctuate by O(1): δs ~ 1. Decompose σ_pq into deterministic sensitivities:

**(i) Count ← injection (the rotation map, measured).** From the K2 rotation map
`a_n ~ a_in·e^{κμ n}`, the count from injection a_in up to the jump amplitude
a_max is `s ≈ ln(a_max/a_in)/(κμ)`. Hence

```
ds / d(ln a_in)  =  −1/(κμ)      ⇒   to change s by 1:  δ(ln a_in) ~ κμ.
```

This is the first μ-handle: **smaller μ ⇒ a smaller log-amplitude shift flips the
count ⇒ easier to dissolve ⇒ lower σ_pq at high q** — the right sign for
high-q-first, from the *measured* κ ≈ 2π².

**(ii) Injection jitter ← noise (the local escape, from the crossover).** The
noise jitters where the global return injects the trajectory, by δ(ln a_in) =
δa_in / a_in. The crossover fixes the noise→displacement transfer as the
**canard-escape** one (μ-independent, prefactor C_q): a physical displacement
δa_in ~ (C_q-transfer)·σ over the injection passage. The injection amplitude is
a_in = **a_min(c)** (the deepest loop). So δ(ln a_in) ~ C_q·σ / a_min(c) (× the
√ε / accumulation bookkeeping to be pinned).

**(iii) Compose.** Setting δ(ln a_in) ~ κμ (dissolution):

```
σ_pq  ~  (1/C_q-transfer) · κ · μ · a_min(c).            [candidate law]
```

So the candidate is **σ_pq ~ μ · a_min(c)** (up to κ ≈ 2π² and the canard-escape
O(1)), i.e. the global μ-dependence is **the rotation-map factor μ times the
injection depth a_min(c)**. This is mechanistically motivated, but **M2 (banner)
shows it is not yet distinguishable from the placeholder μ^{3/2}** at accessible
resolution — the a_min ambiguity dominates. Separating them is the open problem
(§5), not a settled reinterpretation.

---

## 3. The key coupling — a_min(c) is the deterministic Path-A unknown

The c-dependence in the candidate law lives almost entirely in **a_min(c)**, the
smallest-SAO amplitude — which is exactly the quantity the deterministic side
left open (`MMO_K2.md`: a_min(c) measured numerically, 0.80 → 0.008 across the
band; closed form is the Path-A K2 inner solution). **Therefore deriving the
global μ-dependence of σ_pq reduces to deriving a_min(c).** One calculation closes
both:

- deterministic Path-A: closed-form a_min(c) ⇒ closed-form funnel-filling f(c), α;
- noise global μ-power: closed-form a_min(c) ⇒ closed-form σ_pq(μ, c).

This is the plan's central structural claim and the reason it is worth doing now.

---

## 4. Task sequence

### Phase M — measure first (the empirical law, cleanly)

The current σ_pq(μ) rests on 2–3 plateaus; the power is not pinned. This phase
overlaps the deferred Task-β counting study — do them together.

**M1.** Resolve σ_pq on q = 2…6 with the converged counter (dt = 0.01, per-plateau
`a_min`-scaled prominence H; `mmo_noise_staircase.py`). Fit σ_pq(μ) and σ_pq(q)
with error bars.
**M2.** Measure the two sensitivities in (2) *directly* and deterministically (no
noise): (i) `ds/d(ln a_in)` by perturbing the injection amplitude and counting
loops — check it equals −1/(κμ); (ii) `a_min(c)` from the deterministic episodes
(already in `mmo_k2_return_map.py`).
**M3.** Test the candidate `σ_pq ~ μ·a_min(c)` vs the fitted `μ^{3/2}`: which
collapses the q = 2…6 data? Also measure the direct rotation-spread γ
(std(ρ) ~ σ^γ) for the §5.2-of-chapter γ slot.

*Deliverable:* the empirical σ_pq(μ, c) law with enough plateaus to distinguish
`μ·a_min(c)` from `μ^{3/2}` — the target the derivation must hit.

### Phase A — Path-B semi-analytic (close the mechanism)

**A1.** Assemble the candidate law from measured ingredients: κ ≈ 2π² (have it),
a_min(c) (numerical, M2), and the canard-escape transfer C_q (crossover). Predict
σ_pq(μ, c) and compare to Phase M.
**A2.** Pin the √ε / accumulation bookkeeping in step (ii) (the injection passage
time τ_inj and how δa_in accumulates) so the prefactor — not just the power —
matches, and check it lands at C_q ≈ 8–10. This is the noise analogue of the
deterministic Path-B that closed f(c).

*Deliverable:* σ_pq(μ, c) reproduced from measured funnel geometry + the
canard-escape constant — the mechanism closed, with a_min(c) still numerical.

### Phase B — Path-A first-principles (close a_min(c) ⇒ close both)

**The literature ingredient (found — Desroches et al. 2012 SIAM Review §3,
citing Brøns–Krupa–Wechselberger 2006 [Thm 3.1] and Wechselberger 2005).** In
the standard folded-node scaling `x = ε^{1/2}x̄, y = εȳ, z = ε^{1/2}z̄`, the
funnel is **O(ε^{1/2})** in physical (v) units. With `2k+1 < μ^{-1} < 2k+3`
there are k secondary canards ξ_i (the i-th twists 2i+1 times); max SAOs
s ≈ k+1 ≈ (μ^{-1}+1)/2 [matches our s_max ~ 1/2μ]. **Theorem 3.1:** at O(1) from
the fold, all secondary canards lie within an **O(ε^{(1-μ)/2})** band of the
strong canard, and the rotational sectors I_i (i ≤ k) have width **O(ε^{(1-μ)/2})**
(the maximal-rotation sector I_{k+1} is O(1)). The global return injects at a
distance **δ from the strong canard**, and δ selects the sector (Thm 3.3,
subsector widths O(ε^{(1-μ)/2})).

This is the missing ingredient. Two consequences for our σ_pq:
- **a_min and the sector scale ~ ε^{1/2}-ish, weakly μ-dependent** (ε^{(1-μ)/2});
  the *spacing* between adjacent secondary canards is the band width / k ~
  **μ·ε^{(1-μ)/2}** — that is the **δ-shift to change s by 1**, and the source of
  the μ-handle in σ_pq. (So the candidate's "μ" is the sector-packing factor, the
  literature scaling — not a guessed power.)
- **⚠ Regime caveat (decisive).** Thm 3.2/3.3 assume **μ ≫ ε^{1/2}** (A0). Our
  FHR has ε^{1/2} ≈ 0.283 and **μ ∈ [0.03, 0.16] < ε^{1/2}** across the whole
  band — so we are in the **Krupa–Wechselberger 2010 transition regime
  (μ = O(ε^{1/2}) or smaller)**, where k = O(ε^{-1/2}) and the SAO amplitudes are
  "tiny" (consistent with our a_min → 0.008). **Phase B must use Krupa–
  Wechselberger 2010, not the standard Thm 3.1**, and the sector/amplitude
  scalings above are upper bounds that need the transition-regime correction.

**B1.** From Krupa–Wechselberger 2010, get the transition-regime sector spacing
and the smallest-SAO amplitude a_min as a function of (μ, ε, δ_inj). Combine with
the FHR global-return δ_inj(c) (numerical) to get a_min(c).
**B2.** Feed into the σ_pq chain (§2): the δ-shift to change s by 1 is the
sector spacing; the noise→δ transfer is the canard-escape C_q (crossover).
Get the closed-form σ_pq(μ, c), read off the true μ-power, confirm against
Phase M, and retire "μ^{3/2}".

*Deliverable:* first-principles σ_pq(μ, c) **and** closed-form a_min(c)/f(c) — both
open items closed together. (Verify the normal-form convention, as for the
crossover Phases C–D.)

### Phase C — reconcile and write up

**C1.** State the final σ_pq(μ, c) law and retire the placeholder "μ^{3/2}"
(replace with the derived form, e.g. μ·a_min(c)) across `MMO_NOISE.md`,
`MMO_CHAPTER.md` §4.3/§6/§7, and the README rows + structural finding.
**C2.** Update the chapter outline §4.3/§6; banner line on `MMO_FHR_PLAN.md`.

---

## 5. Success criteria / falsifiability

- **Clean success:** Phase M distinguishes the forms and one of them (likely
  μ·a_min(c)) collapses q = 2…6; Phase A reproduces it from measured geometry with
  C_q ≈ 8–10; Phase B derives a_min(c) and hence the σ_pq μ-power from first
  principles. The noise side is then first-principles, and the deterministic
  a_min(c) is closed as a by-product.
- **Partial (still valuable):** Phases M + A close the *mechanism*
  (σ_pq ~ μ·a_min(c), prefactor C_q) with a_min(c) numerical — the same
  universal-mechanism + numerical-global-return status as the deterministic Path B.
- **Informative failure:** if no clean σ_pq(μ, c) law survives better resolution
  (e.g. σ_pq depends on more than μ and a_min — say on a_max or the full injection
  distribution), then the dissolution is not a single-injection-depth effect; that
  is itself a result (the funnel's noise response is richer than the rotation map)
  and reframes the §4 story.

**Watch:** the current "μ^{3/2}" and the candidate "μ·a_min(c)" are **degenerate**
on the present 2–3 plateaus — only Phase M (more plateaus) can separate them. Do
not declare the power before M.

---

## 6. Hazards

1. **Degeneracy of forms (the chief hazard).** Multiple μ/c combinations fit a
   narrow plateau range; the derivation must be tested on q = 2…6, not 3–4.
2. **a_min counting at high q.** The deepest loops sit at the noise-/dt-countability
   floor (the Phase-3 limiting systematic); a_min(c) itself gets uncertain as
   c → −1. The per-plateau a_min-scaled H (Task β) is the prerequisite.
3. **The √ε / accumulation bookkeeping (A2)** is exactly where the canard chapter's
   own prefactor work was subtle; expect the *power* to come easily and the
   *prefactor* (C_q ≈ 8–10) to need care.
4. **Path-A convention.** a_min(c) from the K2 inner solution needs the verified
   normal-form convention (the same caveat flagged for the crossover Phases C–D).
5. **Don't over-claim μ^{3/2}.** It is a placeholder; the honest deliverable may be
   μ·a_min(c) (or another derived form). Revising the label *is* part of the result.

---

## 7. References

- Wechselberger (2005); Krupa & Wechselberger (2010) — the K2 inner solution /
  a_min(c) (Path B of this plan's Phase B).
- Brøns, Krupa, Wechselberger (2006); Desroches et al. (2012) — secondary-canard
  spacing and the rotation structure.
- Internal: `MMO_K2.md` (κ ≈ 2π², a_min(c) numerical, f(c)), `MMO_NOISE.md`
  (σ_pq data, the exponent reconciliation), `MMO_CROSSOVER.md` (the local escape =
  canard, the reason μ^{3/2} is global), `MMO_CROSSOVER_PLAN.md` (the coupled FSN
  scaling).

---

## 8. Payoff

The deterministic and noise sides share one unknown — **a_min(c)**. Closing it
(Phase B) delivers, in one stroke: closed-form funnel-filling f(c) and α
(deterministic Path-A), and closed-form σ_pq(μ, c) (the noise global μ-power). The
chapter then has no "empirical/derivation-open" entries left — every row in
`MMO_CHAPTER.md` §6 becomes derived. Phases M + A alone already upgrade the noise
μ-dependence from "empirical placeholder μ^{3/2}" to "mechanism (μ·a_min(c)) closed
with the shared canard prefactor," which is the defensible, write-up-ready core.
