# Next directions after the Tracy–Widom paper — a tiered, grounded shortlist

*What this is: the next-step map **conditioned on the noisy-folded-cycle paper as it now
stands** (TW exit law proved to exponential order + θ-uniform, recrossing ceiling verified).
It tiers the open ground by how much is already in hand, and is honest about difficulty,
novelty risk, and competition. Unlike `FUTURE_DIRECTIONS_SOURCED.md` (a broad literature-lead
map), this one is anchored in the result you already have.*

*Epistemic calibration, same as the sourced map: the **in-house gaps** (Tier 1–2) I can
certify against your own confirmed work; the **literature-novelty** claims (Tier 3, "no one
has done X") are search-supported but not certifiable from outside the expert network — treat
those as Popović questions, not green lights.*

---

## Landscape note (why the targeting below is what it is)

The fold-plus-noise corner is being actively worked by the people who built your foundations:

- **Bergeot–Berglund–Zogheib, "The dynamic saddle-node bifurcation with noise on the slow
  variable," arXiv:2512.10460 (Dec 2025, rev. Apr 2026).** Noise on the *slow* variable
  through a fold; the contribution comes out "in terms of **Airy functions**"; results are at
  the level of the **mean and variance** of the slow variable after the bifurcation.
- Kuehn's group owns the deterministic blow-up line (JKK 2024) and the critical-transitions /
  variance-up-to-codim-2 framework.

The key observation: **this nearby activity is moment-level / deterministic-delay-level, not
distributional, and not a universality statement.** Your differentiated lane is the
**edge/universality layer** — the full exit *law*, its sub-exponential *prefactor*, the
*second edge*, and the finer KPZ objects. Recommendation: do **not** go head-to-head on
slow-variable noise (Berglund is there); push the distribution/universality side, which the
moment papers do not touch and which is where your TW result is already unique.

---

## Tier 1 — Finish the result you have

### A. The sharp prefactor / inner exit measure  ★ recommended next
**What:** Turn the TW *identification* into the rigorous **sub-exponential prefactor** — the
full Eyring–Kramers constant `A₀` for the folded-cycle canard escape. This is the single named
gap in `FOLDED_CYCLE_NOISYAIRY_CONJECTURES.md`: the inner exit measure `μ_η(Y_peel,R_peel)`
that corrects the quasi-static hazard `H=η²/4π` in the window `Y ≲ η^{4/3}`.

**Why it's reachable (not open-ended):** Conjecture 1 (inner Riccati = stochastic Airy
operator, `η=2/√β`, exit law = TW_β) is **numerically confirmed** — parameter-free in
mean/std/skew/kurt at β=1,2, raw histograms on the exact Painlevé II densities. So the analytic
task is to *read off* a prefactor that is already a known TW-edge quantity, with three classical
tools all aimed at the same object:
- **Frisch–Lloyd current** (C2) — the most computable; the escape rate as the IDOS/rotation-number
  current of `u''=(Y−ηξ)u`, a closed quadrature.
- **RRV diffusion** (C1) — the rigorous backbone; RRV's diffusion *is* your inner Riccati.
- **Gelfand–Yaglom determinant** (C4) — the physicist's one-loop check.

**The one genuinely new piece:** the **uniform treatment at the vanishing-barrier caustic
`Y→0`**, where instanton-merging meets the Airy edge. By the search done for this note, this
precise junction (fold-canard escape read as a TW edge, with the caustic uniformised by the
Airy edge itself) **does not appear in the literature in this form** — it is the lone spot a
new idea may be required, and the most publishable sub-result.

**Payoff:** upgrades `FOLDED_CYCLE_PATHA_PROOFS.md` from "critical-noise law to exponential
order" to "**with the sharp prefactor**," and cements the **KPZ-edge-class** membership as a
theorem, not a numerical match. **Fit:** perfect — your Cole–Hopf / SAO / FW machinery.
**Ceiling:** high and in-wheelhouse. **Concrete first move:** the Frisch–Lloyd falsifier (see
"First move" below).

---

## Tier 2 — Natural next theorems (reuse what you built)

### B. The phase channel — close the second edge
**What:** You have the amplitude (TW) edge in full; the **SNIC phase channel** is the twin you
have only to leading exponent. At a saddle-node-on-invariant-circle the inner phase equation is
the noisy saddle-node-on-circle with `ω, D_φ ∼ σ^{2/3}` and constant `J=∬_{w<u} e^{(w³−u³)/3}`,
which **"lacks a Tracy–Widom-style closed form"** and has slowly-converging tails (paper,
`thm:phase` + discussion). **Open object:** the *universal first-passage law* of the noisy
saddle-node-on-circle — the phase-edge analogue of your amplitude edge.
**Why it matters:** turns "two channels" into a single **two-edge** statement (an amplitude/TW
edge and a phase/SNIC edge, the destroying bifurcation selecting which is singular) — a clean,
self-contained follow-on paper. **Fit:** your inner-equation + FW toolkit, plus stochastic phase
reduction (Goldobin–Teramae–Nakao–Ermentrout 2010). **Difficulty:** medium; **novelty:** clear.

### D. Folded-node universality — one example into a class
**What:** Your Universality subsection already conjectures the **same TW_β exit law** for the
noisy fold and the noisy **folded node**; BGK's "Hunting French ducks" is the immediate neighbour.
First verify TW_β numerically for the noisy folded-node passage, then prove it by the same
Cole–Hopf → SAO route. **Why it matters:** converts "the folded *cycle* is in the TW edge class"
into "this **family** is" — the stronger, more quotable claim, and the natural way to argue the
result is structural rather than incidental. **Fit:** reuses your machinery almost verbatim.
**Risk:** lowest of the lot (it's an extension, not a new mountain), which is exactly why it's a
safe, fast second paper.

> **TESTED (2026-06-07) — universality holds.** The folded node's canard inner equation is the
> **Weber (parabolic-cylinder)** equation `u''=(z²/4−(ν+½)−ηξ)u` — a genuinely *different* global
> operator from the fold's Airy, with an oscillatory well (the rotations) and a simple turning
> point `z*=2√(ν+½)`, local slope `s=√(ν+½)`. Running the same recessive Cole–Hopf / first-node
> machinery and rescaling to local Airy coordinates `ζ=s^{1/3}(z−z*)`, the noise-broadened
> peel-off **converges to TW_β (β=4s/η²) as the fold sharpens (s→∞)**: at `s≈11`, std `0.907`
> (TW₂ `0.902`), peel-off location `−2.35→−2.338`, skew `0.19` and kurt `0.09` rising
> monotonically to TW₂ `0.224, 0.093` from large finite-curvature (Weber) corrections at small
> `s`. The folded *cycle* (Airy) reference matches TW_β exactly across `β=4,2,1`
> (skew `0.178/0.226/0.310` vs `0.166/0.224/0.293`). **Verdict: the folded node is in the same
> TW_β edge universality class** — the location is model-dependent (Weber-shifted), the edge
> *fluctuation law* is universal. So "the folded cycle is in the TW class" → "this family is."
> Scripts `folded_node_universality.py`, `figures/folded_node_universality.png`.

---

## Tier 3 — Ambitious / structural (higher ceiling, higher cost)

### C. The Airy process for successive / coupled peel-offs — the headline
**What:** Your discussion flags it directly: a *single* peel-off is TW (the **edge**, a 0+1-dim
object); the **joint law of successive peel-offs** in a periodically *forced* or *coupled*
folded cycle should be a finer KPZ-class object — the **Airy₂ process**. **Status:**
search-supported **open** — no work connects Airy processes to canard escape or to coupled
FitzHugh–Nagumo; the nearest deterministic neighbour is coupled-FHN MMOs via cusped-singularity
blow-up (arXiv:2202.12027). **Why it matters:** the "**independent interest beyond JKK**"
result — it would place the *temporal* fluctuation structure of a noisy oscillator in the KPZ
fixed-point picture. **Difficulty:** high; PhD-arc-sized; needs the determinantal / multi-point
machinery and a forced-or-coupled model where successive peel-offs are well-defined. Best as a
thesis spine, not a single paper.

### E. Stochastic blow-up of the folded limit cycle — your "real, even if small" target
**What:** From `REAL_CONTRIBUTION_CANDIDATES.md` #1: take JKK's (2024) deterministic blow-up of
the folded limit cycle, add noise, and prove a **leading-order stochastic passage through one
chart** of the blown-up singularity (sample-path concentration / exit law). **Status:** no
stochastic version exists (search-supported). **Bounded first result is feasible.** **Fit:**
blow-up + Freidlin–Wentzell + oscillations — exactly your stack. **Risk — now sharper:** with
**both Berglund and Kuehn** active in fold-plus-noise (see Landscape note), the pipeline risk
your candidates note flagged is real. **This is the one to clear with Popović first** — the exact
question: *"Is a noise extension of the JKK folded-limit-cycle blow-up in Christian's pipeline? If
not, I'd like the leading-order stochastic passage through one chart."*

---

## Tier 4 — Bounded side quests

- **Closed form for `c = 5.4439`** (the Painlevé II action setting the early-escape rate).
  A clean "is there an elementary closed form" question; bounded, low-stakes, satisfying if it
  cracks.
- **SISR / coherence-resonance sharpening on the FHN model.** It *is* your model and your
  FW/Kramers + numerics toolkit, but — as `FUTURE_DIRECTIONS_SOURCED.md` already warns — the
  risk is sliding into "incremental sharpening" of analytical SISR beyond idealised limits. Pin
  the precise open quantity with Popović before committing.

---

## Ranking

1. **A** — finish the prefactor (completes the paper; mostly scaffolded; one genuinely-new crux).
2. **B** — the phase/second edge (natural next theorem, self-contained).
3. **D** — folded-node universality (hardens the class claim; lowest risk).
4. **C** — Airy process (highest ceiling, thesis-sized bet).
5. **E** — stochastic blow-up (structurally novel; gated on the Popović pipeline question).

---

## First move (needs no new theory — already in the conjecture note)

The **Frisch–Lloyd falsifier** for direction A, straight from
`FOLDED_CYCLE_NOISYAIRY_CONJECTURES.md` §7:

1. Integrate the **stationary Riccati current** for `dR=(R²−Y)dT+η dB` (equivalently the IDOS
   current of `u''=(Y−ηξ)u`) at **fixed Y**.
2. Compare to the quasi-static Kramers rate `(√Y/π) e^{−8Y^{3/2}/3η²}`.
3. Locate where the two diverge — predicted at `Y ∼ η^{4/3}` — which **localises the inner
   correction** and gives the prefactor as a quadrature.

Companion check (tail reconciliation): measure the early-escape fraction vs η and fit against
both `η²/4π` (quasi-static) and the SAO left tail `∝ exp(−s³/6η²)` (Dumaz–Virág); report which
tracks and the prefactor ratio. Both are bare-numpy, fast, and decisive about whether A reduces
to "read off a known edge quantity."

---

## References

**Recent / web-verified (June 2026):**
- Bergeot, Berglund, Zogheib — *The dynamic saddle-node bifurcation with noise on the slow
  variable*, arXiv:2512.10460 (2025/26). [Adjacent: slow-variable noise, Airy functions, moments.]
- Jelbart, Kuehn, Kuntz — *Geometric Blow-Up for Folded Limit Cycle Manifolds in Three
  Time-Scale Systems*, arXiv:2208.01361 (2024). [Deterministic scaffold for E.]
- Dumaz, Virág — *The right tail exponent of the Tracy–Widom-β distribution*, arXiv:1102.4818.
  [SAO tail used in A.]
- *Mixed-mode oscillations in coupled FitzHugh–Nagumo oscillators: blow-up analysis of cusped
  singularities*, arXiv:2202.12027. [Nearest deterministic neighbour for C.]

**Standing (own these cold — see `THEORY_ROADMAP.md`):**
- Ramírez, Rider, Virág, *Beta ensembles, stochastic Airy spectrum, and a diffusion*, J. AMS 24
  (2011). [SAO = your inner Riccati; backbone of A.]
- Berglund, Gentz (–Kuehn), *Noise-Induced Phenomena…* (2006); *Hunting French ducks in a noisy
  environment* (2012/15). [Concentration tubes; the folded-node neighbour for D.]
- Goldobin, Teramae, Nakao, Ermentrout (2010), *Dynamics of limit-cycle oscillators subject to
  general noise*. [Stochastic phase reduction for B.]
- Corwin (2012), *The KPZ equation and universality class*; Quastel–Remenik on Airy processes.
  [The Airy-process target in C.]

## How this maps to your files

A → `FOLDED_CYCLE_NOISYAIRY_CONJECTURES.md`, `FOLDED_CYCLE_PATHA_PROOFS.md`,
`inner_exit_measure.py`, `folded_cycle_noisyairy.py` · B → `thm:phase` in
`NoisyFoldedCycle_paper.tex`, the Channel-B chapters · C → paper discussion (Airy process),
`folded_cycle_tracy_widom.py` · D → Universality subsection, `folded_cycle_normal_form_map.py` ·
E → `REAL_CONTRIBUTION_CANDIDATES.md`, the blow-up chart machinery.
