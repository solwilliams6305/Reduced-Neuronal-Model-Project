# Attacking the final 4% of the cusp rung (B) — a strategy

_June 2026. Strategic map (thinking, not closure) for the last 4% of rung B, after the determinantal route was
falsified and 𝒲 was reframed as a first-passage law. One feasibility probe run (`persistence_scaling.py`).
Tags **[DERIVED]/[NUMERIC]/[CONJECTURAL]**; ⚑ = pivotal._

---

## 0. The 4%, honestly decomposed

| piece | what | size | character |
|---|---|---|---|
| **1. Intrinsic process** | the genuine Dyson–Weber multi-point object (own covariance, not a surrogate) | ~2% | deep frontier, shared with rungs A & D; **may resolve as a negative** |
| **2. Persistence exponent** | the β-dependent left-tail exponent of the stochastic Weber field, computed | ~1% | bounded research; **anomalous β-scaling (probe below) ⇒ no naive closed form** |
| **3. Canonical definition** | give 𝒲 a defining equation (renewal/Feynman–Kac), like TW = Painlevé II | the ceiling | **most achievable**; turns "no closed form" into a positive definition |

A blunt up-front truth: part of the 4% is a **proven negative** (𝒲 has no single closed form — not a Fredholm
gap, not a σ-form, not a clean Gamow expansion). That sliver cannot be "calculated"; the right move is **Piece 3**
(define 𝒲 by its equation), not to keep hunting a closed form.

---

## 1. The unifying reframe — attack the *field*, not the marginal ⚑

All three pieces are properties of the **zeros of the stochastic Weber field** $u$ ($u''=(\operatorname{sign}(Y)Y^2-\eta\dot W)u$,
$p=-u'/u$ explodes at a zero):
- 𝒲 = law of the **first zero** of $u$;
- left tail = **persistence** (no zero up to depth $\Theta$) of $u$;
- process = the **multi-zero** structure of $u$.

So the natural "named law" is *the first-zero law of the stochastic Weber field*, and the complete object is the
**field $u$** itself. This reframe is the through-line: stop chasing transforms of the marginal, characterize the
field. **[DERIVED reframe.]**

---

## 2. Piece 3 — the canonical defining equation (do FIRST; highest confidence)

**Idea.** TW isn't "elementary" either — it's *defined* by Painlevé II. Give 𝒲 the same status: write the
**exact equation it already satisfies** as its definition. 𝒲's CDF $F(s)=\mathbb P(Y^\star<s)$ solves the
backward-Kolmogorov / Feynman–Kac PDE for the escape of the inner Riccati,
$$\partial_Y F + (\operatorname{sign}(Y)Y^2-p^2)\,\partial_p F + \tfrac{\eta^2}{2}\partial_{pp}F=0,\quad
F(\text{escaped})=1,$$
plus the recessive initial data at $Y\to+\infty$. Reduce/normalize this to a **self-contained functional/renewal
equation** for $\hat{\mathcal W}$ and *declare it the definition* of $\mathcal W_\beta$.
- **Cheap check:** the equation's numerical solution = the FP-𝒲 we trust (it already does — it *is* the FP
  operator). So this is bookkeeping, not a gamble.
- **Wall risk:** low. The PDE exists; the work is writing it in canonical, reduced, β-explicit form and proving
  self-containedness.
- **Outcome:** 𝒲 becomes a **defined named law** — the honest "completion" given the no-closed-form theorem.
  Likely **+1%**. **[CONJECTURAL-in-form, DERIVED-in-existence.]**

---

## 3. Piece 2 — the persistence exponent (anomalous; needs real machinery)

**Probe result [NUMERIC ⚑].** `persistence_scaling.py`: the left-tail survival rate in the phase-depth
$\Theta=Y^{\star2}/2$ scales **$\sim\beta^{0.6}$** (β=2,4,8 → 3.19, 4.66, 7.49; the β=16 point under-sampled,
discarded). This is **anomalous** — between $\sqrt\beta$ and $\beta$, and crucially **not** the drift-diffusion
$\beta^1$. So the easy route (treat the Prüfer phase as a drifted Brownian motion, persistence = its first-passage,
exponent $\propto\mu^2/D\propto\beta$) **fails**. The persistence is genuinely non-trivial.

**The real route.** The stochastic Weber field is a **non-stationary** Gaussian-type process (instantaneous
frequency $|Y|$ grows). Persistence theory (Bray–Majumdar–Schehr) is built for **stationary** processes, so:
1. **Stationarize.** In the phase variable, rescale to constant frequency; the log-amplitude / Prüfer process
   becomes (approximately) stationary with a computable correlator $C(\tau)$.
2. **Apply the machinery.** Rice formula (zero-crossing rate), the **independent-interval approximation**, or
   exact results for the specific $C(\tau)$, give the persistence exponent $\theta_p(\beta)$.
- **Cheap checks:** measure the zero-crossing correlation of $u$ directly; test IIA against the measured survival;
  the $\beta^{0.6}$ anomaly is the target number.
- **Wall risk:** medium–high. Non-stationary persistence is genuinely hard; realistic outcome is the exponent
  **computed numerically + identified as a (non-universal, β-family) persistence exponent**, possibly not in
  closed form. But the *identification* (left tail = persistence of the stochastic Weber field, anomalous
  $\beta^{0.6}$) is already a real reframe. Likely **+0.5–1%**. **[NUMERIC + CITED route.]**

---

## 4. Piece 1 — the intrinsic process (deepest; run a GO/NO-GO first) ⚑

**The catch.** The unbounded-below operator means this is **not** a standard Dyson-gas / soft-edge process — a
Dyson gas in the cusp potential **falls down the $-Y^2$ side** (no stable top edge). So the "process" is a
**multi-first-passage**, and the open question is whether it carries **intrinsic rigidity** (genuine repulsion,
sub-linear variance growth, like Airy₂) or only **shared-noise correlation** (effectively iid escapes).

**GO/NO-GO (cheap, pivotal).** Simulate the genuine multi-line / successive escapes (no imposed-OU surrogate)
and measure:
- the **variance of the $k$-th escape** vs $k$ — *sub-linear* ⇒ rigidity (GO); *linear* ⇒ iid (NO-GO);
- the **connected 2-point** $\mathrm{Cov}(Y^\star_i,Y^\star_j)$ beyond the shared-noise part.

**Outcomes, both informative:**
- **GO** (rigidity): characterize the covariance — *that* is the intrinsic cusp process (the B/D/A lever). Deep,
  high-value.
- **NO-GO** (iid): a **clean negative** — the cusp has *no* nontrivial intrinsic process; the marginal 𝒲 is the
  whole story (unlike the fold, whose Airy₂ process is genuinely rigid). This would *complete* the process
  question as a theorem-shaped result, and explain why every "process" result needed a surrogate.
- **Wall risk:** the GO/NO-GO itself is cheap; the *characterization* (if GO) is the deep part. **[NUMERIC GO/NO-GO, then DERIVED/CONJECTURAL.]**

---

## 5. Honest assessment — how far the 4% realistically goes

- **Piece 3** (definition): achievable, **~+1%** — the honest completion given the no-closed-form theorem.
- **Piece 2** (persistence): identification + anomalous $\beta^{0.6}$ scaling already in hand (**~+0.5%**); the
  exponent value via non-stationary theory is a hard bounded calc (the rest of the ~1%, possibly numeric-only).
- **Piece 1** (process): GO/NO-GO is cheap and **either answer advances B** — a NO-GO closes the process question
  as a negative (~+1%); a GO opens a deep characterization (the genuinely hard, shared frontier).
- **The proven ceiling** (no closed-form law) is **not closable** — it is a theorem; Piece 3 is its positive
  replacement.

**Realistic ceiling:** B → **~98–99%** with Piece 3 + Piece-2 identification + Piece-1 GO/NO-GO. The final ~1%
is either a deep process characterization (if GO) or is already closed-as-negative (if NO-GO). **100% in the
closed-form sense is provably unreachable** — which is exactly why B is flagged "honest ceiling".

---

## 6. Recommended sequence

1. **Piece 3** — write 𝒲's canonical defining equation (cheap, high-confidence, the honest completion).
2. **Piece 1 GO/NO-GO** — the pivotal cheap simulation: does intrinsic rigidity exist? (Shapes everything; shared with A & D.)
3. **Piece 2** — the non-stationary persistence-exponent calculation (the hard bounded research; the $\beta^{0.6}$ anomaly is the target).
4. **If Piece 1 = GO:** characterize the intrinsic covariance (the deep prize).

**Through-line:** all three are facets of one object — *the zeros of the stochastic Weber field*. Define it
(3), test its multi-zero rigidity (1), compute its persistence (2). That is the complete, honest closure of B —
with the standing caveat that the closed-form law provably does not exist.
