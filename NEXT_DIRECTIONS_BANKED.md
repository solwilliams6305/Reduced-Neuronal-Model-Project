# Banked next directions (post-session, 2026-06-08)

*Snapshot after the session that cleared most of `NEXT_DIRECTIONS_POST_TW.md`, plus the
directions parked for later. #1 (QIF/θ-neuron grounding of the phase edge) is being developed
now (`QIF_PHASE_EDGE.md`); 2–5 are banked here.*

---

## Status vs `NEXT_DIRECTIONS_POST_TW.md`

- **A (prefactor) — ~80% done.** Frisch–Lloyd current + closed `H = K·η²/4π` (`K≈1.039`) + tail
  reconciliation (exponential SAO tail, not the polynomial hazard) — `frisch_lloyd_current.py`,
  `close_H_matching.py`, `tail_reconciliation.py`. **Remaining:** the analytic Eyring–Kramers
  constant via the vanishing-barrier caustic (banked as #2).
- **B (second edge) — done++.** Quartic anharmonic FPT operator; closed-form rate constant
  `J=(√π/3)·12^{1/6}·Γ(1/6)≈4.976`; `σ^{2/3}` verified; the two-edge statement; the off-critical
  one-parameter crossover; the two-channel synthesis. `phase_edge_snic.py`,
  `offcritical_phase_edge.py`, `FOLDED_CYCLE_PHASE_EDGE.md`, `TWO_CHANNEL_SYNTHESIS.md`.
  **Paper-ready.**
- **C (Airy process) — resolved (two-part verdict).** *Yes* for **successive** peel-offs of one
  cycle (= Airy₂ point process, matched to the GUE edge to ~1%, `airy_line_ensemble.py`); *no*
  for **coupled** arrays (local exponential shadow, not the process — `coupled_array_peeloff.py`,
  `airy_covariance.py`). Full arc with proved/numeric/refuted/retracted tags in
  `FOLDED_CYCLE_AIRY_PROCESS.md`.
- **D (folded-node universality) — done++.** Verified (Weber→SAO → TW_β) + the perturbed-SAO
  edge lemma. `folded_node_universality.py`, `FOLDED_NODE_UNIVERSALITY.md`,
  `FOLDED_NODE_EDGE_LEMMA.md`. **Paper-ready.**
- **E — untouched.** Banked as #3.
- **Side quests:** amplitude `c=5.4439` closed form still open (a Γ(1/6) closed form was found for
  the *phase* constant `J`, not this one); SISR untouched.

## The real bottleneck now: two paper-ready write-ups
- **Second-edge / two-edge paper:** B + the off-critical crossover + the two-channel synthesis.
- **Universality-class paper:** D + the perturbed-SAO edge lemma.

---

## Banked directions (pick up later)

### 2. Finish A — the analytic prefactor
Numerical prefactor pinned (`K≈1.039`). Remaining theorem: the rigorous Eyring–Kramers constant
`A₀` via the **uniform treatment at the vanishing-barrier caustic `Y→0`** (instanton-merging
meets the Airy edge) — the one genuinely-new spot
(`FOLDED_CYCLE_NOISYAIRY_CONJECTURES.md` §C4 / §5). Upgrades the main paper from "exponential
order" to "with the sharp prefactor." In-wheelhouse (Cole–Hopf / SAO / Gelfand–Yaglom).

### 3. E — stochastic blow-up of the folded limit cycle
The one **theorem-grade, structurally-novel** target (`REAL_CONTRIBUTION_CANDIDATES.md` #1): noise
on JKK's 2024 deterministic blow-up; prove the leading-order **stochastic passage through one
chart** of the blown-up singularity (sample-path concentration / exit law). Toolkit-perfect
(blow-up + Freidlin–Wentzell + oscillations). **Gated on the Popović pipeline question:** is a
noise extension of the JKK folded-limit-cycle blow-up in Kuehn's pipeline? (Berglund + Kuehn
active in fold-plus-noise.) The one direction that yields a genuine new *theorem*, not an
extension of the edge story.

### 4. The true Airy *process* via a spectral-interaction model
C's lesson: the Airy₂ **process** needs all-to-all **Dyson** eigenvalue repulsion, which coupling
separate oscillators does **not** supply (nearest-neighbour → exponential, local; the
single-operator point process IS Airy₂, `airy_line_ensemble.py`). The right model is a
**many-mode / spatially-extended folded system** — a **folded PDE**, or a folded system whose
linearisation is a large band/random matrix — where the modes genuinely repel. *"Airy process
from a folded PDE"* is genuinely new, PhD-arc-sized, the real home for the headline.

### 5. The two-edge dichotomy as an observable
Fold-of-cycles vs SNIC oscillators carry different **noise fingerprints** — the TW amplitude edge
(`σ_*∼√ε₂`, peel-off *level* law) vs the quartic phase edge (`σ^{2/3}`, ISI *time* law). A clean
classifier / experimental distinction. Ties directly to #1 (the QIF/θ-neuron is the SNIC member;
a fold-of-cycles bursting neuron is the amplitude member).

---

## Done

### 1. QIF / θ-neuron grounding of the phase edge — **done.**
The SNIC normal form **is** the quadratic-integrate-and-fire / theta neuron; the phase-edge
results **are** the near-threshold inter-spike-interval statistics and the noise-smoothed f–I
curve of a noisy Type-I neuron. Three predictions delivered and tagged: (i) noise-smoothed
`r=σ^{2/3}Φ(ν)` f–I curve with `√I/π` (supra) and Arrhenius `e^{−(8/3)|ν|^{3/2}}` (sub) limits;
(ii) ISI law = the **quartic** FPT law — provably **not** Tracy–Widom, with a `Γ(1/6)` closed-form
mean-ISI constant; (iii) `CV: 1→0.57→0` regularity crossover in `ν=I/σ^{4/3}`. The `σ^{2/3}`
collapse and f–I match confirmed across two noise levels to ~2% / ~3–8%. The two-edge **neural
fingerprint** (SNIC/phase = quartic/`σ^{2/3}` vs fold-of-cycles/amplitude = `TW_β`/`√ε₂`) is the
classifier. Write-up + numerics in `QIF_PHASE_EDGE.md` / `qif_firing_statistics.py` /
`figures/qif_firing_statistics.png`. Ties to banked #5 (the two-edge observable).

### 5. The two-edge dichotomy as an observable — **done** (this session, "D").
Noise-only classifier built: SNIC/Type-I rate `~σ^{2/3}` (timing/phase edge) vs Hopf/Type-II
amplitude `~σ^{1/2}` (amplitude edge), both within 1.5% of prediction; plus the Hodgkin onset
dichotomy (continuous `f→0` vs jump to `f_H`) and the fingerprint table. `two_edge_classifier.py`,
`TWO_EDGE_CLASSIFIER.md`, `figures/two_edge_classifier.png`. (TW stays the fold-of-cycles canard
fine structure, not generic Hopf = Rayleigh.)

### Also done this session (neuro round)
- **QIF data validation.** Synthetic ground-truth harness + ready-to-run AllenSDK pull; estimators
  recover the `σ^{2/3}` collapse / `0.57` CV / quartic-FPT ISI, and the two confounds (colored noise,
  adaptation) leave **opposite, separable** fingerprints (threshold CV `>0.57` ⇒ colored, `<0.57` ⇒
  adapting). `qif_validation.py`, `allen_qif_validation.py`, `QIF_DATA_VALIDATION.md`.
- **Coherence resonance — the ¾ law.** Regularity-optimising noise of a SNIC neuron obeys
  `σ* ∝ |I−I_c|^{3/4}` (CV collapse under `ν=I/σ^{4/3}`; log-log slope `0.78 → 0.75` as finite-size
  removed). Type-I **regularises to the `0.57` plateau** rather than showing a sharp CR minimum —
  the minimum is a Type-II feature (ties to #5). `coherence_resonance_snic.py`,
  `COHERENCE_RESONANCE_SNIC.md`.

---

## Banked — new ideas (neuro round, 2026-06-08), filtered for neuroscience connection

### 6. The Pearcey edge of a cusped fold — at the Bogdanov–Takens point
*[most novel; higher risk]* Merge two folds into a **cusp** (A₃): the inner kernel should leave Airy
and become **Pearcey / higher-Painlevé**, with a new noise exponent replacing `σ^{2/3}`. Neuro home:
the **Bogdanov–Takens** point — the organising centre where Type-I (SNIC) meets Type-II (Hopf)
excitability; "the noise statistics exactly at the Type-I/II boundary." Hard part: derive the inner
equation at the cusp of the critical manifold and identify the Pearcey structure + exponent. Turns
the edge story into a **hierarchy** (catastrophe order → edge universality class).

### 7. A BBP detachment transition in canard escape
*[sharp; medium risk]* A deterministic bias at the fold = a finite-rank "spike" on the SAO; past a
critical strength the top eigenvalue **detaches** from the TW edge and goes Gaussian (Baik–Ben
Arous–Péché). Predict a critical detuning `ν_c` where escape stops being TW-governed → Gaussian; the
neuro reading is **signal detection / mode-locking** (coherent input capturing spike timing vs
noise). Uses the off-critical `ν` already parametrised.

### 8. The exact peel-off law via Painlevé II
*[lowest risk; paper-strengthener]* The amplitude edge is integrable: pin the full TW/Painlevé II CDF
and the explicit tails (upper `~exp(−(4/3)s^{3/2})`, lower `~exp(−|s|^3/12)`) for the canard-escape
amplitude. Upgrades "it's TW" to "here is the distribution and both tails." Distinct from #2 (the
Kramers *prefactor*) — this is the amplitude *distribution*. Mostly known structure, just apply it.

### 9. MMO small-oscillation counts = the Airy line ensemble
*[concrete dynamical observable]* C's capstone: successive peel-offs of one cycle = Airy₂ point
process. In a folded **node** those peel-offs are the small oscillations of a **mixed-mode
oscillation**, so the noisy MMO signature `L^s` is predicted by the line ensemble. Recordable:
entorhinal stellate (grid-cell input), dopaminergic, pre-Bötzinger respiratory neurons; ties to
Desroches–Wechselberger MMOs + Berglund–Gentz noise.

### 10. Heavy-tailed (Lévy) noise: where TW breaks
*[universality boundary]* Swap Gaussian inner noise for α-stable; extreme-value theory says light
tails → TW, heavy tails → **Fréchet**. Find the α at which the canard edge stops being TW — "the
edge of the edge." Heavy-tailed synaptic input is real (avalanches, non-Gaussian bombardment).
